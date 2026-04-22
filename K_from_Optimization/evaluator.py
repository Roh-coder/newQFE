"""
evaluator.py — single function from (r1, r2) to (cost, sigma_cost).

Each call:
  1. Find beta_c on the test lattice with k1=r1, k2=r2, k3=1
     (warm-started from the previous call's beta_c if available).
  2. Run a longer production MC at that beta_c.
  3. Compute the L^2 cost vs the cached reference correlator.
  4. Append a record to eval_log.jsonl.
  5. Push a frame to OptimizerPlotter (if attached).

The class is stateful only because it caches the last beta_c (warm start)
and re-uses the reference correlator data.  All other state lives on disk.
"""

from __future__ import annotations

import json
import os
import shutil
import time
from dataclasses import dataclass, asdict
from typing import Optional

import mc_engine
import cost as cost_module


@dataclass
class EvalResult:
    eval_id: int
    r1: float
    r2: float
    beta_c: float
    cost: float
    sigma_cost: float
    snr: float
    snr_status: str
    per_dir: list
    per_dir_sigma: list
    n_traj_prod: int
    n_traj_scan_total: int
    wall_time_s: float
    # β_c-scan diagnostics for the criticality-check visualizer.
    scan_betas: list = None
    scan_chis: list = None
    scan_chi_errs: list = None


class Evaluator:
    """Callable evaluator: ``cost, sigma = ev(r1, r2)``.

    Parameters
    ----------
    exe : str
        Path to the C++ simulator (`bin/ising_tri_twisted_parallelogram`).
    ref_data : dict
        Pre-loaded reference correlator (from `mc_engine.load_all_to_all`).
    ref_geom : tuple[int, int, int, int]
        (Lx, Ly, Tx, Ty) of the reference lattice.
    test_geom : tuple[int, int, int, int]
        (Lx, Ly, Tx, Ty) of the test lattice.
    output_dir : str
        Where to write `eval_log.jsonl` and to host `frames/`.
    n_traj_prod : int
        Production-run trajectories per evaluation.
    n_traj_scan_coarse, n_traj_scan_fine : int
        Trajectories for the coarse / fine GC-scan passes.
    beta_seed : tuple[float, float]
        Initial (beta_lo, beta_hi) bracket for the very first call.
    optimizer_plot : visualization.OptimizerPlotter or None
        Optional 4-panel plotter; updated after each evaluation.
    beta_plot_dir : str or None
        If given, a fresh `BetaScanPlotter` is created per call writing into
        this directory (one sub-set of frames per evaluation).
    keep_mc_subdirs : bool
        If False (default), MC scratch subdirs are wiped after each call so
        disk doesn't fill up.  Set True when debugging.
    """

    def __init__(self, exe, ref_data, ref_geom, test_geom, output_dir,
                 n_traj_prod=30000, n_traj_scan_coarse=4000,
                 n_traj_scan_fine=10000,
                 beta_seed=(0.20, 0.32),
                 optimizer_plot=None, beta_plot_dir=None,
                 keep_mc_subdirs=False, mu=0.0):
        self.exe = exe
        self.ref_data = ref_data
        self.ref_geom = tuple(ref_geom)
        self.test_geom = tuple(test_geom)
        self.output_dir = output_dir
        self.n_traj_prod = int(n_traj_prod)
        self.n_traj_scan_coarse = int(n_traj_scan_coarse)
        self.n_traj_scan_fine = int(n_traj_scan_fine)
        self.beta_seed = tuple(beta_seed)
        self.optimizer_plot = optimizer_plot
        self.beta_plot_dir = beta_plot_dir
        self.keep_mc_subdirs = keep_mc_subdirs
        self.mu = float(mu)

        os.makedirs(output_dir, exist_ok=True)
        self.log_path = os.path.join(output_dir, "eval_log.jsonl")
        # Wipe any previous log so each benchmark run is self-contained.
        if os.path.exists(self.log_path):
            os.remove(self.log_path)

        self._mc_root = os.path.join(output_dir, "mc_scratch")
        os.makedirs(self._mc_root, exist_ok=True)

        self._eval_id = 0
        self._beta_prev: Optional[float] = None
        # Set by run_nelder_mead before each evaluator call so OptimizerPlotter
        # can draw the current simplex in the trajectory panel.
        self.current_simplex: Optional[list] = None

    # ------------------------------------------------------------------
    def __call__(self, r1: float, r2: float) -> EvalResult:
        """Run one full evaluation; return the EvalResult."""
        self._eval_id += 1
        eid = self._eval_id
        Lx, Ly, Tx, Ty = self.test_geom
        ref_Lx, ref_Ly, ref_Tx, ref_Ty = self.ref_geom

        t0 = time.time()
        label = f"eval{eid:04d}_r1_{r1:.4f}_r2_{r2:.4f}"
        scratch = os.path.join(self._mc_root, label)
        os.makedirs(scratch, exist_ok=True)

        # ---- 1. find beta_c (warm-started bracket) ----------------------
        # Pass-0 sweeps a WIDE window so we rarely need to translate.
        # Window translations are pure waste (each one adds n_coarse fresh
        # MC points), so it is far cheaper to cover more β-space up front.
        # eps = 20% of beta_prev, floor 0.05 → typical span ≈ 0.10.
        if self._beta_prev is not None:
            eps = max(0.20 * self._beta_prev, 0.05)
            beta_lo = max(0.01, self._beta_prev - eps)
            beta_hi = self._beta_prev + eps
        else:
            beta_lo, beta_hi = self.beta_seed

        beta_cb = None
        if self.beta_plot_dir is not None:
            from visualization import BetaScanPlotter
            beta_cb = BetaScanPlotter(self.beta_plot_dir, label,
                                      every_point=True)

        print(f"[ev {eid:03d}] r1={r1:.4f} r2={r2:.4f}  "
              f"β∈[{beta_lo:.4f},{beta_hi:.4f}]")
        beta_c, chi_peak, sb, sc, sce = mc_engine.find_beta_c(
            self.exe, Lx, Ly, Tx, Ty, r1, r2, 1.0,
            beta_lo, beta_hi,
            n_traj_coarse=self.n_traj_scan_coarse,
            n_traj_fine=self.n_traj_scan_fine,
            data_dir=os.path.join(scratch, "scan"),
            progress_cb=beta_cb,
        )
        n_scan = len(sb) * (self.n_traj_scan_coarse + self.n_traj_scan_fine) // 2
        print(f"[ev {eid:03d}]  β_c={beta_c:.6f}  χ_peak={chi_peak:.3g}  "
              f"({len(sb)} scan pts)")

        # ---- 2. production MC at beta_c --------------------------------
        prod_dir = os.path.join(scratch, "prod")
        _, subdir = mc_engine.run_simulator(
            self.exe, Lx, Ly, Tx, Ty, r1, r2, 1.0, beta_c,
            n_traj=self.n_traj_prod,
            data_dir=prod_dir,
        )
        if subdir is None:
            raise RuntimeError(f"no MC output subdir for {label}")
        a2a = os.path.join(subdir, "two_point_all_to_all.dat")
        test_data = mc_engine.load_all_to_all(a2a)

        # ---- 3. cost ----------------------------------------------------
        c_val, sig, per_dir, per_dir_sig = cost_module.l2_cost(
            self.ref_data, test_data,
            Lx, Ly, Tx, Ty,
            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
            mu=self.mu,
        )
        snr = cost_module.snr(c_val, sig)
        status = cost_module.snr_status(c_val, sig)

        wall = time.time() - t0
        result = EvalResult(
            eval_id=eid, r1=float(r1), r2=float(r2),
            beta_c=float(beta_c),
            cost=float(c_val), sigma_cost=float(sig),
            snr=float(snr), snr_status=status,
            per_dir=[float(x) for x in per_dir],
            per_dir_sigma=[float(x) for x in per_dir_sig],
            n_traj_prod=self.n_traj_prod,
            n_traj_scan_total=int(n_scan),
            wall_time_s=float(wall),
            scan_betas=[float(x) for x in sb],
            scan_chis=[float(x) for x in sc],
            scan_chi_errs=[float(x) for x in sce],
        )

        # ---- 4. log -----------------------------------------------------
        with open(self.log_path, "a") as f:
            f.write(json.dumps(asdict(result)) + "\n")

        print(f"[ev {eid:03d}]  cost={c_val:.4e}±{sig:.2e}  "
              f"SNR={snr:.2f} ({status})  wall={wall:.1f}s")

        # ---- 5. plotter -------------------------------------------------
        if self.optimizer_plot is not None:
            self.optimizer_plot.update(r1, r2, c_val, sig, beta_c, test_data,
                                       simplex=self.current_simplex)

        # ---- cleanup ----------------------------------------------------
        if not self.keep_mc_subdirs:
            shutil.rmtree(scratch, ignore_errors=True)

        self._beta_prev = beta_c
        return result

    # ------------------------------------------------------------------
    def cost_only(self, r1: float, r2: float) -> float:
        """Convenience: drop the EvalResult and return scalar cost."""
        return self(r1, r2).cost
