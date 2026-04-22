"""
evaluator.py — single function from (r1, r2) to (cost, sigma_cost).

Each call:
  1. Find beta_c on the test lattice with k1=r1, k2=r2, k3=1
     (warm-started from the previous call's beta_c if available).
  2. Run a longer production MC at that beta_c.
  3. Compute the L² cost vs the cached reference correlator.
  4. Append a record to eval_log.jsonl.
  5. Push a frame to OptimizerPlotter (if attached).
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
    scan_betas: list = None
    scan_chis: list = None
    scan_chi_errs: list = None


class Evaluator:
    """Callable evaluator: ``result = ev(r1, r2)``."""

    def __init__(self, exe, ref_data, ref_geom, test_geom, output_dir,
                 n_traj_prod=30000, n_traj_scan_coarse=4000,
                 n_traj_scan_fine=10000,
                 beta_seed=(0.20, 0.32),
                 optimizer_plot=None, beta_plot_dir=None,
                 keep_mc_subdirs=False, dashboard=None):
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
        self.dashboard = dashboard

        os.makedirs(output_dir, exist_ok=True)
        self.log_path = os.path.join(output_dir, "eval_log.jsonl")
        if os.path.exists(self.log_path):
            os.remove(self.log_path)

        self._mc_root = os.path.join(output_dir, "mc_scratch")
        os.makedirs(self._mc_root, exist_ok=True)

        self._eval_id = 0
        self._beta_prev: Optional[float] = None
        # Set by the optimizer before each evaluator call so the visualizer
        # can draw the current simplex.
        self.current_simplex: Optional[list] = None

    def __call__(self, r1: float, r2: float) -> EvalResult:
        self._eval_id += 1
        eid = self._eval_id
        Lx, Ly, Tx, Ty = self.test_geom
        ref_Lx, ref_Ly, ref_Tx, ref_Ty = self.ref_geom

        t0 = time.time()
        label = f"eval{eid:04d}_r1_{r1:.4f}_r2_{r2:.4f}"
        scratch = os.path.join(self._mc_root, label)
        os.makedirs(scratch, exist_ok=True)

        # --- 1. β_c bracket: warm start or seed ---
        if self._beta_prev is not None:
            eps = max(0.20 * self._beta_prev, 0.05)
            beta_lo = max(0.01, self._beta_prev - eps)
            beta_hi = self._beta_prev + eps
        else:
            beta_lo, beta_hi = self.beta_seed

        beta_cb = None
        if self.dashboard is not None:
            self.dashboard.begin_eval(eid, r1, r2, (beta_lo, beta_hi))
            beta_cb = self.dashboard.update_scan
        elif self.beta_plot_dir is not None:
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

        # --- 2. production MC at β_c ---
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

        # --- 3. cost ---
        c_val, sig, per_dir, per_dir_sig = cost_module.l2_cost(
            self.ref_data, test_data,
            Lx, Ly, Tx, Ty,
            ref_Lx, ref_Ly, ref_Tx, ref_Ty,
        )
        snr_val = cost_module.snr(c_val, sig)
        status = cost_module.snr_status(c_val, sig)

        wall = time.time() - t0
        result = EvalResult(
            eval_id=eid, r1=float(r1), r2=float(r2),
            beta_c=float(beta_c),
            cost=float(c_val), sigma_cost=float(sig),
            snr=float(snr_val), snr_status=status,
            per_dir=[float(x) for x in per_dir],
            per_dir_sigma=[float(x) for x in per_dir_sig],
            n_traj_prod=self.n_traj_prod,
            n_traj_scan_total=int(n_scan),
            wall_time_s=float(wall),
            scan_betas=[float(x) for x in sb],
            scan_chis=[float(x) for x in sc],
            scan_chi_errs=[float(x) for x in sce],
        )

        with open(self.log_path, "a") as f:
            f.write(json.dumps(asdict(result)) + "\n")

        print(f"[ev {eid:03d}]  cost={c_val:.4e}±{sig:.2e}  "
              f"SNR={snr_val:.2f} ({status})  wall={wall:.1f}s")

        if self.optimizer_plot is not None:
            self.optimizer_plot.update(r1, r2, c_val, sig, beta_c, test_data,
                                       simplex=self.current_simplex)
        if self.dashboard is not None:
            self.dashboard.update_eval(result)

        if not self.keep_mc_subdirs:
            shutil.rmtree(scratch, ignore_errors=True)

        self._beta_prev = beta_c
        return result
