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
    beta_c_sigma: float
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
                 scan_n_coarse=11, scan_n_refine=5,
                 scan_n_refine2=5, scan_n_refine3=5,
                 scan_max_shifts=4, scan_jackknife=False,
                 beta_seed=(0.20, 0.32),
                 optimizer_plot=None, beta_plot_dir=None,
                 keep_mc_subdirs=False, dashboard=None,
                 betac_cache=None,
                 reweight=False, reweight_n_traj=40000,
                 reweight_n_grid=201, reweight_n_eff_floor=0.10,
                 reweight_max_recenters=3,
                 cost_mode="l4mean_both_interp", cost_power=2):
        self.exe = exe
        self.ref_data = ref_data
        self.ref_geom = tuple(ref_geom)
        self.test_geom = tuple(test_geom)
        self.output_dir = output_dir
        self.n_traj_prod = int(n_traj_prod)
        self.n_traj_scan_coarse = int(n_traj_scan_coarse)
        self.n_traj_scan_fine = int(n_traj_scan_fine)
        self.scan_n_coarse = int(scan_n_coarse)
        self.scan_n_refine = int(scan_n_refine)
        self.scan_n_refine2 = int(scan_n_refine2)
        self.scan_n_refine3 = int(scan_n_refine3)
        self.scan_max_shifts = int(scan_max_shifts)
        self.scan_jackknife = bool(scan_jackknife)
        self.beta_seed = tuple(beta_seed)
        self.optimizer_plot = optimizer_plot
        self.beta_plot_dir = beta_plot_dir
        self.keep_mc_subdirs = keep_mc_subdirs
        self.dashboard = dashboard
        # Optional persistent β_c cache (Speedup 2).  When attached, we
        # try the cache before launching a fresh 3-pass scan; misses fall
        # through to mc_engine.find_beta_c and are then added back.
        self.betac_cache = betac_cache
        # Speedup 4 — Ferrenberg–Swendsen reweighting.  When True, the
        # 3-pass scan is replaced by a single long parent MC at the
        # bracket midpoint plus reweighted χ(β) over a dense grid.
        self.reweight = bool(reweight)
        self.reweight_n_traj = int(reweight_n_traj)
        self.reweight_n_grid = int(reweight_n_grid)
        self.reweight_n_eff_floor = float(reweight_n_eff_floor)
        self.reweight_max_recenters = int(reweight_max_recenters)
        self.cost_mode = str(cost_mode)
        self.cost_power = int(cost_power)

        os.makedirs(output_dir, exist_ok=True)
        self.log_path = os.path.join(output_dir, "eval_log.jsonl")
        if os.path.exists(self.log_path):
            os.remove(self.log_path)

        self._mc_root = os.path.join(output_dir, "mc_scratch")
        os.makedirs(self._mc_root, exist_ok=True)

        self._eval_id = 0
        self._beta_prev: Optional[float] = None
        # Set by the optimizer before each evaluator call so the visualizer
        # can draw the current simplex (Nelder-Mead) or Gaussian search
        # distribution (CMA-ES).
        self.current_simplex: Optional[list] = None
        # CMA-ES writes {"mean": [m1, m2], "cov": [[..],[..]], "sigma": float,
        # "gen": int} here at the start of each generation.
        self.current_gaussian: Optional[dict] = None
        # BO writes {"r1": 2‑D grid, "r2": 2‑D grid, "mean": 2‑D grid,
        # "std": 2‑D grid, "acq": 2‑D grid, "next": (r1, r2), "step": int}
        # here after each GP refit so the visualiser can draw the surface.
        self.current_gp_surface: Optional[dict] = None

    def __call__(self, r1: float, r2: float) -> EvalResult:
        self._eval_id += 1
        eid = self._eval_id
        Lx, Ly, Tx, Ty = self.test_geom
        ref_Lx, ref_Ly, ref_Tx, ref_Ty = self.ref_geom

        t0 = time.time()
        # Keep label short to stay under Windows MAX_PATH (260 chars).
        # The full r1/r2 values are recorded in eval_log.jsonl.
        label = f"ev{eid:04d}"
        scratch = os.path.join(self._mc_root, label)
        os.makedirs(scratch, exist_ok=True)

        # --- 1. β_c bracket: warm start or seed ---
        if self._beta_prev is not None:
            eps = max(0.20 * self._beta_prev, 0.05)
            beta_lo = max(0.01, self._beta_prev - eps)
            beta_hi = self._beta_prev + eps
        else:
            beta_lo, beta_hi = self.beta_seed

        beta_cbs = []
        if self.dashboard is not None:
            self.dashboard.begin_eval(eid, r1, r2, (beta_lo, beta_hi))
            beta_cbs.append(self.dashboard.update_scan)
        if self.optimizer_plot is not None and hasattr(self.optimizer_plot,
                                                       "update_scan"):
            beta_cbs.append(self.optimizer_plot.update_scan)
        if self.beta_plot_dir is not None:
            from visualization import BetaScanPlotter
            beta_cbs.append(BetaScanPlotter(self.beta_plot_dir, label,
                                            every_point=False))

        if not beta_cbs:
            beta_cb = None
        elif len(beta_cbs) == 1:
            beta_cb = beta_cbs[0]
        else:
            def beta_cb(state, _cbs=tuple(beta_cbs)):
                for cb in _cbs:
                    cb(state)

        print(f"[ev {eid:03d}] r1={r1:.4f} r2={r2:.4f}  "
              f"β∈[{beta_lo:.4f},{beta_hi:.4f}]")

        # ---- Speedup 2: try the persistent β_c cache first ----
        cache_hit = None
        if self.betac_cache is not None:
            cache_hit = self.betac_cache.lookup(r1, r2)

        if cache_hit is not None:
            beta_c, beta_c_sigma = cache_hit
            chi_peak = float("nan")
            sb, sc, sce = [], [], []
            n_scan = 0
            print(f"[ev {eid:03d}]  CACHE HIT β_c={beta_c:.6f}±{beta_c_sigma:.2e} "
                  f"(skipping 3-pass scan)")
        elif self.reweight:
            beta_c, beta_c_sigma, chi_peak, sb, sc, sce = mc_engine.find_beta_c_reweight(
                self.exe, Lx, Ly, Tx, Ty, r1, r2, 1.0,
                beta_lo, beta_hi,
                n_traj_parent=self.reweight_n_traj,
                n_grid=self.reweight_n_grid,
                n_eff_floor=self.reweight_n_eff_floor,
                max_recenters=self.reweight_max_recenters,
                data_dir=os.path.join(scratch, "scan_rw"),
                jackknife=self.scan_jackknife,
                progress_cb=beta_cb,
            )
            # Per-eval trajectory budget for accounting:
            # (max_recenters + 1) parent runs of n_traj_parent each.
            n_scan = self.reweight_n_traj * (self.reweight_max_recenters + 1)
            print(f"[ev {eid:03d}]  β_c={beta_c:.6f}"
                  f"{(' ±' + format(beta_c_sigma, '.2e')) if beta_c_sigma > 0 else ''}"
                  f"  χ_peak={chi_peak:.3g}  (reweight, {self.reweight_n_grid} grid pts)")
            if self.betac_cache is not None:
                self.betac_cache.add(
                    r1, r2, beta_c, beta_c_sigma,
                    source_run=os.path.basename(self.output_dir),
                    n_traj_total=n_scan,
                )
        else:
            beta_c, beta_c_sigma, chi_peak, sb, sc, sce = mc_engine.find_beta_c(
                self.exe, Lx, Ly, Tx, Ty, r1, r2, 1.0,
                beta_lo, beta_hi,
                n_coarse=self.scan_n_coarse,
                n_refine=self.scan_n_refine,
                n_refine2=self.scan_n_refine2,
                n_refine3=self.scan_n_refine3,
                n_traj_coarse=self.n_traj_scan_coarse,
                n_traj_fine=self.n_traj_scan_fine,
                max_shifts=self.scan_max_shifts,
                jackknife=self.scan_jackknife,
                data_dir=os.path.join(scratch, "scan"),
                progress_cb=beta_cb,
            )
            n_scan = (len(sb) *
                      (self.n_traj_scan_coarse + self.n_traj_scan_fine) // 2)
            if beta_c_sigma > 0:
                print(f"[ev {eid:03d}]  β_c={beta_c:.6f}±{beta_c_sigma:.2e}  "
                      f"χ_peak={chi_peak:.3g}  ({len(sb)} scan pts)")
            else:
                print(f"[ev {eid:03d}]  β_c={beta_c:.6f}  χ_peak={chi_peak:.3g}  "
                      f"({len(sb)} scan pts)")
            # Add the freshly-computed β_c to the cache for future hits.
            if self.betac_cache is not None:
                self.betac_cache.add(
                    r1, r2, beta_c, beta_c_sigma,
                    source_run=os.path.basename(self.output_dir),
                    n_traj_total=n_scan,
                )

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
            cost_mode=self.cost_mode,
            cost_power=self.cost_power,
        )
        snr_val = cost_module.snr(c_val, sig)
        status = cost_module.snr_status(c_val, sig)

        wall = time.time() - t0
        result = EvalResult(
            eval_id=eid, r1=float(r1), r2=float(r2),
            beta_c=float(beta_c),
            beta_c_sigma=float(beta_c_sigma),
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
                                       simplex=self.current_simplex,
                                       gaussian=self.current_gaussian,
                                       gp_surface=self.current_gp_surface)
        if self.dashboard is not None:
            self.dashboard.update_eval(result)

        if not self.keep_mc_subdirs:
            shutil.rmtree(scratch, ignore_errors=True)

        self._beta_prev = beta_c
        return result
