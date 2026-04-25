# Archive — K_from_Optimization_pnp (pre-speedup snapshot)

Captured **2026-04-24**, immediately before starting the
**S2 → S4 → S3 → S5** speedup implementation in the live
`K_from_Optimization_pnp/` directory.

This is a frozen snapshot of the working pipeline at the point where:

- **S2 (β_c cache)**: code landed (`betac_cache.py`, `tools/preseed_betac.py`,
  unit tests) but `CONFIG["betac_cache"] = False` (off by default).
- **S3 (parallel CMA-ES)**: code landed (`parallel.py`, optimizer patch,
  unit tests) but `CONFIG["n_workers"] = 1` (serial by default).
- **S4 (β reweighting)**: not started.
- **S5 (cost-surface surrogate / BO)**: not started.
- **S1 (C++ RPC backend)**: not started; deferred indefinitely per the
  revised plan.

All run artifacts (`results/j*/`, `results/nelder_mead/`, `__pycache__`,
`.pytest_cache`) were stripped before archiving.  The reference
correlator caches under `results/_reference_*` were kept because they
are expensive to regenerate and are immutable.

To recover this snapshot, copy this entire directory back to
`K_from_Optimization_pnp/`.
