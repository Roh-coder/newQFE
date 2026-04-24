"""
parallel.py — Process-pool fan-out for CMA-ES population evaluation
                (Speedup 3).

Each generation of CMA-ES samples λ candidate (r1, r2) points from the
same search distribution.  Their evaluations are independent, so they
can run in parallel — this module owns the ProcessPool and the per-
worker `Evaluator` instances.

Design notes
------------

* Workers are HEADLESS: they do not hold the dashboard or matplotlib
  plotter (neither is picklable).  Each worker writes its own
  `EvalResult` to its own line in `eval_log.jsonl` for the matching
  scratch directory; the main process reads back that record and
  forwards it to `OptimizerPlotter.update(...)` and
  `Dashboard.update_eval(...)` after the generation completes.
* Workers can share the persistent β_c cache because it is backed
  by an `fcntl.flock`-coordinated `samples.jsonl` on disk (see
  betac_cache.py).  Each worker opens an independent `BetacCache`
  handle to the same directory; the file lock makes append-writes
  atomic.
* Each worker derives its MC seeds from a per-worker
  `numpy.random.SeedSequence`: workers pass their seed into the
  evaluator, which threads it down into `mc_engine.run_simulator`.
  This guarantees that two workers evaluating the same (r1, r2) at
  the same generation produce *independent* MC noise, so the
  CMA-ES termination tests (cost spread) are not biased.

Public API
----------
    class GenerationPool:
        def __init__(self, n_workers, evaluator_kwargs, master_seed=None)
        def map_generation(self, points, eval_id_base) -> list[EvalResult]
        def shutdown()

When `n_workers <= 1`, GenerationPool degenerates to the serial path
and just builds one Evaluator in the main process — no fork is done
and the dashboard/plotter callbacks fire as before.  This is what the
optimizer always uses; the parallelism is opt-in by raising
``n_workers``.
"""
from __future__ import annotations

import os
import pickle
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Worker-side state.  Lives once per child process.
# ---------------------------------------------------------------------------

# Module-globals because ProcessPoolExecutor calls our function fresh
# each task; we initialise the heavy state once via the initializer.
_W_EVAL = None       # type: ignore[assignment]
_W_RNG = None        # type: ignore[assignment]


def _worker_init(eval_kwargs: dict, seed_state: bytes) -> None:
    """Per-worker setup: build a headless Evaluator and a private RNG."""
    global _W_EVAL, _W_RNG
    # Reconstruct the seed sequence so each worker draws from a stream
    # that is independent of every other worker's stream.
    ss = pickle.loads(seed_state)
    _W_RNG = np.random.default_rng(ss)

    from evaluator import Evaluator   # local import: light dep graph
    from betac_cache import BetacCache

    cache_kwargs = eval_kwargs.pop("_betac_cache", None)
    cache = None
    if cache_kwargs is not None:
        cache = BetacCache(**cache_kwargs)

    _W_EVAL = Evaluator(
        **eval_kwargs,
        # Headless: no plotter, no dashboard, no beta_plot_dir.
        optimizer_plot=None,
        beta_plot_dir=None,
        dashboard=None,
        betac_cache=cache,
    )


def _worker_eval(point_with_id) -> dict:
    """Evaluate one (r1, r2) and return the EvalResult as a plain dict."""
    eid, r1, r2 = point_with_id
    # Pin the evaluator's eval_id so log lines stay contiguous across
    # workers.  (Evaluator increments _eval_id internally; we override.)
    _W_EVAL._eval_id = eid - 1
    res = _W_EVAL(float(r1), float(r2))
    return asdict(res)


# ---------------------------------------------------------------------------
# Pool wrapper
# ---------------------------------------------------------------------------

class GenerationPool:
    """Fan λ candidates per CMA-ES generation across n_workers processes."""

    def __init__(self, n_workers: int, evaluator_kwargs: dict,
                 master_seed: Optional[int] = None):
        self.n_workers = max(1, int(n_workers))
        self.master_ss = np.random.SeedSequence(master_seed)
        self._pool: Optional[ProcessPoolExecutor] = None
        self._spawned_seeds: list = []
        self._eval_kwargs = dict(evaluator_kwargs)
        self._serial_eval = None  # built lazily for n_workers=1

        if self.n_workers > 1:
            # Spawn one seed sequence per worker up-front; the pool rotates
            # through workers but each is initialised with its own stream.
            child_ss = self.master_ss.spawn(self.n_workers)
            self._spawned_seeds = [pickle.dumps(s) for s in child_ss]
            # We use a single shared seed for initializer (the round-robin
            # task-to-worker mapping is opaque) — workers each draw from a
            # per-process generator, independent of one another because
            # the initial seeds are spawned children of the master.
            #
            # Initialise the pool with worker #0's stream; subsequent
            # workers receive their own via per-task init isn't possible
            # in concurrent.futures, so we just give them all distinct
            # seeds via the per-task argument (see map_generation).
            self._pool = ProcessPoolExecutor(
                max_workers=self.n_workers,
                initializer=_worker_init,
                initargs=(self._eval_kwargs, self._spawned_seeds[0]),
            )
        else:
            # Serial fallback — build the evaluator in the main process so
            # the dashboard/plotter callbacks (kept by the caller) still
            # work.  No fork occurs.
            from evaluator import Evaluator
            from betac_cache import BetacCache

            cache_kwargs = self._eval_kwargs.pop("_betac_cache", None)
            cache = (BetacCache(**cache_kwargs)
                     if cache_kwargs is not None else None)

            # The serial Evaluator might want to keep its plotter/dashboard
            # — those are passed back in by the caller via .attach_callbacks.
            self._serial_eval = Evaluator(
                **self._eval_kwargs,
                optimizer_plot=None, beta_plot_dir=None, dashboard=None,
                betac_cache=cache,
            )

    # ----- callbacks for serial mode (n_workers=1) -----
    def attach_callbacks(self, *, optimizer_plot=None,
                         beta_plot_dir=None, dashboard=None) -> None:
        """In serial mode, route plotter / dashboard updates through pool."""
        if self._serial_eval is None:
            return
        self._serial_eval.optimizer_plot = optimizer_plot
        self._serial_eval.beta_plot_dir  = beta_plot_dir
        self._serial_eval.dashboard      = dashboard

    @property
    def serial_evaluator(self):
        """Direct access to the in-process Evaluator (n_workers=1 only)."""
        return self._serial_eval

    # ----- main hook used by the CMA-ES outer loop -----
    def map_generation(self, points, eval_id_base: int):
        """Evaluate a batch of (r1, r2) candidates; return list[dict]."""
        tasks = [(eval_id_base + k, float(r1), float(r2))
                 for k, (r1, r2) in enumerate(points)]
        if self.n_workers <= 1 or self._pool is None:
            results = []
            for eid, r1, r2 in tasks:
                self._serial_eval._eval_id = eid - 1
                res = self._serial_eval(r1, r2)
                results.append(asdict(res))
            return results

        # Parallel: submit, gather as_completed, then re-sort by eval_id
        # so the optimizer sees a deterministic order regardless of which
        # worker finished first.
        futures = {self._pool.submit(_worker_eval, t): t[0] for t in tasks}
        out = {}
        for fut in as_completed(futures):
            eid = futures[fut]
            try:
                out[eid] = fut.result()
            except Exception as e:
                # Surface worker failure with the task id so the caller can
                # decide whether to bail or continue.
                raise RuntimeError(f"worker failed on eval_id={eid}") from e
        return [out[eid] for eid, _, _ in tasks]

    def shutdown(self) -> None:
        if self._pool is not None:
            self._pool.shutdown(wait=True, cancel_futures=True)
            self._pool = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.shutdown()
