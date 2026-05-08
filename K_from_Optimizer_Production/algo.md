# algo.md — Twisted ↔ Untwisted 4‑5‑6 Matching in the Continuum Limit

> **Audience.** A PhD‑level lattice / CFT practitioner picking up this
> project cold. This document explains *what* problem we are solving,
> *why* the obvious approaches did not work, and *how* the validated
> matching algorithm runs end‑to‑end. Section anchors map directly to
> source files in [`K_from_Optimizer_Production/`](.).
>
> Companion documents:
> * [`README.md`](README.md) — operational guide for the optimizer machinery.
> * [`plan.md`](plan.md) — historical cost‑function design notes.
> * [`LANDSCAPE.md`](LANDSCAPE.md) — cost‑landscape diagnostics.
>
> §10 (“Build from scratch”) gives the exact commands to reproduce
> every result in this document starting from a fresh clone of the
> [`brower/newQFE`](https://github.com/brower/newQFE) repository.

---

## 1. Project context

### 1.1 Physics goal

We are tuning a **2D triangular‑lattice Ising model on a parallelogram
torus** to realise a chosen continuum CFT geometry — a flat torus whose
modular parameter $\tau$ is fixed by the three boundary‑cycle lengths
of a target Euclidean triangle (the canonical case is the **4‑5‑6
triangle**). The triangular Ising critical surface is

$$
\sinh(2\beta K_1)\sinh(2\beta K_2)+\sinh(2\beta K_2)\sinh(2\beta K_3)
+\sinh(2\beta K_3)\sinh(2\beta K_1)=1.
$$

There are two equivalent ways to embed a 4‑5‑6 triangle into a
triangular lattice:

| Approach | Geometry $(L_x,L_y,T_x,T_y)$ | Couplings $(K_1,K_2,K_3)$ | Critical $\beta_c$ |
|---|---|---|---|
| **Twisted, isotropic** (REF) | $(13\alpha, 16\alpha, -3\alpha, 3\alpha)$ | $(1,1,1)$ | from MC peak of $\chi(\beta)$ |
| **Untwisted, anisotropic** (TEST) | $(L,L,0,0)$ | $(r_1, r_2, 1)$ exact | $\ln(3/\sqrt 7)/2 \approx 0.06283$ |

with the closed‑form anisotropic optima

$$
r_1=\frac{2\ln 5-\ln 7}{2\ln 3-\ln 7}\approx 5.06523,
\qquad
r_2=\frac{\ln 7}{2\ln 3-\ln 7}\approx 7.74293.
$$

Both lattices should flow to the **same continuum CFT** (free Majorana
at $c=1/2$). The matching question is: *can we verify this from MC
two‑point functions alone, without invoking 2D‑Ising specifics
(modular forms, $\Delta_\sigma=1/8$, theta functions)?* Answering yes
is the prerequisite for porting the workflow to $\phi^4$ and other
scalar QFTs.

### 1.2 Why the question matters operationally

The whole optimizer stack in this directory ([`run.py`](run.py),
[`cost.py`](cost.py), [`mc_engine.py`](mc_engine.py), …) tries to find
the anisotropic couplings $(r_1,r_2)$ that make a small **untwisted**
test torus match a high‑statistics **twisted** reference torus. If the
twisted and untwisted lattices are *not* CFT‑equivalent (or if the
observable used to compare them is mis‑specified), the optimizer is
chasing a phantom. The matching algorithm in §4 is the empirical
ground‑truth check that the two simulations realise the same continuum
limit, with the residual difference being a single overall operator
normalisation $Z_\sigma$.

---

## 2. The 4‑5‑6 geometry, in detail

### 2.1 Bond / angle dictionary

The triangular‑lattice unit cell uses bond vectors
$\vec e_1=(1,0)$, $\vec e_2=(0,1)$, $\vec e_3=\vec e_2-\vec e_1=(-1,1)$
in lattice indices $(m,n)$, mapped to triangular Cartesian coordinates
via $(x,y)=(m+\tfrac12 n,\,\tfrac{\sqrt 3}{2}n)$ (see
[`cost.py:_triangular_xy`](cost.py)).

Solving the sinh‑rule with the half‑angle identity gives the
correspondence

| Coupling | Bond direction | Triangle angle | Triangle side opposite |
|---|---|---|---|
| $K_1$ | $\vec e_1$ | $B$ | side **5** |
| $K_2$ | $\vec e_2$ | $A$ | side **4** |
| $K_3$ | $\vec e_3$ | $C$ | side **6** |

For the 4‑5‑6 triangle, $\sin A:\sin B:\sin C = 4:5:6$, so the test
lattice has *physical* bond lengths in ratio $5:4:6$ for the three
bond directions $\vec e_1,\vec e_2,\vec e_3$.

### 2.2 Three boundary cycles

The torus has three fundamental boundary direction vectors
([`cost.py:boundary_paths`](cost.py)):

```python
def boundary_paths(Lx, Ly, Tx, Ty):
    return [(Lx, Ty), (Tx, -Ly), (-Lx - Tx, Ly - Ty)]
```

Pulling the actual cycle lengths through the triangular embedding:

| Cycle index $c$ | REF $(13,16,-3,3)$ vector | REF length | REF physical side | TEST $(L,L,0,0)$ vector | TEST physical side |
|---:|---|---:|---:|---|---:|
| 0 | $(13, 3)$  | 14.731 | **5** | $(L, 0)\propto \vec e_1$        | **5** |
| 1 | $(-3,-16)$ | 17.692 | **6** | $(0,-L)\propto -\vec e_2$       | **4** |
| 2 | $(-10,13)$ | 11.790 | **4** | $(-L, L)\propto \vec e_3$       | **6** |

Both lattices realise the same triangle, but the **cycle index ↔
physical side correspondence is permuted** between them, because the
reference encodes the triangle through twists while the test encodes
it through the bond‑coupling map of §2.1. This permutation is the
single technical insight that the rest of this document depends on.

```text
                cycle 0 cycle 1 cycle 2
   REF   →  side  5       6       4
   TEST  →  side  5       4       6
```

### 2.3 Reference $\beta_c$

For REF $(13,16,-3,3)$ with $K_1=K_2=K_3=1$, the analytic continuum
critical inverse temperature is $\beta_c=\ln(3/\sqrt 7)/2\approx
0.06283$, *far* below the default $[0.20,0.32]$ scan bracket used for
near‑equilateral runs. The bracket‑translation logic in
[`mc_engine.find_beta_c_reweight`](mc_engine.py) /
[`find_beta_c_multidonor_2pass`](mc_engine.py) is what makes the β
finder robust here (see [`README.md` §4](README.md)).

---

## 3. What did **not** work, and why

The history is documented at length in [`plan.md`](plan.md) and
[`LANDSCAPE.md`](LANDSCAPE.md). Headlines:

1. **All 14 standard cost kernels** (raw / log / Huber L_p, χ²,
   relative, drop‑short‑log, effmass, slope‑loglog, affine, cosine)
   place their global minimum at the small‑$r$ corner $(0.5,0.5)$
   when run as `cost(G_\text{test}, G_\text{ref})` on the dense
   precomputed grid (see [`landscape_dozen.py`](landscape_dozen.py),
   [`results/_landscape/.../dozen.png`](results/_landscape/Lx16_Ly16_Tx0_Ty0/dozen.png)).
   The tilt is **structural**: at small $r$ all cycle lengths in
   physical units shrink, so $G\to 1$ and any residual built on
   $\log G$ is trivially small.

2. **Continuum‑extrapolated FSS by cycle index** (the original
   [`plot_fss_continuum_limit_456.py`](plot_fss_continuum_limit_456.py))
   gave per‑cycle test/ref disagreements of $-5\%$, $+15\%$, $-20\%$
   that *did not improve* with better fitting (we tried fixed‑$\omega$,
   free‑$\omega$, polynomial $a/L+b/L^2$, Padé, exponential
   $A\exp(b/L)$). The disagreement was nearly **antisymmetric across
   cycles 1 and 2**, with the asymptote ratio approximately reciprocal
   $(A_1^{\text{test}}/A_1^{\text{ref}}, A_2^{\text{test}}/A_2^{\text{ref}})
   \approx (1.12,\,0.93)$ vs $(0.92,\,1.11)$. That is the fingerprint of
   a **cycle‑label permutation**, not a fitting failure.

3. **Rescaling $L_\text{phys}$ on the test family** by the per‑cycle
   reference base length had no effect on the asymptote because $A$ in
   $G=A\exp(b/L)$ is invariant under $L\to\lambda L$.

The conclusion that drove the fix: the comparison was being done on the
*wrong observable pairing*. Once you pair by physical side rather than
cycle index, and you strip the operator normalisation $Z_\sigma$ before
comparing, the disagreement collapses to $\lesssim 0.3\%$.

---

## 4. The matching algorithm

Reference implementation:
[`match_456_twisted_vs_untwisted.py`](match_456_twisted_vs_untwisted.py).
Output figure:
[`results/_fss_456/match_twisted_vs_untwisted.png`](results/_fss_456/match_twisted_vs_untwisted.png).

### 4.1 Inputs

* REF dataset family: $(L_x,L_y,T_x,T_y)=(13\alpha,16\alpha,-3\alpha,3\alpha)$,
  $K_1=K_2=K_3=1$, $\alpha\in\{1,2,3,4\}$. Generated by
  [`precompute_456_fss.py`](precompute_456_fss.py) and
  [`precompute_456_fss_extend.py`](precompute_456_fss_extend.py).
* TEST dataset family: $(L,L,0,0)$ at exact analytic $(K_1,K_2,K_3)=(r_1,r_2,1)$,
  $L\in\{8,16,24,32,48,64\}$. Same precompute scripts.
* Each dataset is the all‑to‑all connected two‑point function
  $G_{\text{conn}}(m,n)$ stored in `two_point_all_to_all.dat` and loaded by
  [`mc_engine.load_all_to_all`](mc_engine.py).

### 4.2 Algorithm (six steps)

**Step 1 — build a tiled interpolant per dataset.**
Use [`cost._tile_interp`](cost.py) to wrap each lattice's MC data into
a `LinearNDInterpolator` over $\geq 5\times 5$ torus copies in
triangular Cartesian coordinates. Cache one interpolant per `(dataset,
geometry)` pair (the fast version in
[`match_456_twisted_vs_untwisted.py`](match_456_twisted_vs_untwisted.py)
keys on `id(dat)` plus geometry; this gives roughly $200\times$ speedup
versus rebuilding per sample point).

**Step 2 — fix the cycle ↔ side mapping.**
Hard‑code the permutation derived in §2.2:

```python
REF_CYCLE_SIDE  = {0: 5, 1: 6, 2: 4}   # twisted, isotropic k
TEST_CYCLE_SIDE = {0: 5, 1: 4, 2: 6}   # untwisted, anisotropic k
SIDE_PAIR = {s: (ref_cycle_for_side[s], test_cycle_for_side[s]) for s in (4,5,6)}
```

This is the single line that converts the failing analysis into the
working one.

**Step 3 — sample each cycle at common fractional positions.**
Take $t_k = k/8$ for $k=1,\dots,7$ along each cycle. The cycle vector
in lattice indices is $(dm,dn)=$ `boundary_paths(...)[c]`; the
Cartesian position is $(x,y)=(dm+\tfrac12 dn,\tfrac{\sqrt 3}{2}dn)\cdot t$.
Evaluate the cached interpolants at these points to get
$G_s^{\text{test}}(t_k)$ and $G_s^{\text{ref}}(t_k)$ with their MC
errors, where $s\in\{4,5,6\}$ is the *physical side*.

**Step 4 — strip the operator normalisation.**
For each side $s$ and each largest available size (TEST $L=64$,
REF $\alpha=4$), form two complementary representations:

$$
R_s(t) \;=\; \frac{G_s(t)}{G_s(1/2)}\qquad
\text{(}Z_\sigma\text{-invariant ratio)}
$$

$$
\hat A_s \;=\; \arg\min_A \sum_k
\frac{\bigl(G_s(t_k) - A\sin(\pi t_k)^{-1/4}\bigr)^2}{\sigma_s(t_k)^2}
\qquad\text{(CFT‑flat amplitude WLS)}
$$

The exponent $1/4=2\Delta_\sigma$ is the only piece of 2D‑Ising input
in the comparison; everything else (the cycle pairing, the ratio
construction) is theory‑agnostic.

**Step 5 — overlay $R_s(t)$ across families.**
Plot $R_s^{\text{test}}(t)$ and $R_s^{\text{ref}}(t)$ on the same axes
for each $s\in\{4,5,6\}$. Continuum‑CFT equivalence demands the two
curves coincide within MC errors. The bottom row of the output figure
is exactly this overlay, with the analytic
$\sin(\pi t)^{-1/4}/\sin(\pi/2)^{-1/4}$ as a guide.

**Step 6 — check $Z_\sigma$ universality across sides.**
Compute the three side‑wise amplitude ratios

$$
\rho_s \equiv \hat A_s^{\text{test}} / \hat A_s^{\text{ref}},
\qquad s\in\{4,5,6\}.
$$

If the two lattices realise the **same** CFT then $\rho_s$ must be
*independent of $s$* — it can only be the ratio of the two operator
normalisations $Z_\sigma^{\text{test}}/Z_\sigma^{\text{ref}}$. The
quantitative test is

$$
\mathrm{spread}(\rho)/\mathrm{mean}(\rho) \;\ll\; 1\%.
$$

### 4.3 Result on current data

From [`match_456_twisted_vs_untwisted.py`](match_456_twisted_vs_untwisted.py)
at TEST $L=64$, REF $\alpha=4$:

| Side $s$ | $\hat A_s^{\text{ref}}$ | $\hat A_s^{\text{test}}$ | $\rho_s$ |
|---:|---:|---:|---:|
| 4 | 0.2429 | 0.2382 | 0.9807 |
| 5 | 0.2170 | 0.2128 | 0.9805 |
| 6 | 0.2003 | 0.1952 | 0.9747 |

$$
\bar\rho = 0.9787,\quad
\sigma(\rho) = 0.0028,\quad
\sigma/\bar\rho = 0.28\%.
$$

The bottom‑row $R_s(t)$ overlay shows test (filled circles) and ref
(open squares) coincident within errors, both tracking the analytic
$\sin(\pi t)^{-1/4}$ curve. The top row $G_s(t)\sin(\pi t)^{1/4}$ is
nearly flat, with the small residual U‑shape signalling the same
finite‑$L$ correction on both lattices.

**Interpretation.** The two lattices realise the same continuum CFT to
sub‑percent precision. The residual $\sim 2\%$ overall offset
$1-\bar\rho$ is the operator normalisation
$Z_\sigma^{\text{test}}/Z_\sigma^{\text{ref}}$ — a non‑universal,
lattice‑geometry‑specific constant that drops out of any CFT
prediction and out of any Z‑invariant observable. The spread across
sides $0.28\%$ is consistent with MC statistics at the available
trajectory counts.

---

## 5. Implementation walkthrough

Below is the minimal code path. All references resolve in this
directory; no external packages beyond `numpy`, `scipy`, `matplotlib`.

```python
import os, sys, math
import numpy as np
from scipy.optimize import curve_fit

import mc_engine
from cost import boundary_paths, _SQRT3_2, _tile_interp

# 1. Load per‑family datasets
ref_data  = {a: mc_engine.load_all_to_all(f"results/_fss_456/ref/a{a}/two_point_all_to_all.dat")
             for a in [1,2,3,4]}
test_data = {L: mc_engine.load_all_to_all(f"results/_fss_456/test/L{L}/two_point_all_to_all.dat")
             for L in [8,16,24,32,48,64]}

# 2. Cache one tiled interpolant per (dataset, geometry)
_cache = {}
def interps(dat, Lx, Ly, Tx, Ty, copies=2):
    key = (id(dat), Lx, Ly, Tx, Ty)
    if key not in _cache:
        iG  = _tile_interp(dat, Lx, Ly, Tx, Ty, "conn",     copies)
        iGe = _tile_interp(dat, Lx, Ly, Tx, Ty, "conn_err", copies)
        _cache[key] = (iG, iGe, boundary_paths(Lx, Ly, Tx, Ty))
    return _cache[key]

def sample(dat, geom, c, t):
    iG, iGe, paths = interps(dat, *geom)
    dm, dn = paths[c]
    ex, ey = dm + 0.5*dn, _SQRT3_2*dn
    xy = np.column_stack([t*ex, t*ey])
    return np.asarray(iG(xy)), np.asarray(iGe(xy))

# 3. Side ↔ cycle mapping
REF_SIDE_OF  = {0: 5, 1: 6, 2: 4}
TEST_SIDE_OF = {0: 5, 1: 4, 2: 6}
PAIR = {s: (next(c for c,v in REF_SIDE_OF.items()  if v==s),
            next(c for c,v in TEST_SIDE_OF.items() if v==s))
        for s in (4,5,6)}

# 4. Compare on largest sizes
T = np.array([k/8 for k in range(1,8)])
ref_geom  = lambda a: (13*a, 16*a, -3*a, 3*a)
test_geom = lambda L: (L, L, 0, 0)

amp = lambda t, G, Ge: float(
    np.sum(np.sin(np.pi*t)**(-0.25) * G / Ge**2)
  / np.sum(np.sin(np.pi*t)**(-0.5)        / Ge**2))

for s, (rc, tc) in PAIR.items():
    Gr, Ger = sample(ref_data[4],  ref_geom(4),  rc, T)
    Gt, Get = sample(test_data[64], test_geom(64), tc, T)
    A_ref, A_test = amp(T, Gr, Ger), amp(T, Gt, Get)
    print(f"side {s}: A_test/A_ref = {A_test/A_ref:.4f}")
```

The physical content of the algorithm is fully captured by the
`PAIR` dictionary plus the WLS `amp` formula; everything else is
plumbing.

---

## 6. Generalisation rules (porting to other geometries / models)

The matching procedure above is *almost* model‑agnostic. To port:

1. **Re‑derive the bond–side dictionary** for the new triangle target.
   The map $K_i \leftrightarrow$ angle is fixed by the sinh‑rule
   half‑angle identity in §2.1; for any target triangle there is a
   unique permutation of cycle indices to physical sides on each
   embedding. The general procedure is in
   [`cost.boundary_paths`](cost.py); the side label per cycle is the
   length $|p_c|$ of that boundary vector in physical units.

2. **Replace the CFT‑flat exponent.** $1/4 = 2\Delta_\sigma$ is
   2D‑Ising specific. For $\phi^4$ the analogous exponent is
   $2\Delta_\phi$; for any CFT, fit the slope of
   $\log G$ vs $\log\sin(\pi t)$ near $t=1/2$ on the *reference*
   only (it is universal), then use it for both families. Any choice
   that is **the same** for ref and test leaves $\rho_s$ unchanged —
   so for matching purposes the exponent is a knob that does not
   affect the universality test, only the physical interpretation of
   $\hat A$.

3. **Z‑invariant ratio is exponent‑free.** $R_s(t) = G_s(t)/G_s(1/2)$
   is the cleanest universality probe and requires no analytic input
   at all. Treat it as primary; treat $\hat A_s$ as a secondary check
   that pins down the $Z_\sigma$ ratio.

4. **Always sample on $t\in[1/8,7/8]$.** Lattice contact terms
   dominate near $t=0,1$; the FSS analysis in
   [`plot_fss_continuum_limit_456.py`](plot_fss_continuum_limit_456.py)
   shows the contact‑term region pollutes the residual without adding
   information about the bulk CFT.

5. **Use the largest sizes available.** Smaller sizes carry the
   same continuum information modulated by FSS corrections that
   *differ* between the two families (twists vs anisotropy give
   different leading corrections). The Z‑invariant overlay
   $R_s^{\text{test}}\stackrel{?}{=}R_s^{\text{ref}}$ holds
   asymptotically; at finite $L$ both curves carry $1/L^p$ corrections
   that bias the comparison if combined naively.

---

## 7. Where this fits in the optimizer pipeline

* **Cost design.** Any new candidate cost in [`cost.py`](cost.py)
  must, *at the truth couplings*, give an output fingerprint
  consistent with the matching algorithm above. The corollary is the
  workflow rule from [`LANDSCAPE.md`](LANDSCAPE.md): a candidate cost
  must show a basin at truth on the dense MC landscape *and* its
  truth‑evaluated correlator must satisfy $\rho_s$ universality.
* **Validation harness.** `match_456_twisted_vs_untwisted.py` is the
  end‑to‑end check that the matching question has an affirmative
  answer for the current geometry. It should be re‑run whenever:
  * the precompute MC data is regenerated;
  * the bond–side dictionary in §2.1 is changed (e.g. new target
    triangle);
  * the simulator
    [`bin/ising_tri_twisted_parallelogram`](bin/ising_tri_twisted_parallelogram)
    is rebuilt with a different bond convention.
* **CMA‑ES interpretation.** The optimizer's objective is to find
  $(r_1,r_2)$ that *recovers* the matching shown here — not to
  verify the matching itself. Treat the matching algorithm as
  **ground truth**; treat optimizer results as estimates of the
  truth couplings whose correctness is judged by re‑running the
  matching at the optimizer's recovered point.

---

## 8. Lessons learned (re‑usable across the project)

1. **Cycle index is not a physical label.** Whenever two lattices
   encode the same target geometry through different mechanisms
   (twist vs anisotropy), establish the cycle ↔ physical observable
   map *before* writing any cost. The 4‑5‑6 cycle permutation
   (REF $5,6,4$ vs TEST $5,4,6$) is the canonical example; analogous
   permutations will appear in any new geometry.
2. **Z_σ is geometry‑specific even at criticality.** Two lattices can
   realise the same continuum CFT and still differ by a
   non‑universal $\sim O(\text{few \%})$ multiplicative constant on
   $\langle\sigma\sigma\rangle$. Any matching observable must be
   either explicitly $Z_\sigma$‑invariant ($R_s(t)$) or treat the
   global $Z$ ratio as a single fitted parameter that is **shared
   across the three sides**.
3. **Better fitting cannot fix a wrong observable.** The progression
   in [`plot_fss_continuum_limit_456.py`](plot_fss_continuum_limit_456.py)
   exhausted six different FSS ansätze before the cycle‑label bug was
   identified; none of them helped. When residuals show a systematic
   sign pattern that is invariant under the fitting choice, the
   problem is in the construction of the residual, not the model.
4. **Z‑invariant first, amplitudes second.** The bottom row of the
   output figure (the $R_s(t)$ overlay) carries the universality
   verdict; the top row ($\hat A_s$) carries the $Z_\sigma$
   diagnostic. Showing both keeps the qualitative and quantitative
   results visually separated.

---

## 9. Future work

* **FSS extrapolation of $\rho_s$.** Repeat §4.2 step 6 at each
  available $(L,\alpha)$ pair, then extrapolate $\rho_s(L)\to\rho_s^\infty$
  with a free‑exponent FSS fit. A successful run gives a single
  continuum number $\rho^\infty$ valid for all three sides.
* **Replace $\sin(\pi t)^{-1/4}$ with a fitted slope.** Fit the
  short‑distance exponent on the ref alone, then use it for both
  families. This drops the only piece of analytic 2D‑Ising input.
* **Apply the same algorithm at non‑critical β.** At a fixed
  off‑critical β the matching test still works; the curves will not
  be CFT but will collapse onto each other if the lattice geometries
  agree.
* **Port to $\phi^4$ on the same triangular torus.** Same code path,
  same cycle ↔ side dictionary, different MC simulator. The Z‑invariant
  ratio test is the cleanest universality probe and requires no
  analytic input from the φ⁴ side.
* **Wire the matching figure into the optimizer dashboard.** When CMA‑ES
  reports a converged $(r_1,r_2)$, automatically re‑run
  [`match_456_twisted_vs_untwisted.py`](match_456_twisted_vs_untwisted.py)
  at that point and embed the resulting $\rho_s$ table in
  `results/<run>/summary.json`. A converged run with
  $\sigma(\rho)/\bar\rho > 1\%$ should be flagged as untrusted
  regardless of the cost value.

---

## 10. Build the workflow from scratch (newQFE repository)

This section reproduces every figure in this document starting from a
fresh clone. All paths are relative to the repository root unless
otherwise noted.

### 10.1 Repository layout (relevant subset)

```
newQFE/
├── include/                                  ← shared C++ headers
│   ├── ising.h, lattice.h, statistics.h, ...
└── K_from_Optimizer_Production/              ← this directory
    ├── Makefile                              ← builds the simulator
    ├── requirements.txt                      ← Python deps
    ├── src/ising_tri_twisted_parallelogram.cc← MC simulator source
    ├── include/                              ← simulator‑local headers
    ├── bin/                                  ← built binary lives here
    ├── mc_engine.py                          ← Python wrapper (β_c finder + runner)
    ├── cost.py                               ← geometry helpers + cost kernels
    ├── precompute_456_fss.py                 ← MC for ref α∈{1,2,3}, test L∈{8..32}
    ├── precompute_456_fss_extend.py          ← MC for ref α=4, test L∈{48,64}
    ├── match_456_twisted_vs_untwisted.py     ← matching algorithm (§4)
    └── results/_fss_456/                     ← MC outputs land here
```

The simulator binary depends only on the headers in
[`K_from_Optimizer_Production/include/`](include/) (a frozen,
self‑contained snapshot of the relevant pieces of the top‑level
[`include/`](../include/) tree). No top‑level build is required.

### 10.2 System prerequisites

| Component | Tested version | Notes |
|---|---|---|
| Linux / macOS / WSL | — | dev‑container `mcr.microsoft.com/devcontainers/base:ubuntu-24.04` works |
| `g++` | ≥ 9 (C++14) | only standard library; no MPI / OpenMP / CUDA |
| `make` | any | one‑target Makefile |
| Python | ≥ 3.10 | stdlib + the four packages below |
| `numpy`, `scipy`, `matplotlib`, `rich`, `pytest` | per [`requirements.txt`](requirements.txt) | install in a virtualenv |

No GPU, no MPI, no proprietary compilers. A cold build + the full
matching pipeline on a 4‑core laptop completes in roughly half a day,
dominated by REF $\alpha=4$ and TEST $L=64$ MC.

### 10.3 Bootstrap commands (cold start)

```bash
# A. Clone
git clone https://github.com/brower/newQFE.git
cd newQFE/K_from_Optimizer_Production

# B. Python environment (one‑time)
python -m venv ../.venv
source ../.venv/bin/activate
pip install -r requirements.txt

# C. Build the C++ simulator
make                          # produces bin/ising_tri_twisted_parallelogram
# Sanity check
./bin/ising_tri_twisted_parallelogram --help | head
```

The simulator has no external dependencies; if `make` fails the only
likely cause is a too‑old `g++` (need C++14 support).

### 10.4 Generate the MC datasets

The two `precompute_*.py` scripts are **resumable** — if a target
`two_point_all_to_all.dat` already exists they skip it. Both internally
call [`mc_engine.find_beta_c`](mc_engine.py) to locate the
finite‑volume pseudo‑critical $\beta$, then
[`mc_engine.run_simulator`](mc_engine.py) to produce the all‑to‑all
two‑point function.

```bash
# Phase 1: REF α∈{1,2,3} and TEST L∈{8,16,24,32}
#   - α=1 and α=3 are symlinked from any pre‑existing
#     results/_reference_Lx{13,39}_Ly{16,48}_Tx{-3,-9}_Ty{3,9}/ caches
#     if present; otherwise simulated fresh.
#   - 4 workers * ~10–30 min per cell.
python precompute_456_fss.py --n-workers 4

# Phase 2: REF α=4 and TEST L∈{48,64}  (the asymptotic sizes used in §4.3)
python precompute_456_fss_extend.py --n-workers 4
```

When complete:

```
results/_fss_456/
├── ref/{a1,a2,a3,a4}/two_point_all_to_all.dat   (+ meta.json)
└── test/{L8,L16,L24,L32,L48,L64}/two_point_all_to_all.dat  (+ meta.json)
```

Each `meta.json` records the geometry $(L_x,L_y,T_x,T_y)$, the
couplings $(K_1,K_2,K_3)$, the recovered $\beta_c$, and the production
trajectory count. Inspect any one of them after the run to confirm:

* **TEST**: `k1 ≈ 5.0652`, `k2 ≈ 7.7429`, `k3 = 1.0`,
  `beta_c → 0.0628` as $L$ grows.
* **REF**: `k1 = k2 = k3 = 1.0`, `beta_c → 0.06283` as $\alpha$ grows
  (independent of $\alpha$ since $K=1$ is the isotropic critical line).

### 10.5 Run the matching algorithm

```bash
python match_456_twisted_vs_untwisted.py
```

Console output (verbatim from the working pipeline):

```
Ref alphas : [1, 2, 3, 4]
Test sizes : [8, 16, 24, 32, 48, 64]
Using largest TEST L = 64, largest REF α = 4
  side 4: ref_c2 A_r=0.2429   test_c1 A_t=0.2382   A_t/A_r=0.9807
  side 5: ref_c0 A_r=0.2170   test_c0 A_t=0.2128   A_t/A_r=0.9805
  side 6: ref_c1 A_r=0.2003   test_c2 A_t=0.1952   A_t/A_r=0.9747
A_test/A_ref ratios per side : [0.9807, 0.9805, 0.9747]
  mean   = 0.9787
  spread = ±0.0028  (0.28% of mean)
  → if << 1% the two lattices realise the SAME continuum CFT.
Wrote results/_fss_456/match_twisted_vs_untwisted.png
```

Figure: top row = $G_s(t)\sin(\pi t)^{1/4}$ (CFT‑flat amplitude),
bottom row = $R_s(t)=G_s(t)/G_s(1/2)$ (Z‑invariant overlay).
Three columns = sides 4, 5, 6.

### 10.6 Acceptance criteria

A successful end‑to‑end build satisfies all of:

1. `bin/ising_tri_twisted_parallelogram` built without warnings of
   substance (sign‑compare and unknown‑warning are explicitly silenced).
2. All 10 `meta.json` files present under
   `results/_fss_456/{ref,test}/...` with the expected couplings
   and finite $\beta_c$.
3. `match_456_twisted_vs_untwisted.py` reports
   $\sigma(\rho)/\bar\rho < 1\%$ (the current pipeline gives $0.28\%$).
4. The bottom‑row overlay in
   `results/_fss_456/match_twisted_vs_untwisted.png` shows test and
   ref symbols within their joint error bars at every sample $t_k$ and
   every side $s$.

If criterion 3 fails by an order of magnitude, the most likely cause
is a **simulator bond‑convention drift**: re‑derive the bond–side
dictionary in §2.1 from the simulator source
([`src/ising_tri_twisted_parallelogram.cc`](src/ising_tri_twisted_parallelogram.cc),
search for `K1`, `K2`, `K3` in the bond‑sum loop) and update
`REF_CYCLE_SIDE` / `TEST_CYCLE_SIDE` in
[`match_456_twisted_vs_untwisted.py`](match_456_twisted_vs_untwisted.py)
accordingly.

### 10.7 Optional: full optimizer stack

The matching algorithm above is the *validation* step. To exercise the
inverse problem (recover $(r_1,r_2)$ from MC alone via CMA‑ES), see
[`README.md` §1, §7](README.md). The minimal command:

```bash
cat > /tmp/cfg_456.json <<'EOF'
{"ref_Tx": -3, "ref_Ty": 3, "test_Tx": 0, "test_Ty": 0}
EOF

python run.py \
    --config /tmp/cfg_456.json \
    --ref-Lx 13 --ref-Ly 16 --test-Lx 13 --test-Ly 16 \
    --x0 3.0 4.0 --cma-sigma0 3.0 \
    --max-evals 150 --multidonor-2pass \
    --n-workers 2 --save-frames \
    --run-name 4_5_6_triangle --no-dashboard
```

A converged run reports a best $(r_1,r_2)$; verify it by re‑running
[`match_456_twisted_vs_untwisted.py`](match_456_twisted_vs_untwisted.py)
after replacing the analytic truth in
[`precompute_456_fss.py`](precompute_456_fss.py) (`_R1`, `_R2`,
`_K3`) with the recovered values, regenerating the TEST family,
and re‑checking the $\rho_s$ universality.

### 10.8 Provenance check

```bash
git -C $(git rev-parse --show-toplevel) log -1 --format='%H %s'
sha256sum bin/ising_tri_twisted_parallelogram \
          results/_fss_456/ref/a4/two_point_all_to_all.dat \
          results/_fss_456/test/L64/two_point_all_to_all.dat
```

Recording these alongside any matching figure makes the result
exactly reproducible by a future collaborator.


---

## 11.  Cost-landscape preview on the matching cost  (2026-05-04)

Before plugging the §4 matching cost into a CMA-ES production run, we
replay it on a coarse $(r_1, r_2)$ grid to verify it produces a clean
basin at the analytic truth. **Result: it does not — at the current
precompute statistics. The cost design is correct (validated to 0.28%
in §4); the cached landscape is statistically under-resolved for it.**

### 11.1  Procedure

Script: [`landscape_matching_456.py`](landscape_matching_456.py).

For each $(r_1, r_2)$ on the 8×8 sub-grid $\{2, 3, …, 9\}^2$
(containing truth $(5.07, 7.74)$ at distance $\approx 0.27$ from grid
point $(5, 8)$) and each test size $L \in \{8, 16\}$:

1. Load the cached `test_data` from
   `results/_landscape/Lx{L}_Ly{L}_Tx0_Ty0/grid/r1_*.pkl`
   (`n_traj_prod = 5000`, untwisted geometry, anisotropic $K = (r_1, r_2, 1)$).
2. Build a tiled `LinearND` interpolant on $G_{\text{conn}}(m, n)$.
3. Sample along each cycle $c \in \{0, 1, 2\}$ at $t_k = k/8$, $k=1..7$.
4. Form the $Z_\sigma$-invariant log observable

   $$\mathcal{L}(t)=\log|G(t)|-\langle\log|G|\rangle_{t\in 1/8..7/8}\,. \tag{11.1}$$

5. Pair test cycle ↔ ref cycle by **physical side** (§2.2):
   `ref c0 ↔ test c0` (side 5),
   `ref c1 ↔ test c2` (side 6),
   `ref c2 ↔ test c1` (side 4).
6. $\text{cost} = \sum_{\text{sides}} \sum_k (\mathcal{L}_{\text{test}}(t_k) - \mathcal{L}_{\text{ref}}(t_k))^2$
   (un-weighted, $21$ residuals per cell).
7. FSS extrapolate cell-by-cell: $\text{cost}(L) \approx c_\infty + a/L \Rightarrow c_\infty = 2 c_{16} - c_8$.

Reference: high-stats $4{-}5{-}6$ $\alpha = 4$ data,
`results/_fss_456/ref/a4/two_point_all_to_all.dat`.

### 11.2  Why $\mathcal{L}(t)$ instead of $R(t) = G(t)/G(1/2)$?

Both kill the global $Z_\sigma$. Empirically (and theoretically) the
log-mean version is **lower-noise**: dividing by a single sample
$G(1/2)$ injects the noise of *that* sample into every ratio, while
subtracting the geometric mean over 7 samples uses
$\sqrt{7} \times$ less of any one sample's variance. This was the
final cost choice in `_test_cost`.

### 11.3  Findings (figure: `match_cost_grid.png`)

> **Pearson(cost, distance-to-truth):**
> $L=8: \rho = -0.04$,  $L=16: \rho = +0.01$,  $L\to\infty: \rho = -0.01$.

Interpretation: the cost on this precompute is **MC noise**; near-truth
and far-from-truth cells are statistically indistinguishable. The
top-5 lowest-cost cells per panel are scattered across the plane and
include both near-truth points (e.g. $L=16$ cell $(5,6)$ at
$d_{\text{truth}}=1.74$) and far-from-truth points (e.g. $L=16$ argmin
$(9,8)$ at $d_{\text{truth}}=3.94$).

The full root-cause analysis, noise-floor calibration, and production
recommendations (precompute knobs, table of required statistics) are
in [`LANDSCAPE.md` §8](LANDSCAPE.md).

### 11.4  Recommendations for the production run

The matching-cost design is keep-as-is (§4 is final). What changes:

| knob          | diagnostic precompute | production precompute |
|---------------|-----------------------|-----------------------|
| `n_traj_prod` | 5 000                 | **≥ 50 000**          |
| FSS sizes     | 2 (L=8, 16)           | **≥ 4** (e.g. 8, 16, 24, 32) |
| grid step     | 0.5                   | 0.5 near truth        |
| extrapolant   | $c_\infty + a/L$ (clipped at 0) | jackknife multi-$L$ fit, no clip |

The order-of-magnitude estimate
$\text{noise floor} \sim 21 \cdot (1/\sqrt{n_{\text{traj}}})^2$
puts the cost noise floor at $\approx 4 \times 10^{-3}$ for
$n_{\text{traj}} = 5\,000$ and $\approx 4 \times 10^{-4}$ for
$50\,000$ — the latter is below the observed truth-cell value of
$0.002$ at $L=16$, so the basin should appear at one decade more
statistics.

### 11.5  Status as a "production seed"

This experiment is a **negative result with high information value**:

- Confirms the algo.md §4 cost design works in principle (already
  validated separately to 0.28%).
- Quantifies the precompute statistics required to use it as a CMA-ES
  initialiser (§11.4 table).
- Provides a reusable replay script + noise-floor calibration figure
  that future production precomputes can be compared against.

The next concrete step is to **build a dedicated higher-stats
precompute** following the §11.4 recipe, then re-run
`landscape_matching_456.py` against it with `TEST_SIZES` extended to
4 sizes and the 1/L extrapolation upgraded to a multi-L jackknife fit.
