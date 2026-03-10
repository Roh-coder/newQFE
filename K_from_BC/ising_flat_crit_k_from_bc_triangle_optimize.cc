// ising_flat_crit_k_from_bc_triangle_optimize.cc
//
// Adaptive coupling search for triangular-lattice Ising with twisted BC.
// Targets isotropy of connected two-point correlators along 3 lattice directions:
//   base:  (dx,dy) = (1, 0)
//   left:  (dx,dy) = (0, 1)
//   right: (dx,dy) = (1,-1)
//
// Score = RMS of pairwise z-differences after linear interpolation of correlators
// and jackknife errors over one period for each direction.

#include <getopt.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include "ising.h"
#include "statistics.h"

namespace {

struct CorrPoint {
  double mean = 0.0;
  double err = 0.0;
};

struct EvalResult {
  double score = std::numeric_limits<double>::infinity();
  double beta = 0.0;
  double m_mean = 0.0;
  std::vector<CorrPoint> base;
  std::vector<CorrPoint> left;
  std::vector<CorrPoint> right;
};

int Mod(int a, int n) {
  int m = a % n;
  return (m < 0) ? (m + n) : m;
}

int FloorDiv(int a, int b) {
  assert(b > 0);
  int q = a / b;
  int r = a % b;
  if (r != 0 && ((r < 0) != (b < 0))) q--;
  return q;
}

int CanonicalSiteTwistedX(int x, int y, int Nx, int Ny, int twist_shift) {
  int qx = FloorDiv(x, Nx);
  int xr = x - qx * Nx;
  int yr = y + qx * twist_shift;
  yr = Mod(yr, Ny);
  return yr * Nx + xr;
}

int PeriodForDirection(int dx, int dy, int Nx, int Ny, int twist_shift) {
  const int s0 = CanonicalSiteTwistedX(0, 0, Nx, Ny, twist_shift);
  int x = 0;
  int y = 0;
  const int nmax = 4 * Nx * Ny;
  for (int n = 1; n <= nmax; n++) {
    x += dx;
    y += dy;
    int s = CanonicalSiteTwistedX(x, y, Nx, Ny, twist_shift);
    if (s == s0) return n;
  }
  return Nx * Ny;
}

void InitTriangleTwistedX(QfeLattice* lattice, int Nx, int Ny, double wt1,
                          double wt2, double wt3, int twist_shift) {
  lattice->ResizeSites(Nx * Ny);
  lattice->vol = double(Nx * Ny);

  for (int s = 0; s < lattice->n_sites; s++) {
    lattice->sites[s].wt = 1.0;
    lattice->sites[s].nn = 0;
    lattice->sites[s].id = 0;
  }

  lattice->links.clear();
  lattice->faces.clear();
  lattice->cells.clear();
  lattice->n_links = 0;
  lattice->n_faces = 0;
  lattice->n_cells = 0;

  for (int y = 0; y < Ny; y++) {
    for (int x = 0; x < Nx; x++) {
      int s = CanonicalSiteTwistedX(x, y, Nx, Ny, twist_shift);
      int sx = CanonicalSiteTwistedX(x + 1, y, Nx, Ny, twist_shift);
      int sy = CanonicalSiteTwistedX(x, y + 1, Nx, Ny, twist_shift);
      int sxy = CanonicalSiteTwistedX(x + 1, y + 1, Nx, Ny, twist_shift);

      lattice->AddFace(s, sx, sy);
      lattice->AddFace(sx, sxy, sy);

      int l1 = lattice->FindLink(s, sx);
      int l2 = lattice->FindLink(s, sy);
      int l3 = lattice->FindLink(sx, sy);
      assert(l1 >= 0 && l2 >= 0 && l3 >= 0);
      lattice->links[l1].wt = wt1;
      lattice->links[l2].wt = wt2;
      lattice->links[l3].wt = wt3;
    }
  }

  lattice->UpdateDistinct();
}

double TriCritStep(double k1, double k2, double k3, double beta) {
  double p1 = exp(-2.0 * beta * (k2 + k3));
  double p2 = exp(-2.0 * beta * (k3 + k1));
  double p3 = exp(-2.0 * beta * (k1 + k2));

  double f = p1 + p2 + p3 - 1.0;
  double df = -2.0 * (p1 * (k2 + k3) + p2 * (k3 + k1) + p3 * (k1 + k2));
  return f / df;
}

bool SolveK3OnCriticalSurface(double k1, double k2, double beta, double* k3) {
  const double c = std::exp(-2.0 * beta * (k1 + k2));
  const double a = std::exp(-2.0 * beta * k1);
  const double b = std::exp(-2.0 * beta * k2);
  const double num = 1.0 - c;
  const double den = a + b;
  if (beta <= 0.0 || den <= 0.0 || num <= 0.0) return false;
  const double x = num / den;  // x = exp(-2 beta k3)
  if (x <= 0.0) return false;
  *k3 = -0.5 * std::log(x) / beta;
  return std::isfinite(*k3);
}

double FerroSurfaceMargin(double k1, double k2, double beta) {
  const double a = std::exp(-2.0 * beta * k1);
  const double b = std::exp(-2.0 * beta * k2);
  const double c = std::exp(-2.0 * beta * (k1 + k2));
  // k3 > 0 is equivalent to a + b + c - 1 > 0 at fixed beta.
  return a + b + c - 1.0;
}

double MaxK2FerroForK1(double k1, double beta, double k2_lo, double k2_hi) {
  if (k2_hi <= k2_lo) return k2_lo;
  if (FerroSurfaceMargin(k1, k2_lo, beta) <= 0.0) return k2_lo;
  if (FerroSurfaceMargin(k1, k2_hi, beta) > 0.0) return k2_hi;

  double lo = k2_lo;
  double hi = k2_hi;
  for (int it = 0; it < 80; it++) {
    const double mid = 0.5 * (lo + hi);
    if (FerroSurfaceMargin(k1, mid, beta) > 0.0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

double FindCrit(double k1, double k2, double k3) {
  double k_mean = (k1 + k2 + k3) / 3.0;
  k1 /= k_mean;
  k2 /= k_mean;
  k3 /= k_mean;

  double beta = 0.267949192431123;
  for (int i = 0; i < 100; i++) beta -= TriCritStep(k1, k2, k3, beta);
  return beta / k_mean;
}

double MeanFromSamples(const std::vector<double>& v) {
  if (v.empty()) return 0.0;
  double sum = 0.0;
  for (double x : v) sum += x;
  return sum / double(v.size());
}

void JackknifeConnectedCorr(const std::vector<double>& corr_samples,
                            const std::vector<double>& mag_samples,
                            double* mean_out, double* err_out) {
  int n = int(corr_samples.size());
  assert(n == int(mag_samples.size()));
  if (n <= 1) {
    *mean_out = 0.0;
    *err_out = 0.0;
    return;
  }

  double sum_corr = 0.0;
  double sum_mag = 0.0;
  for (int i = 0; i < n; i++) {
    sum_corr += corr_samples[i];
    sum_mag += mag_samples[i];
  }

  double mean_corr = sum_corr / double(n);
  double mean_mag = sum_mag / double(n);
  double connected_mean = mean_corr - mean_mag * mean_mag;

  double sumsq = 0.0;
  for (int i = 0; i < n; i++) {
    double corr_leave = (sum_corr - corr_samples[i]) / double(n - 1);
    double mag_leave = (sum_mag - mag_samples[i]) / double(n - 1);
    double conn_leave = corr_leave - mag_leave * mag_leave;
    double d = conn_leave - connected_mean;
    sumsq += d * d;
  }

  *mean_out = connected_mean;
  *err_out = std::sqrt((double(n) - 1.0) / double(n) * sumsq);
}

double InterpPeriodic(const std::vector<CorrPoint>& c, double t, bool use_err) {
  const int n = int(c.size());
  if (n <= 0) return 0.0;
  double u = t - std::floor(t);
  double pos = u * double(n);
  int i0 = int(std::floor(pos)) % n;
  double a = pos - std::floor(pos);
  int i1 = (i0 + 1) % n;
  double v0 = use_err ? c[i0].err : c[i0].mean;
  double v1 = use_err ? c[i1].err : c[i1].mean;
  return (1.0 - a) * v0 + a * v1;
}

double PairwiseRmsScore(const std::vector<CorrPoint>& a,
                        const std::vector<CorrPoint>& b,
                        const std::vector<CorrPoint>& c, int dense_points) {
  double sum = 0.0;
  int cnt = 0;
  for (int m = 0; m < dense_points; m++) {
    double t = (double(m) + 0.5) / double(dense_points);

    double am = InterpPeriodic(a, t, false);
    double bm = InterpPeriodic(b, t, false);
    double cm = InterpPeriodic(c, t, false);

    double ae = InterpPeriodic(a, t, true);
    double be = InterpPeriodic(b, t, true);
    double ce = InterpPeriodic(c, t, true);

    double eab = std::sqrt(ae * ae + be * be);
    double eac = std::sqrt(ae * ae + ce * ce);
    double ebc = std::sqrt(be * be + ce * ce);

    if (eab <= 1.0e-14) eab = 1.0e-14;
    if (eac <= 1.0e-14) eac = 1.0e-14;
    if (ebc <= 1.0e-14) ebc = 1.0e-14;

    double zab = (am - bm) / eab;
    double zac = (am - cm) / eac;
    double zbc = (bm - cm) / ebc;

    sum += zab * zab;
    sum += zac * zac;
    sum += zbc * zbc;
    cnt += 3;
  }
  return (cnt > 0) ? std::sqrt(sum / double(cnt)) : 0.0;
}

EvalResult RunSingleEvaluation(int Nx, int Ny, int twist_shift, double k1,
                               double k2, double k3, double beta_mult,
                               double fixed_beta,
                               int n_therm, int n_traj, int n_skip,
                               int n_wolff, int n_metropolis,
                               bool verbose_traj) {
  EvalResult out;

  const int dx[3] = {1, 0, 1};
  const int dy[3] = {0, 1, -1};
  int period[3] = {
      PeriodForDirection(dx[0], dy[0], Nx, Ny, twist_shift),
      PeriodForDirection(dx[1], dy[1], Nx, Ny, twist_shift),
      PeriodForDirection(dx[2], dy[2], Nx, Ny, twist_shift),
  };

  if (verbose_traj) {
    std::printf("  periods(base,left,right) = (%d,%d,%d)\n", period[0],
                period[1], period[2]);
  }

  QfeLattice lattice;
  InitTriangleTwistedX(&lattice, Nx, Ny, k1, k2, k3, twist_shift);

  double beta_use = 0.0;
  if (fixed_beta > 0.0) {
    beta_use = fixed_beta;
  } else {
    beta_use = FindCrit(k1, k2, k3) * beta_mult;
  }

  QfeIsing field(&lattice, beta_use);
  field.HotStart();
  out.beta = field.beta;

  std::vector<double> mag;
  std::vector<std::vector<std::vector<double>>> corr_samples(3);
  for (int d = 0; d < 3; d++) corr_samples[d].resize(period[d]);

  for (int n = 0; n < (n_therm + n_traj); n++) {
    for (int j = 0; j < n_wolff; j++) field.WolffUpdate();
    for (int j = 0; j < n_metropolis; j++) field.Metropolis();

    if (n < n_therm || (n % n_skip) != 0) continue;

    std::vector<std::vector<double>> corr_sum(3);
    for (int d = 0; d < 3; d++) corr_sum[d].assign(period[d], 0.0);

    for (int y = 0; y < Ny; y++) {
      for (int x = 0; x < Nx; x++) {
        int s = CanonicalSiteTwistedX(x, y, Nx, Ny, twist_shift);
        double spin_s = field.spin[s];
        for (int d = 0; d < 3; d++) {
          for (int r = 0; r < period[d]; r++) {
            int t = CanonicalSiteTwistedX(x + r * dx[d], y + r * dy[d], Nx, Ny,
                                          twist_shift);
            corr_sum[d][r] += spin_s * field.spin[t];
          }
        }
      }
    }

    const double norm = double(Nx * Ny);
    for (int d = 0; d < 3; d++) {
      for (int r = 0; r < period[d]; r++) {
        corr_samples[d][r].push_back(corr_sum[d][r] / norm);
      }
    }

    mag.push_back(field.MeanSpin());
    if (verbose_traj) {
      std::printf("    meas %04zu m=%+.8f\n", mag.size(), mag.back());
    }
  }

  out.m_mean = MeanFromSamples(mag);
  out.base.resize(period[0]);
  out.left.resize(period[1]);
  out.right.resize(period[2]);

  for (int r = 0; r < period[0]; r++) {
    JackknifeConnectedCorr(corr_samples[0][r], mag, &out.base[r].mean,
                           &out.base[r].err);
  }
  for (int r = 0; r < period[1]; r++) {
    JackknifeConnectedCorr(corr_samples[1][r], mag, &out.left[r].mean,
                           &out.left[r].err);
  }
  for (int r = 0; r < period[2]; r++) {
    JackknifeConnectedCorr(corr_samples[2][r], mag, &out.right[r].mean,
                           &out.right[r].err);
  }

  int dense = 4 * std::max(period[0], std::max(period[1], period[2]));
  if (dense < 64) dense = 64;
  out.score = PairwiseRmsScore(out.base, out.left, out.right, dense);

  return out;
}

std::vector<double> UniformGrid(double lo, double hi, int n) {
  std::vector<double> out;
  out.reserve(n);
  if (n <= 1) {
    out.push_back(0.5 * (lo + hi));
    return out;
  }
  for (int i = 0; i < n; i++) {
    double t = double(i) / double(n - 1);
    out.push_back((1.0 - t) * lo + t * hi);
  }
  return out;
}

void RefineBoundsFromBestIndexBoundaryAware(const std::vector<double>& grid,
                                            int best_idx, double* lo,
                                            double* hi, double hard_lo,
                                            double hard_hi) {
  int n = int(grid.size());
  if (n < 2) {
    *lo = *hi = grid.empty() ? 0.0 : grid[0];
    return;
  }

  const double cur_lo = *lo;
  const double cur_hi = *hi;
  const double width = std::max(cur_hi - cur_lo, 1.0e-12);

  if (best_idx <= 0) {
    // Best at lower boundary: shift window downward and slightly expand.
    const double new_width = 1.2 * width;
    double nlo = cur_lo - 0.6 * width;
    double nhi = nlo + new_width;
    if (nlo < hard_lo) {
      nlo = hard_lo;
      nhi = std::min(hard_hi, nlo + new_width);
    }
    *lo = nlo;
    *hi = std::max(nhi, nlo + 1.0e-12);
    return;
  }

  if (best_idx >= n - 1) {
    // Best at upper boundary: shift window upward and slightly expand.
    const double new_width = 1.2 * width;
    double nhi = cur_hi + 0.6 * width;
    double nlo = nhi - new_width;
    if (nhi > hard_hi) {
      nhi = hard_hi;
      nlo = std::max(hard_lo, nhi - new_width);
    }
    *lo = nlo;
    *hi = std::max(nhi, nlo + 1.0e-12);
    return;
  }

  // Interior best: standard local 2-point refinement around best cell.
  int left = best_idx;
  int right = best_idx + 1;
  *lo = std::min(grid[left], grid[right]);
  *hi = std::max(grid[left], grid[right]);
  if (*lo < hard_lo) *lo = hard_lo;
  if (*hi > hard_hi) *hi = hard_hi;
}

void WriteBestDat(const char* tag, int Nx, int Ny, int shift, double k1,
                  double k2, double k3, const EvalResult& r) {
  mkdir("ising_flat_crit_k_from_bc", 0755);
  char path[512];
  std::snprintf(path, sizeof(path),
                "ising_flat_crit_k_from_bc/triangle_opt_%s_Nx%d_Ny%d_s%d_"
                "k%.6f_%.6f_%.6f.dat",
                tag, Nx, Ny, shift, k1, k2, k3);

  FILE* f = std::fopen(path, "w");
  assert(f != nullptr);
  for (int i = 0; i < int(r.base.size()); i++) {
    std::fprintf(f, "10 %04d %.12e %.12e\n", i, r.base[i].mean, r.base[i].err);
  }
  for (int i = 0; i < int(r.left.size()); i++) {
    std::fprintf(f, "11 %04d %.12e %.12e\n", i, r.left[i].mean, r.left[i].err);
  }
  for (int i = 0; i < int(r.right.size()); i++) {
    std::fprintf(f, "12 %04d %.12e %.12e\n", i, r.right[i].mean,
                 r.right[i].err);
  }
  std::fclose(f);

  std::printf("output_dat: %s\n", path);
}

}  // namespace

int main(int argc, char* argv[]) {
  int Nx = 32;
  int Ny = 32;
  int shift = 16;

  double k1_center = 1.0;
  double k2_center = 1.0;
  double k3_center = 1.0;

  double k1_min = 0.0;
  double k1_max = 0.0;
  double k2_min = 0.0;
  double k2_max = 0.0;
  double k3_min = 0.0;
  double k3_max = 0.0;

  int rounds = 3;
  int grid_n = 4;

  double beta_mult = 1.0;
  double fixed_beta = 0.0;
  int constrain_surface = 0;
  int ferro_only = 0;
  int n_therm = 400;
  int n_traj = 1200;
  int n_skip = 10;
  int n_wolff = 3;
  int n_metropolis = 5;
  int verbose_traj = 0;

  const struct option long_options[] = {
      {"n_x", required_argument, 0, 'X'},
      {"n_y", required_argument, 0, 'Y'},
      {"shift", required_argument, 0, 'S'},
      {"k1", required_argument, 0, 'a'},
      {"k2", required_argument, 0, 'b'},
      {"k3", required_argument, 0, 'c'},
      {"k1_min", required_argument, 0, 'd'},
      {"k1_max", required_argument, 0, 'e'},
      {"k2_min", required_argument, 0, 'f'},
      {"k2_max", required_argument, 0, 'g'},
      {"k3_min", required_argument, 0, 'i'},
      {"k3_max", required_argument, 0, 'j'},
      {"rounds", required_argument, 0, 'R'},
      {"grid_n", required_argument, 0, 'G'},
      {"beta_mult", required_argument, 0, 'm'},
      {"fixed_beta", required_argument, 0, 'B'},
      {"constrain_surface", required_argument, 0, 'C'},
      {"ferro_only", required_argument, 0, 'F'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'p'},
      {"verbose_traj", required_argument, 0, 'v'},
      {0, 0, 0, 0}};

  const char* short_options =
      "X:Y:S:a:b:c:d:e:f:g:i:j:R:G:m:B:C:F:h:t:s:w:p:v:";

  while (true) {
    int opt_index = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &opt_index);
    if (c == -1) break;

    switch (c) {
      case 'X': Nx = std::atoi(optarg); break;
      case 'Y': Ny = std::atoi(optarg); break;
      case 'S': shift = std::atoi(optarg); break;
      case 'a': k1_center = std::stod(optarg); break;
      case 'b': k2_center = std::stod(optarg); break;
      case 'c': k3_center = std::stod(optarg); break;
      case 'd': k1_min = std::stod(optarg); break;
      case 'e': k1_max = std::stod(optarg); break;
      case 'f': k2_min = std::stod(optarg); break;
      case 'g': k2_max = std::stod(optarg); break;
      case 'i': k3_min = std::stod(optarg); break;
      case 'j': k3_max = std::stod(optarg); break;
      case 'R': rounds = std::atoi(optarg); break;
      case 'G': grid_n = std::atoi(optarg); break;
      case 'm': beta_mult = std::stod(optarg); break;
      case 'B': fixed_beta = std::stod(optarg); break;
      case 'C': constrain_surface = std::atoi(optarg); break;
      case 'F': ferro_only = std::atoi(optarg); break;
      case 'h': n_therm = std::atoi(optarg); break;
      case 't': n_traj = std::atoi(optarg); break;
      case 's': n_skip = std::atoi(optarg); break;
      case 'w': n_wolff = std::atoi(optarg); break;
      case 'p': n_metropolis = std::atoi(optarg); break;
      case 'v': verbose_traj = std::atoi(optarg); break;
      default: break;
    }
  }

  if (grid_n < 2) grid_n = 2;
  if (rounds < 1) rounds = 1;

  // Default initial search box if user did not specify explicit bounds.
  if (!(k1_max > k1_min)) {
    k1_min = 0.6 * k1_center;
    k1_max = 1.4 * k1_center;
  }
  if (!(k2_max > k2_min)) {
    k2_min = 0.6 * k2_center;
    k2_max = 1.4 * k2_center;
  }
  if (!(k3_max > k3_min)) {
    k3_min = 0.6 * k3_center;
    k3_max = 1.4 * k3_center;
  }

  const double hard_k1_min = k1_min;
  const double hard_k1_max = k1_max;
  const double hard_k2_min = k2_min;
  const double hard_k2_max = k2_max;
  const double hard_k3_min = k3_min;
  const double hard_k3_max = k3_max;

  std::printf("=== triangle coupling optimizer ===\n");
  std::printf("Nx=%d Ny=%d shift=%d\n", Nx, Ny, shift);
  std::printf("MC params: n_therm=%d n_traj=%d n_skip=%d n_wolff=%d n_metropolis=%d\n",
              n_therm, n_traj, n_skip, n_wolff, n_metropolis);
  if (fixed_beta > 0.0) {
    std::printf("beta mode: fixed beta=%.10f\n", fixed_beta);
  } else {
    std::printf("beta mode: critical beta * beta_mult, beta_mult=%.10f\n",
                beta_mult);
  }
  std::printf("search rounds=%d grid_n=%d (points per round=%d)\n", rounds,
              grid_n, grid_n * grid_n * grid_n);
  if (constrain_surface) {
    if (fixed_beta <= 0.0) fixed_beta = 1.0;
    std::printf("constraint mode: critical surface at fixed beta=%.10f\n",
                fixed_beta);
  }
  if (ferro_only) {
    std::printf("coupling constraint: ferromagnetic-only (k1,k2,k3 > 0)\n");
  }

  EvalResult best_global;
  double best_k1 = 0.0;
  double best_k2 = 0.0;
  double best_k3 = 0.0;

  for (int round = 1; round <= rounds; round++) {
    double rk1_min = k1_min;
    double rk1_max = k1_max;
    double rk2_min = k2_min;
    double rk2_max = k2_max;
    if (constrain_surface && ferro_only && fixed_beta > 0.0) {
      // Restrict round window to the ferro-admissible region k3>0.
      rk2_max = std::min(rk2_max,
                         MaxK2FerroForK1(rk1_max, fixed_beta, rk2_min, rk2_max));
      if (rk2_max <= rk2_min) rk2_max = rk2_min + 1.0e-6;
    }

    std::printf("\n--- round %d/%d ---\n", round, rounds);
    std::printf("range k1=[%.9f, %.9f] k2=[%.9f, %.9f] k3=[%.9f, %.9f]\n",
                rk1_min, rk1_max, rk2_min, rk2_max, k3_min, k3_max);

    std::vector<double> g1 = UniformGrid(rk1_min, rk1_max, grid_n);
    std::vector<double> g2 = UniformGrid(rk2_min, rk2_max, grid_n);
    std::vector<double> g3 = UniformGrid(k3_min, k3_max, grid_n);

    double best_round_score = std::numeric_limits<double>::infinity();
    int bi = 0, bj = 0, bk = 0;
    EvalResult best_round_result;

    int eval_id = 0;
    int total_eval = constrain_surface ? (grid_n * grid_n)
                                       : (grid_n * grid_n * grid_n);
    for (int i = 0; i < grid_n; i++) {
      for (int j = 0; j < grid_n; j++) {
        const int k_end = constrain_surface ? 1 : grid_n;
        for (int k = 0; k < k_end; k++) {
          eval_id++;
          double k1 = g1[i], k2 = g2[j], k3 = constrain_surface ? 0.0 : g3[k];
          if (constrain_surface) {
            if (!SolveK3OnCriticalSurface(k1, k2, fixed_beta, &k3)) {
              std::printf("eval %03d/%03d: k1=%.9f k2=%.9f invalid on surface, skipping\n",
                          eval_id, total_eval, k1, k2);
              continue;
            }
          }
          if (ferro_only && (k1 <= 0.0 || k2 <= 0.0 || k3 <= 0.0)) {
            std::printf("eval %03d/%03d: k=(%.9f, %.9f, %.9f) violates ferro-only, skipping\n",
                        eval_id, total_eval, k1, k2, k3);
            continue;
          }
          std::printf("eval %03d/%03d: k=(%.9f, %.9f, %.9f) ...\n", eval_id,
                      total_eval, k1, k2, k3);
          EvalResult r = RunSingleEvaluation(
              Nx, Ny, shift, k1, k2, k3, beta_mult, fixed_beta, n_therm,
              n_traj, n_skip, n_wolff, n_metropolis, verbose_traj != 0);
          std::printf("  -> score=%.8f beta=%.10f <m>=%+.7f\n", r.score, r.beta,
                      r.m_mean);

          if (r.score < best_round_score) {
            best_round_score = r.score;
            bi = i;
            bj = j;
            bk = k;
            best_round_result = r;
          }
          if (r.score < best_global.score) {
            best_global = r;
            best_k1 = k1;
            best_k2 = k2;
            best_k3 = k3;
            std::printf(
                "  ** new global best: score=%.8f at k=(%.9f, %.9f, %.9f)\n",
                best_global.score, best_k1, best_k2, best_k3);
          }
        }
      }
    }

    std::printf("round %d best score=%.8f at indices (%d,%d,%d)\n", round,
                best_round_score, bi, bj, bk);

    if (round < rounds) {
      RefineBoundsFromBestIndexBoundaryAware(g1, bi, &k1_min, &k1_max,
                                             hard_k1_min, hard_k1_max);
      RefineBoundsFromBestIndexBoundaryAware(g2, bj, &k2_min, &k2_max,
                                             hard_k2_min, hard_k2_max);
      if (!constrain_surface) {
        RefineBoundsFromBestIndexBoundaryAware(g3, bk, &k3_min, &k3_max,
                                               hard_k3_min, hard_k3_max);
      }
    }
  }

  std::printf("\n=== final best ===\n");
  std::printf("best k1=%.10f\n", best_k1);
  std::printf("best k2=%.10f\n", best_k2);
  std::printf("best k3=%.10f\n", best_k3);
  std::printf("best score=%.10f\n", best_global.score);
  std::printf("best beta=%.10f\n", best_global.beta);

  char tag[128];
  std::snprintf(tag, sizeof(tag), "bc_Nx%d_Ny%d_s%d", Nx, Ny, shift);
  WriteBestDat(tag, Nx, Ny, shift, best_k1, best_k2, best_k3, best_global);

  return 0;
}
