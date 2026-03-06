// ising_flat_crit_k_from_bc_triangle_sides.cc
//
// Direction-focused anisotropy test on a triangular lattice with twisted BC.
// Measures only two side directions in the imposed real geometry:
//   1) lattice x direction
//   2) one oblique lattice side (x +/- y)
//
// The x-wrap twist is implemented as:
//   (x + Nx, y) ~ (x, y + twist_shift)
// so the lattice itself is a skewed torus.

#include <getopt.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include "ising.h"
#include "statistics.h"

namespace {

struct RightTriangleRule {
  // Geometry-scale inputs interpreted as right triangle legs.
  double lx = 1.0;
  double ly = 1.0;
  // Derived couplings for triangular-lattice bond directions.
  // We map:
  //   k1 -> x bond
  //   k2 -> +y bond
  //   k3 -> (y-x) oblique bond
  double k1 = 1.0;
  double k2 = 1.0;
  double k3 = 1.0;
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

int GcdInt(int a, int b) {
  a = std::abs(a);
  b = std::abs(b);
  while (b != 0) {
    int t = a % b;
    a = b;
    b = t;
  }
  return (a == 0) ? 1 : a;
}

int CanonicalSiteTwistedX(int x, int y, int Nx, int Ny, int twist_shift) {
  int qx = FloorDiv(x, Nx);
  int xr = x - qx * Nx;
  int yr = y + qx * twist_shift;
  yr = Mod(yr, Ny);
  return yr * Nx + xr;
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

double FindCrit(double k1, double k2, double k3) {
  double k_mean = (k1 + k2 + k3) / 3.0;
  k1 /= k_mean;
  k2 /= k_mean;
  k3 /= k_mean;

  double beta = 0.267949192431123;  // 2 - sqrt(3)
  for (int i = 0; i < 100; i++) beta -= TriCritStep(k1, k2, k3, beta);
  return beta / k_mean;
}

RightTriangleRule BuildRuleFromRightTriangle(double lx, double ly) {
  RightTriangleRule rule;
  rule.lx = lx;
  rule.ly = ly;

  // l*_i / l_i prescription encoded through sinh(2k_i).
  // Chosen mapping for a right triangle with hypotenuse l_h:
  //   sinh(2k_x)          = ly / lx
  //   sinh(2k_oblique)    = lx / ly
  //   sinh(2k_cross-link) = lx*ly / l_h^2
  // This keeps k_x and k_oblique as the two side-direction controls.
  const double lh2 = lx * lx + ly * ly;
  const double lh = std::sqrt(lh2);
  assert(lx > 0.0 && ly > 0.0 && lh > 0.0);

  const double sx = ly / lx;
  const double sobl = lx / ly;
  const double scross = (lx * ly) / lh2;

  rule.k1 = 0.5 * std::asinh(sx);
  rule.k2 = 0.5 * std::asinh(scross);
  rule.k3 = 0.5 * std::asinh(sobl);
  return rule;
}

double MeanFromSamples(const std::vector<double>& v) {
  if (v.empty()) return 0.0;
  double sum = 0.0;
  for (double x : v) sum += x;
  return sum / double(v.size());
}

double JackknifeErrorMean(const std::vector<double>& v) {
  int n = int(v.size());
  if (n <= 1) return 0.0;
  double total = 0.0;
  for (double x : v) total += x;
  double mean = total / double(n);
  double sumsq = 0.0;
  for (double x : v) {
    double leave = (total - x) / double(n - 1);
    double d = leave - mean;
    sumsq += d * d;
  }
  return std::sqrt((double(n) - 1.0) / double(n) * sumsq);
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

}  // namespace

int main(int argc, char* argv[]) {
  int Nx = 48;
  int Ny = 64;

  // Right-triangle geometry for sinh(2k)=l*/l mapping.
  double lx = 1.0;
  double ly = 0.75;
  bool use_sinh_rule = true;

  // Manual override (used when --use_sinh_rule=0).
  double k1 = 1.0;
  double k2 = 1.0;
  double k3 = 1.0;

  // Twisted BC: (x+Nx, y) ~ (x, y+twist_shift).
  int twist_shift = 0;

  // oblique_sign = +1 uses x+y; -1 uses y-x.
  int oblique_sign = -1;

  double beta_mult = 1.0;
  int n_therm = 2000;
  int n_traj = 20000;
  int n_skip = 20;
  int n_wolff = 3;
  int n_metropolis = 5;

  const struct option long_options[] = {
      {"n_x", required_argument, 0, 'X'},
      {"n_y", required_argument, 0, 'Y'},
      {"lx", required_argument, 0, 'u'},
      {"ly", required_argument, 0, 'v'},
      {"use_sinh_rule", required_argument, 0, 'r'},
      {"k1", required_argument, 0, 'a'},
      {"k2", required_argument, 0, 'b'},
      {"k3", required_argument, 0, 'c'},
      {"twist_shift", required_argument, 0, 'T'},
      {"oblique_sign", required_argument, 0, 'g'},
      {"beta_mult", required_argument, 0, 'm'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'e'},
      {0, 0, 0, 0}};

  const char* short_options = "X:Y:u:v:r:a:b:c:T:g:m:h:t:s:w:e:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'X': Nx = atoi(optarg); break;
      case 'Y': Ny = atoi(optarg); break;
      case 'u': lx = std::stod(optarg); break;
      case 'v': ly = std::stod(optarg); break;
      case 'r': use_sinh_rule = (atoi(optarg) != 0); break;
      case 'a': k1 = std::stod(optarg); break;
      case 'b': k2 = std::stod(optarg); break;
      case 'c': k3 = std::stod(optarg); break;
      case 'T': twist_shift = atoi(optarg); break;
      case 'g': oblique_sign = (atoi(optarg) >= 0) ? 1 : -1; break;
      case 'm': beta_mult = std::stod(optarg); break;
      case 'h': n_therm = atoi(optarg); break;
      case 't': n_traj = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'w': n_wolff = atoi(optarg); break;
      case 'e': n_metropolis = atoi(optarg); break;
      default: break;
    }
  }

  if (use_sinh_rule) {
    RightTriangleRule rule = BuildRuleFromRightTriangle(lx, ly);
    k1 = rule.k1;
    k2 = rule.k2;
    k3 = rule.k3;
  }

  printf("Nx Ny: %d %d\n", Nx, Ny);
  printf("lx ly (right-triangle legs): %.12f %.12f\n", lx, ly);
  printf("use_sinh_rule: %d\n", int(use_sinh_rule));
  printf("k1 k2 k3: %.12f %.12f %.12f\n", k1, k2, k3);
  printf("twisted BC: (x+Nx,y)~(x,y+%d)\n", twist_shift);
  printf("real-y direction uses lattice step: (dx,dy)=(1,%d)\n", oblique_sign);

  int effective_shift = Mod(twist_shift, Ny);
  int close_numer = Nx * Ny;
  int close_denom = GcdInt(Ny, Nx - oblique_sign * effective_shift);
  int oblique_period = close_numer / close_denom;
  printf("oblique closure period under twist: %d lattice steps\n", oblique_period);

  QfeLattice lattice;
  InitTriangleTwistedX(&lattice, Nx, Ny, k1, k2, k3, twist_shift);

  QfeIsing field(&lattice, FindCrit(k1, k2, k3) * beta_mult);
  field.HotStart();

  printf("beta: %.12f\n", field.beta);
  printf("initial action: %.12f\n", field.Action());

  int r_bins_side = Nx;

  std::vector<QfeMeasReal> corr_x(r_bins_side);
  std::vector<QfeMeasReal> corr_side(r_bins_side);
  std::vector<std::vector<double>> corr_x_samples(r_bins_side);
  std::vector<std::vector<double>> corr_side_samples(r_bins_side);

  std::vector<double> mag;
  std::vector<double> action;
  QfeMeasReal cluster_size;
  QfeMeasReal accept_metropolis;

  for (int n = 0; n < (n_traj + n_therm); n++) {
    int cluster_size_sum = 0;
    for (int j = 0; j < n_wolff; j++) cluster_size_sum += field.WolffUpdate();

    double metropolis_sum = 0.0;
    for (int j = 0; j < n_metropolis; j++) metropolis_sum += field.Metropolis();

    cluster_size.Measure(double(cluster_size_sum) / double(Nx * Ny));
    accept_metropolis.Measure(metropolis_sum);

    if (n % n_skip || n < n_therm) continue;

    std::vector<double> corr_x_sum(r_bins_side, 0.0);
    std::vector<double> corr_side_sum(r_bins_side, 0.0);

    for (int y = 0; y < Ny; y++) {
      for (int x = 0; x < Nx; x++) {
        int s = CanonicalSiteTwistedX(x, y, Nx, Ny, twist_shift);
        double spin_s = field.spin[s];

        for (int r = 0; r < r_bins_side; r++) {
          int sx = CanonicalSiteTwistedX(x + r, y, Nx, Ny, twist_shift);
          int so = CanonicalSiteTwistedX(x + r, y + oblique_sign * r, Nx, Ny,
                                         twist_shift);
          corr_x_sum[r] += spin_s * field.spin[sx];
          corr_side_sum[r] += spin_s * field.spin[so];
        }
      }
    }

    double norm = double(Nx * Ny);
    for (int i = 0; i < r_bins_side; i++) {
      double vx = corr_x_sum[i] / norm;
      double vo = corr_side_sum[i] / norm;
      corr_x[i].Measure(vx);
      corr_side[i].Measure(vo);
      corr_x_samples[i].push_back(vx);
      corr_side_samples[i].push_back(vo);
    }

    action.push_back(field.Action());
    mag.push_back(field.MeanSpin());
    printf("%06d %.12f %+.12f %.4f %.4f\n", n, action.back(), mag.back(),
           accept_metropolis.last, cluster_size.last);
  }

  std::vector<double> mag_abs(mag.size());
  std::vector<double> mag2(mag.size());
  for (int i = 0; i < int(mag.size()); i++) {
    double m = mag[i];
    mag_abs[i] = std::abs(m);
    mag2[i] = m * m;
  }

  printf("accept_metropolis: %.4f\n", accept_metropolis.Mean());
  printf("cluster_size/V: %.4f\n", cluster_size.Mean());
  printf("action: %.12e (%.12e), %.4f\n", Mean(action), JackknifeMean(action),
         AutocorrTime(action));
  printf("m: %.12e (%.12e), %.4f\n", Mean(mag), JackknifeMean(mag),
         AutocorrTime(mag));
  printf("susceptibility: %.12e (%.12e)\n", Susceptibility(mag2, mag_abs),
         JackknifeSusceptibility(mag2, mag_abs));

  mkdir("ising_flat_crit_k_from_bc", 0755);
  char path[256];
  std::snprintf(path, sizeof(path),
                "ising_flat_crit_k_from_bc/triangle_sides_twist_Nx%d_Ny%d_s%d_"
                "ob%d_lx%.3f_ly%.3f_k%.3f_%.3f_%.3f.dat",
                Nx, Ny, twist_shift, oblique_sign, lx, ly, k1, k2, k3);

  FILE* file = std::fopen(path, "w");
  assert(file != nullptr);

  // type 0: raw x-direction correlator
  for (int i = 0; i < r_bins_side; i++) {
    double m = MeanFromSamples(corr_x_samples[i]);
    double e = JackknifeErrorMean(corr_x_samples[i]);
    std::fprintf(file, "0 %04d %.12e %.12e\n", i, m, e);
  }

  // type 1: raw oblique-side correlator (real-y direction)
  for (int i = 0; i < r_bins_side; i++) {
    double m = MeanFromSamples(corr_side_samples[i]);
    double e = JackknifeErrorMean(corr_side_samples[i]);
    std::fprintf(file, "1 %04d %.12e %.12e\n", i, m, e);
  }

  // type 4: connected x-direction correlator
  for (int i = 0; i < r_bins_side; i++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_x_samples[i], mag, &m_conn, &e_conn);
    std::fprintf(file, "4 %04d %.12e %.12e\n", i, m_conn, e_conn);
  }

  // type 5: connected oblique-side correlator
  for (int i = 0; i < r_bins_side; i++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_side_samples[i], mag, &m_conn, &e_conn);
    std::fprintf(file, "5 %04d %.12e %.12e\n", i, m_conn, e_conn);
  }

  std::fclose(file);

  double sum_z2 = 0.0;
  int n_z = 0;
  for (int r = 1; r < r_bins_side; r++) {
    double mx = 0.0, ex = 0.0;
    double my = 0.0, ey = 0.0;
    JackknifeConnectedCorr(corr_x_samples[r], mag, &mx, &ex);
    JackknifeConnectedCorr(corr_side_samples[r], mag, &my, &ey);
    double e = std::sqrt(ex * ex + ey * ey);
    if (e <= 0.0) continue;
    double z = (my - mx) / e;
    sum_z2 += z * z;
    n_z++;
  }
  double rms_z = (n_z > 0) ? std::sqrt(sum_z2 / double(n_z)) : 0.0;

  printf("output_dat: %s\n", path);
  printf("symmetry_check_rms_z_connected: %.6f over %d bins\n", rms_z, n_z);

  return 0;
}
