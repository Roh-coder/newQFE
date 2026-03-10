#include <getopt.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#include "ising.h"
#include "statistics.h"

namespace {

int Mod(int a, int n) {
  int r = a % n;
  return (r < 0) ? (r + n) : r;
}

int SiteIndex(int x, int y, int Nx) { return y * Nx + x; }

struct XY {
  int x;
  int y;
};

XY StepYPlus(int x, int y, int Nx, int Ny, int shift) {
  if (y + 1 < Ny) return {x, y + 1};
  return {Mod(x + shift, Nx), 0};
}

void AddUniqueLink(QfeLattice* lattice, int a, int b, double wt) {
  if (a == b) return;
  int l = lattice->FindLink(a, b);
  if (l == -1) {
    lattice->AddLink(a, b, wt);
  } else {
    lattice->links[l].wt = wt;
  }
}

void InitTwistedTriangle(QfeLattice* lattice, int Nx, int Ny, int shift,
                        double k1, double k2, double k3) {
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
      int s = SiteIndex(x, y, Nx);

      // Direction 1: +x
      int sx = SiteIndex(Mod(x + 1, Nx), y, Nx);
      AddUniqueLink(lattice, s, sx, k1);

      // Direction 2: +y with twisted wrap (x -> x + shift at y-boundary)
      XY yp = StepYPlus(x, y, Nx, Ny, shift);
      int sy = SiteIndex(yp.x, yp.y, Nx);
      AddUniqueLink(lattice, s, sy, k2);

      // Direction 3: (-x, +y) = (+y) then (-x)
      int sd = SiteIndex(Mod(yp.x - 1, Nx), yp.y, Nx);
      AddUniqueLink(lattice, s, sd, k3);
    }
  }
}

double tri_crit(double k1, double k2, double k3, double beta) {
  double p1 = std::exp(-2.0 * beta * (k2 + k3));
  double p2 = std::exp(-2.0 * beta * (k3 + k1));
  double p3 = std::exp(-2.0 * beta * (k1 + k2));

  double r1 = p1 + p2 + p3 - 1.0;
  double r2 = -2.0 * (p1 * (k2 + k3) + p2 * (k3 + k1) + p3 * (k1 + k2));
  return r1 / r2;
}

double find_crit(double k1, double k2, double k3) {
  double k_mean = (k1 + k2 + k3) / 3.0;
  k1 /= k_mean;
  k2 /= k_mean;
  k3 /= k_mean;

  double beta = 0.267949192431123;  // 2 - sqrt(3)
  for (int i = 0; i < 100; i++) beta -= tri_crit(k1, k2, k3, beta);
  return beta / k_mean;
}

double MeanFromSamples(const std::vector<double>& v) {
  if (v.empty()) return 0.0;
  double s = 0.0;
  for (double x : v) s += x;
  return s / double(v.size());
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

int OrbitPeriod(const std::vector<int>& next) {
  int V = int(next.size());
  int s0 = 0;
  int p = 0;
  int t = s0;
  do {
    t = next[t];
    p++;
  } while (t != s0 && p <= V);
  if (p > V) {
    std::fprintf(stderr, "ERROR: failed to find direction period (bound %d)\n",
                 V);
    std::exit(1);
  }
  return p;
}

}  // namespace

int main(int argc, char* argv[]) {
  int Nx = 32;
  int Ny = 32;
  int shift = 16;

  double k1 = 1.0;
  double k2 = 1.0;
  double k3 = 1.0;
  double beta_mult = 1.0;
  double beta_override = -1.0;

  int n_therm = 2000;
  int n_traj = 10000;
  int n_skip = 10;
  int n_wolff = 3;
  int n_metropolis = 5;
  int seed = 12345;

  std::string out_path;

  const struct option long_options[] = {
      {"n_x", required_argument, 0, 'X'},
      {"n_y", required_argument, 0, 'Y'},
      {"shift", required_argument, 0, 'S'},
      {"k1", required_argument, 0, 'a'},
      {"k2", required_argument, 0, 'b'},
      {"k3", required_argument, 0, 'c'},
      {"beta_mult", required_argument, 0, 'm'},
      {"beta", required_argument, 0, 'B'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'e'},
      {"seed", required_argument, 0, 'r'},
      {"out", required_argument, 0, 'o'},
      {0, 0, 0, 0}};

  const char* short_options = "X:Y:S:a:b:c:m:B:h:t:s:w:e:r:o:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'X':
        Nx = std::atoi(optarg);
        break;
      case 'Y':
        Ny = std::atoi(optarg);
        break;
      case 'S':
        shift = std::atoi(optarg);
        break;
      case 'a':
        k1 = std::atof(optarg);
        break;
      case 'b':
        k2 = std::atof(optarg);
        break;
      case 'c':
        k3 = std::atof(optarg);
        break;
      case 'm':
        beta_mult = std::atof(optarg);
        break;
      case 'B':
        beta_override = std::atof(optarg);
        break;
      case 'h':
        n_therm = std::atoi(optarg);
        break;
      case 't':
        n_traj = std::atoi(optarg);
        break;
      case 's':
        n_skip = std::atoi(optarg);
        break;
      case 'w':
        n_wolff = std::atoi(optarg);
        break;
      case 'e':
        n_metropolis = std::atoi(optarg);
        break;
      case 'r':
        seed = std::atoi(optarg);
        break;
      case 'o':
        out_path = optarg;
        break;
      default:
        break;
    }
  }

  if (Nx <= 0 || Ny <= 0) {
    std::fprintf(stderr, "ERROR: Nx and Ny must be positive.\n");
    return 1;
  }
  shift = Mod(shift, Nx);

  QfeLattice lattice;
  InitTwistedTriangle(&lattice, Nx, Ny, shift, k1, k2, k3);
  lattice.SeedRng(seed);

  double beta_crit = find_crit(k1, k2, k3);
  double beta = (beta_override > 0.0) ? beta_override : (beta_crit * beta_mult);

  QfeIsing field(&lattice, beta);
  field.HotStart();

  int V = Nx * Ny;
  std::vector<int> next_x(V, 0);
  std::vector<int> next_left(V, 0);
  std::vector<int> next_diag(V, 0);

  for (int y = 0; y < Ny; y++) {
    for (int x = 0; x < Nx; x++) {
      int s = SiteIndex(x, y, Nx);

      int x1 = Mod(x + 1, Nx);
      next_x[s] = SiteIndex(x1, y, Nx);

      XY yp = StepYPlus(x, y, Nx, Ny, shift);
      next_left[s] = SiteIndex(yp.x, yp.y, Nx);
      next_diag[s] = SiteIndex(Mod(yp.x - 1, Nx), yp.y, Nx);
    }
  }

  int p_x = OrbitPeriod(next_x);
  int p_left = OrbitPeriod(next_left);
  int p_diag = OrbitPeriod(next_diag);

  std::printf("Nx Ny shift: %d %d %d\n", Nx, Ny, shift);
  std::printf("k1 k2 k3: %.12f %.12f %.12f\n", k1, k2, k3);
  std::printf("beta_crit: %.12f\n", beta_crit);
  std::printf("beta: %.12f\n", beta);
  std::printf("periods x/left/diag: %d %d %d\n", p_x, p_left, p_diag);

  std::vector<std::vector<double>> corr_x_samples(p_x);
  std::vector<std::vector<double>> corr_left_samples(p_left);
  std::vector<std::vector<double>> corr_diag_samples(p_diag);
  std::vector<double> mag;

  QfeMeasReal cluster_size;
  QfeMeasReal accept_metropolis;

  for (int n = 0; n < (n_traj + n_therm); n++) {
    int cluster_size_sum = 0;
    for (int j = 0; j < n_wolff; j++) cluster_size_sum += field.WolffUpdate();

    double metropolis_sum = 0.0;
    for (int j = 0; j < n_metropolis; j++) metropolis_sum += field.Metropolis();

    cluster_size.Measure(double(cluster_size_sum) / double(V));
    accept_metropolis.Measure(metropolis_sum / double(n_metropolis));

    if (n % n_skip || n < n_therm) continue;

    std::vector<double> corr_x_sum(p_x, 0.0);
    std::vector<double> corr_left_sum(p_left, 0.0);
    std::vector<double> corr_diag_sum(p_diag, 0.0);

    for (int s = 0; s < V; s++) {
      double spin_s = field.spin[s];

      int t = s;
      for (int r = 0; r < p_x; r++) {
        corr_x_sum[r] += spin_s * field.spin[t];
        t = next_x[t];
      }

      t = s;
      for (int r = 0; r < p_left; r++) {
        corr_left_sum[r] += spin_s * field.spin[t];
        t = next_left[t];
      }

      t = s;
      for (int r = 0; r < p_diag; r++) {
        corr_diag_sum[r] += spin_s * field.spin[t];
        t = next_diag[t];
      }
    }

    double norm = double(V);
    for (int r = 0; r < p_x; r++) corr_x_samples[r].push_back(corr_x_sum[r] / norm);
    for (int r = 0; r < p_left; r++)
      corr_left_samples[r].push_back(corr_left_sum[r] / norm);
    for (int r = 0; r < p_diag; r++)
      corr_diag_samples[r].push_back(corr_diag_sum[r] / norm);

    mag.push_back(field.MeanSpin());
  }

  std::printf("accept_metropolis: %.4f\n", accept_metropolis.Mean());
  std::printf("cluster_size/V: %.4f\n", cluster_size.Mean());
  std::printf("m: %.12e (%.12e), %.4f\n", Mean(mag), JackknifeMean(mag),
              AutocorrTime(mag));

  mkdir("K_from_BC/results", 0755);
  if (out_path.empty()) {
    char tmp[256];
    std::snprintf(tmp, sizeof(tmp),
                  "K_from_BC/results/ising_twist_two_point_Nx%d_Ny%d_shift%d.dat",
                  Nx, Ny, shift);
    out_path = tmp;
  }

  FILE* file = std::fopen(out_path.c_str(), "w");
  if (!file) {
    std::fprintf(stderr, "ERROR: failed to open output file %s\n",
                 out_path.c_str());
    return 1;
  }

  std::fprintf(file,
               "# Nx %d Ny %d shift %d k1 %.12f k2 %.12f k3 %.12f beta %.12f\n",
               Nx, Ny, shift, k1, k2, k3, beta);
  std::fprintf(file, "# periods x %d left %d diag %d\n", p_x, p_left, p_diag);
  std::fprintf(file, "# dir r mean err\n");

  // dir=0: x
  for (int r = 0; r < p_x; r++) {
    double m = MeanFromSamples(corr_x_samples[r]);
    double e = JackknifeErrorMean(corr_x_samples[r]);
    std::printf("0 %04d %.12e %.12e\n", r, m, e);
    std::fprintf(file, "0 %04d %.12e %.12e\n", r, m, e);
  }

  // dir=1: left (+y with twist)
  for (int r = 0; r < p_left; r++) {
    double m = MeanFromSamples(corr_left_samples[r]);
    double e = JackknifeErrorMean(corr_left_samples[r]);
    std::printf("1 %04d %.12e %.12e\n", r, m, e);
    std::fprintf(file, "1 %04d %.12e %.12e\n", r, m, e);
  }

  // dir=2: diag (-x,+y)
  for (int r = 0; r < p_diag; r++) {
    double m = MeanFromSamples(corr_diag_samples[r]);
    double e = JackknifeErrorMean(corr_diag_samples[r]);
    std::printf("2 %04d %.12e %.12e\n", r, m, e);
    std::fprintf(file, "2 %04d %.12e %.12e\n", r, m, e);
  }

  // Connected correlators via jackknife with matching per-measurement magnetization.
  for (int r = 0; r < p_x; r++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_x_samples[r], mag, &m_conn, &e_conn);
    std::printf("3 %04d %.12e %.12e\n", r, m_conn, e_conn);
    std::fprintf(file, "3 %04d %.12e %.12e\n", r, m_conn, e_conn);
  }

  for (int r = 0; r < p_left; r++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_left_samples[r], mag, &m_conn, &e_conn);
    std::printf("4 %04d %.12e %.12e\n", r, m_conn, e_conn);
    std::fprintf(file, "4 %04d %.12e %.12e\n", r, m_conn, e_conn);
  }

  for (int r = 0; r < p_diag; r++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_diag_samples[r], mag, &m_conn, &e_conn);
    std::printf("5 %04d %.12e %.12e\n", r, m_conn, e_conn);
    std::fprintf(file, "5 %04d %.12e %.12e\n", r, m_conn, e_conn);
  }

  std::fclose(file);
  std::printf("wrote %s\n", out_path.c_str());
  return 0;
}
