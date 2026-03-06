// ising_flat_crit_k_from_bc.cc
//
// Adapted from examples/ising_flat_crit.cc.
// This variant lets the user specify a target triangle shape by angles
// (theta1, theta2) and builds an effective skewed-torus embedding for the
// periodic boundary conditions of an equilateral triangular lattice.

#include <getopt.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include "ising.h"
#include "statistics.h"

namespace {

struct SkewShape {
  double theta1;
  double theta2;
  double theta3;
  double skewness;   // x-offset over base for the triangle apex
  double aspect;     // height over base
  double e1x, e1y;   // primitive vector #1 (equilateral lattice direction)
  double e2x, e2y;   // primitive vector #2 (skewed by target triangle)
  double L1x, L1y;   // torus period vector #1
  double L2x, L2y;   // torus period vector #2
};

SkewShape BuildSkewShapeFromAngles(double theta1, double theta2,
                                   int Nx, int Ny) {
  const double pi = acos(-1.0);
  assert(theta1 > 0.0 && theta2 > 0.0);
  assert(theta1 + theta2 < pi);

  SkewShape shape;
  shape.theta1 = theta1;
  shape.theta2 = theta2;
  shape.theta3 = pi - theta1 - theta2;

  // Triangle with base AB = 1 and apex C determined by (theta1, theta2).
  // AC = sin(theta2)/sin(theta3), C = (AC cos theta1, AC sin theta1).
  double side_ac = sin(theta2) / sin(shape.theta3);
  shape.skewness = side_ac * cos(theta1);
  shape.aspect = side_ac * sin(theta1);

  // Local helper "augmentation" for skewed periodic geometry:
  // the graph is still periodic InitTriangle(Nx, Ny), but we interpret lattice
  // displacements in this oblique basis, so the periodic cell is a parallelogram.
  shape.e1x = 1.0;
  shape.e1y = 0.0;
  shape.e2x = shape.skewness;
  shape.e2y = shape.aspect;

  shape.L1x = Nx * shape.e1x;
  shape.L1y = Nx * shape.e1y;
  shape.L2x = Ny * shape.e2x;
  shape.L2y = Ny * shape.e2y;

  return shape;
}

double MinImageDistanceSkewTorus(int dx, int dy, const SkewShape& shape) {
  // Convert lattice displacement to physical displacement in oblique basis.
  double rx = dx * shape.e1x + dy * shape.e2x;
  double ry = dx * shape.e1y + dy * shape.e2y;

  // Minimum image using neighboring torus copies.
  double best2 = 1.0e300;
  for (int n1 = -1; n1 <= 1; n1++) {
    for (int n2 = -1; n2 <= 1; n2++) {
      double x = rx + n1 * shape.L1x + n2 * shape.L2x;
      double y = ry + n1 * shape.L1y + n2 * shape.L2y;
      double d2 = x * x + y * y;
      if (d2 < best2) best2 = d2;
    }
  }
  return sqrt(best2);
}

double tri_crit(double k1, double k2, double k3, double beta) {
  double p1 = exp(-2.0 * beta * (k2 + k3));
  double p2 = exp(-2.0 * beta * (k3 + k1));
  double p3 = exp(-2.0 * beta * (k1 + k2));

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

int LcmInt(int a, int b) {
  return (a / GcdInt(a, b)) * b;
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
  int Nx = 32;
  int Ny = 32;

  double theta1 = 1.0471975512;  // 60 deg default (equilateral)
  double theta2 = 1.0471975512;  // 60 deg default (equilateral)

  double k1 = 1.0;
  double k2 = 1.0;
  double k3 = 1.0;
  double beta_mult = 1.0;

  int n_therm = 2000;
  int n_traj = 50000;
  int n_skip = 20;
  int n_wolff = 3;
  int n_metropolis = 5;

  const struct option long_options[] = {
    {"n_x", required_argument, 0, 'X'},
    {"n_y", required_argument, 0, 'Y'},
    {"theta1", required_argument, 0, 'p'},
    {"theta2", required_argument, 0, 'q'},
    {"k1", required_argument, 0, 'a'},
    {"k2", required_argument, 0, 'b'},
    {"k3", required_argument, 0, 'c'},
    {"beta_mult", required_argument, 0, 'm'},
    {"n_therm", required_argument, 0, 'h'},
    {"n_traj", required_argument, 0, 't'},
    {"n_skip", required_argument, 0, 's'},
    {"n_wolff", required_argument, 0, 'w'},
    {"n_metropolis", required_argument, 0, 'e'},
    {0, 0, 0, 0}};

  const char* short_options = "X:Y:p:q:a:b:c:m:h:t:s:w:e:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'X': Nx = atoi(optarg); break;
      case 'Y': Ny = atoi(optarg); break;
      case 'p': theta1 = std::stod(optarg); break;
      case 'q': theta2 = std::stod(optarg); break;
      case 'a': k1 = std::stod(optarg); break;
      case 'b': k2 = std::stod(optarg); break;
      case 'c': k3 = std::stod(optarg); break;
      case 'm': beta_mult = std::stod(optarg); break;
      case 'h': n_therm = atoi(optarg); break;
      case 't': n_traj = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'w': n_wolff = atoi(optarg); break;
      case 'e': n_metropolis = atoi(optarg); break;
      default: break;
    }
  }

  SkewShape shape = BuildSkewShapeFromAngles(theta1, theta2, Nx, Ny);

  printf("Nx Ny: %d %d\n", Nx, Ny);
  printf("theta1 theta2 theta3: %.12f %.12f %.12f\n", shape.theta1, shape.theta2,
         shape.theta3);
  printf("target skewness/aspect: %.12f %.12f\n", shape.skewness, shape.aspect);
  printf("period vectors L1=(%.6f, %.6f), L2=(%.6f, %.6f)\n", shape.L1x,
         shape.L1y, shape.L2x, shape.L2y);

  printf("k1 k2 k3: %.12f %.12f %.12f\n", k1, k2, k3);
  printf("beta_mult: %.12f\n", beta_mult);

  QfeLattice lattice;
  lattice.InitTriangle(Nx, Ny, k1, k2, k3);

  QfeIsing field(&lattice, find_crit(k1, k2, k3) * beta_mult);
  field.HotStart();

  printf("beta: %.12f\n", field.beta);
  printf("initial action: %.12f\n", field.Action());

  // All-to-all directional two-point correlators over full periodic lengths.
  // corr_x[r]      = < s(x,y) s(x+r,y) >, r = 0..Nx-1
  // corr_y[r]      = < s(x,y) s(x,y+r) >, r = 0..Ny-1
  // corr_xpy[r]    = < s(x,y) s(x+r,y+r) >, r = 0..Nx-1
  // corr_ymyx[r]   = < s(x,y) s(x+r,y-r) >, r = 0..Nx-1
  int r_bins_x = Nx;
  int r_bins_y = Ny;
  // Oblique (x+y and y-x) periods are lcm(Nx, Ny) on a rectangular torus.
  int r_bins_oblique = LcmInt(Nx, Ny);

  std::vector<QfeMeasReal> corr_x(r_bins_x);
  std::vector<QfeMeasReal> corr_y(r_bins_y);
  std::vector<QfeMeasReal> corr_xpy(r_bins_oblique);
  std::vector<QfeMeasReal> corr_ymyx(r_bins_oblique);
  std::vector<std::vector<double>> corr_x_samples(r_bins_x);
  std::vector<std::vector<double>> corr_y_samples(r_bins_y);
  std::vector<std::vector<double>> corr_xpy_samples(r_bins_oblique);
  std::vector<std::vector<double>> corr_ymyx_samples(r_bins_oblique);

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

    std::vector<double> corr_x_sum(r_bins_x, 0.0);
    std::vector<double> corr_y_sum(r_bins_y, 0.0);
    std::vector<double> corr_xpy_sum(r_bins_oblique, 0.0);
    std::vector<double> corr_ymyx_sum(r_bins_oblique, 0.0);

    // Direct all-to-all estimator along x and y directions.
    for (int y = 0; y < Ny; y++) {
      for (int x = 0; x < Nx; x++) {
        int s = y * Nx + x;
        double spin_s = field.spin[s];

        for (int r = 0; r < r_bins_x; r++) {
          int xp = (x + r) % Nx;
          int sp = y * Nx + xp;
          corr_x_sum[r] += spin_s * field.spin[sp];
        }

        for (int r = 0; r < r_bins_y; r++) {
          int yp = (y + r) % Ny;
          int sp = yp * Nx + x;
          corr_y_sum[r] += spin_s * field.spin[sp];
        }

        for (int r = 0; r < r_bins_oblique; r++) {
          int xp = (x + r) % Nx;
          int yp_plus = (y + r) % Ny;
          int yp_minus = (y - (r % Ny) + Ny) % Ny;

          int sp_plus = yp_plus * Nx + xp;
          int sp_minus = yp_minus * Nx + xp;

          corr_xpy_sum[r] += spin_s * field.spin[sp_plus];
          corr_ymyx_sum[r] += spin_s * field.spin[sp_minus];
        }
      }
    }

    double norm = double(Nx * Ny);
    for (int i = 0; i < r_bins_x; i++) {
      double v = corr_x_sum[i] / norm;
      corr_x[i].Measure(v);
      corr_x_samples[i].push_back(v);
    }

    for (int i = 0; i < r_bins_y; i++) {
      double v = corr_y_sum[i] / norm;
      corr_y[i].Measure(v);
      corr_y_samples[i].push_back(v);
    }

    for (int i = 0; i < r_bins_oblique; i++) {
      double v = corr_xpy_sum[i] / norm;
      corr_xpy[i].Measure(v);
      corr_xpy_samples[i].push_back(v);
    }

    for (int i = 0; i < r_bins_oblique; i++) {
      double v = corr_ymyx_sum[i] / norm;
      corr_ymyx[i].Measure(v);
      corr_ymyx_samples[i].push_back(v);
    }

    action.push_back(field.Action());
    mag.push_back(field.MeanSpin());
    printf("%06d %.12f %+.12f %.4f %.4f\n", n, action.back(), mag.back(),
           accept_metropolis.last, cluster_size.last);
  }

  std::vector<double> mag_abs(mag.size());
  std::vector<double> mag2(mag.size());
  std::vector<double> mag4(mag.size());
  for (int i = 0; i < int(mag.size()); i++) {
    double m = mag[i];
    double m2 = m * m;
    mag_abs[i] = fabs(m);
    mag2[i] = m2;
    mag4[i] = m2 * m2;
  }

  printf("accept_metropolis: %.4f\n", accept_metropolis.Mean());
  printf("cluster_size/V: %.4f\n", cluster_size.Mean());
  printf("action: %.12e (%.12e), %.4f\n", Mean(action), JackknifeMean(action),
         AutocorrTime(action));
  printf("m: %.12e (%.12e), %.4f\n", Mean(mag), JackknifeMean(mag),
         AutocorrTime(mag));
  printf("m^2: %.12e (%.12e), %.4f\n", Mean(mag2), JackknifeMean(mag2),
         AutocorrTime(mag2));
  printf("m^4: %.12e (%.12e), %.4f\n", Mean(mag4), JackknifeMean(mag4),
         AutocorrTime(mag4));
  printf("U4: %.12e (%.12e)\n", U4(mag2, mag4), JackknifeU4(mag2, mag4));
  printf("susceptibility: %.12e (%.12e)\n", Susceptibility(mag2, mag_abs),
         JackknifeSusceptibility(mag2, mag_abs));

  mkdir("ising_flat_crit_k_from_bc", 0755);

  char path[200];
  sprintf(path,
          "ising_flat_crit_k_from_bc/%d_%d_t1_%.3f_t2_%.3f_k_%.3f_%.3f_%.3f.dat",
          Nx, Ny, theta1, theta2, k1, k2, k3);
  FILE* file = fopen(path, "w");
  assert(file != nullptr);

  printf("\ncorr_x:\n");
  for (int i = 0; i < r_bins_x; i++) {
    double m = MeanFromSamples(corr_x_samples[i]);
    double e = JackknifeErrorMean(corr_x_samples[i]);
    printf("0 %04d %.12e %.12e\n", i, m, e);
    fprintf(file, "0 %04d %.12e %.12e\n", i, m, e);
  }

  printf("\ncorr_y:\n");
  for (int i = 0; i < r_bins_y; i++) {
    double m = MeanFromSamples(corr_y_samples[i]);
    double e = JackknifeErrorMean(corr_y_samples[i]);
    printf("1 %04d %.12e %.12e\n", i, m, e);
    fprintf(file, "1 %04d %.12e %.12e\n", i, m, e);
  }

  printf("\ncorr_x_plus_y:\n");
  for (int i = 0; i < r_bins_oblique; i++) {
    double m = MeanFromSamples(corr_xpy_samples[i]);
    double e = JackknifeErrorMean(corr_xpy_samples[i]);
    printf("2 %04d %.12e %.12e\n", i, m, e);
    fprintf(file, "2 %04d %.12e %.12e\n", i, m, e);
  }

  printf("\ncorr_y_minus_x:\n");
  for (int i = 0; i < r_bins_oblique; i++) {
    double m = MeanFromSamples(corr_ymyx_samples[i]);
    double e = JackknifeErrorMean(corr_ymyx_samples[i]);
    printf("3 %04d %.12e %.12e\n", i, m, e);
    fprintf(file, "3 %04d %.12e %.12e\n", i, m, e);
  }

  printf("\nconnected_corr_x (jackknife):\n");
  for (int i = 0; i < r_bins_x; i++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_x_samples[i], mag, &m_conn, &e_conn);
    printf("4 %04d %.12e %.12e\n", i, m_conn, e_conn);
    fprintf(file, "4 %04d %.12e %.12e\n", i, m_conn, e_conn);
  }

  printf("\nconnected_corr_y (jackknife):\n");
  for (int i = 0; i < r_bins_y; i++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_y_samples[i], mag, &m_conn, &e_conn);
    printf("5 %04d %.12e %.12e\n", i, m_conn, e_conn);
    fprintf(file, "5 %04d %.12e %.12e\n", i, m_conn, e_conn);
  }

  printf("\nconnected_corr_x_plus_y (jackknife):\n");
  for (int i = 0; i < r_bins_oblique; i++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_xpy_samples[i], mag, &m_conn, &e_conn);
    printf("6 %04d %.12e %.12e\n", i, m_conn, e_conn);
    fprintf(file, "6 %04d %.12e %.12e\n", i, m_conn, e_conn);
  }

  printf("\nconnected_corr_y_minus_x (jackknife):\n");
  for (int i = 0; i < r_bins_oblique; i++) {
    double m_conn = 0.0;
    double e_conn = 0.0;
    JackknifeConnectedCorr(corr_ymyx_samples[i], mag, &m_conn, &e_conn);
    printf("7 %04d %.12e %.12e\n", i, m_conn, e_conn);
    fprintf(file, "7 %04d %.12e %.12e\n", i, m_conn, e_conn);
  }

  int midpoint = r_bins_x / 2;
  double corr_x_mid = MeanFromSamples(corr_x_samples[midpoint]);
  double corr_y_mid = MeanFromSamples(corr_y_samples[std::min(midpoint, r_bins_y - 1)]);
  double corr_xpy_mid = MeanFromSamples(corr_xpy_samples[midpoint]);
  double corr_ymyx_mid = MeanFromSamples(corr_ymyx_samples[midpoint]);
  printf("\nmidpoint_bin: %d\n", midpoint);
  printf("corr_x(midpoint): %.12e\n", corr_x_mid);
  printf("corr_y(midpoint): %.12e\n", corr_y_mid);
  printf("corr_x_plus_y(midpoint): %.12e\n", corr_xpy_mid);
  printf("corr_y_minus_x(midpoint): %.12e\n", corr_ymyx_mid);
  if (std::abs(corr_x_mid) > 0.0) {
    printf("corr_y/corr_x(midpoint): %.12e\n", corr_y_mid / corr_x_mid);
  }

  fclose(file);
  return 0;
}