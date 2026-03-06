// tri_histogram_reweight.cc
//
// Simple data collection for histogram reweighting analysis.
// Runs equilateral triangular Ising model on square lattice (N x N),
// saves m and m^2 for each configuration to enable later histogram
// reweighting analysis to find the critical point.
//
// As beta is varied, we collect statistics on the magnetization that
// can later be reweighted to construct the free energy landscape.

#include <getopt.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include "ising.h"
#include "statistics.h"

int main(int argc, char* argv[]) {
  int N = 32;
  int n_therm = 1000;
  int n_traj = 5000;
  int n_skip = 10;
  int n_wolff = 3;
  int n_metropolis = 5;
  
  double beta_min = 0.20;
  double beta_max = 0.28;
  int n_beta = 5;

  const struct option long_options[] = {
    {"N", required_argument, 0, 'N'},
    {"n_therm", required_argument, 0, 'h'},
    {"n_traj", required_argument, 0, 't'},
    {"n_skip", required_argument, 0, 's'},
    {"n_wolff", required_argument, 0, 'w'},
    {"n_metropolis", required_argument, 0, 'e'},
    {"beta_min", required_argument, 0, 'a'},
    {"beta_max", required_argument, 0, 'b'},
    {"n_beta", required_argument, 0, 'n'},
    {0, 0, 0, 0}};

  const char* short_options = "N:h:t:s:w:e:a:b:n:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'N': N = atoi(optarg); break;
      case 'h': n_therm = atoi(optarg); break;
      case 't': n_traj = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'w': n_wolff = atoi(optarg); break;
      case 'e': n_metropolis = atoi(optarg); break;
      case 'a': beta_min = std::stod(optarg); break;
      case 'b': beta_max = std::stod(optarg); break;
      case 'n': n_beta = atoi(optarg); break;
      default: break;
    }
  }

  printf("# Equilateral triangular Ising lattice histogram reweighting data\n");
  printf("# N: %d\n", N);
  printf("# n_therm: %d, n_traj: %d, n_skip: %d\n", n_therm, n_traj, n_skip);
  printf("# n_wolff: %d, n_metropolis: %d\n", n_wolff, n_metropolis);
  printf("# beta_min: %.12f, beta_max: %.12f, n_beta: %d\n", beta_min, beta_max,
         n_beta);

  mkdir("K_from_BC/histogram/tri_histogram_reweight", 0755);

  // Loop over beta values
  for (int beta_idx = 0; beta_idx < n_beta; beta_idx++) {
    double beta;
    if (n_beta == 1) {
      beta = beta_min;
    } else {
      beta = beta_min + (beta_max - beta_min) * beta_idx / (n_beta - 1.0);
    }

    printf("\n# Processing beta = %.12f\n", beta);

    QfeLattice lattice;
    lattice.InitTriangle(N, 1.0, 1.0, 1.0);

    QfeIsing field(&lattice, beta);
    field.HotStart();

    char dat_filename[256];
    sprintf(dat_filename, "K_from_BC/histogram/tri_histogram_reweight/N%d_beta_%.4f.dat", N, beta);
    FILE* dat_file = fopen(dat_filename, "w");
    assert(dat_file != nullptr);

    // Write header
    fprintf(dat_file, "# beta %.12f\n", beta);
    fprintf(dat_file, "# N %d\n", N);
    fprintf(dat_file, "# m m2 action\n");

    QfeMeasReal cluster_size;
    QfeMeasReal accept_metropolis;

    for (int n = 0; n < (n_traj + n_therm); n++) {
      int cluster_size_sum = 0;
      for (int j = 0; j < n_wolff; j++) {
        cluster_size_sum += field.WolffUpdate();
      }
      double metropolis_sum = 0.0;
      for (int j = 0; j < n_metropolis; j++) {
        metropolis_sum += field.Metropolis();
      }
      cluster_size.Measure(double(cluster_size_sum) / lattice.vol);
      accept_metropolis.Measure(metropolis_sum);

      // Skip if in thermalization or not at measurement interval
      if (n < n_therm || n % n_skip != 0) continue;

      double m = fabs(field.MeanSpin());
      double m2 = m * m;
      double action = field.Action();

      fprintf(dat_file, "%.12e %.12e %.12e\n", m, m2, action);

      if ((n - n_therm) % (n_skip * 100) == 0) {
        printf("  traj %06d: m=%.6f m2=%.6f action=%.6f accept=%.4f "
               "cluster_frac=%.4f\n",
               n, m, m2, action, accept_metropolis.last, cluster_size.last);
      }
    }

    fclose(dat_file);
    printf("# Wrote %s\n", dat_filename);

    printf("# Final measurements for beta = %.12f:\n", beta);
    printf("# accept_metropolis: %.4f\n", accept_metropolis.Mean());
    printf("# cluster_size/V: %.4f\n", cluster_size.Mean());
  }

  printf("\n# Data collection complete.\n");
  printf(
      "# Files saved in K_from_BC/histogram/tri_histogram_reweight/ directory.\n"
      "# Format: m m^2 action\n"
      "# Ready for histogram reweighting analysis.\n");

  return 0;
}
