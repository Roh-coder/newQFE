// potts_sw4_scan.cc

#include <getopt.h>

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>

#include "lattice.h"

namespace {

constexpr int kQ = 4;

struct DSU {
  explicit DSU(int n = 0) { Reset(n); }

  void Reset(int n) {
    parent.resize(n);
    rank.assign(n, 0);
    for (int i = 0; i < n; ++i) parent[i] = i;
  }

  int Find(int x) {
    if (parent[x] != x) parent[x] = Find(parent[x]);
    return parent[x];
  }

  void Union(int a, int b) {
    int ra = Find(a);
    int rb = Find(b);
    if (ra == rb) return;
    if (rank[ra] < rank[rb]) {
      parent[ra] = rb;
    } else if (rank[ra] > rank[rb]) {
      parent[rb] = ra;
    } else {
      parent[rb] = ra;
      rank[ra]++;
    }
  }

  std::vector<int> parent;
  std::vector<int> rank;
};

struct BetaResult {
  double beta = 0.0;
  double mag_mean = 0.0;
  double mag_err = 0.0;
  double chi = 0.0;
  double chi_err = 0.0;
};

struct PeakResult {
  double beta_peak = 0.0;
  double beta_err = 0.0;
  double chi_peak = 0.0;
  double chi_err = 0.0;
  int peak_idx = 0;
};

bool MakeDirIfMissing(const char* path) {
  if (mkdir(path, 0775) == 0) return true;
  if (errno == EEXIST) return true;
  return false;
}

std::vector<int> ParseLList(const std::string& text) {
  std::vector<int> out;
  std::string clean = text;
  for (char& c : clean) {
    if (c == ';' || c == ':' || c == ' ') c = ',';
  }

  std::stringstream ss(clean);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item.empty()) continue;
    int value = std::atoi(item.c_str());
    if (value > 1) out.push_back(value);
  }

  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

double Mean(const std::vector<double>& a) {
  if (a.empty()) return 0.0;
  double sum = 0.0;
  for (double x : a) sum += x;
  return sum / double(a.size());
}

double MeanError(const std::vector<double>& a) {
  if (a.size() < 2) return 0.0;
  double mean = Mean(a);
  double sum2 = 0.0;
  for (double x : a) {
    double d = x - mean;
    sum2 += d * d;
  }
  return std::sqrt(sum2 / (double(a.size()) * double(a.size() - 1)));
}

void SwendsenWangSweep(QfeLattice& lattice, std::vector<int>& spin, double beta,
                       DSU& dsu) {
  const int n_sites = lattice.n_sites;
  dsu.Reset(n_sites);

  for (int l = 0; l < lattice.n_links; ++l) {
    const QfeLink& link = lattice.links[l];
    const int s0 = link.sites[0];
    const int s1 = link.sites[1];

    if (spin[s0] != spin[s1]) continue;

    const double p = 1.0 - std::exp(-beta * link.wt);
    if (lattice.rng.RandReal() < p) {
      dsu.Union(s0, s1);
    }
  }

  std::unordered_map<int, int> cluster_spin;
  cluster_spin.reserve(n_sites);

  for (int s = 0; s < n_sites; ++s) {
    int root = dsu.Find(s);
    auto it = cluster_spin.find(root);
    if (it == cluster_spin.end()) {
      int new_spin = lattice.rng.RandInt(0, kQ - 1);
      cluster_spin[root] = new_spin;
      spin[s] = new_spin;
    } else {
      spin[s] = it->second;
    }
  }
}

double PottsMagnetization(const std::vector<int>& spin) {
  std::vector<int> count(kQ, 0);
  for (int s : spin) count[s]++;

  int n_sites = int(spin.size());
  int n_max = *std::max_element(count.begin(), count.end());
  return (double(kQ * n_max) / double(n_sites) - 1.0) / double(kQ - 1);
}

BetaResult RunAtBeta(QfeLattice& lattice, int n_therm, int n_meas, int n_skip,
                     double beta) {
  const int n_sites = lattice.n_sites;
  std::vector<int> spin(n_sites, 0);
  for (int s = 0; s < n_sites; ++s) {
    spin[s] = lattice.rng.RandInt(0, kQ - 1);
  }

  DSU dsu(n_sites);

  for (int n = 0; n < n_therm; ++n) {
    SwendsenWangSweep(lattice, spin, beta, dsu);
  }

  std::vector<double> mag_samples;
  mag_samples.reserve(n_meas);

  int sweeps_between = std::max(1, n_skip);
  for (int n = 0; n < n_meas; ++n) {
    for (int k = 0; k < sweeps_between; ++k) {
      SwendsenWangSweep(lattice, spin, beta, dsu);
    }
    mag_samples.push_back(PottsMagnetization(spin));
  }

  BetaResult res;
  res.beta = beta;
  res.mag_mean = Mean(mag_samples);
  res.mag_err = MeanError(mag_samples);

  if (mag_samples.size() < 2) {
    res.chi = 0.0;
    res.chi_err = 0.0;
    return res;
  }

  double sum_m = 0.0;
  double sum_m2 = 0.0;
  for (double m : mag_samples) {
    sum_m += m;
    sum_m2 += m * m;
  }

  const int n = int(mag_samples.size());
  const double mean_m = sum_m / double(n);
  const double mean_m2 = sum_m2 / double(n);
  res.chi = double(n_sites) * (mean_m2 - mean_m * mean_m);

  std::vector<double> chi_jk;
  chi_jk.resize(n);
  for (int i = 0; i < n; ++i) {
    const double m = mag_samples[i];
    const double m2 = m * m;
    const double m_excl = (sum_m - m) / double(n - 1);
    const double m2_excl = (sum_m2 - m2) / double(n - 1);
    chi_jk[i] = double(n_sites) * (m2_excl - m_excl * m_excl);
  }

  const double chi_jk_mean = Mean(chi_jk);
  double chi_var = 0.0;
  for (double x : chi_jk) {
    const double d = x - chi_jk_mean;
    chi_var += d * d;
  }
  res.chi_err = std::sqrt((double(n - 1) / double(n)) * chi_var);

  return res;
}

PeakResult FindSusceptibilityPeak(const std::vector<BetaResult>& data) {
  assert(!data.empty());
  int idx = 0;
  for (int i = 1; i < int(data.size()); ++i) {
    if (data[i].chi > data[idx].chi) idx = i;
  }

  PeakResult peak;
  peak.peak_idx = idx;
  peak.beta_peak = data[idx].beta;
  peak.chi_peak = data[idx].chi;
  peak.chi_err = data[idx].chi_err;

  if (data.size() == 1) {
    peak.beta_err = 0.0;
    return peak;
  }

  double dbeta = std::fabs(data[1].beta - data[0].beta);
  peak.beta_err = 0.5 * dbeta;

  if (idx > 0 && idx + 1 < int(data.size())) {
    const double x1 = data[idx].beta;
    const double x2 = data[idx + 1].beta;
    const double y0 = data[idx - 1].chi;
    const double y1 = data[idx].chi;
    const double y2 = data[idx + 1].chi;

    const double denom = (y0 - 2.0 * y1 + y2);
    if (std::fabs(denom) > 1.0e-14) {
      const double dx = x2 - x1;
      const double shift = 0.5 * (y0 - y2) / denom * dx;
      peak.beta_peak = x1 + shift;

      const double a = 0.5 * denom / (dx * dx);
      const double b = (y2 - y0) / (2.0 * dx);
      const double x_shift = shift;
      peak.chi_peak = y1 + b * x_shift + a * x_shift * x_shift;
    }
  }

  return peak;
}

void PrintUsage(const char* prog) {
  std::printf(
      "Usage: %s [options]\n"
      "Options:\n"
      "  --L_list <csv>      L values (default: 16,24,32,48)\n"
      "  --beta_min <x>      minimum beta (default: 0.95)\n"
      "  --beta_max <x>      maximum beta (default: 1.20)\n"
      "  --n_beta <n>        number of beta points (default: 20)\n"
      "  --n_therm <n>       thermalization SW sweeps (default: 1000)\n"
      "  --n_meas <n>        measured samples per beta (default: 4000)\n"
      "  --n_skip <n>        SW sweeps between samples (default: 1)\n"
      "  --seed <n>          base RNG seed (default: 123456)\n"
      "  --out_dir <path>    output directory (default: Potts/output)\n"
      "\n"
      "For square-lattice q=4 Potts, beta_c(infinite volume)=ln(3)~1.098612.\n",
      prog);
}

}  // namespace

int main(int argc, char* argv[]) {
  std::string l_list_str = "16,24,32,48";
  double beta_min = 0.95;
  double beta_max = 1.20;
  int n_beta = 20;
  int n_therm = 1000;
  int n_meas = 4000;
  int n_skip = 1;
  int seed = 123456;
  std::string out_dir = "Potts/output";

  const struct option long_options[] = {
      {"L_list", required_argument, 0, 'L'},
      {"beta_min", required_argument, 0, 'a'},
      {"beta_max", required_argument, 0, 'b'},
      {"n_beta", required_argument, 0, 'B'},
      {"n_therm", required_argument, 0, 't'},
      {"n_meas", required_argument, 0, 'm'},
      {"n_skip", required_argument, 0, 's'},
      {"seed", required_argument, 0, 'r'},
      {"out_dir", required_argument, 0, 'o'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};

  while (true) {
    int opt_idx = 0;
    int c = getopt_long(argc, argv, "L:a:b:B:t:m:s:r:o:h", long_options,
                        &opt_idx);
    if (c == -1) break;

    switch (c) {
      case 'L': l_list_str = optarg; break;
      case 'a': beta_min = std::atof(optarg); break;
      case 'b': beta_max = std::atof(optarg); break;
      case 'B': n_beta = std::atoi(optarg); break;
      case 't': n_therm = std::atoi(optarg); break;
      case 'm': n_meas = std::atoi(optarg); break;
      case 's': n_skip = std::atoi(optarg); break;
      case 'r': seed = std::atoi(optarg); break;
      case 'o': out_dir = optarg; break;
      case 'h':
        PrintUsage(argv[0]);
        return 0;
      default:
        PrintUsage(argv[0]);
        return 1;
    }
  }

  if (beta_max <= beta_min || n_beta < 2 || n_therm < 0 || n_meas < 2 ||
      n_skip < 1) {
    std::fprintf(stderr,
                 "Invalid parameters: require beta_max>beta_min, n_beta>=2, "
                 "n_meas>=2, n_skip>=1.\n");
    return 1;
  }

  std::vector<int> L_list = ParseLList(l_list_str);
  if (L_list.empty()) {
    std::fprintf(stderr, "No valid L values parsed from --L_list.\n");
    return 1;
  }

  if (!MakeDirIfMissing("Potts")) {
    std::fprintf(stderr, "Warning: could not ensure directory Potts exists.\n");
  }
  if (!MakeDirIfMissing(out_dir.c_str())) {
    std::fprintf(stderr, "Warning: could not create output directory %s\n",
                 out_dir.c_str());
  }

  std::vector<double> beta_grid(n_beta, 0.0);
  for (int i = 0; i < n_beta; ++i) {
    beta_grid[i] = beta_min + (beta_max - beta_min) * double(i) / double(n_beta - 1);
  }

  char summary_path[512];
  std::snprintf(summary_path, sizeof(summary_path), "%s/critical_points.dat",
                out_dir.c_str());
  FILE* summary = std::fopen(summary_path, "w");
  if (summary == nullptr) {
    std::fprintf(stderr, "Failed to open %s for writing.\n", summary_path);
    return 1;
  }

  std::fprintf(summary,
               "# q=4 Potts SW pseudocritical points from susceptibility peaks\n");
  std::fprintf(summary,
               "# L beta_peak beta_err chi_peak chi_err peak_grid_index\n");

  std::printf("q = %d SW Potts scan on square lattices\n", kQ);
  std::printf("beta range: [%.6f, %.6f], n_beta = %d\n", beta_min, beta_max,
              n_beta);
  std::printf("n_therm = %d, n_meas = %d, n_skip = %d\n", n_therm, n_meas,
              n_skip);

  for (int li = 0; li < int(L_list.size()); ++li) {
    int L = L_list[li];
    QfeLattice lattice;
    lattice.InitRect(L, L, 1.0, 1.0);

    std::vector<BetaResult> series;
    series.reserve(beta_grid.size());

    std::printf("\n--- L = %d (V = %d) ---\n", L, lattice.n_sites);

    for (int bi = 0; bi < int(beta_grid.size()); ++bi) {
      double beta = beta_grid[bi];
      int beta_seed = seed + 100000 * li + bi;
      lattice.SeedRng(beta_seed);

      BetaResult r = RunAtBeta(lattice, n_therm, n_meas, n_skip, beta);
      series.push_back(r);

      std::printf("L=%d beta=%.8f m=%.6f +/- %.6f chi=%.6f +/- %.6f\n", L,
                  r.beta, r.mag_mean, r.mag_err, r.chi, r.chi_err);
    }

    PeakResult peak = FindSusceptibilityPeak(series);

    char scan_path[512];
    std::snprintf(scan_path, sizeof(scan_path), "%s/susceptibility_L%d.dat",
                  out_dir.c_str(), L);
    FILE* scan = std::fopen(scan_path, "w");
    if (scan == nullptr) {
      std::fprintf(stderr, "Failed to open %s for writing.\n", scan_path);
      std::fclose(summary);
      return 1;
    }

    std::fprintf(scan,
                 "# L=%d q=4 Potts SW scan\n"
                 "# beta m_mean m_err susceptibility chi_err\n",
                 L);
    for (const BetaResult& r : series) {
      std::fprintf(scan, "%.12f %.12e %.12e %.12e %.12e\n", r.beta, r.mag_mean,
                   r.mag_err, r.chi, r.chi_err);
    }
    std::fclose(scan);

    std::fprintf(summary, "%d %.12f %.12f %.12e %.12e %d\n", L,
                 peak.beta_peak, peak.beta_err, peak.chi_peak, peak.chi_err,
                 peak.peak_idx);

    std::printf("Peak(L=%d): beta_c(L)=%.8f +/- %.8f, chi_peak=%.6f\n", L,
                peak.beta_peak, peak.beta_err, peak.chi_peak);
    std::printf("Wrote: %s\n", scan_path);
  }

  std::fclose(summary);
  std::printf("\nWrote pseudocritical summary: %s\n", summary_path);

  return 0;
}
