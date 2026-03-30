#include <getopt.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#include "ising.h"
#include "statistics.h"
#include "timer.h"

namespace {

struct Coord {
  int m;
  int n;
};

struct TwistedMap {
  std::vector<int> next_e1;
  std::vector<int> next_e2;
  std::vector<int> next_e2me1;
  std::vector<int> disp_m;
  std::vector<int> disp_n;
  std::vector<int> add_table;
};

int Mod(int a, int n) {
  int r = a % n;
  return (r < 0) ? (r + n) : r;
}

int64_t MakeKey(int k1, int k2) {
  return (int64_t(uint32_t(k1)) << 32) | uint32_t(k2);
}

int64_t CosetKey(int m, int n, int Lx, int Ly, int Tx, int Ty, int Ncell) {
  // Invariants for Z^2 / <(Lx,Ty), (Tx,-Ly)>.
  int k1 = Mod(-Ly * m - Tx * n, Ncell);
  int k2 = Mod(-Ty * m + Lx * n, Ncell);
  return MakeKey(k1, k2);
}

bool MakeDir(const std::string& path) {
  if (path.empty()) return true;
  if (::mkdir(path.c_str(), 0755) == 0) return true;
  return errno == EEXIST;
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

void InitTriangleTwistedParallelogram(QfeLattice* lattice,
                                      int Lx, int Ly, int Tx, int Ty,
                                      double K1, double K2, double K3,
                                      TwistedMap* map_out) {
  const int Ncell = Lx * Ly + Tx * Ty;
  assert(Ncell > 0);

  // Representatives are the integer points in one half-open fundamental
  // parallelogram: p = a v + b u, with a,b in [0,1),
  // v=(Lx,Ty), u=(Tx,-Ly).
  std::unordered_map<int64_t, int> key_to_id;
  std::vector<Coord> reps;
  reps.reserve(Ncell);

  const int m0 = 0;
  const int n0 = 0;
  const int m1 = Lx;
  const int n1 = Ty;
  const int m2 = Tx;
  const int n2 = -Ly;
  const int m3 = Lx + Tx;
  const int n3 = Ty - Ly;

  int m_min = std::min(std::min(m0, m1), std::min(m2, m3));
  int m_max = std::max(std::max(m0, m1), std::max(m2, m3));
  int n_min = std::min(std::min(n0, n1), std::min(n2, n3));
  int n_max = std::max(std::max(n0, n1), std::max(n2, n3));

  const double invN = 1.0 / double(Ncell);
  const double eps = 1.0e-12;

  for (int m = m_min; m <= m_max; m++) {
    for (int n = n_min; n <= n_max; n++) {
      // Solve p = a v + b u using inverse of [[Lx,Tx],[Ty,-Ly]].
      double a = (double(Ly) * double(m) + double(Tx) * double(n)) * invN;
      double b = (double(Ty) * double(m) - double(Lx) * double(n)) * invN;

      if (a < -eps || a >= 1.0 - eps) continue;
      if (b < -eps || b >= 1.0 - eps) continue;

      int64_t key = CosetKey(m, n, Lx, Ly, Tx, Ty, Ncell);
      if (key_to_id.find(key) != key_to_id.end()) continue;

      int id = int(reps.size());
      key_to_id[key] = id;
      reps.push_back({m, n});
    }
  }

  assert(int(reps.size()) == Ncell);

  lattice->ResizeSites(Ncell);
  lattice->vol = double(Ncell);
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

  // Triangular directions: +e1, +e2, +(e2-e1).
  const int dm[3] = {1, 0, -1};
  const int dn[3] = {0, 1, 1};
  const double wt[3] = {K1, K2, K3};

  if (map_out != nullptr) {
    map_out->next_e1.assign(Ncell, -1);
    map_out->next_e2.assign(Ncell, -1);
    map_out->next_e2me1.assign(Ncell, -1);
    map_out->disp_m.assign(Ncell, 0);
    map_out->disp_n.assign(Ncell, 0);
    map_out->add_table.assign(Ncell * Ncell, -1);
  }

  for (int s = 0; s < Ncell; s++) {
    Coord c = reps[s];
    for (int d = 0; d < 3; d++) {
      int64_t nk = CosetKey(c.m + dm[d], c.n + dn[d], Lx, Ly, Tx, Ty, Ncell);
      auto it = key_to_id.find(nk);
      assert(it != key_to_id.end());
      AddUniqueLink(lattice, s, it->second, wt[d]);

      if (map_out != nullptr) {
        if (d == 0) map_out->next_e1[s] = it->second;
        if (d == 1) map_out->next_e2[s] = it->second;
        if (d == 2) map_out->next_e2me1[s] = it->second;
      }
    }

    if (map_out != nullptr) {
      map_out->disp_m[s] = c.m;
      map_out->disp_n[s] = c.n;
    }
  }

  if (map_out != nullptr) {
    for (int s = 0; s < Ncell; s++) {
      Coord c = reps[s];
      for (int d = 0; d < Ncell; d++) {
        Coord dd = reps[d];
        int64_t k = CosetKey(c.m + dd.m, c.n + dd.n, Lx, Ly, Tx, Ty, Ncell);
        auto it = key_to_id.find(k);
        assert(it != key_to_id.end());
        map_out->add_table[s * Ncell + d] = it->second;
      }
    }
  }

  lattice->UpdateDistinct();
}

int OrbitPeriod(const std::vector<int>& next) {
  int V = int(next.size());
  if (V == 0) return 0;
  int t = 0;
  int p = 0;
  do {
    t = next[t];
    p++;
  } while (t != 0 && p <= V);
  assert(p <= V);
  return p;
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

void JackknifeBinderU4(const std::vector<double>& mag_samples,
                       double* mean_out, double* err_out) {
  const int n = int(mag_samples.size());
  if (n <= 1) {
    *mean_out = 0.0;
    *err_out = 0.0;
    return;
  }

  double sum_m2 = 0.0;
  double sum_m4 = 0.0;
  for (int i = 0; i < n; i++) {
    const double m2 = mag_samples[i] * mag_samples[i];
    sum_m2 += m2;
    sum_m4 += m2 * m2;
  }

  const double m2_mean = sum_m2 / double(n);
  const double m4_mean = sum_m4 / double(n);
  if (m2_mean <= 0.0) {
    *mean_out = 0.0;
    *err_out = 0.0;
    return;
  }

  const double u4_full = 1.5 * (1.0 - m4_mean / (3.0 * m2_mean * m2_mean));

  double sumsq = 0.0;
  for (int i = 0; i < n; i++) {
    const double m2 = mag_samples[i] * mag_samples[i];
    const double m4 = m2 * m2;
    const double m2_leave = (sum_m2 - m2) / double(n - 1);
    const double m4_leave = (sum_m4 - m4) / double(n - 1);
    if (m2_leave <= 0.0) continue;
    const double u4_leave = 1.5 * (1.0 - m4_leave / (3.0 * m2_leave * m2_leave));
    const double d = u4_leave - u4_full;
    sumsq += d * d;
  }

  *mean_out = u4_full;
  *err_out = std::sqrt((double(n) - 1.0) / double(n) * sumsq);
}

void JackknifeMean(const std::vector<double>& samples,
                   double* mean_out, double* err_out) {
  const int n = int(samples.size());
  if (n <= 1) {
    *mean_out = 0.0;
    *err_out = 0.0;
    return;
  }

  double sum = 0.0;
  for (double x : samples) sum += x;
  const double mean = sum / double(n);

  double sumsq = 0.0;
  for (int i = 0; i < n; i++) {
    const double leave = (sum - samples[i]) / double(n - 1);
    const double d = leave - mean;
    sumsq += d * d;
  }

  *mean_out = mean;
  *err_out = std::sqrt((double(n) - 1.0) / double(n) * sumsq);
}

void JackknifeBinderFromMoments(const std::vector<double>& m2_samples,
                                const std::vector<double>& m4_samples,
                                double* mean_out, double* err_out) {
  const int n = int(m2_samples.size());
  if (n <= 1 || n != int(m4_samples.size())) {
    *mean_out = 0.0;
    *err_out = 0.0;
    return;
  }

  double sum_m2 = 0.0;
  double sum_m4 = 0.0;
  for (int i = 0; i < n; i++) {
    sum_m2 += m2_samples[i];
    sum_m4 += m4_samples[i];
  }

  const double m2 = sum_m2 / double(n);
  const double m4 = sum_m4 / double(n);
  if (m2 <= 0.0) {
    *mean_out = 0.0;
    *err_out = 0.0;
    return;
  }

  const double u4_full = 1.5 * (1.0 - m4 / (3.0 * m2 * m2));

  double sumsq = 0.0;
  for (int i = 0; i < n; i++) {
    const double m2_leave = (sum_m2 - m2_samples[i]) / double(n - 1);
    const double m4_leave = (sum_m4 - m4_samples[i]) / double(n - 1);
    if (m2_leave <= 0.0) continue;
    const double u4_leave = 1.5 * (1.0 - m4_leave / (3.0 * m2_leave * m2_leave));
    const double d = u4_leave - u4_full;
    sumsq += d * d;
  }

  *mean_out = u4_full;
  *err_out = std::sqrt((double(n) - 1.0) / double(n) * sumsq);
}

int UFRoot(std::vector<int>* parent, int x) {
  int r = x;
  while ((*parent)[r] != r) r = (*parent)[r];
  while (x != r) {
    int p = (*parent)[x];
    (*parent)[x] = r;
    x = p;
  }
  return r;
}

void UFJoin(std::vector<int>* parent, std::vector<int>* rank, int a, int b) {
  int ra = UFRoot(parent, a);
  int rb = UFRoot(parent, b);
  if (ra == rb) return;
  if ((*rank)[ra] < (*rank)[rb]) {
    int t = ra;
    ra = rb;
    rb = t;
  }
  (*parent)[rb] = ra;
  if ((*rank)[ra] == (*rank)[rb]) (*rank)[ra]++;
}

// Improved moments from one sampled FK bond graph at fixed spins.
// For cluster sizes |C| = s_i, with V sites:
//   m2_imp = (sum_i s_i^2) / V^2
//   m4_imp = (3 (sum_i s_i^2)^2 - 2 sum_i s_i^4) / V^4
void SampleFKImprovedEstimator(QfeIsing* field,
                               double* m2_imp_out,
                               double* m4_imp_out) {
  const int V = field->lattice->n_sites;
  std::vector<int> parent(V);
  std::vector<int> rank(V, 0);
  for (int s = 0; s < V; s++) parent[s] = s;

  for (int l = 0; l < field->lattice->n_links; l++) {
    int s1 = field->lattice->links[l].sites[0];
    int s2 = field->lattice->links[l].sites[1];
    double w = field->lattice->links[l].wt;

    // Same bond criterion as Wolff/SW implementation.
    const double rate = -2.0 * field->beta * field->spin[s1] * field->spin[s2] * w;
    if (rate >= 0.0) continue;
    if (field->lattice->rng.RandReal() < std::exp(rate)) continue;
    UFJoin(&parent, &rank, s1, s2);
  }

  std::vector<int> csize(V, 0);
  for (int s = 0; s < V; s++) {
    int r = UFRoot(&parent, s);
    csize[r]++;
  }

  double sum2 = 0.0;
  double sum4 = 0.0;
  for (int s = 0; s < V; s++) {
    if (csize[s] == 0) continue;
    const double x = double(csize[s]);
    const double x2 = x * x;
    sum2 += x2;
    sum4 += x2 * x2;
  }

  const double V2 = double(V) * double(V);
  const double V4 = V2 * V2;
  *m2_imp_out = (V2 > 0.0) ? (sum2 / V2) : 0.0;
  *m4_imp_out = (V4 > 0.0) ? ((3.0 * sum2 * sum2 - 2.0 * sum4) / V4) : 0.0;
}

}  // namespace

int main(int argc, char* argv[]) {
  // Twisted parallelogram vectors: v=(L_x,T_y), u=(T_x,-L_y).
  int L_x = 5;
  int L_y = 5;
  int T_x = 0;
  int T_y = 0;

  // Directional couplings on triangular links (+e1, +e2, +(e2-e1)).
  double K1 = 0.5;
  double K2 = 0.5;
  double K3 = 0.5;
  double Kt = 0.5;
  double beta = 1.0;
  int N_t = 1;

  unsigned int seed = 1234u;
  int n_therm = 2000;
  int n_traj = 50000;
  int n_skip = 20;
  int n_wolff = 3;
  int n_metropolis = 5;
  double wall_time = 0.0;
  std::string data_dir = "ising_tri_twisted_parallelogram";
  std::string two_point_name = "two_point_typed.dat";
  std::string full_two_point_name = "two_point_all_to_all.dat";
  std::string time_two_point_name = "two_point_time.dat";
  std::string full_3d_two_point_name = "two_point_3d.dat";

  const struct option long_options[] = {
      {"L_x", required_argument, 0, 'X'},
      {"L_y", required_argument, 0, 'Y'},
      {"T_x", required_argument, 0, 'P'},
      {"T_y", required_argument, 0, 'Q'},
      {"N_t", required_argument, 0, 'Z'},
      {"k_x", required_argument, 0, 'a'},
      {"k_y", required_argument, 0, 'b'},
      {"k_xy", required_argument, 0, 'c'},
      {"k1", required_argument, 0, 'a'},
      {"k2", required_argument, 0, 'b'},
      {"k3", required_argument, 0, 'c'},
      {"k_t", required_argument, 0, 'g'},
      {"beta", required_argument, 0, 'B'},
      {"seed", required_argument, 0, 'S'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'e'},
      {"wall_time", required_argument, 0, 'W'},
      {"data_dir", required_argument, 0, 'd'},
      {"two_point_name", required_argument, 0, 'p'},
      {"full_two_point_name", required_argument, 0, 'f'},
      {"time_two_point_name", required_argument, 0, 'T'},
      {"full_3d_two_point_name", required_argument, 0, 'F'},
      {0, 0, 0, 0}};

  const char* short_options = "X:Y:P:Q:Z:a:b:c:g:B:S:h:t:s:w:e:W:d:p:f:T:F:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'X': L_x = atoi(optarg); break;
      case 'Y': L_y = atoi(optarg); break;
      case 'P': T_x = atoi(optarg); break;
      case 'Q': T_y = atoi(optarg); break;
      case 'Z': N_t = atoi(optarg); break;
      case 'a': K1 = std::stod(optarg); break;
      case 'b': K2 = std::stod(optarg); break;
      case 'c': K3 = std::stod(optarg); break;
      case 'g': Kt = std::stod(optarg); break;
      case 'B': beta = std::stod(optarg); break;
      case 'S': seed = static_cast<unsigned int>(std::strtoul(optarg, nullptr, 10)); break;
      case 'h': n_therm = atoi(optarg); break;
      case 't': n_traj = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'w': n_wolff = atoi(optarg); break;
      case 'e': n_metropolis = atoi(optarg); break;
      case 'W': wall_time = std::stod(optarg); break;
      case 'd': data_dir = optarg; break;
      case 'p': two_point_name = optarg; break;
      case 'f': full_two_point_name = optarg; break;
      case 'T': time_two_point_name = optarg; break;
      case 'F': full_3d_two_point_name = optarg; break;
      default: break;
    }
  }

  const int N_cell = L_x * L_y + T_x * T_y;
  if (L_x <= 0 || L_y <= 0 || N_cell <= 0 || N_t <= 0) {
    fprintf(stderr, "ERROR: require L_x>0, L_y>0, N_t>0, and L_x*L_y + T_x*T_y > 0\n");
    return 1;
  }

  printf("L_x: %d\n", L_x);
  printf("L_y: %d\n", L_y);
  printf("T_x: %d\n", T_x);
  printf("T_y: %d\n", T_y);
  printf("N_cell: %d\n", N_cell);
  printf("N_t: %d\n", N_t);
  printf("K1: %.12f\n", K1);
  printf("K2: %.12f\n", K2);
  printf("K3: %.12f\n", K3);
  printf("Kt: %.12f\n", Kt);
  printf("beta: %.12f\n", beta);
  printf("seed: 0x%08X\n", seed);
  printf("n_therm: %d\n", n_therm);
  printf("n_traj: %d\n", n_traj);
  printf("n_skip: %d\n", n_skip);
  printf("n_wolff: %d\n", n_wolff);
  printf("n_metropolis: %d\n", n_metropolis);
  printf("wall_time: %f\n", wall_time);
  printf("data_dir: %s\n", data_dir.c_str());
  printf("two_point_name: %s\n", two_point_name.c_str());
  printf("full_two_point_name: %s\n", full_two_point_name.c_str());
  printf("time_two_point_name: %s\n", time_two_point_name.c_str());
  printf("full_3d_two_point_name: %s\n", full_3d_two_point_name.c_str());

  QfeLattice lattice;
  TwistedMap twisted_map;
  lattice.SeedRng(seed);
  InitTriangleTwistedParallelogram(&lattice, L_x, L_y, T_x, T_y, K1, K2, K3,
                                   &twisted_map);

  if (N_t > 1) {
    lattice.AddDimension(N_t);
    for (int s = 0; s < lattice.n_sites; s++) {
      int sp1 = (s + N_cell) % lattice.n_sites;
      int l = lattice.FindLink(s, sp1);
      if (l != -1) lattice.links[l].wt = Kt;
    }
  }

  QfeIsing field(&lattice, beta);
  field.HotStart();

  printf("n_sites: %d\n", lattice.n_sites);
  printf("n_links: %d\n", lattice.n_links);
  printf("initial action: %.12f\n", field.Action());

  const int p_e1 = OrbitPeriod(twisted_map.next_e1);
  const int p_e2 = OrbitPeriod(twisted_map.next_e2);
  const int p_e2me1 = OrbitPeriod(twisted_map.next_e2me1);
  printf("periods e1/e2/e2-e1: %d %d %d\n", p_e1, p_e2, p_e2me1);

  QfeMeasReal mag;
  QfeMeasReal mag_2;
  QfeMeasReal mag_4;
  QfeMeasReal energy;
  QfeMeasReal energy_2;
  std::vector<QfeMeasReal> corr_e1(p_e1);
  std::vector<QfeMeasReal> corr_e2(p_e2);
  std::vector<QfeMeasReal> corr_e2me1(p_e2me1);
  std::vector<QfeMeasReal> corr_all(N_cell);
  std::vector<std::vector<double>> corr_e1_samples(p_e1);
  std::vector<std::vector<double>> corr_e2_samples(p_e2);
  std::vector<std::vector<double>> corr_e2me1_samples(p_e2me1);
  std::vector<std::vector<double>> corr_all_samples(N_cell);
  // Temporal correlator: G_time(dt) = (1/V) sum_{t0,x} spin(t0,x)*spin(t0+dt,x)
  std::vector<QfeMeasReal> corr_time(N_t);
  std::vector<std::vector<double>> corr_time_samples(N_t);
  // Full 3D correlator: G_3d[dt*N_cell + d] = (1/V) sum_{t0,x} spin(t0,x)*spin(t0+dt, x+d)
  // Uses QfeMeasReal online stats (no per-sample storage to keep memory bounded).
  std::vector<QfeMeasReal> corr_3d(N_t * N_cell);
  std::vector<double> mag_samples;
  std::vector<double> m2_imp_samples;
  std::vector<double> m4_imp_samples;
  QfeMeasReal cluster_size;
  QfeMeasReal accept_metropolis;

  Timer timer;

    const int total_traj = n_traj + n_therm;
    const int progress_step = (total_traj / 20 > 0) ? (total_traj / 20) : 1;
    int measured_count = 0;
    int next_progress = progress_step;

    printf("progress_step: every %d trajectories (total %d)\n",
      progress_step, total_traj);
    fflush(stdout);

  for (int n = 0; n < (n_traj + n_therm); n++) {
    int cluster_size_sum = 0;
    for (int j = 0; j < n_wolff; j++) {
      cluster_size_sum += field.WolffUpdate();
    }
    double metropolis_sum = 0.0;
    for (int j = 0; j < n_metropolis; j++) {
      metropolis_sum += field.Metropolis();
    }
    cluster_size.Measure(double(cluster_size_sum) / double(lattice.n_sites));
    accept_metropolis.Measure(metropolis_sum);

    if (n % n_skip || n < n_therm) {
      const int traj_done = n + 1;
      if (traj_done >= next_progress || traj_done == total_traj) {
        const double pct = 100.0 * double(traj_done) / double(total_traj);
        printf("progress traj=%d/%d (%.1f%%) measured=%d\n",
               traj_done, total_traj, pct, measured_count);
        fflush(stdout);
        while (next_progress <= traj_done) next_progress += progress_step;
      }
      continue;
    }

    measured_count++;

    double m_sum = 0.0;
    for (int s = 0; s < lattice.n_sites; s++) m_sum += field.spin[s];
    double m_mean = m_sum / double(lattice.n_sites);

    double e_sum = 0.0;
    for (int l = 0; l < lattice.n_links; l++) {
      int s_a = lattice.links[l].sites[0];
      int s_b = lattice.links[l].sites[1];
      e_sum -= field.spin[s_a] * field.spin[s_b] * lattice.links[l].wt;
    }
    double e_mean = e_sum / double(lattice.n_sites);

    double m2 = m_mean * m_mean;
    mag_samples.push_back(m_mean);
    mag.Measure(fabs(m_mean));
    mag_2.Measure(m2);
    mag_4.Measure(m2 * m2);

    double m2_imp_sample = 0.0;
    double m4_imp_sample = 0.0;
    SampleFKImprovedEstimator(&field, &m2_imp_sample, &m4_imp_sample);
    m2_imp_samples.push_back(m2_imp_sample);
    m4_imp_samples.push_back(m4_imp_sample);

    energy.Measure(e_mean);
    energy_2.Measure(e_mean * e_mean);

    std::vector<double> corr_sum_e1(p_e1, 0.0);
    std::vector<double> corr_sum_e2(p_e2, 0.0);
    std::vector<double> corr_sum_e2me1(p_e2me1, 0.0);
    std::vector<double> corr_sum_all(N_cell, 0.0);

    for (int t0 = 0; t0 < N_t; t0++) {
      int t_offset = t0 * N_cell;
      for (int s2 = 0; s2 < N_cell; s2++) {
        int s = t_offset + s2;
        double spin_s = field.spin[s];

        int u = s2;
        for (int r = 0; r < p_e1; r++) {
          corr_sum_e1[r] += spin_s * field.spin[t_offset + u];
          u = twisted_map.next_e1[u];
        }

        u = s2;
        for (int r = 0; r < p_e2; r++) {
          corr_sum_e2[r] += spin_s * field.spin[t_offset + u];
          u = twisted_map.next_e2[u];
        }

        u = s2;
        for (int r = 0; r < p_e2me1; r++) {
          corr_sum_e2me1[r] += spin_s * field.spin[t_offset + u];
          u = twisted_map.next_e2me1[u];
        }

        int add_base = s2 * N_cell;
        for (int d = 0; d < N_cell; d++) {
          int u_all = twisted_map.add_table[add_base + d];
          corr_sum_all[d] += spin_s * field.spin[t_offset + u_all];
        }
      }
    }

    const double norm = double(lattice.n_sites);
    for (int r = 0; r < p_e1; r++) {
      double c = corr_sum_e1[r] / norm;
      corr_e1[r].Measure(c);
      corr_e1_samples[r].push_back(c);
    }
    for (int r = 0; r < p_e2; r++) {
      double c = corr_sum_e2[r] / norm;
      corr_e2[r].Measure(c);
      corr_e2_samples[r].push_back(c);
    }
    for (int r = 0; r < p_e2me1; r++) {
      double c = corr_sum_e2me1[r] / norm;
      corr_e2me1[r].Measure(c);
      corr_e2me1_samples[r].push_back(c);
    }
    for (int d = 0; d < N_cell; d++) {
      double c = corr_sum_all[d] / norm;
      corr_all[d].Measure(c);
      corr_all_samples[d].push_back(c);
    }

    // ── Temporal correlator G_time(dt) ──────────────────────────────────────
    // G_time(dt) = (1/V) sum_{t0,x} spin(t0,x) * spin((t0+dt)%Nt, x)
    {
      std::vector<double> ctime_sum(N_t, 0.0);
      for (int t0 = 0; t0 < N_t; t0++) {
        const int t0_off = t0 * N_cell;
        for (int dt = 0; dt < N_t; dt++) {
          const int t1_off = ((t0 + dt) % N_t) * N_cell;
          for (int s2 = 0; s2 < N_cell; s2++) {
            ctime_sum[dt] += field.spin[t0_off + s2] * field.spin[t1_off + s2];
          }
        }
      }
      for (int dt = 0; dt < N_t; dt++) {
        const double c = ctime_sum[dt] / norm;
        corr_time[dt].Measure(c);
        corr_time_samples[dt].push_back(c);
      }
    }

    // ── Full 3D all-to-all correlator G_3d(d, dt) ───────────────────────────
    // G_3d(d, dt) = (1/V) sum_{t0,x} spin(t0,x) * spin((t0+dt)%Nt, x+d)
    // Loops ordered for cache efficiency: innermost d-loop reads from a
    // contiguous spin slice of size N_cell (fits in L1 cache).
    {
      std::vector<double> c3d_sum(N_t * N_cell, 0.0);
      for (int t0 = 0; t0 < N_t; t0++) {
        const int t0_off = t0 * N_cell;
        for (int dt = 0; dt < N_t; dt++) {
          const int t1_off = ((t0 + dt) % N_t) * N_cell;
          double* acc = c3d_sum.data() + dt * N_cell;
          for (int s2 = 0; s2 < N_cell; s2++) {
            const double spin_s = field.spin[t0_off + s2];
            const int* trow = twisted_map.add_table.data() + s2 * N_cell;
            for (int d = 0; d < N_cell; d++) {
              acc[d] += spin_s * field.spin[t1_off + trow[d]];
            }
          }
        }
      }
      for (int i = 0; i < N_t * N_cell; i++) {
        corr_3d[i].Measure(c3d_sum[i] / norm);
      }
    }

    const int traj_done = n + 1;
    if (traj_done >= next_progress || traj_done == total_traj) {
      const double pct = 100.0 * double(traj_done) / double(total_traj);
      printf("progress traj=%d/%d (%.1f%%) measured=%d\n",
             traj_done, total_traj, pct, measured_count);
      fflush(stdout);
      while (next_progress <= traj_done) next_progress += progress_step;
    }

    if (wall_time > 0.0 && timer.Duration() > wall_time) break;
  }

  timer.Stop();
  printf("duration: %.6f\n", timer.Duration());
  printf("cluster_size/V: %.4f\n", cluster_size.Mean());
  printf("accept_metropolis: %.4f\n", accept_metropolis.Mean());

  double m_mean = mag.Mean();
  double m_err = mag.Error();
  double m2_mean = mag_2.Mean();
  double m2_err = mag_2.Error();

  double U4_spin_mean = 0.0;
  double U4_spin_err = 0.0;
  JackknifeBinderU4(mag_samples, &U4_spin_mean, &U4_spin_err);

  double m2_imp_mean = 0.0;
  double m2_imp_err = 0.0;
  JackknifeMean(m2_imp_samples, &m2_imp_mean, &m2_imp_err);

  double U4_imp_mean = 0.0;
  double U4_imp_err = 0.0;
  JackknifeBinderFromMoments(m2_imp_samples, m4_imp_samples,
                             &U4_imp_mean, &U4_imp_err);

  // Primary reported Binder for scans/plots uses the cluster-improved estimator.
  const double U4_mean = U4_imp_mean;
  const double U4_err = U4_imp_err;
  printf("U4 (cluster-improved): %.12e %.12e\n", U4_mean, U4_err);
  printf("U4 (spin-jackknife): %.12e %.12e\n", U4_spin_mean, U4_spin_err);
  printf("m2 (cluster-improved): %.12e %.12e\n", m2_imp_mean, m2_imp_err);

  double m_susc_mean = m2_mean - m_mean * m_mean;
  double m_susc_err = sqrt(pow(m2_err, 2.0) + pow(2.0 * m_mean * m_err, 2.0));
  printf("m_susc: %.12e %.12e\n", m_susc_mean, m_susc_err);

  if (!MakeDir(data_dir)) {
    fprintf(stderr, "WARNING: could not create data_dir: %s\n", data_dir.c_str());
    return 0;
  }

  char run_id[128];
  char subdir[256];
  char path[512];
    sprintf(run_id, "Lx%d_Ly%d_Tx%d_Ty%d_Nt%d_k%.3f_%.3f_%.3f_kt%.3f",
      L_x, L_y, T_x, T_y, N_t, K1, K2, K3, Kt);
  sprintf(subdir, "%s/%s", data_dir.c_str(), run_id);
  if (!MakeDir(subdir)) {
    fprintf(stderr, "WARNING: could not create run subdir: %s\n", subdir);
    return 0;
  }
  sprintf(path, "%s/%s_%08X.dat", subdir, run_id, seed);

  printf("opening file: %s\n", path);
  FILE* file = fopen(path, "w");
  assert(file != nullptr);

  fprintf(file, "L_x %d\n", L_x);
  fprintf(file, "L_y %d\n", L_y);
  fprintf(file, "T_x %d\n", T_x);
  fprintf(file, "T_y %d\n", T_y);
  fprintf(file, "N_cell %d\n", N_cell);
  fprintf(file, "N_t %d\n", N_t);
  fprintf(file, "K1 %.12e\n", K1);
  fprintf(file, "K2 %.12e\n", K2);
  fprintf(file, "K3 %.12e\n", K3);
  fprintf(file, "Kt %.12e\n", Kt);
  fprintf(file, "beta %.12e\n", beta);
  fprintf(file, "mag %.12e %.12e %d\n", mag.Mean(), mag.Error(), mag.n);
  fprintf(file, "mag^2 %.12e %.12e %d\n", mag_2.Mean(), mag_2.Error(), mag_2.n);
  fprintf(file, "mag^4 %.12e %.12e %d\n", mag_4.Mean(), mag_4.Error(), mag_4.n);
  fprintf(file, "energy %.12e %.12e %d\n", energy.Mean(), energy.Error(), energy.n);
  fprintf(file, "energy^2 %.12e %.12e %d\n", energy_2.Mean(), energy_2.Error(), energy_2.n);
  fprintf(file, "U4 %.12e %.12e\n", U4_mean, U4_err);
  fprintf(file, "U4_spin %.12e %.12e\n", U4_spin_mean, U4_spin_err);
  fprintf(file, "mag^2_imp %.12e %.12e %d\n", m2_imp_mean, m2_imp_err, int(m2_imp_samples.size()));
  fprintf(file, "U4_imp %.12e %.12e %d\n", U4_imp_mean, U4_imp_err, int(m4_imp_samples.size()));
  fprintf(file, "m_susc %.12e %.12e\n", m_susc_mean, m_susc_err);
  fclose(file);

  char two_point_path[512];
  sprintf(two_point_path, "%s/%s", subdir, two_point_name.c_str());
  FILE* corr_file = fopen(two_point_path, "w");
  assert(corr_file != nullptr);
  // Typed format: type r mean err
  // type 0: e1, type 1: e2, type 2: e2-e1.
  // type 4: connected e1, type 5: connected e2, type 6: connected e2-e1.
  for (int r = 0; r < p_e1; r++) {
    fprintf(corr_file, "0 %d %.16e %.16e\n", r, corr_e1[r].Mean(), corr_e1[r].Error());
  }
  for (int r = 0; r < p_e2; r++) {
    fprintf(corr_file, "1 %d %.16e %.16e\n", r, corr_e2[r].Mean(), corr_e2[r].Error());
  }
  for (int r = 0; r < p_e2me1; r++) {
    fprintf(corr_file, "2 %d %.16e %.16e\n", r, corr_e2me1[r].Mean(), corr_e2me1[r].Error());
  }

  for (int r = 0; r < p_e1; r++) {
    double mean_conn = 0.0;
    double err_conn = 0.0;
    JackknifeConnectedCorr(corr_e1_samples[r], mag_samples, &mean_conn, &err_conn);
    fprintf(corr_file, "4 %d %.16e %.16e\n", r, mean_conn, err_conn);
  }
  for (int r = 0; r < p_e2; r++) {
    double mean_conn = 0.0;
    double err_conn = 0.0;
    JackknifeConnectedCorr(corr_e2_samples[r], mag_samples, &mean_conn, &err_conn);
    fprintf(corr_file, "5 %d %.16e %.16e\n", r, mean_conn, err_conn);
  }
  for (int r = 0; r < p_e2me1; r++) {
    double mean_conn = 0.0;
    double err_conn = 0.0;
    JackknifeConnectedCorr(corr_e2me1_samples[r], mag_samples, &mean_conn, &err_conn);
    fprintf(corr_file, "6 %d %.16e %.16e\n", r, mean_conn, err_conn);
  }
  fclose(corr_file);
  printf("wrote: %s\n", two_point_path);

  char full_two_point_path[512];
  sprintf(full_two_point_path, "%s/%s", subdir, full_two_point_name.c_str());
  FILE* full_file = fopen(full_two_point_path, "w");
  assert(full_file != nullptr);
  fprintf(full_file, "# Equal-time spatial all-to-all average over all time slices.\n");
  fprintf(full_file, "# d m n corr err corr_conn err_conn\n");
  for (int d = 0; d < N_cell; d++) {
    double cconn = 0.0;
    double econn = 0.0;
    JackknifeConnectedCorr(corr_all_samples[d], mag_samples, &cconn, &econn);
    fprintf(full_file, "%d %d %d %.16e %.16e %.16e %.16e\n",
            d,
            twisted_map.disp_m[d],
            twisted_map.disp_n[d],
            corr_all[d].Mean(),
            corr_all[d].Error(),
            cconn,
            econn);
  }
  fclose(full_file);
  printf("wrote: %s\n", full_two_point_path);

  return 0;
}
