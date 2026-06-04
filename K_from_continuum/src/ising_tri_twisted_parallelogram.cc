#include <getopt.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _WIN32
#  include <direct.h>
#  define MKDIR(p) ::_mkdir(p)
#else
#  include <sys/stat.h>
#  include <sys/types.h>
#  define MKDIR(p) ::mkdir((p), 0755)
#endif

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
  std::vector<int> add_selected;
  std::vector<int> selected_disp_ids;
  std::vector<int> selected_disp_m;
  std::vector<int> selected_disp_n;
  int single_disp_id = -1;
  int single_disp_m = 0;
  int single_disp_n = 0;
};

std::vector<Coord> ParseSelectedDispList(const std::string& raw) {
  std::vector<Coord> selected;
  if (raw.empty()) return selected;

  std::stringstream stream(raw);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (token.empty()) continue;
    size_t sep = token.find(':');
    if (sep == std::string::npos) {
      fprintf(stderr, "ERROR: invalid selected displacement token '%s' (expected m:n)\n", token.c_str());
      exit(1);
    }
    Coord coord;
    coord.m = std::atoi(token.substr(0, sep).c_str());
    coord.n = std::atoi(token.substr(sep + 1).c_str());
    selected.push_back(coord);
  }
  return selected;
}

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
  if (MKDIR(path.c_str()) == 0) return true;
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
                                      TwistedMap* map_out,
                                      bool selected_disp_only,
                                      const std::vector<Coord>& selected_disps) {
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
    map_out->add_table.clear();
    map_out->add_selected.clear();
    map_out->selected_disp_ids.clear();
    map_out->selected_disp_m.clear();
    map_out->selected_disp_n.clear();
    map_out->single_disp_id = -1;
    map_out->single_disp_m = 0;
    map_out->single_disp_n = 0;
    if (selected_disp_only) {
      map_out->add_selected.assign(Ncell * int(selected_disps.size()), -1);
    } else {
      map_out->add_table.assign(Ncell * Ncell, -1);
    }
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
    if (selected_disp_only) {
      map_out->selected_disp_ids.reserve(selected_disps.size());
      map_out->selected_disp_m.reserve(selected_disps.size());
      map_out->selected_disp_n.reserve(selected_disps.size());

      for (size_t d = 0; d < selected_disps.size(); d++) {
        const Coord& coord = selected_disps[d];
        int64_t selected_key = CosetKey(coord.m, coord.n, Lx, Ly, Tx, Ty, Ncell);
        auto selected_it = key_to_id.find(selected_key);
        assert(selected_it != key_to_id.end());
        const int selected_id = selected_it->second;
        map_out->selected_disp_ids.push_back(selected_id);
        map_out->selected_disp_m.push_back(reps[selected_id].m);
        map_out->selected_disp_n.push_back(reps[selected_id].n);
        if (d == 0) {
          map_out->single_disp_id = selected_id;
          map_out->single_disp_m = reps[selected_id].m;
          map_out->single_disp_n = reps[selected_id].n;
        }

        for (int s = 0; s < Ncell; s++) {
          Coord c = reps[s];
          int64_t k = CosetKey(c.m + coord.m, c.n + coord.n, Lx, Ly, Tx, Ty, Ncell);
          auto it = key_to_id.find(k);
          assert(it != key_to_id.end());
          map_out->add_selected[s * int(selected_disps.size()) + int(d)] = it->second;
        }
      }
    } else {
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
  int single_disp_m = 0;
  int single_disp_n = 0;
  bool have_single_disp_m = false;
  bool have_single_disp_n = false;
  std::string single_disp_samples_name;
  std::string selected_disp_list_raw;
  std::string selected_disp_bundle_name;
  // Speedup 4 (FS reweighting): when nonzero, dump one (E_per_site, |m|, m^2)
  // line per accepted trajectory to <subdir>/traces_<seed>.dat so the Python
  // side can reweight chi(beta) over a window without rerunning MC.
  int dump_traces = 0;

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
      {"single_disp_m", required_argument, 0, 'm'},
      {"single_disp_n", required_argument, 0, 'n'},
      {"single_disp_samples_name", required_argument, 0, 'J'},
        {"selected_disp_list", required_argument, 0, 'M'},
        {"selected_disp_bundle_name", required_argument, 0, 'K'},
      {"dump_traces", required_argument, 0, 'D'},
      {0, 0, 0, 0}};

      const char* short_options = "X:Y:P:Q:Z:a:b:c:g:B:S:h:t:s:w:e:W:d:p:f:m:n:J:M:K:D:";

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
      case 'm': single_disp_m = atoi(optarg); have_single_disp_m = true; break;
      case 'n': single_disp_n = atoi(optarg); have_single_disp_n = true; break;
      case 'J': single_disp_samples_name = optarg; break;
      case 'M': selected_disp_list_raw = optarg; break;
      case 'K': selected_disp_bundle_name = optarg; break;
      case 'D': dump_traces = atoi(optarg); break;
      default: break;
    }
  }

  if (have_single_disp_m != have_single_disp_n) {
    fprintf(stderr, "ERROR: --single_disp_m and --single_disp_n must be provided together\n");
    return 1;
  }
  if (!selected_disp_list_raw.empty() && (have_single_disp_m || have_single_disp_n)) {
    fprintf(stderr, "ERROR: choose either --single_disp_m/--single_disp_n or --selected_disp_list, not both\n");
    return 1;
  }

  std::vector<Coord> selected_disps;
  if (have_single_disp_m && have_single_disp_n) {
    selected_disps.push_back({single_disp_m, single_disp_n});
  } else {
    selected_disps = ParseSelectedDispList(selected_disp_list_raw);
  }

  const bool selected_disp_only = !selected_disps.empty();
  if (!selected_disp_only && !single_disp_samples_name.empty()) {
    fprintf(stderr, "ERROR: --single_disp_samples_name requires --single_disp_m and --single_disp_n\n");
    return 1;
  }
  if (!selected_disp_only && !selected_disp_bundle_name.empty()) {
    fprintf(stderr, "ERROR: --selected_disp_bundle_name requires --selected_disp_list\n");
    return 1;
  }
  if (!single_disp_samples_name.empty() && !selected_disp_bundle_name.empty()) {
    fprintf(stderr, "ERROR: choose either --single_disp_samples_name or --selected_disp_bundle_name\n");
    return 1;
  }
  if (selected_disp_only && selected_disps.size() > 1 && !single_disp_samples_name.empty()) {
    fprintf(stderr, "ERROR: --single_disp_samples_name supports exactly one selected displacement\n");
    return 1;
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
  if (selected_disp_only) {
    for (size_t d = 0; d < selected_disps.size(); d++) {
      printf("selected_disp[%zu]: (%d, %d)\n", d, selected_disps[d].m, selected_disps[d].n);
    }
  }
  if (!single_disp_samples_name.empty()) {
    printf("single_disp_samples_name: %s\n", single_disp_samples_name.c_str());
  }
  if (!selected_disp_bundle_name.empty()) {
    printf("selected_disp_bundle_name: %s\n", selected_disp_bundle_name.c_str());
  }

  QfeLattice lattice;
  TwistedMap twisted_map;
  lattice.SeedRng(seed);
  InitTriangleTwistedParallelogram(&lattice, L_x, L_y, T_x, T_y, K1, K2, K3,
                                   &twisted_map,
                                   selected_disp_only,
                                   selected_disps);

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
  const bool compute_typed_corr = !selected_disp_only;
  const int n_selected_disp = selected_disp_only ? int(selected_disps.size()) : 0;
  printf("periods e1/e2/e2-e1: %d %d %d\n", p_e1, p_e2, p_e2me1);

  QfeMeasReal mag;
  QfeMeasReal mag_2;
  QfeMeasReal mag_4;
  QfeMeasReal energy;
  QfeMeasReal energy_2;
  std::vector<QfeMeasReal> corr_e1(compute_typed_corr ? p_e1 : 0);
  std::vector<QfeMeasReal> corr_e2(compute_typed_corr ? p_e2 : 0);
  std::vector<QfeMeasReal> corr_e2me1(compute_typed_corr ? p_e2me1 : 0);
  std::vector<QfeMeasReal> corr_all(selected_disp_only ? n_selected_disp : N_cell);
  std::vector<std::vector<double>> corr_e1_samples(compute_typed_corr ? p_e1 : 0);
  std::vector<std::vector<double>> corr_e2_samples(compute_typed_corr ? p_e2 : 0);
  std::vector<std::vector<double>> corr_e2me1_samples(compute_typed_corr ? p_e2me1 : 0);
  std::vector<std::vector<double>> corr_all_samples(selected_disp_only ? n_selected_disp : N_cell);
  std::vector<double> mag_samples;
  std::vector<double> abs_mag_samples;
  std::vector<double> m2_samples;
  std::vector<double> e_samples;   // Speedup 4: per-traj energy/site for FS reweighting
  std::vector<double> e1_samples;  // Per-traj directional energy/site for K1 links.
  std::vector<double> e2_samples;  // Per-traj directional energy/site for K2 links.
  std::vector<double> e3_samples;  // Per-traj directional energy/site for K3 links.
  QfeMeasReal cluster_size;
  QfeMeasReal accept_metropolis;

  Timer timer;

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

    if (n % n_skip || n < n_therm) continue;

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
    abs_mag_samples.push_back(fabs(m_mean));
    m2_samples.push_back(m2);
    e_samples.push_back(e_mean);
    mag.Measure(fabs(m_mean));
    mag_2.Measure(m2);
    mag_4.Measure(m2 * m2);
    energy.Measure(e_mean);
    energy_2.Measure(e_mean * e_mean);

    std::vector<double> corr_sum_e1(compute_typed_corr ? p_e1 : 0, 0.0);
    std::vector<double> corr_sum_e2(compute_typed_corr ? p_e2 : 0, 0.0);
    std::vector<double> corr_sum_e2me1(compute_typed_corr ? p_e2me1 : 0, 0.0);
    std::vector<double> corr_sum_all(selected_disp_only ? n_selected_disp : N_cell, 0.0);
    double bond_sum_e1 = 0.0;
    double bond_sum_e2 = 0.0;
    double bond_sum_e3 = 0.0;

    for (int t0 = 0; t0 < N_t; t0++) {
      int t_offset = t0 * N_cell;
      for (int s2 = 0; s2 < N_cell; s2++) {
        int s = t_offset + s2;
        double spin_s = field.spin[s];

        bond_sum_e1 += spin_s * field.spin[t_offset + twisted_map.next_e1[s2]];
        bond_sum_e2 += spin_s * field.spin[t_offset + twisted_map.next_e2[s2]];
        bond_sum_e3 += spin_s * field.spin[t_offset + twisted_map.next_e2me1[s2]];

        if (compute_typed_corr) {
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
        }

        if (selected_disp_only) {
          int add_base = s2 * n_selected_disp;
          for (int d = 0; d < n_selected_disp; d++) {
            int u_all = twisted_map.add_selected[add_base + d];
            corr_sum_all[d] += spin_s * field.spin[t_offset + u_all];
          }
        } else {
          int add_base = s2 * N_cell;
          for (int d = 0; d < N_cell; d++) {
            int u_all = twisted_map.add_table[add_base + d];
            corr_sum_all[d] += spin_s * field.spin[t_offset + u_all];
          }
        }
      }
    }

    const double norm = double(lattice.n_sites);
    double bond_e1 = -bond_sum_e1 / norm;
    double bond_e2 = -bond_sum_e2 / norm;
    double bond_e3 = -bond_sum_e3 / norm;
    e1_samples.push_back(bond_e1);
    e2_samples.push_back(bond_e2);
    e3_samples.push_back(bond_e3);
    if (compute_typed_corr) {
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
    }
    if (selected_disp_only) {
      for (int d = 0; d < n_selected_disp; d++) {
        double c = corr_sum_all[d] / norm;
        corr_all[d].Measure(c);
        corr_all_samples[d].push_back(c);
      }
    } else {
      for (int d = 0; d < N_cell; d++) {
        double c = corr_sum_all[d] / norm;
        corr_all[d].Measure(c);
        corr_all_samples[d].push_back(c);
      }
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
  double m4_mean = mag_4.Mean();
  double m4_err = mag_4.Error();

  double U4_mean = 1.5 * (1.0 - m4_mean / (3.0 * m2_mean * m2_mean));
  double U4_err = 0.5 * U4_mean * sqrt(pow(m4_err / m4_mean, 2.0) + pow(2.0 * m2_err / m2_mean, 2.0));
  printf("U4: %.12e %.12e\n", U4_mean, U4_err);

  double m_susc_mean = m2_mean - m_mean * m_mean;
  // Jackknife error on chi = <m^2> - <|m|>^2
  double m_susc_err = 0.0;
  {
    int n_jack = (int)m2_samples.size();
    double sum_m2 = 0.0, sum_absmag = 0.0;
    for (int i = 0; i < n_jack; i++) {
      sum_m2 += m2_samples[i];
      sum_absmag += abs_mag_samples[i];
    }
    std::vector<double> jack_chi(n_jack);
    for (int i = 0; i < n_jack; i++) {
      double mean_m2_i = (sum_m2 - m2_samples[i]) / double(n_jack - 1);
      double mean_abs_i = (sum_absmag - abs_mag_samples[i]) / double(n_jack - 1);
      jack_chi[i] = mean_m2_i - mean_abs_i * mean_abs_i;
    }
    double jack_mean = 0.0;
    for (int i = 0; i < n_jack; i++) jack_mean += jack_chi[i];
    jack_mean /= double(n_jack);
    double jack_var = 0.0;
    for (int i = 0; i < n_jack; i++) {
      double d = jack_chi[i] - jack_mean;
      jack_var += d * d;
    }
    m_susc_err = sqrt((double(n_jack) - 1.0) / double(n_jack) * jack_var);
  }
  printf("m_susc: %.12e %.12e\n", m_susc_mean, m_susc_err);

  if (!MakeDir(data_dir)) {
    fprintf(stderr, "WARNING: could not create data_dir: %s\n", data_dir.c_str());
    return 0;
  }

  char run_id_buffer[128];
    // Compact format keeps path short (<260 chars) on Windows MAX_PATH.
    // k values stored as integer * 1000 (e.g. 0.850 -> 850).
    snprintf(run_id_buffer, sizeof(run_id_buffer), "%dx%d_t%dx%d_k%d_%d_%d_%d",
      L_x, L_y, T_x, T_y,
      (int)round(K1*1000), (int)round(K2*1000),
      (int)round(K3*1000), (int)round(Kt*1000));
  std::string run_id(run_id_buffer);
  std::string subdir = data_dir + "/" + run_id;
  if (!MakeDir(subdir)) {
    fprintf(stderr, "WARNING: could not create run subdir: %s\n", subdir.c_str());
    return 0;
  }
  char seed_suffix[32];
  snprintf(seed_suffix, sizeof(seed_suffix), "_%08X.dat", seed);
  std::string path = subdir + "/" + run_id + seed_suffix;

  printf("opening file: %s\n", path.c_str());
  FILE* file = fopen(path.c_str(), "w");
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
  fprintf(file, "m_susc %.12e %.12e\n", m_susc_mean, m_susc_err);
  fclose(file);

  std::string two_point_path = subdir + "/" + two_point_name;
  FILE* corr_file = fopen(two_point_path.c_str(), "w");
  assert(corr_file != nullptr);
  // Typed format: type r mean err
  // type 0: e1, type 1: e2, type 2: e2-e1.
  // type 4: connected e1, type 5: connected e2, type 6: connected e2-e1.
  if (compute_typed_corr) {
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
  } else {
    fprintf(corr_file, "# typed correlators omitted in single-displacement mode\n");
  }
  fclose(corr_file);
  printf("wrote: %s\n", two_point_path.c_str());

  std::string full_two_point_path = subdir + "/" + full_two_point_name;
  FILE* full_file = fopen(full_two_point_path.c_str(), "w");
  assert(full_file != nullptr);
  fprintf(full_file, "# Equal-time spatial all-to-all average over all time slices.\n");
  fprintf(full_file, "# d m n corr err corr_conn err_conn\n");
  if (selected_disp_only) {
    for (int idx = 0; idx < n_selected_disp; idx++) {
      double cconn = 0.0;
      double econn = 0.0;
      JackknifeConnectedCorr(corr_all_samples[idx], mag_samples, &cconn, &econn);
      int d = twisted_map.selected_disp_ids[idx];
      fprintf(full_file, "%d %d %d %.16e %.16e %.16e %.16e\n",
              d,
              twisted_map.selected_disp_m[idx],
              twisted_map.selected_disp_n[idx],
              corr_all[idx].Mean(),
              corr_all[idx].Error(),
              cconn,
              econn);
    }
  } else {
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
  }
  fclose(full_file);
  printf("wrote: %s\n", full_two_point_path.c_str());

  if (selected_disp_only && !single_disp_samples_name.empty()) {
    std::string sample_path = subdir + "/" + single_disp_samples_name;
    FILE* sample_file = fopen(sample_path.c_str(), "w");
    if (sample_file == nullptr) {
      fprintf(stderr, "WARNING: could not open single-disp sample file: %s\n", sample_path.c_str());
    } else {
      const int d = twisted_map.single_disp_id;
      const std::vector<double>& corr_samples = corr_all_samples[0];
      assert(corr_samples.size() == mag_samples.size());
      fprintf(sample_file, "# single-displacement corr/mag samples for true ratio jackknife\n");
      fprintf(sample_file, "# seed %u\n", seed);
      fprintf(sample_file, "# L_x %d\n", L_x);
      fprintf(sample_file, "# L_y %d\n", L_y);
      fprintf(sample_file, "# T_x %d\n", T_x);
      fprintf(sample_file, "# T_y %d\n", T_y);
      fprintf(sample_file, "# K1 %.16e\n", K1);
      fprintf(sample_file, "# K2 %.16e\n", K2);
      fprintf(sample_file, "# K3 %.16e\n", K3);
      fprintf(sample_file, "# Kt %.16e\n", Kt);
      fprintf(sample_file, "# beta %.16e\n", beta);
      fprintf(sample_file, "# n_therm %d\n", n_therm);
      fprintf(sample_file, "# n_traj %d\n", n_traj);
      fprintf(sample_file, "# n_skip %d\n", n_skip);
      fprintf(sample_file, "# d %d\n", d);
      fprintf(sample_file, "# m %d\n", twisted_map.disp_m[d]);
      fprintf(sample_file, "# n %d\n", twisted_map.disp_n[d]);
      fprintf(sample_file, "# columns: sample corr mag e1_per_site e2_per_site e3_per_site\n");
      for (size_t i = 0; i < corr_samples.size(); i++) {
        fprintf(sample_file, "%zu %.16e %.16e %.16e %.16e %.16e\n",
                i,
                corr_samples[i],
                mag_samples[i],
                e1_samples[i],
                e2_samples[i],
                e3_samples[i]);
      }
      fclose(sample_file);
      printf("wrote: %s\n", sample_path.c_str());
    }
  }

  if (selected_disp_only && !selected_disp_bundle_name.empty()) {
    std::string bundle_path = subdir + "/" + selected_disp_bundle_name;
    FILE* bundle_file = fopen(bundle_path.c_str(), "w");
    if (bundle_file == nullptr) {
      fprintf(stderr, "WARNING: could not open selected-disp bundle file: %s\n", bundle_path.c_str());
    } else {
      fprintf(bundle_file, "# selected-displacement correlator samples for local geometry-match reweighting\n");
      fprintf(bundle_file, "# seed %u\n", seed);
      fprintf(bundle_file, "# L_x %d\n", L_x);
      fprintf(bundle_file, "# L_y %d\n", L_y);
      fprintf(bundle_file, "# T_x %d\n", T_x);
      fprintf(bundle_file, "# T_y %d\n", T_y);
      fprintf(bundle_file, "# K1 %.16e\n", K1);
      fprintf(bundle_file, "# K2 %.16e\n", K2);
      fprintf(bundle_file, "# K3 %.16e\n", K3);
      fprintf(bundle_file, "# Kt %.16e\n", Kt);
      fprintf(bundle_file, "# beta %.16e\n", beta);
      fprintf(bundle_file, "# n_therm %d\n", n_therm);
      fprintf(bundle_file, "# n_traj %d\n", n_traj);
      fprintf(bundle_file, "# n_skip %d\n", n_skip);
      fprintf(bundle_file, "# n_selected %d\n", n_selected_disp);
      for (int idx = 0; idx < n_selected_disp; idx++) {
        fprintf(bundle_file, "# selected_disp %d %d %d %d\n",
                idx,
                twisted_map.selected_disp_ids[idx],
                twisted_map.selected_disp_m[idx],
                twisted_map.selected_disp_n[idx]);
      }
      fprintf(bundle_file, "# columns: sample");
      for (int idx = 0; idx < n_selected_disp; idx++) {
        fprintf(bundle_file, " corr_%d", idx);
      }
      fprintf(bundle_file, " mag e1_per_site e2_per_site e3_per_site\n");
      const size_t n_samples = mag_samples.size();
      for (size_t i = 0; i < n_samples; i++) {
        fprintf(bundle_file, "%zu", i);
        for (int idx = 0; idx < n_selected_disp; idx++) {
          fprintf(bundle_file, " %.16e", corr_all_samples[idx][i]);
        }
        fprintf(bundle_file, " %.16e %.16e %.16e %.16e\n",
                mag_samples[i],
                e1_samples[i],
                e2_samples[i],
                e3_samples[i]);
      }
      fclose(bundle_file);
      printf("wrote: %s\n", bundle_path.c_str());
    }
  }

  // -------------------------------------------------------------------------
  // Speedup 4 (FS reweighting): per-trajectory traces.
  //
  // Writes one line per measured trajectory:  E_per_site  |m|  m^2
  // The Python side parses this back into NumPy arrays and reweights
  // chi(beta) = V * (<m^2> - <|m|>^2) over a beta window via
  //     w_i(b) = exp( -(b - beta_parent) * E_total_i )
  // with E_total_i = E_per_site_i * V.  Dumping is opt-in to keep the
  // legacy CLI bit-equivalent when not used (--dump_traces 0 = default).
  // -------------------------------------------------------------------------
  if (dump_traces) {
    char traces_name[32];
    snprintf(traces_name, sizeof(traces_name), "traces_%08X.dat", seed);
    std::string traces_path = subdir + "/" + traces_name;
    FILE* tf = fopen(traces_path.c_str(), "w");
    if (tf == nullptr) {
      fprintf(stderr, "WARNING: could not open traces file: %s\n", traces_path.c_str());
    } else {
      fprintf(tf, "# E_per_site  abs_m  m^2  (one line per measured trajectory)\n");
      fprintf(tf, "# n_sites %d  beta %.12e  K1 %.12e  K2 %.12e  K3 %.12e  Kt %.12e\n",
              lattice.n_sites, beta, K1, K2, K3, Kt);
      const int n_meas = (int)m2_samples.size();
      for (int i = 0; i < n_meas; i++) {
        double abs_m_i = (i < (int)abs_mag_samples.size()) ? abs_mag_samples[i] : 0.0;
        double m2_i = m2_samples[i];
        double e_i = (i < (int)e_samples.size()) ? e_samples[i] : 0.0;
        fprintf(tf, "%.12e %.12e %.12e\n", e_i, abs_m_i, m2_i);
      }
      fclose(tf);
      printf("wrote: %s\n", traces_path.c_str());
    }
  }

  return 0;
}
