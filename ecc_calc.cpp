#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <cassert>
#include <map>
#include <tuple>

typedef long long ll;
typedef std::vector<int> vi;
typedef std::vector<vi> vvi;
typedef std::tuple<int, std::string, int> isi;

const vi NUM_UNLABELED_TREES = {
  1,1,1,1,2,3,6,11,23,47,106,235,551,1301,3159,
 7741,19320,48629,123867,317955,823065,2144505,
 5623756,14828074,39299897,104636890,279793450,
 751065460
};

struct tree {
  int n;
  vi dist, par, ecc;
  vvi adj;
  tree(int n) : n(n) {
    dist = vi(n);
    par  = vi(n);
    ecc  = vi(n, -1);
    adj  = vvi(n, vi());
  }
  void from_prufer(vi &code) {
    assert(code.size() + 2 == n);
    for (int u = 0; u < n; ++u) adj[u].clear();
    vi deg(n, 1);
    for (int u : code) deg[u]++;

    int ptr = 0;
    while (deg[ptr] != 1) ptr++;
    int leaf = ptr;

    for (int u : code) {
      add_edge(leaf, u);
      if (--deg[u] == 1 and u < ptr)
        leaf = u;
      else {
        ptr++;
        while (deg[ptr] != 1) ptr++;
        leaf = ptr;
      }
    }
    add_edge(leaf, n-1);
  }
  void add_edge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  int bfs(int u) {
    for (int i = 0; i < n; ++i) par[i] = -1;
    std::queue<int> q; q.push(u);
    par[u] = u; dist[u] = 0;
    while (!q.empty()) {
      u = q.front();  q.pop();
      for (int v : adj[u]) {
        if (par[v] == -1) {
          par[v] = u; dist[v] = dist[u] + 1;
          q.push(v);
        }
      }
    }
    return u;
  }
  void calc_eccs(int r=0) {
    vi path;
    int a = bfs(r);
    int b = bfs(a);
    for (int u = b; u != par[u]; u = par[u])
      path.push_back(u);
    path.push_back(a);

    std::queue<int> q;
    for (int up : path) {
      ecc[up] = std::max(dist[up], int(path.size())-1-dist[up]);
      q.push(up);
    }

    while (!q.empty()) {
      int u = q.front();  q.pop();
      for (int v : adj[u]) {
        if (ecc[v] == -1) {
          ecc[v] = ecc[u] + 1;
          q.push(v);
        }
      }
    }
  }
};

int main(int argc, char * argv[]) {
  int START_N = std::stoi(argv[1]);
  int END_N = std::stoi(argv[2]);

  for (int N = START_N; N <= END_N; ++N) {
    std::string tree_filename = "trees/tree_prufer_" + std::to_string(N) + ".txt";
    std::ifstream tree_file(tree_filename);

    std::string ecc_filename = "trees/ecc_" + std::to_string(N) + ".txt";
    std::ofstream ecc_file(ecc_filename);

    if (!std::filesystem::exists(tree_filename)) {
      std::cout << tree_filename << " does not exist. Continuing...";
      continue;
    }

    int N_, M; tree_file >> N_ >> M;
    assert(N_ == N);

    std::vector<isi> encs;
    for (int m = 0; m < M; ++m) {
      vi code(N-2);
      for (int j = 0; j < N-2; ++j)
        tree_file >> code[j];

      tree t(N);
      t.from_prufer(code);
      t.calc_eccs();

      ll ecc_sum = 0;
      std::sort(t.ecc.begin(), t.ecc.end(), std::greater<int>());
      std::map<int, int> cnt;
      for (int e : t.ecc) {
        ecc_sum += e;
        cnt[e]++;
      }

      int BI;
      for (int i = 0; i < N; ++i)
        if (t.ecc[i]*N < ecc_sum) {
          BI = i;
          break;
        }

      std::string enc = "";
      for (auto &[e, c] : cnt)
        if (c == 1)
          enc += std::to_string(e) + " ";
        else
          enc += std::to_string(e) + "x" + std::to_string(c) + " ";
      enc.pop_back();

      encs.push_back({BI, enc, m});
    }

    std::sort(encs.begin(), encs.end());

    ecc_file << N << " " << M << "\n";
    for (auto &[bi, enc, i] : encs)
      ecc_file << bi << " " << enc << " " << i << "\n";

    tree_file.close();
    ecc_file.close();
    std::cout << N << " | " << M << " " << NUM_UNLABELED_TREES[N] << " | " <<  "\n";
  }

  return 0;
}
