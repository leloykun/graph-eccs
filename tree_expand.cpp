#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <filesystem>
#include <cassert>

typedef long long ll;
typedef std::vector<int> vi;

// labeled trees <-> prufer codes
// trees -> integers
// labeled trees -> integers

// Step 1: Generate prufer codes
// Step 2: Convert prufer code -> labeled tree
// Step 3: Remove redundant trees (up to isomorphism)

// isomorphic trees:
// 0 <-> 1 <-> 2
// 2 <-> 0 <-> 1

const int MAXN = 1000000;
// taken from http://oeis.org/A000055
const vi NUM_UNLABELED_TREES = {
  1,1,1,1,2,3,6,11,23,47,106,235,551,1301,3159,
 7741,19320,48629,123867,317955,823065,2144505,
 5623756,14828074,39299897,104636890,279793450,
 751065460
};
int PRIME_START = 400;

int is_prime[MAXN];
std::vector<ll> primes;

void prep_primes() {
  for (int i = 0; i < MAXN; ++i) is_prime[i] = 1;
  is_prime[0] = is_prime[1] = 0;
  primes.push_back(2);
  for (int i = 4; i < MAXN; i += 2) is_prime[i] = 0;
  for (ll i = 3; i < MAXN; i += 2)
    if (is_prime[i]) {
      primes.push_back({i});
      for (ll j = i*i; j < MAXN; j += i)
        is_prime[j] = 0;
    }
}

struct tree {
  int n;
  vi par;
  std::vector<vi> adj;
  tree() {}
  tree(int n) : n(n) {
    par = vi(n);
    adj = std::vector<vi>(n, vi());
  }
  void expand() {
    n += 1;
    par.push_back(-1);
    adj.push_back({});
  }
  void add_edge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  /*
  void print(std::ofstream &file) {
    file << n << "\n";
    for (int u = 0; u < n; ++u)
      for (int v : adj[u])
        if (u < v)
          file << u << " " << v << "\n";
  }
  void print_prufer(std::ofstream &file) {
    vi code = to_prufer();
    for (int u : code)
      file << u << " ";
    file << "\n";
  }
  */
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
  vi to_prufer() {
    dfs(n-1, -1);

    int ptr = -1;
    vi deg(n);
    for (int u = 0; u < n; ++u) {
      deg[u] = adj[u].size();
      if (deg[u] == 1 and ptr == -1)
        ptr = u;
    }

    vi code(n-2);
    int leaf = ptr;
    for (int i = 0; i < n-2; ++i) {
      int next = par[leaf];
      code[i] = next;
      if (--deg[next] == 1 and next < ptr)
        leaf = next;
      else {
        ptr++;
        while (deg[ptr] != 1) ptr++;
        leaf = ptr;
      }
    }
    return code;
  }
  void dfs(int u, int p) {
    par[u] = p;
    for (int v : adj[u])
      if (v != p) {
        par[v] = u;
        dfs(v, u);
      }
  }
  int bfs(int u) {
    for (int i = 0; i < n; ++i) par[i] = -1;
    std::queue<int> q; q.push(u); par[u] = u;
    while (!q.empty()) {
      u = q.front();  q.pop();
      for (int v : adj[u]) {
        if (par[v] == -1) {
          par[v] = u;
          q.push(v);
        }
      }
    }
    return u;
  }
  vi tree_centers(int r) {
    vi path;
    int u = bfs(bfs(r));
    for (; u != par[u]; u = par[u])
      path.push_back(u);
    path.push_back(u);

    vi med = {path[path.size()/2]};
    if (path.size() % 2 == 0) med.push_back(path[path.size()/2-1]);
    return med;
  }
  ll rootcode(int u, int p=-1, int d=400) {
    std::vector<ll> k; int nd = (d+1) % primes.size();
    for (int v : adj[u])
      if (v != p)
        k.push_back(rootcode(v, u, nd));
    std::sort(k.begin(), k.end());
    ll h = k.size() + 1;
    for (int i = 0; i < k.size(); ++i)
      h = h*primes[d] + k[i];
    return h;
  }
  ll treecode(int r=0) {
    vi c = tree_centers(r);
    if (c.size() == 1)
      return (rootcode(c[0]) << 1) | 1;
    return (rootcode(c[0])*rootcode(c[1])) << 1;
  }
  bool isomorphic(tree &other) {
    return treecode() == other.treecode();
  }
};

//std::vector<tree> past_trees;
//std::vector<tree> new_trees;
std::vector<vi> past_trees_prufer;
std::vector<vi> new_trees_prufer;
std::set<ll> tree_hashes;

void load_prev_trees(std::string prev_filename) {
  std::ifstream tree_file(prev_filename);

  int N, M;  tree_file >> N >> M;
  for (int i = 0; i < M; ++i) {
    vi code(N-2);
    for (int j = 0; j < N-2; ++j)
      tree_file >> code[j];
    past_trees_prufer.push_back(code);
  }

  /*
  std::cout << "checking for isomorphism\n";
  for (int i = 0; i < M; ++i) {
    tree t1(N);
    t1.from_prufer(past_trees_prufer[i]);
    for (int j = i+1; j < M; ++j) {
      tree t2(N);
      t2.from_prufer(past_trees_prufer[j]);
      assert(!t1.isomorphic(t2));
    }
  }
  std::cout << "done check\n";
  */
}

int main(int argc, char * argv[]) {
  int START_N = std::stoi(argv[1]);
  int END_N = std::stoi(argv[2]);

  if (argc >= 4) {
    PRIME_START = std::stoi(argv[3]);
  }

  std::string prev_filename = "trees/tree_prufer_" + std::to_string(START_N-1) + ".txt";
  if (START_N > 3 and !std::filesystem::exists(prev_filename)) {
    std::cout << prev_filename << " does not exist. Can't continue.\n";
    return 0;
  } else if (START_N < 3) {
    std::cout << "Start node count is less than 1. Can't continue.\n";
    return 0;
  } else if (START_N > 3) {
    std::cout << "Loading past trees...\n";
    load_prev_trees(prev_filename);
    std::cout << "Loaded past trees...\n";
  }

  prep_primes();

  //std::cout << "ready!\n";

  for (int n = START_N; n <= END_N; ++n) {
    new_trees_prufer.clear();
    tree_hashes.clear();
    if (n == 3) {
      tree t(3);
      t.add_edge(0, 1);
      t.add_edge(1, 2);
      new_trees_prufer = {t.to_prufer()};
    } else {
      for (auto p : past_trees_prufer) {
        tree t(n-1);
        t.from_prufer(p);
        t.expand();
        for (int u = 0; u < n-1; ++u) {
          t.add_edge(u, n-1);

          vi code = t.to_prufer();
          tree g(n);
          g.from_prufer(code);

          if (!g.isomorphic(t)) {
            vi code_ = g.to_prufer();
            for (int u : code)
              std::cout << u << " ";
            std::cout << "\n";
            for (int u : code_)
              std::cout << u << " ";
            std::cout << "\n";
          }
          assert(g.isomorphic(t));

          ll h = t.treecode();
          if (tree_hashes.find(h) == tree_hashes.end()) {
            new_trees_prufer.push_back(t.to_prufer());
            tree_hashes.insert(h);
          }

          t.adj[u].pop_back();
          t.adj[n-1].pop_back();
        }
      }
    }

    std::string filename = "trees/tree_prufer_" + std::to_string(n) + ".txt";
    std::ofstream tree_file(filename);

    tree_file << n << " " << NUM_UNLABELED_TREES[n] << "\n";
    for (auto p: new_trees_prufer) {
      for (int u : p)
        tree_file << u << " ";
      tree_file << "\n";
    }
    tree_file.close();

    std::cout << n << " | " << filename << " | " << new_trees_prufer.size() << " " << NUM_UNLABELED_TREES[n] << " | " <<  "\n";
    past_trees_prufer = new_trees_prufer;
  }

  return 0;
}
