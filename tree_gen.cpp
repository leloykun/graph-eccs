#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <string>

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
  1, 1, 1, 1, 2, 3,
  6, 11, 23, 47, 106,
  235, 551, 1301, 3159, 7741
};

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
  vi deg_, par;
  std::vector<vi> adj;
  tree() {}
  tree(int n) : n(n) {
    deg_ = vi(n);
    par = vi(n);
    adj = std::vector<vi>(n, vi());
  }/*
  ~tree() {
    delete deg_;
    delete par;
    delete adj;
  }*/
  void from_prufer(vi &code) {
    for (int u = 0; u < n; ++u) deg_[u] = 1;
    for (int u : code) deg_[u]++;

    int ptr = 0;
    while (deg_[ptr] != 1) ptr++;
    int leaf = ptr;

    for (int u : code) {
      add_edge(leaf, u);
      if (--deg_[u] == 1 and u < ptr)
        leaf = u;
      else {
        ptr++;
        while (deg_[ptr] != 1) ptr++;
        leaf = ptr;
      }
    }
    add_edge(leaf, n-1);
  }
  void add_edge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  void print() {
    std::cout << n << "\n";
    for (int u = 0; u < n; ++u)
      for (int v : adj[u])
        if (u < v)
          std::cout << u << " " << v << "\n";
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
  ll rootcode(int u, int p=-1, int d=200) {
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

vi cur_code;
std::set<ll> tree_hashes;
//std::map<ll, vi > uniq_trees;

// Theorem:
// There are N^(N-2) labeled trees of size N
//
// for each labeled tree of size N, there is
// a prufer code of size N-2

void generate_labeled_trees(int n) {
  if (cur_code.size() == n-2) {
    tree g(n);
    g.from_prufer(cur_code);
    ll h = g.treecode();
    if (tree_hashes.find(h) == tree_hashes.end()) {
      tree_hashes.insert(h);
      g.print();
    }
  } else {
    for (int u = 0; u < n; ++u) {
      cur_code.push_back(u);
      generate_labeled_trees(n);
      cur_code.pop_back();
    }
  }
}

int main(int argc, char * argv[]) {
  prep_primes();

  int N = std::stoi(argv[1]);

  std::cout << NUM_UNLABELED_TREES[N] << "\n";

  cur_code.clear();
  generate_labeled_trees(N);

  /*std::cout << uniq_trees.size() << "\n";
  std::vector<int> ecc(N);
  for (auto &[h, c] : uniq_trees) {
    tree g(N);
    g.from_prufer(c);
    g.print();
  }*/

  return 0;
}
