# Sumário

- [**Estruturas de Dados**](#estruturas-de-dados)
  - [Segment Tree](#segment-tree)
    - [Point update](#point-update)
    - [Lazy propagation](#lazy-propagation)
  - [Fenwick Tree / BIT](#bit)
    - [BIT 2D](#bit-2d)
  - [Sparse Table](#sparse-table)
  - [SQRT Decomposition](#sqrt-decomposition)
  - [Trie](#trie)
  - [Union Find](#union-find)
  - [Ordered Set](#ordered-set)
- [**Grafos**](#grafos)
  - [Dijkstra](#dijkstra)
  - [Bellman-Ford](#bellman-ford)
  - [Floyd-Warshall](#floyd-warshall)
  - [Pontes e Pontos de Articulação](#pontes-e-pontos-de-articulação)
  - [Componentes Fortemente Conexos](#componentes-fortemente-conexos)
  - [Ordenação Topológica](#ordenação-topológica)
  - [Minimum Spanning Tree](#minimum-spanning-tree)
  - [Lowest Common Ancestor](#lowest-common-ancestor)
  - [Bipartite Matching](#bipartite-matching)
  - [Fluxo Máximo](#fluxo-máximo)
    - [Edmons-Karp](#edmons-karp)
- [**Strings**](#strings)
  - [KMP](#kmp)
  - [Z Function](#z-function)
  - [Suffix Array](#suffix-array)
- [**Matemática**](#matemática)
  - [MDC e MMC](#mdc-e-mmc)
  - [Euclides Extendido](#euclides-extendido)
  - [Inverso Multiplicativo](#inverso-multiplicativo)
  - [Crivo de Erastóstenes](#crivo-de-erastóstenes)
    - [Crivo segmentado](#crivo-segmentado)
  - [Totiente de Euler](#totiente-de-euler)
  - [Exponenciação de Matrizes](#exponenciação-de-matrizes)
  - [Miller-Rabin + Pollard's Rho](#miller-rabin--pollards-rho)
- [**Geometria Computacional**](#geometria-computacional)
  - [Interseção de Retas](#interseção-de-retas)
  - [Outras Operações com Ponto e Reta](#outras-operações-com-ponto-e-reta)
  - [Operações com Vetores](#operações-com-vetores)
  - [Interseção de Segmentos](#interseção-de-segmentos)
  - [Área de Polígono](#área-de-polígono)
  - [Convex Hull](#convex-hull)
  - [Par de Pontos Mais Próximos](#par-de-pontos-mais-próximos)
  - [Ponto Dentro do Polígono](#ponto-dentro-do-polígono)


# Estruturas de dados

### Segment Tree

#### Point update

https://leetcode.com/problems/range-sum-query-mutable/description/

Range Sum Query

```c
struct SegmentTree {
    
    vector<int> tree;
    int s, e;

    SegmentTree(int s, int e, int val = 0): s(s), e(e) {
        int sz = e - s + 1;
        tree = vector<int>(4 * sz, val);
    }
    
    SegmentTree(vector<int> &v) {
        *this = SegmentTree(0, v.size() - 1);
        build(v, 1, s, e);
    }
    
    void build(vector<int> &v, int i, int a, int b) {
        
        if(a == b)
            tree[i] = v[a];
        else {
            int m = (a + b) >> 1;
            
            build(v, 2*i, a, m);
            build(v, 2*i + 1, m+1, b);
            
            tree[i] = tree[2*i] + tree[2*i + 1];
        }
    }
    
    void update(int p, int val, int i, int a, int b) {
        
        if(a > p || b < p) return;
        
        if(a == b && p == a)
            tree[i] = val;
        else {
            int m = (a + b) >> 1;
            
            update(p, val, 2*i, a, m);
            update(p, val, 2*i + 1, m+1, b);
            
            tree[i] = tree[2*i] + tree[2*i + 1];
        }
    }
    
    int query(int A, int B, int i, int a, int b) {
        
        if(a > B || b < A) return 0;
        
        if(a >= A && b <= B) return tree[i];
        
        int m = (a + b) >> 1;
        return query(A, B, 2*i, a, m) + query(A, B, 2*i + 1, m+1, b);
    }
    
    void update(int p, int val) {
        update(p, val, 1, s, e);
    }
    
    int query(int A, int B) {
        return query(A, B, 1, s, e);
    }
};
```

#### Lazy propagation

https://www.urionlinejudge.com.br/judge/pt/problems/view/1500

Range Sum Query

```c
struct SegmentTree {
    
    vector<ll> tree, lazy;
    int s, e;

    SegmentTree(int s, int e): s(s), e(e) {
        int sz = e - s + 1;
        tree = lazy = vector<ll>(4*sz, 0);
    }

    SegmentTree(vector<ll> &v) {
        *this = SegmentTree(0, v.size() - 1);
        build(v, 1, s, e);
    }

    void build(vector<ll> &v, int i, int a, int b) {
        lazy[i] = 0;
        if(a == b)
            tree[i] = v[i];
        else {
            int m = (a + b) >> 1;
            build(v, 2*i, a, m);
            build(v, 2*i + 1, m+1, b);
            tree[i] = tree[2*i] + tree[2*i + 1];
        }
    }

    void propagate(int i, int a, int b) {

        if(!lazy[i]) return;

        tree[i] += (b - a + 1) * lazy[i];
        if(a != b) {
            lazy[2*i] += lazy[i];
            lazy[2*i + 1] += lazy[i];
        }
        lazy[i] = 0;
    }

    void update(int A, int B, ll val, int i, int a, int b) {

        propagate(i, a, b);
        if(a > B || b < A) return;

        if(a >= A && b <= B) {
            tree[i] += (b - a + 1) * val;
            if(a != b) {
                lazy[2*i] += val;
                lazy[2*i + 1] += val;
            }
        }
        else {
            int m = (a + b) >> 1;
            update(A, B, val, 2*i, a, m);
            update(A, B, val, 2*i + 1, m+1, b);
            tree[i] = tree[2*i] + tree[2*i + 1];
        }
    }

    ll query(int A, int B, int i, int a, int b) {
        
        if(a > B || b < A) return 0;
        propagate(i, a, b);

        if(a >= A && b <= B) return tree[i];

        int m = (a + b) >> 1;
        return query(A, B, 2*i, a, m) + query(A, B, 2*i + 1, m+1, b);
    }

    void update(int A, int B, ll val) {
        update(A, B, val, 1, s, e);
    }

    ll query(int A, int B) {
        return query(A, B, 1, s, e);
    }
};
```

### BIT

https://www.urionlinejudge.com.br/judge/pt/problems/view/2857

```c
struct Bit { //1-indexado

    int n;
    vector<int> arr;

    int lsone(int x) {
        return x & -x;
    }

    Bit(int N, int val = 0): n(N + 1) {
        arr = vector<int>(N + 1, val);
    }

    Bit(vector<int> &v) {

        *this = Bit(v.size());
        for (int i = 1; i <= n; ++i)
            update(v[i], i);
    }

    void update(int pos, int val) {

        for (; pos < n; pos += lsone(pos))
            arr[pos] += val;
    }

    int get(int pos) {

        int sum = 0;
        for (; pos > 0; pos -= lsone(pos))
            sum += arr[pos];
        return sum;
    }

    int get(int a, int b) {
        return get (b) - get(a - 1);
    }
};
```

#### BIT 2D

https://www.urionlinejudge.com.br/judge/pt/problems/view/1112

```c
int bit[MAX][MAX];

void update(int x, int y, int val) {
    for(int i = x; i < MAX; i += i & -i)
        for(int j = y; j < MAX; j += j & -j)
            bit[i][j] += val;
}

int get(int x, int y) {
    int ans = 0;
    for(int i = x; i > 0; i -= i & -i)
        for(int j = y; j > 0; j -= j & -j)
            ans += bit[i][j];
    return ans;
}

int get(int x1, int y1, int x2, int y2) {
    if(x1 > x2) swap(x1, x2);
    if(y1 > y2) swap(y1, y2);
    return get(x2, y2) - get(x1-1, y2) - get(x2, y1-1) + get(x1-1, y1-1);
}
```

### Sparse Table

https://www.spoj.com/problems/RMQSQ/

Range Minimum Query

```c
#define MAX 100100
#define LOG_MAX 20

int v[MAX];
int table[MAX][LOG_MAX];

void build(int n) {
    for(int i = 0; i < n; ++i)
        table[i][0] = v[i];
    for(int i = 1; i < LOG_MAX; ++i)
        for(int j = 0; j < MAX; ++j)
            table[j][i] = min(table[j][i-1], table[min(n-1, j + (1 << (i-1)))][i-1]);
}

int get_min(int i, int j) {
    int d = log2(j - i + 1); 
    return min(table[i][d], table[j - (1 << d) + 1][d]);
}
```

### SQRT Decomposition

https://www.urionlinejudge.com.br/judge/pt/problems/view/2800

```c
const int MAXR = 350; //raiz sempre fixa
const int MAXN = 100010;

int bucket[MAXR][MAXN];
int arr[MAXN];

// Altera o valor da posição id para x
void update(int id, int x){
    int bloco = id / MAXR;
    bucket[bloco][arr[id]]--;
    bucket[bloco][x]++;
    arr[id] = x;
}
// A query retorna o número de valores iguais a W no intervalo A - B
int query(int A, int B, int W){
    int b1 = A/MAXR;
    int b2 = B/MAXR;
    int cnt = 0;
    if(b1 == b2){
        for(int i = A; i <= B; i++){
            if(arr[i] == W) cnt++;
        }
    }
    else{
        for(int i = b1+1; i <= b2-1; i++)
            cnt += bucket[i][W];
        for(int i = A; i < (b1+1) * MAXR; i++)
            if(arr[i] == W) cnt++;
        for(int i = b2 * MAXR; i <= B; i++)
            if(arr[i] == W) cnt++;
    }
    return cnt;
}
```

### Trie

https://www.spoj.com/problems/STRMATCH/

Versão recursiva

```c
struct Trie {

    int cnt = 0; // numero de prefixos que terminam neste no
    Trie *c[26];

    Trie() {
        memset(c, 0, sizeof c);
    }

    // insere s[i..] na trie
    void insert(string &s, int i = 0) {
        
        ++cnt;
        if(i >= s.length()) return;
        if(c[s[i] - 'a'] == NULL)
            c[s[i] - 'a'] = new Trie();
        c[s[i] - 'a']->insert(s, i+1);
    }

    // retorna o numero de prefixos iguais a s
    int count(string &s, int i = 0) {
        
        if(i == s.length()) return cnt;
        if(c[s[i] - 'a'] == NULL) return 0;
        return c[s[i] - 'a']->count(s, i+1);
    }
};
```

Versão iterativa

```c
struct Node{
    Node *children[26];
    int isEnd;
    Node(){
        for(int i = 0; i < 26; i++){
            children[i] = NULL;
        }
        isEnd = 0;
    }
};

struct Trie{
    Node *root = new Node();
    void insert(const string &s, int ini) {
        Node *it = root;
        for(int i = ini; i < (int)s.size(); i++){
            char c = s[i];
            if(!it->children[c-'a']){
                it->children[c-'a'] = new Node();
            }
            it = it->children[c-'a'];
            it->isEnd++;
        } 
    }
    int search(const string &s){
        Node *it = root;
        int vezes=0;
        for(auto c: s){
            if(it->children[c-'a'] == NULL)
                return 0;
            it = it->children[c-'a'];
        }
        return it->isEnd;
    }
};
```

### Union Find

https://www.hackerearth.com/practice/data-structures/disjoint-data-strutures/basics-of-disjoint-data-structures/practice-problems/algorithm/count-friends/

```c

typedef vector<int> vi;

struct UnionFind {
    
    vi p, sz;
    int n;
    
    UnionFind(int n): n(n), p(vi(n, 0)), sz(vi(n, 1)) {
        for(int i = 0; i < n; ++i)
            p[i] = i;
    }
    
    int find(int v) {
        return p[v] == v ? v : v = find(p[v]);
    }
    
    void uni(int u, int v) {
        u = find(u), v = find(v);
        if(u == v) return;
        if(sz[v] > sz[u]) swap(u, v);
        p[v] = u, sz[u] += sz[v];
    }
    
    int size(int v) {
        return sz[find(v)];
    }
};
```

### Ordered Set

http://codeforces.com/blog/entry/11080?locale=en

```c
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;

typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;

// st.find_by_order(k) - iterador para o k-ésimo menor elemento (0-indexado)
// st.order_by_key(x) - número de elementos menores que x
```

# Grafos

### Dijkstra

https://www.hackerrank.com/challenges/dijkstrashortreach/problem

```c
int dist[MAX];

int dijkstra(int orig, int dest) {

    memset(dist, 0x3f, sizeof dist);
    dist[orig] = 0;
    priority_queue<pii> pq;
    pq.push(pii(0, orig));
    while(!pq.empty()) {
        
        auto p = pq.top(); pq.pop();
        int d = -p.first;
        int u = p.second;
        
        if(u == dest) return d;
        if(d > dist[u]) continue;
        
        for(auto v : g[u])
            if(d + v.second < dist[v.first]) {
                dist[v.first] = d + v.second;
                pq.push(pii(-dist[v.first], v.first));
            }
    }
    return -1;
}
```

### Bellman-Ford

https://practice.geeksforgeeks.org/problems/negative-weight-cycle/0

```c
#define INF 0x3f3f3f3f

struct Edge { int u, v, w; } edges[MAXM];

int n, m, d[MAXN];

// retorna true se não houver ciclo negativo
bool bellman_ford(int s) {

    memset(d, 0x3f, sizeof d);
    d[s] = 0;
    bool ok = true; // indica se relaxou alguma aresta
    for(int k = 0; k < n-1 && ok; ++k) {
        ok = false;
        for(int i = 0; i < m; ++i) {
            auto e = edges[i];
            if(d[e.v] > d[e.u] + e.w)
                d[e.v] = d[e.u] + e.w, ok = true;
        }
    }

    for(int i = 0; i < m; ++i) {
        auto e = edges[i];
        if(d[e.v] > d[e.u] + e.w)
            return false;
    }
    return true;
}
```

### Floyd-Warshall

https://practice.geeksforgeeks.org/problems/implementing-floyd-warshall/0

```c
int g[MAX][MAX]; // g[i][i] = INF
int dist[MAX][MAX];

void floyd(int n) {
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            dist[i][j] = g[i][j];
    for(int k = 0; k < n; ++k)
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j)
                dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
}
```

### Pontes e Pontos de Articulação

```c
const int MAXN = 550;

int d[MAXN], low[MAXN], tempo=0, raiz;
/*
d[u] é o tempo de descoberta de u
low[u] é o menor d de um descendente próprio de u
*/
vector<int> adj[MAXN];
vector<int> vertices_corte;
vector<pair<int, int> > pontes;

void dfs(int u, int p){

    int nf = 0;

    bool any = false;

    d[u] = low[u] = ++tempo;
    for(int v: adj[u]){

        if(!d[v]){
            dfs(v, u);
            nf++;

            if(low[v] >= d[u]) any = true;

            low[u] = min(low[u], low[v]);

            if(low[v] > d[u]){
                // u-v é uma ponte
                pontes.push_back({v, u});
            }

        }
        else if(v != p){
            low[u] = min(low[u], d[v]);
        }
    }
    if( (u == raiz && nf >= 2) || (u != raiz && any) ){
        // u é um vértice de corte
        vertices_corte.push_back(u);
    }
}

```

### Componentes Fortemente Conexos

```c
int  d[MAX], low[MAX], tempo, pilha[MAX], topo=-1, cont=0, comp[MAX];
lli memo[MAX];

vector<int> g[MAX];

void dfs(int u) {

    pilha[++topo] = u;
    low[u] = d[u] = ++tempo;

    for(int i = 0; i < (int)g[u].size(); i++) {
        int v = g[u][i];
        if(!d[v]) {
            dfs(v);
            low[u] = min(low[u], low[v]); 
        }
        else
            low[u] = min(low[u], d[v]);
    }

    if(d[u] == low[u]) {
        int x;
        ++cont;
        do {
            x = pilha[topo--];
            comp[x] = cont;
            d[x] = 5*MAX; //equivalente a d[x] = INF
        } while(x != u);
    }
}
```

### Ordenação Topológica

https://olimpiada.ic.unicamp.br/pratique/p2/2011/f2/escalona/

```c
void bfs(vector<int> &ans){

    priority_queue<int> pq;

    for(int i = 0; i < n; i++) 
        if(!grau[i]) pq.push(-i);

    while(!pq.empty()){
        int v = -pq.top(); pq.pop();
        ans.push_back(v);
        for(int u: adj[v]){
            if(grau[u]){
                grau[u]--;
                if(!grau[u]) pq.push(-u);
            }
        }
    }
}
```

### Minimum Spanning Tree

https://www.urionlinejudge.com.br/judge/pt/problems/view/2404

```c
struct Aresta {
    int u, v, peso;
};

bool comp(Aresta a, Aresta b) {
    return a.peso < b.peso;
}

int kruskal(vector<Aresta> &arestas, int n) {

    sort(arestas.begin(), arestas.end(), comp);
    int soma = 0;
    UnionFind uf(n);
    for(Aresta aresta: arestas) {
        int v = uf.find(aresta.v);
        int u = uf.find(aresta.u);
        int peso = aresta.peso;
        if(v != u) {
            soma += peso;
            uf.uni(v, u);
        }
    }
    return soma;
}
```

### Lowest Common Ancestor

https://www.spoj.com/problems/LCA/

```c
#define MAX 100100
#define LOG_MAX 20

vector<int> g[MAX];
int anc[MAX][LOG_MAX], h[MAX];

void dfs(int u, int p, int hu) {
    if(h[u] > -1) return;
    h[u] = hu;
    anc[u][0] = p;
    for(auto v : g[u])
        dfs(v, u, hu + 1);
}

void build(int n) {
    memset(h, -1, sizeof h);
    dfs(1, 0, 0);
    for(int i = 1; i < LOG_MAX; ++i)
        for(int v = 1; v <= n; ++v)
            anc[v][i] = anc[ anc[v][i-1] ][i-1];
}

int lca(int u, int v) {

    if(h[u] > h[v]) swap(u, v);

    for(int i = LOG_MAX-1; i >= 0; --i)
        if(h[v] - h[u] >= (1 << i))
            v = anc[v][i];
    
    if(u == v) return u;

    for(int i = LOG_MAX-1; i >= 0; --i)
        if(anc[u][i] != anc[v][i])
            u = anc[u][i], v = anc[v][i];

    return anc[u][0];
}
```

### Bipartite Matching

### Fluxo Máximo

#### Edmons-Karp

https://practice.geeksforgeeks.org/problems/find-the-maximum-flow/0

```c
int n, g[MAXN][MAXN], gr[MAXN][MAXN], f[MAXN][MAXN], p[MAXN];

void bfs(int s, int t) {
    
    queue<int> q;
    memset(p, -1, sizeof p);
    q.push(s);
    while(!q.empty() && p[t] == -1) {
        int u = q.front(); q.pop();
        for(int v = 0; v < n; ++v)
            if(gr[u][v] && p[v] == -1)
                q.push(v), p[v] = u;
    }
}

int edmons_karp(int s, int t) {
    
    memset(f, 0, sizeof f);
    memcpy(gr, g, sizeof gr);
    
    for(bfs(s, t); p[t] != -1; bfs(s, t)) {
        
        int mn = INT_MAX;
        for(int v = t; v != s; v = p[v])
            mn = min(mn, gr[p[v]][v]);
        
        for(int v = t; v != s; v = p[v]) {
            int u = p[v];
            gr[u][v] -= mn;
            gr[v][u] += mn;
            if(g[u][v] > 0) f[u][v] += mn;
            else f[v][u] -= mn;
        }
    }
    
    int ans = 0;
    for(int i = 0; i < n; ++i)
        ans += f[s][i] - f[i][s];
    return ans;
}
```

# Strings

### KMP

### Z Function

https://practice.geeksforgeeks.org/problems/search-pattern/0

```c
void zf(string &s, vector<int> &z) {
    int L = 0, R = 0, n = s.length();    
    z.resize(n);
    z[0] = n;
    for(int i = 1; i < n; ++i) {
        if(R < i) {
            for(L = R = i; R < n && s[R] == s[R - L]; ++R);
            z[i] = R - L;
            --R;
        }
        else if(z[i - L] <= R - i)
            z[i] = z[i - L];
        else {
            for(L = i, ++R; R < n && s[R] == s[R - L]; ++R);
            z[i] = R - L;
            --R;
        }
    }
}
```

### Suffix Array

# Matemática

### MDC e MMC

https://practice.geeksforgeeks.org/problems/lcm-and-gcd/0

```c
// __gcd(a, b)
int mdc(int a, int b) {
    return b ? mdc(b, a % b) : a;
}

int mmc(int a, int b) {
    return a / mdc(a, b) * b;
}
```

### Euclides Extendido

```c
// x * a + y * b = mdc(a, b)
int mdc(int a, int b, int &x, int &y) {
    
    if(b == 0) {
        x = 1, y = 0;
        return a;
    }

    int x2, y2;
    int m = mdc(b, a % b, x2, y2);
    x = y2;
    y = x2 - (a / b) * y2;
    return m;
}
```

### Inverso Multiplicativo

```c
// m primo => invMod(a, m) = a^(m-2) (mod m)
int invMod(int a, int m) {
    int x, y;
    if(mdc(a, m, x, y) != 1) return -1; // não existe inverso    
    return (x % m + m) % m;
}
```

### Crivo de Erastóstenes

https://practice.geeksforgeeks.org/problems/sieve-of-eratosthenes/0/

```c
int p[MAX], np;
bool isPrime[MAX];

void crivo() {
    p[np++] = 2;
    isPrime[2] = false;
    for (int i = 3; i < MAX; i += 2)
        isPrime[i] = true;
    for (int i = 3; i < MAX; i += 2)
        if (isPrime[i]) {
            p[np++] = i;
            for (int j = 2*i; j < MAX; j += i)
                isPrime[j] = false;
        }
}
```

#### Crivo segmentado

https://www.spoj.com/problems/PRINT/

```c
// Primeiro gerar crivo até sqrt(maior valor)
// Adiciona todos os primos no intervalo [l, u] em ans
void primes(int l, int u, vector<int> &ans) {
    
    vector<bool> prime(u - l + 1, true);
    if (l < 1) prime[0] = false;
    if (l < 2) prime[1 - l] = false;

    for (int i = 0; i < np && p[i] <= u; ++i) {
        int start = (l > 1 ? p[i] * (l / p[i]) : 2*p[i]);
        while (start < l || start == p[i]) start += p[i];
        for (int j = start; j <= u; j += p[i])
            prime[j - l] = false;
    }

    for (int i = 0; i < prime.size(); ++i)
        if (prime[i])
            ans.push_back(i + l);
}
```

### Totiente de Euler

https://practice.geeksforgeeks.org/problems/euler-totient-function/0/

Para um único número

```c
int phi(int n) {
    int ans = 1;
    for(int i = 0; i < np && p[i]*p[i] <= n; ++i) {
        if(n % p[i] == 0) {
            ans *= p[i] - 1;
            for(n /= p[i]; n % p[i] == 0; n /= p[i])
                ans *= p[i];
        }
    }
    if(n > 1) ans *= n-1;
    return ans;
}
```

Para um intervalo

```c
int phi[MAX];

void build_phi() {
    for(int i = 1; i < MAX; ++i)
        phi[i] = i;
    for(int i = 2; i < MAX; ++i)
        if(phi[i] == i) {
            phi[i] = i - 1;
            for(int j = 2*i; j < MAX; j += i)
                phi[j] = (phi[j] / i) * (i - 1);
        }
}
```

### Exponenciação de Matrizes

### Miller-Rabin + Pollard's Rho

https://uva.onlinejudge.org/index.php?option=onlinejudge&page=show_problem&problem=1333

```c
ll randll() {
    ll t = rand();
    return (t << 31) | rand();
}

ll mdc(ll a, ll b) {
    return b ? mdc(b, a % b) : a;
}

ll sum(ll a, ll b, ll m) {
    a += b;
    if(a > m) a -= m;
    return a;
}

ll mul(ll a, ll b, ll m) {
    ll ans = 0;
    while (b) {
        if (b & 1) ans = sum(ans, a, m);
        a = sum(a, a, m);
        b >>= 1;
    }
    return ans;
}

ll exp(ll a, ll b, ll m) {
    ll ans = 1;
    while (b) {
        if (b & 1) ans = mul(ans, a, m);
        a = mul(a, a, m);
        b >>= 1;
    }
    return ans;
}

// Miller-Rabin
bool isPrime(ll n) {
    
    // testes suficientes para garantir corretude para n <= 10^18
    static vector<int> a = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    
    if (n <= 1) return false;
    if (n <= 3) return true;

    ll s = 0, d = n-1;
    while(d % 2 == 0)
        ++s, d >>= 1;

    for(int i = 0; i < a.size() && a[i] < n; ++i) {
        
        ll x = exp(a[i], d, n);
        if(x == 1 || x == n-1)
            continue;
        for(int r = 1; r < s; ++r) {
            x = mul(x, x, n);
            if(x == 1) return false;
            if(x == n-1) break;
        }
        if(x != n-1) return false;
    }
    return true;
}

// usar em numeros impares
ll pollardRho(ll n) {

    ll x = (randll() % (n-1)) + 1;
    ll c = (randll() % (n-1)) + 1;
    ll y = x, d = 1;
    while (d == 1) {
        x = sum(mul(x, x, n), c, n);
        y = sum(mul(y, y, n), c, n);
        y = sum(mul(y, y, n), c, n);

        d = mdc(abs(x - y), n);
        if(d == n) return pollardRho(n);
    }
    return d;
}

// adiciona os fatores primos de n em ans (fora de ordem)
void factors(ll n, vector<ll> &ans) {
    
    if(n <= 1) return;

    if(isPrime(n))
        ans.push_back(n);
    else {
        ll f = (n % 2 == 0 ? 2 : pollardRho(n));
        factors(f, ans);
        factors(n / f, ans);
    }
}
```

# Geometria Computacional

### Interseção de Retas

https://uva.onlinejudge.org/index.php?option=onlinejudge&page=show_problem&problem=314

```c
#define EPS 1e-9
#define same(x, y) (fabs(x - y) < EPS)
#define inRange(a, b, c) (c >= fmin(a, b) - EPS && c <= fmax(a, b) + EPS)

struct Point {
    
    double x, y;

    Point(double x = 0, double y = 0): x(x), y(y) {}
};

struct Line { double a, b, c; };

Line pointsToLine(Point p1, Point p2) {

    Line l;
    if (same(p1.x, p2.x)) {
        l.a = 1.0;
        l.b = 0.0;
        l.c = -p1.x;
    } else {
        l.a = -(p1.y - p2.y) / (p1.x - p2.x);
        l.b = 1.0;
        l.c = -(l.a * p1.x) - p1.y;
    }
    return l;
}

bool parallel(Line l1, Line l2) {
    return same(l1.a, l2.a) && same(l1.b, l2.b);
}

bool sameLine(Line l1, Line l2) {
    return parallel(l1, l2) && same(l1.c, l2.c);
}

bool intersect(Line l1, Line l2) {
    return !parallel(l1, l2);
}

Point pointIntersection(Line l1, Line l2) {

    Point p;
    p.x = (l2.b * l1.c - l1.b * l2.c) / (l2.a * l1.b - l1.a * l2.b);
    if (same(l1.b, 0))
        p.y = -(l2.a * p.x + l2.c);
    else p.y = -(l1.a * p.x + l1.c);
    return p;
}
```


### Outras Operações com Ponto e Reta

```c
//Distancia entre os pontos a e b
double dist(Point a, Point b) {
    return hypot(fabs(b.x - a.x), fabs(b.y - a.y));
}

//Retorna a distancia de p para a reta que contem ab
double distToLine(Point a, Point b, Point p) {

    Vec ap = toVec(a, p), ab = toVec(a, b);
    double u = dot(ap, ab) / norm_sq(ab);
    a = translate(a, scale(ab, u));
    return dist(p, a);
}

double toRad(double t) { return t * M_PI / 180.0; }
double toDeg(double t) { return t * 180.0 / M_PI; }

// Rotaciona p em theta graus anti-horario
Point rotate(Point p, double theta) {

    double rad = toRad(theta);
    return Point(p.x * cos(rad) - p.y * sin(rad),
                 p.x * sin(rad) + p.y * cos(rad));
}
```


### Operações com Vetores

```c
typedef Point Vec;

Vec toVec(Point a, Point b) { return Vec(b.x - a.x, b.y - a.y); }

//Produto escalar
double dot(Vec a, Vec b) { return a.x * b.x + a.y * b.y; }

//Produto vetorial
double cross(Vec a, Vec b) { return a.x * b.y - a.y * b.x; }

Point translate(Point p, Vec v) { return Point(p.x + v.x, p.y + v.y); }

Vec scale(Vec v, double u) { return Vec(v.x * u, v.y * u); }

//Norma(modulo) do vetor ao quadrado
double norm_sq(Vec v) { return v.x * v.x + v.y * v.y; }

//Orientacao anti-horaria
bool ccw(Point a, Point b, Point c) {
    return cross(toVec(a, b), toVec(a, c)) > 0;
}

bool collinear(Point a, Point b, Point c) {
    return same(cross(toVec(a, b), toVec(a, c)), 0);
}
```


### Interseção de Segmentos

https://uva.onlinejudge.org/index.php?option=onlinejudge&page=show_problem&problem=127

```c
struct Segment {

    Point s, e;
    double dist;

    Segment(Point s, Point e): s(s), e(e) {
        dist = hypot(fabs(e.x - s.x), fabs(e.y - s.y));
    }

    Segment(double sx = 0, double sy = 0, double ex = 0, double ey = 0) {
        *this = Segment(Point(sx, sy), Point(ex, ey));
    }

    //Retorna true se p esta no segmento
    //Deve ser usado apos collinear
    bool contains(Point p) {
        return inRange(s.x, e.x, p.x) && inRange(s.y, e.y, p.y);
    }

    bool intersect(Segment seg) {

        bool o1 = ccw(s, e, seg.s);
        bool o2 = ccw(s, e, seg.e);
        bool o3 = ccw(seg.s, seg.e, s);
        bool o4 = ccw(seg.s, seg.e, e);

        return (o1 != o2 && o3 != o4) ||
               (collinear(s, e, seg.s) && contains(seg.s)) ||
               (collinear(s, e, seg.e) && contains(seg.e)) ||
               (collinear(seg.s, seg.e, s) && seg.contains(s)) ||
               (collinear(seg.s, seg.e, e) && seg.contains(e));
    }
};
```


### Área de Polígono

https://www.math10.com/en/geometry/geogebra/fullscreen.html
https://uva.onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&page=show_problem&problem=45

```c
double area(const vector<Point>& p) {

    double a = 0.0;
    for (int i = 1; i < p.size() - 1; ++i)
        a += cross(toVec(p[0], p[i]), toVec(p[0], p[i + 1]));
    return fabs(a * 0.5);
}
```


### Convex Hull

https://practice.geeksforgeeks.org/problems/convex-hull/0
https://uva.onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&page=show_problem&problem=45

```c
Point p0; //Vertice inicial do convex hull

//Ordena pontos no sentido anti-horario
bool cmp(Point p1, Point p2) {

    double ori = cross(toVec(p0, p1), toVec(p0, p2));
    return same(ori, 0) ? dist(p0, p1) < dist(p0, p2):
           atan2(p1.y - p0.y, p1.x - p0.x) < atan2(p2.y - p0.y, p2.x - p0.x);
}

//Retorna vetor de pontos do convex hull
vector<Point> grahamScan(vector<Point>& p) {

    vector<Point> newP;
    int iMin = 0, n = p.size();

    for (int i = 1; i < n; ++i)
        if (p[i].y < p[iMin].y || (p[i].y == p[iMin].y && p[i].x < p[iMin].x))
            iMin = i;
    swap(p[iMin], p[0]);
    p0 = p[0];
    newP.push_back(p0);
    sort(p.begin() + 1, p.end(), cmp);

    for (int i = 1; i < n; ++i) {
        while (i < n - 1 && same(cross(toVec(p0, p[i]), toVec(p0, p[i + 1])), 0))
            i++;
        newP.push_back(p[i]);
    }

    vector<Point> poly(newP.size());
    if (newP.size() > 2) {
        for (int i = 0; i < 3; ++i)
            poly[i] = newP[i];
        int m = 3;
        for (int i = 3; i < newP.size(); ++i) {
            while (cross(toVec(poly[m - 2], poly[m - 1]), toVec(poly[m - 2], newP[i])) < EPS)
                --m;
            poly[m++] = newP[i];
        }
        poly.resize(m);
    }
    return poly;
}
```



### Par de Pontos Mais Próximos

https://www.urionlinejudge.com.br/judge/pt/problems/view/1295

```c
typedef pair<Point, Point> ppp;

bool cmpX(Point p1, Point p2) { return p1.x < p2.x; }

bool cmpY(Point p1, Point p2) { return p1.y < p2.y; }

ppp closestStrip(Point strip[], int m, double minD) {

    ppp res = {{-INF, -INF}, {INF, INF}};
    for (int i = 0; i < m - 1; ++i)
        for (int j = i + 1; j < m && strip[j].y - strip[i].y < minD; ++j)
            if (dist(strip[i], strip[j]) < minD) {
                minD = dist(strip[i], strip[j]);
                res = { strip[i], strip[j] };
            }
    return res;
}

ppp closestByBruteForce(Point vp[], int sz) {

    ppp res = {{-INF, -INF}, {INF, INF}};
    for (int i = 0; i < sz - 1; ++i)
        for (int j = i + 1; j < sz; ++j)
            if (dist(vp[i], vp[j]) < dist(res.first, res.second))
                res = {vp[i], vp[j]};
    return res;
}

ppp closestUtil(Point px[], int sz) {

    if (sz < 4)
        return closestByBruteForce(px, sz);

    int mid = (sz - 1) / 2, l = 0, r = 0;
    Point midP = px[mid];
    Point pxl[sz], pxr[sz];
    ppp res;

    for(int i = 0; i < sz; i++) {
        if(px[i].x < midP.x || (px[i].x == midP.x && r > l))
            pxl[l++] = px[i];
        else
            pxr[r++] = px[i];
    }
    ppp pl = closestUtil(pxl, l);
    ppp pr = closestUtil(pxr, r);
    double plDist = dist(pl.first, pl.second);
    double prDist = dist(pr.first, pr.second);
    double d = fmin(plDist, prDist);
    res = plDist < prDist ? pl : pr;

    sort(px, px + sz, cmpY);

    Point strip[sz];
    int m = 0;
    for (int i = 0; i < sz; ++i)
        if (fabs(px[i].x - midP.x) < d)
            strip[m++] = px[i];

    ppp spDist = closestStrip(strip, m, d);
    if (dist(spDist.first, spDist.second) < d)
        return spDist;
    return res;
}

ppp closestPair(vector<Point>& ps, int n) {

    Point px[n];
    for (int i = 0; i < n; ++i) px[i] = ps[i];
    sort(px, px + n, cmpX);
    return closestUtil(px, n);
}

double minimumDist(vector<Point>& ps) {

    if (ps.size() < 3) return INF;
    ppp close = closestPair(ps, ps.size());
    return dist(close.first, close.second);
}
```

### Ponto Dentro do Polígono

https://uva.onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&page=show_problem&problem=45

```c
//Retorna o angulo aob em radianos
double angle(Point a, Point o, Point b) {
    return acos(dot(toVec(o, a), toVec(o, b)) / (dist(o, a) * dist(o, b)));
}

//Nao considera pontos nos vertices
//como pontos internos
bool inPolygon(vector<Point>& p, Point pt) {

    if (p.size() < 3) return false;
    double sum = 0;
    int j = p.size() - 1;
    for (int i = 0; i < p.size(); ++i) {
        double ang = angle(p[j], pt, p[i]);
        if (ccw(pt, p[j], p[i]))
            sum += ang;
        else sum -= ang;
        j = i;
    }
    return same(fabs(sum) - 2 * M_PI, 0);
}
```
