# Sumário

- [**Estruturas de Dados**](#estruturas-de-dados)
  - [Segment Tree](#segment-tree)
    - [Point update](#point-update)
    - [Lazy propagation](#lazy-propagation)
    - [Persistent](#persistent)
    - [Segment Tree 2D](#segment-tree-2d)
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
  - [MST](#mst)
  - [Lowest Common Ancestor](#lowest-common-ancestor)
  - [Bipartite Matching](#bipartite-matching)
  - [Fluxo](#fluxo)
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
  - [Área de polígono](#área-de-polígono)
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

#### Persistent

#### Segment Tree 2D

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
    
    vi p, rnk, cnt;
    int n;
    
    UnionFind(int n): n(n), p(vi(n)), rnk(vi(n, 0)), cnt(vi(n, 1)) {
        for(int i = 0; i < n; ++i)
            p[i] = i;
    }
    
    int find(int v) {
        return p[v] == v ? v : v = find(p[v]);
    }
    
    void uni(int u, int v) {
        int pu = find(u), pv = find(v);
        rnk[pu] < rnk[pv] ? p[pu] = pv : p[pv] = pu;
        rnk[pu] < rnk[pv] ? cnt[pv] += cnt[pu] : cnt[pu] += cnt[pv];
        if(rnk[pu] == rnk[pv]) ++rnk[pu];
    }
    
    int count(int v) {
        return cnt[find(v)];
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

### Floyd-Warshall

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

### MST

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

### Fluxo

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
    return a * b / mdc(a, b);
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

https://practice.geeksforgeeks.org/problems/check-if-two-line-segments-intersect/0

```c
struct Point {
    
    double x, y;

    Point(double x = 0, double y = 0): x(x), y(y) {}

    int orientation(Point p1, Point p2) {

        double d = (p1.y - this->y) * (p2.x - p1.x) -
                   (p1.x - this->x) * (p2.y - p1.y);
        if (fabs(d) < EPS) return 0; // colineares
        if (d > 0) return 1;         // horario
        return 2;                    // anti-horario
    }
};

struct Line {

    Point s, e;
    double dist;

    Line(Point s, Point e): s(s), e(e) {
        dist = hypot(fabs(e.x - s.x), fabs(e.y - s.y));
    }

    Line(double sx = 0, double sy = 0, double ex = 0, double ey = 0) {
        *this = Line(Point(sx, sy), Point(ex, ey));
    }

    //Retorna true se p e colinear com os
    //extremos do segmento (s e)
    bool collinear(Point p) {
        return !s.orientation(e, p);
    }

    //Retorna true se p esta no segmento
    //Deve ser usado apos collinear
    bool contains(Point p) {
        
        return fmin(s.x, e.x) <= p.x && fmax(s.x, e.x) >= p.x &&
               fmin(s.y, e.y) <= p.y && fmax(s.y, e.y) >= p.y;
    }

    bool intersect(Line other) {

        int o1 = this->s.orientation(this->e, other.s);
        int o2 = this->s.orientation(this->e, other.e);
        int o3 = other.s.orientation(other.e, this->s);
        int o4 = other.s.orientation(other.e, this->e);

        return (o1 != o2 && o3 != o4)  ||
               (!o1 && this->contains(other.s)) ||
               (!o2 && this->contains(other.e)) ||
               (!o3 && other.contains(this->s)) ||
               (!o4 && other.contains(this->e));
    }

    //Produto escalar
    double dot(Line other) {
        
        return (this->e.x - this->s.x) * (other.e.x - other.s.x) +
               (this->e.y - this->s.y) * (other.e.y - other.s.y);
    }

    //Produto vetorial
    double cross(Line other) {

        return (this->e.x - this->s.x) * (other.e.y - other.s.y) -
               (this->e.y - this->s.y) * (other.e.x - other.s.x);
    }
};
```

### Área de polígono

https://www.math10.com/en/geometry/geogebra/fullscreen.html
https://uva.onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&page=show_problem&problem=45

```c
struct PointSet {

    vector<Point> p;
    int n;

    PointSet(): n(0) {}

    PointSet(int n): n(n) {
        p = vector<Point>(n);
    }

    PointSet(vector<Point> points): p(points) {
        n = p.size();
    }

    double area() {

        double a = 0.0;
        for (int i = 1; i < p.size() - 1; ++i)
            a += Line(p[0], p[i]).cross(Line(p[0], p[i + 1]));
        return fabs(a * 0.5);
    }

    //Convex Hull
    vector<Point> grahamScan();

    //Par de pontos mais proximos (forca bruta)
    ppp closestByBruteForce(Point[], int);

    //Par de pontos mais proximos que estao a
    //uma distancia minima do ponto central
    ppp closestStrip(Point[], int, double);

    //Par de pontos mais proximos (Recursivo)
    ppp closestUtil(Point[], int);

    //Retorna o par de pontos mais proximos
    ppp closestPair();

    //Retorna a distancia entre os pontos
    //mais proximos
    double minimumDist();

    //Retorna verdadeiro se o ponto esta dentro
    //do poligono representado pelo Pointset
    bool inPolygon(Point);
};
```


### Convex Hull

https://practice.geeksforgeeks.org/problems/convex-hull/0
https://uva.onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&page=show_problem&problem=45

```c
Point p0; //Vertice inicial do convex hull

//Ordena pontos no sentido anti-horario
bool cmp(Point p1, Point p2) {

    double ori = Line(p0, p1).cross(Line(p0, p2));
    return ori == 0 ? Line(p0, p1).dist < Line(p0, p2).dist :
           atan2(p1.y - p0.y, p1.x - p0.x) < atan2(p2.y - p0.y, p2.x - p0.x);
}

//Retorna vetor de pontos do convex hull
vector<Point> PointSet::grahamScan() {

    vector<Point> newP;
    int iMin = 0;

    for (int i = 1; i < n; ++i)
        if (p[i].y < p[iMin].y || (p[i].y == p[iMin].y && p[i].x < p[iMin].x))
            iMin = i;
    swap(p[iMin], p[0]);
    p0 = p[0];
    newP.push_back(p0);
    sort(p.begin() + 1, p.end(), cmp);

    for (int i = 1; i < n; ++i) {
        while (i < n - 1 && fabs(Line(p0, p[i]).cross(Line(p0, p[i + 1]))) < EPS)
            i++;
        newP.push_back(p[i]);
    }

    vector<Point> poly(newP.size());
    if (newP.size() > 2) {
        for (int i = 0; i < 3; ++i)
            poly[i] = newP[i];
        int m = 3;
        for (int i = 3; i < newP.size(); ++i) {
            while (Line(poly[m - 2], poly[m - 1]).cross(Line(poly[m - 2], newP[i])) < EPS)
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

//Ordena pontos pelo x
bool cmpX(Point p1, Point p2) {
    return p1.x < p2.x;
}

//Ordena pontos pelo y
bool cmpY(Point p1, Point p2) {
    return p1.y < p2.y;
}

ppp PointSet::closestStrip(Point strip[], int m, double minD) {

    ppp res = {{-INF, -INF}, {INF, INF}};
    for (int i = 0; i < m - 1; ++i)
        for (int j = i + 1; j < m && strip[j].y - strip[i].y < minD; ++j)
            if (Line(strip[i], strip[j]).dist < minD) {
                minD = Line(strip[i], strip[j]).dist;
                res = {strip[i], strip[j]};
            }
    return res;
}

ppp PointSet::closestByBruteForce(Point vp[], int sz) {

    ppp res = {{-INF, -INF}, {INF, INF}};
    for (int i = 0; i < sz - 1; ++i)
        for (int j = i + 1; j < sz; ++j)
            if (Line(vp[i], vp[j]).dist < Line(res.first, res.second).dist)
                res = {vp[i], vp[j]};
    return res;
}

ppp PointSet::closestUtil(Point px[], int sz) {

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
    double plDist = Line(pl.first, pl.second).dist;
    double prDist = Line(pr.first, pr.second).dist;
    double d = fmin(plDist, prDist);
    res = plDist < prDist ? pl : pr;

    sort(px, px + sz, cmpY);

    Point strip[sz];
    int m = 0;
    for (int i = 0; i < sz; ++i)
        if (fabs(px[i].x - midP.x) < d)
            strip[m++] = px[i];

    ppp spDist = closestStrip(strip, m, d);
    if (Line(spDist.first, spDist.second).dist < d)
        return spDist;
    return res;
}

ppp PointSet::closestPair() {

    Point px[n];
    for (int i = 0; i < n; ++i) px[i] = p[i];
    sort(px, px + n, cmpX);
    return closestUtil(px, n);
}

double PointSet::minimumDist() {

    if (n < 2) return INF;
    ppp close = closestPair();
    return Line(close.first, close.second).dist;
}
```

### Ponto Dentro do Polígono

https://uva.onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&page=show_problem&problem=45

```c
//Retorna o angulo aob em radianos
double angle(Point a, Point o, Point b) {

    Line oa(o, a), ob(o, b);
    return acos(oa.dot(ob) / (oa.dist * ob.dist));
}

//Nao considera vertices como pontos internos
bool PointSet::inPolygon(Point pt) {

    double sum = 0;
    int j = n - 1;
    for (int i = 0; i < n; ++i) {
        double ang = angle(p[j], pt, p[i]);
        if (Line(pt, p[j]).cross(Line(pt, p[i])) > 0)
            sum += ang;
        else sum -= ang;
        j = i;
    }
    return fabs(fabs(sum) - 2 * M_PI) < EPS;
}
```
