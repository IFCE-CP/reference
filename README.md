# Sumário

- [**Estruturas de Dados**](#estruturas-de-dados)
  - [Segment Tree](#segment-tree)
    - [Point update](#point-update)
    - [Lazy propagation](#lazy-propagation)
    - [Persistent](#persistent)
    - [Segment Tree 2D](#segment-tree-2d)
  - [Fenwick Tree / BIT](#bit)
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
  - [LCA](#lca)
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

# Estruturas de dados

### Segment Tree

#### Point update

https://leetcode.com/problems/range-sum-query-mutable/description/

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

#### Persistent

#### Segment Tree 2D

### BIT

### Sparse Table

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

### Ordered Set

# Grafos

### Dijkstra

### Bellman-Ford

### Floyd-Warshall

### Pontes e Pontos de Articulação

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

### LCA

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

### Euclides Extendido

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
