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

#### Lazy propagation

#### Persistent

#### Segment Tree 2D

### BIT

### Sparse Table

### SQRT Decomposition

### Trie

### Union Find

### Ordered Set

# Grafos

### Dijkstra

### Bellman-Ford

### Floyd-Warshall

### Pontes e Pontos de Articulação

### Componentes Fortemente Conexos

### Ordenação Topológica

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
