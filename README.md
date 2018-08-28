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

# Geometria Computacional
