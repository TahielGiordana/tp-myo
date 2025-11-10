"""
Implementación de las tres heurísticas mencionadas en el paper.
"""
import argparse
import random
import itertools
import math
from collections import defaultdict

def parse_dimacs(path):
    n = None
    edges = []
    max_idx = -1
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('c'):
                continue
            parts = line.split()
            if parts[0] == 'p' and len(parts) >= 4:
                try:
                    n = int(parts[2])
                except:
                    n = None
            elif parts[0] == 'e' and len(parts) >= 3:
                try:
                    u = int(parts[1]) - 1
                    v = int(parts[2]) - 1
                except ValueError:
                    continue
                edges.append((u, v))
                max_idx = max(max_idx, u, v)
    if n is None:
        n = max_idx + 1 if max_idx >= 0 else 0
    # Lista de vecinos (Conveniente ya que se consultan continuamente los vecinos de un vector)
    adj = [set() for _ in range(n)]
    for u, v in edges:
        if 0 <= u < n and 0 <= v < n and u != v:
            adj[u].add(v)
            adj[v].add(u)
    return n, edges, adj

# Peso total de un conjunto 
def pesoSet(S, pesos):
    return sum(pesos[v] for v in S)

# Check conjunto es independiente
def esIndependiente(S, adj):
    Sset = set(S)
    for v in Sset:
        if any((u in Sset) for u in adj[v]):
            return False
    return True

# Heurística 1: Ordenar nodos por peso decreciente y agregar mientras se cumpla la independencia
def greedy_1(n, adj, weights):
    order = sorted(range(n), key=lambda v: weights[v], reverse=True)
    S = []
    elegidos = set()
    for v in order:
        # Si v no tiene vecino en los elegidos, lo agrego
        if not any(u in elegidos for u in adj[v]):
            elegidos.add(v)
            S.append(v)
    return tuple(sorted(S)), pesoSet(S, weights)

# -----------------------
# Heurística 2: En cada iteración calculamos un surplus para seleccionar el siguiente nodo
# -----------------------
def greedy_2(n, adj, weights):
    R = set(range(n))
    S = []
    elegidos = set()
    while R:
        scores = {}
        for v in R:
            costo = 0
            for u in adj[v]:
                if u in R:
                    costo += weights[u]
            scores[v] = weights[v] - costo
        # Elegimos el nodo mas conveniente
        v_star = max(R, key=lambda x: (scores[x], weights[x]))
        # Lo agregamos a S si es posible
        if not any(u in elegidos for u in adj[v_star]):
            elegidos.add(v_star)
            S.append(v_star)
        # Eliminamos el nodo de R
        R.remove(v_star)
    return tuple(sorted(S)), pesoSet(S, weights)

# Heurística 3: Similar a lo anterior pero con un surplus estático.
def greedy_3(n, adj, weights):
    score_static = [0.0] * n
    for v in range(n):
        s = 0.0
        for u in adj[v]:
            s += weights[u]
        score_static[v] = weights[v] - s
    orden = sorted(range(n), key=lambda v: (score_static[v], weights[v]), reverse=True)
    S = []
    elegidos = set()
    for v in orden:
        if not any(u in elegidos for u in adj[v]):
            elegidos.add(v)
            S.append(v)
    return tuple(sorted(S)), pesoSet(S, weights)

# Invoca las 3 heurísticas
def ejecutarHeuristicas(n, adj, weights):
    """
    Ejecuta greedy_1, greedy_2, greedy_3; sobre cada solución intenta local_search_1_2_swap;
    devuelve la mejor solución (set tuple ordenada) y su peso.
    """
    results = []
    s1, w1 = greedy_1(n, adj, weights)
    results.append(('Heurística-1', s1, w1))
    s2, w2 = greedy_2(n, adj, weights)
    results.append(('Heurística-2', s2, w2))
    s3, w3 = greedy_3(n, adj, weights)
    results.append(('Heurística-3', s3, w3))

    best_name, best_S, best_w = None, (), -1.0

    for name, S, w in results:
        if w > best_w:
            best_name = name
            best_S = S
            best_w = w

    # Devolvemos el mejor resultado
    return best_name, best_S, best_w, results

# Testing
def readRandomWeights(n, seed=None, scale=1000.0):
    random.seed(seed)
    return [random.random() * scale for _ in range(n)]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Archivo DIMACS (edge list)")
    parser.add_argument("--weights-file", "-w", default=None, help="Archivo con pesos (una línea por vértice)")
    parser.add_argument("--random-weights", action="store_true", help="Generar pesos aleatorios para test")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--scale", type=float, default=1000.0)
    args = parser.parse_args()

    n, edges, adj = parse_dimacs(args.input)
    if args.random_weights:
        weights = readRandomWeights(n, seed=args.seed, scale=args.scale)
    else:
        weights = [1.0] * n

    best_name, best_S, best_w, all_results = ejecutarHeuristicas(n, adj, weights)

    print(f"Grafo n={n}, m={len(edges)}")

    print("-" * 40)
    print(f"Mejor Heurística: {best_name}")
    print(f"Tamaño máximo de conjunto estable: {len(best_S)}, weight: {best_w:.6f}")
    print("Vertices:", " ".join(str(v+1) for v in best_S))

if __name__ == "__main__":
    main()
