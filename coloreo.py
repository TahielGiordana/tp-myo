import argparse
from pyscipopt import Model

def parserDimacs(path):
    n = None
    aristas = []
    max_index = -1
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line[0] == 'c':
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
                    aristas.append((u,v))
                    max_index = max(max_index, u, v)
                except ValueError:
                    continue
    if n is None:
        if max_index >= 0:
            n = max_index + 1
        else:
            n = 0
    return n, aristas

def getColoreoConjEstables(n, aristas, max_colors=None):
    """
    Modelo de coloreo de grafos basado en conjuntos estables.
    Cada color es un conjunto independiente.
    """
    model = Model("ColoreoConjuntosEstables")

    # Si no se especifica, usar un número máximo de colores igual a n
    if max_colors is None:
        max_colors = n

    # Matriz de adyacencia
    adj = [[False]*n for _ in range(n)]
    for u, v in aristas:
        adj[u][v] = True
        adj[v][u] = True

    # Variables binarias:
    # x[v][c] = 1 si el vértice v tiene el color c
    x = {}
    for v in range(n):
        for c in range(max_colors):
            x[v, c] = model.addVar(vtype="B", name=f"x_{v}_{c}")

    # y[c] = 1 si el color c se usa
    y = {c: model.addVar(vtype="B", name=f"y_{c}") for c in range(max_colors)}

    # Cada vértice debe tener exactamente un color
    for v in range(n):
        model.addCons(sum(x[v, c] for c in range(max_colors)) == 1)

    # Si dos vértices son vecinos, no pueden tener el mismo color
    for (u, v) in aristas:
        for c in range(max_colors):
            model.addCons(x[u, c] + x[v, c] <= 1)

    # Un vértice solo puede tener color c si ese color se usa
    for v in range(n):
        for c in range(max_colors):
            model.addCons(x[v, c] <= y[c])

    # Minimizar cantidad de colores usados
    model.setObjective(sum(y[c] for c in range(max_colors)), "minimize")

    model.optimize()

    sol = model.getBestSol()
    if not sol:
        return None

    # Construir solución
    color_asignado = {}
    for v in range(n):
        for c in range(max_colors):
            if sol[x[v, c]] > 0.5:
                color_asignado[v] = c
                break

    usados = [c for c in range(max_colors) if sol[y[c]] > 0.5]
    return color_asignado, usados, len(usados)


def getColoreoRepresentantes(n, aristas):
    model = Model("Coloreo_Representantes")

    # Matriz de adyacencia
    adj = [[False]*n for _ in range(n)]
    for u, v in aristas:
        adj[u][v] = True
        adj[v][u] = True

    # Ñ[v] = no vecinos de v ∪ {v}
    Ntil = {v: [u for u in range(n) if not adj[v][u] or u == v] for v in range(n)}

    # Variables: x[u,v] = 1 si u representa a v
    x = {}
    for u in range(n):
        for v in Ntil[u]:
            x[u, v] = model.addVar(vtype="B", name=f"x_{u}_{v}")

    # Cada vértice tiene exactamente un representante
    for v in range(n):
        model.addCons(sum(x[u, v] for u in Ntil[v]) == 1)

    # Un vértice u solo puede representar si es representante
    for u in range(n):
        for v in Ntil[u]:
            model.addCons(x[u, v] <= x[u, u])

    # Los vértices representados por u deben formar un conjunto estable
    for u in range(n):
        for (i, j) in aristas:
            if i in Ntil[u] and j in Ntil[u]:
                model.addCons(x[u, i] + x[u, j] <= 1)

    # Minimizar cantidad de representantes (colores)
    model.setObjective(sum(x[u, u] for u in range(n)), "minimize")

    # Resolver
    model.optimize()

    sol = model.getBestSol()
    if not sol:
        return None

    colores = {}
    clases_color = {}

    # Asignar representante a cada vértice
    for v in range(n):
        for u in Ntil[v]:
            if sol[x[u, v]] > 0.5:
                colores[v] = u
                clases_color.setdefault(u, []).append(v)
                break

    k = sum(sol[x[u, u]] > 0.5 for u in range(n))
    return colores, clases_color, k

def main():
    parser = argparse.ArgumentParser(description="Coloreo")
    parser.add_argument("input", help="Grafo en formato DIMACS")
    args = parser.parse_args()

    n, aristas = parserDimacs(args.input)
    print(f"Vertices: {n}, Aristas: {len(aristas)}")

    colores, clases_color, k = getColoreoConjEstables(n, aristas)

    print("Color de cada vértice:")
    for v in range(n):
        print(f"  Vértice {v} → Color {colores[v]}")

    print(f"\nColores usados: {clases_color}")
    print(f"Número mínimo de colores: {k}")

    colores, clases_color, k = getColoreoRepresentantes(n, aristas)

    print("Color (representante) asignado a cada vértice:")
    for v in range(n):
        print(f"  Vértice {v} -> Representante {colores[v]}")

    print("\nClases de color:")
    for rep, grupo in clases_color.items():
        print(f"  Color {rep}: {grupo}")

    print(f"\nNúmero mínimo de colores usados: {k}")

if __name__ == "__main__":
    main()