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

def getConjuntoIndependienteMax(n, aristas):
    model = Model("ConjuntoIndependienteMax")

    # Variable binaria Xi indica si el nodo i pertenece al conjunto independiente.
    x = [model.addVar(vtype="B", name=f"x_{i}") for i in range(n)]

    # Si dos nodos son vecinos no pueden pertenecer a un conjunto independiente.
    for u,v in aristas:
        model.addCons(x[u] + x[v] <= 1)

    # Se busca maximizar el tamaño del conjunto independiente.
    model.setObjective(sum(x), "maximize")

    model.optimize()

    sol = model.getBestSol()
    if sol is None:
        sol = model.getSol()

    conjuntoFinal = []
    for i, var in enumerate(x):
        val = model.getSolVal(sol, var)
        if val is not None and val > 0.5:
            conjuntoFinal.append(i+1)
    
    return conjuntoFinal

def main():
    parser = argparse.ArgumentParser(description="Conjunto Independiente Máximo")
    parser.add_argument("input", help="Grafo en formato DIMACS")
    args = parser.parse_args()

    n, aristas = parserDimacs(args.input)
    print(f"Vertices: {n}, Aristas: {len(aristas)}")

    result = getConjuntoIndependienteMax(n, aristas)

    print(f"Tamaño del conjunto independiente: {len(result)}")
    print("Nodos: ")
    if result:
        print(" ".join(map(str,result)))
    else:
        print("Ninguno")

if __name__ == "__main__":
    main()