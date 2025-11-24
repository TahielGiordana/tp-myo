import auxFuncs as aux
import heuristics as gr

def parserDimacs(path):
    n_nodos = 0
    n_aristas = 0
    adj = None
    try:
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
            
                parts = line.split()

                # Línea de comentario
                if parts[0] == 'c':
                    continue

                # Línea del problema: p edge N M
                elif parts[0] == 'p':
                    if parts[1] == 'edge':
                        n_nodos = int(parts[2])
                        n_aristas = int(parts[3])
                        adj = {i: set() for i in range(1,n_nodos+1)}

                # Línea de arista: e u v
                elif parts[0] == 'e':
                    try:
                        u = int(parts[1])
                        v = int(parts[2])

                        # Agregamos la arista (u,v) y (v,u)
                        adj[u].add(v)
                        adj[v].add(u)
                    except ValueError:
                        # Ignoramos las líneas que no sigan el formato
                        continue

    except FileNotFoundError:
        print(f"No se encuentra el archvo en la ruta {path}")
        return None

    return n_nodos, n_aristas, adj

if __name__ == "__main__":
    n_nodos, n_aristas, adj = parserDimacs("coloreoCG/grafoTest")
    print(f"Cantidad de Nodos={n_nodos}")
    print(f"Cantidad de Aristas={n_aristas}")
    for i in adj:
        print(f"Vecinos de {i}: {adj.get(i)}")

    gr.greedy1(n_nodos, adj, [2,3,7,6,2,3,4,5])
    gr.greedy2(n_nodos, adj, [2,3,7,6,2,3,4,5])
    gr.greedy3(n_nodos, adj, [2,3,7,6,2,3,4,5])
    
