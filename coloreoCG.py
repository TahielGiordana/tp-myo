from pyscipopt import Model
from collections import defaultdict
import time
import random
import struct

def decode_dimacs_binary_graph(file_path):
    """
    Decodifica un archivo binario de grafos DIMACS (.col.b)
    y devuelve la lista de aristas en el formato ASCII estándar.
    """
    
    # Parsing de la cabecera (ASCII)
    try:
        with open(file_path, 'r', encoding='latin-1') as f:
            header_lines = []
            n_nodes = 0
            n_edges = 0
            
            # Leer la cabecera hasta encontrar la línea 'p'
            for line in f:
                header_lines.append(line.strip())
                parts = line.split()
                if parts and parts[0] == 'p' and parts[1] == 'edge':
                    n_nodes = int(parts[2])
                    n_edges = int(parts[3])
                    break
            
            if n_nodes == 0:
                raise ValueError("No se encontró la línea 'p edge N M' en la cabecera.")
            
            # Calcular la posición donde termina la cabecera y empieza el binario
            #binary_start_position = f.tell()

    except FileNotFoundError:
        print(f"Error: Archivo no encontrado en la ruta {file_path}")
        return None
    except Exception as e:
        print(f"Error al leer la cabecera: {e}")
        return None


    # Procesar la Bit Matrix
    
    # Abrir el archivo en modo binario ('rb') para leer el cuerpo
    with open(file_path, 'rb') as f:
        #f.seek(binary_start_position)
        binary_data = f.read()

    decoded_edges = []
    bit_index = 0
    
    # Iterar sobre las coordenadas (i, j) en el triángulo superior
    # Los nodos DIMACS son 1-indexed (1 a N_nodes)
    for i in range(1, n_nodes + 1):
        for j in range(i + 1, n_nodes + 1):
            
            # Calcular el byte y el bit dentro de ese byte
            byte_offset = bit_index // 8
            bit_in_byte = bit_index % 8
            
            # Si nos quedamos sin datos antes de terminar el triángulo superior, algo salió mal o el archivo está truncado.
            if byte_offset >= len(binary_data):
                # Esto es padding o fin inesperado
                break 

            # Obtener el byte actual
            current_byte = binary_data[byte_offset]
            
            # Extraer el bit (asumiendo MSB primero: el primer bit es el más significativo)
            # El bit se encuentra en la posición 7 - bit_in_byte
            bit_value = (current_byte >> (7 - bit_in_byte)) & 1
            
            # Si el bit es 1, existe una arista (i, j)
            if bit_value == 1:
                decoded_edges.append((i, j))
            
            bit_index += 1

    # Generar DIMACS ASCII
    
    # El número M_decoded de aristas debe coincidir con el M del header (21375)
    M_decoded = len(decoded_edges)
    
    output_lines = header_lines[:] # Copiamos la cabecera original
    
    # Actualizar la línea 'p' si es necesario (el original ya tiene el valor correcto)
    output_lines[len(header_lines)-1] = f"p edge {n_nodes} {n_edges}"
    
    # Agregar las aristas decodificadas
    for u, v in decoded_edges:
        output_lines.append(f"e {u} {v}")

    # Devolver la lista de líneas ASCII o el texto completo
    return "\n".join(output_lines)


def parse_dimacs_file(file_path):
    """
    Lee un archivo de grafos en formato DIMACS (.col) y devuelve un 
    diccionario de adyacencia.
    """
    adj = defaultdict(set)
    n_nodes = 0
    n_edges = 0
    
    with open(file_path, 'r', encoding='latin-1') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            
            # Línea de comentario
            if parts[0] == 'c':
                continue
            
            # Línea de problema: p edge N M
            elif parts[0] == 'p':
                if parts[1] == 'edge':
                    n_nodes = int(parts[2])
                    n_edges = int(parts[3])
                    # Inicializar el diccionario con todos los nodos (0 a N-1)
                    # El formato DIMACS generalmente usa 1-based indexing, ajustamos a 0-based.
                    adj = {i: set() for i in range(1, n_nodes + 1)}
            
            # Línea de arista: e U V
            elif parts[0] == 'e':
                try:
                    u = int(parts[1])
                    v = int(parts[2])
                    
                    # Añadir la arista (u, v) y (v, u)
                    adj[u].add(v)
                    adj[v].add(u)
                except ValueError:
                    # Ignorar líneas mal formadas
                    continue
            
    # Devolver el grafo en formato 1-based (como lo lee el parser)
    return adj

class GraphColoringCG:
    def __init__(self, adj_list):
        """
        Inicializa el problema de coloreo mediante generación de columnas.
        :param adj_list: Diccionario donde key=vertice, value=set(vecinos)
        """
        self.adj = adj_list
        self.nodes = list(adj_list.keys())
        self.n = len(self.nodes)
        
        # Inicializar PySCIPOpt Model (Master Problem)
        self.model = Model("FractionalColoring_LP")
        
        # Desactivar Presolving para evitar que modifique la estructura del problema
        self.model.setPresolve(0)
        self.model.setParam("presolving/maxrounds", 0)
        self.model.setParam("presolving/maxrestarts", 0)
        
        # Desactivar Propagación: Evita que SCIP fije variables "obvias" 
        # y elimine las restricciones asociadas, lo cual causa el error de duales NULL.
        self.model.setParam("propagating/maxrounds", 0)
        self.model.setParam("propagating/maxroundsroot", 0)

        self.model.setParam("separating/maxrounds", 0)
        self.model.setParam("separating/maxroundsroot", 0)

        self.model.setParam("lp/presolving", False)
        
        self.model.hideOutput()
        
        self.conss = {}
        self.vars = []
        
        self._init_master_problem()

    def _init_master_problem(self):
        """
        Inicializa variables y restricciones simultáneamente.
        """
        for i, v in enumerate(self.nodes):
            # Al no poner cota superior estricta de 1.0, evitamos que el propagador
            # detecte la solución trivial inmediata y elimine la restricción.
            var = self.model.addVar(name=f"S_init_{v}", vtype="C", lb=0.0, ub=None, obj=1.0)
            self.vars.append(var)
            
            # Crear restricción: var >= 1
            # modifiable=True es vital para poder agregar coeficientes luego
            cons = self.model.addCons(var >= 1, name=f"cover_{v}", separate=False, modifiable=True, removable=False)
            self.conss[v] = cons

    def add_column(self, stable_set):
        """Agrega una nueva variable (columna) al modelo."""
        col_idx = len(self.vars)
        # También usamos ub=None para las nuevas columnas por consistencia
        var = self.model.addVar(name=f"S_{col_idx}", vtype="C", lb=0.0, ub=None, obj=1.0)
        self.vars.append(var)
        
        for v in stable_set:
            self.model.addConsCoeff(self.conss[v], var, 1.0)

    def get_dual_values(self):
        """Obtiene los valores duales (pi_v)."""
        duals = {}
        for v in self.nodes:
            # Ahora que la propagación está desactivada y la restricción es 'modifiable',
            # SCIP mantendrá la fila en la matriz LP y podremos leer su dual.
            pi = self.model.getDualsolLinear(self.conss[v])
            if pi > 1e15:
                duals[v] = 0.0
            else:
                duals[v] = pi
        return duals

    def solve(self, max_iter=100):
        """Ejecuta el bucle de generación de columnas."""
        print(f"{'Iter':<5} | {'LP Obj':<10} | {'Heuristic':<15} | {'Weight':<10} | {'Size':<5}")
        print("-" * 60)

        for it in range(max_iter):
            # 1. Optimizar
            self.model.optimize()
            lp_obj = self.model.getObjVal()
            
            # 2. Obtener duales (Antes de liberar transformación)
            duals = self.get_dual_values()
            
            # 3. Resolver Pricing (Heurísticas MWSS)
            stable_set, weight, method = self.run_mwss_heuristics(duals)
            
            # 4. Criterio de parada: Peso <= 1 indica que no hay columnas con costo reducido negativo
            # (Nota: costo reducido = 1 - weight. Si weight <= 1, costo reducido >= 0 -> óptimo)
            if weight <= 1.0 + 1e-6:
                print("-" * 60)
                print(f"Terminado: No se encontraron conjuntos con peso > 1 (Heurísticas agotadas).")
                break
            
            print(f"{it:<5} | {lp_obj:<10.4f} | {method:<15} | {weight:<10.4f} | {len(stable_set):<5}")
            
            # 5. Liberar transformación para modificar el modelo
            self.model.freeTransform()
            
            # 6. Agregar columna
            self.add_column(stable_set)

        self.model.optimize()
        return self.model.getObjVal()

    # Heurísticas
    def run_mwss_heuristics(self, duals):

        S = self.greedy_strategy_2(duals)
        w = sum(duals[v] for v in S)
        if w > 1.00001: return S, w, "DynSurplus"

        S = self.greedy_strategy_3(duals)
        w = sum(duals[v] for v in S)
        if w > 1.00001: return S, w, "StatSurplus"
        
        S = self.greedy_strategy_1(duals)
        w = sum(duals[v] for v in S)
        if w > 1.00001: return S, w, "MaxWeight"
            
        return [], 0.0, "None"

    def _build_greedy_stable_set(self, sorted_candidates, duals):
        S = set()
        forbidden = set()
        for v in sorted_candidates:
            if v not in forbidden:
                S.add(v)
                forbidden.update(self.adj[v])
        return list(S)

    def greedy_strategy_1(self, duals):
        candidates = [v for v in self.nodes if duals[v] > 1e-6]
        candidates.sort(key=lambda v: duals[v], reverse=True)
        return self._build_greedy_stable_set(candidates, duals)

    def greedy_strategy_3(self, duals):
        candidates = [v for v in self.nodes if duals[v] > 1e-6]
        def static_score(v):
            return duals[v] - sum(duals[w] for w in self.adj[v])
        candidates.sort(key=static_score, reverse=True)
        return self._build_greedy_stable_set(candidates, duals)

    def greedy_strategy_2(self, duals):
        S = set()
        N_S = set()
        candidates = set(v for v in self.nodes if duals[v] > 1e-6)
        while candidates:
            best_v = None
            best_score = -float('inf')
            for v in candidates:
                pi_v = duals[v]
                penalty = 0.0
                for w in self.adj[v]:
                    if w not in N_S:
                        penalty += duals[w]
                score = pi_v - penalty
                if score > best_score:
                    best_score = score
                    best_v = v
            if best_v is None: break
            S.add(best_v)
            new_neighbors = self.adj[best_v]
            N_S.update(new_neighbors)
            candidates.remove(best_v)
            candidates = {v for v in candidates if v not in new_neighbors}
        return list(S)
    
# Función auxiliar para testear grafos aleatorios
def generate_random_graph(n, p):
    """Genera un grafo Erdos-Renyi G(n, p)."""
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < p:
                adj[i].add(j)
                adj[j].add(i)
    return adj

if __name__ == "__main__":

    #FILE_PATH = "DSJC125.1.col.b"  # Mejor k conseguido = 16
    #FILE_PATH = "DSJC250.5.col.b"  # Mejor k conseguido = 52
    #FILE_PATH = "DSJC250.9.col.b"  # Mejor k conseguido = 165
    FILE_PATH = "flat300_20_0.col"  # Mejor k conseguido = 40 

    print(f"Iniciando decodificación binaria de {FILE_PATH}...")
    
    ascii_graph_content = decode_dimacs_binary_graph(FILE_PATH)
    output_file_path = FILE_PATH.replace(".col.b", ".col")

    if ascii_graph_content:
        # Imprimir las primeras 10 líneas de la salida
        print("\n--- Cabecera y Primeras Aristas Decodificadas (DIMACS ASCII) ---")
        print("\n".join(ascii_graph_content.split('\n')[:15]))
        
        # Opcional: Guardar la salida decodificada en un nuevo archivo .col
        
        with open(output_file_path, 'w') as f:
             f.write(ascii_graph_content)
        
        print(f"\nDecodificación completa. {len(ascii_graph_content.splitlines()) - 10} aristas encontradas.")
        print(f"El grafo decodificado se ha guardado en: {output_file_path}")

    print(f"1. Parseando archivo DIMACS: {FILE_PATH}...")
    try:
        adj_list = parse_dimacs_file(output_file_path)
        N = len(adj_list)
        E = sum(len(neighbors) for neighbors in adj_list.values()) // 2
        print(f"   Grafo cargado: {N} nodos, {E} aristas.")
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo '{FILE_PATH}'. Asegúrate de que esté en el directorio correcto.")
        exit()

    #N_NODES = 150
    #DENSITY = 0.5
    
    #print(f"Generando grafo aleatorio ({N_NODES} nodos, densidad {DENSITY})...")
    #adj_list = generate_random_graph(N_NODES, DENSITY)

    print("Iniciando Generación de Columnas...")
    cg_solver = GraphColoringCG(adj_list)
    final_obj = cg_solver.solve()
    print(f"\nResultado: {final_obj:.4f}")