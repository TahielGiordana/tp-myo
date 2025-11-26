"""Conjunto de Heurísticas para encontrar un conjunto estable"""
import parserDimacs

import auxFuncs as aux

def greedy1(n, adj, weights):
    """
    Ordena los nodos por peso decreciente y 
    agrega mientras se cumpla la independencia
    """
    order = sorted(range(n),key= lambda v:weights[v], reverse=True)
    #print(order)
    S = []
    for v in order:
        if not any (u in S for u in adj[v+1]):
            S.append(v+1)
    print(f"Greddy1: {sorted(S)}")
    return tuple(sorted(S)), aux.weightOfSet(S, weights)

def greedy1(nodes_weights, adj):
    """
    Ordena los nodos por peso decreciente y 
    agrega mientras se cumpla la independencia
    """
    order = sorted(nodes_weights.items(),key= lambda nodes_weights:nodes_weights[1], reverse=True)
    #print(order)
    order = dict(order)
    relevant_nodes = [v for v in order if order[v]>0.0]
    print(relevant_nodes)
    #print(order)
    S = []
    weight = 0
    for v in relevant_nodes:
        if not any (u in S for u in adj[v]):
            S.append(v)
            weight += order[v]
    print(f"Greddy1: {sorted(S)}")
    return tuple(sorted(S)), weight

def greedy2(n,adj,weights):
    """
    En cada iteración calculamos un surplus para seleccionar el siguiente
    nodo a agregar
    """
    R = set(range(1,n+1))
    S = []
    while R:
        scores = {}
        for v in R:
            costo = 0
            for u in adj[v]:
                if u in R:
                    costo += weights[u-1]
            scores[v] = weights[v-1] - costo
        
        best_v = max(R,key=lambda x:(scores[x], weights[x-1]))

        if not any(u in S for u in adj[best_v]):
            S.append(best_v)

        R.remove(best_v)
    print(f"Greddy2: {sorted(S)}")
    return tuple(sorted(S)), aux.weightOfSet(S,weights)

def greedy2(nodes_weights,adj):
    """
    En cada iteración calculamos un surplus para seleccionar el siguiente
    nodo a agregar
    """
    relevant_nodes = [v for v in nodes_weights if nodes_weights[v]>0.0]
    R = set(relevant_nodes)
    #print(R)
    S = []
    weight = 0
    while R:
        scores = {}
        for v in R:
            costo = 0
            for u in adj[v]:
                if u in R:
                    costo += nodes_weights[u]
            scores[v] = nodes_weights[v] - costo
        
        best_v = max(R,key=lambda x:(scores[x], nodes_weights[x]))

        if not any(u in S for u in adj[best_v]):
            S.append(best_v)
            weight += nodes_weights[best_v]

        R.remove(best_v)
    print(f"Greddy2: {sorted(S)}")
    return tuple(sorted(S)), weight

def greedy3(n,adj,weights):
    """
    Calculamos un surplus estático al comienzo de la ejecución
    """
    score_static = [0.0] * n
    for v in range(n):
        score = 0.0
        for u in adj[v+1]:
            score += weights[u-1]
        score_static[v] = weights[v] - score
    order = sorted(range(1,n+1),key=lambda v: (score_static[v-1], weights[v-1]),reverse=True)
    S = []
    for v in order:
        if not any (u in S for u in adj[v]):
            S.append(v)
    print(f"Greddy3: {sorted(S)}")
    return tuple(sorted(S)), aux.weightOfSet(S,weights)

def greedy3(nodes_weights,adj):
    """
    Calculamos un surplus estático al comienzo de la ejecución
    """
    relevant_nodes = [v for v in nodes_weights if nodes_weights[v]>0.0]
    score_static = {}
    for v in relevant_nodes:
        score = 0.0
        for u in adj[v]:
            score += nodes_weights[u]
        score_static[v] = nodes_weights[v] - score

    order = sorted(score_static.items(),key=lambda score_static:score_static[1],reverse=True)
    order = dict(order)
    #print(order)
    S = []
    weight = 0
    for v in order:
        if not any (u in S for u in adj[v]):
            S.append(v)
            weight += nodes_weights[v]
    print(f"Greddy3: {sorted(S)}")
    return tuple(sorted(S)), weight

def improveStableSet(S,nodes_weights,adj):
    S_actual = set(S)
    #print(f"S Actual: {S_actual}")
    improved = True

    while improved:
        improved = False

        peso_actual = sum(nodes_weights.get(v,0) for v in S_actual)
        #print(f"Peso Actual: {peso_actual}")

        # Sacamos un nodo de S y agregamos sus vecinos
        list_S = list(S_actual)
        for u in list_S:
            #print(f"U: {u}")
            peso_u = nodes_weights.get(u,0)
            #print(f"Peso de U: {peso_u}")

            vecinos_u = set(adj.get(u,[]))
            #print(f"Vecinos de U: {vecinos_u}")

            candidatos = []
            S_minus_u = S_actual - {u}

            for v in vecinos_u:
                compatible = True
                for s_node in S_minus_u:
                    if s_node in adj.get(v,[]):
                        compatible = False
                        break
                if compatible:
                    candidatos.append(v)

            candidatos.sort(key=lambda x: nodes_weights.get(x,0), reverse=True)

            added_nodes = []
            added_weight = 0

            temp_added = set()
            for cand in candidatos:
                is_indep = True
                for i in temp_added:
                    if i in adj.get(cand, []):
                        is_indep = False
                        break
                if is_indep:
                    temp_added.add(cand)
                    added_weight+=nodes_weights.get(cand,0)

            if added_weight > peso_u:
                S_actual.remove(u)
                S_actual.update(temp_added)
                improved = True
                break
            
            #print(f"Candidatos: {candidatos}")
    
    return tuple(sorted(S_actual)), sum(nodes_weights.get(v,0) for v in S_actual)

if __name__ == "__main__":
    n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/grafoTest")
    print(f"Cantidad de Nodos={n_nodos}")
    print(f"Cantidad de Aristas={n_aristas}")
    for i in adj:
        print(f"Vecinos de {i}: {adj.get(i)}")

    S = {}
    S[1]=1.0
    S[2]=1.0
    S[3]=1.0
    S[4]=1.0
    S[5]=1.0
    S[6]=1.0
    S[7]=1.0
    S[8]=1.0
    #print(greedy2(S, adj))
    #S2, w2 = greedy2(S,adj)
    S2 = list([1,3,5])
    print(f"Set Mejorado: {improveStableSet(S2,S,adj)}")

    #gr.greedy2(n_nodos, adj, [2,3,7,6,2,3,4,5])
    #gr.greedy3(n_nodos, adj, [2,3,7,6,2,3,4,5])