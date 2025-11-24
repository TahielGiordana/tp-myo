import auxFuncs
import parserDimacs

#S=Conj.Estable, F=Resto de vertices, X=Vertices excluidos
def mwssRecursion(S,F,X, adj, maxIt):

    print(F)
    relevant_nodes = dict({key:value for key,value in F.items() if value > 0.0})
    print(f"RELEVANTE NODES: {relevant_nodes}")
    bestS = []
    bestW = 0.0
    n_it = 0
    pi_S = 0
    pi_F = auxFuncs.weightOfSet(F)
    #print(S)

    def mwssBucle(S,relevant_nodes,X,adj):

        nonlocal bestS, bestW, n_it, pi_S

        if n_it > maxIt:
            return False
        n_it += 1

        pi_S = auxFuncs.weightOfSet(S)
        pi_F = auxFuncs.weightOfSet(relevant_nodes)
        print(f"S: {S}")
        print(f"F: {relevant_nodes}")
        print(f"PI_S: {pi_S}")
        print(f"PI_F: {pi_F}")

        if pi_S > 1.0:
            bestS = list(S)
            bestW = pi_S
            return True

        # Si f ya está vacío retornamos
        if not relevant_nodes:
            print("Ya no hay restantes")
            return False
    
        # Elegimos el vértice de mayor peso
        v = max(relevant_nodes.items(), key=lambda item:item[1])
        print(f"V de mayor Peso: {v[0]}")
        print(f"Peso de V: {v[1]}")

        S2 = dict(S)
        S2[v[0]] = v[1]
        #print(f"S2: {S2}")
        #weights2 = pi_S + weights[v-1]

        # Nuevo F2 = F - {v} - N(v)
        F2 = dict(relevant_nodes)
        print(f"F2 Actual{F2}")
        F2.pop(v[0])
        print(f"F2 Despues{F2}")
        for u in adj[v[0]]:
            print(f"Vecino de {v[0]} : {u}")
            F2.pop(u,0)
        
        #print(f"F2: {F2}")

        X2 = X | adj[v[0]]
        #print(f"X2: {X2}")

        if mwssBucle(S2,F2,X2,adj):
            return True
    
        F3 = dict(relevant_nodes)
        F3.pop(v[0])
        X3 = X | {v[0]}
        if mwssBucle(S,F3,X3,adj):
            return True
    
        #print(f"S:{S} - F:{F} - X:{X}")

        return False
    
    mwssBucle(S,relevant_nodes,X,adj)
    print(f"Mejor S: {bestS}")
    print(f"Mejor w: {bestW}")
            
    return tuple(sorted(bestS)),bestW

if __name__ == "__main__":
    #n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/grafoTest")
    n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/DSJC125.1.col")
    #weights=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    weights=[-0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, 1.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0]
    #Vals: [1.0, 1.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0]
    S= set()
    F = set(range(1,n_nodos+1))
    X = set()

    bestS = []
    bestW = 0.0

    print(mwssRecursion(S,F,X,weights,adj))