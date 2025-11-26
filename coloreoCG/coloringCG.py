import parserDimacs
import heuristics
import pyscipopt
import mwssRecursion
from pyscipopt import Model, SCIP_PARAMSETTING

if __name__ == "__main__":

    #n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/grafoTest")
    #n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/DSJC125.1.col") # Mejor K Heuristicas: 8 - Mejor k MWSS = 7 - Mejor conocido = ?
    #n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/DSJC250.5.col") # Mejor k Heuristicas = 41 - Mejor k MWSS = 33 - Mejor conocido = 26
    #n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/DSJC250.9.col") # Mejor k Heuristicas = 95 - Mejor k MWSS = 86 - Mejor conocido = 71
    #n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/flat300_20_0.col") # Mejor k Heuristicas = 43 - Mejor k MWSS = 20 - Mejor conocido = 20
    n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/le450_5c.col") # Mejor k heuristicas = 11 - Mejor k MWSS = 10 - Mejor conocido = 5
    n_nodos, n_aristas, adj = parserDimacs.parserDimacs("coloreoCG/queen12_12.col") # Mejor k heuristicas = 19 - Mejor k MWSS = 15 - Mejor conocido = 12

    g_densidad = (2*n_aristas)/(n_nodos*(n_nodos-1))

    print(f"Cantidad de Nodos={n_nodos}")
    print(f"Cantidad de Aristas={n_aristas}")
    print(f"Densidad = {g_densidad}")
    #for i in adj:
        #print(f"Vecinos de {i}: {adj.get(i)}")

    model = pyscipopt.Model("ColoringCG")

    model.setPresolve(0)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    model.disablePropagation()
    model.setParam("presolving/maxrounds",0)
    model.setParam("propagating/maxrounds",0)
    model.setParam("separating/maxrounds",0)
    model.setParam("lp/presolving",0)

    #model.setParam("limits/time", 20) # Limitamos el tiempo de ejecución

    model.hideOutput()

    model.setMinimize()

    constraints = {}
    variables = []
    nodes_weights = {}
    color_asign = {}

    # Problema maestro restringido
    for v in range(1,n_nodos+1):
        var = model.addVar(name=f"{v}_", vtype="C", lb=0.0, ub=1.0, obj=1.0)
        variables.append(var)

        cons = model.addCons(var >= 1, name=f"cover_{v}", separate=False,modifiable=True,removable=False)
        constraints[v] = cons

        nodes_weights[v] = 0
        color_asign[v] = v

    model.optimize()


    max_it = 100
    i=0
    while(i<=max_it):
        print(f"Iteración {i}")
        model.optimize()
        if model.getStatus() != 'optimal':
            print(f"ADVERTENCIA: EL LP NO ÓPTIMO. ESTADO:{model.getStatus()}")
        print(f"{model.getStatus()}")
        for v in range(n_nodos):
            pi = model.getDualsolLinear(model.getConss()[v])
            if pi > 1e15: nodes_weights[v+1] = 0.0
            else : nodes_weights[v+1] = pi

        bestW = 0.0
        bestS = []
        bestStrat = ""
        
        #Ejecutamos las 3 heurísticas
        '''
        #S,w = heuristics.greedy1(n_nodos,adj,vals)
        S,w = heuristics.greedy1(nodes_weights,adj)

        bestS = S
        bestW = w
        bestStrat = "Greedy1"

        #S2,w2 = heuristics.greedy2(n_nodos,adj,vals)
        S2,w2 = heuristics.greedy2(nodes_weights,adj)
        if w2 > bestW : 
            bestW = w2
            bestS = S2
            bestStrat = "Greedy2"
        #S3,w3 = heuristics.greedy3(n_nodos,adj,vals)
        S3,w3 = heuristics.greedy3(nodes_weights,adj)
        if w3 > bestW : 
            bestW = w3
            bestS = S3
            bestStrat = "Greedy3"

        #Si el peso es mayor a 1 se encontró una columna

        print(f"BestS: {bestS,bestW}")
        bestS, bestW = heuristics.improveStableSet(bestS, nodes_weights, adj)
        print(f"Nuevo BestS: {bestS, bestW}")
        '''
        if bestW > 1:
            print(f"Se encontró una columna S:{bestS} - w:{bestW} con {bestStrat}")
            varName = ''.join(str(val)+"_" for val in bestS)
            model.freeTransform()
            new_col_var = model.addVar( name=varName,vtype="C", obj=0.0, lb=0.0, ub=1.0)
            for v in bestS:
                model.addConsCoeff(model.getConss()[v-1], new_col_var, 1.0)
        else:      
                 
            print("Ejecutando MWSS Exacto")
            S = {}
            F = dict(nodes_weights) # {1,1,....,1} en la primer iteración
            X = set()
            mwssSol,mwssW = mwssRecursion.mwssRecursion(S=S,F=F,X=X,adj=adj,maxIt=200000)
            print(f"MWSS S{mwssSol}")
            print(f"MWSW S{mwssW}")

            model.freeTransform()
            if mwssW > 1.0:
                print("Agrego Columna por MWSS")
                varName = ''.join(str(val)+"_" for val in mwssSol)         
                new_mwss = model.addVar( name=varName,vtype="C", obj=0.0, lb=0.0, ub=1.0)
                for v in mwssSol:
                    model.addConsCoeff(model.getConss()[v-1], new_mwss, 1.0)
            else:
                break
            
            #break
        i=i+1

    model.optimize()

    print(f"ITERACIONES: {i}")

    solucion = model.getBestSol()
    solution_dict = {}
    variables = model.getVars()
    for var in variables:
        solution_dict[var.name] = model.getSolVal(solucion,var)
    
    stableSets = dict({key:value for key,value in solution_dict.items() if value > 0.0})
    
    color = 0
    for sol in stableSets:
        color += 1
        nodes = sol.split('_')
        nodes.pop()
        for node in nodes:
            intNode = int(node)
            color_asign[intNode] = color

    #Salida
    # Exito: s optimal <k> o s feasible <k>
    # No existe: s unsatisfiable
    # v <node_id> <color_id>

    output_comment = 'c Coloreo generado por coloreo HCS'
    output_status = 's ' + model.getStatus() + " " + str(color)
    output_coloring = ''
    for node in color_asign:
        output_coloring += "v "+ str(node) + " " + str(color_asign[node]) + "\n"
    #print(output_comment)
    print(output_status)
    #print(output_coloring)

    print(stableSets)
    #print(solucion)

    with open("oputputHCS", 'w') as f:
        f.write(output_comment+"\n")
        f.write(output_status+"\n")
        f.write(output_coloring)