def isStableWithSet(adj, u, set):
    """True si u puede ser agregado al set"""
    for v in set:
        if u in adj[v]:
            return False
    return True

def weightOfSet(S):
    """"Devuelve el peso total de un set"""
    return sum(S[v] for v in S)

