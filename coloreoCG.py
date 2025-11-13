import argparse
import time
from pyscipopt import Model

# ============================================================
#  DIMACS parser
# ============================================================
def parse_dimacs(path):
    n = None
    edges = []
    max_idx = -1

    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('c'):
                continue
            parts = s.split()

            if parts[0] == 'p':
                n = int(parts[2])

            elif parts[0] == 'e':
                u = int(parts[1]) - 1
                v = int(parts[2]) - 1
                edges.append((u, v))
                max_idx = max(max_idx, u, v)

    if n is None:
        n = max_idx + 1 if max_idx >= 0 else 0

    adj = [set() for _ in range(n)]
    for u, v in edges:
        if u != v:
            adj[u].add(v)
            adj[v].add(u)

    return n, edges, adj

# ============================================================
#  Initial Greedy Coloring
# ============================================================
def initial_coloring_columns(n, adj):
    order = sorted(range(n), key=lambda v: len(adj[v]), reverse=True)
    color_sets = []
    color_of = [-1]*n

    for v in order:
        assigned = False
        for cidx, cset in enumerate(color_sets):
            if all((u not in adj[v]) for u in cset):
                cset.append(v)
                color_of[v] = cidx
                assigned = True
                break
        if not assigned:
            color_sets.append([v])
            color_of[v] = len(color_sets)-1

    return [tuple(sorted(c)) for c in color_sets]

# ============================================================
#  Greedy heuristics (no local search)
# ============================================================
def weight_of_set(S, w):
    return sum(w[v] for v in S)

def greedy_1(n, adj, w):
    order = sorted(range(n), key=lambda v: w[v], reverse=True)
    S = []
    chosen = set()

    for v in order:
        if not any(u in chosen for u in adj[v]):
            chosen.add(v)
            S.append(v)

    return tuple(sorted(S)), weight_of_set(S, w)

def greedy_2(n, adj, w):
    R = set(range(n))
    chosen = set()
    S = []

    while R:
        best_v = None
        best_score = None

        for v in R:
            penalty = sum(w[u] for u in adj[v] if u in R)
            score = w[v] - penalty
            if best_v is None or (score, w[v]) > (best_score, w[best_v]):
                best_v = v
                best_score = score

        v_star = best_v
        if not any(u in chosen for u in adj[v_star]):
            chosen.add(v_star)
            S.append(v_star)

        R.remove(v_star)

    return tuple(sorted(S)), weight_of_set(S, w)

def greedy_3(n, adj, w):
    static = [0]*n
    for v in range(n):
        static[v] = w[v] - sum(w[u] for u in adj[v])

    order = sorted(range(n), key=lambda v: (static[v], w[v]), reverse=True)
    S = []
    chosen = set()

    for v in order:
        if not any(u in chosen for u in adj[v]):
            chosen.add(v)
            S.append(v)

    return tuple(sorted(S)), weight_of_set(S, w)

def run_greedies(n, adj, w):
    g1 = greedy_1(n, adj, w)
    g2 = greedy_2(n, adj, w)
    g3 = greedy_3(n, adj, w)
    return max([g1, g2, g3], key=lambda x: x[1])

# ============================================================
#  Restricted Master Problem
# ============================================================
class RestrictedMaster:
    def __init__(self, n):
        self.n = n
        self.columns = []
        self.model = None
        self.xvars = []
        self.cons = []
        self.build()

    def build(self):
        self.model = Model("CLP-r")
        self.model.setMinimize()
        self.model.setParam("display/verblevel", 0)
        self.dummy = self.model.addVar(vtype="C", lb=0.0, ub=0.0, name="dummy")

        self.cons = []
        for v in range(self.n):
            cons = self.model.addCons(self.dummy >= 1, name=f"cov_{v}")
            self.cons.append(cons)

        self.columns = []
        self.xvars = []

    def add_column(self, S):
        idx = len(self.columns)
        self.columns.append(tuple(S))
        var = self.model.addVar(vtype="C", name=f"x_{idx}")
        self.model.chgVarLb(var, 1.0)
        self.xvars.append(var)

        for v in S:
            self.model.addCoefLinear(self.cons[v], var, 1.0)

    def solve_lp(self):
        self.model.optimize()
        duals = []
        for c in self.cons:
            try:
                d = self.model.getDualsolLinear(c)
            except:
                try:
                    d = self.model.getDualsol(c)
                except:
                    d = 0.0
            duals.append(d)
        return [float(d) for d in duals]

# ============================================================
#  Column Generation Loop (NO SCALING)
# ============================================================
def column_generation(n, adj, time_limit=600, max_iter=200):

    master = RestrictedMaster(n)
    init_cols = initial_coloring_columns(n, adj)
    for S in init_cols:
        master.add_column(S)

    start = time.time()

    for it in range(1, max_iter+1):

        if time.time() - start > time_limit:
            print("Time limit reached.")
            break

        duals = master.solve_lp()
        duals = [max(0.0, min(1.0, p)) for p in duals]
        sum_pi = sum(duals)

        print(f"\nIteration {it}: LP dual sum = {sum_pi:.6f}, columns = {len(master.columns)}")

        # pricing weights directly from duals (no scaling)
        weights = {v: duals[v] for v in range(n)}

        # run heuristics
        S, wS = run_greedies(n, adj, weights)
        print(f"  Heuristic best set: |S|={len(S)}, weight={wS:.6f}")

        # violation test (no scaling):
        # column is violated when sum(pi_v for v in S) > 1
        if wS > 1.0 + 1e-12:
            print("  Violated column found → Adding.")
            master.add_column(S)
        else:
            print("  No violated column found → STOP.")
            break

    return sum_pi, master

# ============================================================
# CLI
# ============================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--time-limit", type=int, default=600)
    parser.add_argument("--max-iter", type=int, default=200)
    args = parser.parse_args()

    n, edges, adj = parse_dimacs(args.input)
    print(f"Graph loaded: n={n}, m={len(edges)}")

    dual_sum, master = column_generation(
        n, adj,
        time_limit=args.time_limit,
        max_iter=args.max_iter
    )

    print("\nFinal dual sum:", dual_sum)
    print("Total columns generated:", len(master.columns))

if __name__ == "__main__":
    main()
