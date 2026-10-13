# 05_flow_matching.py -- network flow, min cut, matching, matroid intersection.
# Python 3 stdlib only, PyPy-friendly, all iterative.
from heapq import heappush, heappop

INF = 1 << 60

# ==== 6.8 Network Flow: max flow / min cut theorem ====
# max s-t flow = min s-t cut capacity. After max flow, S-side of min cut = vertices reachable
# from s in the residual graph (Dinic.min_cut); cut edges = original edges S-side -> T-side.

# ==== 6.8 Ford-Fulkerson / Edmonds-Karp (note) ====
# FF: augment any path, O(E * maxflow). EK: BFS shortest path, O(V E^2). Dinic below is
# strictly better; never write EK. Unit-capacity graphs: Dinic O(E sqrt E); bipartite O(E sqrt V).


# ==== 6.8 Dinic's algorithm ====
class Dinic:
    """O(V^2 E) (much faster in practice). d=Dinic(n); e=d.add_edge(u,v,c[,rev_c]);
    d.max_flow(s,t); d.flow(e) flow on edge e; d.min_cut(s) -> bool list (True = S side).
    Undirected edge: add_edge(u,v,c,c). d.reset() restores original capacities."""
    def __init__(self, n):
        self.n = n
        self.g = [[] for _ in range(n)]
        self.to = []
        self.cap = []
        self.orig = []

    def add_edge(self, u, v, c, rc=0):
        e = len(self.to)
        self.g[u].append(e); self.to.append(v); self.cap.append(c)
        self.g[v].append(e + 1); self.to.append(u); self.cap.append(rc)
        self.orig += [c, rc]
        return e                                     # reverse edge is e^1

    def flow(self, e):
        return self.orig[e] - self.cap[e]

    def reset(self):
        self.cap = self.orig[:]

    def _bfs(self, s, t):
        g, to, cap = self.g, self.to, self.cap
        lv = [-1] * self.n
        lv[s] = 0
        q = [s]
        for u in q:
            for e in g[u]:
                v = to[e]
                if cap[e] > 0 and lv[v] < 0:
                    lv[v] = lv[u] + 1
                    q.append(v)
        self.lv = lv
        return lv[t] >= 0

    def max_flow(self, s, t):
        g, to, cap = self.g, self.to, self.cap
        total = 0
        while self._bfs(s, t):
            lv = self.lv
            it = [0] * self.n
            while True:                              # find one blocking-flow path at a time
                path = []
                u = s
                while u != t:
                    adj = g[u]; i = it[u]; L = len(adj); nl = lv[u] + 1
                    while i < L and not (cap[adj[i]] > 0 and lv[to[adj[i]]] == nl):
                        i += 1
                    it[u] = i
                    if i < L:
                        e = adj[i]; path.append(e); u = to[e]
                    else:                            # dead end: retreat
                        lv[u] = -1
                        if not path:
                            break
                        e = path.pop(); u = to[e ^ 1]; it[u] += 1
                if u != t:
                    break
                f = min(cap[e] for e in path)
                for e in path:
                    cap[e] -= f
                    cap[e ^ 1] += f
                total += f
        return total

    def min_cut(self, s):
        self._bfs(s, s)
        return [x >= 0 for x in self.lv]


# ==== 6.8 Push-relabel (note) ====
# Highest-label push-relabel with gap heuristic: O(V^2 sqrt E), beats Dinic on dense graphs
# in C++. In Python, Dinic is simpler and usually enough; not included.


# ==== 6.8 Min cost max flow (SPFA potentials + Dijkstra with potentials) ====
class MCMF:
    """O(F * E log V) (+ one SPFA if negative costs). Negative costs OK, no negative cycles.
    m=MCMF(n); m.add_edge(u,v,cap,cost); flow,cost = m.flow(s,t[,maxf]).
    Min cost for ANY flow amount: stop when path cost h[t] >= 0 (see comment)."""
    def __init__(self, n):
        self.n = n
        self.g = [[] for _ in range(n)]
        self.to = []; self.cap = []; self.cost = []

    def add_edge(self, u, v, c, w):
        e = len(self.to)
        self.g[u].append(e); self.to += [v]; self.cap += [c]; self.cost += [w]
        self.g[v].append(e + 1); self.to += [u]; self.cap += [0]; self.cost += [-w]
        return e

    def flow(self, s, t, maxf=INF):
        n, g, to, cap, cost = self.n, self.g, self.to, self.cap, self.cost
        h = [0] * n
        if any(cost[e] < 0 for e in range(0, len(cost), 2)):   # SPFA for initial potentials
            d = [INF] * n; d[s] = 0; inq = [False] * n; q = [s]; inq[s] = True; qi = 0
            while qi < len(q):
                u = q[qi]; qi += 1; inq[u] = False
                for e in g[u]:
                    if cap[e] > 0 and d[u] + cost[e] < d[to[e]]:
                        v = to[e]; d[v] = d[u] + cost[e]
                        if not inq[v]:
                            inq[v] = True; q.append(v)
            h = [x if x < INF else 0 for x in d]
        fl = tc = 0
        while fl < maxf:
            dist = [INF] * n; dist[s] = 0; pe = [-1] * n; pq = [(0, s)]
            while pq:
                du, u = heappop(pq)
                if du > dist[u]:
                    continue
                hu = h[u]
                for e in g[u]:
                    if cap[e] > 0:
                        v = to[e]; nd = du + cost[e] + hu - h[v]
                        if nd < dist[v]:
                            dist[v] = nd; pe[v] = e; heappush(pq, (nd, v))
            if dist[t] == INF:
                break
            for v in range(n):
                if dist[v] < INF:
                    h[v] += dist[v]
            # if h[t] - h[s] >= 0: break     # <- for min-cost (not max) flow
            f = maxf - fl; v = t
            while v != s:
                e = pe[v]; f = min(f, cap[e]); v = to[e ^ 1]
            v = t
            while v != s:
                e = pe[v]; cap[e] -= f; cap[e ^ 1] += f; v = to[e ^ 1]
            fl += f
            tc += f * (h[t] - h[s])
        return fl, tc


# ==== 6.8 Flows with lower bounds (circulations) ====
def flow_lower_bounds(n, edges, s=-1, t=-1):
    """edges: (u,v,lo,hi). s=t=-1: feasible circulation; else MAX s-t flow (net out of s).
    Returns (value, per-edge flow list) or None if infeasible. One/two Dinic runs.
    Min flow instead: after feasibility, value - max_flow(t,s) (with t-s/s-t edges removed)."""
    S, T = n, n + 1
    D = Dinic(n + 2)
    exc = [0] * n
    ids = []
    for u, v, lo, hi in edges:
        ids.append(D.add_edge(u, v, hi - lo)); exc[v] += lo; exc[u] -= lo
    need = 0
    for v in range(n):
        if exc[v] > 0:
            D.add_edge(S, v, exc[v]); need += exc[v]
        elif exc[v] < 0:
            D.add_edge(v, T, -exc[v])
    if s >= 0:
        ts = D.add_edge(t, s, INF); st = D.add_edge(s, t, INF)   # both: net may be negative
    if D.max_flow(S, T) != need:
        return None
    val = 0
    if s >= 0:
        val = D.flow(ts) - D.flow(st)
        for e in (ts, st):
            D.cap[e] = D.cap[e ^ 1] = 0
        val += D.max_flow(s, t)
    return val, [hi - D.cap[e] for e, (u, v, lo, hi) in zip(ids, edges)]


# ==== 6.8 Modeling problems as flows (notes) ====
# - Vertex capacity c(v): split v -> v_in, v_out with edge (v_in, v_out, c); in-edges go to v_in,
#   out-edges leave v_out.
# - Multiple sources/sinks: super source S->s_i, t_i->T.
# - Undirected edge: add_edge(u, v, c, c) (one pair, both directions capacity c).
# - Bipartite with capacities (assignment, b-matching): S->left (cap), left->right, right->T.
# - "At least" constraints -> lower bounds (above). Min cost with time layers -> time-expanded graph.
# - Min cut applications: cut = choose S/T side per vertex. Pay a if v in S, b if v in T:
#   s->v cap b, v->t cap a. Pay c if u in S and v in T: u->v cap c. Forbid u in S, v in T: cap INF.
#   (only "submodular" pairwise costs representable; flip one bipartite side for "same side" costs.)
# - Edge connectivity of (u,v): max flow with unit caps (undirected: cap 1 both ways).
#   Vertex connectivity (u,v non-adjacent): split each vertex with cap 1 (u,v cap INF), edges INF.
#   Global edge connectivity: Stoer-Wagner below, or n-1 flows from vertex 0.


# ==== 6.8 Project selection / maximum weight closure ====
def max_closure(w, deps):
    """w[i]: profit of taking i (may be negative). deps (u,v): taking u requires taking v.
    Returns (best total, chosen bool list). One min cut: S->i (w>0), i->T (-w), u->v INF."""
    n = len(w)
    S, T = n, n + 1
    D = Dinic(n + 2)
    pos = 0
    for i, x in enumerate(w):
        if x > 0:
            D.add_edge(S, i, x); pos += x
        elif x < 0:
            D.add_edge(i, T, -x)
    for u, v in deps:
        D.add_edge(u, v, INF)
    f = D.max_flow(S, T)
    return pos - f, D.min_cut(S)[:n]


# ==== 6.8 Gomory-Hu tree (Gusfield) ====
def gomory_hu(n, edges):
    """Undirected weighted graph -> tree [(u, parent, w)] with n-1 edges; min cut(u,v) =
    min weight on tree path. n-1 max flows. edges: (u,v,c)."""
    D = Dinic(n)
    for u, v, c in edges:
        D.add_edge(u, v, c, c)
    par = [0] * n
    tree = []
    for i in range(1, n):
        D.reset()
        f = D.max_flow(i, par[i])
        side = D.min_cut(i)
        tree.append((i, par[i], f))
        for j in range(i + 1, n):
            if par[j] == par[i] and side[j]:
                par[j] = i
    return tree


# ==== 6.8 Global min cut (Stoer-Wagner) ====
def stoer_wagner(mat):
    """O(n^3). mat: symmetric n x n weight matrix (n >= 2). Returns (cut value, one side list)."""
    n = len(mat)
    w = [r[:] for r in mat]
    alive = list(range(n))
    co = [[i] for i in range(n)]
    best, best_set = INF, []
    while len(alive) > 1:
        key = [0] * n
        inA = [False] * n
        prev = last = -1
        for _ in range(len(alive)):          # maximum adjacency ordering
            sel = -1
            for v in alive:
                if not inA[v] and (sel < 0 or key[v] > key[sel]):
                    sel = v
            inA[sel] = True
            prev, last = last, sel
            ws = w[sel]
            for v in alive:
                if not inA[v]:
                    key[v] += ws[v]
        if key[last] < best:                 # cut of the phase: {last group} vs rest
            best, best_set = key[last], co[last][:]
        co[prev] += co[last]
        for v in alive:
            w[prev][v] += w[last][v]
            w[v][prev] = w[prev][v]
        alive.remove(last)
    return best, best_set


# ==== 6.9 Bipartite matching: Kuhn (note) ====
# Kuhn: for each left u, DFS for an augmenting path, O(VE). Hopcroft-Karp below is O(E sqrt V)
# and also iterative; use it. (Greedy initial matching speeds Kuhn a lot if you do write it.)


# ==== 6.9 Hopcroft-Karp ====
def hopcroft_karp(nl, nr, adj):
    """O(E sqrt V). adj[u] = right neighbours of left u. Returns (size, ml, mr); -1 = free."""
    ml = [-1] * nl
    mr = [-1] * nr
    res = 0
    while True:
        dist = [-1] * nl
        q = [u for u in range(nl) if ml[u] < 0]
        for u in q:
            dist[u] = 0
        found = False
        for x in q:
            for y in adj[x]:
                w = mr[y]
                if w < 0:
                    found = True
                elif dist[w] < 0:
                    dist[w] = dist[x] + 1
                    q.append(w)
        if not found:
            break
        it = [0] * nl
        for u in range(nl):
            if ml[u] >= 0:
                continue
            st = [u]
            while st:
                x = st[-1]
                if it[x] < len(adj[x]):
                    y = adj[x][it[x]]; it[x] += 1
                    w = mr[y]
                    if w < 0:                        # augment along stack
                        for z in st:
                            yy = adj[z][it[z] - 1]; ml[z] = yy; mr[yy] = z
                        res += 1
                        break
                    if dist[w] == dist[x] + 1:
                        st.append(w)
                else:
                    dist[x] = -1                     # dead
                    st.pop()
    return res, ml, mr


# ==== 6.9 Konig's theorem: min vertex cover / max independent set (bipartite) ====
def konig(nl, nr, adj, ml, mr):
    """Given a MAX matching: returns (coverL, coverR, indepL, indepR) as index lists.
    |cover| = matching size; independent set = complement of cover (size n - matching)."""
    visL = [False] * nl
    visR = [False] * nr
    q = [u for u in range(nl) if ml[u] < 0]
    for u in q:
        visL[u] = True
    for x in q:                                       # alternating BFS from free left vertices
        for y in adj[x]:
            if not visR[y]:
                visR[y] = True
                w = mr[y]
                if w >= 0 and not visL[w]:
                    visL[w] = True
                    q.append(w)
    coverL = [u for u in range(nl) if not visL[u]]
    coverR = [y for y in range(nr) if visR[y]]
    indepL = [u for u in range(nl) if visL[u]]
    indepR = [y for y in range(nr) if not visR[y]]
    return coverL, coverR, indepL, indepR


# ==== 6.9 Hall's marriage theorem (note) ====
# Left side L has a perfect matching into R iff |N(S)| >= |S| for every S subset of L.
# Deficiency version: max matching = |L| - max_S(|S| - |N(S)|). A violating S = free left
# vertices + left vertices reached by Konig's alternating BFS (visL). Regular bipartite graphs
# always have perfect matchings. Use for "is it always possible" / greedy-proof arguments.


# ==== 6.9 Weighted bipartite matching: Hungarian algorithm ====
def hungarian(a):
    """O(n^2 m). a: n x m cost matrix, n <= m. Min-cost assignment of every row to a distinct
    column. Returns (cost, col) with col[i] = column of row i. Max weight: negate costs."""
    n, m = len(a), len(a[0])
    u = [0] * (n + 1); v = [0] * (m + 1); p = [0] * (m + 1); way = [0] * (m + 1)
    for i in range(1, n + 1):
        p[0] = i; j0 = 0
        minv = [INF] * (m + 1); used = [False] * (m + 1)
        while True:
            used[j0] = True
            i0 = p[j0]; ai = a[i0 - 1]; ui = u[i0]
            delta = INF; j1 = 0
            for j in range(1, m + 1):
                if not used[j]:
                    cur = ai[j - 1] - ui - v[j]
                    if cur < minv[j]:
                        minv[j] = cur; way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]; j1 = j
            for j in range(m + 1):
                if used[j]:
                    u[p[j]] += delta; v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]; p[j0] = p[j1]; j0 = j1
    col = [0] * n
    for j in range(1, m + 1):
        if p[j]:
            col[p[j] - 1] = j - 1
    return sum(a[i][col[i]] for i in range(n)), col


# ==== 6.9 General matching: Edmonds' blossom ====
def blossom(n, adj):
    """O(n^3) max cardinality matching in a general graph. adj: undirected adjacency lists.
    Returns (size, match) with match[v] = partner or -1."""
    match = [-1] * n
    for v in range(n):                               # greedy init
        if match[v] < 0:
            for u in adj[v]:
                if match[u] < 0 and u != v:
                    match[u] = v; match[v] = u
                    break

    def find_path(root):
        used = [False] * n; p = [-1] * n; base = list(range(n))
        used[root] = True
        q = [root]; qi = 0

        def lca(a, b):
            seen = [False] * n
            while True:
                a = base[a]; seen[a] = True
                if match[a] < 0:
                    break
                a = p[match[a]]
            while True:
                b = base[b]
                if seen[b]:
                    return b
                b = p[match[b]]

        def mark(v, b, child, bl):
            while base[v] != b:
                bl[base[v]] = bl[base[match[v]]] = True
                p[v] = child; child = match[v]; v = p[match[v]]

        while qi < len(q):
            v = q[qi]; qi += 1
            for to in adj[v]:
                if base[v] == base[to] or match[v] == to:
                    continue
                if to == root or (match[to] >= 0 and p[match[to]] >= 0):   # odd cycle
                    cur = lca(v, to)
                    bl = [False] * n
                    mark(v, cur, to, bl); mark(to, cur, v, bl)
                    for i in range(n):
                        if bl[base[i]]:
                            base[i] = cur
                            if not used[i]:
                                used[i] = True; q.append(i)
                elif p[to] < 0:
                    p[to] = v
                    if match[to] < 0:
                        return to, p
                    used[match[to]] = True; q.append(match[to])
        return -1, p

    for r in range(n):
        if match[r] < 0:
            v, p = find_path(r)
            while v >= 0:                            # augment
                pv = p[v]; ppv = match[pv]
                match[v] = pv; match[pv] = v; v = ppv
    return sum(1 for x in match if x >= 0) // 2, match
# Weighted general matching (O(n^3) weighted blossom) is very long -- bitmask DP for n<=20.


# ==== 6.9 Stable marriage (Gale-Shapley) ====
def gale_shapley(mpref, wpref):
    """O(n^2). mpref[m] = women in m's order, wpref[w] = men in w's order. Proposers (men)
    get their best stable partner; women their worst. Returns wife[m]."""
    n = len(mpref)
    rank = [[0] * n for _ in range(n)]
    for w in range(n):
        for r, m in enumerate(wpref[w]):
            rank[w][m] = r
    nxt = [0] * n; husband = [-1] * n; wife = [-1] * n
    free = list(range(n))
    while free:
        m = free.pop()
        w = mpref[m][nxt[m]]; nxt[m] += 1
        h = husband[w]
        if h < 0:
            husband[w] = m; wife[m] = w
        elif rank[w][m] < rank[w][h]:
            husband[w] = m; wife[m] = w; wife[h] = -1; free.append(h)
        else:
            free.append(m)
    return wife


# ==== 6.9 Minimum path cover in a DAG ====
def min_path_cover(n, edges, disjoint=True):
    """Vertex-disjoint: n - maxmatching(u_out -> v_in). disjoint=False (paths may share
    vertices): same on the transitive closure (O(n^3/64) bitsets). Returns (count, paths);
    in non-disjoint mode consecutive path vertices are reachable, not necessarily adjacent."""
    if disjoint:
        adj = [[] for _ in range(n)]
        for u, v in edges:
            adj[u].append(v)
    else:
        reach = [0] * n
        for u, v in edges:
            reach[u] |= 1 << v
        for k in range(n):
            bk = 1 << k; rk = reach[k]
            for i in range(n):
                if reach[i] & bk:
                    reach[i] |= rk
        adj = [[v for v in range(n) if reach[u] >> v & 1] for u in range(n)]
    res, ml, mr = hopcroft_karp(n, n, adj)
    paths = []
    for v in range(n):
        if mr[v] < 0:
            path = [v]
            while ml[path[-1]] >= 0:
                path.append(ml[path[-1]])
            paths.append(path)
    return n - res, paths


# ==== 6.9 Dilworth's theorem ====
# In a finite poset, min number of chains covering it = max antichain size. Min chain cover =
# non-disjoint path cover of the comparability DAG = n - matching on (i < j) pairs.
# Mirsky (dual): min antichain partition = longest chain (layer by longest-path depth).
# 2D special case: points (a,b), order a<a',b<b' -> chains via LIS-type greedy, O(n log n).
def dilworth(n, less):
    """O(n^2 sqrt n). less(i,j): STRICT partial order (must be transitive).
    Returns (chains, max_antichain); len(chains) == len(max_antichain)."""
    adj = [[j for j in range(n) if less(i, j)] for i in range(n)]
    res, ml, mr = hopcroft_karp(n, n, adj)
    chains = []
    for v in range(n):
        if mr[v] < 0:
            c = [v]
            while ml[c[-1]] >= 0:
                c.append(ml[c[-1]])
            chains.append(c)
    cL, cR, _, _ = konig(n, n, adj, ml, mr)
    inc = [False] * n
    for x in cL + cR:
        inc[x] = True
    return chains, [x for x in range(n) if not inc[x]]


# ==== 12.4 Matroid intersection ====
# Weighted (note): max weight common independent set. Same exchange graph; vertex weight
# -w(y) for y not in I, +w(x) for x in I; take the shortest path from X1 to X2 by
# (total weight, #edges) with Bellman-Ford/SPFA (no negative cycles). After k augmentations
# I is a max-weight common independent set of size k; take the best over k. O(r^2 n) oracles + BF.
def matroid_intersection(n, ind1, ind2):
    """Max-cardinality common independent set of two matroids on elements 0..n-1.
    ind1/ind2(list of elements) -> bool. O(r^2 * n) oracle calls (r = answer size).
    Examples: bipartite matching = partition(left) x partition(right);
    colorful spanning forest = graphic x partition(color)."""
    inI = [False] * n
    while True:
        I = [i for i in range(n) if inI[i]]
        out = [i for i in range(n) if not inI[i]]
        X2 = [False] * n
        for y in out:
            X2[y] = ind2(I + [y])
        adj = [[] for _ in range(n)]
        for x in I:
            J = [z for z in I if z != x]
            for y in out:
                if ind1(J + [y]):
                    adj[x].append(y)
                if ind2(J + [y]):
                    adj[y].append(x)
        par = [-2] * n
        q = []
        for y in out:
            if ind1(I + [y]):
                par[y] = -1; q.append(y)
        end = -1
        for u in q:                                  # BFS: shortest path X1 -> X2
            if X2[u]:
                end = u
                break
            for v in adj[u]:
                if par[v] == -2:
                    par[v] = u; q.append(v)
        if end < 0:
            return I
        while end >= 0:
            inI[end] = not inI[end]; end = par[end]


# =============================== tests ===============================
if __name__ == "__main__":
    import random
    from itertools import permutations, product
    rnd = random.Random(12345)

    # --- Dinic vs brute min cut ---
    for _ in range(1000):
        n = rnd.randint(2, 7)
        E = [(rnd.randrange(n), rnd.randrange(n), rnd.randint(0, 9)) for _ in range(rnd.randint(0, 14))]
        s, t = rnd.sample(range(n), 2)
        D = Dinic(n)
        ids = [D.add_edge(u, v, c) for u, v, c in E]
        f = D.max_flow(s, t)
        best = min(sum(c for u, v, c in E if m >> u & 1 and not m >> v & 1)
                   for m in range(1 << n) if m >> s & 1 and not m >> t & 1)
        assert f == best
        side = D.min_cut(s)
        assert side[s] and not side[t]
        assert sum(c for u, v, c in E if side[u] and not side[v]) == f
        for e, (u, v, c) in zip(ids, E):
            assert 0 <= D.flow(e) <= c

    # --- brute flow enumerator: edges (u,v,lo,hi[,cost]) ---
    def brute(n, E, s, t):
        res = []
        for fl in product(*[range(e[2], e[3] + 1) for e in E]):
            bal = [0] * n
            for e, x in zip(E, fl):
                bal[e[0]] -= x; bal[e[1]] += x
            if all(bal[v] == 0 for v in range(n) if v != s and v != t):
                res.append((-bal[s] if s >= 0 else 0, sum(x * (e[4] if len(e) > 4 else 0) for e, x in zip(E, fl))))
        return res

    # --- MCMF vs brute (DAG, negative costs) ---
    for _ in range(200):
        n = rnd.randint(2, 5)
        E = []
        for _ in range(rnd.randint(0, 6)):
            u, v = sorted(rnd.sample(range(n), 2))
            E.append((u, v, 0, rnd.randint(0, 2), rnd.randint(-5, 5)))
        M = MCMF(n)
        for u, v, lo, hi, c in E:
            M.add_edge(u, v, hi, c)
        fl, cost = M.flow(0, n - 1)
        b = brute(n, E, 0, n - 1)
        mf = max(x for x, _ in b)
        assert fl == mf and cost == min(c for x, c in b if x == mf)
        D = Dinic(n)
        for u, v, lo, hi, c in E:
            D.add_edge(u, v, hi)
        assert D.max_flow(0, n - 1) == fl

    # --- lower bounds vs brute (general graphs) ---
    for _ in range(300):
        n = rnd.randint(2, 4)
        E = []
        for _ in range(rnd.randint(1, 6)):
            u, v = rnd.sample(range(n), 2)
            lo = rnd.randint(0, 2)
            E.append((u, v, lo, lo + rnd.randint(0, 2)))
        circ = rnd.random() < 0.3
        s, t = (-1, -1) if circ else rnd.sample(range(n), 2)
        r = flow_lower_bounds(n, E, s, t)
        b = brute(n, E, s, t)
        if not b:
            assert r is None
        else:
            assert r is not None and r[0] == max(x for x, _ in b)
            val, fls = r
            bal = [0] * n
            for (u, v, lo, hi), x in zip(E, fls):
                assert lo <= x <= hi
                bal[u] -= x; bal[v] += x
            assert all(bal[v] == 0 for v in range(n) if v != s and v != t)
            if s >= 0:
                assert -bal[s] == val

    # --- max closure vs brute ---
    for _ in range(200):
        n = rnd.randint(1, 8)
        w = [rnd.randint(-10, 10) for _ in range(n)]
        deps = [tuple(rnd.sample(range(n), 2)) for _ in range(rnd.randint(0, 10))] if n > 1 else []
        val, ch = max_closure(w, deps)
        ok = lambda m: all(not (m >> u & 1) or (m >> v & 1) for u, v in deps)
        best = max(sum(w[i] for i in range(n) if m >> i & 1) for m in range(1 << n) if ok(m))
        assert val == best
        assert all(not ch[u] or ch[v] for u, v in deps) and sum(w[i] for i in range(n) if ch[i]) == val

    # --- Gomory-Hu vs pairwise max flow; Stoer-Wagner vs brute ---
    for _ in range(100):
        n = rnd.randint(2, 7)
        E = [(*rnd.sample(range(n), 2), rnd.randint(1, 9)) for _ in range(rnd.randint(0, 12))]
        tree = gomory_hu(n, E)
        tadj = [[] for _ in range(n)]
        for u, v, c in tree:
            tadj[u].append((v, c)); tadj[v].append((u, c))
        for a in range(n):
            mn = [None] * n; mn[a] = INF; st = [a]
            while st:
                x = st.pop()
                for y, c in tadj[x]:
                    if mn[y] is None:
                        mn[y] = min(mn[x], c); st.append(y)
            for b in range(a + 1, n):
                D = Dinic(n)
                for u, v, c in E:
                    D.add_edge(u, v, c, c)
                assert mn[b] == D.max_flow(a, b)
        mat = [[0] * n for _ in range(n)]
        for u, v, c in E:
            mat[u][v] += c; mat[v][u] += c
        cut = lambda m: sum(c for u, v, c in E if (m >> u & 1) != (m >> v & 1))
        val, side = stoer_wagner(mat)
        assert val == min(cut(m) for m in range(1, (1 << n) - 1))
        assert 0 < len(side) < n and cut(sum(1 << x for x in side)) == val

    # --- Hopcroft-Karp + Konig vs brute vertex cover ---
    for _ in range(300):
        nl, nr = rnd.randint(1, 6), rnd.randint(1, 6)
        adj = [[y for y in range(nr) if rnd.random() < 0.4] for _ in range(nl)]
        res, ml, mr = hopcroft_karp(nl, nr, adj)
        assert all(ml[u] < 0 or (ml[u] in adj[u] and mr[ml[u]] == u) for u in range(nl))
        assert sum(1 for u in range(nl) if ml[u] >= 0) == res
        E = [(u, y) for u in range(nl) for y in adj[u]]
        bestvc = min(bin(m).count("1") for m in range(1 << (nl + nr))
                     if all(m >> u & 1 or m >> (nl + y) & 1 for u, y in E))
        assert res == bestvc
        cL, cR, iL, iR = konig(nl, nr, adj, ml, mr)
        assert len(cL) + len(cR) == res
        sL, sR = set(cL), set(cR)
        assert all(u in sL or y in sR for u, y in E)
        jL, jR = set(iL), set(iR)
        assert len(iL) + len(iR) == nl + nr - res and not any(u in jL and y in jR for u, y in E)

    # --- Hungarian vs brute ---
    for _ in range(300):
        n = rnd.randint(1, 5); m = rnd.randint(n, 6)
        a = [[rnd.randint(-20, 20) for _ in range(m)] for _ in range(n)]
        cost, col = hungarian(a)
        assert len(set(col)) == n and sum(a[i][col[i]] for i in range(n)) == cost
        assert cost == min(sum(a[i][p[i]] for i in range(n)) for p in permutations(range(m), n))

    # --- Blossom vs bitmask DP ---
    for _ in range(300):
        n = rnd.randint(1, 13)
        p = rnd.random()
        E = [(u, v) for u in range(n) for v in range(u + 1, n) if rnd.random() < p]
        adj = [[] for _ in range(n)]
        for u, v in E:
            adj[u].append(v); adj[v].append(u)
        sz, mt = blossom(n, adj)
        assert all(mt[v] < 0 or (mt[mt[v]] == v and mt[v] in adj[v]) for v in range(n))
        assert sum(1 for x in mt if x >= 0) == 2 * sz
        am = [0] * n
        for u, v in E:
            am[u] |= 1 << v; am[v] |= 1 << u
        dp = [0] * (1 << n)
        for m in range(1, 1 << n):
            i = (m & -m).bit_length() - 1
            r = m & ~(1 << i)
            best = dp[r]
            c = am[i] & r
            while c:
                j = (c & -c).bit_length() - 1; c &= c - 1
                best = max(best, 1 + dp[r & ~(1 << j)])
            dp[m] = best
        assert sz == dp[(1 << n) - 1]

    # --- Gale-Shapley: stable + man-optimal vs brute ---
    for _ in range(150):
        n = rnd.randint(1, 5)
        mp = [rnd.sample(range(n), n) for _ in range(n)]
        wp = [rnd.sample(range(n), n) for _ in range(n)]
        mr_ = [[0] * n for _ in range(n)]; wr = [[0] * n for _ in range(n)]
        for i in range(n):
            for r, x in enumerate(mp[i]):
                mr_[i][x] = r
            for r, x in enumerate(wp[i]):
                wr[i][x] = r

        def stable(wife):
            hus = [0] * n
            for m, w in enumerate(wife):
                hus[w] = m
            return not any(mr_[m][w] < mr_[m][wife[m]] and wr[w][m] < wr[w][hus[w]]
                           for m in range(n) for w in range(n))
        wife = gale_shapley(mp, wp)
        assert sorted(wife) == list(range(n)) and stable(wife)
        st = [p for p in permutations(range(n)) if stable(p)]
        for m in range(n):
            assert mr_[m][wife[m]] == min(mr_[m][p[m]] for p in st)

    # --- path cover (disjoint vs brute edge subsets; non-disjoint vs brute antichain) ---
    for _ in range(200):
        n = rnd.randint(1, 7)
        E = list({tuple(sorted(rnd.sample(range(n), 2))) for _ in range(rnd.randint(0, 10))}) if n > 1 else []
        cnt, paths = min_path_cover(n, E)
        Es = set(E)
        assert sorted(v for p in paths for v in p) == list(range(n)) and len(paths) == cnt
        assert all((p[i], p[i + 1]) in Es for p in paths for i in range(len(p) - 1))
        best = 0
        for m in range(1 << len(E)):
            sel = [E[i] for i in range(len(E)) if m >> i & 1]
            if len({u for u, v in sel}) == len(sel) and len({v for u, v in sel}) == len(sel):
                best = max(best, len(sel))
        assert cnt == n - best
        reach = [[False] * n for _ in range(n)]
        for u, v in E:
            reach[u][v] = True
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    reach[i][j] |= reach[i][k] and reach[k][j]
        cnt2, paths2 = min_path_cover(n, E, disjoint=False)
        anti = max(bin(m).count("1") for m in range(1 << n)
                   if not any(m >> i & 1 and m >> j & 1 and reach[i][j] for i in range(n) for j in range(n)))
        assert cnt2 == anti
        assert all(reach[p[i]][p[i + 1]] for p in paths2 for i in range(len(p) - 1))

    # --- Dilworth on 2D dominance posets ---
    for _ in range(200):
        n = rnd.randint(1, 9)
        P = [(rnd.randint(0, 5), rnd.randint(0, 5)) for _ in range(n)]
        less = lambda i, j: P[i][0] < P[j][0] and P[i][1] < P[j][1]
        chains, anti = dilworth(n, less)
        assert len(chains) == len(anti)
        assert sorted(v for c in chains for v in c) == list(range(n))
        assert all(less(c[i], c[i + 1]) for c in chains for i in range(len(c) - 1))
        assert not any(less(i, j) for i in anti for j in anti)
        best = max(bin(m).count("1") for m in range(1 << n)
                   if not any(m >> i & 1 and m >> j & 1 and less(i, j) for i in range(n) for j in range(n)))
        assert len(anti) == best

    # --- matroid intersection: colorful spanning forest vs brute; bipartite matching vs HK ---
    def graphic(V, E):
        def ind(S):
            p = list(range(V))
            def f(x):
                while p[x] != x:
                    p[x] = p[p[x]]; x = p[x]
                return x
            for i in S:
                a, b = f(E[i][0]), f(E[i][1])
                if a == b:
                    return False
                p[a] = b
            return True
        return ind

    def partition(lab):
        return lambda S: len({lab[i] for i in S}) == len(S)

    for _ in range(150):
        V = rnd.randint(1, 6); m = rnd.randint(0, 10)
        E = [(rnd.randrange(V), rnd.randrange(V)) for _ in range(m)]
        col = [rnd.randrange(4) for _ in range(m)]
        g1, p2 = graphic(V, E), partition(col)
        I = matroid_intersection(m, g1, p2)
        assert g1(I) and p2(I)
        best = max(bin(s).count("1") for s in range(1 << m)
                   if g1([i for i in range(m) if s >> i & 1]) and p2([i for i in range(m) if s >> i & 1]))
        assert len(I) == best
    for _ in range(100):
        nl, nr = rnd.randint(1, 6), rnd.randint(1, 6)
        adj = [[y for y in range(nr) if rnd.random() < 0.4] for _ in range(nl)]
        E = [(u, y) for u in range(nl) for y in adj[u]]
        I = matroid_intersection(len(E), partition([u for u, y in E]), partition([y for u, y in E]))
        assert len(I) == hopcroft_karp(nl, nr, adj)[0]

    print("ALL OK")
