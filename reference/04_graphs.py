# 04_graphs.py -- Section 6.1-6.7, 6.10 graph algorithms (flow/matching: 05, trees/LCA: 06).
# Python 3 stdlib only, PyPy-friendly, NO recursion. Vertices 0..n-1.
# Conventions: g[u] = [v, ...] (unweighted) or [(v, w), ...] (weighted); edges = [(u, v)] or [(u, v, w)].
import random
from heapq import heappush, heappop
from collections import deque
from itertools import combinations, permutations, product

INF = float('inf')  # used by snippets below; ints compare fine against it

# ==== 6.1 Graph basics ====
# Adjacency list: g = [[] for _ in range(n)]; g[u].append(v) (+ g[v].append(u) if undirected).
# Matrix: a = [[0]*n for _ in range(n)] (n<=~2000); bitmask rows adj[u] |= 1<<v for n<=~60 tricks.
# Edge list [(u,v,w)] for Kruskal/Bellman-Ford. Implicit graphs (grid, states): generate
# neighbours on the fly; encode state as int (r*W+c, mask*n+v) for flat dist arrays.
# Degrees: sum(deg) = 2m (undirected). Reading fast: data = sys.stdin.buffer.read().split().

def csr(n, edges):  # O(n+m) compressed adjacency: neighbours of u are to[start[u]:start[u+1]]
    start = [0] * (n + 1)  # faster/leaner than list-of-lists for n,m ~ 1e6 (PyPy)
    for u, v in edges:
        start[u + 1] += 1
    for i in range(n):
        start[i + 1] += start[i]
    pos = start[:]; to = [0] * len(edges)
    for u, v in edges:
        to[pos[u]] = v; pos[u] += 1
    return start, to

def bipartite(g):  # O(n+m); returns color list (0/1) or None if odd cycle
    n = len(g); c = [-1] * n
    for s in range(n):
        if c[s] >= 0: continue
        c[s] = 0; st = [s]
        while st:
            u = st.pop()
            for v in g[u]:
                if c[v] < 0: c[v] = c[u] ^ 1; st.append(v)
                elif c[v] == c[u]: return None
    return c

# ==== 6.2 Graph traversal ====
def dfs_times(g):  # O(n+m) iterative DFS on directed g. kind[u][i] for edge g[u][i]:
    n = len(g)     # 0 tree, 1 back (incl. self-loop), 2 forward, 3 cross. Undirected: only 0/1.
    tin = [-1] * n; tout = [-1] * n; par = [-1] * n; it = [0] * n; col = [0] * n
    kind = [[0] * len(g[u]) for u in range(n)]; t = 0
    for r in range(n):
        if col[r]: continue
        tin[r] = t; t += 1; col[r] = 1; st = [r]
        while st:
            u = st[-1]
            if it[u] < len(g[u]):
                i = it[u]; v = g[u][i]; it[u] += 1
                if col[v] == 0:
                    par[v] = u; tin[v] = t; t += 1; col[v] = 1; st.append(v)
                else:
                    kind[u][i] = 1 if col[v] == 1 else (2 if tin[v] > tin[u] else 3)
            else:
                st.pop(); col[u] = 2; tout[u] = t; t += 1
    return tin, tout, par, kind   # v in subtree(u) <=> tin[u] <= tin[v] and tout[v] <= tout[u]
# DFS topological order = vertices sorted by tout descending (DAG).

def bfs(g, srcs):  # O(n+m) unweighted dist from a list of sources (multi-source); -1 = unreachable
    n = len(g); d = [-1] * n; par = [-1] * n; q = list(srcs)
    for s in q: d[s] = 0
    for u in q:  # list as queue (appending while iterating is fine)
        for v in g[u]:
            if d[v] < 0: d[v] = d[u] + 1; par[v] = u; q.append(v)
    return d, par

def bfs01(g, s):  # O(n+m); g[u]=[(v,w)] with w in {0,1}
    d = [INF] * len(g); d[s] = 0; dq = deque([s])
    while dq:
        u = dq.popleft(); du = d[u]
        for v, w in g[u]:
            if du + w < d[v]:
                d[v] = du + w
                if w: dq.append(v)
                else: dq.appendleft(v)
    return d

def flood_fill(grid, free='.'):  # O(RC) 4-dir labelling; returns (lab, cnt), lab=-1 for walls
    R = len(grid); C = len(grid[0]); lab = [[-1] * C for _ in range(R)]; cnt = 0
    for i in range(R):
        for j in range(C):
            if grid[i][j] != free or lab[i][j] >= 0: continue
            lab[i][j] = cnt; st = [(i, j)]
            while st:
                r, c = st.pop()
                for nr, nc in ((r + 1, c), (r - 1, c), (r, c + 1), (r, c - 1)):
                    if 0 <= nr < R and 0 <= nc < C and lab[nr][nc] < 0 and grid[nr][nc] == free:
                        lab[nr][nc] = cnt; st.append((nr, nc))
            cnt += 1
    return lab, cnt

def components(g):  # O(n+m) undirected; returns (comp id list, count)
    n = len(g); comp = [-1] * n; k = 0
    for s in range(n):
        if comp[s] >= 0: continue
        comp[s] = k; st = [s]
        while st:
            u = st.pop()
            for v in g[u]:
                if comp[v] < 0: comp[v] = k; st.append(v)
        k += 1
    return comp, k

def find_cycle_directed(g):  # O(n+m); returns [v0..vk-1] with v_i->v_{i+1}, v_{k-1}->v0, or None
    n = len(g); col = [0] * n; par = [-1] * n; it = [0] * n
    for r in range(n):
        if col[r]: continue
        col[r] = 1; st = [r]
        while st:
            u = st[-1]
            if it[u] < len(g[u]):
                v = g[u][it[u]]; it[u] += 1
                if col[v] == 0: col[v] = 1; par[v] = u; st.append(v)
                elif col[v] == 1:
                    cyc = [u]
                    while cyc[-1] != v: cyc.append(par[cyc[-1]])
                    return cyc[::-1]
            else:
                col[u] = 2; st.pop()
    return None

def find_cycle_undirected(n, edges):  # O(n+m); multi-edges/self-loops ok (edge ids); vertex list or None
    g = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        g[u].append((v, i)); g[v].append((u, i))
    col = [0] * n; par = [-1] * n; pe = [-1] * n; it = [0] * n
    for r in range(n):
        if col[r]: continue
        col[r] = 1; st = [r]
        while st:
            u = st[-1]
            if it[u] < len(g[u]):
                v, e = g[u][it[u]]; it[u] += 1
                if e == pe[u]: continue
                if col[v] == 0: col[v] = 1; par[v] = u; pe[v] = e; st.append(v)
                elif col[v] == 1:  # back edge to ancestor
                    cyc = [u]
                    while cyc[-1] != v: cyc.append(par[cyc[-1]])
                    return cyc
            else:
                col[u] = 2; st.pop()
    return None

# ==== 6.3 Shortest paths ====
def dijkstra(g, s):  # O(m log m) lazy heap; g[u]=[(v,w)], w>=0; returns (dist, parent)
    n = len(g); d = [INF] * n; p = [-1] * n; d[s] = 0; pq = [(0, s)]
    while pq:  # speed: push single int du*n+u instead of tuples (if dist ints bounded)
        du, u = heappop(pq)
        if du > d[u]: continue
        for v, w in g[u]:
            nd = du + w
            if nd < d[v]: d[v] = nd; p[v] = u; heappush(pq, (nd, v))
    return d, p

def get_path(p, t):  # path root..t from parent array (check t reachable first)
    path = []
    while t != -1: path.append(t); t = p[t]
    return path[::-1]

def dijkstra_dense(W, s):  # O(n^2) for dense graphs; W[u][v] = weight or INF
    n = len(W); d = [INF] * n; d[s] = 0; done = [False] * n
    for _ in range(n):
        u = -1; best = INF
        for v in range(n):
            if not done[v] and d[v] < best: best = d[v]; u = v
        if u < 0: break
        done[u] = True; row = W[u]
        for v in range(n):
            if best + row[v] < d[v]: d[v] = best + row[v]
    return d

def bellman_ford(n, edges, s):  # O(nm); edges (u,v,w). returns (dist, neg_cycle or None)
    d = [INF] * n; p = [-1] * n; d[s] = 0  # cycle only if reachable from s; it's a vertex list
    for _ in range(n):                     # with edges cyc[i]->cyc[i+1] (cyclic)
        x = -1
        for u, v, w in edges:
            if d[u] + w < d[v]: d[v] = d[u] + w; p[v] = u; x = v
        if x < 0: return d, None
    for _ in range(n): x = p[x]  # step back n times -> surely on the cycle
    cyc = [x]; y = p[x]
    while y != x: cyc.append(y); y = p[y]
    return d, cyc[::-1]
# Vertices with dist -inf: everything reachable from any vertex relaxed in the n-th pass.
# Detect ANY negative cycle (not only reachable): start with d=[0]*n (virtual source).
# SPFA = Bellman-Ford with a queue of changed vertices. Avg fast, worst O(nm) and setters
# hack it (grid + random weights). Prefer Dijkstra (+potentials) whenever weights allow.

def floyd(d):  # O(n^3) in place; d[i][i]=0, INF if no edge. Neg cycle through i iff d[i][i] < 0
    n = len(d)
    for k in range(n):
        dk = d[k]
        for i in range(n):
            di = d[i]; dik = di[k]
            if dik == INF: continue
            for j in range(n):
                if dik + dk[j] < di[j]: di[j] = dik + dk[j]
    return d

def transitive_closure(adj):  # O(n^3/64) bitsets; adj[i] mask; r[i] bit j <=> path i->j (len>=1)
    r = adj[:]; n = len(r)
    for k in range(n):
        bk = 1 << k; rk = r[k]
        for i in range(n):
            if r[i] & bk: r[i] |= rk
    return r

def johnson(n, edges):  # O(nm log m) APSP with negative edges; matrix (INF unreachable) or None
    h, cyc = bellman_ford(n + 1, edges + [(n, v, 0) for v in range(n)], n)  # uses bellman_ford, dijkstra
    if cyc: return None
    g = [[] for _ in range(n)]
    for u, v, w in edges: g[u].append((v, w + h[u] - h[v]))  # reduced weights >= 0
    res = []
    for s in range(n):
        d, _ = dijkstra(g, s)
        res.append([d[v] - h[s] + h[v] if d[v] < INF else INF for v in range(n)])
    return res

def dag_shortest(g, s):  # O(n+m) any weights on DAG; g[u]=[(v,w)]; longest path: negate weights
    n = len(g); indeg = [0] * n
    for u in range(n):
        for v, _ in g[u]: indeg[v] += 1
    order = [u for u in range(n) if indeg[u] == 0]
    for u in order:
        for v, _ in g[u]:
            indeg[v] -= 1
            if indeg[v] == 0: order.append(v)
    d = [INF] * n; d[s] = 0
    for u in order:
        if d[u] < INF:
            for v, w in g[u]:
                if d[u] + w < d[v]: d[v] = d[u] + w
    return d

def k_shortest_walks(g, s, t, k):  # O(k m log(km)); k smallest s->t WALK lengths (repeats allowed)
    cnt = [0] * len(g); res = []; pq = [(0, s)]  # each vertex popped at most k times; w>=0
    while pq and len(res) < k:
        d, u = heappop(pq)
        if cnt[u] >= k: continue
        cnt[u] += 1
        if u == t: res.append(d)
        for v, w in g[u]:
            if cnt[v] < k: heappush(pq, (d + w, v))
    return res
# Simple (loopless) k-shortest paths: Yen = for each prefix of previous path, ban its next edges
#   + prefix vertices, run Dijkstra for spur; O(k n (m log n)). Eppstein: O(m + n log n + k log k)
#   walks via sidetrack-edge heaps on the shortest-path tree (persistent leftist heaps).
# Layered/state graph: vertex = (v, extra) e.g. (v, coupons used), (v, parity), (v, time mod k);
#   index id = v*K + extra, run Dijkstra/BFS on the product graph (n*K vertices).

def diff_constraints(n, cons):  # cons (u,v,c) meaning x_v - x_u <= c; returns x or None. O(nm)
    d, cyc = bellman_ford(n + 1, list(cons) + [(n, i, 0) for i in range(n)], n)
    return None if cyc else d[:n]  # solution has max(x)<=0; shift freely. x_v>=x_u+c: (v,u,-c)
# Tight-ness: d gives the max x with x_i <= 0; min of x_t-x_s subject to constraints: dist s->t.

# ==== 6.4 Minimum spanning trees ====
def kruskal(n, edges):  # O(m log m); edges (u,v,w); returns (weight, used edge ids); forest if disconnected
    p = list(range(n))  # minimal inline DSU (full DSU in 02)
    def find(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    tot = 0; used = []
    for i in sorted(range(len(edges)), key=lambda i: edges[i][2]):
        u, v, w = edges[i]; a, b = find(u), find(v)
        if a != b: p[a] = b; tot += w; used.append(i)
    return tot, used

def prim(g):  # O(m log m); g[u]=[(v,w)] undirected connected; returns weight or None if disconnected
    n = len(g); seen = [False] * n; tot = 0; cnt = 0; pq = [(0, 0)]
    while pq:
        w, u = heappop(pq)
        if seen[u]: continue
        seen[u] = True; tot += w; cnt += 1
        for v, wv in g[u]:
            if not seen[v]: heappush(pq, (wv, v))
    return tot if cnt == n else None  # dense graph: O(n^2) array version like dijkstra_dense

def boruvka(n, edges):  # O(m log n); edges (u,v,w); returns (weight, used ids). ties broken by id
    p = list(range(n))
    def find(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    tot = 0; used = []
    while True:
        best = [-1] * n
        for i, (u, v, w) in enumerate(edges):
            a, b = find(u), find(v)
            if a == b: continue
            for c in (a, b):
                j = best[c]
                if j < 0 or (w, i) < (edges[j][2], j): best[c] = i
        merged = False
        for i in best:
            if i < 0: continue
            u, v, w = edges[i]; a, b = find(u), find(v)
            if a != b: p[a] = b; tot += w; used.append(i); merged = True
        if not merged: return tot, used

def second_mst(n, edges):  # O(n^2 + m log m) STRICTLY second-best spanning tree weight, or None
    tot, used = kruskal(n, edges)  # uses kruskal. large n: LCA lifting storing (max, strict 2nd max)
    if len(used) < n - 1: return None
    inT = [False] * len(edges); t = [[] for _ in range(n)]; other = [[] for _ in range(n)]
    for i in used:
        u, v, w = edges[i]; inT[i] = True; t[u].append((v, w)); t[v].append((u, w))
    for i, (u, v, w) in enumerate(edges):
        if not inT[i] and u != v: other[u].append((v, w))
    best = INF
    for s in range(n):
        if not other[s]: continue
        m1 = [-INF] * n; m2 = [-INF] * n; seen = [False] * n; seen[s] = True; st = [s]
        while st:  # m1/m2 = max / strictly smaller second max on tree path s..x
            u = st.pop()
            for v, w in t[u]:
                if not seen[v]:
                    seen[v] = True; a, b = m1[u], m2[u]
                    if w > a: a, b = w, a
                    elif a > w > b: b = w
                    m1[v], m2[v] = a, b; st.append(v)
        for x, w in other[s]:
            if w > m1[x]: best = min(best, w - m1[x])
            elif m2[x] > -INF: best = min(best, w - m2[x])
    return None if best == INF else tot + best
# Minimum bottleneck spanning tree: any MST is one. Minimax path u-v (min possible max edge)
#   = max edge on MST path = val[lca] in the Kruskal reconstruction tree below.

def kruskal_tree(n, edges):  # O(m log m) Kruskal reconstruction tree. leaves 0..n-1, internal n..
    p = list(range(2 * n)); par = [-1] * (2 * n - 1); val = [0] * (2 * n - 1); x = n
    def find(a):
        while p[a] != a: p[a] = p[p[a]]; a = p[a]
        return a
    for u, v, w in sorted(edges, key=lambda e: e[2]):
        a, b = find(u), find(v)
        if a != b:
            p[a] = p[b] = x; par[a] = par[b] = x; val[x] = w; x += 1
    return par[:x], val[:x]  # val increases toward root (a heap). minimax(u,v) = val[lca(u,v)];
# vertices reachable from u via edges <= X = leaves under highest ancestor a of u with val[a] <= X
# (binary lifting, see 06). Max-spanning variant: sort descending.

def directed_mst(n, edges, root):  # Chu-Liu/Edmonds O(nm); edges (u,v,w) u->v; min arborescence or None
    res = 0
    while True:
        inw = [INF] * n; pre = [-1] * n
        for u, v, w in edges:
            if u != v and w < inw[v]: inw[v] = w; pre[v] = u
        for v in range(n):
            if v != root and inw[v] == INF: return None
        inw[root] = 0; cid = [-1] * n; vis = [-1] * n; cnt = 0
        for v in range(n):
            res += inw[v]; x = v
            while vis[x] != v and cid[x] == -1 and x != root: vis[x] = v; x = pre[x]
            if x != root and cid[x] == -1:  # found a new cycle through x
                y = pre[x]
                while y != x: cid[y] = cnt; y = pre[y]
                cid[x] = cnt; cnt += 1
        if cnt == 0: return res
        for v in range(n):
            if cid[v] == -1: cid[v] = cnt; cnt += 1
        edges = [(cid[u], cid[v], w - inw[v]) for u, v, w in edges if cid[u] != cid[v]]
        n = cnt; root = cid[root]
# Edge recovery needs bookkeeping per contraction; O(m log n) version uses mergeable heaps + DSU.

def steiner_tree(g, terms):  # O(3^k n + 2^k m log m); min weight tree connecting terminals
    n = len(g); k = len(terms); full = (1 << k) - 1
    dp = [[INF] * n for _ in range(1 << k)]
    for i, t in enumerate(terms): dp[1 << i][t] = 0
    for S in range(1, full + 1):
        d = dp[S]; low = S & -S; sub = (S - 1) & S
        while sub:
            if sub & low:  # each split once
                a = dp[sub]; b = dp[S ^ sub]
                for v in range(n):
                    if a[v] + b[v] < d[v]: d[v] = a[v] + b[v]
            sub = (sub - 1) & S
        pq = [(d[v], v) for v in range(n) if d[v] < INF]; pq.sort()
        while pq:
            dv, v = heappop(pq)
            if dv > d[v]: continue
            for u, w in g[v]:
                if dv + w < d[u]: d[u] = dv + w; heappush(pq, (d[u], u))
    return min(dp[full]) if k else 0

# ==== 6.5 Directed graphs ====
def toposort(g, lex=False):  # Kahn O(n+m) (lex=True: smallest-first heap, O(m log n)); None if cycle
    n = len(g); indeg = [0] * n
    for u in range(n):
        for v in g[u]: indeg[v] += 1
    q = [u for u in range(n) if indeg[u] == 0]; order = []
    if lex:
        while q:
            u = heappop(q); order.append(u)  # q already sorted -> valid heap
            for v in g[u]:
                indeg[v] -= 1
                if indeg[v] == 0: heappush(q, v)
    else:
        order = q
        for u in order:
            for v in g[u]:
                indeg[v] -= 1
                if indeg[v] == 0: order.append(v)
    return order if len(order) == n else None
# Lexicographically smallest order where small LABELS come early (not same!): reverse graph,
# pop largest first via max-heap, reverse result.

def scc(g):  # iterative Tarjan O(n+m); returns (k, comp); comp ids in TOPOLOGICAL order (u->v => comp[u]<=comp[v])
    n = len(g); idx = [-1] * n; low = [0] * n; onst = [False] * n; st = []
    comp = [-1] * n; it = [0] * n; t = 0; k = 0
    for s in range(n):
        if idx[s] >= 0: continue
        idx[s] = low[s] = t; t += 1; st.append(s); onst[s] = True; cs = [s]
        while cs:
            u = cs[-1]
            if it[u] < len(g[u]):
                v = g[u][it[u]]; it[u] += 1
                if idx[v] < 0:
                    idx[v] = low[v] = t; t += 1; st.append(v); onst[v] = True; cs.append(v)
                elif onst[v] and idx[v] < low[u]: low[u] = idx[v]
            else:
                cs.pop()
                if cs and low[u] < low[cs[-1]]: low[cs[-1]] = low[u]
                if low[u] == idx[u]:
                    while True:
                        x = st.pop(); onst[x] = False; comp[x] = k
                        if x == u: break
                    k += 1
    return k, [k - 1 - c for c in comp]
# Kosaraju alternative: order by DFS exit time, then DFS on reversed graph in decreasing exit.

def condensation(g, k, comp):  # DAG on k SCC nodes (deduplicated edges), ids already topo-sorted
    seen = set(); dag = [[] for _ in range(k)]
    for u in range(len(g)):
        for v in g[u]:
            a, b = comp[u], comp[v]
            if a != b and (a, b) not in seen: seen.add((a, b)); dag[a].append(b)
    return dag

class TwoSat:  # O(n+m). literal: 2*i = x_i, 2*i+1 = not x_i (neg = lit^1). uses scc above
    def __init__(self, n): self.n = n; self.g = [[] for _ in range(2 * n)]
    def either(self, a, b):  # clause (a or b)
        self.g[a ^ 1].append(b); self.g[b ^ 1].append(a)
    def implies(self, a, b): self.either(a ^ 1, b)
    def force(self, a): self.either(a, a)
    def solve(self):  # list of bools or None
        _, comp = scc(self.g); res = []
        for i in range(self.n):
            if comp[2 * i] == comp[2 * i + 1]: return None
            res.append(comp[2 * i] > comp[2 * i + 1])  # later in topo order => true
        return res
# at-most-one of many literals: prefix-variable trick (3n clauses) instead of O(n^2) pairs.

def functional_graph(f):  # O(n); f[v] = successor. returns (cycles, cyc_id, dist_to_cycle, entry)
    n = len(f); state = [0] * n; cid = [-1] * n; dist = [0] * n; entry = [-1] * n; cycles = []
    for s in range(n):
        if state[s]: continue
        path = []; x = s
        while state[x] == 0: state[x] = 1; path.append(x); x = f[x]
        if state[x] == 1:  # closed a new cycle starting at x
            i = path.index(x); cyc = path[i:]
            for v in cyc: cid[v] = len(cycles); entry[v] = v
            cycles.append(cyc); tail = path[:i]
        else: tail = path
        for v in reversed(tail): dist[v] = dist[f[v]] + 1; entry[v] = entry[f[v]]
        for v in path: state[v] = 2
    return cycles, cid, dist, entry  # k-th successor: binary lifting, or via tail+cycle index.
# Floyd (tortoise-hare) / Brent: cycle of a single sequence x->f(x) in O(1) memory.

def dominator_tree(g, root):  # Lengauer-Tarjan O(m log n), iterative; idom[v] (root: itself, -1 unreachable)
    n = len(g); rg = [[] for _ in range(n)]
    for u in range(n):
        for v in g[u]: rg[v].append(u)
    num = [-1] * n; order = [root]; par = [0]; num[root] = 0; it = [0] * n; st = [root]
    while st:  # true DFS preorder; vertices renumbered 0..N-1
        u = st[-1]
        if it[u] < len(g[u]):
            v = g[u][it[u]]; it[u] += 1
            if num[v] < 0: num[v] = len(order); order.append(v); par.append(num[u]); st.append(v)
        else: st.pop()
    N = len(order); sdom = list(range(N)); label = list(range(N)); anc = [-1] * N
    idom = [0] * N; bucket = [[] for _ in range(N)]
    def ev(v):
        if anc[v] < 0: return v
        path = []
        while anc[anc[v]] >= 0: path.append(v); v = anc[v]
        for x in reversed(path):
            a = anc[x]
            if sdom[label[a]] < sdom[label[x]]: label[x] = label[a]
            anc[x] = anc[a]
        return label[path[0]] if path else label[v]
    for w in range(N - 1, 0, -1):
        for pv in rg[order[w]]:
            if num[pv] >= 0:
                u = ev(num[pv])
                if sdom[u] < sdom[w]: sdom[w] = sdom[u]
        bucket[sdom[w]].append(w); p = par[w]; anc[w] = p
        for v in bucket[p]:
            u = ev(v); idom[v] = u if sdom[u] < sdom[v] else p
        bucket[p] = []
    for w in range(1, N):
        if idom[w] != sdom[w]: idom[w] = idom[idom[w]]
    res = [-1] * n
    for w in range(N): res[order[w]] = order[idom[w]]
    return res

# ==== 6.6 Connectivity ====
def lowlink(n, edges):  # O(n+m) iterative; undirected, multi-edges ok (edge ids), self-loops ignored.
    g = [[] for _ in range(n)]  # returns (is_art[v], is_bridge[e], bccs = vertex lists of blocks)
    for i, (u, v) in enumerate(edges):
        if u != v: g[u].append((v, i)); g[v].append((u, i))
    tin = [-1] * n; low = [0] * n; it = [0] * n; pe = [-1] * n
    is_art = [False] * n; is_bridge = [False] * len(edges); bccs = []; t = 0
    for r in range(n):
        if tin[r] >= 0: continue
        tin[r] = low[r] = t; t += 1
        if not g[r]: bccs.append([r]); continue  # isolated vertex = own block
        st = [r]; vs = [r]; ch = 0
        while st:
            u = st[-1]
            if it[u] < len(g[u]):
                v, e = g[u][it[u]]; it[u] += 1
                if e == pe[u]: continue
                if tin[v] < 0:
                    tin[v] = low[v] = t; t += 1; pe[v] = e; st.append(v); vs.append(v)
                elif tin[v] < low[u]: low[u] = tin[v]
            else:
                st.pop()
                if not st: break
                p = st[-1]
                if low[u] < low[p]: low[p] = low[u]
                if low[u] > tin[p]: is_bridge[pe[u]] = True
                if low[u] >= tin[p]:
                    if p == r: ch += 1
                    else: is_art[p] = True
                    comp = [p]
                    while True:
                        x = vs.pop(); comp.append(x)
                        if x == u: break
                    bccs.append(comp)
        if ch > 1: is_art[r] = True
    return is_art, is_bridge, bccs

def two_edge_cc(n, edges):  # O(n+m) 2-edge-connected component ids (remove bridges, flood). uses lowlink
    _, br, _ = lowlink(n, edges); g = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        if not br[i]: g[u].append(v); g[v].append(u)
    return components(g)  # (comp, k); bridge tree: comps as nodes, bridges as edges

def block_cut_tree(n, bccs):  # nodes 0..n-1 vertices, n+i block i; edge vertex--block (forest)
    t = [[] for _ in range(n + len(bccs))]  # variant: keep only articulation vertices as nodes
    for i, b in enumerate(bccs):
        for v in b: t[v].append(n + i); t[n + i].append(v)
    return t
# Edge connectivity lambda(s,t) = max flow unit caps; vertex connectivity: split v -> v_in->v_out
# cap 1 (see 05). Global min cut: Stoer-Wagner (05). Online dynamic connectivity: Holm et al.
# (HDT, levels of Euler-tour trees, O(log^2) amortized) or link-cut trees for forests - rarely needed.

def offline_connectivity(n, queries):  # O(q log q log n) segment tree over time + rollback DSU (inline; also in 02)
    # queries: (0,u,v) add edge, (1,u,v) remove existing edge, (2,u,v) ask connected. multi-edges ok
    Q = len(queries); open_ = {}; iv = []
    for t, (tp, u, v) in enumerate(queries):
        if u > v: u, v = v, u
        if tp == 0: open_.setdefault((u, v), []).append(t)
        elif tp == 1: iv.append((open_[(u, v)].pop(), t, u, v))
    for (u, v), ts in open_.items():
        for s in ts: iv.append((s, Q, u, v))
    size = 1
    while size < max(Q, 1): size *= 2
    seg = [[] for _ in range(2 * size)]
    for l, r, u, v in iv:  # edge alive on time [l, r)
        l += size; r += size
        while l < r:
            if l & 1: seg[l].append((u, v)); l += 1
            if r & 1: r -= 1; seg[r].append((u, v))
            l >>= 1; r >>= 1
    p = list(range(n)); sz = [1] * n; hist = []; mark = [0] * (2 * size); ans = []
    def find(x):  # no path compression (rollback)
        while p[x] != x: x = p[x]
        return x
    st = [1]
    while st:
        x = st.pop()
        if x < 0:  # exit node: rollback
            x = ~x
            while len(hist) > mark[x]:
                a, b = hist.pop(); p[a] = a; sz[b] -= sz[a]
            continue
        mark[x] = len(hist)
        for u, v in seg[x]:
            a, b = find(u), find(v)
            if a != b:
                if sz[a] > sz[b]: a, b = b, a
                p[a] = b; sz[b] += sz[a]; hist.append((a, b))
        st.append(~x)
        if x >= size:
            t = x - size
            if t < Q and queries[t][0] == 2: ans.append(find(queries[t][1]) == find(queries[t][2]))
        else:
            st.append(2 * x + 1); st.append(2 * x)
    return ans  # component count: keep comps -= 1 per union, += 1 per rollback

# ==== 6.7 Eulerian and Hamiltonian paths ====
def euler_directed(n, edges):  # Hierholzer O(n+m); edges (u,v). returns (verts, edge ids) or None
    g = [[] for _ in range(n)]; deg = [0] * n
    for i, (u, v) in enumerate(edges): g[u].append(i); deg[u] += 1; deg[v] -= 1
    m = len(edges)
    if m == 0: return [0], []
    st_ = [v for v in range(n) if deg[v] == 1]
    if any(abs(d) > 1 for d in deg) or len(st_) > 1: return None
    start = st_[0] if st_ else edges[0][0]  # circuit wanted: require all deg==0
    it = [0] * n; st = [(start, -1)]; vs = []; es = []
    while st:
        u, e = st[-1]
        if it[u] < len(g[u]):
            i = g[u][it[u]]; it[u] += 1; st.append((edges[i][1], i))
        else:
            st.pop(); vs.append(u); es.append(e)
    if len(vs) != m + 1: return None  # disconnected
    return vs[::-1], es[::-1][1:]

def euler_undirected(n, edges):  # Hierholzer O(n+m); (verts, edge ids) or None; self-loops ok
    g = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges): g[u].append((v, i)); g[v].append((u, i))
    m = len(edges)
    if m == 0: return [0], []
    odd = [v for v in range(n) if len(g[v]) & 1]
    if len(odd) > 2: return None
    start = odd[0] if odd else edges[0][0]
    used = [False] * m; it = [0] * n; st = [(start, -1)]; vs = []; es = []
    while st:
        u, e = st[-1]
        while it[u] < len(g[u]) and used[g[u][it[u]][1]]: it[u] += 1
        if it[u] < len(g[u]):
            v, i = g[u][it[u]]; it[u] += 1; used[i] = True; st.append((v, i))
        else:
            st.pop(); vs.append(u); es.append(e)
    if len(vs) != m + 1: return None
    return vs[::-1], es[::-1][1:]

def de_bruijn(k, n):  # O(k^n) cyclic sequence with every length-n word over [0,k) once (FKM, Lyndon words)
    seq = []; w = [-1]  # linear version: seq + seq[:n-1]
    while w:
        w[-1] += 1; m = len(w)
        if n % m == 0: seq.extend(w)
        while len(w) < n: w.append(w[-m])
        while w and w[-1] == k - 1: w.pop()
    return seq  # lexicographically smallest. Alternative: Euler circuit on (n-1)-words graph.

def ham_path(W):  # O(2^n n^2) min cost Hamiltonian path (any endpoints); W[u][v] or INF. (cost, path)
    n = len(W); N = 1 << n; dp = [INF] * (N * n); par = [-1] * (N * n)
    for v in range(n): dp[(1 << v) * n + v] = 0
    for S in range(1, N):
        base = S * n
        for v in range(n):
            cur = dp[base + v]
            if cur == INF: continue
            row = W[v]
            for u in range(n):
                if S >> u & 1: continue
                j = (S | 1 << u) * n + u
                if cur + row[u] < dp[j]: dp[j] = cur + row[u]; par[j] = v
    full = N - 1; v = min(range(n), key=lambda v: dp[full * n + v]); best = dp[full * n + v]
    if best == INF: return INF, None
    path = []; S = full
    while v != -1: path.append(v); pv = par[S * n + v]; S ^= 1 << v; v = pv
    return best, path[::-1]  # cycle: fix start 0, close with W[v][0]. count paths: same DP with +=
# Chinese postman (undirected): sum(w) + min-weight perfect matching of odd-degree vertices
#   under shortest-path distances (bitmask DP over odd set, <=~20). Directed: min cost flow
#   from vertices with in>out to out>in. Must be connected.

# ==== 6.8 / 6.9 Network flow, matching: see 05_flow_matching.py ====

# ==== 6.10 Special graph topics ====
def chromatic_number(adj):  # O(2^n n); adj[v] bitmask. inclusion-exclusion over independent sets
    n = len(adj); N = 1 << n; MOD = (1 << 61) - 1  # mod: tiny false-negative risk, k-coloring count
    if n == 0: return 0
    f = [0] * N; f[0] = 1; sgn = [0] * N  # f[S] = #independent subsets of S (incl. empty)
    for S in range(1, N):
        v = (S & -S).bit_length() - 1
        f[S] = f[S ^ (1 << v)] + f[S & ~adj[v] & ~(1 << v)]; sgn[S] = sgn[S >> 1] ^ (S & 1)
    odd = n & 1; pw = [1] * N
    for k in range(1, n + 1):
        tot = 0
        for S in range(N):
            pw[S] = pw[S] * f[S] % MOD
            tot += -pw[S] if sgn[S] ^ odd else pw[S]  # sign (-1)^(n-|S|)
        if tot % MOD: return k
    return n
# Greedy/DSatur give upper bounds; 2-colorable = bipartite; planar <= 4; interval/chordal: greedy on PEO.

def max_clique(adj):  # Bron-Kerbosch with pivot, iterative, bitmasks; n <= ~50-60. returns clique mask
    n = len(adj); best = 0; bc = 0; st = [(0, (1 << n) - 1, 0, 0)]  # max independent set: complement graph
    while st:
        R, P, X, rc = st.pop()
        if not P:
            if rc > bc: best, bc = R, rc
            continue
        if rc + bin(P).count('1') <= bc: continue  # bound
        PX = P | X; piv = 0; pc = -1
        while PX:
            b = PX & -PX; u = b.bit_length() - 1; PX ^= b
            c = bin(P & adj[u]).count('1')
            if c > pc: pc = c; piv = u
        cand = P & ~adj[piv]
        while cand:
            b = cand & -cand; v = b.bit_length() - 1; cand ^= b
            st.append((R | b, P & adj[v], X & adj[v], rc + 1)); P ^= b; X |= b
    return best
# Meet in the middle for MIS n<=40: clique-check masks of the left half, for each right-half
# clique C take dp[common neighbours in left] (dp = best clique size inside a left subset).
# Planar: m <= 3n-6 (n>=3), bipartite planar m <= 2n-4; Euler V - E + F = 1 + C (C components).
#   K5 / K3,3 minors forbidden; planar graphs are 5-degenerate (some vertex deg <= 5).

def tournament_path(n, beats):  # O(n log n) beats() calls; beats(u,v) True iff edge u->v
    path = []  # Hamiltonian path always exists (binary insertion)
    for v in range(n):
        if not path or beats(v, path[0]): path.insert(0, v)
        elif beats(path[-1], v): path.append(v)
        else:
            lo, hi = 0, len(path) - 1  # invariant: path[lo]->v, v->path[hi]
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if beats(path[mid], v): lo = mid
                else: hi = mid
            path.insert(hi, v)
    return path
# Tournaments: SCC condensation is a total order (a path); strongly connected => Hamiltonian
#   cycle (Camion). Max out-degree vertex is a king (reaches all in <= 2). Landau: sorted scores
#   s valid iff prefix sums >= C(k,2) and total = C(n,2).
# Cactus (each edge in <= 1 cycle): blocks are cycles or bridges -> use lowlink bccs / block-cut
#   tree, then tree DP over it; vertex cactus = each vertex in <= 1 cycle.
# Chordal (every cycle >= 4 has a chord): perfect elimination ordering via max cardinality
#   search (reverse MCS order); then max clique / chromatic number / MIS are greedy on PEO.

def count_spanning_trees(n, edges, mod):  # Kirchhoff O(n^3): det of Laplacian minor mod prime p
    L = [[0] * n for _ in range(n)]  # multigraph ok, self-loops ignored. det also in 09
    for u, v in edges:
        if u != v: L[u][u] += 1; L[v][v] += 1; L[u][v] -= 1; L[v][u] -= 1
    a = [[x % mod for x in row[1:]] for row in L[1:]]; m = n - 1; det = 1
    for c in range(m):
        piv = next((r for r in range(c, m) if a[r][c]), -1)
        if piv < 0: return 0
        if piv != c: a[c], a[piv] = a[piv], a[c]; det = -det
        det = det * a[c][c] % mod; inv = pow(a[c][c], mod - 2, mod)
        for r in range(c + 1, m):
            if a[r][c]:
                f = a[r][c] * inv % mod; ar = a[r]; ac = a[c]
                for j in range(c, m): ar[j] = (ar[j] - f * ac[j]) % mod
    return det % mod
# Directed (arborescences rooted at r): L = diag(indeg) - A (A[u][v] = #edges u->v), delete row/col r.
# Weighted: use weights instead of 1 (sum over trees of product of weights).

def prufer_encode(n, edges):  # O(n) tree on 0..n-1 (n>=2) -> code of length n-2
    g = [[] for _ in range(n)]
    for u, v in edges: g[u].append(v); g[v].append(u)
    par = [-1] * n; st = [n - 1]; seen = [False] * n; seen[n - 1] = True
    while st:
        u = st.pop()
        for v in g[u]:
            if not seen[v]: seen[v] = True; par[v] = u; st.append(v)
    deg = [len(x) for x in g]; ptr = 0
    while deg[ptr] != 1: ptr += 1
    leaf = ptr; code = []
    for _ in range(n - 2):
        nx = par[leaf]; code.append(nx); deg[nx] -= 1
        if deg[nx] == 1 and nx < ptr: leaf = nx
        else:
            ptr += 1
            while deg[ptr] != 1: ptr += 1
            leaf = ptr
    return code

def prufer_decode(code):  # O(n) code -> edge list; #labeled trees = n^(n-2); deg(v) = 1 + count(v)
    n = len(code) + 2; deg = [1] * n
    for x in code: deg[x] += 1
    ptr = 0
    while deg[ptr] != 1: ptr += 1
    leaf = ptr; edges = []
    for v in code:
        edges.append((leaf, v)); deg[v] -= 1
        if deg[v] == 1 and v < ptr: leaf = v
        else:
            ptr += 1
            while deg[ptr] != 1: ptr += 1
            leaf = ptr
    edges.append((leaf, n - 1))
    return edges
# Forests: #forests on n labeled vertices with k trees rooted at given k vertices = k n^(n-k-1).
# Tree isomorphism (AHU / hashing): see 06_trees.py. General graph isomorphism heuristic:
#   Weisfeiler-Lehman colour refinement (iteratively hash (colour, sorted neighbour colours)),
#   compare colour histograms; equal histograms != isomorphic (regular graphs fool it).


# ==================================== tests ====================================
if __name__ == "__main__":
    R = random.Random(12345)

    def rand_edges(n, m, loops=False):
        es = []
        while len(es) < m:
            u, v = R.randrange(n), R.randrange(n)
            if u != v or loops: es.append((u, v))
        return es

    def adj_und(n, es):
        g = [[] for _ in range(n)]
        for u, v in es: g[u].append(v); g[v].append(u)
        return g

    def adj_dir(n, es):
        g = [[] for _ in range(n)]
        for u, v in es: g[u].append(v)
        return g

    def reach(g, s, banned=-1):
        if s == banned: return set()
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for v in g[u]:
                if v not in seen and v != banned: seen.add(v); st.append(v)
        return seen

    def n_comps(n, es, skip_v=-1, skip_e=-1):
        p = list(range(n))
        def f(x):
            while p[x] != x: x = p[x]
            return x
        c = n - (1 if skip_v >= 0 else 0)
        for i, (u, v) in enumerate(es):
            if i == skip_e or u == skip_v or v == skip_v: continue
            a, b = f(u), f(v)
            if a != b: p[a] = b; c -= 1
        return c

    def is_dag(n, g):
        return not any(u in reach(g, v) for u in range(n) for v in g[u])

    # csr, bipartite, bfs, components
    for _ in range(200):
        n = R.randint(1, 7); es = rand_edges(n, R.randint(0, 10)) if n > 1 else []
        st, to = csr(n, es)
        assert all(sorted(to[st[u]:st[u + 1]]) == sorted(v for a, v in es if a == u) for u in range(n))
        g = adj_und(n, es)
        c = bipartite(g)
        ok = any(all(m >> u & 1 != m >> v & 1 for u, v in es) for m in range(1 << n))
        assert (c is not None) == ok
        if c: assert all(c[u] != c[v] for u, v in es)
        comp, k = components(g)
        assert k == n_comps(n, es)
        assert all(comp[u] == comp[v] for u, v in es)
        srcs = R.sample(range(n), R.randint(1, n))
        d, par = bfs(g, srcs)
        bd = [0 if v in srcs else -1 for v in range(n)]
        for _ in range(n):
            for u, v in es + [(b, a) for a, b in es]:
                if bd[u] >= 0 and (bd[v] < 0 or bd[u] + 1 < bd[v]): bd[v] = bd[u] + 1
        assert d == bd
    # dfs_times / edge classification, directed cycles, toposort
    for _ in range(300):
        n = R.randint(1, 7); es = rand_edges(n, R.randint(0, 10), loops=True)
        g = adj_dir(n, es)
        tin, tout, par, kind = dfs_times(g)
        def anc(a, b):  # a ancestor-or-self of b in DFS forest
            while b != -1:
                if a == b: return True
                b = par[b]
            return False
        for u in range(n):
            assert anc(u, u) == (tin[u] <= tin[u] and tout[u] <= tout[u])
            for i, v in enumerate(g[u]):
                kd = kind[u][i]
                assert anc(u, v) == (tin[u] <= tin[v] and tout[v] <= tout[u])
                if kd == 0: assert par[v] == u
                elif kd == 1: assert anc(v, u)
                elif kd == 2: assert anc(u, v) and u != v
                else: assert not anc(u, v) and not anc(v, u) and tout[v] < tin[u]
        has_cycle = any(kind[u][i] == 1 for u in range(n) for i in range(len(g[u])))
        brute_cyc = not is_dag(n, g)
        assert has_cycle == brute_cyc
        cyc = find_cycle_directed(g)
        assert (cyc is not None) == brute_cyc
        if cyc:
            assert len(set(cyc)) == len(cyc)
            assert all(cyc[(i + 1) % len(cyc)] in g[cyc[i]] for i in range(len(cyc)))
        for lex in (False, True):
            o = toposort(g, lex)
            assert (o is None) == brute_cyc
            if o:
                pos = {v: i for i, v in enumerate(o)}
                assert sorted(o) == list(range(n)) and all(pos[u] < pos[v] for u, v in es)
        if not brute_cyc and n <= 6:
            allt = [p for p in permutations(range(n)) if all(p.index(u) < p.index(v) for u, v in es)]
            assert list(min(allt)) == toposort(g, True)
    # undirected cycle
    for _ in range(300):
        n = R.randint(1, 7); es = rand_edges(n, R.randint(0, 8), loops=True)
        cyc = find_cycle_undirected(n, es)
        exists = len(es) > n - n_comps(n, es)
        assert (cyc is not None) == exists
        if cyc:
            L = len(cyc); assert len(set(cyc)) == L
            from collections import Counter
            cnt = Counter(tuple(sorted(e)) for e in es)
            if L == 1: assert cnt[(cyc[0], cyc[0])] >= 1
            elif L == 2: assert cnt[tuple(sorted(cyc))] >= 2
            else: assert all(cnt[tuple(sorted((cyc[i], cyc[(i + 1) % L])))] for i in range(L))
    # flood fill
    for _ in range(100):
        Rr, C = R.randint(1, 6), R.randint(1, 6)
        grid = [''.join(R.choice('.#') for _ in range(C)) for _ in range(Rr)]
        lab, cnt = flood_fill(grid)
        cells = [(i, j) for i in range(Rr) for j in range(C) if grid[i][j] == '.']
        idx = {c: k for k, c in enumerate(cells)}
        es = [(idx[(i, j)], idx[(i + di, j + dj)]) for (i, j) in cells for di, dj in ((1, 0), (0, 1))
              if (i + di, j + dj) in idx]
        assert cnt == n_comps(len(cells), es)
        assert all(lab[i][j] == lab[i + di][j + dj] for (i, j) in cells for di, dj in ((1, 0), (0, 1)) if (i + di, j + dj) in idx)
        assert all(lab[i][j] == -1 for i in range(Rr) for j in range(C) if grid[i][j] == '#')
    # shortest paths
    for _ in range(200):
        n = R.randint(1, 7); es = [(u, v, R.randint(0, 9)) for u, v in rand_edges(n, R.randint(0, 14), True)]
        g = [[] for _ in range(n)]
        for u, v, w in es: g[u].append((v, w))
        M = [[0 if i == j else INF for j in range(n)] for i in range(n)]
        for u, v, w in es: M[u][v] = min(M[u][v], w)
        F = floyd([r[:] for r in M])
        s = R.randrange(n)
        d, p = dijkstra(g, s)
        assert d == F[s]
        for t in range(n):
            if d[t] < INF:
                P = get_path(p, t); assert P[0] == s and P[-1] == t
                assert sum(M[P[i]][P[i + 1]] for i in range(len(P) - 1)) == d[t]
        assert dijkstra_dense(M, s) == F[s]
        assert bellman_ford(n, es, s) == (F[s], None)
        g01 = [[(v, w & 1) for v, w in g[u]] for u in range(n)]
        assert bfs01(g01, s) == dijkstra(g01, s)[0]
        adj = [sum(1 << v for v in set(v for v, _ in g[u])) for u in range(n)]
        tc = transitive_closure(adj)
        for i in range(n):
            ri = set()
            for v, _ in g[i]: ri |= reach(adj_dir(n, [(a, b) for a, b, _ in es]), v)
            assert tc[i] == sum(1 << j for j in ri)
        # k shortest walks (weights >= 1)
        gk = [[(v, w + 1) for v, w in g[u]] for u in range(n)]
        t = R.randrange(n); k = R.randint(1, 8); B = 25
        res = k_shortest_walks(gk, s, t, k)
        ways = [[0] * n for _ in range(B + 1)]; ways[0][s] = 1
        for c in range(B + 1):
            for u in range(n):
                if ways[c][u]:
                    for v, w in gk[u]:
                        if c + w <= B: ways[c + w][v] += ways[c][u]
        brute = []
        for c in range(B + 1): brute += [c] * min(ways[c][t], k + 1)
        Rr = [x for x in res if x <= B]
        assert Rr == brute[:len(Rr)] and (len(Rr) == k or len(brute) == len(Rr))
    # negative weights: Bellman-Ford cycles, Johnson, difference constraints, DAG
    for _ in range(300):
        n = R.randint(1, 6); es = [(u, v, R.randint(-4, 9)) for u, v in rand_edges(n, R.randint(0, 10), True)]
        M = [[0 if i == j else INF for j in range(n)] for i in range(n)]
        for u, v, w in es: M[u][v] = min(M[u][v], w)
        F = floyd([r[:] for r in M]); s = R.randrange(n)
        negany = any(F[i][i] < 0 for i in range(n))
        negs = any(F[s][v] < INF and F[v][v] < 0 for v in range(n))
        d, cyc = bellman_ford(n, es, s)
        assert (cyc is not None) == negs
        if cyc:
            L = len(cyc); assert sum(M[cyc[i]][cyc[(i + 1) % L]] for i in range(L)) < 0
        else: assert d == F[s]
        J = johnson(n, es)
        assert (J is None) == negany
        if J: assert J == F
        x = diff_constraints(n, es)
        assert (x is None) == negany
        if x: assert all(x[v] - x[u] <= w for u, v, w in es)
        dag = [(u, v, w) for u, v, w in es if u < v]
        gd = [[] for _ in range(n)]
        for u, v, w in dag: gd[u].append((v, w))
        Md = [[0 if i == j else INF for j in range(n)] for i in range(n)]
        for u, v, w in dag: Md[u][v] = min(Md[u][v], w)
        assert dag_shortest(gd, s) == floyd(Md)[s]
    # MST family
    def spanning_trees(n, es):
        out = []
        for c in combinations(range(len(es)), n - 1):
            if n_comps(n, [es[i][:2] for i in c]) == 1: out.append(sum(es[i][2] for i in c))
        return sorted(out)
    for _ in range(200):
        n = R.randint(1, 6); es = [(u, v, R.randint(1, 6)) for u, v in rand_edges(n, R.randint(0, 9), True)]
        tr = spanning_trees(n, es)
        tot, used = kruskal(n, es)
        assert len(used) == n - n_comps(n, [e[:2] for e in es])
        tb, ub = boruvka(n, es)
        assert tb == tot and len(ub) == len(used) and n_comps(n, [es[i][:2] for i in ub]) == n_comps(n, [e[:2] for e in es])
        g = [[] for _ in range(n)]
        for u, v, w in es: g[u].append((v, w)); g[v].append((u, w))
        if tr:
            assert tot == tr[0] and prim(g) == tot
            gt = [x for x in tr if x > tr[0]]
            assert second_mst(n, es) == (gt[0] if gt else None)
        else:
            assert prim(g) is None and second_mst(n, es) is None
        par, val = kruskal_tree(n, es)
        for u in range(n):
            for v in range(n):
                A = set(); x = u
                while x != -1: A.add(x); x = par[x] if x < len(par) else -1
                y = v
                while y != -1 and y not in A: y = par[y]
                bw = None
                for W in sorted(set(e[2] for e in es)):
                    if v in reach(adj_und(n, [e[:2] for e in es if e[2] <= W]), u): bw = W; break
                if u == v: assert y == u
                elif bw is None: assert y == -1
                else: assert val[y] == bw
        # Steiner
        k = R.randint(0, min(n, 4)); terms = R.sample(range(n), k)
        if len(es) <= 9:
            best = INF
            for mask in range(1 << len(es)):
                sub = [es[i] for i in range(len(es)) if mask >> i & 1]
                if k == 0 or all(t in reach(adj_und(n, [e[:2] for e in sub]), terms[0]) for t in terms):
                    best = min(best, sum(e[2] for e in sub))
            assert steiner_tree(g, terms) == best
    # directed MST
    for _ in range(300):
        n = R.randint(1, 5); es = [(u, v, R.randint(-3, 9)) for u, v in rand_edges(n, R.randint(0, 10), True)]
        root = R.randrange(n); best = None
        choices = [[(a, w_) for (a, b, w_) in es if b == v and a != v] if v != root else [(-1, 0)] for v in range(n)]
        for ch in product(*choices):
            ok = True
            for v in range(n):
                x = v; steps = 0
                while x != root and steps <= n: x = ch[x][0]; steps += 1
                if x != root: ok = False; break
            if ok:
                s_ = sum(c[1] for c in ch)
                best = s_ if best is None else min(best, s_)
        assert directed_mst(n, es, root) == best
    # SCC, condensation, 2-SAT
    for _ in range(300):
        n = R.randint(1, 8); es = rand_edges(n, R.randint(0, 14), True); g = adj_dir(n, es)
        k, comp = scc(g)
        rs = [reach(g, u) for u in range(n)]
        for u in range(n):
            for v in range(n):
                assert (comp[u] == comp[v]) == (v in rs[u] and u in rs[v])
        assert all(comp[u] <= comp[v] for u, v in es) and len(set(comp)) == k
        dag = condensation(g, k, comp)
        assert all(a < b for a in range(k) for b in dag[a])
    for _ in range(300):
        n = R.randint(1, 6); ts = TwoSat(n); cl = []
        for _ in range(R.randint(0, 12)):
            a, b = R.randrange(2 * n), R.randrange(2 * n); cl.append((a, b)); ts.either(a, b)
        res = ts.solve()
        lit = lambda asg, l: asg[l >> 1] ^ (l & 1)
        sat = any(all(lit(a_, x) or lit(a_, y) for x, y in cl) for a_ in product((False, True), repeat=n))
        assert (res is not None) == sat
        if res: assert all(lit(res, x) or lit(res, y) for x, y in cl)
    # functional graph
    for _ in range(300):
        n = R.randint(1, 9); f = [R.randrange(n) for _ in range(n)]
        cycles, cid, dist, entry = functional_graph(f)
        for v in range(n):
            x = v; seq = []
            for _ in range(2 * n + 2): seq.append(x); x = f[x]
            oncyc = any(seq[i] == v for i in range(1, n + 1))
            assert (cid[v] >= 0) == oncyc
            d0 = next(i for i in range(len(seq)) if any(seq[j] == seq[i] for j in range(i + 1, i + n + 1)))
            assert dist[v] == d0 and entry[v] == seq[d0]
        assert sorted(v for c in cycles for v in c) == sorted(v for v in range(n) if cid[v] >= 0)
        for c in cycles: assert all(f[c[i]] == c[(i + 1) % len(c)] for i in range(len(c)))
    # dominator tree
    for _ in range(300):
        n = R.randint(1, 8); es = rand_edges(n, R.randint(0, 16), True); g = adj_dir(n, es)
        root = R.randrange(n); idom = dominator_tree(g, root); rr = reach(g, root)
        dom = {v: {x for x in rr if v not in reach(g, root, x)} | {v} for v in rr}
        for v in range(n):
            if v not in rr: assert idom[v] == -1
            elif v == root: assert idom[v] == root
            else:
                sd = dom[v] - {v}
                assert idom[v] == max(sd, key=lambda x: len(dom[x]))
    # lowlink: bridges, articulation, bcc, 2ecc, block-cut tree
    from collections import Counter
    for _ in range(300):
        n = R.randint(1, 8); es = rand_edges(n, R.randint(0, 11), True)
        art, br, bccs = lowlink(n, es)
        c0 = n_comps(n, es)
        for i in range(len(es)): assert br[i] == (n_comps(n, es, skip_e=i) > c0)
        for v in range(n): assert art[v] == (n_comps(n, es, skip_v=v) > c0 - (1 if all(v not in e for e in es) else 0))
        sets = [set(b) for b in bccs]
        assert all(any(v in s for s in sets) for v in range(n))
        for u, v in es:
            if u != v: assert sum(1 for s in sets if u in s and v in s) == 1
        for i in range(len(sets)):
            for j in range(i):
                common = sets[i] & sets[j]; assert len(common) <= 1
                for v in common: assert art[v]
        for s in sets:
            sub = [(u, v) for u, v in es if u in s and v in s and u != v]; idx = {v: i for i, v in enumerate(s)}
            se = [(idx[u], idx[v]) for u, v in sub]
            assert n_comps(len(s), se) == 1
            if len(s) >= 3:
                for x in range(len(s)): assert n_comps(len(s), se, skip_v=x) == 1
        comp, k = two_edge_cc(n, es)
        for u in range(n):
            for v in range(n):
                same = all(v in reach(adj_und(n, [e for j, e in enumerate(es) if j != i]), u) for i in range(len(es))) and v in reach(adj_und(n, es), u)
                assert (comp[u] == comp[v]) == same
        t = block_cut_tree(n, bccs)
        te = [(a, b) for a in range(len(t)) for b in t[a] if a < b]
        assert n_comps(len(t), te) == len(t) - len(te) == c0
    # offline dynamic connectivity
    for _ in range(100):
        n = R.randint(1, 7); qs = []; cur = []
        for _ in range(R.randint(0, 40)):
            r = R.random()
            if r < 0.4: u, v = R.randrange(n), R.randrange(n); qs.append((0, u, v)); cur.append((u, v))
            elif r < 0.6 and cur: e = cur.pop(R.randrange(len(cur))); qs.append((1, e[1], e[0]) if R.random() < .5 else (1, e[0], e[1]))
            else: qs.append((2, R.randrange(n), R.randrange(n)))
        exp = []; cur = []
        for tp, u, v in qs:
            if tp == 0: cur.append((u, v))
            elif tp == 1:
                cur.remove((u, v) if (u, v) in cur else (v, u))
            else: exp.append(v in reach(adj_und(n, cur), u))
        assert offline_connectivity(n, qs) == exp
    # Euler
    def euler_exists(n, es, directed):
        m = len(es)
        if m == 0: return True
        used = [False] * m
        def dfs(u, cnt):
            if cnt == m: return True
            for i, (a, b) in enumerate(es):
                if used[i]: continue
                if a == u: nxt = b
                elif not directed and b == u: nxt = a
                else: continue
                used[i] = True
                if dfs(nxt, cnt + 1): used[i] = False; return True
                used[i] = False
            return False
        return any(dfs(v, 0) for v in range(n))  # test-only recursion, depth <= m
    for _ in range(300):
        n = R.randint(1, 5); es = rand_edges(n, R.randint(0, 7), True)
        for directed in (True, False):
            r = (euler_directed if directed else euler_undirected)(n, es)
            assert (r is not None) == euler_exists(n, es, directed)
            if r and es:
                vs, ids = r; assert sorted(ids) == list(range(len(es))) and len(vs) == len(es) + 1
                for i, e in enumerate(ids):
                    a, b = es[e]
                    assert (vs[i], vs[i + 1]) == (a, b) or (not directed and (vs[i], vs[i + 1]) == (b, a))
    for k in range(1, 4):
        for nn in range(1, 5):
            s = de_bruijn(k, nn); L = k ** nn; assert len(s) == L
            ss = s + s[:nn - 1]
            assert len({tuple(ss[i:i + nn]) for i in range(L)}) == L
    # Hamiltonian path DP
    for _ in range(100):
        n = R.randint(1, 6)
        W = [[INF if i == j or R.random() < 0.3 else R.randint(1, 9) for j in range(n)] for i in range(n)]
        best = min(sum(W[p[i]][p[i + 1]] for i in range(n - 1)) for p in permutations(range(n)))
        c, path = ham_path(W); assert c == best
        if path: assert sorted(path) == list(range(n)) and sum(W[path[i]][path[i + 1]] for i in range(n - 1)) == c
    # chromatic number, max clique
    for _ in range(150):
        n = R.randint(0, 7); es = rand_edges(n, R.randint(0, 14)) if n > 1 else []
        adj = [0] * n
        for u, v in es: adj[u] |= 1 << v; adj[v] |= 1 << u
        chi = next((k for k in range(n + 1) if any(all(c[u] != c[v] for u, v in es) for c in product(range(k), repeat=n))), n)
        assert chromatic_number(adj) == chi
        cl = max_clique(adj)
        assert all(adj[u] >> v & 1 for u in range(n) for v in range(n) if u != v and cl >> u & 1 and cl >> v & 1)
        bestc = max(bin(m).count('1') for m in range(1 << n) if all(adj[u] >> v & 1 for u in range(n) for v in range(n) if u != v and m >> u & 1 and m >> v & 1))
        assert bin(cl).count('1') == bestc
    # tournaments
    for _ in range(100):
        n = R.randint(1, 30); T = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(i):
                if R.random() < .5: T[i][j] = True
                else: T[j][i] = True
        p = tournament_path(n, lambda u, v: T[u][v])
        assert sorted(p) == list(range(n)) and all(T[p[i]][p[i + 1]] for i in range(n - 1))
    # Kirchhoff
    for _ in range(150):
        n = R.randint(1, 6); es = rand_edges(n, R.randint(0, 9), True)
        cnt = sum(1 for c in combinations(range(len(es)), n - 1) if n_comps(n, [es[i] for i in c]) == 1)
        assert count_spanning_trees(n, es, 998244353) == cnt % 998244353
    # Prufer
    for n in range(2, 7):
        codes = set()
        for code in product(range(n), repeat=n - 2):
            es = prufer_decode(list(code))
            assert len(es) == n - 1 and n_comps(n, es) == 1
            assert prufer_encode(n, es) == list(code); codes.add(tuple(sorted(tuple(sorted(e)) for e in es)))
        assert len(codes) == n ** (n - 2)
    # speed smoke test
    nB = 200000; esB = rand_edges(nB, 400000)
    gB = adj_dir(nB, esB); scc(gB); dfs_times(gB)
    lowlink(nB, esB)
    gw = [[] for _ in range(nB)]
    for u, v in esB: gw[u].append((v, R.randint(1, 100)))
    dijkstra(gw, 0)
    print("ALL OK")
