# 06_trees.py -- trees: rooting, diameter/center/centroid, LCA (3 ways), Euler tour, sack,
# HLD, centroid decomposition, virtual tree, isomorphism, long-path decomposition.
# Python 3 stdlib only, PyPy-friendly, everything iterative. Trees: adj = list of lists, nodes 0..n-1.

# ==== 7.1 Tree basics: rooting (parent, depth, BFS order, subtree sizes) ====
# O(n). order is BFS order: parents before children; reversed(order) = children before parents.
def root_tree(adj, root=0):
    n = len(adj)
    par = [-1] * n
    dep = [0] * n
    order = [root]
    for v in order:                     # list grows while iterating: that's the BFS queue
        for u in adj[v]:
            if u != par[v]:
                par[u] = v
                dep[u] = dep[v] + 1
                order.append(u)
    sz = [1] * n
    for v in reversed(order):
        if par[v] >= 0:
            sz[par[v]] += sz[v]
    return par, dep, order, sz
# height of node v (longest downward path) = same reversed loop with h[p] = max(h[p], h[v]+1).


# ==== 7.1 Tree traversals (note) ====
# Preorder = order nodes are popped in a stack DFS (euler_tin_tout). Postorder = reverse of a
# preorder that visits children in reversed order. For bottom-up DP any "children before parent"
# order works: reversed(BFS order). Level order = BFS order. Inorder (binary trees only): walk
# left pushing nodes, pop + visit, go right.


# ==== 7.1 Diameter (two BFS) and center(s) ====
def bfs_dist(adj, s):                   # O(n) unweighted distances from s
    d = [-1] * len(adj)
    d[s] = 0
    q = [s]
    for v in q:
        for u in adj[v]:
            if d[u] < 0:
                d[u] = d[v] + 1
                q.append(u)
    return d


def diameter(adj):                      # O(n) -> (length in edges, endpoint a, endpoint b)
    d = bfs_dist(adj, 0)
    a = d.index(max(d))
    d = bfs_dist(adj, a)
    b = d.index(max(d))
    return d[b], a, b


def tree_centers(adj):                  # O(n) -> [c] or [c1, c2]: middle of a diameter path
    D, a, b = diameter(adj)
    par = root_tree(adj, a)[0]
    path = [b]
    while path[-1] != a:
        path.append(par[path[-1]])
    return [path[D // 2]] if D % 2 == 0 else [path[D // 2], path[D // 2 + 1]]
# weighted diameter: same two-pass idea with weighted distances (or DP: best two child depths).


# ==== 7.1 Centroid(s) ====
def tree_centroids(adj):                # O(n) -> 1 or 2 nodes whose removal leaves parts <= n/2
    n = len(adj)
    par, _, _, sz = root_tree(adj)
    res = []
    for v in range(n):
        mx = n - sz[v]
        for u in adj[v]:
            if u != par[v] and sz[u] > mx:
                mx = sz[u]
        if 2 * mx <= n:
            res.append(v)
    return res


# ==== 7.2 LCA via binary lifting (+ k-th ancestor, distance) ====
# Build O(n log n), query O(log n).  L = LCA(adj); L.lca(u,v); L.kth(v,k); L.dist(u,v)
class LCA:
    def __init__(self, adj, root=0):
        par, self.dep, _, _ = root_tree(adj, root)
        n = len(adj)
        self.up = up = [[root if p < 0 else p for p in par]]
        for _ in range(max(1, (n - 1).bit_length()) - 1):
            prv = up[-1]
            up.append([prv[prv[v]] for v in range(n)])

    def kth(self, v, k):                # k-th ancestor of v, -1 if k > depth
        if k > self.dep[v]:
            return -1
        i = 0
        while k:
            if k & 1:
                v = self.up[i][v]
            k >>= 1
            i += 1
        return v

    def lca(self, u, v):
        dep, up = self.dep, self.up
        if dep[u] < dep[v]:
            u, v = v, u
        u = self.kth(u, dep[u] - dep[v])
        if u == v:
            return u
        for k in range(len(up) - 1, -1, -1):
            a = up[k]
            if a[u] != a[v]:
                u, v = a[u], a[v]
        return up[0][u]

    def dist(self, u, v):
        return self.dep[u] + self.dep[v] - 2 * self.dep[self.lca(u, v)]
# k-th node on path u->v: w=lca; if k<=dep[u]-dep[w]: kth(u,k) else kth(v, dist-k).
# weighted dist: keep wdep[] (root distance) and use wdep[u]+wdep[v]-2*wdep[lca].


# ==== 7.2 LCA via Euler tour + sparse table: O(n log n) build, O(1) query ====
class LCAEuler:
    def __init__(self, adj, root=0):
        n = len(adj)
        dep = [0] * n
        par = [-1] * n
        first = [0] * n
        it = [0] * n
        euler = [root]
        stack = [root]
        while stack:
            v = stack[-1]
            if it[v] < len(adj[v]):
                u = adj[v][it[v]]
                it[v] += 1
                if u != par[v]:
                    par[u] = v
                    dep[u] = dep[v] + 1
                    first[u] = len(euler)
                    euler.append(u)
                    stack.append(u)
            else:
                stack.pop()
                if stack:
                    euler.append(stack[-1])
        self.dep, self.first = dep, first
        sp = [euler]                    # sp[k][i] = min-depth node in euler[i : i+2^k]
        j = 1
        while 2 * j <= len(euler):
            p = sp[-1]
            sp.append([a if dep[a] < dep[b] else b for a, b in zip(p, p[j:])])
            j *= 2
        self.sp = sp

    def lca(self, u, v):
        l, r = self.first[u], self.first[v]
        if l > r:
            l, r = r, l
        k = (r - l + 1).bit_length() - 1
        a, b = self.sp[k][l], self.sp[k][r - (1 << k) + 1]
        return a if self.dep[a] < self.dep[b] else b
# Memory trick: sparse table over preorder (n entries): for tin[u]<tin[v], lca = parent of the
# min-depth node in preorder positions (tin[u], tin[v]].


# ==== 7.2 Tarjan offline LCA (iterative DFS + DSU), O((n+q) alpha) ====
def tarjan_lca(adj, queries, root=0):   # queries: list of (u, v) -> list of lca
    n = len(adj)
    qs = [[] for _ in range(n)]
    for i, (u, v) in enumerate(queries):
        qs[u].append((v, i))
        qs[v].append((u, i))
    dsu = list(range(n))
    anc = list(range(n))
    done = [False] * n
    par = [-1] * n
    it = [0] * n
    ans = [-1] * len(queries)

    def find(x):
        while dsu[x] != x:
            dsu[x] = dsu[dsu[x]]
            x = dsu[x]
        return x
    stack = [root]
    while stack:
        v = stack[-1]
        if it[v] < len(adj[v]):
            u = adj[v][it[v]]
            it[v] += 1
            if u != par[v]:
                par[u] = v
                stack.append(u)
            continue
        stack.pop()                     # v finished: its children's sets are already merged into v
        done[v] = True
        for w, i in qs[v]:
            if done[w]:
                ans[i] = anc[find(w)]
        p = par[v]
        if p >= 0:
            rp = find(p)
            dsu[find(v)] = rp
            anc[rp] = p
    return ans


# ==== 7.2 Farach-Colton-Bender (conceptual) ====
# O(n) build / O(1) LCA: Euler tour depths differ by +-1; split into blocks of (log n)/2,
# sparse table over block minima, and precompute in-block answers per block "shape" (2^(b-1)
# shapes). Never needed in contests: LCAEuler / binary lifting are fast enough.


# ==== 7.3 Euler tour flattening (tin/tout for subtree queries) ====
# O(n). Subtree of v = positions [tin[v], tout[v]) in preorder; ord[tin[v]] = v.
# u in subtree(v)  <=>  tin[v] <= tin[u] < tout[v].
def euler_tin_tout(adj, root=0):
    n = len(adj)
    tin = [0] * n
    par = [-1] * n
    ord_ = []
    stack = [root]
    while stack:
        v = stack.pop()
        tin[v] = len(ord_)
        ord_.append(v)
        for u in adj[v]:
            if u != par[v]:
                par[u] = v
                stack.append(u)
    sz = [1] * n
    for v in reversed(ord_):
        if par[v] >= 0:
            sz[par[v]] += sz[v]
    tout = [tin[v] + sz[v] for v in range(n)]
    return tin, tout, ord_
# Example: subtree-sum with point updates -> Fenwick over positions: add(tin[v], x),
# query sum(tin[v], tout[v]). Path-to-root sums with point updates: add +x at tin[v], -x at tout[v],
# then prefix(tin[u]) = sum on root..u. (Fenwick: see 02_data_structures.py)


# ==== 7.3 Small-to-large / DSU on tree (sack): distinct colors in every subtree ====
# O(n log n). Light children computed and wiped; heavy child's counts kept and extended.
def subtree_distinct_colors(adj, col, root=0):
    n = len(adj)
    par, _, _, sz = root_tree(adj, root)
    tin, tout, ord_ = euler_tin_tout(adj, root)
    heavy = [-1] * n
    for v in range(n):
        for u in adj[v]:
            if u != par[v] and (heavy[v] < 0 or sz[u] > sz[heavy[v]]):
                heavy[v] = u
    cnt = [0] * (max(col) + 1)
    ans = [0] * n
    distinct = 0
    stack = [(root, 1, 0)]
    while stack:
        v, keep, phase = stack.pop()
        h = heavy[v]
        if phase == 0:
            stack.append((v, keep, 1))
            if h >= 0:
                stack.append((h, 1, 0))     # heavy runs after all light ones (LIFO)
            for c in adj[v]:
                if c != par[v] and c != h:
                    stack.append((c, 0, 0))
            continue
        cv = col[v]
        if cnt[cv] == 0:
            distinct += 1
        cnt[cv] += 1
        for c in adj[v]:
            if c != par[v] and c != h:
                for u in ord_[tin[c]:tout[c]]:
                    cu = col[u]
                    if cnt[cu] == 0:
                        distinct += 1
                    cnt[cu] += 1
        ans[v] = distinct               # <- answer subtree queries of v here
        if not keep:
            for u in ord_[tin[v]:tout[v]]:
                cnt[col[u]] -= 1
            distinct = 0
    return ans
# Alternative small-to-large with sets: for v in reversed(BFS order): S[v] = largest child set,
# merge the others into it (S[v] |= S[c]); O(n log n) set ops; easiest when values are hashable.


# ==== 7.3 Heavy-light decomposition (path query / point update) ====
# Build O(n). Each path = O(log n) contiguous segments of pos[]. Subtree v = [pos[v], pos[v]+sz[v]).
class HLD:
    def __init__(self, adj, root=0):
        n = len(adj)
        self.par, self.dep, _, self.sz = par, dep, _, sz = root_tree(adj, root)
        heavy = [-1] * n
        for v in range(n):
            for u in adj[v]:
                if u != par[v] and (heavy[v] < 0 or sz[u] > sz[heavy[v]]):
                    heavy[v] = u
        self.head = head = [root] * n
        self.pos = pos = [0] * n
        cur = 0
        stack = [root]
        while stack:
            v = stack.pop()
            pos[v] = cur
            cur += 1
            for u in adj[v]:
                if u != par[v] and u != heavy[v]:
                    head[u] = u
                    stack.append(u)
            if heavy[v] >= 0:                # pushed last -> popped next -> chain contiguous
                head[heavy[v]] = head[v]
                stack.append(heavy[v])

    def segments(self, u, v, edge=False):   # inclusive (l, r) pos-ranges covering path u-v
        head, pos, par, dep = self.head, self.pos, self.par, self.dep
        res = []
        while head[u] != head[v]:
            if dep[head[u]] < dep[head[v]]:
                u, v = v, u
            res.append((pos[head[u]], pos[u]))
            u = par[head[u]]
        if dep[u] > dep[v]:
            u, v = v, u
        if pos[u] + edge <= pos[v]:         # edge=True: values stored on child endpoint, skip lca
            res.append((pos[u] + edge, pos[v]))
        return res

    def lca(self, u, v):
        head, par, dep = self.head, self.par, self.dep
        while head[u] != head[v]:
            if dep[head[u]] < dep[head[v]]:
                u, v = v, u
            u = par[head[u]]
        return u if dep[u] < dep[v] else v


class _BIT:                             # tiny Fenwick for the HLD example (full one: 02)
    def __init__(self, n):
        self.n, self.t = n, [0] * (n + 1)

    def add(self, i, x):
        i += 1
        while i <= self.n:
            self.t[i] += x
            i += i & -i

    def pre(self, i):                   # sum of [0, i)
        s = 0
        while i:
            s += self.t[i]
            i &= i - 1
        return s


class HLDPathSum:                       # point set / path sum, O(log^2 n) per query
    def __init__(self, adj, vals, root=0):
        self.h = HLD(adj, root)
        self.val = list(vals)
        self.bit = _BIT(len(adj))
        for v, x in enumerate(vals):
            self.bit.add(self.h.pos[v], x)

    def set(self, v, x):
        self.bit.add(self.h.pos[v], x - self.val[v])
        self.val[v] = x

    def query(self, u, v):
        b = self.bit
        return sum(b.pre(r + 1) - b.pre(l) for l, r in self.h.segments(u, v))
# For max/min or range updates on paths use a segment tree / lazy segtree (02) over pos instead.


# ==== 7.3 Centroid decomposition ====
# O(n log n). Returns cpar (parent in centroid tree, -1 for top) and order (ancestors first).
# Depth of centroid tree <= log2 n. Any path u-v passes through the centroid-tree LCA of u, v.
def centroid_decomp(adj):
    n = len(adj)
    removed = [False] * n
    cpar = [-1] * n
    order = []
    par = [-1] * n
    sz = [0] * n
    stack = [(0, -1)]
    while stack:
        s, p = stack.pop()
        par[s] = -1
        comp = [s]
        for v in comp:
            for u in adj[v]:
                if u != par[v] and not removed[u]:
                    par[u] = v
                    comp.append(u)
        for v in comp:
            sz[v] = 1
        for v in reversed(comp):
            if par[v] >= 0:
                sz[par[v]] += sz[v]
        tot = len(comp)
        c = s
        while True:                     # walk toward the heavy side
            for u in adj[c]:
                if u != par[c] and not removed[u] and 2 * sz[u] > tot:
                    c = u
                    break
            else:
                break
        removed[c] = True
        cpar[c] = p
        order.append(c)
        for u in adj[c]:
            if not removed[u]:
                stack.append((u, c))
    return cpar, order


# Example: number of unordered pairs (u != v) with dist(u, v) <= K.  O(n log^2 n).
# Process centroids in `order`; when c is processed, its component = reachable without removed.
def count_pairs_dist_le(adj, K):
    n = len(adj)
    _, order = centroid_decomp(adj)
    removed = [False] * n

    def pairs(ds):                      # pairs i<j with ds[i]+ds[j] <= K
        ds.sort()
        i, j, r = 0, len(ds) - 1, 0
        while i < j:
            if ds[i] + ds[j] <= K:
                r += j - i
                i += 1
            else:
                j -= 1
        return r
    total = 0
    for c in order:
        alld = [0]
        for s in adj[c]:
            if removed[s]:
                continue
            ds = []
            q = [(s, c, 1)]
            for v, p, d in q:
                ds.append(d)
                for u in adj[v]:
                    if u != p and not removed[u]:
                        q.append((u, v, d + 1))
            total -= pairs(ds)          # remove pairs inside the same child subtree
            alld += ds
        total += pairs(alld)
        removed[c] = True
    return total
# Nearest marked node queries: best[c] = min dist to a marked node in c's component;
# mark v: for a in centroid ancestors of v (follow cpar): best[a] = min(best[a], dist(a, v));
# query v: min over ancestors a of best[a] + dist(a, v). Use LCA.dist, or store per-level dists.


# ==== 7.3 Rerooting DP ====
# see 03_dp.py: rerooting


# ==== 7.3 Virtual (auxiliary) tree ====
# O(k log k) given O(1)/O(log) lca. Returns (nodes sorted by tin, vpar dict), vpar[top] = -1.
# Contains the keys plus pairwise LCAs; vpar[x] = nearest proper ancestor of x in that set.
def virtual_tree(keys, tin, lca):
    vs = sorted(set(keys), key=tin.__getitem__)
    vs += [lca(vs[i], vs[i + 1]) for i in range(len(vs) - 1)]
    vs = sorted(set(vs), key=tin.__getitem__)
    vpar = {vs[0]: -1}
    for i in range(1, len(vs)):         # parent = lca of tin-adjacent nodes (set is LCA-closed)
        vpar[vs[i]] = lca(vs[i - 1], vs[i])
    return vs, vpar
# Edge weights = depth differences. Do the DP on vs in reverse tin order (children first).


# ==== 7.3 Tree hashing / isomorphism (AHU canonical labels) ====
# O(n log n). ids: dict shared by all trees being compared (tuple of child labels -> int).
# Rooted trees isomorphic <=> equal labels. Unrooted: root at center(s), take min label.
def rooted_label(adj, root, ids):
    par, _, order, _ = root_tree(adj, root)
    ch = [[] for _ in adj]
    lab = 0
    for v in reversed(order):
        key = tuple(sorted(ch[v]))
        lab = ids.get(key)
        if lab is None:
            lab = ids[key] = len(ids)
        if par[v] >= 0:
            ch[par[v]].append(lab)
    return lab                          # last processed = root


def unrooted_label(adj, ids):
    return min(rooted_label(adj, c, ids) for c in tree_centers(adj))
# Random hashing alternative (one pass, mergeable): h(v) = prod(r[depth(v)] + h(c)) mod P with
# random r per depth, or h(v) = 1 + sum f(h(c)) with f a random-ish mixing function (mod 2^64).
# Subtree-isomorphism counting: number of distinct labels among all subtrees = len(set(labels)).


# ==== 7.3 Link-cut trees ====
# see 02_data_structures.py: LinkCutTree


# ==== 7.3 Mo's algorithm on trees (note) ====
# Path queries: Euler tour where each node appears twice (in at st[v], out at en[v]), length 2n.
# For path u-v with st[u] <= st[v]: if lca == u -> range [st[u], st[v]], else [en[u], st[v]] plus
# the lca handled separately. Nodes appearing twice inside the range cancel (toggle on each add).
# Then run plain Mo over these ranges: see 01_techniques.py: Mo's algorithm.
# Subtree queries need no Mo: they're contiguous [tin, tout) ranges (euler_tin_tout).


# ==== 7.3 Long-path decomposition: depth-count DP in O(n) ====
# Example (CF 1009F "Dominant Indices"): for every v, smallest d maximizing #nodes at depth d in
# subtree(v). Each node's depth array lives in one shared buffer; the longest child's array starts
# at off[v]+1 so it is inherited for free; short children merged in O(their height).
def dominant_indices(adj, root=0):
    n = len(adj)
    par, _, order, _ = root_tree(adj, root)
    h = [1] * n
    lng = [-1] * n
    for v in reversed(order):
        p = par[v]
        if p >= 0 and h[v] + 1 > h[p]:
            h[p] = h[v] + 1
            lng[p] = v
    off = [0] * n
    nxt = h[root]
    for v in order:                     # parents first
        for c in adj[v]:
            if c == par[v]:
                continue
            if c == lng[v]:
                off[c] = off[v] + 1
            else:
                off[c] = nxt
                nxt += h[c]
    buf = [0] * n                       # buf[off[v] + d] = #nodes at depth d below v
    ans = [0] * n
    for v in reversed(order):
        o = off[v]
        buf[o] = 1
        b, bv = 0, 1
        if lng[v] >= 0:
            b2 = ans[lng[v]] + 1
            if buf[o + b2] > 1:
                b, bv = b2, buf[o + b2]
        for c in adj[v]:
            if c == par[v] or c == lng[v]:
                continue
            oc = off[c]
            for i in range(h[c]):
                x = buf[o + i + 1] + buf[oc + i]
                buf[o + i + 1] = x
                if x > bv or (x == bv and i + 1 < b):
                    b, bv = i + 1, x
        ans[v] = b
    return ans


# ==================================== tests ====================================
if __name__ == "__main__":
    import random
    from itertools import permutations
    random.seed(12345)

    def rand_tree(n, shape=None):
        perm = list(range(n))
        random.shuffle(perm)
        adj = [[] for _ in range(n)]
        edges = []
        for i in range(1, n):
            if shape == "path":
                p = i - 1
            elif shape == "star":
                p = 0
            else:
                p = random.randrange(max(0, i - 3) if random.random() < 0.3 else 0, i)
            a, b = perm[i], perm[p]
            adj[a].append(b)
            adj[b].append(a)
            edges.append((a, b))
        for l in adj:
            random.shuffle(l)
        return adj, edges

    def brute_par(adj, root):           # parent array via naive DFS
        n = len(adj)
        par = [-1] * n
        seen = [False] * n
        st = [root]
        seen[root] = True
        while st:
            v = st.pop()
            for u in adj[v]:
                if not seen[u]:
                    seen[u] = True
                    par[u] = v
                    st.append(u)
        return par

    def anc_list(par, v):
        res = [v]
        while par[res[-1]] >= 0:
            res.append(par[res[-1]])
        return res

    def all_dist(adj):
        return [bfs_dist(adj, s) for s in range(len(adj))]

    shapes = [None, None, "path", "star"]
    for it in range(300):
        n = random.randint(1, 40)
        adj, edges = rand_tree(n, shapes[it % 4])
        root = random.randrange(n)
        bp = brute_par(adj, root)
        dist = all_dist(adj)
        par, dep, order, sz = root_tree(adj, root)
        assert par == bp
        for v in range(n):
            assert dep[v] == dist[root][v]
            assert sz[v] == sum(1 for u in range(n) if v in anc_list(bp, u))
        # diameter / centers / centroids
        D, a, b = diameter(adj)
        assert D == max(max(r) for r in dist) and dist[a][b] == D
        ecc = [max(r) for r in dist]
        assert sorted(tree_centers(adj)) == [v for v in range(n) if ecc[v] == min(ecc)]
        cen = []
        for v in range(n):              # brute centroid: remove v, BFS components
            seen = [False] * n
            seen[v] = True
            mx = 0
            for s in adj[v]:
                q = [s]
                seen[s] = True
                for x in q:
                    for y in adj[x]:
                        if not seen[y]:
                            seen[y] = True
                            q.append(y)
                mx = max(mx, len(q))
            if 2 * mx <= n:
                cen.append(v)
        assert tree_centroids(adj) == cen
        # LCA x3
        L, E = LCA(adj, root), LCAEuler(adj, root)
        H = HLD(adj, root)
        qs = [(random.randrange(n), random.randrange(n)) for _ in range(30)]
        tl = tarjan_lca(adj, qs, root)
        for (u, v), t in zip(qs, tl):
            au = anc_list(bp, u)
            w = next(x for x in anc_list(bp, v) if x in au)
            assert L.lca(u, v) == w == E.lca(u, v) == t == H.lca(u, v)
            assert L.dist(u, v) == dist[u][v]
            k = random.randint(0, n)
            assert L.kth(u, k) == (au[k] if k < len(au) else -1)
        # Euler tin/tout
        tin, tout, ord_ = euler_tin_tout(adj, root)
        for v in range(n):
            assert ord_[tin[v]] == v
            for u in range(n):
                assert (tin[v] <= tin[u] < tout[v]) == (v in anc_list(bp, u))
        # sack
        col = [random.randrange(5) for _ in range(n)]
        sd = subtree_distinct_colors(adj, col, root)
        for v in range(n):
            assert sd[v] == len({col[u] for u in range(n) if v in anc_list(bp, u)})
        # HLD path sum
        vals = [random.randint(-9, 9) for _ in range(n)]
        P = HLDPathSum(adj, vals, root)
        for _ in range(30):
            if random.random() < 0.4:
                v = random.randrange(n)
                vals[v] = random.randint(-9, 9)
                P.set(v, vals[v])
            else:
                u, v = random.randrange(n), random.randrange(n)
                au, av = anc_list(bp, u), anc_list(bp, v)
                w = L.lca(u, v)
                path = au[:au.index(w) + 1] + av[:av.index(w)]
                assert P.query(u, v) == sum(vals[x] for x in path)
                es = H.segments(u, v, edge=True)
                assert sum(r - l + 1 for l, r in es) == len(path) - 1
        # centroid decomposition
        cpar, corder = centroid_decomp(adj)
        assert sorted(corder) == list(range(n)) and cpar.count(-1) == 1
        cdep = [len(anc_list(cpar, v)) for v in range(n)]
        assert max(cdep) <= n.bit_length()
        for u in range(n):              # path u-v must contain their centroid-tree LCA
            for v in range(n):
                cu = anc_list(cpar, u)
                c = next(x for x in anc_list(cpar, v) if x in cu)
                assert dist[u][c] + dist[c][v] == dist[u][v]
        K = random.randint(0, 6)
        assert count_pairs_dist_le(adj, K) == sum(1 for u in range(n) for v in range(u + 1, n)
                                                  if dist[u][v] <= K)
        # virtual tree
        keys = random.sample(range(n), random.randint(1, n))
        vs, vpar = virtual_tree(keys, tin, E.lca)
        closure = set(keys)
        for u in keys:
            for v in keys:
                closure.add(L.lca(u, v))
        assert set(vs) == closure
        for x in vs:
            ax = anc_list(bp, x)[1:]
            exp = next((y for y in ax if y in closure), -1)
            assert vpar[x] == exp
        # long-path dominant indices
        di = dominant_indices(adj, root)
        for v in range(n):
            cnt = [0] * n
            for u in range(n):
                if v in anc_list(bp, u):
                    cnt[dep[u] - dep[v]] += 1
            assert di[v] == cnt.index(max(cnt))

    # isomorphism vs brute (permutations) on small trees
    for it in range(150):
        n = random.randint(1, 7)
        a1, e1 = rand_tree(n)
        if random.random() < 0.5:
            pm = list(range(n))
            random.shuffle(pm)
            a2 = [[] for _ in range(n)]
            for x, y in e1:
                a2[pm[x]].append(pm[y])
                a2[pm[y]].append(pm[x])
            e2 = [(pm[x], pm[y]) for x, y in e1]
        else:
            a2, e2 = rand_tree(n)
        s1 = {frozenset(e) for e in e1}
        s2 = {frozenset(e) for e in e2}
        iso = any({frozenset((p[x], p[y])) for x, y in s1} == s2 for p in permutations(range(n)))
        ids = {}
        assert (unrooted_label(a1, ids) == unrooted_label(a2, ids)) == iso
        r1, r2 = random.randrange(n), random.randrange(n)
        dd1, dd2 = bfs_dist(a1, r1), bfs_dist(a2, r2)
        riso = any({frozenset((p[x], p[y])) for x, y in s1} == s2 and p[r1] == r2
                   for p in permutations(range(n)))
        assert (rooted_label(a1, r1, ids) == rooted_label(a2, r2, ids)) == riso

    # big deep trees: no recursion
    n = 200000
    adj, _ = rand_tree(n, "path")
    assert diameter(adj)[0] == n - 1
    L = LCA(adj)
    H = HLD(adj)
    assert dominant_indices(adj)[0] == (0 if len(adj[0]) == 1 else 1)
    assert len(centroid_decomp(adj)[1]) == n
    print("ALL OK")
