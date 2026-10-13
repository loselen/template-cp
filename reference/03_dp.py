# 03_dp.py -- Dynamic programming: fundamentals, classics, variants, optimizations.
# Python 3 stdlib only, PyPy-friendly. Snippets are self-contained (imports listed at top).
from bisect import bisect_left, bisect_right
from collections import deque
import heapq

INF = float('inf')

# ==== 5.1 Fundamentals: memoization vs tabulation ====
# Memo: @functools.cache (or lru_cache(None)) on f(state). Easy, but each call ~1us+, hashing
#   tuples is slow, and depth > ~1e4 needs sys.setrecursionlimit + bigger stack
#   (threading.stack_size(1<<27); run in a Thread) -- on PyPy deep recursion may still crash.
#   Call f.cache_clear() between test cases. Encode state as one int (i*M+j) for speed.
# Tabulation: loop states in topological order (dependencies first). Preferred in Python:
#   flat lists, no call overhead, enables rolling arrays. Use memo only when reachable
#   states are sparse or the order is awkward.
# State design: "what minimal info about the prefix/subproblem determines the future?"

# ==== 5.1 Pull vs push DP ====
# Pull: dp[i] = combine(dp[j] for predecessors j)    -- natural when predecessors are easy to list.
# Push: for each i (already final), dp[j] op= f(dp[i]) for successors j -- natural for
#   "from state i I can go to ...", probability/expectation forward, and skipping dp[i]==0/INF.

# ==== 5.1 Space optimization (rolling arrays) ====
# dp[i][*] depends only on dp[i-1][*] -> keep prev/cur (cur = [0]*m each row, or swap).
# 1D in-place: iterate j DESCENDING to use old values (0/1 knapsack), ASCENDING to reuse the
#   updated row (unbounded knapsack). Reconstruction needs full table (or Hirschberg).


# ==== 5.1 Reconstructing solutions: LCS string ====
# O(nm) time & memory. Returns one longest common subsequence (works for str or list).
def lcs_string(a, b):
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n):
        ai, row, nrow = a[i], dp[i], dp[i + 1]
        for j in range(m):
            nrow[j + 1] = row[j] + 1 if ai == b[j] else max(row[j + 1], nrow[j])
    res, i, j = [], n, m
    while i and j:                       # walk back from (n,m)
        if a[i - 1] == b[j - 1]:
            res.append(a[i - 1]); i -= 1; j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1
        else:
            j -= 1
    res.reverse()
    return ''.join(res) if isinstance(a, str) else res


# ==== 5.2 LCS length (rolling) ====
# O(nm) time, O(m) memory.
def lcs_len(a, b):
    m = len(b)
    prev = [0] * (m + 1)
    for x in a:
        cur = [0] * (m + 1)
        for j in range(m):
            cur[j + 1] = prev[j] + 1 if x == b[j] else max(prev[j + 1], cur[j])
        prev = cur
    return prev[m]


# ==== 5.2 Fibonacci-style recurrences / counting paths / grid DP ====
# Count right/down paths from (0,0) to (n-1,m-1) avoiding '#'. O(nm). Fibonacci-like
# linear recurrences with huge n -> see 09_algebra_poly.py: mat_pow.
def grid_paths(grid, mod=None):
    m = len(grid[0])
    dp = [0] * m
    dp[0] = 1 if grid[0][0] != '#' else 0
    for r, row in enumerate(grid):
        for c in range(m):
            if row[c] == '#':
                dp[c] = 0
            elif c:
                dp[c] += dp[c - 1]       # dp[c] (from above) + left
                if mod: dp[c] %= mod
    return dp[m - 1]


# Min-cost right/down path sum. O(nm).
def grid_min_path(cost):
    m = len(cost[0])
    dp = [INF] * m
    dp[0] = 0
    for row in cost:
        dp[0] += row[0]
        for c in range(1, m):
            dp[c] = min(dp[c], dp[c - 1]) + row[c]
    return dp[m - 1]


# ==== 5.2 Coin change ====
# Number of ways (order doesn't matter: coins outer loop). O(len(coins)*T).
# For ordered sequences (compositions) swap loops: for t: for c.
def coin_ways(coins, T, mod=None):
    dp = [0] * (T + 1)
    dp[0] = 1
    for c in coins:
        for t in range(c, T + 1):
            dp[t] += dp[t - c]
            if mod: dp[t] %= mod
    return dp[T]


# Min number of coins to make T (unbounded); -1 if impossible. O(len(coins)*T).
def coin_min(coins, T):
    dp = [0] + [INF] * T
    for c in coins:
        for t in range(c, T + 1):
            if dp[t - c] + 1 < dp[t]:
                dp[t] = dp[t - c] + 1
    return dp[T] if dp[T] != INF else -1


# ==== 5.2 Knapsack ====
# 0/1 knapsack, items=[(w,v)], max value with weight <= W. O(nW). j descending!
# For "exactly W" init dp=[0]+[-INF]*W.
def knap01(items, W):
    dp = [0] * (W + 1)
    for w, v in items:
        for j in range(W, w - 1, -1):
            if dp[j - w] + v > dp[j]:
                dp[j] = dp[j - w] + v
    return dp[W]


# Unbounded knapsack. O(nW). j ascending.
def knap_unbounded(items, W):
    dp = [0] * (W + 1)
    for w, v in items:
        for j in range(w, W + 1):
            if dp[j - w] + v > dp[j]:
                dp[j] = dp[j - w] + v
    return dp[W]


# Bounded knapsack via binary splitting, items=[(w,v,cnt)]. O(W * sum log cnt).
def knap_bounded_binary(items, W):
    dp = [0] * (W + 1)
    for w, v, c in items:
        k = 1
        while c > 0:
            t = min(k, c); c -= t; k <<= 1
            tw, tv = t * w, t * v
            for j in range(W, tw - 1, -1):
                if dp[j - tw] + tv > dp[j]:
                    dp[j] = dp[j - tw] + tv
    return dp[W]


# Bounded knapsack via monotonic queue per residue mod w. O(nW).
# new[r+k*w] = k*v + max_{k-c<=t<=k} (dp[r+t*w] - t*v)
def knap_bounded_mq(items, W):
    dp = [0] * (W + 1)
    for w, v, c in items:
        new = dp[:]
        for r in range(min(w, W + 1)):
            q = deque()                  # (t, value) with decreasing value
            t = 0
            for j in range(r, W + 1, w):
                val = dp[j] - t * v
                while q and q[-1][1] <= val:
                    q.pop()
                q.append((t, val))
                if q[0][0] < t - c:
                    q.popleft()
                new[j] = q[0][1] + t * v
                t += 1
        dp = new
    return dp[W]


# ==== 5.2 Longest increasing subsequence ====
# O(n log n). Returns list of indices of one LIS. strict=False -> non-decreasing.
# O(n^2) version: dp[i] = 1 + max(dp[j] for j<i if a[j]<a[i]).
def lis(a, strict=True):
    find = bisect_left if strict else bisect_right
    tails, tidx, par = [], [], [-1] * len(a)
    for i, x in enumerate(a):
        p = find(tails, x)
        if p == len(tails):
            tails.append(x); tidx.append(i)
        else:
            tails[p] = x; tidx[p] = i
        par[i] = tidx[p - 1] if p else -1
    res, i = [], tidx[-1] if tidx else -1
    while i != -1:
        res.append(i); i = par[i]
    return res[::-1]


# ==== 5.2 Edit distance (Levenshtein) ====
# O(nm) time, O(m) memory. Insert/delete/substitute cost 1.
def edit_distance(a, b):
    m = len(b)
    prev = list(range(m + 1))
    for i, x in enumerate(a, 1):
        cur = [i] + [0] * m
        for j in range(1, m + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != b[j - 1]))
        prev = cur
    return prev[m]


# ==== 5.2 Matrix chain multiplication ====
# dims: matrix i is dims[i] x dims[i+1]. O(n^3). Min scalar multiplications.
def matrix_chain(dims):
    n = len(dims) - 1
    dp = [[0] * n for _ in range(n)]
    for length in range(2, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1
            di, dj = dims[i], dims[j + 1]
            dp[i][j] = min(dp[i][k] + dp[k + 1][j] + di * dims[k + 1] * dj for k in range(i, j))
    return dp[0][n - 1] if n else 0


# ==== 5.2 Subset sum with Python-int bitset ====
# Bit s of result set <=> some subset sums to s. O(n * S / 64) -- very fast.
# Mask with (1 << (S+1)) - 1 if only sums <= S matter (keeps ints small).
def subset_sums(a, S=None):
    bits = 1
    mask = (1 << (S + 1)) - 1 if S is not None else -1
    for x in a:
        bits |= bits << x
        if S is not None: bits &= mask
    return bits                          # test: bits >> s & 1


# ==== 5.2 Partition problems ====
# Min |sum(A) - sum(B)| over 2-partitions. O(n*sum/64).
# Many equal values: group counts + binary splitting (as knap_bounded_binary) first.
def min_partition_diff(a):
    tot = sum(a)
    bits = subset_sums(a, tot // 2)
    return tot - 2 * (bits.bit_length() - 1)   # highest reachable sum <= tot//2


# ==== 5.3 Interval DP ====
# Merge adjacent piles, cost = size of merged pile; min total. O(n^3). (Knuth -> O(n^2) below.)
# Pattern: for length: for i: j=i+len-1: dp[i][j] = min over split k.
def merge_stones(a):
    n = len(a)
    P = [0] * (n + 1)
    for i, x in enumerate(a): P[i + 1] = P[i] + x
    dp = [[0] * n for _ in range(n)]
    for length in range(2, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1
            dp[i][j] = min(dp[i][k] + dp[k + 1][j] for k in range(i, j)) + P[j + 1] - P[i]
    return dp[0][n - 1] if n else 0


# Burst balloons: bursting k gains l*k*r of current neighbours. O(n^3).
# Trick: choose the LAST balloon k burst in open interval (l, r).
def burst_balloons(a):
    v = [1] + a + [1]
    n = len(v)
    dp = [[0] * n for _ in range(n)]
    for length in range(2, n):
        for l in range(n - length):
            r = l + length
            vl, vr, best = v[l], v[r], 0
            for k in range(l + 1, r):
                t = dp[l][k] + dp[k][r] + vl * v[k] * vr
                if t > best: best = t
            dp[l][r] = best
    return dp[0][n - 1]


# ==== 5.3 Tree DP (iterative): max weight independent set ====
# O(n). adj: adjacency lists, w: weights (default 1). Process vertices in reverse BFS order.
def tree_mis(n, adj, w=None, root=0):
    w = w or [1] * n
    par = [-1] * n
    order = [root]
    par[root] = root
    for v in order:
        for c in adj[v]:
            if par[c] == -1 and c != root:
                par[c] = v; order.append(c)
    take, skip = w[:], [0] * n
    for v in reversed(order):
        if v != root:
            p = par[v]
            take[p] += skip[v]
            skip[p] += max(take[v], skip[v])
    return max(take[root], skip[root])


# ==== 5.3 Rerooting (all-roots DP), generic iterative ====
# O(n * cost(merge)). merge must be associative with identity e.
#   lift(acc, v): acc = merge of neighbour values around v (excluding the parent side)
#                 -> value of "subtree hanging at v" as seen from its parent.
# Returns res[v] = merge of ALL neighbour values of v (finalize yourself, e.g. lift(res[v], v)).
# Example (sum of distances): e=(0,0), merge=add pairwise, lift=(s+1, d+s+1), ans=res[v][1].
def reroot(n, adj, e, merge, lift, root=0):
    par = [-1] * n
    par[root] = root
    order = [root]
    for v in order:
        for c in adj[v]:
            if par[c] == -1:
                par[c] = v; order.append(c)
    par[root] = -1
    down = [e] * n
    for v in reversed(order):
        acc = e
        for c in adj[v]:
            if c != par[v]:
                acc = merge(acc, down[c])
        down[v] = lift(acc, v)
    up = [e] * n                          # up[v]: value of the parent-side part, seen from v
    res = [e] * n
    for v in order:
        nb = adj[v]
        pv = par[v]
        vals = [up[v] if c == pv else down[c] for c in nb]
        k = len(nb)
        pre = [e] * (k + 1)
        for i in range(k):
            pre[i + 1] = merge(pre[i], vals[i])
        res[v] = pre[k]
        suf = e
        for i in range(k - 1, -1, -1):
            c = nb[i]
            if c != pv:
                up[c] = lift(merge(pre[i], suf), v)
            suf = merge(vals[i], suf)
    return res


# ==== 5.3 DAG DP ====
# Any DP is a DAG over states. Explicit DAG: topological order (Kahn, see 04_graphs.py),
# then dp[v] = f(dp[u] for u->v) (longest path, #paths mod p). Cycles -> SCC condense first.

# ==== 5.3 Bitmask DP ====
# TSP: min Hamiltonian cycle from 0. O(2^n n^2). dist: n x n matrix.
# Hamiltonian path variant: init dp[1<<v][v]=0 for all v, answer min dp[full][v].
def tsp(dist):
    n = len(dist)
    full = 1 << n
    dp = [[INF] * n for _ in range(full)]
    dp[1][0] = 0
    for mask in range(1, full, 2):       # masks containing 0
        row = dp[mask]
        for v in range(n):
            cur = row[v]
            if cur == INF or not (mask >> v & 1): continue
            dv = dist[v]
            for u in range(n):
                if not (mask >> u & 1):
                    nm = mask | 1 << u
                    if cur + dv[u] < dp[nm][u]:
                        dp[nm][u] = cur + dv[u]
    return min(dp[full - 1][v] + dist[v][0] for v in range(n)) if n > 1 else 0


# Assignment: worker popcount(mask) takes a job in mask. O(2^n n). (Hungarian: 05_flow_matching.py)
def assignment_bitmask(cost):
    n = len(cost)
    dp = [INF] * (1 << n)
    dp[0] = 0
    for mask in range(1 << n):
        cur = dp[mask]
        i = bin(mask).count('1')
        if cur == INF or i == n: continue
        row = cost[i]
        for j in range(n):
            if not mask >> j & 1:
                nm = mask | 1 << j
                if cur + row[j] < dp[nm]:
                    dp[nm] = cur + row[j]
    return dp[(1 << n) - 1]


# ==== 5.3 Digit DP ====
# Count x in [0, N] with digit sum == S. Iterative with tight flag: O(len(N) * S * 10).
# Range [L,R]: f(R) - f(L-1). Generic: state = (pos, tight, extra); non-tight states pooled.
def count_digit_sum(N, S):
    if N < 0: return 0
    cnt = [0] * (S + 1)                  # non-tight prefixes by digit sum
    tsum = 0                             # digit sum of the tight prefix
    for d in map(int, str(N)):
        new = [0] * (S + 1)
        for s in range(S + 1):
            c = cnt[s]
            if c:
                for x in range(min(9, S - s) + 1):
                    new[s + x] += c
        for x in range(d):               # leave tight by placing x < d
            if tsum + x <= S: new[tsum + x] += 1
        tsum += d
        cnt = new
    return cnt[S] + (tsum == S)


# ==== 5.3 Probability / expected value DP ====
# Expected number of fair-die rolls until position >= n (start 0). O(6n).
# E[i] = 1 + sum_d E[i+d]/6 computed backwards (states with no future first).
# Mod p version: replace /6 with * pow(6, p-2, p). Cyclic dependencies -> Gaussian elimination.
def expected_rolls(n, frac=False):
    from fractions import Fraction
    one = Fraction(1) if frac else 1.0
    E = [one * 0] * (n + 6)
    for i in range(n - 1, -1, -1):
        E[i] = one + sum(E[i + 1:i + 7]) / 6
    return E[0]


# ==== 5.3 Game DP ====
# Players alternately take from either end of a; returns (first - second) under optimal play.
# O(n^2). General: win[state] = any(not win[next]) over moves (Sprague-Grundy: see 08_math_nt_comb.py 9.6).
def ends_game(a):
    n = len(a)
    dp = a[:]                            # dp[i] for intervals of current length starting at i
    for length in range(2, n + 1):
        dp = [max(a[i] - dp[i + 1], a[i + length - 1] - dp[i]) for i in range(n - length + 1)]
    return dp[0] if n else 0


# ==== 5.3 DP on subsequences and strings ====
# Count distinct subsequences: dp = 2*dp - dp_before_last_occurrence. O(n). (includes empty)
def count_distinct_subseq(s, mod=None):
    dp, last = 1, {}
    for ch in s:
        prev = dp
        dp = 2 * dp - last.get(ch, 0)
        if mod: dp %= mod
        last[ch] = prev
    return dp


# ==== 5.3 Broken profile DP: domino tilings of n x m ====
# Cell-by-cell, mask bit j = cell in column j of the frontier already covered. O(n m 2^m).
def domino_tilings(n, m):
    if m > n: n, m = m, n
    dp = [0] * (1 << m)
    dp[0] = 1
    for i in range(n):
        for j in range(m):
            new = [0] * (1 << m)
            bj = 1 << j
            for mask in range(1 << m):
                c = dp[mask]
                if not c: continue
                if mask & bj:
                    new[mask ^ bj] += c                      # already covered from above
                else:
                    if i + 1 < n: new[mask | bj] += c          # vertical, covers (i+1,j)
                    if j + 1 < m and not mask >> (j + 1) & 1:
                        new[mask | bj << 1] += c               # horizontal, covers (i,j+1)
            dp = new
    return dp[0]


# ==== 5.3 SOS DP (zeta / Mobius over subsets) ====
# In place, O(n 2^n). len(f) == 1<<n.
def sos_subset(f, n):                    # f[S] = sum_{T subset S} f[T]
    for i in range(n):
        b = 1 << i
        for S in range(1 << n):
            if S & b: f[S] += f[S ^ b]
    return f


def sos_superset(f, n):                  # f[S] = sum_{T superset S} f[T]
    for i in range(n):
        b = 1 << i
        for S in range(1 << n):
            if not S & b: f[S] += f[S | b]
    return f


def mobius_subset(f, n):                 # inverse of sos_subset
    for i in range(n):
        b = 1 << i
        for S in range(1 << n):
            if S & b: f[S] -= f[S ^ b]
    return f
# Subset convolution / FWHT: see 09_algebra_poly.py. Enumerate submasks: s=m; while s: ...; s=(s-1)&m (O(3^n) total).


# ==== 5.3 DP on permutations: insertion / "connected components" DP ====
# Insertion DP: build permutation by inserting values 1..n (or positions) one at a time;
#   state tracks only what matters (e.g. # descents, last relative rank -> dp[i][j] with prefix sums).
# Components DP: insert values in increasing order; state = #components (blocks to be glued).
#   New value: open new component, append to an end of a component, or merge two.
# Example (CEOI16 kangaroo): count zigzag permutations of 1..n with p[0]=s, p[-1]=t. O(n^2).
def zigzag_perms(n, s, t, mod=None):
    dp = [1] + [0] * (n + 1)            # dp[j]: j components
    for i in range(1, n + 1):
        new = [0] * (n + 2)
        for j in range(1, i + 1):
            if i == s or i == t:        # endpoint: new comp at its end, or attach to comp at that end
                v = dp[j - 1] + dp[j]
            else:                       # new comp in free gaps, or merge two (zigzag forbids append)
                v = dp[j - 1] * (j - (i > s) - (i > t)) + dp[j + 1] * j
            new[j] = v % mod if mod else v
        dp = new
    return dp[1]


# ==== 5.3 DP with data structures (BIT transitions) ====
# Count strictly increasing subsequences (non-empty). dp[i] = 1 + sum dp[j], j<i, a[j]<a[i];
# Fenwick over compressed values -> O(n log n). Same pattern: max-weight LIS with a max-BIT,
# or range-of-j transitions with a segment tree (02_data_structures.py).
def count_increasing_subseq(a, mod=None):
    vals = sorted(set(a))
    rk = {v: i + 1 for i, v in enumerate(vals)}
    m = len(vals)
    bit = [0] * (m + 1)
    total = 0
    for x in a:
        i = rk[x] - 1                    # prefix sum over ranks < x
        s = 1
        while i > 0:
            s += bit[i]; i &= i - 1
        if mod: s %= mod
        total += s
        i = rk[x]
        while i <= m:
            bit[i] += s
            if mod: bit[i] %= mod
            i += i & -i
    return total % mod if mod else total


# ==== 5.3 Counting DP with inclusion-exclusion ====
# #objects avoiding all "bad" properties = sum_{S subset bad} (-1)^|S| #(objects having all of S).
# DP form: sort bad points, g[i] = ways(start->p_i) - sum_{j before i} g[j]*ways(p_j->p_i)
#   (first-bad-point decomposition, e.g. grid paths avoiding k blocked cells in O(k^2)).
# "Exactly k" from "at least k": e_k = sum_{j>=k} (-1)^(j-k) C(j,k) a_j (binomial inversion).

# ==== 5.3 Plug DP (connectivity profile DP) ====
# Broken profile where the mask stores, per frontier edge, a plug id / bracket type
# (0 none, 1 '(' , 2 ')') encoding which plugs connect -> Hamiltonian cycles, single loop
# tiling on grids. States kept in dict {encoded_profile: count}; normalize labels each step.
# Cases at each cell: (none,none) open pair; one plug -> extend; two plugs -> merge (fix match).


# ==== 5.4 Prefix sum optimization ====
# dp[i][j] = sum_{k=j-L..j} dp[i-1][k] -> keep prefix sums of row i-1: O(1) per state.
# Same for min/max over sliding ranges (monotonic queue) and for "sum over all j<i" (running total).

# ==== 5.4 Monotonic queue optimization ====
# dp[i] = c[i] + min(dp[i-k..i-1]): choose positions so every window of k consecutive has
# a chosen one; min total cost. O(n). Answer = min over last k dp's.
def min_cost_cover_windows(c, k):
    n = len(c)
    dp = [0] * n
    q = deque()                          # indices with increasing dp
    for i in range(n):
        while q and q[0] < i - k:
            q.popleft()
        dp[i] = c[i] + (dp[q[0]] if i >= k else 0)   # first pick must be < k
        while q and dp[q[-1]] >= dp[i]:
            q.pop()
        q.append(i)
    return min(dp[max(0, n - k):])


# ==== 5.4 Convex Hull Trick (monotone slopes) ====
# Min of lines y = m*x + b. Add lines in DECREASING slope order (for max: negate m,b and result).
# add amortized O(1); query: qmono O(1) amortized for non-decreasing x, query O(log n) any x.
class CHT:
    def __init__(self):
        self.M, self.B, self.p = [], [], 0

    def add(self, m, b):
        M, B = self.M, self.B
        if M and M[-1] == m:
            if b >= B[-1]: return
            M.pop(); B.pop()
        while len(M) >= 2 and (b - B[-2]) * (M[-2] - M[-1]) <= (B[-1] - B[-2]) * (M[-2] - m):
            M.pop(); B.pop()
        M.append(m); B.append(b)

    def query(self, x):                  # binary search, any x
        M, B = self.M, self.B
        lo, hi = 0, len(M) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if M[mid] * x + B[mid] >= M[mid + 1] * x + B[mid + 1]: lo = mid + 1
            else: hi = mid
        return M[lo] * x + B[lo]

    def qmono(self, x):                  # x non-decreasing across calls
        M, B = self.M, self.B
        p = min(self.p, len(M) - 1)
        while p + 1 < len(M) and M[p + 1] * x + B[p + 1] <= M[p] * x + B[p]:
            p += 1
        self.p = p
        return M[p] * x + B[p]
# Arbitrary slope order / online: Li Chao tree -> see 02_data_structures.py: LiChao


# ==== 5.4 Divide and conquer optimization ====
# new[i] = min_{j<i} prev[j] + cost(j, i), when opt(i) is monotone (cost satisfies QI).
# O(n log n) per layer with explicit stack (no recursion). cost(j,i): segment (j, i].
def dc_layer(prev, cost, n):
    new = [INF] * (n + 1)
    st = [(1, n, 0, n - 1)]
    while st:
        lo, hi, ol, orr = st.pop()
        if lo > hi: continue
        mid = (lo + hi) // 2
        best, bj = INF, ol
        for j in range(ol, min(mid - 1, orr) + 1):
            v = prev[j] + cost(j, mid)
            if v < best: best, bj = v, j
        new[mid] = best
        st.append((lo, mid - 1, ol, bj))
        st.append((mid + 1, hi, bj, orr))
    return new


# Example: split a into k contiguous groups minimizing sum of (group sum)^2. O(k n log n).
def partition_k_squares(a, k):
    n = len(a)
    P = [0] * (n + 1)
    for i, x in enumerate(a): P[i + 1] = P[i] + x
    cost = lambda j, i: (P[i] - P[j]) ** 2
    dp = [0] + [INF] * n
    for _ in range(k):
        dp = dc_layer(dp, cost, n)
        dp[0] = INF
    return dp[n]


# ==== 5.4 Knuth optimization ====
# dp[i][j] = min_{i<=k<j} dp[i][k]+dp[k+1][j] + C(i,j), with opt[i][j-1] <= opt[i][j] <= opt[i+1][j]
# (C satisfies QI + monotone). O(n^2). Example: merging adjacent piles (cf. merge_stones).
# Optimal BST is identical in shape (C(i,j) = sum of frequencies in [i,j]).
def merge_stones_knuth(a):
    n = len(a)
    if n == 0: return 0
    P = [0] * (n + 1)
    for i, x in enumerate(a): P[i + 1] = P[i] + x
    dp = [[0] * n for _ in range(n)]
    opt = [[0] * n for _ in range(n)]
    for i in range(n): opt[i][i] = i
    for length in range(2, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1
            best, bk = INF, i
            dpi = dp[i]
            for k in range(opt[i][j - 1], min(opt[i + 1][j], j - 1) + 1):
                v = dpi[k] + dp[k + 1][j]
                if v < best: best, bk = v, k
            dp[i][j] = best + P[j + 1] - P[i]
            opt[i][j] = bk
    return dp[0][n - 1]


# ==== 5.4 Aliens trick / WQS binary search / Lagrangian relaxation ====
# Need g(k) = best with EXACTLY k parts, g convex in k. Solve h(lam) = min (cost + lam*#parts)
# without k (tie-break: fewest parts); find smallest lam with cnt <= k; answer h - lam*k.
# Integer costs -> integer lam suffices. Inner DP here O(n^2) (use CHT for O(n)).
# Example: same as partition_k_squares (a positive), total O(n^2 log C).
def aliens_partition(a, k):
    n = len(a)
    P = [0] * (n + 1)
    for i, x in enumerate(a): P[i + 1] = P[i] + x

    def solve(lam):                      # returns (value, #parts) minimal value, min parts
        dp = [0] * (n + 1); cnt = [0] * (n + 1)
        for i in range(1, n + 1):
            best, bc = INF, 0
            for j in range(i):
                v = dp[j] + (P[i] - P[j]) ** 2 + lam
                if v < best or (v == best and cnt[j] + 1 < bc):
                    best, bc = v, cnt[j] + 1
            dp[i], cnt[i] = best, bc
        return dp[n], cnt[n]

    lo, hi = 0, P[n] ** 2 + 1            # lam=hi forces 1 part
    while lo < hi:
        mid = (lo + hi) // 2
        if solve(mid)[1] <= k: hi = mid
        else: lo = mid + 1
    v, _ = solve(lo)
    return v - lo * k


# ==== 5.4 Slope trick ====
# Min sum |a_i - b_i| with b non-decreasing. Max-heap of breakpoints. O(n log n).
# (Strictly increasing: apply to a_i - i.) Generalize: keep L (max-heap) / R (min-heap) + lazy add.
def min_cost_nondecreasing(a):
    h, ans = [], 0
    for x in a:
        heapq.heappush(h, -x)
        if -h[0] > x:
            ans += -h[0] - x
            heapq.heapreplace(h, -x)
    return ans


# ==== 5.4 Matrix exponentiation for linear recurrences ====
# see 09_algebra_poly.py: mat_pow (and Berlekamp-Massey / Kitamasa there).

# ==== 5.4 SMAWK / monotone matrices (conceptual) ====
# Row minima of a totally monotone n x m matrix in O(n+m) (REDUCE columns to n, recurse on
# even rows, interpolate odd rows). Rarely needed: D&C opt gives O((n+m) log n) for the same
# monotone-argmin property and is far simpler. LARSCH handles online dp[i] = min_j dp[j]+w(j,i).

# ==== 5.4 Bitset DP speedups ====
# Python int = bitset: shifts/or/and are C-speed over whole words (~n/64 ops).
#   subset sum: bits |= bits << x;  reachability/transitive closure: reach[v] |= reach[u];
#   LCS-style bit-parallel (Hyyro): O(n m / 64). Count bits: int.bit_count() (3.10+) or bin(x).count('1').
#   Keep ints bounded by masking, otherwise they grow without limit.

# ==== 5.4 Open and close interval trick ====
# Sort items; process left to right, state = #currently "open" groups (+ accumulated cost).
# Cost of a group = max - min = sum of gaps it spans, so when moving from a[i-1] to a[i]
# add (a[i]-a[i-1]) * open_count. Transitions per item: open new group / close / add to
# an open one (x open choices) / singleton. Makes "sum over groups of (max-min) <= K" O(n^2 K).


# =====================================================================
if __name__ == "__main__":
    import random
    from itertools import combinations, permutations, product
    from functools import lru_cache
    rnd = random.Random(12345)

    def is_subseq(s, t):
        it = iter(t)
        return all(c in it for c in s)

    # LCS
    for _ in range(300):
        a = ''.join(rnd.choice('abc') for _ in range(rnd.randint(0, 7)))
        b = ''.join(rnd.choice('abc') for _ in range(rnd.randint(0, 7)))
        best = max((k for k in range(len(a) + 1) for c in combinations(a, k)
                    if is_subseq(c, b)), default=0)
        s = lcs_string(a, b)
        assert len(s) == best == lcs_len(a, b) and is_subseq(s, a) and is_subseq(s, b)

    # grid paths
    for _ in range(200):
        n, m = rnd.randint(1, 5), rnd.randint(1, 5)
        g = [''.join(rnd.choice('..#') for _ in range(m)) for _ in range(n)]
        cst = [[rnd.randint(0, 9) for _ in range(m)] for _ in range(n)]

        @lru_cache(None)
        def cnt(r, c):
            if r >= n or c >= m or g[r][c] == '#': return 0
            if (r, c) == (n - 1, m - 1): return 1
            return cnt(r + 1, c) + cnt(r, c + 1)

        @lru_cache(None)
        def mc(r, c):
            if r >= n or c >= m: return INF
            if (r, c) == (n - 1, m - 1): return cst[r][c]
            return cst[r][c] + min(mc(r + 1, c), mc(r, c + 1))
        assert grid_paths(g) == cnt(0, 0) and grid_min_path(cst) == mc(0, 0)

    # coins
    for _ in range(200):
        coins = list(set(rnd.randint(1, 8) for _ in range(rnd.randint(1, 4))))
        T = rnd.randint(0, 25)
        ways, mn = 0, INF
        for combo in product(*[range(T // c + 1) for c in coins]):
            if sum(x * c for x, c in zip(combo, coins)) == T:
                ways += 1; mn = min(mn, sum(combo))
        assert coin_ways(coins, T) == ways
        assert coin_min(coins, T) == (mn if mn != INF else -1)

    # knapsacks
    for _ in range(200):
        n = rnd.randint(0, 6); W = rnd.randint(0, 20)
        items = [(rnd.randint(1, 7), rnd.randint(0, 10), rnd.randint(0, 4)) for _ in range(n)]
        best = 0
        for mask in range(1 << n):
            w = sum(items[i][0] for i in range(n) if mask >> i & 1)
            if w <= W: best = max(best, sum(items[i][1] for i in range(n) if mask >> i & 1))
        assert knap01([(w, v) for w, v, _ in items], W) == best
        bb = 0
        for combo in product(*[range(c + 1) for _, _, c in items]):
            if sum(x * it[0] for x, it in zip(combo, items)) <= W:
                bb = max(bb, sum(x * it[1] for x, it in zip(combo, items)))
        assert knap_bounded_binary(items, W) == bb == knap_bounded_mq(items, W)
        ub = knap_bounded_binary([(w, v, W) for w, v, _ in items], W)
        assert knap_unbounded([(w, v) for w, v, _ in items], W) == ub

    # LIS
    for _ in range(300):
        a = [rnd.randint(0, 6) for _ in range(rnd.randint(0, 12))]
        for strict in (True, False):
            dp = []
            for i in range(len(a)):
                dp.append(1 + max([dp[j] for j in range(i)
                                   if (a[j] < a[i] if strict else a[j] <= a[i])], default=0))
            r = lis(a, strict)
            assert len(r) == max(dp, default=0)
            assert all(r[i] < r[i + 1] and (a[r[i]] < a[r[i + 1]] if strict else a[r[i]] <= a[r[i + 1]])
                       for i in range(len(r) - 1))

    # edit distance (vs memo brute)
    for _ in range(200):
        a = ''.join(rnd.choice('ab') for _ in range(rnd.randint(0, 6)))
        b = ''.join(rnd.choice('ab') for _ in range(rnd.randint(0, 6)))

        @lru_cache(None)
        def ed(i, j):
            if i == len(a): return len(b) - j
            if j == len(b): return len(a) - i
            return min(ed(i + 1, j) + 1, ed(i, j + 1) + 1, ed(i + 1, j + 1) + (a[i] != b[j]))
        assert edit_distance(a, b) == ed(0, 0)

    # matrix chain (vs recursion over all parenthesizations)
    for _ in range(100):
        d = [rnd.randint(1, 10) for _ in range(rnd.randint(2, 7))]

        def mcb(i, j):
            if i == j: return 0
            return min(mcb(i, k) + mcb(k + 1, j) + d[i] * d[k + 1] * d[j + 1] for k in range(i, j))
        assert matrix_chain(d) == mcb(0, len(d) - 2)

    # subset sum / partition
    for _ in range(200):
        a = [rnd.randint(1, 15) for _ in range(rnd.randint(0, 8))]
        sums = {sum(c) for k in range(len(a) + 1) for c in combinations(a, k)}
        bits = subset_sums(a)
        assert {s for s in range(sum(a) + 1) if bits >> s & 1} == sums
        assert min_partition_diff(a) == min(abs(sum(a) - 2 * s) for s in sums)
        S = rnd.randint(0, 30)
        b2 = subset_sums(a, S)
        assert b2.bit_length() <= S + 1 and all((b2 >> s & 1) == (s in sums) for s in range(S + 1))

    # interval DP + Knuth
    for _ in range(200):
        a = [rnd.randint(1, 20) for _ in range(rnd.randint(0, 9))]

        def mb(l):
            if len(l) <= 1: return 0
            return min(mb(l[:i] + [l[i] + l[i + 1]] + l[i + 2:]) + l[i] + l[i + 1]
                       for i in range(len(l) - 1))
        if len(a) <= 6: assert merge_stones(a) == mb(a)
        assert merge_stones(a) == merge_stones_knuth(a)

        b = a[:6]

        def bb_(l):
            if not l: return 0
            return max(bb_(l[:i] + l[i + 1:]) + (l[i - 1] if i else 1) * l[i] *
                       (l[i + 1] if i + 1 < len(l) else 1) for i in range(len(l)))
        assert burst_balloons(b) == bb_(b)

    # tree MIS + rerooting
    for _ in range(150):
        n = rnd.randint(1, 11)
        adj = [[] for _ in range(n)]
        for v in range(1, n):
            p = rnd.randrange(v); adj[p].append(v); adj[v].append(p)
        w = [rnd.randint(1, 9) for _ in range(n)]
        best = 0
        for mask in range(1 << n):
            if all(not (mask >> v & 1 and mask >> u & 1) for v in range(n) for u in adj[v]):
                best = max(best, sum(w[v] for v in range(n) if mask >> v & 1))
        assert tree_mis(n, adj, w, rnd.randrange(n)) == best
        res = reroot(n, adj, (0, 0), lambda x, y: (x[0] + y[0], x[1] + y[1]),
                     lambda acc, v: (acc[0] + 1, acc[1] + acc[0] + 1), rnd.randrange(n))
        for s in range(n):
            dist = [-1] * n; dist[s] = 0; q = [s]
            for v in q:
                for u in adj[v]:
                    if dist[u] < 0: dist[u] = dist[v] + 1; q.append(u)
            assert res[s] == (n - 1, sum(dist))

    # bitmask: TSP, assignment
    for _ in range(60):
        n = rnd.randint(1, 6)
        d = [[rnd.randint(1, 20) if i != j else 0 for j in range(n)] for i in range(n)]
        bt = min((sum(d[p[i]][p[(i + 1) % n]] for i in range(n)) for p in permutations(range(n))
                  if p[0] == 0), default=0) if n > 1 else 0
        assert tsp(d) == bt
        assert assignment_bitmask(d) == min(sum(d[i][p[i]] for i in range(n))
                                            for p in permutations(range(n)))

    # digit DP
    for _ in range(200):
        N = rnd.randint(-1, 3000); S = rnd.randint(0, 30)
        assert count_digit_sum(N, S) == sum(1 for x in range(N + 1) if sum(map(int, str(x))) == S)

    # expected value (vs explicit enumeration over roll sequences)
    from fractions import Fraction
    for n in range(1, 7):
        def ex(pos, depth):
            if pos >= n: return Fraction(depth)
            return sum(ex(pos + d, depth + 1) for d in range(1, 7)) / 6
        assert expected_rolls(n, True) == ex(0, 0)
        assert abs(expected_rolls(n) - float(ex(0, 0))) < 1e-9

    # game DP (vs explicit minimax on scores)
    for _ in range(200):
        a = [rnd.randint(0, 20) for _ in range(rnd.randint(0, 9))]

        def play(l, r, turn):           # returns (score0, score1)
            if l > r: return (0, 0)
            opts = []
            for take, nl, nr in ((a[l], l + 1, r), (a[r], l, r - 1)):
                s = list(play(nl, nr, 1 - turn)); s[turn] += take; opts.append(s)
            return max(opts, key=lambda s: s[turn] - s[1 - turn])
        s = play(0, len(a) - 1, 0)
        assert ends_game(a) == s[0] - s[1]

    # distinct subsequences
    for _ in range(200):
        s = ''.join(rnd.choice('abc') for _ in range(rnd.randint(0, 8)))
        allsub = {''.join(c) for k in range(len(s) + 1) for c in combinations(s, k)}
        assert count_distinct_subseq(s) == len(allsub)

    # domino tilings
    def brute_tile(n, m):
        g = [[0] * m for _ in range(n)]

        def go(k):
            if k == n * m: return 1
            r, c = divmod(k, m)
            if g[r][c]: return go(k + 1)
            t = 0
            if c + 1 < m and not g[r][c + 1]:
                g[r][c] = g[r][c + 1] = 1; t += go(k + 1); g[r][c] = g[r][c + 1] = 0
            if r + 1 < n:
                g[r][c] = g[r + 1][c] = 1; t += go(k + 1); g[r][c] = g[r + 1][c] = 0
            return t
        return go(0)
    for n in range(1, 6):
        for m in range(1, 6):
            assert domino_tilings(n, m) == brute_tile(n, m)
    assert domino_tilings(8, 8) == 12988816

    # SOS
    for _ in range(50):
        n = rnd.randint(0, 6)
        f = [rnd.randint(-9, 9) for _ in range(1 << n)]
        sub = sos_subset(f[:], n); sup = sos_superset(f[:], n)
        for S in range(1 << n):
            assert sub[S] == sum(f[T] for T in range(1 << n) if T & S == T)
            assert sup[S] == sum(f[T] for T in range(1 << n) if T & S == S)
        assert mobius_subset(sub, n) == f

    # components DP
    for n in range(2, 8):
        for s in range(1, n + 1):
            for t in range(1, n + 1):
                if s == t: continue
                c = 0
                for p in permutations(range(1, n + 1)):
                    if p[0] == s and p[-1] == t and all(
                            (p[i - 1] < p[i]) != (p[i] < p[i + 1]) for i in range(1, n - 1)):
                        c += 1
                assert zigzag_perms(n, s, t) == c, (n, s, t)

    # BIT transitions
    for _ in range(200):
        a = [rnd.randint(0, 5) for _ in range(rnd.randint(0, 10))]
        dp = []
        for i in range(len(a)):
            dp.append(1 + sum(dp[j] for j in range(i) if a[j] < a[i]))
        assert count_increasing_subseq(a) == sum(dp)
        assert count_increasing_subseq(a, 7) == sum(dp) % 7

    # monotonic queue
    for _ in range(200):
        n = rnd.randint(1, 10); k = rnd.randint(1, n)
        c = [rnd.randint(0, 9) for _ in range(n)]
        best = INF
        for mask in range(1 << n):
            if all(any(mask >> j & 1 for j in range(i, i + k)) for i in range(n - k + 1)):
                best = min(best, sum(c[i] for i in range(n) if mask >> i & 1))
        assert min_cost_cover_windows(c, k) == best, (c, k)

    # CHT
    for _ in range(200):
        lines = sorted({(rnd.randint(-20, 20), rnd.randint(-50, 50)) for _ in range(rnd.randint(1, 15))},
                       key=lambda t: -t[0])
        h = CHT()
        for m, b in lines: h.add(m, b)
        xs = sorted(rnd.randint(-30, 30) for _ in range(30))
        for x in xs:
            want = min(m * x + b for m, b in lines)
            assert h.query(x) == want
        for x in xs:
            assert h.qmono(x) == min(m * x + b for m, b in lines)
    # interleaved add / qmono (typical DP usage)
    for _ in range(100):
        h = CHT(); lines = []
        m = 50
        for x in sorted(rnd.randint(-30, 30) for _ in range(20)):
            m -= rnd.randint(0, 3); b = rnd.randint(-50, 50)
            h.add(m, b); lines.append((m, b))
            assert h.qmono(x) == min(mm * x + bb for mm, bb in lines)

    # D&C opt + Aliens vs O(n^2 k)
    for _ in range(150):
        n = rnd.randint(1, 12); k = rnd.randint(1, n)
        a = [rnd.randint(1, 10) for _ in range(n)]
        P = [0]
        for x in a: P.append(P[-1] + x)
        dp = [0] + [INF] * n
        for _ in range(k):
            nd = [INF] * (n + 1)
            for i in range(1, n + 1):
                nd[i] = min(dp[j] + (P[i] - P[j]) ** 2 for j in range(i))
            dp = nd
        assert partition_k_squares(a, k) == dp[n] == aliens_partition(a, k), (a, k)

    # slope trick vs DP over candidate values
    for _ in range(200):
        a = [rnd.randint(0, 10) for _ in range(rnd.randint(0, 10))]
        vals = sorted(set(a)) or [0]
        dp = [0] * len(vals)
        for x in a:
            run, nd = INF, []
            for i, v in enumerate(vals):
                run = min(run, dp[i]); nd.append(run + abs(x - v))
            dp = nd
        assert min_cost_nondecreasing(a) == min(dp)

    print("ALL OK")
