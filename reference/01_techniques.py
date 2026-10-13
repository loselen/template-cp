# 01_techniques.py -- basic & advanced techniques, offline query tricks, search, heuristics.
# Python 3 stdlib only, PyPy-friendly.
import random, math, heapq
from bisect import bisect_right
from collections import deque
from itertools import combinations, permutations

# ==== 3.1 Brute force: subsets / permutations ====
# Subsets of n items: for m in range(1 << n): items = [i for i in range(n) if m >> i & 1]
# itertools: combinations(a, k), permutations(a[, k]), product(range(3), repeat=n) (base-3 enum).
# Recursion-free k-subsets in increasing mask order: Gosper's hack (see 3.8).


def next_permutation(a):
    # O(n) in place, lexicographic next; returns False if a was the last (and leaves it). Handles dups.
    i = len(a) - 2
    while i >= 0 and a[i] >= a[i + 1]:
        i -= 1
    if i < 0:
        return False
    j = len(a) - 1
    while a[j] <= a[i]:
        j -= 1
    a[i], a[j] = a[j], a[i]
    a[i + 1:] = a[:i:-1]
    return True


# ==== 3.1 Backtracking with pruning: N-queens count ====
def n_queens(n):
    # O(#solutions-ish), bitmask backtracking, depth n (fine recursively). n=12 -> 14200.
    full = (1 << n) - 1

    def go(cols, d1, d2):
        if cols == full:
            return 1
        cnt = 0
        free = full & ~(cols | d1 | d2)
        while free:
            b = free & -free
            free ^= b
            cnt += go(cols | b, ((d1 | b) << 1) & full, (d2 | b) >> 1)
        return cnt
    return go(0, 0, 0)


# ==== 3.1 / 12.3 Meet in the middle: count subsets with sum <= S ====
def subset_sums(a):
    # O(2^n): all subset sums (list doubling, no recursion)
    s = [0]
    for x in a:
        s += [v + x for v in s]
    return s


def mitm_count_leq(a, S):
    # O(2^(n/2) * n). n <= 40. For "== S" use Counter on one half.
    h = len(a) // 2
    L, R = subset_sums(a[:h]), sorted(subset_sums(a[h:]))
    return sum(bisect_right(R, S - x) for x in L)


# ==== 3.2 Grid problems: direction arrays ====
DIR4 = ((1, 0), (-1, 0), (0, 1), (0, -1))
DIR8 = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1))
KNIGHT = ((1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1))
# 'URDL' -> index: {c: i for i, c in enumerate('DURL')} matching DIR4 order (D=+row).
# Fast grid BFS: flatten id = r*W + c, pad grid with a wall border to skip bounds checks.


def grid_bfs(grid, sr, sc, wall="#"):
    # O(HW). Returns dist matrix (-1 unreachable).
    H, W = len(grid), len(grid[0])
    dist = [[-1] * W for _ in range(H)]
    dist[sr][sc] = 0
    q = deque([(sr, sc)])
    while q:
        r, c = q.popleft()
        for dr, dc in DIR4:
            nr, nc = r + dr, c + dc
            if 0 <= nr < H and 0 <= nc < W and dist[nr][nc] < 0 and grid[nr][nc] != wall:
                dist[nr][nc] = dist[r][c] + 1
                q.append((nr, nc))
    return dist


# ==== 3.2 Date/time and calendar ====
# Leap: y%4==0 and (y%100!=0 or y%400==0). Gregorian cycle = 400 years = 146097 days = 20871 weeks.
# datetime.date(y,m,d) works for 1..9999: .weekday() (Mon=0), (d2-d1).days, date.fromordinal.
def is_leap(y):
    return y % 4 == 0 and (y % 100 != 0 or y % 400 == 0)


def days_in_month(y, m):
    return 29 if m == 2 and is_leap(y) else (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[m - 1]


def day_of_week(y, m, d):
    # Sakamoto. 0 = Sunday ... 6 = Saturday. Gregorian.
    t = (0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4)
    if m < 3:
        y -= 1
    return (y + y // 4 - y // 100 + y // 400 + t[m - 1] + d) % 7


def days_from_civil(y, m, d):
    # Days since 1970-01-01 (any year, proleptic Gregorian). Diff of two = days between.
    y -= m <= 2
    era = y // 400
    yoe = y - era * 400
    doy = (153 * (m + (-3 if m > 2 else 9)) + 2) // 5 + d - 1
    doe = yoe * 365 + yoe // 4 - yoe // 100 + doy
    return era * 146097 + doe - 719468
# Time: convert hh:mm:ss to seconds, work mod 86400.


# ==== 3.2 Parsing expressions (shunting-yard, iterative) ====
def eval_expr(s):
    # O(len). Ints, + - * /, parentheses, unary +/-. '/' truncates toward zero (C-style).
    # Add '^' (right-assoc): pop only while prec[top] > prec[c].
    prec = {"+": 1, "-": 1, "*": 2, "/": 2, "u": 3}
    nums, ops = [], []

    def apply():
        op = ops.pop()
        if op == "u":
            nums.append(-nums.pop())
            return
        b = nums.pop()
        a = nums.pop()
        if op == "+":
            nums.append(a + b)
        elif op == "-":
            nums.append(a - b)
        elif op == "*":
            nums.append(a * b)
        else:
            q = abs(a) // abs(b)
            nums.append(q if (a >= 0) == (b >= 0) else -q)
    i, n, prev = 0, len(s), False   # prev: last token was an operand
    while i < n:
        c = s[i]
        if c.isdigit():
            j = i
            while j < n and s[j].isdigit():
                j += 1
            nums.append(int(s[i:j]))
            i, prev = j, True
            continue
        if c == "(":
            ops.append(c)
            prev = False
        elif c == ")":
            while ops[-1] != "(":
                apply()
            ops.pop()
            prev = True
        elif c in prec:
            if not prev:            # unary (prefix ops never pop)
                if c == "-":
                    ops.append("u")
            else:
                while ops and ops[-1] != "(" and prec[ops[-1]] >= prec[c]:
                    apply()
                ops.append(c)
                prev = False
        i += 1
    while ops:
        apply()
    return nums[-1]
# Edge cases checklist (3.2): empty input, n=1, inclusive/exclusive bounds, 0- vs 1-indexing,
# overflow irrelevant in Python but mod arithmetic still needed; re-read output format.


# ==== 3.3 Sorting: exchange arguments ====
# To order items so that cost is optimal, compare adjacent pair a,b: a before b iff
# cost(a,b) < cost(b,a). If it is a strict weak order, sort by key, e.g. concatenation
# "largest number": sort by cmp_to_key(lambda x, y: (y+x > x+y) - (y+x < x+y)).
# Weighted jobs min sum w*C: sort by t/w (key t_i*w_j < t_j*w_i -> use cmp_to_key or Fraction).


# ==== 3.3 Coordinate compression ====
def compress(a):
    # O(n log n). Returns (ranks, sorted unique values). Value of rank r: vals[r].
    vals = sorted(set(a))
    idx = {v: i for i, v in enumerate(vals)}
    return [idx[x] for x in a], vals
# Lookups of non-present values: bisect_left(vals, x).


# ==== 3.3 Counting sort / radix sort ====
def counting_sort(a, maxv):
    # O(n + maxv), values in [0, maxv].
    cnt = [0] * (maxv + 1)
    for x in a:
        cnt[x] += 1
    out = []
    for v in range(maxv + 1):
        if cnt[v]:
            out += [v] * cnt[v]
    return out


def radix_sort(a, bits=16):
    # O(n * w/bits), nonneg ints, stable LSD. Usually sorted() is faster in Python; useful in PyPy
    # for sorting pairs packed as ints: key = (x << 20) | i.
    if not a:
        return a
    mask = (1 << bits) - 1
    mx, sh = max(a), 0
    while (mx >> sh) or sh == 0:
        buckets = [[] for _ in range(mask + 1)]
        for x in a:
            buckets[(x >> sh) & mask].append(x)
        a = [x for b in buckets for x in b]
        sh += bits
    return a


# ==== 3.3 Binary search on sorted arrays ====
# bisect_left(a, x): first idx with a[i] >= x; bisect_right: first a[i] > x.
# count in [lo, hi]: bisect_right(a, hi) - bisect_left(a, lo). Py3.10+: key= param.


# ==== 3.3 Binary search on the answer (monotonic predicate) ====
def first_true(lo, hi, pred):
    # O(log(hi-lo)) calls. Smallest x in [lo, hi) with pred(x) true (pred F..FT..T); hi if none.
    while lo < hi:
        mid = (lo + hi) // 2
        if pred(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo
# Last true for T..TF..F: first_true(lo, hi, lambda x: not pred(x)) - 1.


def bsearch_real(lo, hi, pred, iters=100):
    # Fixed iterations (no eps loop -> no infinite loop). Smallest x with pred(x) true. 100 iters
    # halves 2^100 times: enough for any double range. ~60 is enough for 1e18-range to 1e-9 rel.
    for _ in range(iters):
        mid = (lo + hi) / 2
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return hi


# ==== 3.3 Ternary search on unimodal functions ====
def ternary_int_max(lo, hi, f):
    # argmax of f on [lo, hi], f strictly increasing then strictly decreasing (plateaus break it!).
    while hi - lo > 2:
        m1 = lo + (hi - lo) // 3
        m2 = hi - (hi - lo) // 3
        if f(m1) < f(m2):
            lo = m1 + 1
        else:
            hi = m2
    return max(range(lo, hi + 1), key=f)
# Alternative for ints: first_true(lo, hi, lambda x: f(x) >= f(x+1)) -- binary search on slope.


def ternary_real_min(lo, hi, f, iters=200):
    # min of unimodal (convex) f on [lo, hi]; returns x. Each iter keeps 2/3.
    for _ in range(iters):
        m1 = lo + (hi - lo) / 3
        m2 = hi - (hi - lo) / 3
        if f(m1) < f(m2):
            hi = m2
        else:
            lo = m1
    return (lo + hi) / 2


# ==== 3.3 / 11 Parallel binary search ====
def parallel_binary_search(Q, T, reset, apply, check):
    # For each query q: smallest t in [1, T] such that check(q) holds after apply(0..t-1); T+1 if never.
    # O((T + Q) log T) * cost(apply/check). reset() clears structure. check must be monotone in t.
    lo, hi = [1] * Q, [T + 1] * Q
    while True:
        act = [q for q in range(Q) if lo[q] < hi[q]]
        if not act:
            return lo
        act.sort(key=lambda q: (lo[q] + hi[q]) // 2)
        reset()
        t = 0
        for q in act:
            mid = (lo[q] + hi[q]) // 2
            while t < mid:
                apply(t)
                t += 1
            if check(q):
                hi[q] = mid
            else:
                lo[q] = mid + 1


# ==== 3.4 Prefix sums 1D / 2D ====
# 1D: P = [0]; for x in a: P.append(P[-1] + x)  -> sum a[l..r) = P[r] - P[l]. itertools.accumulate(a, initial=0).
# Prefix XOR: same with ^ (xor of [l,r) = X[r]^X[l]). Prefix products: only invertible (mod p, no zeros)
# -> count zeros separately; or use sparse table / segment tree.
def prefix2d(g):
    # O(HW). P[i][j] = sum g[0..i)[0..j). rect sum [r1,r2) x [c1,c2): see rect().
    H, W = len(g), len(g[0])
    P = [[0] * (W + 1) for _ in range(H + 1)]
    for i in range(H):
        row, pr, cur = g[i], P[i], P[i + 1]
        s = 0
        for j in range(W):
            s += row[j]
            cur[j + 1] = pr[j + 1] + s
    return P


def rect(P, r1, c1, r2, c2):
    return P[r2][c2] - P[r1][c2] - P[r2][c1] + P[r1][c1]


# ==== 3.4 Difference arrays (1D / 2D range add) ====
# 1D: D=[0]*(n+1); add v on [l,r): D[l]+=v; D[r]-=v; then a = prefix sums of D.
def diff2d_apply(H, W, updates):
    # updates: (r1, c1, r2, c2, v) half-open rectangles. O(HW + U).
    D = [[0] * (W + 1) for _ in range(H + 1)]
    for r1, c1, r2, c2, v in updates:
        D[r1][c1] += v
        D[r1][c2] -= v
        D[r2][c1] -= v
        D[r2][c2] += v
    for i in range(H + 1):
        for j in range(1, W + 1):
            D[i][j] += D[i][j - 1]
    for i in range(1, H + 1):
        for j in range(W + 1):
            D[i][j] += D[i - 1][j]
    return [row[:W] for row in D[:H]]


# ==== 3.4 Higher-order difference: add polynomial to ranges ====
def add_polys(n, updates, deg):
    # updates: (l, r, p): a[i] += p(i - l) for l <= i <= r; p: int -> int, degree <= deg.
    # O(n*K + U*K^2), K = deg+1. AP: deg=1, p = lambda x: a0 + d*x.
    K = deg + 1
    sgn = [math.comb(K, j) * (-1) ** j for j in range(K)]
    D = [0] * (n + 2 * K + 1)

    def put(start, f, s):
        vals = [f(x) for x in range(K)]
        for j in range(K):
            D[start + j] += s * sum(sgn[j - i] * vals[i] for i in range(j + 1))
    for l, r, p in updates:
        put(l, p, 1)
        sh = r + 1 - l
        put(r + 1, lambda x, p=p, sh=sh: p(x + sh), -1)
    for _ in range(K):           # K prefix-sum passes undo K-th difference
        s = 0
        for i in range(len(D)):
            s += D[i]
            D[i] = s
    return D[:n]


# ==== 3.5 Two pointers on sorted arrays ====
def count_pairs_leq(a, S):
    # a sorted. #pairs i<j with a[i]+a[j] <= S. O(n).
    i, j, cnt = 0, len(a) - 1, 0
    while i < j:
        if a[i] + a[j] <= S:
            cnt += j - i
            i += 1
        else:
            j -= 1
    return cnt


# ==== 3.5 Variable-size sliding window ====
def longest_at_most_k_distinct(a, k):
    # O(n). Longest subarray with <= k distinct values. Pattern: extend r, shrink l while invalid.
    cnt, l, best, distinct = {}, 0, 0, 0
    for r, x in enumerate(a):
        c = cnt.get(x, 0)
        if c == 0:
            distinct += 1
        cnt[x] = c + 1
        while distinct > k:
            y = a[l]
            cnt[y] -= 1
            if cnt[y] == 0:
                distinct -= 1
            l += 1
        best = max(best, r - l + 1)
    return best


# ==== 3.5 Sliding window min/max (monotonic deque) ====
def sliding_min(a, k):
    # O(n). min of each window a[i-k+1..i] for i >= k-1. For max flip comparison (<=  -> >=).
    dq, out = deque(), []
    for i, x in enumerate(a):
        while dq and a[dq[-1]] >= x:
            dq.pop()
        dq.append(i)
        if dq[0] <= i - k:
            dq.popleft()
        if i >= k - 1:
            out.append(a[dq[0]])
    return out


# ==== 3.5 Fast and slow pointers (Floyd cycle detection) ====
def floyd(f, x0):
    # O(mu + lam) time, O(1) memory. Returns (mu, lam): first index on cycle, cycle length.
    t, h = f(x0), f(f(x0))
    while t != h:
        t, h = f(t), f(f(h))
    mu, t = 0, x0
    while t != h:
        t, h = f(t), f(h)
        mu += 1
    lam, h = 1, f(t)
    while t != h:
        h = f(h)
        lam += 1
    return mu, lam


# ==== 3.6 Greedy: proving correctness / matroids ====
# Prove by exchange argument (swap an optimal solution's choice with greedy's, never worse) or
# "greedy stays ahead" induction. Greedy max-weight independent set is optimal iff the family is a
# matroid (hereditary + exchange): e.g. forests (Kruskal), linearly independent vectors (xor basis),
# scheduling unit jobs with deadlines. Matroid intersection (2 matroids) -> see 05_flow_matching.py.


# ==== 3.6 Interval scheduling (max non-overlapping) ====
def max_disjoint(iv):
    # half-open [l, r). Sort by right end. O(n log n).
    last, cnt = -float("inf"), 0
    for l, r in sorted(iv, key=lambda x: x[1]):
        if l >= last:
            cnt += 1
            last = r
    return cnt


# ==== 3.6 Interval covering (min intervals to cover [L, R]) ====
def min_cover(iv, L, R):
    # closed intervals; returns min count or -1. O(n log n).
    iv = sorted(iv)
    i, n, cur, cnt = 0, len(iv), L, 0
    while True:
        best = None
        while i < n and iv[i][0] <= cur:
            if best is None or iv[i][1] > best:
                best = iv[i][1]
            i += 1
        if best is None or best < cur:
            return -1
        cnt += 1
        if best >= R:
            return cnt
        if best == cur:          # can't extend (only zero-progress intervals)
            return -1
        cur = best


# ==== 3.6 Interval partitioning (min rooms, heap) ====
def min_rooms(iv):
    # half-open [l, r). Returns (k, room id per interval). O(n log n). k = max overlap.
    order = sorted(range(len(iv)), key=lambda i: iv[i][0])
    heap, free, room, k = [], [], [0] * len(iv), 0
    for i in order:
        l, r = iv[i]
        while heap and heap[0][0] <= l:
            heapq.heappush(free, heapq.heappop(heap)[1])
        if free:
            rid = heapq.heappop(free)
        else:
            rid, k = k, k + 1
        room[i] = rid
        heapq.heappush(heap, (r, rid))
    return k, room


# ==== 3.6 Job sequencing with deadlines (unit time) ====
def job_sequencing(jobs):
    # jobs: (deadline, profit); each takes 1 unit, must finish by deadline. Max profit. O(n log n).
    h = []
    for d, p in sorted(jobs):
        heapq.heappush(h, p)
        if len(h) > d:
            heapq.heappop(h)     # drop least profitable
    return sum(h)
# Variant with durations (max #jobs, deadlines): sort by deadline, keep max-heap of durations,
# pop the longest when total time exceeds deadline (Moore-Hodgson).


# ==== 3.6 Huffman coding (optimal merge cost) ====
def huffman_cost(w):
    # O(n log n). Sum of merge costs = total weighted code length.
    h = list(w)
    heapq.heapify(h)
    cost = 0
    while len(h) > 1:
        a = heapq.heappop(h)
        b = heapq.heappop(h)
        cost += a + b
        heapq.heappush(h, a + b)
    return cost
# k-ary: pad with zeros until (n-1) % (k-1) == 0, merge k at a time.


# ==== 3.6 Regret greedy (greedy with undo) ====
def stock_max_profit(prices):
    # Each day: buy 1, sell 1, or nothing (CF 865D). O(n log n). Push p twice: one is the
    # "undo" option (selling at p now can later be revoked and resold at higher q, gain q-p).
    h, ans = [], 0
    for p in prices:
        if h and h[0] < p:
            ans += p - heapq.heappop(h)
            heapq.heappush(h, p)
        heapq.heappush(h, p)
    return ans
# General pattern: when taking an item, push a compensating "regret" item (value -gain) so a later
# choice can undo it. Also: choose k non-adjacent max sum (linked list + heap, merge neighbors a+c-b).


# ==== 3.7 Merge sort and counting inversions ====
def count_inversions(a):
    # O(n log n), iterative bottom-up. #pairs i<j with a[i] > a[j]. (Alt: BIT over compressed values.)
    a = list(a)
    n, inv, w = len(a), 0, 1
    buf = [0] * n
    while w < n:
        for lo in range(0, n, 2 * w):
            mid, hi = min(lo + w, n), min(lo + 2 * w, n)
            i, j, k = lo, mid, lo
            while i < mid and j < hi:
                if a[j] < a[i]:
                    buf[k] = a[j]
                    inv += mid - i
                    j += 1
                else:
                    buf[k] = a[i]
                    i += 1
                k += 1
            buf[k:hi] = a[i:mid] + a[j:hi]
        a, buf = buf, a
        w *= 2
    return inv


# ==== 3.7 Quickselect ====
def quickselect(a, k, rng=random):
    # k-th smallest (0-based). Expected O(n). 3-way partition handles duplicates.
    a = list(a)
    while True:
        p = a[rng.randrange(len(a))]
        lt = [x for x in a if x < p]
        if k < len(lt):
            a = lt
            continue
        eq = sum(1 for x in a if x == p)
        if k < len(lt) + eq:
            return p
        k -= len(lt) + eq
        a = [x for x in a if x > p]
# In practice: sorted(a)[k] or heapq.nsmallest(k+1, a)[-1] is often faster in CPython.


# ==== 3.7 Closest pair ====
# see 10_geometry.py: closest pair


# ==== 3.7 Master theorem / Karatsuba ====
# T(n) = a T(n/b) + O(n^d): d > log_b a -> O(n^d); d == log_b a -> O(n^d log n); d < -> O(n^log_b a).
# Examples: merge sort 2T(n/2)+n = n log n; binary search T(n/2)+1 = log n; Karatsuba 3T(n/2)+n = n^1.585.
# Karatsuba: Python ints already use it for big multiplication -- just multiply ints. For polynomial
# multiplication pack coefficients into one big int (Kronecker) or use NTT: see 09_algebra_poly.py.


# ==== 3.8 Bit manipulation ====
# lowbit: x & -x; clear lowest: x & (x-1); test: x >> i & 1; set: x | 1<<i; clear: x & ~(1<<i);
# toggle: x ^ 1<<i; popcount: x.bit_count() (3.10+) or bin(x).count("1"); highest bit: x.bit_length()-1;
# trailing zeros: (x & -x).bit_length() - 1; is power of 2: x and not x & (x-1); all n bits: (1<<n)-1.
# Python ints are unbounded: ~x is negative (-x-1) -> mask with & ((1<<n)-1).
def submasks(m):
    # all submasks of m incl. 0, decreasing. Over all m: total O(3^n).
    s = m
    while True:
        yield s
        if s == 0:
            return
        s = (s - 1) & m


def supersets(m, n):
    # all supersets of m within n bits, increasing.
    s = m
    full = 1 << n
    while s < full:
        yield s
        s = (s + 1) | m


def k_subsets(n, k):
    # Gosper's hack: all n-bit masks with k bits set, increasing. k >= 1.
    s = (1 << k) - 1
    while s < 1 << n:
        yield s
        c = s & -s
        r = s + c
        s = (((r ^ s) >> 2) // c) | r


def gray(i):
    return i ^ (i >> 1)


def gray_inv(g):
    i = 0
    while g:
        i ^= g
        g >>= 1
    return i
# Gray: consecutive codes differ in exactly 1 bit (bit = trailing zeros of i+1) -> enumerate subsets
# updating state by one add/remove per step.
# XOR facts: a^a=0, xor is its own inverse (prefix xor), xor of 0..n = [n,1,n+1,0][n%4];
# a+b = (a^b) + 2(a&b); a^b <= a+b; min pairwise xor of a set is among adjacent elements after sorting;
# max xor pair / subset -> binary trie or linear basis (see 07_strings.py: BinaryTrie / 09_algebra_poly.py).
# SOS DP (sum over subsets) -> see 03_dp.py.


# ==== 11 Offline queries: overview ====
# Offline: read all queries, reorder (sort by r, by value, by time), answer, output in original order.
# Sweep line + BIT: sort queries by right end, add elements as r advances (e.g. distinct count:
#   keep only last occurrence of each value in BIT).
# Persistent structures (online version of the sweep) -> see 02_data_structures.py.
# Answering in reverse: deletions are hard, insertions easy -> process ops backwards
#   (deleting edges -> adding edges with DSU), then reverse the answers.
# Sqrt decomposition on queries: buffer up to ~sqrt(q) pending updates, answer queries with the
#   static structure + brute over the buffer, rebuild when buffer is full.
# Segment tree over time (offline dynamic connectivity, DSU with rollback) -> see 04_graphs.py.


# ==== 11 Mo's algorithm (Hilbert order) ====
def hilbert_order(x, y, k=21):
    # Position on Hilbert curve of order k (coords < 2^k). Sorting by this beats block order.
    d, s, n = 0, 1 << (k - 1), 1 << k
    while s:
        rx = 1 if x & s else 0
        ry = 1 if y & s else 0
        d += s * s * ((3 * rx) ^ ry)
        if ry == 0:
            if rx == 1:
                x, y = n - 1 - x, n - 1 - y
            x, y = y, x
        s >>= 1
    return d


def mo_distinct(a, queries):
    # O((n + q) sqrt n). Distinct count on inclusive [l, r]. Template: replace add/remove bodies.
    # Classic order alternative: key (l // B, r if (l//B) even else -r), B = n / sqrt(q).
    vals, _ = compress(a)
    cnt = [0] * (len(a) + 1)
    order = sorted(range(len(queries)), key=lambda i: hilbert_order(*queries[i]))
    ans = [0] * len(queries)
    L, R, cur = 0, -1, 0
    for qi in order:
        l, r = queries[qi]
        while R < r:
            R += 1
            v = vals[R]; cnt[v] += 1
            if cnt[v] == 1: cur += 1
        while L > l:
            L -= 1
            v = vals[L]; cnt[v] += 1
            if cnt[v] == 1: cur += 1
        while R > r:
            v = vals[R]; cnt[v] -= 1
            if cnt[v] == 0: cur -= 1
            R -= 1
        while L < l:
            v = vals[L]; cnt[v] -= 1
            if cnt[v] == 0: cur -= 1
            L += 1
        ans[qi] = cur
    return ans


# ==== 11 Mo's algorithm with updates ====
def mo_updates_distinct(a, ops):
    # ops: ('Q', l, r) inclusive or ('U', pos, val). O(n^(5/3)). Block B ~ n^(2/3).
    # Sort queries by (l//B, r//B, t); move time pointer applying/undoing updates by swapping values.
    a = list(a)
    allv = list(a) + [o[2] for o in ops if o[0] == "U"]
    _, vs = compress(allv)
    idx = {v: i for i, v in enumerate(vs)}
    a = [idx[x] for x in a]
    ups, qs = [], []
    for o in ops:
        if o[0] == "U":
            ups.append([o[1], idx[o[2]]])
        else:
            qs.append((o[1], o[2], len(ups), len(qs)))
    n = len(a)
    B = max(1, round(n ** (2 / 3)))
    qs.sort(key=lambda q: (q[0] // B, q[1] // B, q[2]))
    cnt = [0] * len(vs)
    ans = [0] * len(qs)
    L, R, T, cur = 0, -1, 0, 0
    for l, r, t, qi in qs:
        while T < t or T > t:
            if T < t:
                p = ups[T]; T += 1
            else:
                T -= 1; p = ups[T]
            pos, v = p
            if L <= pos <= R:
                cnt[a[pos]] -= 1
                if cnt[a[pos]] == 0: cur -= 1
                cnt[v] += 1
                if cnt[v] == 1: cur += 1
            a[pos], p[1] = v, a[pos]      # swap: applying again undoes it
        while R < r:
            R += 1; cnt[a[R]] += 1
            if cnt[a[R]] == 1: cur += 1
        while L > l:
            L -= 1; cnt[a[L]] += 1
            if cnt[a[L]] == 1: cur += 1
        while R > r:
            cnt[a[R]] -= 1
            if cnt[a[R]] == 0: cur -= 1
            R -= 1
        while L < l:
            cnt[a[L]] -= 1
            if cnt[a[L]] == 0: cur -= 1
            L += 1
        ans[qi] = cur
    return ans


# ==== 11 CDQ divide and conquer (3D partial order) ====
def cdq_3d(pts):
    # For each i: #j != i with a_j<=a_i, b_j<=b_i, c_j<=c_i. O(n log^2 n). Recursion depth log n.
    # Idea: split by a; left half contributes to right half via sweep on b + BIT on c.
    uniq = sorted(set(pts))
    from collections import Counter
    mult = Counter(pts)
    m = len(uniq)
    w = [mult[p] for p in uniq]
    _, cs = compress([p[2] for p in uniq])
    cidx = {v: i + 1 for i, v in enumerate(cs)}
    C = [cidx[p[2]] for p in uniq]
    Bv = [p[1] for p in uniq]
    res = [0] * m
    bit = [0] * (len(cs) + 1)
    K = len(cs)

    def solve(lo, hi):
        if hi - lo <= 1:
            return
        mid = (lo + hi) // 2
        solve(lo, mid)
        solve(mid, hi)
        left = sorted(range(lo, mid), key=lambda i: Bv[i])
        right = sorted(range(mid, hi), key=lambda i: Bv[i])
        p = 0
        for i in right:
            while p < len(left) and Bv[left[p]] <= Bv[i]:
                j = C[left[p]]
                while j <= K:
                    bit[j] += w[left[p]]; j += j & -j
                p += 1
            j, s = C[i], 0
            while j:
                s += bit[j]; j -= j & -j
            res[i] += s
        for q in range(p):        # undo BIT
            j = C[left[q]]
            while j <= K:
                bit[j] -= w[left[q]]; j += j & -j
    solve(0, m)
    pos = {p: i for i, p in enumerate(uniq)}
    return [res[pos[p]] + w[pos[p]] - 1 for p in pts]
# Uniques sorted by (a,b,c) => any j dominating-before i with equal a lies left of i. Duplicates merged.
# CDQ also: offline dynamic problems where updates in left half affect queries in right half.


# ==== 12.1 Sqrt techniques ====
# Sqrt decomposition on arrays: blocks of size B~sqrt(n), keep per-block aggregate (+ lazy tag);
#   query = partial blocks brute + full blocks aggregate. (Block structures: see 02_data_structures.py.)
# Heavy/light by degree: vertices with deg > sqrt(m) are heavy (<= 2 sqrt(m) of them); handle heavy
#   specially (precompute / push updates to them), light ones by brute over neighbors.
#   Triangle counting: orient edges low->high (deg, id), O(m sqrt m).
# By frequency: values occurring > sqrt(n) times are few -> precompute per heavy value.
# Python: sqrt algorithms are often too slow in CPython; prefer PyPy and flat lists.


# ==== 12.2 Randomization ====
# Las Vegas (always correct, random time: quicksort/quickselect) vs Monte Carlo (fixed time, may err:
#   Miller-Rabin, hashing, Freivalds' matrix check A*B==C in O(n^2) per trial).
# Shuffle input (random.shuffle) to defeat anti-quicksort / anti-incremental tests (e.g. randomized
#   incremental min enclosing circle -> see 10_geometry.py).
# Hashing with random base/mod (string hashing -> 07_strings.py); dict keys salt -> 00_core.py.
# Random sampling: element occurring > n/2 times: pick random, check (P(fail) 2^-k after k tries).
# Sum/xor hashing of sets: assign random 64-bit value per element; multiset equality = sum/xor equality.
# Probabilistic counting: estimate #distinct by min hash.
# Color coding: find path of k distinct vertices: color randomly with k colors, DP over (v, colorset)
#   O(2^k m); success prob k!/k^k ~ e^-k -> repeat ~e^k times.


# ==== 12.2 / 13.4 Simulated annealing / hill climbing ====
def simulated_annealing(init, cost, neighbor, iters, t0, t1, rng):
    # Minimize cost. neighbor(state, rng) -> new state (don't mutate input). Geometric cooling t0 -> t1.
    # Contest tips: time-based schedule (time.perf_counter), O(1) delta evaluation instead of full cost,
    # t0 ~ typical |delta| at start, t1 ~ smallest meaningful delta.
    cur, cc = init, cost(init)
    best, bc = cur, cc
    for i in range(iters):
        T = t0 * (t1 / t0) ** (i / iters)
        nxt = neighbor(cur, rng)
        nc = cost(nxt)
        if nc <= cc or rng.random() < math.exp((cc - nc) / T):
            cur, cc = nxt, nc
            if cc < bc:
                best, bc = cur, cc
    return best, bc


def hill_climb(init, cost, neighbor, iters, rng, restarts=1, rand_init=None):
    # Accept only non-worsening moves; random restarts via rand_init(rng).
    best, bc = None, float("inf")
    for r in range(restarts):
        cur = init if r == 0 or rand_init is None else rand_init(rng)
        cc = cost(cur)
        for _ in range(iters):
            nxt = neighbor(cur, rng)
            nc = cost(nxt)
            if nc <= cc:
                cur, cc = nxt, nc
        if cc < bc:
            best, bc = cur, cc
    return best, bc
# Beam search: keep the best W partial states per depth, expand each, keep top W (heapq.nlargest).
# Genetic algorithms: rarely worth it in contests; SA usually beats them.


# ==== 12.3 Iterative deepening (IDDFS / IDA*) ====
def ida_star(start, is_goal, nbrs, h, limit=10 ** 9):
    # Returns optimal cost or -1. Memory O(depth). h admissible. nbrs(u) -> iterable of (v, w).
    # Recursion depth = solution length (small by nature). Order nbrs by h for speed.
    path, found = {start}, [None]

    def dfs(u, g, bound):
        f = g + h(u)
        if f > bound:
            return f
        if is_goal(u):
            found[0] = g
            return -1
        mn = float("inf")
        for v, w in nbrs(u):
            if v in path:
                continue
            path.add(v)
            t = dfs(v, g + w, bound)
            path.discard(v)
            if t == -1:
                return -1
            mn = min(mn, t)
        return mn
    bound = h(start)
    while bound <= limit:
        t = dfs(start, 0, bound)
        if t == -1:
            return found[0]
        if t == float("inf"):
            return -1
        bound = t
    return -1


# IDDFS = IDA* with h = 0: DFS with depth limit 0,1,2,... (BFS-optimal, DFS memory).


# ==== 12.3 A* search ====
def astar(start, goal, nbrs, h):
    # Dijkstra with priority g + h(u). h consistent (h(u) <= w + h(v)) => first pop of goal optimal.
    # O(E log V) worst. Returns dist or -1. nbrs(u) -> (v, w).
    dist = {start: 0}
    pq = [(h(start), 0, start)]
    while pq:
        f, g, u = heapq.heappop(pq)
        if g > dist[u]:
            continue
        if u == goal:
            return g
        for v, w in nbrs(u):
            ng = g + w
            if ng < dist.get(v, float("inf")):
                dist[v] = ng
                heapq.heappush(pq, (ng + h(v), ng, v))
    return -1
# States must be hashable & comparable (tuples) for the heap tie-breaks; else add a counter.


# ==== 12.3 Bidirectional BFS ====
def bidir_bfs(s, t, nbrs):
    # Unweighted, undirected (or give reverse nbrs for t side). O(b^(d/2)). Returns dist or -1.
    if s == t:
        return 0
    ds, dt = {s: 0}, {t: 0}
    fs, ft = [s], [t]
    while fs and ft:
        if len(fs) > len(ft):
            fs, ft, ds, dt = ft, fs, dt, ds
        best = None
        nf = []
        for u in fs:
            du = ds[u] + 1
            for v in nbrs(u):
                if v in dt:
                    c = du + dt[v]
                    if best is None or c < best:
                        best = c
                if v not in ds:
                    ds[v] = du
                    nf.append(v)
        if best is not None:
            return best
        fs = nf
    return -1
# Expand a whole level of the smaller frontier, take min over all meetings in that level.


# ==== 12.3 Branch and bound / pruning ====
# Keep global best; prune a branch when (current + optimistic bound of remaining) can't beat best.
# Good bounds: fractional relaxation (knapsack: greedy by ratio), sum of remaining positives.
# Order branches best-first (greedy choice first) so a good best is found early.
# Pruning strategies: symmetry breaking (fix first element), feasibility checks, memo of dead states,
# sort items descending so infeasibility shows early.


# ==== 12.4 Binary lifting / doubling on functional graph ====
class FuncGraph:
    # nxt[v] = successor. Build O(n log K), kth(v, k) O(log K). Also carry sums/min in parallel tables.
    def __init__(self, nxt, maxk):
        self.up = [list(nxt)]
        for _ in range(max(1, maxk.bit_length()) - 1):
            p = self.up[-1]
            self.up.append([p[p[v]] for v in range(len(nxt))])

    def kth(self, v, k):
        j = 0
        while k:
            if k & 1:
                v = self.up[j][v]
            k >>= 1
            j += 1
        return v
# Huge k (1e18) with one start: walk until cycle found (dict of first visit), then k mod cycle length.
# Small-to-large merging, binary lifting LCA -> see 06_trees.py.


# ==== 12.4 Other advanced techniques (notes) ====
# Bitset speedup: Python big ints as bitsets: subset sum reach |= reach << w (mask to limit);
#   adjacency rows as ints, (row_u & row_v).bit_count() for common neighbors. Huge constant win.
# Four Russians: split into blocks of size ~log n, precompute all 2^blk block answers/transitions
#   (e.g. boolean matrix mult n^3/log n, LCS bit-parallel, +-1 RMQ). In Python use the bitset trick.
# FPT (conceptual): exponential only in a parameter k, e.g. vertex cover O(2^k * m) by branching on
#   an uncovered edge's endpoints; kernelization shrinks the instance first.
# Matroid intersection -> 05_flow_matching.py.
# Offline-to-online: if queries are forced online (xor with last answer), use persistent DS or
#   sqrt rebuilding; binary grouping (keep log n static structures of sizes 2^i, merge like a counter).


# ==== 13 Problem types and paradigms (notes) ====
# 13.1 Constructive: brute force small n and print answers to spot the pattern; try simple families
#   (sorted, alternating, 1..n shifted, identity + swaps, star/path graphs); prove impossibility via
#   invariant (parity, sum mod k, count bound). Invariant = never changes; monovariant = strictly
#   monotone -> termination.
# 13.2 Interactive: bound queries by information theory (need >= log2(#answers) bits; ternary replies
#   -> log3). Binary search via queries; adaptive judges may change the hidden answer consistently,
#   so strategies must be worst-case. Flush after each query (00_core.py: ask/answer).
# 13.3 Ad hoc: look for invariants, parity, symmetry; reformulate (reduce to known problem); count by
#   contribution (sum over elements of #subarrays/pairs where it matters) instead of over structures.
# 13.4 Output-only / heuristic: run long offline searches; scoring-based -> SA / beam search / local
#   search; always keep a valid greedy baseline answer first.
# 13.5 Communication (two-run): encode info in output of run 1 for run 2. Grid puzzles / games:
#   state search (BFS/IDA*) or game theory (09). Subtasks: grab partial points with brute force first.


# ============================== tests ==============================
if __name__ == "__main__":
    rng = random.Random(12345)

    # next_permutation
    for n in range(1, 6):
        a = sorted(rng.randrange(3) for _ in range(n))
        seen = [tuple(a)]
        while next_permutation(a):
            seen.append(tuple(a))
        assert seen == sorted(set(permutations(a))), (a, seen)

    # n queens
    assert [n_queens(n) for n in range(1, 10)] == [1, 0, 0, 2, 10, 4, 40, 92, 352]

    # MITM
    for _ in range(200):
        n = rng.randint(0, 12)
        a = [rng.randint(-10, 20) for _ in range(n)]
        S = rng.randint(-10, 60)
        br = sum(1 for m in range(1 << n) if sum(a[i] for i in range(n) if m >> i & 1) <= S)
        assert mitm_count_leq(a, S) == br

    # grid bfs
    g = ["..#.", "...#", "#...", "...."]
    d = grid_bfs(g, 0, 0)
    assert d[3][3] == 6 and d[0][2] == -1 and d[0][3] == -1 and d[2][3] == 5

    # dates
    import datetime
    base = datetime.date(1970, 1, 1)
    for _ in range(2000):
        y, m = rng.randint(1, 9999), rng.randint(1, 12)
        dd = rng.randint(1, days_in_month(y, m))
        dt = datetime.date(y, m, dd)
        assert day_of_week(y, m, dd) == (dt.weekday() + 1) % 7
        assert days_from_civil(y, m, dd) == (dt - base).days
        assert is_leap(y) == ((datetime.date(y, 3, 1) - datetime.date(y, 2, 1)).days == 29)

    # expression parsing vs eval
    def gen_expr(depth):
        r = rng.random()
        if depth == 0 or r < 0.3:
            return str(rng.randint(0, 20))
        if r < 0.4:
            return "-" + gen_expr(depth - 1)
        if r < 0.55:
            return "(" + gen_expr(depth - 1) + ")"
        return gen_expr(depth - 1) + rng.choice(" + - * ".split()) + gen_expr(depth - 1)
    for _ in range(2000):
        e = gen_expr(5)
        assert eval_expr(e) == eval(e), e
    assert eval_expr("7/2") == 3 and eval_expr("-7/2") == -3 and eval_expr("7/-2") == -3
    assert eval_expr("2*(3+4)*5-10/3") == 67 and eval_expr(" 1 - - 1 ") == 2

    # compression / sorts
    for _ in range(200):
        a = [rng.randint(0, 1000) for _ in range(rng.randint(0, 50))]
        r, vals = compress(a)
        assert [vals[x] for x in r] == a and vals == sorted(set(a))
        assert counting_sort(a, 1000) == sorted(a)
        b = [rng.getrandbits(rng.randint(1, 60)) for _ in range(rng.randint(0, 50))]
        assert radix_sort(b, 8) == sorted(b)

    # binary search on answer
    for _ in range(500):
        lo, hi = rng.randint(-50, 50), rng.randint(-50, 60)
        if lo > hi:
            lo, hi = hi, lo
        th = rng.randint(lo - 5, hi + 5)
        exp = min([x for x in range(lo, hi) if x >= th], default=hi)
        assert first_true(lo, hi, lambda x: x >= th) == exp
    assert abs(bsearch_real(0, 10, lambda x: x * x >= 2) - math.sqrt(2)) < 1e-12
    for _ in range(300):
        lo, hi = rng.randint(-30, 0), rng.randint(0, 30)
        pk = rng.randint(lo, hi)
        f = lambda x: -abs(x - pk) * (2 if x < pk else 3)
        assert ternary_int_max(lo, hi, f) == pk
    assert abs(ternary_real_min(-10, 10, lambda x: (x - 1.234) ** 2 + 5) - 1.234) < 1e-7

    # parallel binary search: range-add updates, query earliest time a[p] >= need
    for _ in range(50):
        n, T, Q = rng.randint(1, 10), rng.randint(1, 15), rng.randint(1, 15)
        U = [(rng.randrange(n), rng.randrange(n), rng.randint(1, 5)) for _ in range(T)]
        U = [(min(l, r), max(l, r), v) for l, r, v in U]
        qs = [(rng.randrange(n), rng.randint(1, 30)) for _ in range(Q)]
        bit = [0] * (n + 2)

        def badd(i, v):
            i += 1
            while i <= n + 1:
                bit[i] += v; i += i & -i

        def reset():
            for i in range(n + 2):
                bit[i] = 0

        def apply(t):
            l, r, v = U[t]
            badd(l, v); badd(r + 1, -v)

        def check(q):
            i, s = qs[q][0] + 1, 0
            while i:
                s += bit[i]; i -= i & -i
            return s >= qs[q][1]
        got = parallel_binary_search(Q, T, reset, apply, check)
        for q, (p, need) in enumerate(qs):
            s, ans = 0, T + 1
            for t in range(T):
                l, r, v = U[t]
                if l <= p <= r:
                    s += v
                if s >= need:
                    ans = t + 1
                    break
            assert got[q] == ans

    # prefix 2D / diff 2D
    for _ in range(100):
        H, W = rng.randint(1, 6), rng.randint(1, 6)
        gr = [[rng.randint(-5, 5) for _ in range(W)] for _ in range(H)]
        P = prefix2d(gr)
        r1, r2 = sorted(rng.randint(0, H) for _ in range(2))
        c1, c2 = sorted(rng.randint(0, W) for _ in range(2))
        assert rect(P, r1, c1, r2, c2) == sum(gr[i][j] for i in range(r1, r2) for j in range(c1, c2))
        ups, br = [], [[0] * W for _ in range(H)]
        for _ in range(rng.randint(0, 5)):
            r1, r2 = sorted(rng.randint(0, H) for _ in range(2))
            c1, c2 = sorted(rng.randint(0, W) for _ in range(2))
            v = rng.randint(-3, 3)
            ups.append((r1, c1, r2, c2, v))
            for i in range(r1, r2):
                for j in range(c1, c2):
                    br[i][j] += v
        assert diff2d_apply(H, W, ups) == br

    # polynomial range add
    for _ in range(200):
        n, deg = rng.randint(1, 15), rng.randint(0, 3)
        br, ups = [0] * n, []
        for _ in range(rng.randint(0, 5)):
            l, r = sorted(rng.randrange(n) for _ in range(2))
            co = [rng.randint(-3, 3) for _ in range(deg + 1)]
            p = lambda x, co=co: sum(c * x ** e for e, c in enumerate(co))
            ups.append((l, r, p))
            for i in range(l, r + 1):
                br[i] += p(i - l)
        assert add_polys(n, ups, deg) == br

    # two pointers / window / deque
    for _ in range(300):
        a = sorted(rng.randint(-10, 10) for _ in range(rng.randint(0, 12)))
        S = rng.randint(-15, 15)
        assert count_pairs_leq(a, S) == sum(1 for i, j in combinations(range(len(a)), 2) if a[i] + a[j] <= S)
        b = [rng.randint(0, 4) for _ in range(rng.randint(0, 12))]
        k = rng.randint(1, 4)
        assert longest_at_most_k_distinct(b, k) == max(
            [0] + [j - i for i in range(len(b)) for j in range(i + 1, len(b) + 1) if len(set(b[i:j])) <= k])
        if b:
            w = rng.randint(1, len(b))
            assert sliding_min(b, w) == [min(b[i - w + 1:i + 1]) for i in range(w - 1, len(b))]

    # floyd
    for _ in range(200):
        n = rng.randint(1, 20)
        nx = [rng.randrange(n) for _ in range(n)]
        x0 = rng.randrange(n)
        seq, pos, x = [], {}, x0
        while x not in pos:
            pos[x] = len(seq); seq.append(x); x = nx[x]
        assert floyd(lambda v: nx[v], x0) == (pos[x], len(seq) - pos[x])

    # greedy: intervals
    for _ in range(300):
        n = rng.randint(0, 7)
        iv = []
        for _ in range(n):
            l = rng.randint(0, 10)
            iv.append((l, l + rng.randint(1, 5)))
        best = 0
        for m in range(1 << n):
            sel = sorted(iv[i] for i in range(n) if m >> i & 1)
            if all(sel[i][1] <= sel[i + 1][0] for i in range(len(sel) - 1)):
                best = max(best, len(sel))
        assert max_disjoint(iv) == best
        k, room = min_rooms(iv)
        ov = max([0] + [sum(1 for l, r in iv if l <= t < r) for t in range(20)])
        assert k == ov
        for i in range(n):
            for j in range(i):
                if room[i] == room[j]:
                    assert iv[i][1] <= iv[j][0] or iv[j][1] <= iv[i][0]
        # covering (closed intervals, integer endpoints, coverage of real segment)
        civ = [(l, r) for l, r in iv] + [(x, x) for x in [rng.randint(0, 10)] if rng.random() < 0.3]
        L, R = sorted(rng.randint(0, 12) for _ in range(2))

        def covers(sel):
            sel = sorted(sel)
            cur, ok = L, False
            for l, r in sel:
                if l <= cur <= r:
                    cur, ok = r, True
            return ok and cur >= R
        m_ = len(civ)
        bestc = min([bin(m).count("1") for m in range(1 << m_) if covers([civ[i] for i in range(m_) if m >> i & 1])],
                    default=-1)
        assert min_cover(civ, L, R) == bestc, (civ, L, R)

    # job sequencing, huffman, regret greedy
    for _ in range(300):
        n = rng.randint(0, 7)
        jobs = [(rng.randint(1, 4), rng.randint(1, 20)) for _ in range(n)]
        best = 0
        for m in range(1 << n):
            sel = sorted(jobs[i] for i in range(n) if m >> i & 1)
            if all(sel[i][0] >= i + 1 for i in range(len(sel))):
                best = max(best, sum(p for _, p in sel))
        assert job_sequencing(jobs) == best

        w = [rng.randint(1, 10) for _ in range(rng.randint(1, 6))]

        def hbrute(ws):
            if len(ws) == 1:
                return 0
            res = float("inf")
            for i, j in combinations(range(len(ws)), 2):
                rest = [ws[k] for k in range(len(ws)) if k != i and k != j] + [ws[i] + ws[j]]
                res = min(res, ws[i] + ws[j] + hbrute(rest))
            return res
        assert huffman_cost(w) == hbrute(w)

        pr = [rng.randint(1, 10) for _ in range(rng.randint(0, 8))]
        dp = {0: 0}                     # holdings -> best cash
        for p in pr:
            nd = {}
            for hld, c in dp.items():
                for nh, nc in ((hld, c), (hld + 1, c - p), (hld - 1, c + p)):
                    if nh >= 0 and nd.get(nh, -10 ** 9) < nc:
                        nd[nh] = nc
            dp = nd
        assert stock_max_profit(pr) == dp[0]

    # inversions, quickselect
    for _ in range(300):
        a = [rng.randint(0, 6) for _ in range(rng.randint(0, 20))]
        assert count_inversions(a) == sum(1 for i, j in combinations(range(len(a)), 2) if a[i] > a[j])
        if a:
            k = rng.randrange(len(a))
            assert quickselect(a, k, rng) == sorted(a)[k]

    # bits
    for _ in range(200):
        n = rng.randint(1, 8)
        m = rng.randrange(1 << n)
        assert list(submasks(m)) == [s for s in range((1 << n) - 1, -1, -1) if s & m == s]
        assert list(supersets(m, n)) == [s for s in range(1 << n) if s & m == m]
        k = rng.randint(1, n)
        assert list(k_subsets(n, k)) == [s for s in range(1 << n) if bin(s).count("1") == k]
    for i in range(1 << 10):
        assert gray_inv(gray(i)) == i
        assert bin(gray(i) ^ gray(i + 1)).count("1") == 1

    # Mo
    for _ in range(100):
        n = rng.randint(1, 30)
        a = [rng.randint(0, 8) for _ in range(n)]
        qs = [tuple(sorted((rng.randrange(n), rng.randrange(n)))) for _ in range(rng.randint(1, 30))]
        assert mo_distinct(a, qs) == [len(set(a[l:r + 1])) for l, r in qs]
        ops = []
        for _ in range(rng.randint(1, 30)):
            if rng.random() < 0.4:
                ops.append(("U", rng.randrange(n), rng.randint(0, 8)))
            else:
                l, r = sorted((rng.randrange(n), rng.randrange(n)))
                ops.append(("Q", l, r))
        b, br = list(a), []
        for o in ops:
            if o[0] == "U":
                b[o[1]] = o[2]
            else:
                br.append(len(set(b[o[1]:o[2] + 1])))
        assert mo_updates_distinct(a, ops) == br

    # CDQ
    for _ in range(100):
        n = rng.randint(1, 40)
        pts = [(rng.randint(0, 4), rng.randint(0, 4), rng.randint(0, 4)) for _ in range(n)]
        br = [sum(1 for j in range(n) if j != i and all(pts[j][t] <= pts[i][t] for t in range(3)))
              for i in range(n)]
        assert cdq_3d(pts) == br

    # SA / hill climbing on small TSP
    n = 7
    P = [(rng.random(), rng.random()) for _ in range(n)]
    D = [[math.dist(p, q) for q in P] for p in P]
    tour_cost = lambda t: sum(D[t[i]][t[i - 1]] for i in range(n))
    opt = min(tour_cost((0,) + p) for p in permutations(range(1, n)))

    def two_opt(t, r):
        i, j = sorted(r.sample(range(n), 2))
        return t[:i] + t[i:j + 1][::-1] + t[j + 1:]
    _, c = simulated_annealing(tuple(range(n)), tour_cost, two_opt, 20000, 1.0, 1e-3, random.Random(7))
    assert abs(c - opt) < 1e-9
    _, c2 = hill_climb(tuple(range(n)), tour_cost, two_opt, 2000, random.Random(3), restarts=10,
                       rand_init=lambda r: tuple(r.sample(range(n), n)))
    assert abs(c2 - opt) < 1e-9

    # 8-puzzle: BFS from goal gives exact distances; compare A*, IDA*, bidir BFS
    goal = (1, 2, 3, 4, 5, 6, 7, 8, 0)

    def pnbrs(s):
        z = s.index(0)
        r, c = divmod(z, 3)
        for dr, dc in DIR4:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 3 and 0 <= nc < 3:
                y = nr * 3 + nc
                t = list(s)
                t[z], t[y] = t[y], 0
                yield tuple(t)
    GP = {v: divmod((v - 1) % 9, 3) for v in range(1, 9)}

    def manh(s):
        h = 0
        for i, v in enumerate(s):
            if v:
                r, c = GP[v]
                h += abs(r - i // 3) + abs(c - i % 3)
        return h
    for _ in range(15):
        s = goal
        for _ in range(rng.randint(0, 40)):
            s = rng.choice(list(pnbrs(s)))
        # exact via plain BFS
        dist = {s: 0}
        q = deque([s])
        while goal not in dist:
            u = q.popleft()
            for v in pnbrs(u):
                if v not in dist:
                    dist[v] = dist[u] + 1
                    q.append(v)
        exact = dist[goal]
        wn = lambda u: ((v, 1) for v in pnbrs(u))
        assert astar(s, goal, wn, manh) == exact
        assert ida_star(s, lambda u: u == goal, wn, manh) == exact
        assert bidir_bfs(s, goal, pnbrs) == exact
    unsolv = (2, 1, 3, 4, 5, 6, 7, 8, 0)
    assert bidir_bfs(unsolv, goal, pnbrs) == -1

    # A* weighted grid vs Dijkstra (h = 0)
    for _ in range(20):
        H, W = rng.randint(1, 8), rng.randint(1, 8)
        cost = [[rng.randint(1, 9) for _ in range(W)] for _ in range(H)]

        def gn(u):
            r, c = u
            for dr, dc in DIR4:
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W:
                    yield (nr, nc), cost[nr][nc]
        tgt = (H - 1, W - 1)
        assert astar((0, 0), tgt, gn, lambda u: abs(u[0] - tgt[0]) + abs(u[1] - tgt[1])) == \
            astar((0, 0), tgt, gn, lambda u: 0)

    # functional graph k-th
    for _ in range(100):
        n = rng.randint(1, 30)
        nx = [rng.randrange(n) for _ in range(n)]
        fg = FuncGraph(nx, 1000)
        v, k = rng.randrange(n), rng.randint(0, 1000)
        x = v
        for _ in range(k):
            x = nx[x]
        assert fg.kth(v, k) == x

    print("ALL OK")
