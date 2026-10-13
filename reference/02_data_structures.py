# 02_data_structures.py -- Section 4 data structures + Python replacements for STL (2.2/2.3).
# Python 3 stdlib only, PyPy-friendly. Indices 0-based, ranges half-open [l, r) unless noted.
import math
import random
from bisect import bisect_left, bisect_right, insort
from heapq import heappush, heappop
from collections import deque, Counter

# ==== 2.2 STL containers -> Python ====
# vector/array -> list (append/pop O(1); pop(0)/insert(0) O(n)!). 2D: [[0]*m for _ in range(n)].
# deque/queue/stack -> collections.deque (append/appendleft/pop/popleft O(1); d[i] is O(n)).
# priority_queue -> heapq (MIN-heap on a list): heappush, heappop, h[0]=top, heapify O(n).
#   max-heap: push -x; key-ordered: push (key, id, payload); tie id avoids comparing payloads.
#   heappushpop(h,x) / heapreplace(h,x) are faster combined ops; nlargest/nsmallest(k, it).
# set/multiset/map + pbds ordered set -> SortedList below (no stdlib balanced BST).
#   SortedDict-ish: SortedList of keys + plain dict key->value; iterate keys in order.
# unordered_set/map -> set/dict (int keys are hackable on CF: salt keys, see 00_core.py 4.3/19);
#   Counter(it) multiset counts, defaultdict(list) adjacency, dict keeps insertion order.
# bitset -> Python int (see 4.8 bitset note). pair/tuple -> tuple (compare lexicographic).
# string -> str (immutable: build with list + ''.join); bytes for fast byte ops.

# ==== 2.3 STL algorithms -> Python ====
# sort/stable_sort -> list.sort()/sorted(): Timsort, ALWAYS stable; key=... (fast) or
#   key=functools.cmp_to_key(cmp) with cmp(a,b)-> neg/0/pos (slow, needs strict weak order).
#   multi-key: sort by secondary first, then primary (stability), or key=lambda t:(t[0],-t[1]).
# lower_bound/upper_bound -> bisect_left/bisect_right(a, x, lo, hi, key=...) (key= py3.10+).
#   equal_range = (bisect_left, bisect_right); binary_search = i<len(a) and a[i]==x.
# next_permutation -> itertools.permutations(sorted(a)) (lexicographic, duplicates repeated);
#   combinations, combinations_with_replacement, product(range(2), repeat=n) for masks.
# unique (sorted) -> sorted(set(a)); keep order: list(dict.fromkeys(a)). reverse -> a[::-1].
# rotate left k -> a[k:]+a[:k]. nth_element/partial_sort -> heapq.nsmallest(k,a) or sort.
# accumulate/partial_sum -> itertools.accumulate(a, initial=0) (prefix sums), op= any func.
# iota -> list(range(n)); max_element -> max(range(n), key=a.__getitem__); min/max/sum builtins.


# SortedList: sqrt-bucket sorted multiset. add/discard O(sqrt n) (memmove, fast), index/kth O(n^(1/2)).
# Replaces std::multiset and pbds order-statistics tree. For a set, check `x in s` before add.
class SortedList:
    BUCKET_RATIO, SPLIT_RATIO = 16, 24

    def __init__(self, a=()):
        a = sorted(a)
        n = self.size = len(a)
        k = int(math.ceil(math.sqrt(n / self.BUCKET_RATIO))) or 1
        self.a = [a[n * i // k: n * (i + 1) // k] for i in range(k)]
        if n == 0: self.a = []

    def __iter__(self):
        for b in self.a: yield from b

    def __len__(self): return self.size

    def __contains__(self, x):
        if self.size == 0: return False
        b, _, i = self._pos(x)
        return i != len(b) and b[i] == x

    def _pos(self, x):  # bucket, bucket idx, bisect_left inside bucket
        for j, b in enumerate(self.a):
            if x <= b[-1]: return b, j, bisect_left(b, x)
        return b, j, len(b)

    def add(self, x):
        if self.size == 0:
            self.a = [[x]]; self.size = 1; return
        b, j, i = self._pos(x)
        b.insert(i, x); self.size += 1
        if len(b) > len(self.a) * self.SPLIT_RATIO:
            m = len(b) >> 1
            self.a[j:j + 1] = [b[:m], b[m:]]

    def _pop(self, b, j, i):
        x = b.pop(i); self.size -= 1
        if not b: del self.a[j]
        return x

    def discard(self, x):  # remove one copy; returns False if absent
        if self.size == 0: return False
        b, j, i = self._pos(x)
        if i == len(b) or b[i] != x: return False
        self._pop(b, j, i); return True

    def __getitem__(self, k):  # k-th smallest (0-indexed, negative ok)
        if k < 0: k += self.size
        for b in self.a:
            if k < len(b): return b[k]
            k -= len(b)
        raise IndexError

    def pop(self, k=-1):
        if k < 0: k += self.size
        for j, b in enumerate(self.a):
            if k < len(b): return self._pop(b, j, k)
            k -= len(b)
        raise IndexError

    def index(self, x):  # number of elements < x  (= bisect_left)
        c = 0
        for b in self.a:
            if b[-1] >= x: return c + bisect_left(b, x)
            c += len(b)
        return c

    def index_right(self, x):  # number of elements <= x (= bisect_right)
        c = 0
        for b in self.a:
            if b[-1] > x: return c + bisect_right(b, x)
            c += len(b)
        return c

    def count(self, x): return self.index_right(x) - self.index(x)

    def lt(self, x):  # largest < x, else None (also le/gt/ge via index)
        i = self.index(x); return self[i - 1] if i else None

    def ge(self, x):  # smallest >= x, else None
        i = self.index(x); return self[i] if i < self.size else None


# ==== 4.1 Linear structures ====
# next_greater: O(n). res[i] = smallest j>i with a[j]>a[i], else n. (>= : pop on <=; prev: scan reversed)
def next_greater(a):
    n = len(a); res = [n] * n; st = []
    for i, x in enumerate(a):
        while st and a[st[-1]] < x: res[st.pop()] = i
        st.append(i)
    return res


# largest_rectangle: O(n) max area in histogram h.
def largest_rectangle(h):
    st = []; best = 0  # stack of (start, height), heights increasing
    for i, x in enumerate(h + [0]):
        start = i
        while st and st[-1][1] >= x:
            s, y = st.pop(); best = max(best, y * (i - s)); start = s
        st.append((start, x))
    return best


# sliding_min: O(n) monotonic deque, min of each window a[i-k+1..i] for i>=k-1. (max: flip compare)
def sliding_min(a, k):
    dq = deque(); res = []
    for i, x in enumerate(a):
        while dq and a[dq[-1]] >= x: dq.pop()
        dq.append(i)
        if dq[0] <= i - k: dq.popleft()
        if i >= k - 1: res.append(a[dq[0]])
    return res

# Linked list via arrays: nxt=[i+1..], prv=[i-1..]; delete i: nxt[prv[i]]=nxt[i]; prv[nxt[i]]=prv[i].
#   Use sentinels 0 and n+1. Undo in reverse order (dancing links). O(1) per op, no objects.
# XOR linked list: memory trick only, pointless in Python.


# ==== 4.2 Heaps and priority queues ====
# Binary heap: array, children 2i+1,2i+2; push/pop O(log n), heapify O(n) bottom-up -> use heapq.
# Pairing / binomial heaps: mergeable heaps, rarely needed; leftist heap below covers merging.

# LazyHeap: min-heap with arbitrary deletion O(log n) amortized. remove(x) assumes x is present.
class LazyHeap:
    def __init__(self): self.h = []; self.d = Counter(); self.n = 0

    def push(self, x): heappush(self.h, x); self.n += 1

    def remove(self, x): self.d[x] += 1; self.n -= 1

    def _clean(self):
        h, d = self.h, self.d
        while h and d[h[0]]: d[h[0]] -= 1; heappop(h)

    def top(self): self._clean(); return self.h[0]

    def pop(self): self._clean(); self.n -= 1; return heappop(self.h)

    def __len__(self): return self.n


# RunningMedian: two heaps, add O(log n), median() = lower median. lo is max-heap (negated).
class RunningMedian:
    def __init__(self): self.lo = []; self.hi = []

    def add(self, x):
        if self.lo and x > -self.lo[0]: heappush(self.hi, x)
        else: heappush(self.lo, -x)
        if len(self.lo) > len(self.hi) + 1: heappush(self.hi, -heappop(self.lo))
        elif len(self.hi) > len(self.lo): heappush(self.lo, -heappop(self.hi))

    def median(self): return -self.lo[0]


# LeftistHeap: pool of mergeable min-heaps; merge/pop O(log n). Heap = root id, 0 = empty.
# Recursion depth <= 2 log n (right spines). Use: r = H.merge(r, H.new(x)); r = H.pop(r).
class LeftistHeap:
    def __init__(self): self.key = [0]; self.L = [0]; self.R = [0]; self.D = [0]

    def new(self, x):
        self.key.append(x); self.L.append(0); self.R.append(0); self.D.append(1)
        return len(self.key) - 1

    def merge(self, a, b):
        if not a or not b: return a or b
        K, L, R, D = self.key, self.L, self.R, self.D
        if K[b] < K[a]: a, b = b, a
        R[a] = self.merge(R[a], b)
        if D[L[a]] < D[R[a]]: L[a], R[a] = R[a], L[a]
        D[a] = D[R[a]] + 1
        return a

    def top(self, a): return self.key[a]

    def pop(self, a): return self.merge(self.L[a], self.R[a])
# Skew heap: same merge but always swap children, no D (amortized log n, recursion can go deep).


# ==== 4.3 Hashing ====
# Hash tables: dict/set (open addressing). Anti-hash for int keys: see 00_core.py 4.3/19.
# Rolling / polynomial hash, randomized bases: see 07_strings.py.
# Zobrist: give each (position, state) a random 64-bit r; hash(state) = XOR of r's of its parts;
#   updates are O(1) (xor out old, xor in new). Used for board states / game memo / sets in DP.
# Set hashing: XOR of random values (elements of even multiplicity cancel -> parity sets).
# Multiset hashing: SUM of random values mod 2^64 (or 2^61-1); order-independent, O(1) update.

# multiset_prefix_hash: O(n). a[l:r] and b[l2:r2] equal as multisets (whp) iff
#   P[r]-P[l] == Q[r2]-Q[l2] (mod 2^64), using the SAME rnd dict for both arrays.
def multiset_prefix_hash(a, rnd=None, M=(1 << 64) - 1):
    if rnd is None: rnd = {}
    P = [0]; s = 0
    for x in a:
        v = rnd.get(x)
        if v is None: v = rnd[x] = random.getrandbits(64)
        s = (s + v) & M; P.append(s)
    return P, rnd


# ==== 4.4 Disjoint Set Union ====
# Path compression + union by size: O(alpha(n)) amortized (inverse Ackermann, <= 4 in practice).
class DSU:
    def __init__(self, n): self.p = [-1] * n  # p[x] < 0: root, -size

    def find(self, x):
        p = self.p; r = x
        while p[r] >= 0: r = p[r]
        while p[x] >= 0 and p[x] != r: p[x], x = r, p[x]
        return r

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b: return False
        p = self.p
        if p[a] > p[b]: a, b = b, a
        p[a] += p[b]; p[b] = a
        return True

    def same(self, a, b): return self.find(a) == self.find(b)

    def size(self, x): return -self.p[self.find(x)]


# WeightedDSU (potentials): maintains val[b]-val[a]. union(a,b,w) asserts val[b]-val[a]=w,
# returns False on contradiction. diff(a,b) -> val[b]-val[a] or None. O(alpha) amortized.
# Parity DSU ("same/different color"): same thing with w in {0,1} and ^ instead of +/-.
class WeightedDSU:
    def __init__(self, n): self.p = list(range(n)); self.sz = [1] * n; self.pot = [0] * n  # val[x]-val[p[x]]

    def find(self, x):
        p, pot = self.p, self.pot; path = []
        while p[x] != x: path.append(x); x = p[x]
        for y in reversed(path):  # nearest-to-root first; pot[root] == 0
            pot[y] += pot[p[y]]; p[y] = x
        return x

    def union(self, a, b, w):
        ra, rb = self.find(a), self.find(b)
        pa, pb = self.pot[a], self.pot[b]
        if ra == rb: return pb - pa == w
        if self.sz[ra] < self.sz[rb]: self.p[ra] = rb; self.pot[ra] = pb - pa - w; self.sz[rb] += self.sz[ra]
        else: self.p[rb] = ra; self.pot[rb] = w + pa - pb; self.sz[ra] += self.sz[rb]
        return True

    def diff(self, a, b):
        if self.find(a) != self.find(b): return None
        return self.pot[b] - self.pot[a]


# RollbackDSU: union by size, NO path compression; find O(log n). snapshot() -> token, rollback(token).
# Used for offline dynamic connectivity (segment tree over time, DFS with undo).
class RollbackDSU:
    def __init__(self, n): self.p = [-1] * n; self.hist = []; self.comps = n

    def find(self, x):
        p = self.p
        while p[x] >= 0: x = p[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b: return False
        p = self.p
        if p[a] > p[b]: a, b = b, a
        self.hist.append((a, b, p[b])); p[a] += p[b]; p[b] = a; self.comps -= 1
        return True

    def snapshot(self): return len(self.hist)

    def rollback(self, t):
        p, h = self.p, self.hist
        while len(h) > t:
            a, b, pb = h.pop(); p[a] -= pb; p[b] = pb; self.comps += 1


# Persistent DSU: union by size (no compression) on a persistent array (persistent segment tree)
#   -> O(log^2 n) per op with versions. Offline alternative: DFS over version tree + RollbackDSU.

# ==== 4.5 Range query structures ====
# SparseTable: build O(n log n), query O(1) for IDEMPOTENT op (min, max, gcd, and, or). q(l, r) on [l, r).
class SparseTable:
    def __init__(self, a, op=min):
        self.op = op; t = [list(a)]; j = 1
        while 2 * j <= len(a):
            p = t[-1]; t.append([op(p[i], p[i + j]) for i in range(len(a) - 2 * j + 1)]); j *= 2
        self.t = t

    def query(self, l, r):  # requires l < r
        k = (r - l).bit_length() - 1
        return self.op(self.t[k][l], self.t[k][r - (1 << k)])


# DisjointSparseTable: build O(n log n), O(1) query for ANY associative op (sum, product mod, matrix).
class DisjointSparseTable:
    def __init__(self, a, op):
        self.op = op; n = len(a); a = list(a); self.t = [a]
        for h in range(1, max(1, (n - 1).bit_length())):
            row = a[:]; half = 1 << h
            for mid in range(half, n, 2 * half):
                for i in range(mid - 2, mid - half - 1, -1): row[i] = op(a[i], row[i + 1])
                for i in range(mid + 1, min(mid + half, n)): row[i] = op(row[i - 1], a[i])
            self.t.append(row)

    def query(self, l, r):  # [l, r), l < r
        r -= 1
        if l == r: return self.t[0][l]
        h = (l ^ r).bit_length() - 1
        return self.op(self.t[h][l], self.t[h][r])


# SqrtDecomp: range add + range sum, O(sqrt n) per op. Template for any "block + lazy" trick.
class SqrtDecomp:
    def __init__(self, a):
        n = len(a); self.B = B = max(1, int(n ** 0.5)); self.a = list(a)
        self.s = [sum(a[i:i + B]) for i in range(0, n, B)]; self.lz = [0] * len(self.s)

    def add(self, l, r, x):
        B, a = self.B, self.a
        while l < r and l % B: a[l] += x; self.s[l // B] += x; l += 1
        while l + B <= r: self.lz[l // B] += x; l += B
        while l < r: a[l] += x; self.s[l // B] += x; l += 1

    def sum(self, l, r):
        B, a, res = self.B, self.a, 0
        while l < r and l % B: res += a[l] + self.lz[l // B]; l += 1
        while l + B <= r: res += self.s[l // B] + self.lz[l // B] * B; l += B
        while l < r: res += a[l] + self.lz[l // B]; l += 1
        return res


# Fenwick: point add / prefix sum O(log n), build O(n). kth via binary lifting.
# Range add + point query: Fenwick on differences: add(l,x), add(r,-x); point i = sum(i+1).
class Fenwick:
    def __init__(self, n_or_a):
        if isinstance(n_or_a, int): self.n = n_or_a; self.t = [0] * (n_or_a + 1); return
        self.n = n = len(n_or_a); t = self.t = [0] + list(n_or_a)
        for i in range(1, n + 1):
            j = i + (i & -i)
            if j <= n: t[j] += t[i]

    def add(self, i, x):
        i += 1; t = self.t
        while i <= self.n: t[i] += x; i += i & -i

    def sum(self, i):  # sum of [0, i)
        t = self.t; s = 0
        while i: s += t[i]; i &= i - 1
        return s

    def range_sum(self, l, r): return self.sum(r) - self.sum(l)

    def lower_bound(self, k):  # min i with sum(i+1) >= k (values >= 0); = k-th element (k>=1) of counts
        t = self.t; pos = 0; step = 1 << self.n.bit_length()
        while step:
            if pos + step <= self.n and t[pos + step] < k: pos += step; k -= t[pos]
            step >>= 1
        return pos  # == n if total < k


# RangeFenwick: range add + range sum, O(log n). prefix(i) = sum1(i)*i - sum2(i).
class RangeFenwick:
    def __init__(self, n): self.n = n; self.b1 = [0] * (n + 1); self.b2 = [0] * (n + 1)

    def _add(self, t, i, x):
        i += 1
        while i <= self.n: t[i] += x; i += i & -i

    def add(self, l, r, x):  # a[l:r] += x
        self._add(self.b1, l, x); self._add(self.b1, r, -x)
        self._add(self.b2, l, x * l); self._add(self.b2, r, -x * r)

    def prefix(self, i):  # sum of [0, i)
        s1 = s2 = 0; j = i
        while j: s1 += self.b1[j]; s2 += self.b2[j]; j &= j - 1
        return s1 * i - s2

    def sum(self, l, r): return self.prefix(r) - self.prefix(l)


# Fenwick2D: point add, rectangle sum, O(log n log m). Flattened array. Sparse coords: compress.
class Fenwick2D:
    def __init__(self, n, m): self.n, self.m = n, m; self.t = [0] * ((n + 1) * (m + 1))

    def add(self, x, y, v):
        m1 = self.m + 1; i = x + 1
        while i <= self.n:
            j = y + 1
            while j <= self.m: self.t[i * m1 + j] += v; j += j & -j
            i += i & -i

    def prefix(self, x, y):  # sum over [0,x) x [0,y)
        m1 = self.m + 1; s = 0; i = x
        while i:
            j = y
            while j: s += self.t[i * m1 + j]; j &= j - 1
            i &= i - 1
        return s

    def rect(self, x1, y1, x2, y2):  # [x1,x2) x [y1,y2)
        return self.prefix(x2, y2) - self.prefix(x1, y2) - self.prefix(x2, y1) + self.prefix(x1, y1)


# ==== 4.6 Segment trees ====
# SegTree (ACL-style, iterative bottom-up) over monoid (op, e). set/get/prod O(log n).
# max_right(l, f): largest r with f(prod(l, r)) True (f(e) must be True, f monotone).
# min_left(r, f): smallest l with f(prod(l, r)) True. Custom node: use tuples, e.g. (sum, maxpref).
class SegTree:
    def __init__(self, op, e, v):
        if isinstance(v, int): v = [e] * v
        self.op, self.e, self.n = op, e, len(v)
        self.log = max(0, (self.n - 1).bit_length()); self.size = 1 << self.log
        d = self.d = [e] * (2 * self.size); d[self.size:self.size + self.n] = v
        for i in range(self.size - 1, 0, -1): d[i] = op(d[2 * i], d[2 * i + 1])

    def set(self, p, x):
        p += self.size; d, op = self.d, self.op; d[p] = x
        while p > 1: p >>= 1; d[p] = op(d[2 * p], d[2 * p + 1])

    def get(self, p): return self.d[p + self.size]

    def prod(self, l, r):
        op, d = self.op, self.d; sml = smr = self.e; l += self.size; r += self.size
        while l < r:
            if l & 1: sml = op(sml, d[l]); l += 1
            if r & 1: r -= 1; smr = op(d[r], smr)
            l >>= 1; r >>= 1
        return op(sml, smr)

    def all_prod(self): return self.d[1]

    def max_right(self, l, f):
        if l == self.n: return self.n
        op, d, size = self.op, self.d, self.size; l += size; sm = self.e
        while True:
            while l % 2 == 0: l >>= 1
            if not f(op(sm, d[l])):
                while l < size:
                    l *= 2
                    if f(op(sm, d[l])): sm = op(sm, d[l]); l += 1
                return l - size
            sm = op(sm, d[l]); l += 1
            if l & -l == l: return self.n

    def min_left(self, r, f):
        if r == 0: return 0
        op, d, size = self.op, self.d, self.size; r += size; sm = self.e
        while True:
            r -= 1
            while r > 1 and r % 2: r >>= 1
            if not f(op(d[r], sm)):
                while r < size:
                    r = 2 * r + 1
                    if f(op(d[r], sm)): sm = op(d[r], sm); r -= 1
                return r + 1 - size
            sm = op(d[r], sm)
            if r & -r == r: return 0


# LazySegTree (ACL-style): monoid (op, e), maps F with mapping(f, x), composition(f, g) = f after g,
# identity idf. apply(l, r, f), prod(l, r), set/get, max_right/min_left: all O(log n).
# Ex. range add/range min: op=min, e=INF, mapping=lambda f,x: f+x, comp=lambda f,g: f+g, idf=0.
# Ex. range affine/range sum: x=(sum,len), f=(b,c): (b*s+c*len, len); comp(f,g)=(fb*gb, fb*gc+fc).
class LazySegTree:
    def __init__(self, op, e, mapping, composition, idf, v):
        if isinstance(v, int): v = [e] * v
        self.op, self.e, self.mp, self.cp, self.idf, self.n = op, e, mapping, composition, idf, len(v)
        self.log = max(0, (self.n - 1).bit_length()); self.size = 1 << self.log
        self.d = [e] * (2 * self.size); self.lz = [idf] * self.size
        self.d[self.size:self.size + self.n] = v
        for i in range(self.size - 1, 0, -1): self._up(i)

    def _up(self, k): self.d[k] = self.op(self.d[2 * k], self.d[2 * k + 1])

    def _ap(self, k, f):
        self.d[k] = self.mp(f, self.d[k])
        if k < self.size: self.lz[k] = self.cp(f, self.lz[k])

    def _push(self, k):
        if self.lz[k] != self.idf:
            self._ap(2 * k, self.lz[k]); self._ap(2 * k + 1, self.lz[k]); self.lz[k] = self.idf

    def set(self, p, x):
        p += self.size
        for i in range(self.log, 0, -1): self._push(p >> i)
        self.d[p] = x
        for i in range(1, self.log + 1): self._up(p >> i)

    def get(self, p):
        p += self.size
        for i in range(self.log, 0, -1): self._push(p >> i)
        return self.d[p]

    def prod(self, l, r):
        if l == r: return self.e
        l += self.size; r += self.size
        for i in range(self.log, 0, -1):
            if ((l >> i) << i) != l: self._push(l >> i)
            if ((r >> i) << i) != r: self._push((r - 1) >> i)
        op, d = self.op, self.d; sml = smr = self.e
        while l < r:
            if l & 1: sml = op(sml, d[l]); l += 1
            if r & 1: r -= 1; smr = op(d[r], smr)
            l >>= 1; r >>= 1
        return op(sml, smr)

    def apply(self, l, r, f):
        if l == r: return
        l += self.size; r += self.size
        for i in range(self.log, 0, -1):
            if ((l >> i) << i) != l: self._push(l >> i)
            if ((r >> i) << i) != r: self._push((r - 1) >> i)
        l2, r2 = l, r
        while l < r:
            if l & 1: self._ap(l, f); l += 1
            if r & 1: r -= 1; self._ap(r, f)
            l >>= 1; r >>= 1
        l, r = l2, r2
        for i in range(1, self.log + 1):
            if ((l >> i) << i) != l: self._up(l >> i)
            if ((r >> i) << i) != r: self._up((r - 1) >> i)

    def max_right(self, l, g):
        if l == self.n: return self.n
        l += self.size
        for i in range(self.log, 0, -1): self._push(l >> i)
        op, d, sm = self.op, self.d, self.e
        while True:
            while l % 2 == 0: l >>= 1
            if not g(op(sm, d[l])):
                while l < self.size:
                    self._push(l); l *= 2
                    if g(op(sm, d[l])): sm = op(sm, d[l]); l += 1
                return l - self.size
            sm = op(sm, d[l]); l += 1
            if l & -l == l: return self.n

    def min_left(self, r, g):
        if r == 0: return 0
        r += self.size
        for i in range(self.log, 0, -1): self._push((r - 1) >> i)
        op, d, sm = self.op, self.d, self.e
        while True:
            r -= 1
            while r > 1 and r % 2: r >>= 1
            if not g(op(d[r], sm)):
                while r < self.size:
                    self._push(r); r = 2 * r + 1
                    if g(op(d[r], sm)): sm = op(d[r], sm); r -= 1
                return r + 1 - self.size
            sm = op(d[r], sm)
            if r & -r == r: return 0


# MergeSortTree: build O(n log n) memory O(n log n); count of a[i] <= x in [l, r): O(log^2 n).
class MergeSortTree:
    def __init__(self, a):
        n = len(a); self.size = s = 1 << max(0, (n - 1).bit_length())
        t = self.t = [[] for _ in range(2 * s)]
        for i, x in enumerate(a): t[s + i] = [x]
        for i in range(s - 1, 0, -1): t[i] = sorted(t[2 * i] + t[2 * i + 1])  # Timsort merges runs

    def count_le(self, l, r, x):
        t = self.t; l += self.size; r += self.size; c = 0
        while l < r:
            if l & 1: c += bisect_right(t[l], x); l += 1
            if r & 1: r -= 1; c += bisect_right(t[r], x)
            l >>= 1; r >>= 1
        return c


# PersistentSegTree: k-th smallest in a[l:r] in O(log n); build O(n log n) nodes. Arrays, iterative.
# roots[i] = count tree of a[:i] over compressed values; node 0 = shared null.
# Same pattern gives versioned arrays: copy the root-to-leaf path on each update.
class PersistentKth:
    def __init__(self, a):
        self.vals = vals = sorted(set(a)); idx = {v: i for i, v in enumerate(vals)}; self.m = m = len(vals)
        L, R, C = self.L, self.R, self.C = [0], [0], [0]
        self.roots = roots = [0]
        for x in a:
            pos = idx[x]; prev = roots[-1]; lo, hi = 0, m; roots.append(len(C))
            while True:
                nid = len(C); L.append(L[prev]); R.append(R[prev]); C.append(C[prev] + 1)
                if hi - lo == 1: break
                mid = (lo + hi) >> 1
                if pos < mid: L[nid] = nid + 1; prev = L[prev]; hi = mid
                else: R[nid] = nid + 1; prev = R[prev]; lo = mid

    def kth(self, l, r, k):  # k-th smallest (0-indexed) of a[l:r]
        L, R, C = self.L, self.R, self.C; u, v = self.roots[r], self.roots[l]; lo, hi = 0, self.m
        while hi - lo > 1:
            mid = (lo + hi) >> 1; c = C[L[u]] - C[L[v]]
            if k < c: u, v, hi = L[u], L[v], mid
            else: k -= c; u, v, lo = R[u], R[v], mid
        return self.vals[lo]


# DynamicSegTree: point add / range sum on huge coordinate range [lo, hi), nodes created lazily.
# O(log(hi-lo)) per op, memory O(q log C). Online alternative to coordinate compression.
class DynamicSegTree:
    def __init__(self, lo, hi): self.lo, self.hi = lo, hi; self.L = [0, 0]; self.R = [0, 0]; self.S = [0, 0]

    def _new(self): self.L.append(0); self.R.append(0); self.S.append(0); return len(self.S) - 1

    def add(self, p, x):
        u, lo, hi = 1, self.lo, self.hi
        while True:
            self.S[u] += x
            if hi - lo == 1: return
            mid = (lo + hi) >> 1
            if p < mid:
                if not self.L[u]: self.L[u] = self._new()
                u, hi = self.L[u], mid
            else:
                if not self.R[u]: self.R[u] = self._new()
                u, lo = self.R[u], mid

    def query(self, l, r):
        res = 0; st = [(1, self.lo, self.hi)]; L, R, S = self.L, self.R, self.S
        while st:
            u, a, b = st.pop()
            if not u or b <= l or r <= a: continue
            if l <= a and b <= r: res += S[u]; continue
            mid = (a + b) >> 1; st.append((L[u], a, mid)); st.append((R[u], mid, b))
        return res

# 2D segment tree / segment tree of BITs: outer tree over x, each node holds a structure over y
#   (sorted ys + Fenwick, built offline) -> O(log^2) query/update. In Python prefer offline
#   sweep over x + Fenwick over y (see 01_techniques.py / 11 offline), or Fenwick2D above.


# SegTreeBeats (Ji Driver): chmin(l,r,x), sum(l,r), max(l,r). O((n+q) log n) amortized.
# Recursion depth = log n. Extend with chmax/add by tracking min/second-min the same way.
class SegTreeBeats:
    NEG = -(1 << 62)

    def __init__(self, a):
        n = len(a); s = self.size = 1 << max(0, (n - 1).bit_length()); NEG = self.NEG
        self.mx = [NEG] * (2 * s); self.smx = [NEG] * (2 * s); self.cnt = [0] * (2 * s); self.sm = [0] * (2 * s)
        for i, x in enumerate(a): self.mx[s + i] = x; self.cnt[s + i] = 1; self.sm[s + i] = x
        for k in range(s - 1, 0, -1): self._pull(k)

    def _pull(self, k):
        mx, smx, cnt = self.mx, self.smx, self.cnt; l, r = 2 * k, 2 * k + 1
        self.sm[k] = self.sm[l] + self.sm[r]
        if mx[l] > mx[r]: mx[k] = mx[l]; cnt[k] = cnt[l]; smx[k] = max(smx[l], mx[r])
        elif mx[l] < mx[r]: mx[k] = mx[r]; cnt[k] = cnt[r]; smx[k] = max(mx[l], smx[r])
        else: mx[k] = mx[l]; cnt[k] = cnt[l] + cnt[r]; smx[k] = max(smx[l], smx[r])

    def _tag(self, k, x):  # requires smx[k] < x
        if self.mx[k] > x: self.sm[k] -= (self.mx[k] - x) * self.cnt[k]; self.mx[k] = x

    def _push(self, k):
        self._tag(2 * k, self.mx[k]); self._tag(2 * k + 1, self.mx[k])

    def chmin(self, l, r, x, k=1, a=0, b=None):
        if b is None: b = self.size
        if r <= a or b <= l or self.mx[k] <= x: return
        if l <= a and b <= r and self.smx[k] < x: self._tag(k, x); return
        self._push(k); m = (a + b) >> 1
        self.chmin(l, r, x, 2 * k, a, m); self.chmin(l, r, x, 2 * k + 1, m, b); self._pull(k)

    def _query(self, l, r, k, a, b, is_sum):
        if r <= a or b <= l: return 0 if is_sum else self.NEG
        if l <= a and b <= r: return self.sm[k] if is_sum else self.mx[k]
        self._push(k); m = (a + b) >> 1
        x, y = self._query(l, r, 2 * k, a, m, is_sum), self._query(l, r, 2 * k + 1, m, b, is_sum)
        return x + y if is_sum else max(x, y)

    def sum(self, l, r): return self._query(l, r, 1, 0, self.size, True)

    def max(self, l, r): return self._query(l, r, 1, 0, self.size, False)


# LiChao: min of lines y=a*x+b at query points xs (known in advance, any order). add O(log n),
# query O(log n), iterative. For max: insert (-a,-b) and negate. Used by 03_dp.py CHT.
class LiChao:
    INF = float('inf')

    def __init__(self, xs):
        xs = sorted(set(xs)); self.idx = {x: i for i, x in enumerate(xs)}
        self.size = s = 1 << max(0, (len(xs) - 1).bit_length())
        self.xs = xs + [xs[-1]] * (s - len(xs)); self.A = [0] * (2 * s); self.B = [self.INF] * (2 * s)

    def add_line(self, a, b):
        A, B, xs = self.A, self.B, self.xs; k, lo, hi = 1, 0, self.size
        while True:
            mid = (lo + hi) >> 1; xm = xs[mid]
            if a * xm + b < A[k] * xm + B[k]: A[k], a = a, A[k]; B[k], b = b, B[k]
            if hi - lo == 1 or b == self.INF: return
            xl, xr = xs[lo], xs[hi - 1]
            if a * xl + b < A[k] * xl + B[k]: k, hi = 2 * k, mid
            elif a * xr + b < A[k] * xr + B[k]: k, lo = 2 * k + 1, mid
            else: return

    def query(self, x):
        k = self.idx[x] + self.size; res = self.INF
        while k: res = min(res, self.A[k] * x + self.B[k]); k >>= 1
        return res
# Segment insertion (line only on [l,r)): decompose [l,r) into O(log n) nodes like SegTree.prod
#   and run add_line's descent from each node -> O(log^2 n).


# SegMerge: pool of dynamic count-segment-trees over values [0, M); merge(a,b) destroys both,
# total merge cost O(n log M) amortized. Typical: small-to-large on trees -> merge children's trees.
class SegMerge:
    def __init__(self, M): self.M = M; self.L = [0]; self.R = [0]; self.C = [0]

    def _node(self, c):
        self.L.append(0); self.R.append(0); self.C.append(c); return len(self.C) - 1

    def single(self, v):  # new tree containing value v once
        root = u = self._node(1); lo, hi = 0, self.M
        while hi - lo > 1:
            mid = (lo + hi) >> 1; w = self._node(1)
            if v < mid: self.L[u] = w; hi = mid
            else: self.R[u] = w; lo = mid
            u = w
        return root

    def merge(self, a, b):  # recursion depth log M
        if not a or not b: return a or b
        self.C[a] += self.C[b]
        self.L[a] = self.merge(self.L[a], self.L[b]); self.R[a] = self.merge(self.R[a], self.R[b])
        return a

    def kth(self, u, k):  # k-th smallest (0-indexed)
        lo, hi = 0, self.M
        while hi - lo > 1:
            mid = (lo + hi) >> 1; c = self.C[self.L[u]]
            if k < c: u, hi = self.L[u], mid
            else: k -= c; u, lo = self.R[u], mid
        return lo
# Split (by value k): walk down, detaching right parts into new nodes; same O(log M).

# Kinetic segment tree: each node stores current max line plus "melt time" when the winner changes;
#   supports heaten(l,r,dt) (x += dt for all) with amortized O(log^2 n). Rare; Li Chao/CHT usually suffice.


# ==== 4.7 Balanced trees ====
# BST (unbalanced): O(depth) ops, degenerates on sorted input -> never use plain in contests.
# Order statistics tree / ordered set: use SortedList (2.2 above), or Fenwick.lower_bound on
#   compressed values (offline), or the treap below with size fields.

# ImplicitTreap: sequence with insert/erase at position, reverse(l,r), sum(l,r), O(log n) expected.
# Arrays-based; split/merge recurse to depth O(log n) expected (~3 log n), safe. Also = rope.
# Add more lazies (range add, assign) in _push/_pull the same way as rev.
class ImplicitTreap:
    def __init__(self, a=()):
        self.val = [0]; self.pri = [0]; self.L = [0]; self.R = [0]; self.sz = [0]; self.sm = [0]; self.rev = [0]
        self.root = 0
        for x in a: self.root = self.merge(self.root, self._new(x))

    def _new(self, x):
        self.val.append(x); self.pri.append(random.random()); self.L.append(0); self.R.append(0)
        self.sz.append(1); self.sm.append(x); self.rev.append(0)
        return len(self.val) - 1

    def _pull(self, t):
        l, r = self.L[t], self.R[t]
        self.sz[t] = self.sz[l] + 1 + self.sz[r]; self.sm[t] = self.sm[l] + self.val[t] + self.sm[r]

    def _push(self, t):
        if self.rev[t]:
            l, r = self.L[t], self.R[t]; self.L[t], self.R[t] = r, l
            if l: self.rev[l] ^= 1
            if r: self.rev[r] ^= 1
            self.rev[t] = 0

    def split(self, t, k):  # -> (first k elements, rest)
        if not t: return 0, 0
        self._push(t)
        if self.sz[self.L[t]] >= k:
            a, b = self.split(self.L[t], k); self.L[t] = b; self._pull(t); return a, t
        a, b = self.split(self.R[t], k - self.sz[self.L[t]] - 1); self.R[t] = a; self._pull(t); return t, b

    def merge(self, a, b):
        if not a or not b: return a or b
        if self.pri[a] > self.pri[b]:
            self._push(a); self.R[a] = self.merge(self.R[a], b); self._pull(a); return a
        self._push(b); self.L[b] = self.merge(a, self.L[b]); self._pull(b); return b

    def __len__(self): return self.sz[self.root]

    def insert(self, pos, x):
        a, b = self.split(self.root, pos); self.root = self.merge(self.merge(a, self._new(x)), b)

    def erase(self, pos):
        a, b = self.split(self.root, pos); _, c = self.split(b, 1); self.root = self.merge(a, c)

    def reverse(self, l, r):
        a, b = self.split(self.root, l); b, c = self.split(b, r - l)
        self.rev[b] ^= 1; self.root = self.merge(self.merge(a, b), c)

    def sum(self, l, r):
        a, b = self.split(self.root, l); b, c = self.split(b, r - l)
        s = self.sm[b]; self.root = self.merge(self.merge(a, b), c); return s

    def to_list(self):  # iterative in-order
        res, st, t = [], [], self.root
        while st or t:
            while t: self._push(t); st.append(t); t = self.L[t]
            t = st.pop(); res.append(self.val[t]); t = self.R[t]
        return res
# Keyed treap (ordered set): split by key instead of size: go left if key[t] >= k.

# Splay tree: amortized O(log n) by rotating accessed node to root (zig-zig/zig-zag); used inside
#   link-cut trees (4.8). As a sequence tree it is equivalent to the implicit treap.
# AVL (|h(l)-h(r)|<=1) / red-black (black-height equal, no red-red): worst-case O(log n), conceptual.
# Scapegoat: rebuild the subtree whose size violates alpha-balance; amortized O(log n), no rotations.
# Skip list: randomized levels, expected O(log n) search; no advantage over SortedList in Python.


# ==== 4.8 Advanced and specialized structures ====
# Sqrt tree: O(1) query for any associative op with O(n log log n) build; DisjointSparseTable
#   (4.5) gives O(1) with O(n log n) build, enough in practice.

# WaveletMatrix: nonneg ints. kth(l,r,k) k-th smallest in a[l:r]; count_less(l,r,x) #{a[i]<x};
# rank(x, r) = #x in a[:r]. Build O(n log V), queries O(log V). Range freq [lo,hi) = cl(hi)-cl(lo).
class WaveletMatrix:
    def __init__(self, a):
        self.lg = lg = max(1, max(a, default=0).bit_length()); n = len(a)
        self.pre = []; self.z = []; cur = list(a)
        for b in range(lg - 1, -1, -1):
            p = [0] * (n + 1); s = 0; zs = []; os_ = []
            for i, x in enumerate(cur):
                if x >> b & 1: s += 1; os_.append(x)
                else: zs.append(x)
                p[i + 1] = s
            self.pre.append(p); self.z.append(len(zs)); cur = zs + os_

    def kth(self, l, r, k):
        res = 0
        for lv in range(self.lg):
            p = self.pre[lv]; ol, orr = p[l], p[r]; zc = (r - l) - (orr - ol)
            if k < zc: l -= ol; r -= orr
            else: k -= zc; res |= 1 << (self.lg - 1 - lv); l = self.z[lv] + ol; r = self.z[lv] + orr
        return res

    def count_less(self, l, r, x):
        if x >= 1 << self.lg: return r - l
        if x <= 0: return 0
        res = 0
        for lv in range(self.lg):
            p = self.pre[lv]; ol, orr = p[l], p[r]
            if x >> (self.lg - 1 - lv) & 1: res += (r - l) - (orr - ol); l = self.z[lv] + ol; r = self.z[lv] + orr
            else: l -= ol; r -= orr
        return res

    def rank(self, x, r): return self.count_less(0, r, x + 1) - self.count_less(0, r, x)
# Wavelet tree: same idea as a tree over value ranges; matrix is simpler & faster.
# van Emde Boas: O(log log U) pred/succ via sqrt(U) clusters + summary (conceptual). In Python a
#   64-ary bitset trie (levels of ints, use bit_length) gives fast pred/succ for U ~ 1e6..1e7.


# LinkCutTree: dynamic forest, nodes 1..n. link(u,v) (different trees), cut(u,v) (edge exists),
# connected, path_sum(u,v), set_val. O(log n) amortized; iterative splay. Non-commutative
# aggregate: also keep reversed aggregate and swap it in _push.
class LinkCutTree:
    def __init__(self, n, vals=None):
        self.L = [0] * (n + 1); self.R = [0] * (n + 1); self.P = [0] * (n + 1); self.rev = [0] * (n + 1)
        self.val = [0] + (list(vals) if vals else [0] * n); self.agg = self.val[:]

    def _isroot(self, x):
        p = self.P[x]; return p == 0 or (self.L[p] != x and self.R[p] != x)

    def _pull(self, x): self.agg[x] = self.agg[self.L[x]] + self.val[x] + self.agg[self.R[x]]

    def _push(self, x):
        if self.rev[x]:
            l, r = self.L[x], self.R[x]; self.L[x], self.R[x] = r, l
            self.rev[l] ^= 1; self.rev[r] ^= 1; self.rev[x] = 0; self.rev[0] = 0

    def _rot(self, x):
        L, R, P = self.L, self.R, self.P; p = P[x]; g = P[p]
        if not self._isroot(p):
            if L[g] == p: L[g] = x
            else: R[g] = x
        P[x] = g
        if L[p] == x: b = R[x]; L[p] = b; R[x] = p
        else: b = L[x]; R[p] = b; L[x] = p
        if b: P[b] = p
        P[p] = x; self._pull(p); self._pull(x)

    def _splay(self, x):
        st = [x]; y = x
        while not self._isroot(y): y = self.P[y]; st.append(y)
        for y in reversed(st): self._push(y)
        P, L = self.P, self.L
        while not self._isroot(x):
            p = P[x]
            if not self._isroot(p):
                g = P[p]
                self._rot(p if (L[g] == p) == (L[p] == x) else x)
            self._rot(x)

    def _access(self, x):
        last, y = 0, x
        while y: self._splay(y); self.R[y] = last; self._pull(y); last = y; y = self.P[y]
        self._splay(x)

    def evert(self, x): self._access(x); self.rev[x] ^= 1  # make x the root of its tree

    def find_root(self, x):
        self._access(x)
        while True:
            self._push(x)
            if not self.L[x]: break
            x = self.L[x]
        self._splay(x); return x

    def connected(self, u, v): return u == v or self.find_root(u) == self.find_root(v)

    def link(self, u, v): self.evert(u); self.P[u] = v

    def cut(self, u, v):
        self.evert(u); self._access(v)
        self.L[v] = 0; self.P[u] = 0; self._pull(v)

    def path_sum(self, u, v): self.evert(u); self._access(v); return self.agg[v]

    def set_val(self, x, v): self._access(x); self.val[x] = v; self._pull(x)

# Euler tour tree: forest as Euler tour sequences in a balanced BST (treap keyed by implicit
#   position); link/cut = O(1) splits/merges; subtree aggregates & connectivity in O(log n).
# Top trees: cluster tree over contractions (rake/compress); path + subtree queries on dynamic
#   forests in O(log n). Conceptual; LCT with virtual-subtree sums covers most uses.
# Persistent structures: path copying (copy the O(log n) nodes on the update path, keep old root;
#   see PersistentKth) or fat nodes (each field stores (version, value) list, bisect by version).
#   Persistent array = persistent segment tree with leaves as cells. Treap persistence: copy in split/merge.
# Retroactive structures: insert/delete operations in the past; partial (query present) vs full.
#   Conceptual; in contests do offline: segment tree over time + rollback (RollbackDSU).
# Rope: sequence with O(log n) cut/concat/insert of substrings = ImplicitTreap (4.7).


# cartesian_tree: O(n) stack build of min-Cartesian tree; returns parent (-1 for root).
# Root = leftmost minimum (ties: earlier index is ancestor). RMQ(l,r) = LCA(l,r) in this tree,
# and LCA reduces to RMQ on Euler tour (see 06_trees.py).
def cartesian_tree(a):
    n = len(a); par = [-1] * n; st = []
    for i in range(n):
        last = -1
        while st and a[st[-1]] > a[i]: last = st.pop()
        if last != -1: par[last] = i
        if st: par[i] = st[-1]
        st.append(i)
    return par

# Bitset via Python int: b |= b << w (subset sum), b & mask, popcount int.bit_count() (3.10+) or
#   bin(b).count('1'), lowest bit (b & -b).bit_length()-1, highest b.bit_length()-1. Shifts/ands
#   are C-speed (~n/64 words): n^2/64-style algorithms (reachability rows, triangle counting).
# Fractional cascading: k sorted lists searched for x in O(log n + k) by storing every 2nd element
#   of list i+1 inside list i with bridge pointers. Mostly conceptual; merge sort tree with it
#   gives O(log n) per query. In Python, bisect on each list is usually fast enough.


# ==== tests ====
if __name__ == "__main__":
    import functools, operator
    rng = random.Random(12345); random.seed(12345)
    R = rng.randint

    # SortedList
    for _ in range(30):
        s = SortedList([R(0, 50) for _ in range(R(0, 40))]); b = sorted(s)
        for _ in range(400):
            t, x = R(0, 5), R(0, 50)
            if t <= 1: s.add(x); insort(b, x)
            elif t == 2: assert s.discard(x) == (x in b); (b.remove(x) if x in b else None)
            else:
                assert s.index(x) == bisect_left(b, x) and s.index_right(x) == bisect_right(b, x)
                assert (x in s) == (x in b) and s.count(x) == b.count(x)
                if b:
                    k = R(-len(b), len(b) - 1); assert s[k] == b[k]
                    i = bisect_left(b, x); assert s.lt(x) == (b[i - 1] if i else None)
                    assert s.ge(x) == (b[i] if i < len(b) else None)
            assert len(s) == len(b)
        assert list(s) == b
        while b: k = R(0, len(b) - 1); assert s.pop(k) == b.pop(k)
    big = SortedList()
    for _ in range(20000): big.add(R(0, 10 ** 9))
    assert list(big) == sorted(big) and len(big.a) > 1

    # 4.1
    for _ in range(200):
        a = [R(0, 6) for _ in range(R(1, 15))]; n = len(a)
        assert next_greater(a) == [next((j for j in range(i + 1, n) if a[j] > a[i]), n) for i in range(n)]
        assert largest_rectangle(a) == max(min(a[i:j]) * (j - i) for i in range(n) for j in range(i + 1, n + 1))
        k = R(1, n); assert sliding_min(a, k) == [min(a[i - k + 1:i + 1]) for i in range(k - 1, n)]

    # 4.2
    h = LazyHeap(); ref = []
    for _ in range(2000):
        t = R(0, 3)
        if t <= 1 or not ref: x = R(0, 30); h.push(x); ref.append(x)
        elif t == 2: x = rng.choice(ref); h.remove(x); ref.remove(x)
        else: assert h.top() == min(ref); assert h.pop() == min(ref); ref.remove(min(ref))
        assert len(h) == len(ref)
    rm = RunningMedian(); ref = []
    for _ in range(500):
        x = R(-100, 100); rm.add(x); insort(ref, x); assert rm.median() == ref[(len(ref) - 1) // 2]
    H = LeftistHeap(); roots = []; refs = []
    for _ in range(1000):
        t = R(0, 2)
        if t == 0 or len(roots) < 2: x = R(0, 99); roots.append(H.new(x)); refs.append([x])
        elif t == 1:
            i, j = rng.sample(range(len(roots)), 2); roots[i] = H.merge(roots[i], roots[j]); refs[i] += refs[j]
            roots.pop(j); refs.pop(j)
        else:
            i = R(0, len(roots) - 1); assert H.top(roots[i]) == min(refs[i])
            refs[i].remove(min(refs[i])); roots[i] = H.pop(roots[i])
            if not refs[i]: assert roots[i] == 0; roots.pop(i); refs.pop(i)

    # 4.3
    for _ in range(100):
        a = [R(0, 5) for _ in range(12)]; b = [R(0, 5) for _ in range(12)]
        P, rnd = multiset_prefix_hash(a); Q, _ = multiset_prefix_hash(b, rnd)
        for _ in range(20):
            l, r = sorted(rng.sample(range(13), 2)); l2 = R(0, 12 - (r - l)); r2 = l2 + r - l
            assert (P[r] - P[l] - Q[r2] + Q[l2]) % (1 << 64) == 0 or sorted(a[l:r]) != sorted(b[l2:r2])
            assert (P[r] - P[l] - Q[r2] + Q[l2]) % (1 << 64) != 0 or sorted(a[l:r]) == sorted(b[l2:r2])

    # 4.4
    def comps_ref(n, edges):
        lab = list(range(n))
        for _ in range(n):
            for u, v in edges: m = min(lab[u], lab[v]); lab[u] = lab[v] = m
        return lab
    for _ in range(100):
        n = R(1, 15); d = DSU(n); rb = RollbackDSU(n); edges = []; snaps = []
        for _ in range(30):
            u, v = R(0, n - 1), R(0, n - 1); t = R(0, 3)
            if t == 0: snaps.append((rb.snapshot(), len(edges)))
            elif t == 1 and snaps: tk, ne = snaps.pop(); rb.rollback(tk); del edges[ne:]
            else: rb.union(u, v); edges.append((u, v))
            lab = comps_ref(n, edges)
            assert all((rb.find(i) == rb.find(j)) == (lab[i] == lab[j]) for i in range(n) for j in range(n))
            assert rb.comps == len(set(lab))
        allE = [(R(0, n - 1), R(0, n - 1)) for _ in range(20)]
        for u, v in allE: d.union(u, v)
        lab = comps_ref(n, allE)
        assert all(d.same(i, j) == (lab[i] == lab[j]) for i in range(n) for j in range(n))
        assert all(d.size(i) == lab.count(lab[i]) for i in range(n))
        val = [R(-50, 50) for _ in range(n)]; w = WeightedDSU(n); E = []
        for _ in range(25):
            u, v = R(0, n - 1), R(0, n - 1)
            if R(0, 4) == 0:
                ok = w.union(u, v, val[v] - val[u] + 1)
                lab = comps_ref(n, E); assert ok == (lab[u] != lab[v])
                if ok: val = None; break
            else: assert w.union(u, v, val[v] - val[u]); E.append((u, v))
            lab = comps_ref(n, E)
            for i in range(n):
                for j in range(n): assert w.diff(i, j) == (val[j] - val[i] if lab[i] == lab[j] else None)

    # 4.5
    for _ in range(100):
        n = R(1, 40); a = [R(-20, 20) for _ in range(n)]
        st = SparseTable(a); sg = SparseTable([abs(x) for x in a], math.gcd)
        ds = DisjointSparseTable(a, lambda x, y: x + y); dm = DisjointSparseTable([(x,) for x in a], operator.add)  # non-commutative
        for _ in range(50):
            l = R(0, n - 1); r = R(l + 1, n)
            assert st.query(l, r) == min(a[l:r]) and ds.query(l, r) == sum(a[l:r])
            assert sg.query(l, r) == functools.reduce(math.gcd, [abs(x) for x in a[l:r]])
            assert dm.query(l, r) == tuple(a[l:r])
        sq = SqrtDecomp(a); rf = RangeFenwick(n); fw = Fenwick(a); b = a[:]
        for i, x in enumerate(a): rf.add(i, i + 1, x)
        for _ in range(100):
            l = R(0, n - 1); r = R(l, n); x = R(-5, 5)
            if R(0, 1): sq.add(l, r, x); rf.add(l, r, x); b[l:r] = [y + x for y in b[l:r]]
            else: assert sq.sum(l, r) == rf.sum(l, r) == sum(b[l:r])
            i = R(0, n - 1); fw.add(i, x); a[i] += x
            assert fw.range_sum(l, r) == sum(a[l:r])
        cnt = [R(0, 3) for _ in range(n)]; f = Fenwick(cnt)
        for k in range(1, sum(cnt) + 2):
            ref = next((i for i in range(n) if sum(cnt[:i + 1]) >= k), n); assert f.lower_bound(k) == ref
    for _ in range(30):
        n, m = R(1, 8), R(1, 8); f2 = Fenwick2D(n, m); g = [[0] * m for _ in range(n)]
        for _ in range(30):
            x, y, v = R(0, n - 1), R(0, m - 1), R(-9, 9); f2.add(x, y, v); g[x][y] += v
            x1, x2 = sorted((R(0, n), R(0, n))); y1, y2 = sorted((R(0, m), R(0, m)))
            assert f2.rect(x1, y1, x2, y2) == sum(g[i][j] for i in range(x1, x2) for j in range(y1, y2))

    # 4.6
    MOD = 998244353
    for _ in range(60):
        n = R(1, 33); a = [R(0, 10) for _ in range(n)]
        sg = SegTree(operator.add, 0, a)
        for _ in range(60):
            t = R(0, 3)
            if t == 0: i = R(0, n - 1); a[i] = R(0, 10); sg.set(i, a[i])
            elif t == 1: l = R(0, n); r = R(l, n); assert sg.prod(l, r) == sum(a[l:r])
            else:
                K = R(0, 40); l = R(0, n); r = R(0, n)
                assert sg.max_right(l, lambda s: s <= K) == max(j for j in range(l, n + 1) if sum(a[l:j]) <= K)
                assert sg.min_left(r, lambda s: s <= K) == min(j for j in range(0, r + 1) if sum(a[j:r]) <= K)
        # affine / sum
        lz = LazySegTree(lambda x, y: ((x[0] + y[0]) % MOD, x[1] + y[1]), (0, 0),
                         lambda f, x: ((f[0] * x[0] + f[1] * x[1]) % MOD, x[1]),
                         lambda f, g: (f[0] * g[0] % MOD, (f[0] * g[1] + f[1]) % MOD), (1, 0), [(x, 1) for x in a])
        # add / min, plus max_right
        lm = LazySegTree(min, 10 ** 18, lambda f, x: f + x, lambda f, g: f + g, 0, a[:])
        b = a[:]; c = a[:]
        for _ in range(60):
            l = R(0, n); r = R(l, n); t = R(0, 4)
            if t == 0: p, q = R(0, 3), R(0, 5); lz.apply(l, r, (p, q)); b[l:r] = [(p * x + q) % MOD for x in b[l:r]]
            elif t == 1: assert lz.prod(l, r)[0] == sum(b[l:r]) % MOD
            elif t == 2: x = R(-3, 3); lm.apply(l, r, x); c[l:r] = [y + x for y in c[l:r]]
            elif t == 3 and l < r: assert lm.prod(l, r) == min(c[l:r]); i = R(0, n - 1); assert lm.get(i) == c[i]
            else:
                K = R(-5, 10)
                assert lm.max_right(l, lambda v: v >= K) == max(j for j in range(l, n + 1) if all(y >= K for y in c[l:j]))
                assert lm.min_left(r, lambda v: v >= K) == min(j for j in range(r + 1) if all(y >= K for y in c[j:r]))
                i = R(0, n - 1); v = R(-5, 5); lm.set(i, v); c[i] = v
        ms = MergeSortTree(a); pk = PersistentKth(a); wm = WaveletMatrix(a); bt = SegTreeBeats(a); d = a[:]
        for _ in range(60):
            l = R(0, n - 1); r = R(l + 1, n); x = R(-1, 11); srt = sorted(a[l:r])
            assert ms.count_le(l, r, x) == sum(y <= x for y in a[l:r])
            k = R(0, r - l - 1); assert pk.kth(l, r, k) == srt[k] == wm.kth(l, r, k)
            assert wm.count_less(l, r, x) == sum(y < x for y in a[l:r]) and wm.rank(x, r) == a[:r].count(x)
            if R(0, 1): bt.chmin(l, r, x); d[l:r] = [min(y, x) for y in d[l:r]]
            else: assert bt.sum(l, r) == sum(d[l:r]) and bt.max(l, r) == max(d[l:r])
    for _ in range(30):
        ds_ = DynamicSegTree(-10 ** 9, 10 ** 9); pts = Counter()
        for _ in range(100):
            p, x = R(-30, 30) * 10 ** 7, R(-5, 5); ds_.add(p, x); pts[p] += x
            l, r = sorted((R(-31, 31) * 10 ** 7, R(-31, 31) * 10 ** 7))
            assert ds_.query(l, r) == sum(v for q, v in pts.items() if l <= q < r)
        xs = [R(-50, 50) for _ in range(R(1, 20))]; lc = LiChao(xs); lines = []
        for _ in range(30):
            ln = (R(-10, 10), R(-100, 100)); lc.add_line(*ln); lines.append(ln)
            x = rng.choice(xs); assert lc.query(x) == min(p * x + q for p, q in lines)
        M = 40; sm_ = SegMerge(M); trees = []; refs = []
        for _ in range(100):
            if len(trees) < 2 or R(0, 2) == 0: v = R(0, M - 1); trees.append(sm_.single(v)); refs.append([v])
            else:
                i, j = rng.sample(range(len(trees)), 2); trees[i] = sm_.merge(trees[i], trees[j]); refs[i] += refs[j]
                trees.pop(j); refs.pop(j)
            i = R(0, len(trees) - 1); srt = sorted(refs[i]); k = R(0, len(srt) - 1)
            assert sm_.kth(trees[i], k) == srt[k]

    # 4.7
    for _ in range(50):
        a = [R(-9, 9) for _ in range(R(0, 20))]; tr = ImplicitTreap(a)
        for _ in range(100):
            t = R(0, 3); n = len(a)
            if t == 0: p, x = R(0, n), R(-9, 9); tr.insert(p, x); a.insert(p, x)
            elif t == 1 and n: p = R(0, n - 1); tr.erase(p); a.pop(p)
            else:
                l = R(0, n); r = R(l, n)
                if t == 2: tr.reverse(l, r); a[l:r] = a[l:r][::-1]
                else: assert tr.sum(l, r) == sum(a[l:r])
            assert len(tr) == len(a)
        assert tr.to_list() == a
    big = ImplicitTreap(range(50000)); big.reverse(0, 50000); assert big.to_list()[:3] == [49999, 49998, 49997]

    # 4.8 link-cut vs brute forest
    for _ in range(60):
        n = R(1, 12); vals = [R(-9, 9) for _ in range(n)]; lct = LinkCutTree(n, vals)
        adj = [set() for _ in range(n + 1)]; val = [0] + vals

        def path(u, v):
            prev = {u: 0}; dq = deque([u])
            while dq:
                x = dq.popleft()
                for y in adj[x]:
                    if y not in prev: prev[y] = x; dq.append(y)
            if v not in prev: return None
            res = []
            while v: res.append(v); v = prev[v]
            return res
        for _ in range(150):
            u, v = R(1, n), R(1, n); t = R(0, 4); pth = path(u, v)
            assert lct.connected(u, v) == (pth is not None)
            if t == 0 and pth is None: lct.link(u, v); adj[u].add(v); adj[v].add(u)
            elif t == 1 and v in adj[u]: lct.cut(u, v); adj[u].discard(v); adj[v].discard(u)
            elif t == 2: x = R(-9, 9); lct.set_val(u, x); val[u] = x
            elif pth is not None: assert lct.path_sum(u, v) == sum(val[x] for x in pth)
    for _ in range(200):
        a = [R(0, 5) for _ in range(R(1, 12))]; par = [-1] * len(a)

        def build(l, r, p):
            st = [(l, r, p)]
            while st:
                l, r, p = st.pop()
                if l >= r: continue
                m = min(range(l, r), key=lambda i: (a[i], i)); par[m] = p
                st.append((l, m, m)); st.append((m + 1, r, m))
        build(0, len(a), -1); assert cartesian_tree(a) == par
    print("ALL OK")
