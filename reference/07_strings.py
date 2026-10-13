# 07_strings.py -- strings: KMP, Z, rolling hash, Aho-Corasick, bitset matching, tries, suffix array
# (SA-IS + doubling) + LCP, suffix automaton, eertree, Manacher, Lyndon/min rotation, classic problems.
# Python 3 stdlib only, PyPy-friendly. Functions take str (most also work on lists of ints).
import random

# ==== 8.1 String basics (notes) ====
# - str is immutable: build with list + "".join(); s[::-1] reverses; bytes ops (translate, count,
#   find) run in C and are much faster than Python loops. s.find/s.count are fine for one pattern
#   (count is non-overlapping!). Map chars: [ord(c) - 97 for c in s].
# - Lexicographic: Python compares str/list/tuple lexicographically; a proper prefix is smaller.
#   Sorting suffixes/rotations with key=lambda i: s[i:] is O(n^2 log n) -- only for brute force.
# - Palindromes: s == s[::-1]; all palindromic substrings -> Manacher / eertree below;
#   interval DP pal[i][j] = s[i]==s[j] and pal[i+1][j-1] for O(n^2).
# - Subsequence (not contiguous) vs substring (contiguous). n(n+1)/2 substrings, 2^n subsequences.
#   Is s a subsequence of t: it = iter(t); all(c in it for c in s)   (O(|t|)).
#   Many queries: nxt[i][c] = first position >= i of c in t (O(26|t|) table), then greedy jumps.
#   Distinct subsequences DP: dp[i+1] = 2*dp[i] - dp[last[c]] (last = index before prev occurrence).


# ==== 8.2 Naive matching ====
# O(nm) worst; in Python prefer str.find in a loop: i = t.find(p, i + 1) (C speed, overlapping).


# ==== 8.2 Prefix function (KMP) + search + automaton ====
def prefix_function(s):                 # O(n). pi[i] = longest proper border of s[:i+1]
    n = len(s)
    pi = [0] * n
    for i in range(1, n):
        k = pi[i - 1]
        while k and s[i] != s[k]:
            k = pi[k - 1]
        if s[i] == s[k]:
            k += 1
        pi[i] = k
    return pi


def kmp_search(text, pat):              # O(n+m) -> start indices of (overlapping) occurrences
    m = len(pat)
    if m == 0:
        return list(range(len(text) + 1))
    pi = prefix_function(pat)
    res = []
    k = 0
    for i, c in enumerate(text):
        while k and c != pat[k]:
            k = pi[k - 1]
        if c == pat[k]:
            k += 1
        if k == m:
            res.append(i - m + 1)
            k = pi[k - 1]
    return res


def kmp_automaton(s, alpha="abcdefghijklmnopqrstuvwxyz"):
    # O(n * |alpha|). aut[i][c] = matched length after reading char index c in state i (0..n).
    # Use for DP over "strings avoiding s as substring" etc.
    n = len(s)
    pi = prefix_function(s)
    A = len(alpha)
    idx = {ch: j for j, ch in enumerate(alpha)}
    aut = [[0] * A for _ in range(n + 1)]
    for i in range(n + 1):
        for c in range(A):
            if i < n and idx[s[i]] == c:
                aut[i][c] = i + 1
            else:
                aut[i][c] = aut[pi[i - 1]][c] if i else 0
    return aut


# ==== 8.2 Z-function ====
def z_function(s):                      # O(n). z[i] = lcp(s, s[i:]), z[0] = n
    n = len(s)
    z = [0] * n
    if n:
        z[0] = n
    l = r = 0
    for i in range(1, n):
        if i < r:
            z[i] = min(r - i, z[i - l])
        while i + z[i] < n and s[z[i]] == s[i + z[i]]:
            z[i] += 1
        if i + z[i] > r:
            l, r = i, i + z[i]
    return z
# Matching: z_function(p + "\x00" + t), occurrence at i where z[m+1+i] >= m.


# ==== 4.3 / 8.4 Polynomial rolling hash mod 2^61-1 (THE library rolling hash) ====
# Build O(n), substring hash O(1). Random base -> safe against anti-hash tests. Collision prob per
# comparison ~ n / 2^61. Hashes from different RollingHash objects are comparable (shared BASE).
class RollingHash:
    MOD = (1 << 61) - 1
    BASE = random.randrange(1 << 40, (1 << 61) - 2)

    def __init__(self, s):
        M, B = self.MOD, self.BASE
        self.s = s
        n = len(s)
        h = [0] * (n + 1)
        p = [1] * (n + 1)
        for i, c in enumerate(s):
            h[i + 1] = (h[i] * B + (ord(c) if type(c) is str else c + 1)) % M
            p[i + 1] = p[i] * B % M
        self.h, self.p = h, p

    def get(self, l, r):                # hash of s[l:r]
        return (self.h[r] - self.h[l] * self.p[r - l]) % self.MOD

    @staticmethod
    def combine(h1, h2, len2):          # hash(A + B) from hash(A), hash(B), |B| (needs pow)
        M = RollingHash.MOD
        return (h1 * pow(RollingHash.BASE, len2, M) + h2) % M


def hash_lcp(A, i, B, j):               # O(log n): lcp of A.s[i:] and B.s[j:] (A, B RollingHash)
    lo, hi = 0, min(len(A.s) - i, len(B.s) - j)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if A.get(i, i + mid) == B.get(j, j + mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


def hash_compare(A, l1, r1, B, l2, r2):  # O(log n): -1/0/1 comparing A.s[l1:r1] vs B.s[l2:r2]
    k = min(hash_lcp(A, l1, B, l2), r1 - l1, r2 - l2)
    if k == r1 - l1 or k == r2 - l2:
        a, b = r1 - l1, r2 - l2
    else:
        a, b = A.s[l1 + k], B.s[l2 + k]
    return (a > b) - (a < b)
# Equal substrings test: A.get(l1, r1) == B.get(l2, r2). Palindrome test: hash s and s[::-1].
# Hashing sets/multisets / Zobrist: see 02_data_structures.py 4.3 (multiset_prefix_hash; anti-hash salt in 00_core.py).


# ==== 8.2 Rabin-Karp (via RollingHash) ====
def rabin_karp(text, pat):              # O(n+m) expected -> start indices
    m = len(pat)
    H = RollingHash(text)
    hp = RollingHash(pat).get(0, m)
    return [i for i in range(len(text) - m + 1) if H.get(i, i + m) == hp]
# Many patterns of the same length: put their hashes in a set, slide once.


# ==== 8.2 Aho-Corasick (goto table, BFS, count occurrences of every pattern) ====
# Build O(total_len * A), scan O(|text|). Alphabet = lowercase by default.
class AhoCorasick:
    def __init__(self, pats, A=26, base=97):
        self.A, self.base = A, base
        go = [[-1] * A]
        term = []
        for p in pats:
            v = 0
            for ch in p:
                c = ord(ch) - base
                if go[v][c] < 0:
                    go[v][c] = len(go)
                    go.append([-1] * A)
                v = go[v][c]
            term.append(v)              # duplicates share a node
        fail = [0] * len(go)
        q = []
        for c in range(A):
            u = go[0][c]
            if u < 0:
                go[0][c] = 0
            else:
                q.append(u)
        for v in q:                     # BFS; q grows
            gv, gf = go[v], go[fail[v]]
            for c in range(A):
                u = gv[c]
                if u < 0:
                    gv[c] = gf[c]
                else:
                    fail[u] = gf[c]
                    q.append(u)
        self.go, self.fail, self.term, self.order = go, fail, term, q

    def count(self, text):              # occurrences (overlapping) of each pattern in text
        go, base = self.go, self.base
        cnt = [0] * len(go)
        v = 0
        for ch in text:
            v = go[v][ord(ch) - base]
            cnt[v] += 1
        fail = self.fail
        for v in reversed(self.order):  # push counts up the fail tree
            cnt[fail[v]] += cnt[v]
        return [cnt[t] for t in self.term]
# "Does any pattern end here": precompute bad[v] |= bad[fail[v]] in BFS order.
# DP over automaton states (go table) for counting strings avoiding/containing patterns.


# ==== 8.2 Bitset matching (Shift-And) with Python ints ====
def shift_and(text, pat):               # O(n * m/64). Classic: state bit i = pat[:i+1] matches here
    m = len(pat)
    mask = {}
    for i, c in enumerate(pat):
        mask[c] = mask.get(c, 0) | 1 << i
    D, hi, res = 0, 1 << (m - 1), []
    for i, c in enumerate(text):
        D = ((D << 1) | 1) & mask.get(c, 0)
        if D & hi:
            res.append(i - m + 1)
    return res


def bitset_match(text, pat):            # O(n*m/64) but only m big-int ops: fastest in Python
    pos = {}
    for i, c in enumerate(text):        # pos[c] bit i set iff text[i] == c
        pos[c] = pos.get(c, 0) | 1 << i
    occ = (1 << len(text)) - 1
    for j, c in enumerate(pat):
        if c != "?":                    # '?' in pattern = wildcard
            occ &= pos.get(c, 0) >> j   # bit i survives iff text[i+j] == pat[j] for all j
    occ &= (1 << max(0, len(text) - len(pat) + 1)) - 1
    return [i for i in range(len(text)) if occ >> i & 1]
# Dynamic text: keep pos[c] updated per point update.


# ==== 8.2 Wildcard matching with FFT (note) ====
# Wildcards (value 0) in text and/or pattern, letters 1..26. Position i matches iff
#   sum_j p_j * t_{i+j} * (p_j - t_{i+j})^2 = sum p^3 t - 2 p^2 t^2 + p t^3 = 0.
# Three convolutions of reversed p with t (see 09_algebra_poly.py: NTT convolution; use 2
# moduli or random letter values to avoid false zeros). Text-only/pattern-only wildcards: sum p*(p-t)^2.


# ==== 8.3 Trie (insert, count words with prefix) ====
class Trie:                             # O(len) per op; nodes are dicts
    def __init__(self):
        self.ch = [{}]
        self.pre = [0]                  # #words passing through node
        self.end = [0]                  # #words ending at node

    def insert(self, w, k=1):           # k=-1 to erase (word must exist)
        v = 0
        self.pre[0] += k
        for c in w:
            nx = self.ch[v].get(c)
            if nx is None:
                nx = self.ch[v][c] = len(self.pre)
                self.ch.append({})
                self.pre.append(0)
                self.end.append(0)
            v = nx
            self.pre[v] += k
        self.end[v] += k

    def count_prefix(self, w):          # #inserted words having prefix w
        v = 0
        for c in w:
            v = self.ch[v].get(c)
            if v is None:
                return 0
        return self.pre[v]


# ==== 8.3 Binary trie (max XOR with multiset of ints) ====
class BinaryTrie:                       # O(B) per op. Flat child array: ch[2v], ch[2v+1]
    def __init__(self, B=30):
        self.B = B
        self.ch = [0, 0]
        self.cnt = [0]

    def add(self, x, k=1):              # k=-1 removes one copy
        v = 0
        self.cnt[0] += k
        ch, cnt = self.ch, self.cnt
        for b in range(self.B - 1, -1, -1):
            d = x >> b & 1
            if not ch[2 * v + d]:
                ch[2 * v + d] = len(cnt)
                ch += (0, 0)
                cnt.append(0)
            v = ch[2 * v + d]
            cnt[v] += k

    def max_xor(self, x):               # max x ^ y over stored y (set must be non-empty)
        v, r = 0, 0
        ch, cnt = self.ch, self.cnt
        for b in range(self.B - 1, -1, -1):
            d = (x >> b & 1) ^ 1
            u = ch[2 * v + d]
            if u and cnt[u] > 0:
                r |= 1 << b
                v = u
            else:
                v = ch[2 * v + (d ^ 1)]
        return r
# min_xor: prefer the same bit. k-th smallest / count of y with x^y < K: walk using cnt.


# ==== 8.3 Suffix array: SA-IS O(n) (ACL port) ====
def sa_is(s, upper):                    # s: list of ints in [0, upper]. Recursion depth O(log n).
    n = len(s)
    if n == 0:
        return []
    if n == 1:
        return [0]
    if n == 2:
        return [0, 1] if s[0] < s[1] else [1, 0]
    sa = [0] * n
    ls = [False] * n
    for i in range(n - 2, -1, -1):
        ls[i] = ls[i + 1] if s[i] == s[i + 1] else s[i] < s[i + 1]
    sum_l = [0] * (upper + 1)
    sum_s = [0] * (upper + 1)
    for i in range(n):
        if not ls[i]:
            sum_s[s[i]] += 1
        else:
            sum_l[s[i] + 1] += 1
    for i in range(upper + 1):
        sum_s[i] += sum_l[i]
        if i < upper:
            sum_l[i + 1] += sum_s[i]

    def induce(lms):
        for i in range(n):
            sa[i] = -1
        buf = sum_s[:]
        for d in lms:
            if d != n:
                sa[buf[s[d]]] = d
                buf[s[d]] += 1
        buf = sum_l[:]
        sa[buf[s[n - 1]]] = n - 1
        buf[s[n - 1]] += 1
        for i in range(n):
            v = sa[i]
            if v >= 1 and not ls[v - 1]:
                sa[buf[s[v - 1]]] = v - 1
                buf[s[v - 1]] += 1
        buf = sum_l[:]
        for i in range(n - 1, -1, -1):
            v = sa[i]
            if v >= 1 and ls[v - 1]:
                buf[s[v - 1] + 1] -= 1
                sa[buf[s[v - 1] + 1]] = v - 1

    lms_map = [-1] * (n + 1)
    lms = []
    for i in range(1, n):
        if not ls[i - 1] and ls[i]:
            lms_map[i] = len(lms)
            lms.append(i)
    m = len(lms)
    induce(lms)
    if m:
        sorted_lms = [v for v in sa if lms_map[v] != -1]
        rec_s = [0] * m
        rec_upper = 0
        for i in range(1, m):
            l, r = sorted_lms[i - 1], sorted_lms[i]
            end_l = lms[lms_map[l] + 1] if lms_map[l] + 1 < m else n
            end_r = lms[lms_map[r] + 1] if lms_map[r] + 1 < m else n
            same = True
            if end_l - l != end_r - r:
                same = False
            else:
                while l < end_l:
                    if s[l] != s[r]:
                        break
                    l += 1
                    r += 1
                if l == n or s[l] != s[r]:
                    same = False
            if not same:
                rec_upper += 1
            rec_s[lms_map[sorted_lms[i]]] = rec_upper
        rec_sa = sa_is(rec_s, rec_upper)
        for i in range(m):
            sorted_lms[i] = lms[rec_sa[i]]
        induce(sorted_lms)
    return sa


def suffix_array(s):                    # str or list of ints -> SA (start indices, sorted)
    a = [ord(c) for c in s] if isinstance(s, str) else list(s)
    if not a:
        return []
    comp = {v: i for i, v in enumerate(sorted(set(a)))}
    return sa_is([comp[v] for v in a], len(comp) - 1)


# ==== 8.3 Suffix array: prefix doubling O(n log^2 n) with sort (short, memorable) ====
def suffix_array_doubling(s):
    n = len(s)
    comp = {v: i for i, v in enumerate(sorted(set(s)))}
    rank = [comp[c] for c in s]
    sa = list(range(n))
    k = 1
    while n > 1:
        key = [rank[i] * (n + 1) + (rank[i + k] + 1 if i + k < n else 0) for i in range(n)]
        sa.sort(key=key.__getitem__)
        new = [0] * n
        for j in range(1, n):
            new[sa[j]] = new[sa[j - 1]] + (key[sa[j]] != key[sa[j - 1]])
        rank = new
        if rank[sa[-1]] == n - 1:
            break
        k <<= 1
    return sa
# SA-IS (conceptual): classify suffixes S/L, sort LMS substrings by induced sorting, recurse on
# their names if not unique, then induce L then S suffixes from sorted LMS suffixes. O(n).


# ==== 8.3 LCP array (Kasai) ====
def lcp_kasai(s, sa):                   # O(n). lcp[i] = lcp(s[sa[i]:], s[sa[i+1]:]), len n-1
    n = len(s)
    rank = [0] * n
    for i, v in enumerate(sa):
        rank[v] = i
    lcp = [0] * max(0, n - 1)
    h = 0
    for i in range(n):
        if h:
            h -= 1
        r = rank[i]
        if r == n - 1:
            h = 0
            continue
        j = sa[r + 1]
        while i + h < n and j + h < n and s[i + h] == s[j + h]:
            h += 1
        lcp[r] = h
    return lcp
# lcp of arbitrary suffixes i, j = min(lcp[rank[i]..rank[j]-1]) -> sparse table (02).
# Pattern search: binary search on SA comparing s[sa[mid]:sa[mid]+m] (O(m log n)).


# ==== 8.3 Suffix automaton (distinct substrings, occurrences, LCS) ====
# Build O(n) states/transitions (dict transitions). <= 2n-1 states.
# cnt[v] = |endpos(v)| = occurrences of every string in state v.
class SuffixAutomaton:
    def __init__(self, s):
        nxt, link, ln, cnt = [{}], [-1], [0], [0]
        last = 0
        for ch in s:
            cur = len(ln)
            nxt.append({})
            ln.append(ln[last] + 1)
            link.append(0)
            cnt.append(1)
            p = last
            while p != -1 and ch not in nxt[p]:
                nxt[p][ch] = cur
                p = link[p]
            if p != -1:
                q = nxt[p][ch]
                if ln[p] + 1 == ln[q]:
                    link[cur] = q
                else:
                    cl = len(ln)
                    nxt.append(nxt[q].copy())
                    ln.append(ln[p] + 1)
                    link.append(link[q])
                    cnt.append(0)
                    while p != -1 and nxt[p].get(ch) == q:
                        nxt[p][ch] = cl
                        p = link[p]
                    link[q] = link[cur] = cl
            last = cur
        # states by decreasing len (counting sort) -> children before suffix-link parents
        buckets = [[] for _ in range(len(s) + 1)]
        for v in range(len(ln)):
            buckets[ln[v]].append(v)
        order = [v for b in reversed(buckets) for v in b]
        for v in order:
            if link[v] >= 0:
                cnt[link[v]] += cnt[v]
        self.nxt, self.link, self.ln, self.cnt, self.order, self.last = nxt, link, ln, cnt, order, last

    def distinct_substrings(self):      # non-empty
        ln, link = self.ln, self.link
        return sum(ln[v] - ln[link[v]] for v in range(1, len(ln)))

    def occurrences(self, p):           # #occurrences of non-empty p
        v = 0
        for ch in p:
            v = self.nxt[v].get(ch)
            if v is None:
                return 0
        return self.cnt[v]
# Also: k-th distinct substring (count paths dp over DAG in reverse order), first occurrence
# (store firstpos at creation; clone copies q's), shortest non-occurring string (BFS on DAG).


def lcs_multi(strs):                    # longest common substring length of all strs. O(sum |s|)
    S = SuffixAutomaton(strs[0])
    nxt, link, ln, order = S.nxt, S.link, S.ln, S.order
    best = ln[:]
    for t in strs[1:]:
        cur = [0] * len(ln)
        v = l = 0
        for ch in t:                    # standard matching walk: longest suffix of t[:i] in SAM
            while v and ch not in nxt[v]:
                v = link[v]
                l = ln[v]
            if ch in nxt[v]:
                v = nxt[v][ch]
                l += 1
            if l > cur[v]:
                cur[v] = l
        for v in order:                 # a match in v fully matches its link state
            if cur[v] and link[v] >= 0:
                cur[link[v]] = ln[link[v]]
        for v in range(len(ln)):
            if cur[v] < best[v]:
                best[v] = cur[v]
        # (cur[v] <= ln[v] always holds)
    return max(best)
# Two strings: just the walk above, tracking max l (and end position i for the substring).
# SA alternative: concat with unique separators, sliding window over SA covering all colors,
# answer = max over windows of min LCP inside.


# ==== 8.3 Suffix tree / Ukkonen (note) ====
# Ukkonen builds the suffix tree online in O(n) but is long and bug-prone. Use instead:
# suffix tree of s = suffix-link tree of SAM of reversed s; or SA + LCP (an LCP interval = node).


# ==== 8.3 Palindromic tree (eertree) ====
# O(n). Node 0: len -1 root, node 1: len 0 root. Distinct palindromes = len(ln) - 2.
# cnt[v] after build = occurrences of palindrome v in s.
class Eertree:
    def __init__(self, s):
        ln, link, nxt, cnt = [-1, 0], [0, 0], [{}, {}], [0, 0]
        last = 1
        for i, ch in enumerate(s):
            v = last
            while not (i - 1 - ln[v] >= 0 and s[i - 1 - ln[v]] == ch):
                v = link[v]
            if ch in nxt[v]:
                last = nxt[v][ch]
                cnt[last] += 1
                continue
            w = len(ln)
            ln.append(ln[v] + 2)
            nxt.append({})
            cnt.append(1)
            if ln[w] == 1:
                link.append(1)
            else:
                u = link[v]
                while not (i - 1 - ln[u] >= 0 and s[i - 1 - ln[u]] == ch):
                    u = link[u]
                link.append(nxt[u][ch])
            nxt[v][ch] = w
            last = w
        for v in range(len(ln) - 1, 1, -1):     # links point to shorter (earlier) nodes
            cnt[link[v]] += cnt[v]
        self.ln, self.link, self.nxt, self.cnt = ln, link, nxt, cnt
# Distinct palindromes so far after each prefix = number of nodes - 2 at that time.


# ==== 8.3 Manacher (odd / even radii) ====
def manacher(s):
    # O(n). d1[i]: #odd palindromes centered at i (longest = s[i-d1+1 : i+d1]).
    #       d2[i]: #even palindromes centered between i-1 and i (longest = s[i-d2 : i+d2]).
    n = len(s)
    d1 = [0] * n
    l, r = 0, -1
    for i in range(n):
        k = 1 if i > r else min(d1[l + r - i], r - i + 1)
        while i - k >= 0 and i + k < n and s[i - k] == s[i + k]:
            k += 1
        d1[i] = k
        if i + k - 1 > r:
            l, r = i - k + 1, i + k - 1
    d2 = [0] * n
    l, r = 0, -1
    for i in range(n):
        k = 0 if i > r else min(d2[l + r - i + 1], r - i + 1)
        while i - k - 1 >= 0 and i + k < n and s[i - k - 1] == s[i + k]:
            k += 1
        d2[i] = k
        if i + k - 1 > r:
            l, r = i - k, i + k - 1
    return d1, d2
# #palindromic substrings = sum(d1) + sum(d2). Longest: max(2*d1-1, 2*d2).


# ==== 8.3 Lyndon factorization (Duval) ====
def duval(s):                           # O(n) -> s = w1 w2 ... wk, Lyndon words, w1 >= w2 >= ...
    n = len(s)
    i = 0
    res = []
    while i < n:
        j, k = i + 1, i
        while j < n and s[k] <= s[j]:
            k = i if s[k] < s[j] else k + 1
            j += 1
        while i <= k:
            res.append(s[i:i + j - k])
            i += j - k
    return res


# ==== 8.3 Minimal string rotation ====
def min_rotation(s):                    # O(n) -> smallest i with s[i:]+s[:i] lexicographically min
    n = len(s)
    s2 = s + s
    i, j, k = 0, 1, 0
    while i < n and j < n and k < n:
        a, b = s2[i + k], s2[j + k]
        if a == b:
            k += 1
            continue
        if a > b:
            i += k + 1
        else:
            j += k + 1
        if i == j:
            j += 1
        k = 0
    return min(i, j)
# (Booth's algorithm / Duval on s+s give the same; this two-pointer form is shortest.)


# ==== 8.3 Run enumeration (conceptual) ====
# A run = maximal periodic substring s[l:r] with smallest period p, (r-l) >= 2p. There are < n runs
# (Runs theorem). Find them: for each Lyndon root w.r.t. both orders (<, >), candidate from the
# Lyndon array (computed with a stack + LCP/LCE queries via hashing or SA), extend left/right with
# LCE. Simpler O(n log n): Main-Lorentz divide & conquer with Z-functions. Rare in contests.


# ==== 8.4 String problems ====
def distinct_substrings_sa(s):          # O(n): n(n+1)/2 - sum(lcp)
    n = len(s)
    return n * (n + 1) // 2 - sum(lcp_kasai(s, suffix_array(s)))
# Via SAM: SuffixAutomaton(s).distinct_substrings().  Per prefix (online): add ln[cur]-ln[link[cur]].


def longest_repeated_substring(s):      # O(n): longest substring occurring >= 2 times (may overlap)
    if len(s) < 2:
        return ""
    sa = suffix_array(s)
    lcp = lcp_kasai(s, sa)
    i = max(range(len(lcp)), key=lcp.__getitem__)
    return s[sa[i]:sa[i] + lcp[i]]
# Occurring >= k times: max over windows of k-1 consecutive lcp values of their min.
# Non-overlapping: binary search length + hashing (track min/max start per hash).


def borders(s):                         # all border lengths (proper prefix == suffix), decreasing
    if not s:
        return []
    pi = prefix_function(s)
    res = []
    b = pi[-1]
    while b:
        res.append(b)
        b = pi[b - 1]
    return res


def periods(s):                         # all periods p (s[i] == s[i+p]), increasing; includes n
    n = len(s)
    return [n - b for b in borders(s)] + [n]
# Minimal period p = n - pi[-1]; s = t^k for some t iff n % p == 0. Count of each prefix as a
# substring: cnt[pi[i]]++ then propagate cnt[pi[k-1]] += cnt[k] for k decreasing (+1 each).


# ==================================== tests ====================================
if __name__ == "__main__":
    random.seed(2024)

    def rs(n, a="ab"):
        return "".join(random.choice(a) for _ in range(n))

    def occ_brute(t, p):
        return [i for i in range(len(t) - len(p) + 1) if t[i:i + len(p)] == p]

    for it in range(400):
        alpha = "ab" if it % 3 else "abc"
        t = rs(random.randint(0, 30), alpha)
        p = rs(random.randint(1, 4), alpha)
        exp = occ_brute(t, p)
        assert kmp_search(t, p) == exp == rabin_karp(t, p) == shift_and(t, p) == bitset_match(t, p)
        z = z_function(p + "#" + t)
        assert [i for i in range(len(t)) if z[len(p) + 1 + i] >= len(p)] == exp
        pw = "".join(c if random.random() < 0.7 else "?" for c in p)
        assert bitset_match(t, pw) == [i for i in range(len(t) - len(p) + 1)
                                       if all(c in "?" + t[i + j] for j, c in enumerate(pw))]
        n = len(t)
        pi = prefix_function(t)
        for i in range(n):
            assert pi[i] == max(k for k in range(i + 1) if t[:k] == t[i + 1 - k:i + 1])
        z = z_function(t)
        for i in range(n):
            k = 0
            while i + k < n and t[k] == t[i + k]:
                k += 1
            assert z[i] == k
        # automaton = KMP state after each char
        aut = kmp_automaton(p, alpha)
        st = 0
        for i, ch in enumerate(t):
            st = aut[st][alpha.index(ch)]
            assert st == max(k for k in range(len(p) + 1) if k <= i + 1 and t[i + 1 - k:i + 1] == p[:k])
        # rolling hash
        H = RollingHash(t)
        G = RollingHash(p)
        for _ in range(10):
            if n == 0:
                break
            l1, r1 = sorted(random.sample(range(n + 1), 2))
            l2, r2 = sorted((random.randint(0, len(p)), random.randint(0, len(p))))
            a, b = t[l1:r1], p[l2:r2]
            assert hash_compare(H, l1, r1, G, l2, r2) == (a > b) - (a < b)
            assert (H.get(l1, r1) == G.get(l2, r2)) == (a == b)
            k = 0
            while l1 + k < n and l2 + k < len(p) and t[l1 + k] == p[l2 + k]:
                k += 1
            assert hash_lcp(H, l1, G, l2) == k
            m = random.randint(l1, r1)
            assert RollingHash.combine(H.get(l1, m), H.get(m, r1), r1 - m) == H.get(l1, r1)
        # suffix arrays + LCP
        sa_b = sorted(range(n), key=lambda i: t[i:])
        assert suffix_array(t) == sa_b == suffix_array_doubling(t)
        arr = [random.randrange(4) for _ in range(n)]
        assert suffix_array(arr) == sorted(range(n), key=lambda i: arr[i:])
        lcp = lcp_kasai(t, sa_b)
        for i in range(n - 1):
            a, b = t[sa_b[i]:], t[sa_b[i + 1]:]
            k = 0
            while k < min(len(a), len(b)) and a[k] == b[k]:
                k += 1
            assert lcp[i] == k
        # distinct substrings: SA vs SAM vs brute; occurrences; LRS
        subs = {t[i:j] for i in range(n) for j in range(i + 1, n + 1)}
        S = SuffixAutomaton(t)
        assert distinct_substrings_sa(t) == len(subs) == S.distinct_substrings()
        for q in list(subs)[:10] + [p]:
            assert S.occurrences(q) == len(occ_brute(t, q))
        rep = [x for x in subs if len(occ_brute(t, x)) >= 2]
        assert len(longest_repeated_substring(t)) == max(map(len, rep), default=0)
        assert len(occ_brute(t, longest_repeated_substring(t))) >= 2 or not rep
        # LCS of several strings
        ss = [rs(random.randint(1, 12), alpha) for _ in range(random.randint(1, 4))]
        cs = {ss[0][i:j] for i in range(len(ss[0])) for j in range(i + 1, len(ss[0]) + 1)}
        cs = [x for x in cs if all(x in y for y in ss)]
        assert lcs_multi(ss) == max(map(len, cs), default=0)
        # palindromes
        pals = {x for x in subs if x == x[::-1]}
        E = Eertree(t)
        assert len(E.ln) - 2 == len(pals)
        d1, d2 = manacher(t)
        for i in range(n):
            k = 0
            while i - k >= 0 and i + k < n and t[i - k] == t[i + k]:
                k += 1
            assert d1[i] == k
            k = 0
            while i - k - 1 >= 0 and i + k < n and t[i - k - 1] == t[i + k]:
                k += 1
            assert d2[i] == k
        assert sum(d1) + sum(d2) == sum(len(occ_brute(t, x)) for x in pals)
        assert sum(E.cnt[2:]) == sum(d1) + sum(d2)
        # Lyndon / rotation / borders / periods
        f = duval(t)
        assert "".join(f) == t and all(f[i] >= f[i + 1] for i in range(len(f) - 1))
        for w in f:
            assert all(w < w[i:] + w[:i] for i in range(1, len(w)))
        if n:
            rots = [t[i:] + t[:i] for i in range(n)]
            assert min_rotation(t) == rots.index(min(rots))
        assert borders(t) == [k for k in range(n - 1, 0, -1) if t[:k] == t[n - k:]]
        assert periods(t) == [q for q in range(1, n + 1) if all(t[i] == t[i + q] for i in range(n - q))] or n == 0

    # Aho-Corasick
    for it in range(200):
        pats = [rs(random.randint(1, 4), "abc") for _ in range(random.randint(1, 8))]
        t = rs(random.randint(0, 40), "abc")
        assert AhoCorasick(pats).count(t) == [len(occ_brute(t, p)) for p in pats]
    # tries
    T = Trie()
    words = []
    BT = BinaryTrie(10)
    nums = []
    for it in range(2000):
        r = random.random()
        if r < 0.5 or not words:
            w = rs(random.randint(0, 5), "abc")
            T.insert(w)
            words.append(w)
        elif r < 0.7:
            w = words.pop(random.randrange(len(words)))
            T.insert(w, -1)
        q = rs(random.randint(0, 3), "abc")
        assert T.count_prefix(q) == sum(w.startswith(q) for w in words)
        if random.random() < 0.6 or not nums:
            x = random.randrange(1024)
            BT.add(x)
            nums.append(x)
        else:
            BT.add(nums.pop(random.randrange(len(nums))), -1)
        if nums:
            x = random.randrange(1024)
            assert BT.max_xor(x) == max(x ^ y for y in nums)
    # eertree occurrence counts
    for it in range(100):
        t = rs(random.randint(1, 25))
        E = Eertree(t)
        cntmap = {}
        for i in range(len(t)):
            for j in range(i + 1, len(t) + 1):
                x = t[i:j]
                if x == x[::-1]:
                    cntmap[x] = cntmap.get(x, 0) + 1
        assert sorted(E.cnt[2:]) == sorted(cntmap.values())
        assert sorted(E.ln[2:]) == sorted(map(len, cntmap))
    # larger sanity / speed
    t = rs(100000, "abcd")
    sa = suffix_array(t)
    assert sa == suffix_array_doubling(t)
    assert distinct_substrings_sa(t) == SuffixAutomaton(t).distinct_substrings()
    t = "a" * 50000
    assert suffix_array(t) == list(range(49999, -1, -1))
    assert len(Eertree(t).ln) - 2 == 50000
    print("ALL OK")
