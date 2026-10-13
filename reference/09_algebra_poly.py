# 09_algebra_poly.py -- Linear algebra (9.4) and polynomials / transforms (9.5)
# Python 3 stdlib only. All poly code goes through ONE convolve().
from operator import mul
import decimal

MOD = 998244353

# ==== 9.4 Linear Algebra: matrix multiply / power ====
# O(n^3). C = A*B mod p. Transpose B once, inner loop = C-level sum(map(mul)), one % per entry.
# PyPy gotcha: the un-reduced sum is a bigint (> 2^63) -> slow on PyPy; there use a loop with % per term.
def mat_mul(A, B, mod):
    Bt = list(zip(*B))
    return [[sum(map(mul, row, col)) % mod for col in Bt] for row in A]

# O(n^3 log e). A^e mod p.
def mat_pow(A, e, mod):
    n = len(A)
    R = [[int(i == j) for j in range(n)] for i in range(n)]
    while e:
        if e & 1:
            R = mat_mul(R, A, mod)
        A = mat_mul(A, A, mod)
        e >>= 1
    return R

# O(d^3 log k). a[n] = sum_{i<d} c[i]*a[n-1-i]; a = first d terms. Returns a[k].
# (Bostan-Mori below is O(d log d log k) -- prefer it for big d.)  Path counting: (adj^k)[u][v].
def linrec_matpow(c, a, k, mod):
    d = len(c)
    if k < d:
        return a[k] % mod
    M = [list(c)] + [[int(j == i) for j in range(d)] for i in range(d - 1)]
    P = mat_pow(M, k - d + 1, mod)
    return sum(P[0][j] * a[d - 1 - j] for j in range(d)) % mod

# ==== 9.4 Gaussian elimination / systems of linear equations ====
# O(n*m*min(n,m)). Over reals, partial pivoting. Returns (status, x): 0 = no solution (x None),
# 1 = unique, 2 = infinitely many (x = one solution, free vars = 0). Gotcha: tune eps to input scale.
def gauss_real(A, b, eps=1e-9):
    n, m = len(A), len(A[0])
    M = [[float(v) for v in A[i]] + [float(b[i])] for i in range(n)]
    where = [-1] * m
    r = 0
    for c in range(m):
        if r == n:
            break
        p = max(range(r, n), key=lambda i: abs(M[i][c]))
        if abs(M[p][c]) < eps:
            continue
        M[r], M[p] = M[p], M[r]
        iv = 1.0 / M[r][c]
        Mr = M[r] = [v * iv for v in M[r]]
        for i in range(n):
            f = M[i][c]
            if i != r and f != 0.0:
                M[i] = [x - f * y for x, y in zip(M[i], Mr)]
        where[c] = r
        r += 1
    for i in range(r, n):
        if abs(M[i][m]) > eps:
            return 0, None
    x = [M[where[c]][m] if where[c] >= 0 else 0.0 for c in range(m)]
    return (1 if r == m else 2), x

# O(n*m*min(n,m)). Gauss-Jordan mod prime p on matrix M (n x ncols+extra), only first ncols pivoted.
# Returns (M in RREF, where[c] = pivot row of column c or -1, rank).
def rref_mod(M, ncols, p):
    M = [[v % p for v in row] for row in M]
    n = len(M)
    where = [-1] * ncols
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, n) if M[i][c]), -1)
        if piv < 0:
            continue
        M[r], M[piv] = M[piv], M[r]
        iv = pow(M[r][c], p - 2, p)
        Mr = M[r] = [v * iv % p for v in M[r]]
        for i in range(n):
            f = M[i][c]
            if i != r and f:
                M[i] = [(x - f * y) % p for x, y in zip(M[i], Mr)]
        where[c] = r
        r += 1
        if r == n:
            break
    return M, where, r

# Solve A x = b mod prime p. Returns (status, x) like gauss_real.
def solve_mod(A, b, p):
    m = len(A[0])
    M, where, r = rref_mod([A[i] + [b[i]] for i in range(len(A))], m, p)
    if any(M[i][m] for i in range(r, len(M))):
        return 0, None
    return (1 if r == m else 2), [M[where[c]][m] if where[c] >= 0 else 0 for c in range(m)]

def rank_mod(A, p):
    return rref_mod(A, len(A[0]), p)[2]

# Inverse mod prime p, or None if singular. O(n^3).
def mat_inv_mod(A, p):
    n = len(A)
    M, where, r = rref_mod([A[i] + [int(i == j) for j in range(n)] for i in range(n)], n, p)
    return None if r < n else [row[n:] for row in M]

# ==== 9.4 Determinants ====
# O(n^3). det mod prime p.
def det_mod(A, p):
    M = [[v % p for v in row] for row in A]
    n, det = len(M), 1
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c]), -1)
        if piv < 0:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det = det * M[c][c] % p
        iv = pow(M[c][c], p - 2, p)
        Mc = M[c]
        for i in range(c + 1, n):
            f = M[i][c] * iv % p
            if f:
                M[i] = [(x - f * y) % p for x, y in zip(M[i], Mc)]
    return det % p

# O(n^3 log m). det mod ANY m (composite ok): Euclid-style row reduction, no inverses needed.
def det_any_mod(A, m):
    M = [[v % m for v in row] for row in A]
    n, det = len(M), 1
    for c in range(n):
        for r in range(c + 1, n):
            while M[r][c]:
                q = M[c][c] // M[r][c]
                if q:
                    M[c] = [(x - q * y) % m for x, y in zip(M[c], M[r])]
                M[c], M[r] = M[r], M[c]
                det = -det
        if M[c][c] == 0:
            return 0
        det = det * M[c][c] % m
    return det % m

# ==== 9.4 Gaussian elimination over GF(2) ====
# O(n*m*n/64). rows: ints, bit c (c<m) = coeff of x_c, bit m = rhs. Returns (rank, x bitmask or None).
# For plain rank pass m = number of columns and rhs bits 0.
def gauss_gf2(rows, m):
    rows = list(rows)
    n = len(rows)
    where = [-1] * m
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if rows[i] >> c & 1), -1)
        if piv < 0:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pr = rows[r]
        for i in range(n):
            if i != r and rows[i] >> c & 1:
                rows[i] ^= pr
        where[c] = r
        r += 1
    if any(rows[i] >> m & 1 for i in range(r, n)):
        return r, None
    x = 0
    for c in range(m):
        if where[c] >= 0 and rows[where[c]] >> m & 1:
            x |= 1 << c
    return r, x

# ==== 9.4 XOR basis (linear basis over GF(2)) ====
# insert/can/max O(B). kth(k): k-th smallest (0-indexed, k=0 -> 0) DISTINCT value in span, -1 if none.
# Gotcha: if the multiset has len(vals) > rank, 0 is also a nonempty-subset xor (adjust k for "nonempty").
class XorBasis:
    def __init__(self, bits=62):
        self.B = bits
        self.b = [0] * bits
        self.rank = 0

    def insert(self, x):
        b = self.b
        for i in range(self.B - 1, -1, -1):
            if x >> i & 1:
                if not b[i]:
                    b[i] = x
                    self.rank += 1
                    return True
                x ^= b[i]
        return False

    def can(self, x):
        for i in range(self.B - 1, -1, -1):
            if x >> i & 1:
                x ^= self.b[i]
        return x == 0

    def max_xor(self, x=0):
        for i in range(self.B - 1, -1, -1):
            if x ^ self.b[i] > x:
                x ^= self.b[i]
        return x

    def kth(self, k):
        b = self.b[:]
        for i in range(self.B):
            for j in range(i - 1, -1, -1):
                if b[j] and b[i] >> j & 1:
                    b[i] ^= b[j]
        basis = [v for v in b if v]
        if k >= 1 << len(basis):
            return -1
        res = 0
        for j, v in enumerate(basis):
            if k >> j & 1:
                res ^= v
        return res

# ==== 9.4 Berlekamp-Massey ====
# O(n^2). Shortest c with s[i] = sum_j c[j]*s[i-1-j] (mod prime). Needs 2*L terms to be sure.
def berlekamp_massey(s, mod=MOD):
    C, B = [1], [1]
    L, m, b = 0, 1, 1
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1):
            d = (d + C[i] * s[n - i]) % mod
        if d == 0:
            m += 1
            continue
        T = C[:]
        coef = d * pow(b, mod - 2, mod) % mod
        if len(C) < len(B) + m:
            C += [0] * (len(B) + m - len(C))
        for i in range(len(B)):
            C[i + m] = (C[i + m] - coef * B[i]) % mod
        if 2 * L <= n:
            L, B, b, m = n + 1 - L, T, d, 1
        else:
            m += 1
    C += [0] * (L + 1 - len(C))
    return [(-c) % mod for c in C[1:L + 1]]

# k-th term of linear recurrence: see linrec_kth (Bostan-Mori) in 9.5. Kitamasa (x^k mod charpoly)
# has the same complexity; Bostan-Mori is shorter.


# ==== 9.5 Polynomial multiplication: Kronecker substitution (any modulus) ====
# convolve(a, b, mod=None): exact. Pack coeffs into one huge number with w digits per slot (w large
# enough that no slot overflows), multiply once, unpack. Works for ANY mod (or exact ints, negatives ok).
# BENCHMARK (CPython 3.14, mod 998244353, len n each):
#   int (Karatsuba)          n=1e4: 0.16s  n=1e5: 5.9s  n=2e5: 17s   <- CPython int mul is NOT FFT!
#   decimal (libmpdec NTT)   n=1e4: 0.06s  n=1e5: 0.66s n=2e5: 1.2s  <- used for big inputs
#   pure-Python ntt_convolve n=1e5: 5.9s.   Built on convolve: poly_inv n=1e5 2.2s, poly_exp n=1e5 5.8s.
# So: decimal-based Kronecker for big sizes, int-based for small. On PyPy use ntt_convolve (decimal
# is slow there). Alternatives for arbitrary mod with NTT: 3 NTT primes + CRT, or split-FFT (15-bit).
_DCTX = decimal.Context(prec=decimal.MAX_PREC, Emax=decimal.MAX_EMAX, Emin=decimal.MIN_EMIN)

def _kron(a, b):  # nonneg ints
    n = len(a) + len(b) - 1
    ma, mb = max(a), max(b)
    if ma == 0 or mb == 0:
        return [0] * n
    bound = ma * mb * min(len(a), len(b))
    if min(len(a), len(b)) <= 300:  # int path (bytes packing)
        W = (bound.bit_length() + 7) >> 3
        fa = int.from_bytes(b''.join(x.to_bytes(W, 'little') for x in a), 'little')
        fb = fa if a is b else int.from_bytes(b''.join(x.to_bytes(W, 'little') for x in b), 'little')
        B = (fa * fb).to_bytes(n * W, 'little')
        return [int.from_bytes(B[i:i + W], 'little') for i in range(0, n * W, W)]
    w = len(str(bound))
    f = '%0' + str(w) + 'd'
    fa = _DCTX.create_decimal(''.join([f % x for x in reversed(a)]))
    fb = fa if a is b else _DCTX.create_decimal(''.join([f % x for x in reversed(b)]))
    s = str(_DCTX.multiply(fa, fb)).zfill(n * w)
    return [int(s[i:i + w]) for i in range(n * w - w, -1, -w)]

def convolve(a, b, mod=None):
    if not a or not b:
        return []
    if mod is not None:
        if len(a) <= 32 or len(b) <= 32:  # naive is faster for tiny inputs
            res = [0] * (len(a) + len(b) - 1)
            for i, x in enumerate(a):
                if x:
                    for j, y in enumerate(b):
                        res[i + j] += x * y
            return [v % mod for v in res]
        sq = a is b
        a = [x % mod for x in a]
        b = a if sq else [x % mod for x in b]
        return [v % mod for v in _kron(a, b)]
    if min(a) >= 0 and min(b) >= 0:
        return _kron(a, b)
    ap, an = [max(x, 0) for x in a], [max(-x, 0) for x in a]
    bp, bn = [max(x, 0) for x in b], [max(-x, 0) for x in b]
    def k(u, v):
        return _kron(u, v) if any(u) and any(v) else [0] * (len(a) + len(b) - 1)
    return [w + x - y - z for w, x, y, z in zip(k(ap, bp), k(an, bn), k(ap, bn), k(an, bp))]

# ==== 9.5 NTT (mod 998244353, primitive root 3) ====
# O(n log n). Iterative in-place; len(a) power of 2. Reference / PyPy (in CPython convolve is faster).
def ntt(a, invert=False):
    n = len(a)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        w = pow(3, (MOD - 1) // length, MOD)
        if invert:
            w = pow(w, MOD - 2, MOD)
        half = length >> 1
        ws = [1] * half
        for i in range(1, half):
            ws[i] = ws[i - 1] * w % MOD
        for s in range(0, n, length):
            for k in range(half):
                u = a[s + k]
                v = a[s + k + half] * ws[k] % MOD
                a[s + k] = (u + v) % MOD
                a[s + k + half] = (u - v) % MOD
        length <<= 1
    if invert:
        ni = pow(n, MOD - 2, MOD)
        for i in range(n):
            a[i] = a[i] * ni % MOD

def ntt_convolve(a, b):
    n = len(a) + len(b) - 1
    sz = 1 << (n - 1).bit_length()
    fa = [x % MOD for x in a] + [0] * (sz - len(a))
    fb = [x % MOD for x in b] + [0] * (sz - len(b))
    ntt(fa)
    ntt(fb)
    fa = [x * y % MOD for x, y in zip(fa, fb)]
    ntt(fa, True)
    return fa[:n]

# ==== 9.5 FFT (complex, cmath) ====
# O(n log n) floats. Exact only while result coeffs < ~1e14 (double precision); rarely needed in
# Python since convolve() is exact. Kept short for reference.
import cmath

def fft(a, invert=False):
    n = len(a)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        ang = (-2 if invert else 2) * cmath.pi / length
        half = length >> 1
        ws = [cmath.exp(1j * ang * k) for k in range(half)]
        for s in range(0, n, length):
            for k in range(half):
                u, v = a[s + k], a[s + k + half] * ws[k]
                a[s + k], a[s + k + half] = u + v, u - v
        length <<= 1
    if invert:
        for i in range(n):
            a[i] /= n

def fft_multiply(a, b):
    n = len(a) + len(b) - 1
    sz = 1 << (n - 1).bit_length()
    fa = [complex(x) for x in a] + [0j] * (sz - len(a))
    fb = [complex(x) for x in b] + [0j] * (sz - len(b))
    fft(fa)
    fft(fb)
    fa = [x * y for x, y in zip(fa, fb)]
    fft(fa, True)
    return [round(x.real) for x in fa[:n]]

# ==== 9.5 Power series mod 998244353: inverse, log, exp, sqrt, pow (Newton) ====
# All take f (list) and n = number of output terms; O(n log n)-ish (convolve-bound).
def _conv(a, b):
    return convolve(a, b, MOD)

def poly_inv(f, n):  # needs f[0] != 0. g = g*(2 - f*g)
    g = [pow(f[0], MOD - 2, MOD)]
    k = 1
    while k < n:
        k <<= 1
        h = _conv(f[:k], g)[:k]
        h = [-x % MOD for x in h]
        h[0] = (h[0] + 2) % MOD
        g = _conv(g, h)[:k]
    return g[:n]

def poly_deriv(f):
    return [f[i] * i % MOD for i in range(1, len(f))]

def poly_integ(f):
    n = len(f) + 1
    inv = [0, 1] + [0] * (n - 1)
    for i in range(2, n + 1):
        inv[i] = (MOD - (MOD // i) * inv[MOD % i] % MOD) % MOD
    return [0] + [f[i] * inv[i + 1] % MOD for i in range(len(f))]

def poly_log(f, n):  # needs f[0] == 1
    if n == 1:
        return [0]
    r = _conv(poly_deriv(f[:n]), poly_inv(f, n))[:n - 1]
    r += [0] * (n - 1 - len(r))
    return poly_integ(r)

def poly_exp(f, n):  # needs f[0] == 0. g = g*(1 - log g + f)
    g = [1]
    k = 1
    while k < n:
        k <<= 1
        l = poly_log(g, k)
        h = [((f[i] if i < len(f) else 0) - l[i]) % MOD for i in range(k)]
        h[0] = (h[0] + 1) % MOD
        g = _conv(g, h)[:k]
    return g[:n]

def poly_sqrt(f, n):  # needs f[0] == 1 (else: sqrt of f[0] via Tonelli, see 08; shift even zeros)
    g = [1]
    inv2 = (MOD + 1) // 2
    k = 1
    while k < n:
        k <<= 1
        h = _conv(f[:k], poly_inv(g, k))[:k]
        h += [0] * (k - len(h))
        g = [(h[i] + (g[i] if i < len(g) else 0)) * inv2 % MOD for i in range(k)]
    return g[:n]

def poly_pow(f, e, n):  # f^e mod x^n, any f, e >= 0 (e may be huge)
    if e == 0:
        return [1] + [0] * (n - 1)
    t = next((i for i, c in enumerate(f) if c % MOD), -1)
    if t < 0 or t * e >= n:
        return [0] * n
    c = f[t] % MOD
    ic = pow(c, MOD - 2, MOD)
    g = [x * ic % MOD for x in f[t:t + n - t * e]]
    m = n - t * e
    l = poly_log(g, m)
    em = e % MOD
    g = poly_exp([x * em % MOD for x in l], m)
    ce = pow(c, e, MOD)
    return [0] * (t * e) + [x * ce % MOD for x in g]

# ==== 9.5 Polynomial division / mod ====
# O(n log n). a = b*q + r, deg r < deg b. Needs b[-1] != 0. Returns (q, r) with len(r) = len(b)-1.
def poly_divmod(a, b):
    n, m = len(a), len(b)
    if n < m:
        return [0], a + [0] * (m - 1 - n)
    k = n - m + 1
    if m <= 64 or k <= 64:  # naive long division
        r = [x % MOD for x in a]
        il = pow(b[-1], MOD - 2, MOD)
        q = [0] * k
        for i in range(k - 1, -1, -1):
            c = q[i] = r[i + m - 1] * il % MOD
            if c:
                for j in range(m):
                    r[i + j] = (r[i + j] - c * b[j]) % MOD
        return q, r[:m - 1]
    q = _conv(a[::-1][:k], poly_inv(b[::-1], k))[:k][::-1]
    bq = _conv(b, q)
    return q, [(a[i] - bq[i]) % MOD for i in range(m - 1)]

def poly_mod(a, b):
    return a if len(a) < len(b) else poly_divmod(a, b)[1]

# ==== 9.5 Multipoint evaluation (subproduct tree) ====
# O(n log^2 n). Values f(x) for x in xs. Small subtrees fall back to Horner.
def multipoint_eval(f, xs):
    m = len(xs)
    if m == 0:
        return []
    sz = 1
    while sz < m:
        sz <<= 1
    T = [[1]] * (2 * sz)
    for i, x in enumerate(xs):
        T[sz + i] = [-x % MOD, 1]
    for i in range(sz - 1, 0, -1):
        T[i] = _conv(T[2 * i], T[2 * i + 1])
    res = [0] * m
    stack = [(1, poly_mod(f, T[1]), 0, sz)]
    while stack:
        v, r, lo, hi = stack.pop()
        if lo >= m:
            continue
        if hi - lo <= 32 or len(r) <= 32:
            for i in range(lo, min(hi, m)):
                x, acc = xs[i], 0
                for c in reversed(r):
                    acc = (acc * x + c) % MOD
                res[i] = acc
            continue
        mid = (lo + hi) >> 1
        stack.append((2 * v, poly_mod(r, T[2 * v]), lo, mid))
        stack.append((2 * v + 1, poly_mod(r, T[2 * v + 1]), mid, hi))
    return res
# (Fast interpolation O(n log^2 n) = same tree + derivative trick; rarely needed, O(n^2) below.)

# ==== 9.5 Lagrange interpolation ====
# O(n^2). Coefficients of the unique poly of deg < n through (xs[i], ys[i]); xs distinct, mod prime.
def interpolate(xs, ys, mod=MOD):
    n = len(xs)
    P = [1]  # prod (X - xs[i])
    for x in xs:
        P = [(u - x * v) % mod for u, v in zip([0] + P, P + [0])]
    res = [0] * n
    for i in range(n):
        xi, d = xs[i], 1
        for j in range(n):
            if j != i:
                d = d * (xi - xs[j]) % mod
        coef = ys[i] * pow(d, mod - 2, mod) % mod
        acc = 0
        for k in range(n, 0, -1):  # synthetic division P / (X - xi)
            acc = (P[k] + acc * xi) % mod
            res[k - 1] += coef * acc
    return [v % mod for v in res]

# O(n). ys = f(0), ..., f(n) (deg f <= n); returns f(x) mod prime (prime > n). E.g. sum of k-th powers.
def lagrange_consecutive(ys, x, mod=MOD):
    n = len(ys) - 1
    x %= mod
    pre = [1] * (n + 2)
    suf = [1] * (n + 2)
    for i in range(n + 1):
        pre[i + 1] = pre[i] * (x - i) % mod
    for i in range(n, -1, -1):
        suf[i] = suf[i + 1] * (x - i) % mod
    fact = [1] * (n + 1)
    for i in range(1, n + 1):
        fact[i] = fact[i - 1] * i % mod
    ifact = [1] * (n + 1)
    ifact[n] = pow(fact[n], mod - 2, mod)
    for i in range(n, 0, -1):
        ifact[i - 1] = ifact[i] * i % mod
    res = 0
    for i in range(n + 1):
        t = ys[i] * pre[i] % mod * suf[i + 1] % mod * ifact[i] % mod * ifact[n - i] % mod
        res += -t if (n - i) & 1 else t
    return res % mod

# ==== 9.5 Walsh-Hadamard transform (XOR convolution), AND / OR convolutions ====
# O(n log n), len power of 2, in place. Slice+zip butterflies (fast in CPython).
def fwht_xor(a, invert=False, mod=MOD):
    n, h = len(a), 1
    while h < n:
        for i in range(0, n, 2 * h):
            L, R = a[i:i + h], a[i + h:i + 2 * h]
            a[i:i + h] = [(x + y) % mod for x, y in zip(L, R)]
            a[i + h:i + 2 * h] = [(x - y) % mod for x, y in zip(L, R)]
        h <<= 1
    if invert:
        ni = pow(n, mod - 2, mod)
        a[:] = [x * ni % mod for x in a]

def fwht_or(a, invert=False, mod=MOD):  # zeta over subsets: A[S] = sum_{T subset S} a[T]
    n, h, s = len(a), 1, (-1 if invert else 1)
    while h < n:
        for i in range(0, n, 2 * h):
            a[i + h:i + 2 * h] = [(y + s * x) % mod for x, y in zip(a[i:i + h], a[i + h:i + 2 * h])]
        h <<= 1

def fwht_and(a, invert=False, mod=MOD):  # zeta over supersets
    n, h, s = len(a), 1, (-1 if invert else 1)
    while h < n:
        for i in range(0, n, 2 * h):
            a[i:i + h] = [(x + s * y) % mod for x, y in zip(a[i:i + h], a[i + h:i + 2 * h])]
        h <<= 1

def bit_convolve(a, b, kind, mod=MOD):  # kind in 'xor','or','and'; c[i op j] += a[i]*b[j]
    tr = {'xor': fwht_xor, 'or': fwht_or, 'and': fwht_and}[kind]
    a, b = a[:], b[:]
    tr(a, mod=mod)
    tr(b, mod=mod)
    c = [x * y % mod for x, y in zip(a, b)]
    tr(c, True, mod=mod)
    return c

# ==== 9.5 Subset convolution ====
# O(2^n n^2). c[S] = sum_{T subset S} a[T]*b[S\T]. Ranked zeta. CPython: n<=~15 ok; n=20 needs PyPy.
def subset_convolve(a, b, mod=MOD):
    N = len(a)
    n = N.bit_length() - 1
    pc = [0] * N
    for i in range(1, N):
        pc[i] = pc[i >> 1] + (i & 1)
    def ranked(f):
        F = [[0] * N for _ in range(n + 1)]
        for i in range(N):
            F[pc[i]][i] = f[i] % mod
        for row in F:
            fwht_or(row, mod=mod)
        return F
    A, B = ranked(a), ranked(b)
    c = [0] * N
    for k in range(n + 1):
        h = [0] * N
        for i in range(k + 1):
            h = [(z + x * y) for z, x, y in zip(h, A[i], B[k - i])]
        h = [z % mod for z in h]
        fwht_or(h, True, mod=mod)
        for s in range(N):
            if pc[s] == k:
                c[s] = h[s]
    return c

# ==== 9.5 Online convolution (CDQ / relaxed) ====
# O(n log^2 n). Computes f[0..n-1] where f[i] = calc(i, s_i), s_i = sum_{j<i} f[j]*g[i-j] (mod),
# g fully known in advance (g[0] unused). E.g. calc = lambda i, s: (h[i] + s) % mod.
# If g also depends on f (e.g. f = 1 + x f^2), use the fully relaxed variant (extra term for l==0 block).
def online_convolve(n, g, calc, mod=MOD):
    g = list(g[:n]) + [0] * (n - len(g))
    f = [0] * n
    acc = [0] * n
    stack = [(0, n, False)]
    while stack:  # iterative CDQ: (l, r, left_done)
        l, r, done = stack.pop()
        if r - l <= 32:
            for i in range(l, r):
                fi = f[i] = calc(i, acc[i] % mod) % mod
                if fi:
                    for k in range(i + 1, r):
                        acc[k] += fi * g[k - i]
            continue
        m = (l + r) >> 1
        if not done:
            stack.append((l, r, True))
            stack.append((l, m, False))
        else:
            c = convolve(f[l:m], g[1:r - l], mod)
            for i in range(m, r):
                acc[i] += c[i - l - 1]
            stack.append((m, r, False))
    return f

# ==== 9.5 Bostan-Mori: [x^k] P(x)/Q(x), and k-th term of linear recurrence ====
# O(d log d log k). Q[0] invertible. Any prime mod (convolve is any-mod).
def bostan_mori(P, Q, k, mod=MOD):
    P = [x % mod for x in P]
    Q = [x % mod for x in Q]
    while k:
        Qm = [c if i % 2 == 0 else (mod - c) % mod for i, c in enumerate(Q)]
        U = convolve(P, Qm, mod)
        V = convolve(Q, Qm, mod)
        P, Q = U[k & 1::2], V[::2]
        k >>= 1
    return P[0] * pow(Q[0], mod - 2, mod) % mod if P else 0

# a[n] = sum_{i<d} c[i]*a[n-1-i], a = first d terms. Returns a[k]. Pair with berlekamp_massey.
def linrec_kth(c, a, k, mod=MOD):
    d = len(c)
    if d == 0:
        return 0
    Q = [1] + [(-x) % mod for x in c]
    P = convolve(a[:d], Q, mod)[:d]
    return bostan_mori(P, Q, k, mod)


# ======================================================================================
if __name__ == "__main__":
    import random
    from fractions import Fraction
    from itertools import permutations
    rnd = random.Random(12345)
    R = lambda: rnd.randrange(MOD)

    def naive_conv(a, b, mod=None):
        if not a or not b:
            return []
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] += x * y
        return [v % mod for v in r] if mod else r

    def naive_mm(A, B, p):
        return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % p for j in range(len(B[0]))]
                for i in range(len(A))]

    # --- matrices ---
    p = 10**9 + 7
    for _ in range(30):
        n, m, k = rnd.randint(1, 6), rnd.randint(1, 6), rnd.randint(1, 6)
        A = [[rnd.randrange(p) for _ in range(m)] for _ in range(n)]
        B = [[rnd.randrange(p) for _ in range(k)] for _ in range(m)]
        assert mat_mul(A, B, p) == naive_mm(A, B, p)
        S = [[rnd.randrange(p) for _ in range(n)] for _ in range(n)]
        e = rnd.randint(0, 20)
        Pw = [[int(i == j) for j in range(n)] for i in range(n)]
        for _ in range(e):
            Pw = naive_mm(Pw, S, p)
        assert mat_pow(S, e, p) == Pw
    for _ in range(30):
        d = rnd.randint(1, 5)
        c = [rnd.randrange(p) for _ in range(d)]
        a = [rnd.randrange(p) for _ in range(d)]
        seq = a[:]
        while len(seq) < 300:
            seq.append(sum(c[j] * seq[-1 - j] for j in range(d)) % p)
        for k in [0, d - 1, d, 57, 299]:
            assert linrec_matpow(c, a, k, p) == seq[k]
            assert linrec_kth(c, a, k, p) == seq[k]

    # --- exact helpers for tests ---
    def frac_rank(A):
        M = [[Fraction(v) for v in row] for row in A]
        r = 0
        for c in range(len(M[0]) if M else 0):
            piv = next((i for i in range(r, len(M)) if M[i][c] != 0), -1)
            if piv < 0:
                continue
            M[r], M[piv] = M[piv], M[r]
            for i in range(len(M)):
                if i != r and M[i][c] != 0:
                    f = M[i][c] / M[r][c]
                    M[i] = [x - f * y for x, y in zip(M[i], M[r])]
            r += 1
        return r

    def det_perm(A, m=None):
        n = len(A)
        tot = 0
        for perm in permutations(range(n)):
            inv = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
            t = 1
            for i in range(n):
                t *= A[i][perm[i]]
            tot += -t if inv & 1 else t
        return tot % m if m else tot

    def low_rank(n, m, r, lo, hi):
        L = [[rnd.randint(lo, hi) for _ in range(r)] for _ in range(n)]
        Rm = [[rnd.randint(lo, hi) for _ in range(m)] for _ in range(r)]
        return [[sum(L[i][t] * Rm[t][j] for t in range(r)) for j in range(m)] for i in range(n)]

    # --- gauss real ---
    for _ in range(200):
        n, m = rnd.randint(1, 5), rnd.randint(1, 5)
        A = low_rank(n, m, rnd.randint(0, min(n, m)), -3, 3)
        if rnd.random() < 0.5:
            x0 = [rnd.randint(-5, 5) for _ in range(m)]
            b = [sum(A[i][j] * x0[j] for j in range(m)) for i in range(n)]
        else:
            b = [rnd.randint(-5, 5) for _ in range(n)]
        rA = frac_rank(A)
        rAb = frac_rank([A[i] + [b[i]] for i in range(n)])
        expect = 0 if rAb > rA else (1 if rA == m else 2)
        st, x = gauss_real(A, b)
        assert st == expect, (A, b, st, expect)
        if st:
            assert all(abs(sum(A[i][j] * x[j] for j in range(m)) - b[i]) < 1e-6 for i in range(n))
        # mod version on same data
        q = 10**9 + 7
        st2, x2 = solve_mod(A, b, q)
        assert st2 == expect
        if st2:
            assert all((sum(A[i][j] * x2[j] for j in range(m)) - b[i]) % q == 0 for i in range(n))
        assert rank_mod(A, q) == rA

    # --- det / inverse ---
    for _ in range(100):
        n = rnd.randint(1, 6)
        A = [[rnd.randint(-10, 10) for _ in range(n)] for _ in range(n)]
        if rnd.random() < 0.2 and n > 1:
            A[0] = A[1][:]
        q = 1000003
        assert det_mod(A, q) == det_perm(A, q)
        for mm in (12, 10**9, 2**20, 1, 7):
            assert det_any_mod(A, mm) == det_perm(A, mm)
        I = mat_inv_mod(A, q)
        if det_perm(A, q) == 0:
            assert I is None
        else:
            assert mat_mul(A, I, q) == [[int(i == j) for j in range(n)] for i in range(n)]

    # --- GF(2) ---
    for _ in range(300):
        n, m = rnd.randint(1, 7), rnd.randint(1, 7)
        rows = [rnd.getrandbits(m + 1) for _ in range(n)]
        r, x = gauss_gf2(rows, m)
        sols = [y for y in range(1 << m)
                if all(bin(row & y & ((1 << m) - 1)).count('1') % 2 == row >> m & 1 for row in rows)]
        if sols:
            assert x in sols
            assert len(sols) == 1 << (m - r)
        else:
            assert x is None
        # rank vs span size
        span = {0}
        for row in rows:
            v = row & ((1 << m) - 1)
            span |= {s ^ v for s in span}
        assert len(span) == 1 << gauss_gf2([row & ((1 << m) - 1) for row in rows], m)[0]

    # --- XOR basis ---
    for _ in range(200):
        vals = [rnd.getrandbits(rnd.randint(1, 8)) for _ in range(rnd.randint(0, 8))]
        xb = XorBasis(10)
        for v in vals:
            xb.insert(v)
        span = {0}
        for v in vals:
            span |= {s ^ v for s in span}
        sp = sorted(span)
        assert xb.rank == len(sp).bit_length() - 1
        for k in range(len(sp) + 2):
            assert xb.kth(k) == (sp[k] if k < len(sp) else -1)
        for y in range(1 << 9):
            assert xb.can(y) == (y in span)
        x0 = rnd.getrandbits(9)
        assert xb.max_xor(x0) == max(x0 ^ s for s in span)
        assert xb.max_xor() == sp[-1]

    # --- Berlekamp-Massey ---
    for _ in range(50):
        d = rnd.randint(1, 10)
        c = [R() for _ in range(d)]
        a = [R() for _ in range(d)]
        seq = a[:]
        while len(seq) < 2 * d + 10:
            seq.append(sum(c[j] * seq[-1 - j] for j in range(d)) % MOD)
        cc = berlekamp_massey(seq)
        assert len(cc) <= d
        for i in range(len(cc), len(seq)):
            assert seq[i] == sum(cc[j] * seq[i - 1 - j] for j in range(len(cc))) % MOD
        full = seq[:]
        while len(full) < 2000:
            full.append(sum(c[j] * full[-1 - j] for j in range(d)) % MOD)
        assert linrec_kth(cc, seq[:len(cc)], 1999) == full[1999]
    assert berlekamp_massey([1, 1, 2, 3, 5, 8, 13]) == [1, 1]

    # --- convolve / NTT / FFT ---
    for _ in range(100):
        a = [rnd.randrange(-10**12, 10**12) for _ in range(rnd.randint(0, 60))]
        b = [rnd.randrange(-10**12, 10**12) for _ in range(rnd.randint(0, 60))]
        assert convolve(a, b) == naive_conv(a, b)
        for md in (MOD, 10**9 + 7, 2**64, 6):
            assert convolve(a, b, md) == naive_conv(a, b, md)
        a2 = [rnd.randrange(MOD) for _ in range(rnd.randint(1, 60))]
        b2 = [rnd.randrange(MOD) for _ in range(rnd.randint(1, 60))]
        assert ntt_convolve(a2, b2) == naive_conv(a2, b2, MOD)
        sa = [rnd.randrange(1000) for _ in range(rnd.randint(1, 60))]
        sb = [rnd.randrange(1000) for _ in range(rnd.randint(1, 60))]
        assert fft_multiply(sa, sb) == naive_conv(sa, sb)
    for n in (500, 3000):  # exercise decimal path + big-mod
        a = [R() for _ in range(n)]
        b = [R() for _ in range(n + 7)]
        assert convolve(a, b, MOD) == ntt_convolve(a, b)
        assert convolve(a, a, MOD) == ntt_convolve(a, a)
        A1 = [rnd.randrange(2**64) for _ in range(400)]
        B1 = [rnd.randrange(2**64) for _ in range(350)]
        assert convolve(A1, B1, 2**64 - 59) == naive_conv(A1, B1, 2**64 - 59)

    # --- power series ---
    def naive_inv(f, n):
        g = [0] * n
        i0 = pow(f[0], MOD - 2, MOD)
        for i in range(n):
            s = (1 if i == 0 else 0) - sum(f[j] * g[i - j] for j in range(1, min(i, len(f) - 1) + 1))
            g[i] = s * i0 % MOD
        return g

    def naive_exp(f, n):  # n g_n = sum_{k>=1} k f_k g_{n-k}
        g = [1] + [0] * (n - 1)
        for i in range(1, n):
            s = sum(k * f[k] * g[i - k] for k in range(1, min(i, len(f) - 1) + 1))
            g[i] = s % MOD * pow(i, MOD - 2, MOD) % MOD
        return g

    for n in (1, 2, 3, 5, 17, 64, 100, 777):
        f = [R() for _ in range(n)]
        f[0] = R() or 1
        assert poly_inv(f, n) == naive_inv(f, n)
        f0 = [0] + f[1:]
        e = poly_exp(f0, n)
        assert e == naive_exp(f0, n)
        assert poly_log(e, n) == f0
        f1 = [1] + f[1:]
        s = poly_sqrt(f1, n)
        assert _conv(s, s)[:n] == f1
        for k in (0, 1, 2, 3, 7):
            for shift in (0, 1, 3):
                g = [0] * shift + f[:max(1, n // 2)]
                ref = [1]
                for _ in range(k):
                    ref = naive_conv(ref, g, MOD)
                ref = (ref + [0] * n)[:n]
                assert poly_pow(g, k, n) == ref
        assert poly_pow([0, 1], 10**18, 5) == [0] * 5
        assert poly_pow([2, 1], 10**18, 3) == [pow(2, 10**18, MOD), 10**18 % MOD * pow(2, 10**18 - 1, MOD) % MOD,
                                                (10**18 * (10**18 - 1) // 2) % MOD * pow(2, 10**18 - 2, MOD) % MOD]

    for _ in range(40):
        a = [R() for _ in range(rnd.randint(1, 400))]
        b = [R() for _ in range(rnd.randint(1, 300))]
        b[-1] = b[-1] or 1
        q, r = poly_divmod(a, b)
        assert len(r) == len(b) - 1
        bq = convolve(b, q, MOD)
        lhs = [((bq[i] if i < len(bq) else 0) + (r[i] if i < len(r) else 0)) % MOD
               for i in range(max(len(bq), len(r), len(a)))]
        assert lhs == [x % MOD for x in a] + [0] * (len(lhs) - len(a))

    for m in (1, 5, 40, 300):
        f = [R() for _ in range(rnd.randint(1, 500))]
        xs = [R() for _ in range(m)]
        ref = []
        for x in xs:
            acc = 0
            for c in reversed(f):
                acc = (acc * x + c) % MOD
            ref.append(acc)
        assert multipoint_eval(f, xs) == ref

    # --- Lagrange ---
    for _ in range(30):
        n = rnd.randint(1, 30)
        coef = [R() for _ in range(n)]
        ev = lambda x: sum(c * pow(x, i, MOD) for i, c in enumerate(coef)) % MOD
        xs = rnd.sample(range(MOD), n)
        assert interpolate(xs, [ev(x) for x in xs]) == coef
        ys = [ev(i) for i in range(n)]
        for x in [0, n - 1, n, 10**18, R()]:
            assert lagrange_consecutive(ys, x) == ev(x % MOD)
    q = 10**9 + 7  # sum_{i=1}^{N} i^3 via consecutive interpolation
    ys = [sum(i**3 for i in range(1, t + 1)) % q for t in range(5)]
    N = 10**12
    assert lagrange_consecutive(ys, N, q) == (N * (N + 1) // 2) ** 2 % q

    # --- bitwise convolutions / subset convolution ---
    for n in range(0, 7):
        N = 1 << n
        a = [rnd.randrange(MOD) for _ in range(N)]
        b = [rnd.randrange(MOD) for _ in range(N)]
        for kind, op in (('xor', lambda i, j: i ^ j), ('or', lambda i, j: i | j), ('and', lambda i, j: i & j)):
            ref = [0] * N
            for i in range(N):
                for j in range(N):
                    ref[op(i, j)] += a[i] * b[j]
            assert bit_convolve(a, b, kind) == [v % MOD for v in ref]
        ref = [0] * N
        for s in range(N):
            t = s
            while True:
                ref[s] += a[t] * b[s ^ t]
                if t == 0:
                    break
                t = (t - 1) & s
        assert subset_convolve(a, b) == [v % MOD for v in ref]

    # --- online convolution ---
    for n in (1, 2, 31, 33, 100, 700):
        g = [R() for _ in range(n)]
        h = [R() for _ in range(n)]
        f = online_convolve(n, g, lambda i, s: h[i] + s * 3)
        ref = []
        for i in range(n):
            s = sum(ref[j] * g[i - j] for j in range(i)) % MOD
            ref.append((h[i] + 3 * s) % MOD)
        assert f == ref

    # --- Bostan-Mori ---
    for _ in range(30):
        dq = rnd.randint(1, 15)
        Q = [R() for _ in range(dq + 1)]
        Q[0] = Q[0] or 1
        P = [R() for _ in range(rnd.randint(0, 25))]
        n = 300
        ser = convolve(P, poly_inv(Q, n), MOD)[:n] if P else [0] * n
        ser += [0] * (n - len(ser))
        for k in (0, 1, 2, 50, 299):
            assert bostan_mori(P, Q, k) == ser[k]
        q = 10**9 + 7
        assert bostan_mori([1], [1, q - 1, q - 1], 90, q) == 4660046610375530309 % q  # Fib(91)

    print("ALL OK")
