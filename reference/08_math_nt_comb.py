# 08_math_nt_comb.py -- Number theory, combinatorics, probability, game theory, misc math.
# Python 3 stdlib only, PyPy-friendly. Matrices/Gauss/NTT/GF ops/Lagrange/BM -> 09_algebra_poly.py.
import random
from math import gcd, isqrt, comb, factorial
from functools import lru_cache
from fractions import Fraction
from itertools import compress

# ==== 9.1 Number theory: divisibility, trial-division factorization, divisors ====
# Builtins: math.gcd(*args) (3.9+ variadic), math.lcm(*args) (3.9+), math.isqrt, math.comb, math.perm.


def factorize(n):  # O(sqrt n). {prime: exponent}
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def divisors(n):  # O(sqrt n) factorize + O(d(n)). Unsorted -> sort if needed. d(n)<=1344 for n<=1e9
    ds = [1]
    for p, e in factorize(n).items():
        ds = [d * p ** k for d in ds for k in range(e + 1)]
    return ds


# ==== 9.1 Sieve of Eratosthenes (bytearray slicing) ====
def sieve(n):  # O(n log log n). is_p[i] in {0,1}; 1e7 in ~1s CPython. primes = list(compress(range(n+1), is_p))
    is_p = bytearray([1]) * (n + 1)
    is_p[0:2] = b'\x00\x00'[:n + 1]
    for i in range(2, isqrt(n) + 1):
        if is_p[i]:
            is_p[i * i::i] = bytes(len(range(i * i, n + 1, i)))
    return is_p


def primes_upto(n):
    return list(compress(range(n + 1), sieve(n))) if n >= 2 else []


# ==== 9.1 Linear sieve: spf, primes, phi, mu ====
def linear_sieve(n):  # O(n) (slow-ish in CPython: ~3s for 1e7; use sieve() if only primes needed)
    spf = [0] * (n + 1); phi = [0] * (n + 1); mu = [0] * (n + 1); primes = []
    if n >= 1:
        phi[1] = mu[1] = 1
    for i in range(2, n + 1):
        if spf[i] == 0:
            spf[i] = i; primes.append(i); phi[i] = i - 1; mu[i] = -1
        si = spf[i]
        for p in primes:
            ip = i * p
            if p > si or ip > n:
                break
            spf[ip] = p
            if p == si:
                phi[ip] = phi[i] * p   # mu[ip] stays 0
            else:
                phi[ip] = phi[i] * (p - 1); mu[ip] = -mu[i]
    return spf, primes, phi, mu
# factor x<=n with spf: while x>1: p=spf[x]; x//=p  -> O(log x)


# ==== 9.1 Segmented sieve ====
def segmented_sieve(L, R):  # primes in [L, R]; O((R-L) log log R + sqrt R). R-L ~1e7 ok
    L = max(L, 2)
    if L > R:
        return []
    seg = bytearray([1]) * (R - L + 1)
    for p in primes_upto(isqrt(R)):
        s = max(p * p, (L + p - 1) // p * p)
        if s <= R:
            seg[s - L::p] = bytes(len(range(s - L, R - L + 1, p)))
    return [L + i for i in compress(range(R - L + 1), seg)]


# ==== 9.1 GCD/LCM, extended Euclid, Bezout ====
# gcd(a,b)=gcd(b,a%b); lcm=a//gcd(a,b)*b. Bezout: ax+by=c solvable iff gcd(a,b)|c.
def ext_gcd(a, b):  # O(log). returns (g, x, y): a*x + b*y = g = gcd(a,b) >= 0. |x|<=|b|, |y|<=|a|
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0


# ==== 9.1 Linear Diophantine equations ====
def diophantine(a, b, c):  # ax+by=c, (a,b)!=(0,0). None or (x0,y0,dx,dy): x=x0+k*dx, y=y0-k*dy, k in Z
    g, x, y = ext_gcd(a, b)
    if c % g:
        return None
    return x * (c // g), y * (c // g), b // g, a // g


def _krange(v0, step, lo, hi):  # k with lo <= v0 + k*step <= hi
    INF = float('inf')
    if step == 0:
        return (-INF, INF) if lo <= v0 <= hi else (1, 0)
    if step < 0:
        v0, step, lo, hi = -v0, -step, -hi, -lo
    return -((v0 - lo) // step), (hi - v0) // step   # ceil((lo-v0)/s), floor((hi-v0)/s)


def count_diophantine(a, b, c, x1, x2, y1, y2):  # #solutions with x in [x1,x2], y in [y1,y2]
    s = diophantine(a, b, c)
    if s is None:
        return 0
    x0, y0, dx, dy = s
    l1, r1 = _krange(x0, dx, x1, x2)
    l2, r2 = _krange(y0, -dy, y1, y2)
    lo, hi = max(l1, l2), min(r1, r2)
    return 0 if lo > hi else int(hi - lo + 1)   # (infinite only if a=0 or b=0 with unbounded k: n/a here)


# ==== 9.1 Modular arithmetic, inverse, fast exponentiation ====
# - Python % is always non-negative for positive mod: (-3) % 7 == 4. Bigints never overflow.
# - Fast exp: built-in pow(a, e, m) (C speed). Inverse: pow(a, -1, m) (3.8+, any m with gcd=1,
#   raises ValueError otherwise). Fermat: pow(a, m-2, m) for prime m. Euler: a^phi(m) = 1 (gcd=1).
# - a^e mod m with huge e, gcd(a,m)>1: if e >= log2(m): a^e = a^(e mod phi + phi) (mod m).
# - Keep values reduced; in PyPy, products < 2^63 stay machine ints (mod ~1e9 ok, mod 2^61 not).
def inv_table(n, m):  # O(n). inv[i] = i^{-1} mod m for 1<=i<=n, m prime > n
    inv = [0, 1] + [0] * (n - 1)
    for i in range(2, n + 1):
        inv[i] = (m - m // i) * inv[m % i] % m
    return inv[:n + 1]


# ==== 9.1 Euler's totient ====
def phi(n):  # O(sqrt n)
    r = n
    for p in factorize(n):
        r -= r // p
    return r
# sum_{d|n} phi(d) = n. phi multiplicative. #{1<=k<=n: gcd(k,n)=1} = phi(n); their sum = n*phi(n)/2 (n>1)


# ==== 9.1 Chinese Remainder Theorem (generalized, non-coprime moduli) ====
def crt(rs, ms):  # x = r_i (mod m_i). returns (x, lcm) with 0<=x<lcm, or None if inconsistent
    x, M = 0, 1
    for r, m in zip(rs, ms):
        g, p, _ = ext_gcd(M, m)
        if (r - x) % g:
            return None
        t = (r - x) // g * p % (m // g)
        x += M * t
        M = M // g * m
        x %= M
    return x, M


# ==== 9.1 Wilson's theorem, Lucas's theorem ====
# Wilson: (p-1)! = -1 (mod p) iff p prime. n! mod p for n>=p is 0; for n near p use (p-1)! = -1 and divide.
def lucas(n, k, p):  # C(n,k) mod small prime p; O(log_p n * cost). Precompute fact mod p if p small & many queries
    res = 1
    while n or k:
        ni, ki = n % p, k % p
        if ki > ni:
            return 0
        res = res * comb(ni, ki) % p
        n //= p; k //= p
    return res
# C(n,k) mod prime power / composite: Granville / CRT over prime powers (rare).


# ==== 9.1 Primitive roots ====
def primitive_root(m):  # smallest g; exists iff m in {1,2,4,p^k,2p^k}. else -1. O(m^{?}) small in practice
    if m <= 4:
        return m - 1 if m > 1 else 0
    ph = phi(m)
    fs = list(factorize(ph))
    for g in range(2, m):
        if gcd(g, m) == 1 and all(pow(g, ph // q, m) != 1 for q in fs):
            return g
    return -1


# ==== 9.1 Discrete logarithm (BSGS, gcd(a,m) != 1 allowed) ====
def discrete_log(a, b, m):  # smallest x>=0 with a^x = b (mod m), or -1. O(sqrt m)
    a %= m; b %= m
    if m == 1:
        return 0
    k, add = 1, 0
    g = gcd(a, m)
    while g > 1:
        if b == k:
            return add
        if b % g:
            return -1
        b //= g; m //= g; add += 1
        k = k * (a // g) % m
        g = gcd(a, m)
    n = isqrt(m) + 1
    an = pow(a, n, m)
    tbl = {}
    cur = b
    for q in range(n + 1):        # baby: b*a^q (later q overwrites -> smallest x)
        tbl[cur] = q
        cur = cur * a % m
    cur = k
    for p in range(1, n + 1):     # giant: k*a^(np)
        cur = cur * an % m
        if cur in tbl:
            return n * p - tbl[cur] + add
    return -1


# ==== 9.1 Discrete root ====
def discrete_root(k, a, p):  # all x in [0,p) with x^k = a (mod p), p prime. sorted list. O(sqrt p)
    a %= p
    if a == 0:
        return [0]
    if p == 2:
        return [1]
    g = primitive_root(p)
    y = discrete_log(pow(g, k, p), a, p)   # (g^y)^k = a
    if y < 0:
        return []
    step = (p - 1) // gcd(k, p - 1)
    return sorted(pow(g, y + i * step, p) for i in range(gcd(k, p - 1)))


# ==== 9.1 Quadratic residues, Tonelli-Shanks ====
# Euler criterion: a^((p-1)/2) = 1 iff a is QR mod odd prime p. (p-1)/2 nonzero QRs.
def sqrt_mod(a, p):  # x with x*x = a (mod p), p prime; -1 if none. Other root: p-x. O(log^2 p)
    a %= p
    if a == 0 or p == 2:
        return a
    if pow(a, (p - 1) // 2, p) != 1:
        return -1
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


# ==== 9.1 Miller-Rabin (deterministic for n < 3.18e23, so all 64-bit) ====
# Bases = first 12 primes; 318665857834031151167461 is the first strong pseudoprime to all of them.
# Arbitrary size: same code is a probable-prime test; add random bases for more confidence.
_MR = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n):  # O(12 log n) modmuls
    if n < 2:
        return False
    for p in _MR:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in _MR:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


# ==== 9.1 Pollard's rho (Brent) full factorization ====
def pollard_rho(n):  # nontrivial factor of composite n, expected O(n^{1/4})
    if n % 2 == 0:
        return 2
    while True:
        y, c, m = random.randrange(1, n), random.randrange(1, n), 128
        g = r = q = 1
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = gcd(q, n); k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = gcd(abs(x - ys), n)
        if g != n:
            return g


def factor(n):  # sorted prime factors with multiplicity; n ~1e18 in ms. Counter(factor(n)) for dict
    res = []
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        while n % p == 0:
            res.append(p); n //= p
    st = [n]
    while st:
        n = st.pop()
        if n == 1:
            continue
        if is_prime(n):
            res.append(n); continue
        d = pollard_rho(n)
        st += [d, n // d]
    return sorted(res)


# ==== 9.1 Multiplicative functions, Mobius, Dirichlet convolution (notes) ====
# f multiplicative: f(ab)=f(a)f(b) for gcd(a,b)=1 -> determined by f(p^k). Linear sieve computes any
#   such f if you also track the power of spf (cnt[i], or i's spf-free part).
# Examples: phi, mu, d(n)=prod(e+1), sigma(n)=prod (p^(e+1)-1)/(p-1), id, 1, e(n)=[n==1].
# Dirichlet: (f*g)(n) = sum_{d|n} f(d) g(n/d). mu*1=e, phi*1=id, mu*id=phi, 1*1=d, id*1=sigma.
#   Computing f*g up to n by "for d: for multiples" is O(n log n).
# Mobius inversion: g = f*1  <=>  f = g*mu. Typical: [gcd(a,b)==1] = sum_{d|gcd} mu(d).
def coprime_pairs(n, mu_prefix=None):  # #{(a,b) in [1,n]^2 : gcd(a,b)=1} = sum_d mu(d)*(n//d)^2. O(n) sieve + O(sqrt n)
    if mu_prefix is None:
        mu = linear_sieve(n)[3]
        mu_prefix = [0] * (n + 1)
        for i in range(1, n + 1):
            mu_prefix[i] = mu_prefix[i - 1] + mu[i]
    res = 0
    for l, r, v in floor_blocks(n):
        res += (mu_prefix[r] - mu_prefix[l - 1]) * v * v
    return res


# ==== 9.1 Floor-division blocks (sqrt(n) distinct values of n//i) ====
def floor_blocks(n):  # yields (l, r, v): n//i == v for all l<=i<=r; ~2sqrt(n) blocks
    l = 1
    while l <= n:
        v = n // l
        r = n // v
        yield l, r, v
        l = r + 1
# Divisor-sum trick: sum_{i<=n} d(i) = sum_i n//i ; sum_{i<=n} sigma(i) = sum_i i*(n//i).
def divisor_count_sum(n):  # sum_{i=1}^n d(i), O(sqrt n)
    return sum((r - l + 1) * v for l, r, v in floor_blocks(n))


# ==== 9.1 Prefix sums of multiplicative functions: Lucy_Hedgehog prime counting ====
def prime_pi(n):  # pi(n) in O(n^{3/4}); 1e10 ~ seconds in PyPy, 1e9 ~ few s CPython
    # Sum of primes: init with v(v+1)/2-1 and subtract p*(S[v/p]-S[p-1]).
    if n < 2:
        return 0
    r = isqrt(n)
    lo = [i - 1 for i in range(r + 1)]                    # lo[v] = S(v), v<=r
    hi = [0] + [n // i - 1 for i in range(1, r + 1)]      # hi[i] = S(n//i)
    for p in range(2, r + 1):
        if lo[p] == lo[p - 1]:
            continue
        sp = lo[p - 1]
        p2 = p * p
        for i in range(1, min(r, n // p2) + 1):
            d = i * p
            hi[i] -= (hi[d] if d <= r else lo[n // d]) - sp
        for v in range(r, p2 - 1, -1):
            lo[v] -= lo[v // p] - sp
    return hi[1]


# ==== 9.1 Du's sieve (prefix sum of phi / mu, n ~ 1e9..1e10) ====
def du_sieve(n, pre, total):
    # S(x) = sum_{i<=x} f(i), given g=1 (so f*1 known): S(x) = total(x) - sum_{d=2}^x S(x//d).
    # pre[x] = S(x) for x < len(pre) (from linear sieve, size ~ n^{2/3}). O(n^{2/3}).
    # phi: total = x(x+1)/2 ; mu: total = 1.
    memo = {}
    lim = len(pre)

    def S(x):
        if x < lim:
            return pre[x]
        if x in memo:
            return memo[x]
        res = total(x)
        l = 2
        while l <= x:
            v = x // l
            r = x // v
            res -= (r - l + 1) * S(v)
            l = r + 1
        memo[x] = res
        return res
    return S(n)   # recursion depth only O(log n)


def phi_sum(n, lim=None):  # sum_{i<=n} phi(i)
    lim = lim or max(10, int(n ** (2 / 3)))
    ph = linear_sieve(lim)[2]
    pre = [0] * (lim + 1)
    for i in range(1, lim + 1):
        pre[i] = pre[i - 1] + ph[i]
    return du_sieve(n, pre, lambda x: x * (x + 1) // 2)


def mertens(n, lim=None):  # sum_{i<=n} mu(i)
    lim = lim or max(10, int(n ** (2 / 3)))
    mu = linear_sieve(lim)[3]
    pre = [0] * (lim + 1)
    for i in range(1, lim + 1):
        pre[i] = pre[i - 1] + mu[i]
    return du_sieve(n, pre, lambda x: 1)
# Min_25 sieve (note): prefix sums of multiplicative f with f(p) a low-degree poly in p, n ~1e10.
#   Part 1 = Lucy-style sums of p^k over primes at all n//i; part 2 = recursive
#   S(x, j) = g(x) - g(p_{j-1}) + sum_{k>=j, p_k^2<=x} sum_e f(p_k^e)(S(x/p_k^e, k+1) + [e>1]). O(n^{3/4}/log n).


# ==== 9.1 Legendre's formula ====
def legendre(n, p):  # v_p(n!) = sum n//p^i = (n - digitsum_p(n))/(p-1). Kummer: v_p(C(m+n,m)) = #carries adding m+n base p
    e = 0
    while n:
        n //= p; e += n
    return e


# ==== 9.1 Continued fractions, Stern-Brocot, best rational approximation ====
def cont_frac(a, b):  # a/b -> [q0; q1, q2, ...]
    res = []
    while b:
        q = a // b
        res.append(q)
        a, b = b, a - q * b
    return res


def convergents(cf):  # list of (p, q) convergents p_k/q_k
    p0, q0, p1, q1 = 0, 1, 1, 0
    out = []
    for a in cf:
        p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
        out.append((p1, q1))
    return out
# Best rational approx with denominator <= N: Fraction(x).limit_denominator(N) (exact, uses CF+semiconvergents).
# Stern-Brocot: mediant (a+c)/(b+d) of a/b < c/d; path to x is run-length = its CF. Neighbours
#   a/b, c/d in SB/Farey satisfy bc - ad = 1. Smallest-denominator fraction in (L,R): walk SB with
#   CF steps (jump k at once) -> O(log).


# ==== 9.1 Pell's equation ====
def pell(D):  # fundamental (x, y) of x^2 - D y^2 = 1, D>0 non-square. Via CF of sqrt(D) (period <~ sqrt(D) log D)
    # All solutions: x_k + y_k sqrt(D) = (x1 + y1 sqrt(D))^k. x^2-Dy^2=-1 solvable iff period odd.
    a0 = isqrt(D)
    m, d, a = 0, 1, a0
    p0, p1, q0, q1 = 1, a0, 0, 1
    while p1 * p1 - D * q1 * q1 != 1:
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d
        p0, p1 = p1, a * p1 + p0
        q0, q1 = q1, a * q1 + q0
    return p1, q1


# ==== 9.2 Combinatorics: factorials, inverse factorials, C, P, multinomial ====
# Sum rule / product rule. P(n,k)=n!/(n-k)!, C(n,k)=n!/(k!(n-k)!), Pascal: C(n,k)=C(n-1,k-1)+C(n-1,k).
# Exact bigints: math.comb / math.perm / math.factorial. Pascal table O(n^2) works for any modulus.
class Comb:  # O(n) precompute, O(1) query. mod must be prime > n
    def __init__(self, n, mod=998244353):
        self.mod = mod
        f = [1] * (n + 1)
        for i in range(1, n + 1):
            f[i] = f[i - 1] * i % mod
        fi = [1] * (n + 1)
        fi[n] = pow(f[n], mod - 2, mod)
        for i in range(n, 0, -1):
            fi[i - 1] = fi[i] * i % mod
        self.f, self.fi = f, fi

    def C(self, n, k):
        if k < 0 or k > n or n < 0:
            return 0
        return self.f[n] * self.fi[k] % self.mod * self.fi[n - k] % self.mod

    def P(self, n, k):
        return 0 if k < 0 or k > n else self.f[n] * self.fi[n - k] % self.mod

    def multinomial(self, ks):  # (sum ks)! / prod(k!)
        r = self.f[sum(ks)]
        for k in ks:
            r = r * self.fi[k] % self.mod
        return r

    def H(self, n, k):  # stars & bars: #ways to put k identical balls into n boxes = C(n+k-1, k)
        return self.C(n + k - 1, k) if n else int(k == 0)


# ==== 9.2 Stars and bars, pigeonhole (notes) ====
# x1+..+xn = k, xi>=0: C(k+n-1, n-1); xi>=1: C(k-1, n-1); xi>=li: subtract sum li first.
# xi<=ui: inclusion-exclusion over the set of violated bounds.
# Lattice paths (0,0)->(a,b) right/up: C(a+b,a). Not crossing y=x (reflection): C(2n,n)-C(2n,n+1).
# Pigeonhole: n+1 items in n boxes -> some box has >=2. Prefix sums mod n: among n+1 prefixes two
#   are equal -> nonempty subarray with sum = 0 mod n. Erdos-Szekeres: seq of length rs+1 has
#   increasing subseq of length r+1 or decreasing of length s+1.


# ==== 9.2 Inclusion-exclusion ====
def count_coprime_to(n, primes):  # #{1<=x<=n : no p in primes divides x}. O(2^k)
    k = len(primes)
    res = 0
    for mask in range(1 << k):
        d, bits = 1, 0
        for i in range(k):
            if mask >> i & 1:
                d *= primes[i]; bits += 1
        res += (-1) ** bits * (n // d)
    return res


def surjections(n, k):  # #onto maps [n]->[k] = sum_i (-1)^i C(k,i) (k-i)^n = k! S(n,k)
    return sum((-1) ** i * comb(k, i) * (k - i) ** n for i in range(k + 1))


# ==== 9.2 Catalan, Stirling (1st/2nd), Bell, derangements ====
def catalan(n):  # C(2n,n)/(n+1); mod: C(2n,n)-C(2n,n+1). 1,1,2,5,14,42.. (parens, BSTs, triangulations, Dyck paths)
    return comb(2 * n, n) // (n + 1)


def stirling1(n, mod=None):  # unsigned c[i][k]: perms of i with k cycles. c(i,k)=c(i-1,k-1)+(i-1)c(i-1,k). O(n^2)
    c = [[0] * (n + 1) for _ in range(n + 1)]
    c[0][0] = 1
    for i in range(1, n + 1):
        for k in range(1, i + 1):
            v = c[i - 1][k - 1] + (i - 1) * c[i - 1][k]
            c[i][k] = v % mod if mod else v
    return c


def stirling2(n, mod=None):  # S[i][k]: partitions of i-set into k nonempty blocks. S(i,k)=S(i-1,k-1)+k S(i-1,k)
    S = [[0] * (n + 1) for _ in range(n + 1)]
    S[0][0] = 1
    for i in range(1, n + 1):
        for k in range(1, i + 1):
            v = S[i - 1][k - 1] + k * S[i - 1][k]
            S[i][k] = v % mod if mod else v
    return S


def bell(n, mod=None):  # B_0..B_n via Bell triangle, O(n^2). B = 1,1,2,5,15,52. B_{n+1} = sum C(n,k) B_k
    res = [1]
    row = [1]
    for _ in range(n):
        nr = [row[-1]]
        for x in row:
            v = nr[-1] + x
            nr.append(v % mod if mod else v)
        row = nr
        res.append(row[0])
    return res


def derangements(n, mod=None):  # D_0..D_n. D(i) = (i-1)(D(i-1)+D(i-2)); D(n) = round(n!/e)
    D = [1, 0] + [0] * max(0, n - 1)
    for i in range(2, n + 1):
        v = (i - 1) * (D[i - 1] + D[i - 2])
        D[i] = v % mod if mod else v
    return D[:n + 1]


# ==== 9.2 Partition numbers (pentagonal number theorem) ====
def partitions(n, mod=None):  # p(0..n), O(n sqrt n). p(i) = sum_k (-1)^{k+1} [p(i-k(3k-1)/2) + p(i-k(3k+1)/2)]
    p = [1] + [0] * n
    for i in range(1, n + 1):
        s, k = 0, 1
        while True:
            g1 = k * (3 * k - 1) // 2
            if g1 > i:
                break
            g2 = g1 + k
            t = p[i - g1] + (p[i - g2] if g2 <= i else 0)
            s += t if k & 1 else -t
            k += 1
        p[i] = s % mod if mod else s
    return p
# Partitions into distinct parts / bounded parts: knapsack DP (O(n^2)) or GF ops in 09.


# ==== 9.2 Burnside / Polya: necklaces and bracelets ====
# Burnside: #orbits = (1/|G|) sum_g |Fix(g)|. Polya: coloring with k colors fixed by g: k^{#cycles(g)}.
# Rotation by i of n beads has gcd(i,n) cycles.
def necklaces(n, k):  # colorings of n-cycle up to rotation: (1/n) sum_{d|n} phi(d) k^{n/d}
    return sum(phi(d) * k ** (n // d) for d in divisors(n)) // n


def bracelets(n, k):  # up to rotation + reflection (dihedral group, size 2n)
    rot = sum(phi(d) * k ** (n // d) for d in divisors(n))
    if n % 2:
        refl = n * k ** ((n + 1) // 2)
    else:
        refl = n // 2 * (k ** (n // 2 + 1) + k ** (n // 2))
    return (rot + refl) // (2 * n)


# ==== 9.2 Generating functions, LGV, bijections, Cayley (notes) ====
# OGF: sum a_n x^n (unlabeled; product = convolution). EGF: sum a_n x^n/n! (labeled; product =
#   binomial convolution; exp(F) = sets of F-structures). 1/(1-x)^k = sum C(n+k-1,k-1) x^n.
#   Series inverse/log/exp/pow and NTT -> see 09_algebra_poly.py.
# LGV: non-intersecting paths a_i->b_i in a DAG (only identity permutation possible) =
#   det(M), M[i][j] = #paths a_i->b_j. Determinant mod p -> see 09_algebra_poly.py.
# Cayley: n^{n-2} labeled trees; forests on n vertices with k given roots: k n^{n-k-1};
#   spanning trees of K_{a,b}: a^{b-1} b^{a-1}; general graph: Kirchhoff matrix-tree (det of Laplacian minor).
#   Prufer code <-> labeled tree; vertex degree = 1 + occurrences. #trees with degrees d_i: (n-2)!/prod(d_i-1)!
# Labeled graphs on n vertices: 2^{C(n,2)}; connected ones via c_n = g_n - sum_{k<n} C(n-1,k-1) c_k g_{n-k}.


# ==== 9.2 Young tableaux, hook length formula ====
def hook_length(shape):  # #standard Young tableaux of partition shape (rows non-increasing) = n!/prod hooks
    cols = [sum(1 for r in shape if r > j) for j in range(shape[0])] if shape else []
    prod = 1
    for i, r in enumerate(shape):
        for j in range(r):
            prod *= (r - j - 1) + (cols[j] - i - 1) + 1
    return factorial(sum(shape)) // prod
# RSK: permutations of n <-> pairs of SYT of same shape; sum f_lambda^2 = n!. LIS = first row length.


# ==== 9.3 Probability and expected value (notes) ====
# P(A|B) = P(A and B)/P(B). Bayes: P(A|B) = P(B|A)P(A)/P(B).
# Linearity: E[sum X_i] = sum E[X_i] even if dependent -> split into indicator variables.
# Contribution technique: answer = sum over elements/pairs of (value * #configs/prob it counts).
# Geometric: success prob p per trial -> E[#trials] = 1/p. Coupon collector: n*H_n.
#   "Expected time until": E[T] = sum_{t>=0} P(T > t).
# Expected value DP / game DP -> see 03_dp.py. Markov chains / random walks: E[steps from i] gives
#   linear system E_i = 1 + sum p_ij E_j -> Gaussian elimination (09_algebra_poly.py); on a path/tree
#   solve in O(n) by expressing E_i = a_i E_parent + b_i. Random walk on line +-1: E[hit 0 or n from k] = k(n-k).
# Answers mod p: probability P/Q -> P * pow(Q, -1, mod) % mod. Fractions add as a/b + c/d = (ad+bc)/bd.
def frac_mod(p, q, mod=998244353):
    return p % mod * pow(q, -1, mod) % mod


# ==== 9.6 Game theory: Nim, Sprague-Grundy ====
# Impartial, normal play: P-position (prev player wins) iff all moves go to N-positions.
# Nim: first player wins iff xor of piles != 0. Winning move: pile a with a ^ x < a, reduce to a ^ x.
# Sprague-Grundy: g(s) = mex{g(t)}; sum of games -> xor of Grundy values; lose iff 0.
def mex(vals):
    s = set(vals)
    m = 0
    while m in s:
        m += 1
    return m


def grundy_table(n, moves):  # states 0..n, moves(s) -> next states all < s. O(n * moves)
    g = [0] * (n + 1)
    for s in range(n + 1):
        g[s] = mex(g[t] for t in moves(s))
    return g
# General hashable states: @lru_cache def G(s): return mex(G(t) for t in moves(s)) (watch recursion depth).
# Tip: print grundy table for small n and look for a period / pattern.


# ==== 9.6 Nim variants: staircase, misere, Moore's, Wythoff ====
def staircase_nim_win(c):  # c[i] coins on stair i (c[0] = ground); move coins i -> i-1. win iff xor odd stairs
    x = 0
    for i in range(1, len(c), 2):
        x ^= c[i]
    return x != 0


def misere_nim_win(piles):  # last to move loses
    x = 0
    for a in piles:
        x ^= a
    if max(piles, default=0) <= 1:
        return x == 0          # even number of 1-piles -> win
    return x != 0


def moore_nim_win(piles, k):  # may reduce up to k piles per move. lose iff every bit column count = 0 mod (k+1)
    for b in range(max(piles, default=0).bit_length()):
        if sum(a >> b & 1 for a in piles) % (k + 1):
            return True
    return False


def wythoff_losing(a, b):  # take any from one pile or equal from both. P-positions: (floor(k phi), floor(k phi^2))
    if a > b:
        a, b = b, a
    k = b - a
    return a == (k + isqrt(5 * k * k)) // 2   # exact integer floor(k*phi)
# Partisan games (conceptual): Conway surreal numbers; Hackenbush values; positions = {L | R}.


# ==== 9.6 Minimax with alpha-beta ====
def alphabeta(s, depth, children, evaluate, alpha=float('-inf'), beta=float('inf'), maxi=True):
    # children(s) -> list of states (order best-first for pruning), evaluate(s) -> score for max player
    ch = children(s) if depth else []
    if not ch:
        return evaluate(s)
    if maxi:
        v = float('-inf')
        for t in ch:
            v = max(v, alphabeta(t, depth - 1, children, evaluate, alpha, beta, False))
            alpha = max(alpha, v)
            if alpha >= beta:
                break
    else:
        v = float('inf')
        for t in ch:
            v = min(v, alphabeta(t, depth - 1, children, evaluate, alpha, beta, True))
            beta = min(beta, v)
            if alpha >= beta:
                break
    return v


# ==== 9.7 Misc: series formulas (notes) ====
# sum i = n(n+1)/2, sum i^2 = n(n+1)(2n+1)/6, sum i^3 = (n(n+1)/2)^2.
# Geometric: sum_{i<n} r^i = (r^n - 1)/(r - 1); mod p with r-1 = 0 mod p -> n; else use inverse,
#   or divide-and-conquer / matrix [[r,1],[0,1]]^n for any modulus. Arithmetic: n(a1+an)/2.
# sum i^k: Faulhaber / Lagrange interpolation (09). H_n ~ ln n + 0.5772.


# ==== 9.7 Fibonacci: fast doubling, Zeckendorf, Pisano period ====
# F(2k) = F(k)(2F(k+1) - F(k)), F(2k+1) = F(k)^2 + F(k+1)^2. gcd(F_m,F_n) = F_gcd(m,n).
# Cassini F(n-1)F(n+1) - F(n)^2 = (-1)^n. sum_{i<=n} F_i = F_{n+2} - 1.
def fib(n, mod=None):  # (F(n), F(n+1)), O(log n)
    a, b = 0, 1
    for bit in bin(n)[2:]:
        c = a * (2 * b - a)
        d = a * a + b * b
        if mod:
            c %= mod; d %= mod
        if bit == '1':
            a, b = d, (c + d) % mod if mod else c + d
        else:
            a, b = c, d
    return a, b


def zeckendorf(n):  # n = sum of non-consecutive Fibonacci numbers (greedy), descending
    fs = [1, 2]
    while fs[-1] <= n:
        fs.append(fs[-1] + fs[-2])
    res = []
    for f in reversed(fs):
        if f <= n:
            res.append(f); n -= f
    return res


def pisano(m):  # period of F mod m, O(pi(m)) <= 6m. pi(lcm) = lcm of pi(prime powers)
    if m == 1:
        return 1
    a, b, i = 0, 1, 0
    while True:
        a, b = b, (a + b) % m
        i += 1
        if a == 0 and b == 1:
            return i


# ==== 9.7 Big integers, fractions (notes) ====
# Python int is arbitrary precision: + - * fast (Karatsuba), // and str() of 1e5+ digit numbers slow
#   (str(int) limited to 4300 digits in 3.11+: sys.set_int_max_str_digits(0)). Use int.bit_length,
#   int.to_bytes, divmod, pow(a, e, m). isqrt exact for bigints; use it, never int(x**0.5).
# Fraction(a, b): exact rationals, auto-reduced, slow (gcd each op). Fraction(0.1) is exact binary
#   value -> use Fraction('0.1'). limit_denominator(N) = best approximation. Compare a/b < c/d as a*d < c*b.
# Decimal with getcontext().prec = 50 for high-precision decimals.


# ==== 9.7 Base conversion (arbitrary base, negative base) ====
def to_base(n, b):  # digits MSB-first. b >= 2 needs n >= 0; b <= -2 works for any n. int(s, b) for b<=36 strings
    if n == 0:
        return [0]
    d = []
    while n:
        n, r = divmod(n, b)
        if r < 0:            # only for negative base
            n += 1; r -= b
        d.append(r)
    return d[::-1]


def from_base(d, b):
    v = 0
    for x in d:
        v = v * b + x
    return v
# bin(n), oct(n), hex(n), format(n, 'b'); int('ff', 16). Gray code: g = n ^ (n >> 1).


# ==== 9.7 Numerical methods: Simpson, Newton, bisection ====
def simpson(f, a, b, eps=1e-10):  # adaptive Simpson integral of f over [a,b]
    def S(a, b, fa, fm, fb):
        return (b - a) / 6 * (fa + 4 * fm + fb)

    def rec(a, b, fa, fm, fb, whole, eps, depth):
        m = (a + b) / 2
        lm, rm = (a + m) / 2, (m + b) / 2
        flm, frm = f(lm), f(rm)
        left, right = S(a, m, fa, flm, fm), S(m, b, fm, frm, fb)
        if depth <= 0 or abs(left + right - whole) <= 15 * eps:
            return left + right + (left + right - whole) / 15
        return rec(a, m, fa, flm, fm, left, eps / 2, depth - 1) + rec(m, b, fm, frm, fb, right, eps / 2, depth - 1)
    fa, fm, fb = f(a), f((a + b) / 2), f(b)
    return rec(a, b, fa, fm, fb, S(a, b, fa, fm, fb), eps, 50)
# Bisection: for _ in range(100): m=(lo+hi)/2; (lo,hi) = (m,hi) if pred(m) else (lo,m). Fixed
#   iteration count beats eps checks. Ternary search for unimodal (or golden-section).
# Newton: x -= f(x)/f'(x), quadratic convergence near simple root. Integer sqrt: math.isqrt;
#   integer k-th root: r = int(round(n ** (1/k))) then adjust r +-1 while r**k > n / (r+1)**k <= n.


# ==== 9.7 Linear programming: simplex ====
def simplex(A, b, c, eps=1e-9):
    # maximize c^T x s.t. A x <= b, x >= 0. returns (value, x); value = -inf infeasible, inf unbounded.
    # Two-phase, Bland's rule (no cycling). O(m n * pivots), fine for ~100x100. Use Fraction for exact.
    m, n = len(b), len(c)
    N = list(range(n)) + [-1]
    B = list(range(n, n + m))
    D = [[0.0] * (n + 2) for _ in range(m + 2)]
    for i in range(m):
        for j in range(n):
            D[i][j] = A[i][j]
        D[i][n] = -1; D[i][n + 1] = b[i]
    for j in range(n):
        D[m][j] = -c[j]
    D[m + 1][n] = 1

    def pivot(r, s):
        Dr = D[r]
        inv = 1 / Dr[s]
        for i in range(m + 2):
            if i != r and abs(D[i][s]) > eps:
                Di = D[i]
                f = Di[s] * inv
                for j in range(n + 2):
                    Di[j] -= Dr[j] * f
                Di[s] = Dr[s] * f
        for j in range(n + 2):
            if j != s:
                Dr[j] *= inv
        for i in range(m + 2):
            if i != r:
                D[i][s] *= -inv
        Dr[s] = inv
        B[r], N[s] = N[s], B[r]

    def run(phase):
        x = m + phase - 1
        while True:
            s = -1
            for j in range(n + 1):
                if N[j] != -phase and (s == -1 or (D[x][j], N[j]) < (D[x][s], N[s])):
                    s = j
            if D[x][s] >= -eps:
                return True
            r = -1
            for i in range(m):
                if D[i][s] <= eps:
                    continue
                if r == -1 or (D[i][n + 1] / D[i][s], B[i]) < (D[r][n + 1] / D[r][s], B[r]):
                    r = i
            if r == -1:
                return False
            pivot(r, s)

    r = min(range(m), key=lambda i: D[i][n + 1]) if m else 0
    if m and D[r][n + 1] < -eps:
        pivot(r, n)
        if not run(2) or D[m + 1][n + 1] < -eps:
            return float('-inf'), None
        for i in range(m):
            if B[i] == -1:
                s = 0
                for j in range(1, n + 1):
                    if s == -1 or (D[i][j], N[j]) < (D[i][s], N[s]):
                        s = j
                pivot(i, s)
    ok = run(1)
    x = [0.0] * n
    for i in range(m):
        if B[i] < n:
            x[B[i]] = D[i][n + 1]
    return (D[m][n + 1] if ok else float('inf')), x
# Duality: max c^T x, Ax<=b, x>=0  ==  min b^T y, A^T y >= c, y>=0. Equality -> two inequalities;
#   free var x = x+ - x-; minimize -> maximize -c.


# ==== 9.7 Floor sum ====
def floor_sum(n, m, a, b):  # sum_{i=0}^{n-1} floor((a*i+b)/m), m>0, any sign a,b. O(log m)
    ans = 0
    while True:
        if a >= m or a < 0:
            ans += n * (n - 1) // 2 * (a // m); a %= m
        if b >= m or b < 0:
            ans += n * (b // m); b %= m
        y = a * n + b
        if y < m:
            return ans
        n, b, m, a = y // m, y % m, a, m


# ==================================== TESTS ====================================
if __name__ == "__main__":
    from itertools import permutations, product
    from collections import Counter
    random.seed(12345)

    # factorize / divisors / phi
    for n in range(1, 2000):
        f = factorize(n)
        t = 1
        for p, e in f.items():
            t *= p ** e
        assert t == n and all(is_prime(p) for p in f)
        assert sorted(divisors(n)) == [d for d in range(1, n + 1) if n % d == 0]
        assert phi(n) == sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)
    # sieves
    N = 5000
    brute = [i for i in range(N + 1) if i >= 2 and all(i % d for d in range(2, isqrt(i) + 1))]
    assert primes_upto(N) == brute and primes_upto(1) == [] and primes_upto(2) == [2]
    spf, pr, ph, mu = linear_sieve(N)
    assert pr == brute
    for n in range(1, N + 1):
        assert ph[n] == phi(n)
        f = factorize(n)
        assert mu[n] == (0 if any(e > 1 for e in f.values()) else (-1) ** len(f))
        if n > 1:
            assert spf[n] == min(f)
    assert len(primes_upto(10 ** 6)) == 78498
    for _ in range(50):
        L = random.randint(0, 4000); R = L + random.randint(-5, 1000)
        assert segmented_sieve(L, R) == [p for p in brute if L <= p <= R]
    assert len(segmented_sieve(10 ** 9, 10 ** 9 + 10 ** 5)) == sum(is_prime(x) for x in range(10 ** 9, 10 ** 9 + 10 ** 5 + 1))

    # ext_gcd / diophantine
    for _ in range(2000):
        a, b = random.randint(-50, 50), random.randint(-50, 50)
        g, x, y = ext_gcd(a, b)
        assert g == gcd(a, b) and a * x + b * y == g
        if a == 0 or b == 0:
            continue
        c = random.randint(-100, 100)
        x1, y1 = random.randint(-30, 0), random.randint(-30, 0)
        x2, y2 = x1 + random.randint(0, 40), y1 + random.randint(0, 40)
        bf = sum(1 for x in range(x1, x2 + 1) for y in range(y1, y2 + 1) if a * x + b * y == c)
        assert count_diophantine(a, b, c, x1, x2, y1, y2) == bf, (a, b, c)
        s = diophantine(a, b, c)
        if s:
            x0, y0, dx, dy = s
            for k in range(-3, 4):
                assert a * (x0 + k * dx) + b * (y0 - k * dy) == c
    # inverses
    MOD = 10 ** 9 + 7
    iv = inv_table(1000, MOD)
    assert all(i * iv[i] % MOD == 1 for i in range(1, 1001))
    assert all(iv[i] == pow(i, -1, MOD) for i in range(1, 1001))
    # CRT
    for _ in range(3000):
        k = random.randint(1, 4)
        ms = [random.randint(1, 20) for _ in range(k)]
        rs = [random.randint(0, m - 1) for m in ms]
        res = crt(rs, ms)
        M = 1
        for m in ms:
            M = M * m // gcd(M, m)
        sols = [x for x in range(M) if all(x % m == r for r, m in zip(rs, ms))]
        assert (res is None) == (not sols)
        if res:
            assert res == (sols[0], M)
    # Wilson, Lucas
    for p in (2, 3, 5, 7, 11, 13):
        assert factorial(p - 1) % p == p - 1
        for n in range(60):
            for k in range(n + 1):
                assert lucas(n, k, p) == comb(n, k) % p
    # primitive root, discrete log, discrete root, sqrt mod
    for m in range(1, 200):
        g = primitive_root(m)
        units = [x for x in range(1, m) if gcd(x, m) == 1] if m > 1 else [0]
        ph_m = phi(m)
        has = any(len({pow(x, e, m) for e in range(1, ph_m + 1)}) == ph_m for x in units) if m > 1 else True
        if not has:
            assert g == -1
        elif m > 1:
            assert len({pow(g, e, m) for e in range(1, ph_m + 1)}) == ph_m
            assert g == min(x for x in units if len({pow(x, e, m) for e in range(1, ph_m + 1)}) == ph_m)
    for m in range(1, 60):
        for a in range(m):
            for b in range(m):
                bf = next((x for x in range(2 * m + 2) if pow(a, x, m) == b % m), -1)
                assert discrete_log(a, b, m) == bf, (a, b, m)
    assert pow(3, discrete_log(3, 123456, 10 ** 9 + 7), 10 ** 9 + 7) == 123456
    for p in (2, 3, 5, 7, 13, 17, 31, 41, 97):
        for k in range(1, 8):
            for a in range(p):
                assert discrete_root(k, a, p) == [x for x in range(p) if pow(x, k, p) == a]
    for p in brute[:200] + [998244353, 10 ** 9 + 7]:
        for _ in range(20):
            a = random.randrange(p)
            r = sqrt_mod(a, p)
            if r == -1:
                assert pow(a, (p - 1) // 2, p) == p - 1
            else:
                assert r * r % p == a
    # Miller-Rabin, Pollard
    isp = sieve(10 ** 5)
    assert all(is_prime(i) == bool(isp[i]) for i in range(10 ** 5 + 1))
    for n in (2 ** 61 - 1, 10 ** 18 + 9, 3215031751, 3825123056546413051):
        assert is_prime(n) == (n in (2 ** 61 - 1, 10 ** 18 + 9))
    for _ in range(30):
        n = random.randint(2, 10 ** 18)
        fs = factor(n)
        t = 1
        for p in fs:
            t *= p
        assert t == n and all(is_prime(p) for p in fs)
    assert factor(1000000007 * 998244353) == [998244353, 1000000007]
    assert factor((2 ** 31 - 1) ** 2 * 4) == [2, 2, 2 ** 31 - 1, 2 ** 31 - 1]
    assert factor(1) == []
    # Mobius / floor blocks / divisor sums
    for n in range(1, 150):
        assert coprime_pairs(n) == sum(1 for a in range(1, n + 1) for b in range(1, n + 1) if gcd(a, b) == 1)
        blocks = list(floor_blocks(n))
        assert [v for l, r, v in blocks for _ in range(l, r + 1)] == [n // i for i in range(1, n + 1)]
        assert divisor_count_sum(n) == sum(len(divisors(i)) for i in range(1, n + 1))
    # Lucy, Du
    P = sieve(2 * 10 ** 6)
    acc = 0
    pref = [0] * (2 * 10 ** 6 + 1)
    for i in range(2 * 10 ** 6 + 1):
        acc += P[i]; pref[i] = acc
    for n in list(range(0, 300)) + [random.randint(1, 2 * 10 ** 6) for _ in range(20)]:
        assert prime_pi(n) == pref[n]
    assert prime_pi(10 ** 9) == 50847534
    LIM = 10 ** 6
    _, _, ph, mu = linear_sieve(LIM)
    sp, sm = [0] * (LIM + 1), [0] * (LIM + 1)
    for i in range(1, LIM + 1):
        sp[i] = sp[i - 1] + ph[i]; sm[i] = sm[i - 1] + mu[i]
    for n in [1, 2, 10, 999, 12345, 654321, 10 ** 6]:
        assert phi_sum(n, lim=50) == sp[n] and mertens(n, lim=50) == sm[n]
        assert phi_sum(n) == sp[n] and mertens(n) == sm[n]
    assert mertens(10 ** 9) == -222 and phi_sum(10 ** 8) == 3039635516365908
    # Legendre
    for n in range(200):
        for p in (2, 3, 5, 7):
            f = factorial(n); e = 0
            while f % p == 0:
                f //= p; e += 1
            assert legendre(n, p) == e
    # continued fractions, Pell
    for _ in range(200):
        a, b = random.randint(1, 10 ** 6), random.randint(1, 10 ** 6)
        cv = convergents(cont_frac(a, b))
        assert Fraction(*cv[-1]) == Fraction(a, b)
        for (p0, q0), (p1, q1) in zip(cv, cv[1:]):
            assert abs(p1 * q0 - p0 * q1) == 1
    assert Fraction(3.141592653589793).limit_denominator(1000) == Fraction(355, 113)
    for D in range(2, 40):
        if isqrt(D) ** 2 == D:
            continue
        y = 1
        while isqrt(D * y * y + 1) ** 2 != D * y * y + 1:
            y += 1
        assert pell(D) == (isqrt(D * y * y + 1), y)
    x, y = pell(61)
    assert (x, y) == (1766319049, 226153980)

    # Comb
    MD = 998244353
    cb = Comb(500, MD)
    for n in range(60):
        for k in range(-1, n + 2):
            assert cb.C(n, k) == (comb(n, k) % MD if 0 <= k <= n else 0)
            if 0 <= k <= n:
                assert cb.P(n, k) == factorial(n) // factorial(n - k) % MD
    assert cb.multinomial([2, 3, 1]) == 60
    for n in range(0, 5):
        for k in range(0, 6):
            assert cb.H(n, k) == sum(1 for t in product(range(k + 1), repeat=n) if sum(t) == k)
    # inclusion-exclusion
    for _ in range(100):
        ps = random.sample(brute[:10], random.randint(0, 4)); n = random.randint(0, 3000)
        assert count_coprime_to(n, ps) == sum(1 for x in range(1, n + 1) if all(x % p for p in ps))
    for n in range(7):
        for k in range(5):
            assert surjections(n, k) == sum(1 for f in product(range(k), repeat=n) if len(set(f)) == k)
    # Catalan / Stirling / Bell / derangements
    def dyck(n):
        dp = {0: 1}
        for _ in range(2 * n):
            nd = Counter()
            for h, c in dp.items():
                nd[h + 1] += c
                if h:
                    nd[h - 1] += c
            dp = nd
        return dp.get(0, 0)
    assert [catalan(n) for n in range(10)] == [dyck(n) for n in range(10)]
    s1, s2 = stirling1(8), stirling2(8)
    for n in range(8):
        cyc = Counter()
        for pm in permutations(range(n)):
            seen, c = [0] * n, 0
            for i in range(n):
                if not seen[i]:
                    c += 1; j = i
                    while not seen[j]:
                        seen[j] = 1; j = pm[j]
            cyc[c] += 1
        for k in range(n + 1):
            assert s1[n][k] == cyc[k]
            assert s2[n][k] * factorial(k) == surjections(n, k)
    B = bell(12)
    assert B == [sum(stirling2(12)[n]) for n in range(13)]
    assert B[:7] == [1, 1, 2, 5, 15, 52, 203]
    D = derangements(9)
    assert D == [sum(1 for pm in permutations(range(n)) if all(pm[i] != i for i in range(n))) for n in range(10)]
    assert derangements(9, 7) == [d % 7 for d in D]
    # partitions
    PN = 300
    dp = [1] + [0] * PN
    for part in range(1, PN + 1):
        for s in range(part, PN + 1):
            dp[s] += dp[s - part]
    assert partitions(PN) == dp and partitions(PN, MD) == [v % MD for v in dp]
    # Burnside
    for n in range(1, 9):
        for k in range(1, 4):
            neck, brac = set(), set()
            for col in product(range(k), repeat=n):
                rots = [col[i:] + col[:i] for i in range(n)]
                neck.add(min(rots))
                brac.add(min(rots + [r[::-1] for r in rots]))
            assert necklaces(n, k) == len(neck) and bracelets(n, k) == len(brac), (n, k)
    # hook length
    @lru_cache(None)
    def syt(shape):
        if sum(shape) == 0:
            return 1
        res = 0
        for i, r in enumerate(shape):
            if r and (i + 1 == len(shape) or shape[i + 1] < r):
                s = list(shape); s[i] -= 1
                res += syt(tuple(s))
        return res
    def parts(n, mx):
        if n == 0:
            yield ()
            return
        for f in range(min(n, mx), 0, -1):
            for rest in parts(n - f, f):
                yield (f,) + rest
    for n in range(1, 11):
        tot = 0
        for sh in parts(n, n):
            assert hook_length(list(sh)) == syt(sh)
            tot += syt(sh) ** 2
        assert tot == factorial(n)
    assert frac_mod(1, 2) * 2 % MD == 1

    # games
    g = grundy_table(30, lambda s: [s - d for d in (1, 3, 4) if s >= d])
    assert g[:14] == [0, 1, 0, 1, 2, 3, 2, 0, 1, 0, 1, 2, 3, 2]
    @lru_cache(None)
    def win_sub(state):  # sum of subtraction games
        for i, s in enumerate(state):
            for d in (1, 3, 4):
                if s >= d and not win_sub(state[:i] + (s - d,) + state[i + 1:]):
                    return True
        return False
    for st in product(range(12), repeat=3):
        assert win_sub(st) == (g[st[0]] ^ g[st[1]] ^ g[st[2]] != 0)
    @lru_cache(None)
    def win_stair(c):
        for i in range(1, len(c)):
            for t in range(1, c[i] + 1):
                nc = list(c); nc[i] -= t; nc[i - 1] += t
                if not win_stair(tuple(nc)):
                    return True
        return False
    for c in product(range(3), repeat=4):
        assert win_stair(c) == staircase_nim_win(list(c))
    @lru_cache(None)
    def win_misere(pl):
        if sum(pl) == 0:
            return True
        for i, a in enumerate(pl):
            for t in range(a):
                if not win_misere(tuple(sorted(pl[:i] + (t,) + pl[i + 1:]))):
                    return True
        return False
    for pl in product(range(5), repeat=3):
        assert win_misere(tuple(sorted(pl))) == misere_nim_win(list(pl))
    @lru_cache(None)
    def win_moore(pl, k):
        n = len(pl)
        for mask in range(1, 1 << n):
            idx = [i for i in range(n) if mask >> i & 1]
            if len(idx) > k:
                continue
            for dec in product(*[range(1, pl[i] + 1) for i in idx]):
                q = list(pl)
                for i, d in zip(idx, dec):
                    q[i] -= d
                if not win_moore(tuple(q), k):
                    return True
        return False
    for pl in product(range(4), repeat=3):
        for k in (1, 2):
            assert win_moore(pl, k) == moore_nim_win(list(pl), k)
    W = 40
    lose = [[False] * W for _ in range(W)]
    for a in range(W):
        for b in range(W):
            mv = [(a - t, b) for t in range(1, a + 1)] + [(a, b - t) for t in range(1, b + 1)] + \
                 [(a - t, b - t) for t in range(1, min(a, b) + 1)]
            lose[a][b] = all(not lose[x][y] for x, y in mv)
            assert lose[a][b] == wythoff_losing(a, b)
    for _ in range(100):  # random game trees
        tree = {}
        def build(node, d):
            if d == 0 or random.random() < 0.2:
                tree[node] = []
                return
            tree[node] = [node * 10 + i for i in range(random.randint(1, 4))]
            for ch in tree[node]:
                build(ch, d - 1)
        build(1, 6)
        vals = {nd: random.randint(-50, 50) for nd in tree}
        def minimax(s, mx):
            if not tree[s]:
                return vals[s]
            r = [minimax(t, not mx) for t in tree[s]]
            return max(r) if mx else min(r)
        assert alphabeta(1, 10, lambda s: tree[s], lambda s: vals[s]) == minimax(1, True)

    # Fibonacci
    F = [0, 1]
    for _ in range(300):
        F.append(F[-1] + F[-2])
    for n in range(300):
        assert fib(n) == (F[n], F[n + 1]) and fib(n, 1000) == (F[n] % 1000, F[n + 1] % 1000)
    for n in range(1, 2000):
        z = zeckendorf(n)
        assert sum(z) == n
        idx = [F.index(f, 2) for f in z]
        assert all(idx[i] - idx[i + 1] >= 2 for i in range(len(idx) - 1))
    for m in range(1, 200):
        pm = pisano(m)
        assert fib(pm, m) == (0, 1 % m)
        assert all(fib(t, m) != (0, 1 % m) for t in range(1, pm)) or m == 1
    assert fib(10 ** 18 + 5, 1000)[0] == fib((10 ** 18 + 5) % pisano(1000), 1000)[0]
    # base conversion
    for n in range(0, 500):
        for b in (2, 3, 7, 10, 16):
            d = to_base(n, b)
            assert from_base(d, b) == n and all(0 <= x < b for x in d)
            assert int(''.join('0123456789abcdef'[x] for x in d), b) == n
    for n in range(-500, 500):
        for b in (-2, -3, -10):
            d = to_base(n, b)
            assert from_base(d, b) == n and all(0 <= x < -b for x in d) and (d[0] != 0 or n == 0)
    # Simpson
    import math
    assert abs(simpson(math.sin, 0, math.pi) - 2) < 1e-8
    assert abs(simpson(lambda x: math.sqrt(max(0.0, 1 - x * x)), -1, 1) - math.pi / 2) < 1e-7
    assert abs(simpson(lambda x: x ** 3, 0, 2) - 4) < 1e-12
    # simplex vs vertex enumeration (2-3 vars, box-bounded so never unbounded)
    def lp_brute(A, b, c):
        n = len(c)
        rows = [(A[i], b[i]) for i in range(len(b))] + [([-1 if j == k else 0 for j in range(n)], 0) for k in range(n)]
        best = None
        from itertools import combinations
        for sub in combinations(rows, n):
            M = [[Fraction(v) for v in r] + [Fraction(rb)] for r, rb in sub]
            ok = True
            for col in range(n):  # Gauss
                piv = next((i for i in range(col, n) if M[i][col] != 0), None)
                if piv is None:
                    ok = False; break
                M[col], M[piv] = M[piv], M[col]
                for i in range(n):
                    if i != col and M[i][col] != 0:
                        f = M[i][col] / M[col][col]
                        M[i] = [x - f * y for x, y in zip(M[i], M[col])]
            if not ok:
                continue
            x = [M[i][n] / M[i][i] for i in range(n)]
            if all(sum(r[j] * x[j] for j in range(n)) <= rb for r, rb in rows):
                v = sum(c[j] * x[j] for j in range(n))
                best = v if best is None else max(best, v)
        return best
    for it in range(400):
        n = random.randint(1, 3); m = random.randint(1, 4)
        A = [[random.randint(-5, 5) for _ in range(n)] for _ in range(m)]
        b = [random.randint(-6, 10) for _ in range(m)]
        A += [[1 if j == k else 0 for j in range(n)] for k in range(n)]
        b += [random.randint(1, 10) for _ in range(n)]
        c = [random.randint(-5, 5) for _ in range(n)]
        val, x = simplex(A, b, c)
        bf = lp_brute(A, b, c)
        if bf is None:
            assert val == float('-inf')
        else:
            assert abs(val - float(bf)) < 1e-6, (A, b, c, val, bf)
            assert all(sum(A[i][j] * x[j] for j in range(n)) <= b[i] + 1e-6 for i in range(len(b)))
            assert all(v >= -1e-9 for v in x)
    assert simplex([[1, -1]], [1], [1, 1])[0] == float('inf')        # unbounded
    assert abs(simplex([[1, 1], [1, 3]], [4, 6], [3, 2])[0] - 12) < 1e-9
    # floor_sum
    for _ in range(3000):
        n, m = random.randint(0, 30), random.randint(1, 30)
        a, b = random.randint(-60, 60), random.randint(-60, 60)
        assert floor_sum(n, m, a, b) == sum((a * i + b) // m for i in range(n))
    assert floor_sum(10 ** 9, 10 ** 9 - 7, 10 ** 9 + 9, 123) >= 0
    print("ALL OK")
