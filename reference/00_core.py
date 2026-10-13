# 00_core.py -- Python contest core: template, I/O, recursion, interactive, hashing, debugging.
# Python 3 stdlib only, PyPy-friendly.

# ==== 1.2 Choosing a language: Python performance notes ====
# - Submit with PyPy when available (5-50x faster on loops); CPython for heavy bigint / string builtins.
# - Rough budget: CPython ~1e7 simple ops/s, PyPy ~1e8. Push work into builtins
#   (sorted, sum, map, str.join, bytes.translate, list slicing) -- they run in C.
# - Local variables are faster than globals: put the solution inside main().
# - PyPy: lists of ints are unboxed and fast; avoid lists of tuples in hot loops (use parallel
#   arrays or pack (a<<20|b) into one int). Recursion is slow & stack-hungry on PyPy too.
# - PyPy: `sys.setrecursionlimit` big + deep recursion can segfault -> use bootstrap or iterative.
# - CPython: `x in set` O(1), but dict/set of tuples is slow & memory heavy -> encode as int.
# - Avoid `deepcopy`, `list.pop(0)` (O(n): use collections.deque), `str +=` in loops.
# - Fractions/Decimal are slow; prefer integer arithmetic (cross-multiply) where possible.

# ==== 1.3 Complexity: constraints -> intended complexity ====
# n<=10: n!   n<=20: 2^n*n   n<=40: 2^(n/2) MITM   n<=500: n^3   n<=5000: n^2
# n<=1e5..2e5: n log n (n sqrt n risky in Python)   n<=1e6: O(n)   n<=1e18: log n / O(1)
# In Python divide those budgets by ~10 (CPython) or ~2 (PyPy) vs C++.
# Memory: Python int ~28 bytes, list slot 8 bytes -> 1e7 ints list ~ 80MB+ (use array('i') to cut).

# ==== 1.4 Input/Output: contest template (fast input, fast output, multi-test) ====
import sys, os


def contest_template():
    # Read everything as bytes tokens; int(b'123') works directly. O(input size).
    data = sys.stdin.buffer.read().split()
    pos = 0

    def nxt():
        nonlocal pos
        pos += 1
        return data[pos - 1]

    out = []
    t = int(nxt())                      # drop if single test case
    for _ in range(t):
        n = int(nxt())
        a = [int(nxt()) for _ in range(n)]
        s = nxt().decode()              # string token
        out.append(str(sum(a)))
    sys.stdout.write("\n".join(out) + "\n")
# Gotchas: read until EOF -> just loop `while pos < len(data)`. Line-based input
# (strings with spaces): lines = sys.stdin.buffer.read().decode().split('\n').
# Avoid input()/print() per line for >1e5 lines; if needed: input = sys.stdin.readline (keeps '\n').
# Bulk ints from one line fast: a = list(map(int, data[pos:pos+n])); pos += n


class FastIO:
    """Token reader over a bytes buffer. io=FastIO(); n=io.int(); a=io.ints(n)."""
    def __init__(self, buf=None):
        self.d = (buf if buf is not None else sys.stdin.buffer.read()).split()
        self.p = 0

    def tok(self):
        self.p += 1
        return self.d[self.p - 1]

    def int(self):
        self.p += 1
        return int(self.d[self.p - 1])

    def ints(self, n):
        self.p += n
        return list(map(int, self.d[self.p - n:self.p]))

    def str(self):
        return self.tok().decode()

    def eof(self):
        return self.p >= len(self.d)


# Output precision: f"{x:.10f}" or "%.10f" % x. Never print(x) raw for floats (may give 1e-05 form).
def fmt_float(x, k=10):
    s = "%.*f" % (k, x)
    return "0." + "0" * k if s.startswith("-") and float(s) == 0 else s   # avoid "-0.000"
# Exact rationals: fractions.Fraction; modular p/q -> p * pow(q, MOD-2, MOD) % MOD.
# Special judge (multiple answers): any valid answer accepted -- don't over-think tie-breaking.


# ==== 1.4 Interactive problems ====
# Every query MUST be flushed; read the reply with readline (NOT the whole-buffer reader!).
def ask(*args, inp=None, outp=None):
    inp = inp or sys.stdin
    outp = outp or sys.stdout
    outp.write("? " + " ".join(map(str, args)) + "\n")
    outp.flush()
    r = inp.readline().split()
    if not r or r[0] == "-1":           # judge signals error: exit immediately (avoid TLE verdict)
        sys.exit(0)
    return list(map(int, r))


def answer(*args, outp=None):
    outp = outp or sys.stdout
    outp.write("! " + " ".join(map(str, args)) + "\n")
    outp.flush()
# Typical: binary search via queries -- count queries vs limit (ceil(log2 n)).
# Local test: write a judge script and connect with a FIFO:
#   mkfifo f; python3 judge.py < f | python3 sol.py > f


# ==== 1.2 Recursion strategies ====
# Default recursion limit is 1000. Option A: raise limit + run in a thread with a big stack.
import threading


def run_with_big_stack(main, stack_mb=256, limit=10 ** 6):
    # Usage: run_with_big_stack(main). CPython only reliable; on PyPy prefer bootstrap/iterative.
    sys.setrecursionlimit(limit)
    threading.stack_size(stack_mb << 20)
    th = threading.Thread(target=main)
    th.start()
    th.join()


# Option B: generator trampoline. Recursive calls become `yield f(x)`, returns become `yield val`
# (or plain `return val` in py3). Depth limited only by memory. Works on PyPy.
from types import GeneratorType


def bootstrap(f, stack=[]):
    def wrappedfunc(*args, **kwargs):
        if stack:
            return f(*args, **kwargs)
        to = f(*args, **kwargs)
        while True:
            if type(to) is GeneratorType:
                stack.append(to)
                to = next(to)
            else:
                stack.pop()
                if not stack:
                    break
                to = stack[-1].send(to)
        return to
    return wrappedfunc
# Example:
#   @bootstrap
#   def dfs(u, p):
#       s = 1
#       for v in g[u]:
#           if v != p: s += yield dfs(v, u)
#       yield s                      # MUST yield (or return) at end, even if None: `yield None`
# Gotcha: every path of the function must end with `yield value`; a non-generator call breaks it.
# Option C (best): manual stack -- iterative DFS with order list, process reversed order for post-order.


# ==== 4.3/19 Anti-hash protection (random salt for dict/set keys) ====
# Hack tests can make CPython/PyPy dicts degrade to O(n) per op with crafted int keys
# (and PyPy sets of ints). XOR keys with a random salt; same salt for insert & lookup.
import random

SALT = random.getrandbits(62) | 1


def salted_counter(arr):
    # O(n) expected. Returns dict keyed by x ^ SALT; look up with d.get(x ^ SALT, 0).
    d = {}
    for x in arr:
        k = x ^ SALT
        d[k] = d.get(k, 0) + 1
    return d
# Alternatives: str keys (string hash is randomized per process); or sort + groupby instead of dict.
# Counter/defaultdict have same weakness -- salt the keys before inserting.


# ==== 14.3 Implementation skills ====
# - Wrap per-test logic in solve() returning the answer; build fresh arrays per test (no stale globals).
# - Allocate arrays sized by this test's n, not max-N, when sum of n is bounded but t is large
#   (t=1e4 tests * [0]*2e5 init = TLE).
# - Python ints are arbitrary precision: no overflow, but huge ints get slow -- reduce mod early.
# - Keep MOD = 998244353 / 10**9+7 as a constant; `%` the result of every multiplication in loops.


# ==== 1.2/15 Python pitfalls ====
# - [[0]*m]*n aliases rows!  Use [[0]*m for _ in range(n)].
# - s += ch in a loop can be O(n^2); collect pieces in a list and "".join(parts).
# - 7 // 2 == 3, -7 // 2 == -4 (floors!), int(-7/2) == -3; -7 % 3 == 2. Use // for ints; / gives float.
# - Float loses precision past 2^53: int(1e18/3) is wrong; use integer //.  Use math.isqrt, not int(sqrt).
# - round() is banker's rounding: round(2.5) == 2. For half-up: int(x + 0.5) (x >= 0) or Decimal.
# - Default recursion limit 1000 -> RecursionError (RE); see Recursion strategies above.
# - sort(key=...) is fast; custom comparator needs functools.cmp_to_key (slow; prefer a key tuple).
#   Comparator must return negative/0/positive, not bool.
# - Default mutable args (def f(a=[])) persist across calls (bootstrap uses this on purpose).
# - Closures in loops capture the variable, not value: lambda x, i=i: ...
# - list.copy()/a[:] is shallow; nested lists need [row[:] for row in a].
# - max()/min() of empty sequence raises: use default=.
# - str.split() vs split(' '): the latter keeps empty strings for repeated spaces.
# - heapq is a min-heap; push -x for max-heap; tuples compare element-wise (ties compare next field!).
# - `is` vs `==`: never use `is` to compare ints/strings.
# - Modifying a list/dict while iterating over it -> bugs / RuntimeError.
# - input() strips '\n', sys.stdin.readline() does not -> .strip() / int() handles it.
# - Exponentiation: pow(a, b, m) is fast modular; a**b % m builds the huge number first.


# ==== 15 Debugging: debug print, asserts ====
DEBUG = os.environ.get("LOCAL") is not None   # run locally with LOCAL=1 python3 sol.py


def dbg(*args):
    # Prints to stderr (judges ignore stderr, but it still costs time -- guarded by DEBUG).
    if DEBUG:
        print(*args, file=sys.stderr)
# assert cond, msg -- crashes with RE instead of WA; good for checking invariants. `python3 -O` disables.
# Profile: python3 -m cProfile -s cumtime sol.py < in.txt ; time: /usr/bin/time -v (max RSS = memory).
# TLE checklist: I/O (use buffer read + join), recursion -> iterative, tuples in hot loops,
#   attribute lookups in loops (bind local: ap = a.append), unnecessary copies/slices, CPython vs PyPy.
# MLE checklist: lists of lists/tuples -> flat arrays (idx = i*m + j), array('i'), del big temporaries.
# Edge cases: n=1, n=0/empty, all equal, sorted/reversed, max values, negatives, duplicates, t large.


# ==== 15 Stress testing ====
def stress(solve, brute, gen, iters=1000, seed=0):
    # solve/brute: input -> output; gen(rng) -> random small input. Returns failing case or None.
    rng = random.Random(seed)
    for it in range(iters):
        inp = gen(rng)
        a, b = solve(inp), brute(inp)
        if a != b:
            print(f"MISMATCH at iter {it}\ninput: {inp!r}\nsolve: {a!r}\nbrute: {b!r}", file=sys.stderr)
            return inp
    return None
# Shell version (files sol.py, brute.py, gen.py taking a seed argument):
#   for i in $(seq 1 1000); do
#     python3 gen.py $i > in; python3 sol.py < in > o1; python3 brute.py < in > o2
#     cmp -s o1 o2 || { echo "diff on seed $i"; cat in; break; }
#   done
# Generators: keep sizes tiny (n<=8, values<=10) so failing cases are easy to read.
# Special judge: replace equality with a checker(inp, out) -> bool.


if __name__ == "__main__":
    import io

    # FastIO
    f = FastIO(b"3\n1 2 3\nhello -5\n")
    assert f.int() == 3 and f.ints(3) == [1, 2, 3] and f.str() == "hello" and f.int() == -5 and f.eof()

    # contest_template via redirected stdin/stdout
    class _B:
        def __init__(self, b): self.buffer = io.BytesIO(b)
    old_in, old_out = sys.stdin, sys.stdout
    sys.stdin, sys.stdout = _B(b"2\n3\n1 2 3\nx\n1\n-4\ny\n"), io.StringIO()
    contest_template()
    res = sys.stdout.getvalue()
    sys.stdin, sys.stdout = old_in, old_out
    assert res == "6\n-4\n", res

    # float formatting
    assert fmt_float(1 / 3, 4) == "0.3333" and fmt_float(-1e-12, 3) == "0.000" and fmt_float(-2.5, 1) == "-2.5"

    # interactive: guess hidden number with binary search, judge simulated via precomputed replies
    hidden, n = 37, 100
    # simulate: the judge answers "1" if hidden >= x else "0"; we drive it with a fake stdin per query
    lo, hi, q = 1, n, 0
    while lo < hi:
        mid = (lo + hi + 1) // 2
        o = io.StringIO()
        r = ask(mid, inp=io.StringIO("1\n" if hidden >= mid else "0\n"), outp=o)
        assert o.getvalue() == f"? {mid}\n"
        q += 1
        if r[0] == 1:
            lo = mid
        else:
            hi = mid - 1
    o = io.StringIO()
    answer(lo, outp=o)
    assert lo == hidden and o.getvalue() == "! 37\n" and q <= 7

    # bootstrap: depth 2e5 recursion
    N = 200000
    g = [[] for _ in range(N)]
    for i in range(1, N):
        g[i - 1].append(i)
        g[i].append(i - 1)

    @bootstrap
    def dfs(u, p):
        s = 1
        for v in g[u]:
            if v != p:
                s += yield dfs(v, u)
        yield s
    assert dfs(0, -1) == N

    @bootstrap
    def fact_mod(k):
        if k == 0:
            yield 1
        yield k * (yield fact_mod(k - 1)) % 1000003
    import math
    assert fact_mod(5000) == math.factorial(5000) % 1000003

    # big stack thread
    box = []

    def deep():
        def rec(k):
            return 0 if k == 0 else 1 + rec(k - 1)
        box.append(rec(50000))
    run_with_big_stack(deep)
    assert box == [50000]
    sys.setrecursionlimit(10000)

    # salted counter
    rng = random.Random(1)
    arr = [rng.randrange(10) for _ in range(1000)]
    d = salted_counter(arr)
    for x in range(10):
        assert d.get(x ^ SALT, 0) == arr.count(x)

    # stress harness: correct pair -> None; buggy pair -> failing case found
    gen = lambda r: [r.randint(-5, 5) for _ in range(r.randint(0, 6))]
    def kadane(a):
        best = cur = 0
        for x in a:
            cur = max(0, cur + x)
            best = max(best, cur)
        return best
    def brute_ms(a):
        return max([0] + [sum(a[i:j]) for i in range(len(a)) for j in range(i + 1, len(a) + 1)])
    assert stress(kadane, brute_ms, gen, 500) is None
    old_err, sys.stderr = sys.stderr, io.StringIO()
    bad = stress(lambda a: max(a, default=0) if a else 0, brute_ms, gen, 500)
    sys.stderr = old_err
    assert bad is not None and max(bad, default=0) != brute_ms(bad)

    dbg("not printed unless LOCAL set")
    print("ALL OK")
