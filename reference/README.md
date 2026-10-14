# Python CP Reference

A Python 3 / PyPy team notebook covering every topic in `../competitive-programming-topics.md`.
Each technique is implemented **once**. Other places that need it point to it with `# see <file>: <name>`.
Topics that aren't code (concepts, strategy, theorems) are in this README or in short comment notes inside the files.

Every `.py` file is self-testing. `python3 <file>` checks each snippet against a brute force and prints `ALL OK`.
Run all of them with `./run_tests.sh`.

The printable ICPC team reference document is `trd.pdf` (25 pages, A4). It is built from these files with `python3 build_trd.py --university "Your University"`, which needs Quarto and a LaTeX distribution (TinyTeX works). The same build writes `trd.md`, a plain Markdown edition with no page limit, for AI assistants.

## File map (topic section → file)

| Topics md section | File |
|---|---|
| 1.2–1.4 language, I/O, interactive · 14.3 implementation · 15 debugging/stress testing | `00_core.py` |
| 3 basic techniques (search, sorting, prefix sums, two pointers, greedy, D&C, bits) · 11 offline/queries (Mo, CDQ, parallel BS) · 12 sqrt, randomization, search optimizations, doubling · 13 paradigms | `01_techniques.py` |
| 2 STL → Python mapping · 4 data structures (SortedList, DSU, Fenwick, segment trees, treap, LCT, wavelet…) | `02_data_structures.py` |
| 5 dynamic programming + all DP optimizations | `03_dp.py` |
| 6.1–6.7, 6.10 graphs (traversal, shortest paths, MST, SCC/2-SAT, bridges/BCC, Euler, special) | `04_graphs.py` |
| 6.8 flows · 6.9 matching · 12.4 matroid intersection | `05_flow_matching.py` |
| 7 trees (LCA, HLD, centroid, small-to-large, virtual tree, isomorphism) | `06_trees.py` |
| 8 strings · 4.3 rolling hash | `07_strings.py` |
| 9.1 number theory · 9.2 combinatorics · 9.3 probability · 9.6 game theory · 9.7 misc math | `08_math_nt_comb.py` |
| 9.4 linear algebra · 9.5 polynomials & transforms | `09_algebra_poly.py` |
| 10 geometry · 3.7 closest pair | `10_geometry.py` |
| 1.1, 1.3, 13, 14, 16–21 (non-code) | this README |

To find a topic, search a file for its number: `grep -n "==== 4.6" reference/*.py`.

---

## 1. Foundations

### Formats & verdicts (1.1)
- **ICPC**: team of 3, one computer, rank by solved count then penalty (minutes + 20 per wrong attempt on solved problems).
- **IOI**: individual, subtasks with partial points, no penalty. Go for subtasks first.
- **Codeforces**: pretests during the round, system tests after. You can hack others (Div1/2). Score decays with time.
- **AtCoder**: ABC/ARC/AGC, penalty 5 min per WA. Library Checker / ACL are reference implementations.
- Verdicts: AC, WA, TLE, MLE, RE (in Python: exceptions, RecursionError, MemoryError), CE, PE (presentation; usually just whitespace).
- Ratings are Elo-like. CF colors: gray < 1200, green, cyan 1400, blue 1600, violet 1900, orange 2100, red 2400+.

### Constraints → intended complexity (1.3)
CPython runs roughly **10^7 simple ops/sec**; PyPy is about **10^8** (closer to C++ / 3–10). Budget accordingly.

| n | Intended | Python note |
|---|---|---|
| ≤ 10–11 | O(n!) | ok |
| ≤ 20 | O(2^n · n) | bitmask DP needs PyPy for n=20 |
| ≤ 40 | O(2^(n/2)) meet in the middle | |
| ≤ 500 | O(n³) | 1.25·10^8: PyPy only, or vectorize with ints/bitsets |
| ≤ 5000 | O(n²) | PyPy; CPython only with list-comprehension inner loops |
| ≤ 10^5–2·10^5 | O(n log n), O(n√n) | √ solutions are tight in Python; prefer log |
| ≤ 10^6 | O(n), O(n log n) w/ small constant | avoid per-element function calls |
| ≤ 10^9 | O(√n) | |
| ≤ 10^18 | O(log n), O(1), O(n^{1/4}) (Pollard) | Python big ints are free here |

- **Amortized analysis** (potential/accounting): monotonic stack/deque, DSU, dynamic arrays, splay trees. Total cost is bounded even when one operation is slow.
- **Expected vs. worst case**: treaps, quickselect, and hashing are expected-time. Shuffle the input or randomize seeds so adversarial tests can't hit the worst case.
- Memory: a Python int is about 28 B and a list slot 8 B. 10^7 ints in a list is roughly 300 MB. Use `array('i')`, `bytearray`, or flat lists instead of nested lists.

---

## 2. Python performance checklist (1.2, 15)
1. Submit with **PyPy** when the judge offers it (except big-int-heavy or recursion-heavy code: CPython is faster there).
2. Read all input at once: `sys.stdin.buffer.read().split()`. Write once: `sys.stdout.write("\n".join(map(str, out)))`.
3. Put the code in a `main()` function. Local variable lookup is faster than global.
4. No recursion deeper than ~10^4. Convert to an iterative stack, or use the `bootstrap` decorator in `00_core.py`.
5. Flatten 2D arrays: `a[i*m+j]`. Precompute `append = lst.append` in hot loops.
6. Avoid `dict` keys that are tuples in hot loops. Encode them as `i*K+j`.
7. Python ints are bitsets: `x |= x << a`, `bin(x).count("1")` / `x.bit_count()` gives an n/64-style speedup for free.
8. Anti-hash: CF hackers break `dict`/`set` of ints. XOR keys with a random salt (see `00_core.py`).
9. `sort` with a `key` is fast. `cmp_to_key` is about 10x slower; pack a tuple key or an int key instead.
10. `heapq` is a min-heap. Push `-x` for a max-heap, and tuples `(priority, id)` to avoid comparing payloads.
11. Big-int multiplication in CPython is Karatsuba, not FFT. For large convolutions use `convolve` in `09_algebra_poly.py` (it switches to `decimal`, whose libmpdec uses NTT: ~0.7s at n=1e5). On PyPy use `ntt_convolve`.

### Python equivalents of C++ language essentials (2.1)
- **Integers**: arbitrary precision, so there is no overflow and no `__int128`. Overflow bugs become *speed* bugs: reduce `% MOD` early so numbers stay small. To emulate unsigned 64-bit: `x & ((1 << 64) - 1)`.
- **Floats**: `float` is a C double. For more precision use `fractions.Fraction` (exact) or `decimal.Decimal` with `getcontext().prec = 50`. Use `math.isclose`. Prefer `math.isqrt(n)` over `int(n**0.5)` (the latter is wrong for n > 2^52).
- **Bit builtins**: `x.bit_count()` (popcount, 3.10+), `x.bit_length()` (= 64 − clz; `__lg(x) = x.bit_length()-1`), `(x & -x).bit_length()-1` (ctz).
- **Lambdas / recursive lambdas**: nested `def` closures. `nonlocal` for mutation. Use `functools.cache` for memoization.
- **Template**: see `00_core.py`. No macros or pragmas; the equivalent of `-O2` is choosing PyPy.
- **C++-only debug tools** (sanitizers, `_GLIBCXX_DEBUG`, gdb) map to: Python already bounds-checks lists (watch out: negative indices wrap silently!), `python3 -X dev`, `faulthandler`, `pdb` / `breakpoint()`.

---

## 3. Theorems & facts cheat sheet (no code needed)
- **Master theorem**: T(n)=aT(n/b)+O(n^d). If d > log_b a → O(n^d). If d = log_b a → O(n^d log n). Otherwise O(n^{log_b a}).
- **Matroid / greedy**: greedy by weight is optimal for every weight function ⇔ the feasible sets form a matroid (graphic → Kruskal, partition, linear → XOR basis).
- **Exchange argument**: assume an optimal order. Swapping an adjacent pair that violates your order doesn't make it worse. So sorted order is optimal.
- **König**: in a bipartite graph, max matching = min vertex cover. Max independent set = V − max matching.
- **Hall**: a perfect matching of side L exists ⇔ |N(S)| ≥ |S| for every S ⊆ L.
- **Dilworth**: min chains covering a poset = max antichain. Min path cover of a DAG = V − max matching (split graph).
- **Max-flow = min-cut**. Edge connectivity = min over t of maxflow(s,t). Vertex connectivity needs vertex splitting.
- **Euler path**: undirected needs 0 or 2 odd vertices. Directed needs in=out everywhere except start (+1 out) and end (+1 in). The graph must also be connected.
- **Euler's formula (planar)**: V − E + F = 1 + C. Simple planar ⇒ E ≤ 3V − 6 (V ≥ 3).
- **Tournament**: always has a Hamiltonian path. Strongly connected ⇒ Hamiltonian cycle.
- **Cayley**: n^{n−2} labeled trees. A forest of k trees rooted at given k vertices: k·n^{n−k−1}.
- **Matrix-tree**: #spanning trees = any cofactor of the Laplacian.
- **Pick**: A = I + B/2 − 1.
- **Sprague–Grundy**: game sum ⇒ XOR of Grundy values. Position is losing ⇔ value 0.
- **Burnside**: #orbits = (1/|G|) Σ_g |Fix(g)|.
- **Lucas**: C(n,k) mod p = Π C(n_i, k_i) over base-p digits. **Kummer**: the exponent of p in C(m+n, m) equals the number of carries when adding m and n in base p.
- **Wilson**: (p−1)! ≡ −1 (mod p). **Fermat**: a^{p−1} ≡ 1. **Euler**: a^φ(m) ≡ 1 when gcd(a,m)=1. Generalized: a^k ≡ a^{k mod φ + φ} for k ≥ log2 m.
- **Erdős–Szekeres / LIS duality**: min number of non-increasing subsequences covering the sequence = LIS length.
- **Linearity of expectation** holds without independence. **Geometric distribution**: with success prob p, expected trials = 1/p.
- **Catalan** C_n = C(2n,n)/(n+1): bracket sequences, binary trees, triangulations, monotone lattice paths under the diagonal.
- Number of divisors of n ≤ 10^18 is at most 103680. For n ≤ 10^9 it's at most 1344.
- Prime gaps below 10^18 are at most ~1500. π(10^6)=78498, π(10^7)=664579.
- Useful primes: 998244353 = 119·2^23+1 (g=3), 10^9+7, 2^61−1 (hashing).

---

## 4. Problem-solving playbook (13, 14)
**Approach (14.1)**
1. Read the constraints first. They reveal the complexity (table above).
2. Hand-solve the samples. Then make up your own tiny cases.
3. Brute force small n and print the answers. Look for a pattern or check OEIS.
4. Simplify: special cases (all equal, n tiny, a tree instead of a graph, a single query, no updates).
5. Reduce to a known problem: shortest path on a state graph, flow, matching, 2-SAT, LIS, knapsack.
6. Change perspective: count the complement, sum contributions per element/pair/edge, think backward (reverse the process or the deletions), fix one variable and iterate over the rest, sort by one key and sweep.
7. Binary search the answer whenever "is answer ≤ X?" is monotone.

**Constructive (13.1)**: try small n by brute force and look for a pattern. Use induction (build n from n−1) and invariants to prove impossibility. Common building blocks are zig-zags, alternating parity, powers of 2, and stars/paths.
**Ad hoc (13.3)**: look for invariants (sum mod k, parity of inversions, XOR), monovariants (something that strictly decreases, so the process terminates), and symmetry/pairing strategies.
**Interactive (13.2)**: flush after every query. The query limit hints at the method (log2 n → binary search; ⌈log3⌉ → ternary weighing). Adaptive judges mean you can't rely on randomness alone. Use information-theoretic lower bounds to sanity-check.
**Heuristic / output-only (13.4)**: get a greedy baseline first, then local search / simulated annealing (template in `01_techniques.py`), then beam search. Watch the clock with `time.perf_counter()`.
**Communication/two-run (13.5)**: design an encoding. Count bits with `log2(#states)`.

**Proofs (14.2)**: use contradiction, induction, an exchange argument, or an invariant. Bound the answer from both sides (construct something achieving X, then prove nothing does better).

**Implementation (14.3)**: reset all per-test state, but only up to the size n of that test (don't reallocate 10^6 arrays per test). Write small functions. Don't mutate globals across tests.

---

## 5. Debugging & testing (15)
- **Common bugs**: off-by-one in binary search bounds; forgetting to reset between tests; using `[[0]*m]*n` (aliasing); integer division `//` vs `/` with negatives (`-7//2 == -4`); `%` of a negative is non-negative in Python (unlike C++); a float `1e18` used as an int; recursion depth; mutating a list while iterating over it.
- **Edge cases**: n=1, n=0/empty, all equal, sorted/reverse sorted, max values, negatives, duplicates, disconnected graphs, self-loops/multi-edges, a star tree / a path tree (depth!).
- **Stress test**: the `stress()` harness in `00_core.py`, or a shell loop: `for i in $(seq 1000); do python3 gen.py $i > in; diff <(python3 sol.py < in) <(python3 brute.py < in) || break; done`.
- **Debug output**: `print(..., file=sys.stderr)`, `assert`s for invariants.
- **Profile**: `python3 -m cProfile -s cumtime sol.py < in`. Time it with `time python3 sol.py < big.in`.
- **TLE**: I/O, recursion, attribute lookups in loops, tuple allocation, dict vs list indexing. Switch to PyPy.
- **MLE**: nested lists of ints, storing all states, keeping big strings. Use `array`, `bytearray`, generators.
- **RecursionError / segfault**: the default limit is 1000. Raising it without `threading.stack_size` crashes CPython at around 10^5 depth.

---

## 6. Contest strategy (16)
**Individual**: skim all problems in the first minutes. Solve in order of (confidence/time). Check the scoreboard for which problems are actually easy. If you're stuck 20–30 min with no new idea, switch problems. Before submitting, test the max case for TLE and edge cases for WA. A WA costs ~10–20 min of penalty, so 1 minute of testing is cheap. CF hacks: look for unsafe hashes, missing overflow handling (rare in Python), and greedy flaws. IOI: grab every subtask first, then go for full solutions.

**Team (ICPC)**: one person types while the others think on paper. The fastest typist starts on the easiest problem immediately. Print code and debug on paper rather than blocking the machine. Hand off a problem with a 1-minute verbal summary + constraints. Follow the scoreboard. The team notebook is this folder (printed): keep it short, indexed, and tested.

**Mental (16.3)**: do timed virtual contests weekly and always upsolve. When a problem tilts you, move on and come back with a fresh head. Consistency beats cramming; sleep before contests.

---

## 7. Tools & environment (17)
- **Judges**: Codeforces, AtCoder, CSES, Kattis, CodeChef, LeetCode, SPOJ, UVa, Timus, LibreOJ, Luogu, yukicoder, DMOJ, Project Euler. **Library Checker** (judge.yosupo.jp) is ideal for verifying the snippets here.
- **Helpers**: Competitive Companion (browser → parse samples), `oj` (online-judge-tools: `oj d URL`, `oj t -c "python3 a.py"`), cf-tool, CP Editor.
- **Libraries to compare against**: AtCoder Library (ACL, and its Python port `ac-library-python`), KACTL, ei1333, Nyaan, PyRival (Python-specific CP library).
- **Problem-setting tools**: Polygon (CF) and testlib (generators, validators, checkers).

## 8. Competitions (18)
- **High school**: IOI, regional olympiads (APIO, BOI, CEOI, EGOI, JOI, NOI, EJOI), national olympiads (USACO, CCO, COCI…), selection camps.
- **University**: ICPC (regionals → continental → World Finals).
- **Open**: CF (Div 1–4, Educational, Global), AtCoder (ABC/ARC/AGC/AHC heuristic), CodeChef Starters, Meta Hacker Cup, LeetCode weekly, Advent of Code, Universal Cup. Google Code Jam / Kick Start / Hash Code and TCO are discontinued; their archives remain.

## 9. Problem setting (19)
Write an unambiguous statement with formal constraints. Subtasks should map to distinct ideas. Tests should include random cases, max cases, edge cases, anti-hash cases (Thue–Morse strings against 2^64 hashes), anti-greedy cases, and anti-SPFA cases (grid graphs). Write a validator (testlib), a checker for multiple answers, the model solution, and the wrong/slow solutions that must fail. Testers should solve it independently. Write the editorial with intuition first, then the proof, then the complexity.

## 10. Practice & learning (20)
- **Problem sets**: CSES Problem Set (classic coverage), AtCoder Educational DP Contest, CF EDU section, USACO Guide, Kattis. For interview-adjacent practice: NeetCode / Blind 75.
- **Method**: upsolve every contest. Read the editorial only after a real attempt, then implement it without looking. Practice ~200 rating points above your current rating.
- **Books**: *Competitive Programmer's Handbook* (Laaksonen), *CP4* (Halim), CLRS, Sedgewick, *Guide to CP*, Kleinberg–Tardos.
- **Websites**: cp-algorithms.com, USACO Guide, OI Wiki, CF blogs. **Video**: Errichto, SecondThread, William Lin, Colin Galen.

## 11. Career (21)
CP skills transfer directly to coding interviews (LeetCode-style is a subset of Div2 A–D). They also carry over to quant trading/research (probability, speed under pressure), performance engineering, and algorithms research. Coaching and problem setting are a common path for staying involved. Balance CP with system design and real-world engineering (testing, code review, large codebases).
