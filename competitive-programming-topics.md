# Competitive Programming Topics: A Comprehensive List

A structured map of everything you can learn in competitive programming, roughly ordered from fundamentals to advanced and specialized topics.

---

## 1. Foundations

### 1.1 What Competitive Programming Is
- Definition: solving well-defined algorithmic problems under time and memory limits
- History: ACM ICPC, IOI, TopCoder, the rise of Codeforces and AtCoder
- Contest formats: ICPC-style (team, penalty time), IOI-style (individual, subtasks, partial scoring), Codeforces rounds, AtCoder ABC/ARC/AGC, TopCoder SRM, Google Code Jam / Kick Start (legacy), Meta Hacker Cup
- Judging: verdicts (AC, WA, TLE, MLE, RE, CE, PE), hidden tests, pretests and system tests, hacking
- Ratings and rating systems (Elo-like, Codeforces colors, AtCoder colors)
- Competitive programming vs. interview prep vs. research algorithms

### 1.2 Choosing a Language
- C++ (the dominant choice): speed, STL, templates
- Java, Python, Kotlin, Rust, and their trade-offs
- Python performance pitfalls (PyPy, recursion limits, fast I/O)
- Compiler flags and judge environments (`-O2`, `-std=c++20`, 32 vs. 64-bit)
- Language-specific undefined behavior and gotchas

### 1.3 Complexity Analysis
- Big-O, Big-Θ, Big-Ω
- Estimating operations per second (~10^8 simple ops/sec rule of thumb)
- Reading constraints to infer the intended complexity (n ≤ 20 → 2^n, n ≤ 500 → n^3, n ≤ 10^5 → n log n, n ≤ 10^18 → log n or O(1))
- Amortized analysis (potential method, accounting method)
- Space complexity and memory limits
- Constant factors and cache effects
- Expected vs. worst-case complexity for randomized algorithms

### 1.4 Input/Output
- Fast I/O (`ios::sync_with_stdio(false)`, `cin.tie(nullptr)`, `scanf`/`printf`, custom readers)
- Reading until EOF, multi-test-case formats
- Interactive problems and flushing output
- Output precision for floating point
- Special judges (multiple valid answers)

---

## 2. C++ and the STL for Contests

### 2.1 Language Essentials
- Integer types and overflow (`int`, `long long`, `__int128`, unsigned)
- Floating types (`double`, `long double`) and precision
- Bit manipulation operators and builtins (`__builtin_popcount`, `__builtin_clz`, `__builtin_ctz`, `__lg`)
- Lambdas, recursive lambdas, `auto`, structured bindings
- Macros, templates, and a personal contest template
- Pragmas and optimization hints (`#pragma GCC optimize`, `target("avx2")`)

### 2.2 STL Containers
- `vector`, `array`, `deque`, `list`
- `stack`, `queue`, `priority_queue` (min-heap via `greater<>`)
- `set`, `multiset`, `map`, `multimap`
- `unordered_set`, `unordered_map` and anti-hash tests (custom hashes)
- `bitset` and bitset-optimized algorithms
- `pair`, `tuple`, `string`, `string_view`

### 2.3 STL Algorithms
- `sort`, `stable_sort`, custom comparators and strict weak ordering
- `lower_bound`, `upper_bound`, `equal_range`, `binary_search`
- `next_permutation`, `prev_permutation`
- `unique`, `reverse`, `rotate`, `nth_element`, `partial_sort`
- `accumulate`, `partial_sum`, `iota`, `max_element`, `min_element`
- Policy-based data structures (`__gnu_pbds` ordered set / order statistics tree)

---

## 3. Basic Techniques

### 3.1 Brute Force and Complete Search
- Iterating over all possibilities
- Generating subsets (bitmasks, recursion)
- Generating permutations
- Backtracking and pruning
- Meet in the middle

### 3.2 Simulation and Implementation
- Carefully translating statements into code
- Grid problems and direction arrays
- Handling edge cases and off-by-one errors
- Date/time and calendar problems
- Parsing strings and expressions

### 3.3 Sorting and Searching
- Sorting-based reasoning and exchange arguments
- Coordinate compression
- Counting sort, radix sort
- Binary search on sorted arrays
- Binary search on the answer (monotonic predicates)
- Parallel binary search
- Ternary search on unimodal functions
- Fractional binary search and fixed-iteration loops

### 3.4 Prefix Sums and Difference Arrays
- 1D and 2D prefix sums
- Difference arrays and range updates
- Prefix XOR, prefix products
- Higher-order difference arrays (adding polynomials to ranges)

### 3.5 Two Pointers and Sliding Window
- Two pointers on sorted arrays
- Variable-size sliding window
- Sliding window minimum/maximum with a monotonic deque
- Fast and slow pointers (cycle detection, Floyd's algorithm)

### 3.6 Greedy Algorithms
- Greedy choice property and proving greedy correctness
- Exchange arguments
- Interval scheduling, interval covering, interval partitioning
- Activity selection, job sequencing with deadlines
- Huffman coding
- Greedy with priority queues
- Regret greedy (greedy with "undo")
- Matroids and when greedy is optimal

### 3.7 Recursion and Divide and Conquer
- Recursive thinking and base cases
- Merge sort and counting inversions
- Quickselect
- Closest pair of points
- Master theorem
- Karatsuba multiplication

### 3.8 Bit Manipulation
- Bit tricks (lowest set bit, toggling, checking, clearing)
- Iterating over submasks
- XOR properties and tricks
- Gray code
- Bitmask enumeration

---

## 4. Data Structures

### 4.1 Linear Structures
- Stacks and monotonic stacks (next greater element, largest rectangle in histogram)
- Queues and monotonic queues
- Deques
- Linked lists (and XOR linked lists)

### 4.2 Heaps and Priority Queues
- Binary heaps
- Heap operations and heapify
- Two-heap technique (running median)
- Mergeable heaps: leftist heaps, skew heaps, pairing heaps, binomial heaps
- Lazy deletion in heaps

### 4.3 Hashing
- Hash tables and collision handling
- Rolling hash / polynomial hashing
- Anti-hash tests and randomized bases
- Zobrist hashing
- Hashing sets and multisets (XOR hashing, sum hashing)

### 4.4 Disjoint Set Union (Union-Find)
- Path compression and union by rank/size
- Inverse Ackermann complexity
- DSU with additional data (sizes, parity, potentials/weights)
- DSU with rollback (for offline dynamic connectivity)
- Persistent DSU

### 4.5 Range Query Structures
- Sparse table (RMQ, idempotent operations)
- Disjoint sparse table
- Square root decomposition
- Fenwick tree (Binary Indexed Tree): point update/range query, range update/point query, range update/range query
- 2D Fenwick tree
- Fenwick tree binary lifting (k-th element)

### 4.6 Segment Trees
- Basic segment tree (point update, range query)
- Lazy propagation (range update, range query)
- Iterative (bottom-up) segment tree
- Segment tree on custom monoids (merging nodes with structs)
- Segment tree descent / binary search on segment tree
- Merge sort tree
- Persistent segment tree (k-th smallest in range, versioning)
- Dynamic / implicit (sparse) segment tree
- 2D segment trees and segment tree of BITs
- Segment tree beats (Ji Driver segment tree)
- Li Chao tree
- Segment tree merging and splitting
- Kinetic segment tree

### 4.7 Balanced Trees
- Binary search trees
- Treaps (Cartesian trees) and implicit treaps (array operations: split, merge, reverse)
- Splay trees
- AVL and red-black trees (conceptual)
- Scapegoat trees
- Skip lists
- Order statistics trees

### 4.8 Advanced and Specialized Structures
- Sqrt tree
- Wavelet tree / wavelet matrix
- Van Emde Boas tree (conceptual)
- Link-cut trees
- Euler tour trees
- Top trees (conceptual)
- Persistent data structures (path copying, fat nodes)
- Retroactive data structures (conceptual)
- Rope data structure
- Cartesian tree and its relation to RMQ/LCA
- Bitset-based structures
- Fractional cascading

---

## 5. Dynamic Programming

### 5.1 Fundamentals
- Overlapping subproblems and optimal substructure
- Memoization (top-down) vs. tabulation (bottom-up)
- Defining states and transitions
- Reconstructing solutions
- Space optimization (rolling arrays)
- Pull vs. push DP

### 5.2 Classic DP Problems
- Fibonacci-style recurrences and counting paths
- Coin change (counting and minimizing)
- 0/1 knapsack, unbounded knapsack, bounded knapsack (binary splitting, monotonic queue)
- Longest increasing subsequence (O(n²) and O(n log n))
- Longest common subsequence
- Edit distance
- Matrix chain multiplication
- Subset sum (with bitset optimization)
- Partition problems
- Grid DP (paths, obstacles)

### 5.3 DP Variants
- Interval DP
- Tree DP (subtree DP, rerooting / all-roots DP)
- DAG DP
- Bitmask DP (TSP, assignment, Hamiltonian paths)
- Digit DP
- Probability and expected value DP
- Game DP
- DP on subsequences and strings
- Profile DP / broken profile DP (tiling)
- DP over subsets: Sum over Subsets (SOS DP), zeta/Möbius transforms
- DP on permutations (insertion DP, "connected components" DP)
- DP with data structures (segment tree / BIT transitions)
- Counting DP with inclusion–exclusion
- Plug DP (connectivity-based profile DP)

### 5.4 DP Optimizations
- Prefix sum optimization
- Monotonic queue optimization
- Convex Hull Trick (CHT)
- Li Chao tree for DP
- Divide and conquer optimization
- Knuth's optimization (Knuth–Yao quadrangle inequality)
- Aliens trick / WQS binary search / Lagrangian relaxation
- Slope trick
- Matrix exponentiation for linear recurrences
- SMAWK algorithm and monotone matrices (conceptual)
- Bitset DP speedups
- Open and close interval trick

---

## 6. Graph Theory

### 6.1 Graph Basics
- Representations: adjacency list, adjacency matrix, edge list
- Directed, undirected, weighted, unweighted graphs
- Degrees, paths, cycles, connectivity
- Implicit graphs (grids, state spaces)
- Bipartite graphs and bipartiteness check

### 6.2 Graph Traversal
- Depth-first search (DFS) and DFS tree
- Breadth-first search (BFS)
- Multi-source BFS
- 0-1 BFS
- Flood fill
- Connected components
- Cycle detection (directed and undirected)
- Entry/exit times and edge classification

### 6.3 Shortest Paths
- BFS for unweighted graphs
- Dijkstra's algorithm (with heap, O(n²) dense version)
- Bellman–Ford and negative cycle detection
- SPFA (and why it can be hacked)
- Floyd–Warshall (all pairs, transitive closure)
- Johnson's algorithm
- Shortest paths on DAGs
- K-th shortest path (Yen's algorithm, Eppstein's algorithm)
- Shortest path with state augmentation (layered graphs)
- Difference constraints systems

### 6.4 Minimum Spanning Trees
- Kruskal's algorithm
- Prim's algorithm
- Borůvka's algorithm
- Second-best MST
- Minimum bottleneck spanning tree
- Kruskal reconstruction tree
- Directed MST / minimum arborescence (Chu–Liu/Edmonds)
- Steiner tree (bitmask DP)

### 6.5 Directed Graphs
- Topological sorting (Kahn's algorithm, DFS-based)
- Strongly connected components (Tarjan, Kosaraju)
- Condensation graph
- 2-SAT
- Functional graphs and cycle finding
- Dominator trees

### 6.6 Connectivity
- Bridges and articulation points
- 2-edge-connected and 2-vertex-connected (biconnected) components
- Block-cut tree
- Edge and vertex connectivity
- Offline dynamic connectivity
- Online dynamic connectivity (conceptual)

### 6.7 Eulerian and Hamiltonian Paths
- Eulerian paths and circuits (Hierholzer's algorithm)
- De Bruijn sequences
- Hamiltonian paths (bitmask DP)
- Chinese postman problem

### 6.8 Network Flow
- Max flow / min cut theorem
- Ford–Fulkerson and Edmonds–Karp
- Dinic's algorithm
- Push–relabel
- Min cost max flow (SPFA/Bellman–Ford, Dijkstra with potentials)
- Flows with lower bounds (circulations)
- Modeling problems as flows: vertex splitting, project selection / closure problem
- Min cut applications
- Gomory–Hu trees
- Global min cut (Stoer–Wagner)

### 6.9 Matching
- Bipartite matching (Kuhn's algorithm, Hopcroft–Karp)
- König's theorem, minimum vertex cover, maximum independent set in bipartite graphs
- Hall's marriage theorem
- Weighted bipartite matching (Hungarian algorithm)
- General matching (Edmonds' blossom algorithm)
- Stable marriage (Gale–Shapley)
- Minimum path cover in DAGs
- Dilworth's theorem

### 6.10 Special Graph Topics
- Graph coloring (and chromatic number for small n)
- Maximum clique / independent set (Bron–Kerbosch, meet in the middle)
- Planar graphs and Euler's formula
- Tournament graphs
- Cactus graphs
- Chordal graphs (conceptual)
- Counting spanning trees (Kirchhoff's matrix-tree theorem)
- Prüfer sequences
- Graph isomorphism heuristics and tree isomorphism

---

## 7. Trees

### 7.1 Tree Basics
- Rooted vs. unrooted trees, parent/child, depth, height
- Tree traversals (preorder, inorder, postorder, level order)
- Subtree sizes
- Diameter and center of a tree
- Centroid of a tree

### 7.2 Lowest Common Ancestor
- Binary lifting
- Euler tour + RMQ
- Tarjan's offline LCA
- Farach-Colton–Bender (conceptual)
- Distance queries and k-th ancestor queries

### 7.3 Tree Techniques
- Euler tour / flattening a tree into an array
- Small-to-large merging (DSU on tree, sack technique)
- Heavy-light decomposition (HLD)
- Centroid decomposition
- Rerooting DP
- Virtual trees (auxiliary trees)
- Tree hashing and isomorphism
- Link-cut trees for dynamic trees
- Mo's algorithm on trees
- Long path decomposition

---

## 8. Strings

### 8.1 String Basics
- String representation and manipulation
- Lexicographic comparison
- Palindromes
- Subsequences vs. substrings

### 8.2 Pattern Matching
- Naive matching
- Prefix function (KMP)
- Z-function
- Rabin–Karp (rolling hash)
- Aho–Corasick automaton (multiple patterns)
- Bitset-based matching (Shift-And)
- Wildcard matching with FFT

### 8.3 Advanced String Structures
- Tries (and binary tries for XOR problems)
- Suffix arrays (O(n log n) construction, SA-IS conceptual)
- LCP array (Kasai's algorithm)
- Suffix automaton
- Suffix trees (Ukkonen's algorithm)
- Palindromic tree (eertree)
- Manacher's algorithm
- Lyndon factorization (Duval's algorithm)
- Minimal string rotation (Booth's algorithm)
- Run enumeration (conceptual)

### 8.4 String Problems
- Counting distinct substrings
- Longest repeated substring
- Longest common substring of multiple strings
- String periods and borders
- Hashing-based string comparison and binary search on LCP

---

## 9. Mathematics

### 9.1 Number Theory
- Divisibility, divisors, and prime factorization
- Sieve of Eratosthenes, linear sieve, segmented sieve
- GCD and LCM, Euclidean algorithm
- Extended Euclidean algorithm and Bézout's identity
- Linear Diophantine equations
- Modular arithmetic
- Modular inverse (Fermat, extended Euclid)
- Fast (binary) exponentiation
- Euler's totient function and Euler's theorem
- Chinese Remainder Theorem (and generalized CRT)
- Wilson's theorem, Lucas's theorem
- Primitive roots and discrete logarithm (baby-step giant-step)
- Discrete root
- Quadratic residues and Tonelli–Shanks
- Miller–Rabin primality test
- Pollard's rho factorization
- Multiplicative functions
- Möbius function and Möbius inversion
- Dirichlet convolution
- Prefix sums of multiplicative functions (Du's sieve, Min_25 sieve, Lucy_Hedgehog)
- Divisor sum tricks and floor division blocks (√n distinct values of ⌊n/i⌋)
- Legendre's formula
- Continued fractions and Stern–Brocot tree
- Pell's equation (conceptual)

### 9.2 Combinatorics
- Counting principles (sum and product rules)
- Permutations and combinations
- Binomial coefficients and Pascal's triangle
- Precomputing factorials and inverse factorials
- Stars and bars
- Inclusion–exclusion principle
- Pigeonhole principle
- Catalan numbers
- Stirling numbers (first and second kind)
- Bell numbers
- Derangements
- Partition numbers
- Burnside's lemma and Pólya enumeration theorem
- Generating functions (ordinary and exponential)
- Lindström–Gessel–Viennot lemma
- Bijective proofs
- Counting labeled graphs and trees (Cayley's formula)
- Young tableaux and the hook length formula

### 9.3 Probability and Expected Value
- Probability basics and conditional probability
- Linearity of expectation
- Expected value DP
- Markov chains
- Random walks
- Geometric distributions and "expected time until" problems
- Contribution technique

### 9.4 Linear Algebra
- Matrices and matrix multiplication
- Matrix exponentiation (linear recurrences, path counting)
- Gaussian elimination (real, modular, and over GF(2))
- Determinants
- Rank and inverse
- XOR basis (linear basis over GF(2))
- Berlekamp–Massey algorithm
- Systems of linear equations

### 9.5 Polynomials and Transforms
- Fast Fourier Transform (FFT)
- Number Theoretic Transform (NTT)
- Polynomial multiplication with arbitrary modulus (three-modulus NTT, splitting FFT)
- Polynomial inverse, log, exp, sqrt, power (Newton's method)
- Polynomial division and multipoint evaluation
- Lagrange interpolation
- Walsh–Hadamard transform (XOR convolution)
- AND/OR convolutions
- Subset convolution
- Online convolution (relaxed / CDQ-style)
- Bostan–Mori algorithm for linear recurrences

### 9.6 Game Theory
- Impartial games
- Winning and losing positions
- Nim and Nim-sum
- Sprague–Grundy theorem and Grundy numbers
- Variants: Misère Nim, staircase Nim, Moore's Nim, Wythoff's game
- Partisan games (conceptual)
- Minimax and alpha-beta pruning

### 9.7 Miscellaneous Math
- Arithmetic and geometric series
- Fibonacci properties (Zeckendorf representation, Pisano periods)
- Big integer arithmetic
- Numerical methods (bisection, Newton's method, Simpson's rule)
- Linear programming and the simplex method
- Fraction handling and exact arithmetic
- Base conversions and numeral systems
- Floor sum (sum of ⌊(ai+b)/c⌋)

---

## 10. Computational Geometry

### 10.1 Basics
- Points, vectors, and lines
- Dot product and cross product
- Orientation tests
- Precision and epsilon comparisons; integer geometry where possible
- Distance computations (point–point, point–line, point–segment)
- Angles and polar sorting

### 10.2 Lines and Segments
- Line representation forms
- Line–line intersection
- Segment intersection
- Checking if a point is on a segment
- Sweep line for segment intersections (Bentley–Ottmann)

### 10.3 Polygons
- Polygon area (shoelace formula)
- Pick's theorem
- Point in polygon (ray casting, winding number)
- Convex polygons and point-in-convex-polygon in O(log n)
- Polygon cutting and clipping
- Polygon triangulation (conceptual)

### 10.4 Convex Hull
- Graham scan
- Andrew's monotone chain
- Rotating calipers (diameter, width, farthest pair)
- Minkowski sum
- Dynamic convex hull
- Half-plane intersection
- 3D convex hull

### 10.5 Circles
- Circle–line and circle–circle intersection
- Tangents
- Minimum enclosing circle (Welzl's algorithm)
- Area of union of circles

### 10.6 Advanced Geometry
- Closest pair of points
- Sweep line techniques
- Voronoi diagrams and Delaunay triangulation
- Point location
- Manhattan MST
- Nearest neighbor and k-d trees
- Pick's theorem and lattice point counting
- Rotational sweep
- Geometric duality

---

## 11. Query Processing and Offline Techniques

- Offline vs. online queries
- Sorting queries
- Mo's algorithm
- Mo's algorithm with updates
- Hilbert order for Mo's algorithm
- Sqrt decomposition on queries
- Divide and conquer on queries (CDQ divide and conquer)
- Offline divide and conquer on time (segment tree over time)
- Parallel binary search
- Sweep line with data structures
- Persistent structures for online answering
- Answering queries in reverse (deleting → adding)

---

## 12. Advanced Algorithmic Techniques

### 12.1 Sqrt Techniques
- Sqrt decomposition on arrays
- Heavy/light splitting by degree or frequency
- Sqrt decomposition on queries (rebuilding every √q)
- Block-based data structures

### 12.2 Randomization
- Randomized algorithms (Las Vegas vs. Monte Carlo)
- Random shuffling for expected complexity
- Hashing with random bases
- Randomized incremental construction
- Color coding
- Random sampling and probabilistic counting
- Simulated annealing and hill climbing (for heuristic problems)

### 12.3 Search Optimizations
- Meet in the middle
- Iterative deepening (IDDFS, IDA*)
- A* search
- Branch and bound
- Bidirectional BFS
- Pruning strategies

### 12.4 Other Advanced Techniques
- Small-to-large merging
- Binary lifting (beyond LCA)
- Doubling techniques
- Matroid intersection
- Fixed-parameter tractable algorithms (conceptual)
- Bitset speedups (n³/64)
- Four Russians technique
- Offline-to-online conversion
- Constructive algorithms

---

## 13. Problem Types and Paradigms

### 13.1 Constructive Problems
- Building arrays, permutations, graphs, or strings that satisfy conditions
- Proving impossibility
- Small-case brute force to discover patterns
- Invariants and monovariants

### 13.2 Interactive Problems
- Query-limited reasoning
- Binary search via queries
- Adaptive vs. non-adaptive judges
- Information-theoretic lower bounds

### 13.3 Ad Hoc and Observation-Based Problems
- Finding invariants
- Parity arguments
- Symmetry
- Reformulation and reduction to known problems
- Thinking in terms of contribution

### 13.4 Output-Only, Approximation, and Heuristic Problems
- Output-only tasks (IOI)
- Optimization contests (AtCoder Heuristic Contest, Topcoder Marathon)
- Local search, simulated annealing, beam search, genetic algorithms
- Scoring-based evaluation

### 13.5 Other Formats
- Communication / two-run problems (IOI-style)
- Reactive problems
- Grid puzzles and board games
- Partial-score / subtask strategy

---

## 14. Problem-Solving Skills

### 14.1 Approaching a Problem
- Reading statements carefully and identifying constraints
- Simplifying the problem (special cases, smaller constraints)
- Working through examples by hand
- Brute-forcing small cases and looking for patterns (OEIS)
- Reducing to a known problem
- Thinking backward
- Fixing one variable and iterating over others
- Changing perspective (complement counting, contribution, transform the graph)

### 14.2 Proofs and Correctness
- Proof by contradiction, induction, exchange argument
- Invariants and monovariants
- Proving greedy and DP correctness
- Upper and lower bounds

### 14.3 Implementation Skills
- Writing clean, bug-resistant code quickly
- Choosing appropriate data types
- Modular code and reusable snippets
- Avoiding global state bugs across test cases
- Handling multiple test cases with resets

---

## 15. Debugging and Testing

- Common bugs: overflow, off-by-one, uninitialized variables, wrong bounds, wrong reset between tests
- Debug prints and debug macros
- Assertions
- Stress testing (brute force vs. solution with random generators)
- Writing test generators
- Edge-case checklists (n = 1, all equal, max values, negative numbers, empty input)
- Sanitizers (`-fsanitize=address,undefined`) and `-D_GLIBCXX_DEBUG`
- Using a debugger (gdb)
- Time and memory profiling
- Diagnosing TLE: constant factors, unnecessary copies, I/O, recursion depth
- Diagnosing MLE and stack overflow

---

## 16. Contest Strategy

### 16.1 Individual Contests
- Reading all problems first vs. solving in order
- Problem difficulty estimation
- Time management
- When to skip a problem
- Penalty time and wrong-answer penalties
- Hacking strategy (Codeforces)
- Partial scoring and subtask strategy (IOI-style)

### 16.2 Team Contests (ICPC)
- Team roles and division of labor
- Using a single computer efficiently
- Printing code and debugging on paper
- Team reference document (TRD / team notebook)
- Communication and problem handoff
- Scoreboard reading and following the crowd

### 16.3 Mental and Physical Preparation
- Handling contest pressure and tilt
- Practicing under timed conditions
- Virtual contests
- Rest and consistency

---

## 17. Tools and Environment

- Online judges: Codeforces, AtCoder, LeetCode, CSES, SPOJ, UVa, Kattis, CodeChef, HackerRank, Timus, LibreOJ, Luogu, yukicoder, Library Checker, DMOJ, Project Euler
- Editors and IDEs (VS Code, Vim, CLion, Sublime, Emacs)
- Competitive programming helpers (Competitive Companion, CP Editor, cf-tool, oj / online-judge-tools)
- Building and maintaining a personal library (KACTL, AtCoder Library, ei1333 library, Nyaan's library)
- Snippets and templates
- Polygon (problem setting platform)
- Testlib (generators, validators, checkers)

---

## 18. Major Competitions

### 18.1 Olympiads (High School)
- IOI (International Olympiad in Informatics)
- Regional: APIO, BOI, CEOI, EGOI, JOI, NOI, EJOI
- National olympiads (USACO, CCO, COCI, etc.)
- Selection processes and training camps

### 18.2 University
- ICPC: regionals, continental finals, World Finals
- Regional university contests and training camps

### 18.3 Open Online Contests
- Codeforces rounds (Div. 1–4, Educational, Global)
- AtCoder (ABC, ARC, AGC, AHC)
- CodeChef (Starters, Cook-Off, Lunchtime legacy)
- TopCoder SRMs and TCO (legacy)
- Meta Hacker Cup
- Google Code Jam, Kick Start, Hash Code (discontinued, archives remain)
- LeetCode weekly and biweekly contests
- Advent of Code
- Universal Cup and other team-based online series

---

## 19. Problem Setting

- Writing clear problem statements
- Designing subtasks and constraints
- Generating strong test data
- Writing validators and checkers
- Anti-hash and anti-heuristic tests
- Writing reference and wrong solutions
- Testing and coordinating rounds
- Editorial writing

---

## 20. Practice, Learning, and Community

- Structured problem sets: CSES Problem Set, AtCoder Educational DP Contest, Codeforces EDU section, USACO Guide, NeetCode / Blind 75 (for interview-adjacent practice), Kattis problem archive
- Upsolving after contests
- Reading editorials effectively
- Rating-based practice ladders (A2OJ legacy, Codeforces problemset filters)
- Books: *Competitive Programmer's Handbook* (Laaksonen), *Competitive Programming 4* (Halim), *Introduction to Algorithms* (CLRS), *Algorithms* (Sedgewick), *Guide to Competitive Programming*, *Algorithm Design* (Kleinberg & Tardos)
- Websites: cp-algorithms.com, USACO Guide, OI Wiki, Codeforces blogs
- Video content and streams (Errichto, SecondThread, William Lin, Colin Galen)
- Training camps and clubs
- Communities: Codeforces blogs, Discord servers, university teams

---

## 21. Career and Beyond

- Competitive programming for technical interviews
- Transitioning skills to software engineering, research, and quantitative finance
- Algorithm engineering and performance-critical systems
- Coaching and mentoring
- Problem setting and contest organization as a hobby or career
- Balancing CP with other skills (system design, real-world engineering)

---

## Suggested Learning Path

1. **Foundations** → a language (C++), complexity, fast I/O (Section 1)
2. **STL and basic techniques** → sorting, binary search, prefix sums, two pointers, greedy (Sections 2–3)
3. **Basic data structures** → stacks, queues, heaps, hashing, DSU (Sections 4.1–4.4)
4. **Introductory DP** → classic problems, knapsack, LIS, grid DP (Sections 5.1–5.2)
5. **Graph basics** → BFS/DFS, shortest paths, MST, topological sort (Sections 6.1–6.5)
6. **Basic math** → number theory, modular arithmetic, combinatorics (Sections 9.1–9.2)
7. **Range queries** → prefix sums, sparse tables, Fenwick trees, segment trees (Sections 4.5–4.6)
8. **Trees** → LCA, Euler tour, tree DP (Sections 5.3, 7)
9. **Strings** → hashing, KMP, Z-function, tries (Sections 8.1–8.3)
10. **Intermediate DP and graphs** → bitmask DP, digit DP, SCC, bridges, 2-SAT, flows, matching (Sections 5.3, 6.5–6.9)
11. **Geometry basics** → cross product, convex hull, polygon area (Section 10)
12. **Advanced data structures and techniques** → lazy/persistent segment trees, HLD, centroid decomposition, Mo's algorithm, CDQ (Sections 4.6–4.8, 7.3, 11)
13. **Advanced math** → FFT/NTT, generating functions, Möbius inversion, game theory (Sections 9.5–9.6)
14. **DP optimizations and research-level topics** → CHT, D&C optimization, Aliens trick, slope trick (Section 5.4, 12)
15. **Throughout** → contest strategy, stress testing, upsolving, and regular contest participation (Sections 14–16, 20)
