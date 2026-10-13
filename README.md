# template-cp

[![tests](https://github.com/loselen/template-cp/actions/workflows/tests.yml/badge.svg)](https://github.com/loselen/template-cp/actions/workflows/tests.yml)

A competitive programming reference library in **Python 3 / PyPy**, plus a printable 25-page ICPC team reference document (TRD).

- 📄 **[`reference/trd.pdf`](reference/trd.pdf)**: the TRD. It is A4, has 25 pages, and puts the university name top-left and the page number top-right on every page.
- 🗺️ **[`reference/README.md`](reference/README.md)**: maps each topic to the file that covers it, plus a cheat sheet of theorems, Python performance tips, debugging, and contest strategy.
- 📋 **[`competitive-programming-topics.md`](competitive-programming-topics.md)**: the topic list this library covers.

## Layout

| File | Contents |
|---|---|
| `00_core.py` | Contest template, fast I/O, deep recursion, interactive problems, anti-hash, stress testing |
| `01_techniques.py` | Binary/ternary search, prefix sums, two pointers, greedy, bits, Mo's algorithm, CDQ, parallel binary search, A*/IDA*, annealing |
| `02_data_structures.py` | SortedList, DSU variants, Fenwick, segment trees (lazy, persistent, beats, Li Chao), treap, link-cut tree, wavelet matrix |
| `03_dp.py` | Classic DP, tree/bitmask/digit/SOS/profile DP, CHT, D&C, Knuth, Aliens, slope trick |
| `04_graphs.py` | Shortest paths, MST, SCC, 2-SAT, dominators, bridges/BCC, offline connectivity, Euler paths, cliques, coloring |
| `05_flow_matching.py` | Dinic, min-cost flow, lower bounds, Gomory–Hu, Stoer–Wagner, Hopcroft–Karp, Hungarian, blossom, matroid intersection |
| `06_trees.py` | LCA, HLD, centroid decomposition, small-to-large, virtual trees, tree isomorphism |
| `07_strings.py` | KMP, Z, rolling hash, Aho–Corasick, suffix array (SA-IS), suffix automaton, eertree, Manacher |
| `08_math_nt_comb.py` | Sieves, CRT, Miller–Rabin, Pollard rho, discrete log, prime counting, combinatorics, game theory, simplex |
| `09_algebra_poly.py` | Gaussian elimination, XOR basis, Berlekamp–Massey, NTT, power series, transforms, Bostan–Mori |
| `10_geometry.py` | Convex hull (2D/3D), half-plane intersection, circles, minimum enclosing circle, Manhattan MST, k-d tree |

Each technique is implemented once. Every snippet is self-contained and copy-pasteable, uses only the stdlib, and is iterative where depth could matter.

## Testing

Every file tests itself against brute force on random inputs:

```sh
./reference/run_tests.sh          # all files, about 40 s on CPython
python3 reference/04_graphs.py    # one file; prints "ALL OK"
```

CI runs the full suite on every push.

## Building the PDF

Requires [Quarto](https://quarto.org) and a LaTeX distribution (`quarto install tinytex` works).

```sh
cd reference
python3 build_trd.py --university "Your University" --team "Team Name"
```

The script generates `trd.qmd` from the `.py` files, leaving out the test blocks. It renders the PDF with xelatex, using 4 columns and syntax highlighting, and fails if the result exceeds 25 pages. Edit the `.py` files, not `trd.qmd`, which is overwritten on each build.

## License

[MIT](LICENSE)
