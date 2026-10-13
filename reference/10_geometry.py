# 10_geometry.py -- computational geometry (Section 10) + closest pair (3.7)
# Points are tuples (x, y). Integer inputs -> exact predicates (cross/dot/area2).
# Alternative: complex numbers (p.real, p.imag; cross = (a.conjugate()*b).imag) -- not used here.
# For exact intersection points pass fractions.Fraction coordinates (all '/' below then stay exact).
import math
import random
from bisect import bisect_left, bisect_right, insort
from collections import deque
from functools import cmp_to_key


# ==== 10.1 Basics ====
EPS = 1e-9
def sgn(x): return (x > EPS) - (x < -EPS)          # float sign with tolerance
def add(a, b): return (a[0] + b[0], a[1] + b[1])
def sub(a, b): return (a[0] - b[0], a[1] - b[1])
def mul(a, k): return (a[0] * k, a[1] * k)
def dot(a, b): return a[0] * b[0] + a[1] * b[1]
def cross(a, b): return a[0] * b[1] - a[1] * b[0]
def cross3(o, a, b): return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
def norm2(a): return a[0] * a[0] + a[1] * a[1]
def norm(a): return math.hypot(a[0], a[1])
def dist2(a, b): return norm2(sub(a, b))
def dist(a, b): return math.hypot(a[0] - b[0], a[1] - b[1])
def perp(a): return (-a[1], a[0])                   # rotate +90 deg (exact)
def rotate(a, t):                                   # rotate CCW by t radians
    c, s = math.cos(t), math.sin(t)
    return (a[0] * c - a[1] * s, a[0] * s + a[1] * c)
def orient(a, b, c):                                # +1 ccw (left turn), -1 cw, 0 collinear (exact for ints)
    v = cross3(a, b, c)
    return (v > 0) - (v < 0)

def dist_point_line(p, a, b):                       # distance p to infinite line ab (a != b)
    return abs(cross3(a, b, p)) / dist(a, b)

def dist_point_seg(p, a, b):                        # distance p to segment ab
    d = sub(b, a)
    if d == (0, 0): return dist(p, a)
    t = dot(sub(p, a), d)
    if t <= 0: return dist(p, a)
    if t >= norm2(d): return dist(p, b)
    return abs(cross(d, sub(p, a))) / norm(d)

# Polar sort by angle in [0, 2pi) around origin, exact (no atan2). Points != (0,0).
# Same direction: tie (here by distance). Sort around c: key on sub(p, c).
def _half(p): return 0 if (p[1] > 0 or (p[1] == 0 and p[0] > 0)) else 1
def _polar_cmp(a, b):
    ha, hb = _half(a), _half(b)
    if ha != hb: return ha - hb
    c = cross(a, b)
    if c: return -1 if c > 0 else 1
    return norm2(a) - norm2(b)
def polar_sort(pts): return sorted(pts, key=cmp_to_key(_polar_cmp))


# ==== 10.2 Lines and Segments ====
# Line forms: two points (a, b) [used here]; point+direction (p, d); ax+by=c with
# a=y2-y1, b=x1-x2, c=a*x1+b*y1 (integer-friendly; normalize by gcd & sign to hash lines);
# slope/intercept -- avoid (vertical lines, precision).

def line_inter(a, b, c, d):                         # lines ab, cd -> point or None if parallel
    den = cross(sub(b, a), sub(d, c))
    if den == 0: return None
    t = cross(sub(c, a), sub(d, c)) / den           # Fraction inputs -> exact
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)

def on_segment(p, a, b):                            # p on closed segment ab (exact for ints)
    return cross3(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) \
        and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])

def seg_intersect(a, b, c, d):                      # closed segments ab, cd share a point? exact; handles collinear/degenerate
    d1, d2 = orient(a, b, c), orient(a, b, d)
    d3, d4 = orient(c, d, a), orient(c, d, b)
    if d1 * d2 < 0 and d3 * d4 < 0: return True
    return (d1 == 0 and on_segment(c, a, b)) or (d2 == 0 and on_segment(d, a, b)) or \
           (d3 == 0 and on_segment(a, c, d)) or (d4 == 0 and on_segment(b, c, d))

def seg_inter_points(a, b, c, d):                   # [] none, [p] single point, [p, q] overlap segment (p<q)
    if not seg_intersect(a, b, c, d): return []
    if cross(sub(b, a), sub(d, c)) != 0: return [line_inter(a, b, c, d)]
    s = sorted(set(p for p in (a, b) if on_segment(p, c, d)) | set(p for p in (c, d) if on_segment(p, a, b)))
    return [s[0]] if len(s) == 1 else [s[0], s[-1]]

# Bentley-Ottmann (report all k intersections, O((n+k) log n)): sweep x, status = BST of
# segments ordered by y at sweep x, event queue = endpoints + found intersections; check only
# neighbours on insert/delete/swap. Painful in Python (needs ordered BST with changing key) --
# in contests usually n small -> O(n^2) pairwise seg_intersect.
# Shamos-Hoey ("does ANY pair intersect?", O(n log n)): same sweep without intersection events;
# on insert check new seg vs its above/below neighbours, on delete check the two neighbours
# that become adjacent; stop at first hit. Status order compares y at current x (exact via cross).


# ==== 10.3 Polygons ====
def area2(P):                                       # 2 * signed area (ccw > 0), exact integer
    s = 0
    for i in range(len(P)):
        x1, y1 = P[i - 1]; x2, y2 = P[i]
        s += x1 * y2 - x2 * y1
    return s

def boundary_points(P):                             # lattice points on boundary of integer polygon
    return sum(math.gcd(abs(P[i][0] - P[i - 1][0]), abs(P[i][1] - P[i - 1][1])) for i in range(len(P)))

def pick_interior(P):                               # Pick: A = I + B/2 - 1  ->  I = (2A - B + 2) / 2
    return (abs(area2(P)) - boundary_points(P) + 2) // 2

def point_in_polygon(P, q):                         # any simple polygon: 1 inside, 0 boundary, -1 outside. O(n), exact
    wn = 0
    n = len(P)
    for i in range(n):
        a, b = P[i], P[(i + 1) % n]
        if on_segment(q, a, b): return 0
        if a[1] <= q[1]:
            if b[1] > q[1] and cross3(a, b, q) > 0: wn += 1
        elif b[1] <= q[1] and cross3(a, b, q) < 0: wn -= 1
    return 1 if wn else -1                          # winding number (also works for self-intersecting, nonzero rule)

def point_in_convex(P, q):                          # P ccw strictly convex (n>=3): 1/0/-1, O(log n)
    n, a = len(P), P[0]
    if cross3(a, P[1], q) < 0 or cross3(a, P[-1], q) > 0: return -1
    lo, hi = 1, n - 1
    while hi - lo > 1:
        m = (lo + hi) // 2
        if cross3(a, P[m], q) >= 0: lo = m
        else: hi = m
    c = cross3(P[lo], P[lo + 1], q)
    if c < 0: return -1
    if c == 0: return 0
    if (lo == 1 and cross3(a, P[1], q) == 0) or (lo == n - 2 and cross3(a, P[-1], q) == 0): return 0
    return 1

def convex_cut(P, a, b):                            # part of convex polygon P left of (on) directed line a->b. O(n)
    res, n = [], len(P)
    for i in range(n):
        p, q = P[i], P[(i + 1) % n]
        cp, cq = cross3(a, b, p), cross3(a, b, q)
        if cp >= 0: res.append(p)
        if cp * cq < 0: res.append(line_inter(p, q, a, b))
    return res                                      # may contain duplicate/collinear pts in degenerate cases

# Sutherland-Hodgman: clip polygon S (any) by CONVEX clipper C = apply convex_cut(S, C[i], C[i+1])
# for each ccw edge of C. O(|S||C|). Result for concave S may have zero-width "bridges".
# Polygon triangulation (conceptual): ear clipping O(n^2): repeatedly find convex vertex i whose
# triangle (i-1,i,i+1) contains no other vertex, cut it. Fan triangulation for convex polygons.
# O(n log n): split into monotone pieces (sweep) then triangulate each in O(n).


# ==== 10.4 Convex Hull ====
def convex_hull(pts, keep_collinear=False):         # Andrew monotone chain O(n log n); ccw, starts at min (x,y)
    pts = sorted(set(pts))                          # dedupe!
    if len(pts) <= 2: return pts
    def chain(seq):
        h = []
        for p in seq:
            while len(h) >= 2:
                c = cross3(h[-2], h[-1], p)
                if c < 0 or (c == 0 and not keep_collinear): h.pop()
                else: break
            h.append(p)
        return h
    lo, up = chain(pts), chain(reversed(pts))
    if keep_collinear and len(lo) == len(pts): return lo   # all collinear
    return lo[:-1] + up[:-1]                        # strict mode: all collinear -> 2 endpoints

# Graham scan: pick lowest point, polar_sort others around it (ties by distance), stack-pop on
# non-left turns. Same O(n log n); monotone chain is simpler (no angle sort).

def hull_diameter2(H):                              # rotating calipers: max squared dist, H strict ccw hull. O(n)
    n = len(H)
    if n == 1: return 0
    if n == 2: return dist2(H[0], H[1])
    best, j = 0, 1
    for i in range(n):
        ni = (i + 1) % n
        while cross(sub(H[ni], H[i]), sub(H[(j + 1) % n], H[j])) > 0: j = (j + 1) % n
        best = max(best, dist2(H[i], H[j]), dist2(H[ni], H[j]))
    return best                                     # farthest pair: record (i,j) at the max

def hull_width(H):                                  # min width (min over edges of max dist), H strict ccw, n>=3. O(n)
    n, j, best = len(H), 1, float('inf')
    for i in range(n):
        a, b = H[i], H[(i + 1) % n]
        while cross3(a, b, H[(j + 1) % n]) > cross3(a, b, H[j]): j = (j + 1) % n
        best = min(best, cross3(a, b, H[j]) / dist(a, b))
    return best

def minkowski(P, Q):                                # convex ccw P, Q -> ccw P+Q, O(n+m) (may keep collinear pts)
    def reorder(R):
        k = min(range(len(R)), key=lambda t: (R[t][1], R[t][0]))
        R = R[k:] + R[:k]
        return R + R[:2]
    P, Q = reorder(P), reorder(Q)
    res, i, j = [], 0, 0
    while i < len(P) - 2 or j < len(Q) - 2:
        res.append(add(P[i], Q[j]))
        c = cross(sub(P[i + 1], P[i]), sub(Q[j + 1], Q[j]))
        if c >= 0 and i < len(P) - 2: i += 1
        if c <= 0 and j < len(Q) - 2: j += 1
    return res or [add(P[0], Q[0])]
# Uses: collision of convex P and Q <=> origin in P + (-Q); distance between polygons = dist(0, P-Q).

# Dynamic convex hull (online point insertion): keep upper & lower hulls separately in a
# sorted container keyed by x (SortedList / treap). Insert p: if under upper hull (check
# neighbours) skip; else insert and pop neighbours on both sides while they make non-right turns.
# Amortized O(log n). For max-of-lines queries (CHT) use Li Chao tree (02_data_structures.py).

def half_plane_intersection(H):                     # H: list of (p, q) keep LEFT of p->q. O(n log n)
    # Include a big bounding box if the region may be unbounded. Returns ccw vertices (floats), [] if empty.
    # Integer p, q -> exact angle sort/parallel tests; Fraction coords -> fully exact (set EPS=0 semantics).
    def dcmp(h1, h2):
        a, b = sub(h1[1], h1[0]), sub(h2[1], h2[0])
        ha, hb = _half(a), _half(b)
        if ha != hb: return ha - hb
        c = cross(a, b)
        return -1 if c > 0 else (1 if c < 0 else 0)
    H = sorted(H, key=cmp_to_key(dcmp))
    def out(h, r): return cross(sub(h[1], h[0]), sub(r, h[0])) < -EPS
    def inter(h1, h2): return line_inter(h1[0], h1[1], h2[0], h2[1])
    dq = deque()
    for h in H:
        while len(dq) > 1 and out(h, inter(dq[-1], dq[-2])): dq.pop()
        while len(dq) > 1 and out(h, inter(dq[0], dq[1])): dq.popleft()
        if dq:
            d1, d2 = sub(h[1], h[0]), sub(dq[-1][1], dq[-1][0])
            if cross(d1, d2) == 0:
                if dot(d1, d2) < 0: return []
                if out(h, dq[-1][0]): dq.pop()
                else: continue
        dq.append(h)
    while len(dq) > 2 and out(dq[0], inter(dq[-1], dq[-2])): dq.pop()
    while len(dq) > 2 and out(dq[-1], inter(dq[0], dq[1])): dq.popleft()
    if len(dq) < 3: return []
    res = []
    for i in range(len(dq)):
        p = inter(dq[i], dq[(i + 1) % len(dq)])
        if p is None: return []
        res.append(p)
    return res

def convex_hull_3d(P):                              # P: list of int (x,y,z). Incremental O(n^2), exact.
    # Returns faces as index triples (i,j,k), ccw seen from outside (outward normal). Coplanar
    # hull faces come out triangulated. Returns [] if all points coplanar (handle as 2D).
    def sub3(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    def vol(a, b, c, d):                            # 6 * signed volume; > 0: d above plane abc (ccw)
        u, v, w = sub3(P[b], P[a]), sub3(P[c], P[a]), sub3(P[d], P[a])
        return (u[1] * v[2] - u[2] * v[1]) * w[0] + (u[2] * v[0] - u[0] * v[2]) * w[1] + (u[0] * v[1] - u[1] * v[0]) * w[2]
    n = len(P)
    i0 = 0
    i1 = next((i for i in range(n) if P[i] != P[i0]), None)
    if i1 is None: return []
    def colin(i):
        u, v = sub3(P[i1], P[i0]), sub3(P[i], P[i0])
        return u[1] * v[2] - u[2] * v[1] == 0 and u[2] * v[0] - u[0] * v[2] == 0 and u[0] * v[1] - u[1] * v[0] == 0
    i2 = next((i for i in range(n) if not colin(i)), None)
    if i2 is None: return []
    i3 = next((i for i in range(n) if vol(i0, i1, i2, i) != 0), None)
    if i3 is None: return []
    if vol(i0, i1, i2, i3) > 0: i1, i2 = i2, i1     # make i3 below face (i0,i1,i2)
    faces = [(i0, i1, i2), (i0, i3, i1), (i1, i3, i2), (i2, i3, i0)]
    used = {i0, i1, i2, i3}
    for p in range(n):
        if p in used: continue
        vis = [vol(a, b, c, p) > 0 for (a, b, c) in faces]
        if not any(vis): continue                   # inside or on hull
        edges = set()
        for f, v in zip(faces, vis):
            if v: edges.update(((f[0], f[1]), (f[1], f[2]), (f[2], f[0])))
        nf = [f for f, v in zip(faces, vis) if not v]
        for (a, b) in edges:
            if (b, a) not in edges: nf.append((a, b, p))   # horizon edge
        faces = nf
    return faces
# Volume = sum over faces of dot(Pa, cross(Pb, Pc)) / 6 ; surface area = sum |cross(Pb-Pa, Pc-Pa)|/2.


# ==== 10.5 Circles ====  (circle = center tuple c, radius r; floats)
def circle_line(c, r, a, b):                        # intersections of circle with line ab: 0/1/2 points
    d = sub(b, a)
    t = dot(sub(c, a), d) / norm2(d)
    h = (a[0] + d[0] * t, a[1] + d[1] * t)          # foot of perpendicular
    s = r * r - dist2(c, h)
    if s < -EPS: return []
    if s <= EPS: return [h]
    k = math.sqrt(s / norm2(d))
    return [(h[0] - d[0] * k, h[1] - d[1] * k), (h[0] + d[0] * k, h[1] + d[1] * k)]

def circle_circle(c1, r1, c2, r2):                  # 0/1/2 points; identical circles -> [] (infinite)
    d2 = dist2(c1, c2)
    if d2 == 0: return []
    d = math.sqrt(d2)
    if d > r1 + r2 + EPS or d < abs(r1 - r2) - EPS: return []
    a = (d2 + r1 * r1 - r2 * r2) / (2 * d)
    h = math.sqrt(max(0.0, r1 * r1 - a * a))
    v = ((c2[0] - c1[0]) / d, (c2[1] - c1[1]) / d)
    m = (c1[0] + v[0] * a, c1[1] + v[1] * a)
    if h <= EPS: return [m]
    return [(m[0] - v[1] * h, m[1] + v[0] * h), (m[0] + v[1] * h, m[1] - v[0] * h)]

def tangents_point(c, r, p):                        # tangent points on circle from external point p
    d2 = dist2(c, p)
    if d2 < r * r - EPS: return []
    if d2 <= r * r + EPS: return [p]
    return circle_circle(c, r, p, math.sqrt(d2 - r * r))

def tangents_circles(c1, r1, c2, r2):               # external tangents; call with -r2 for internal ones
    # returns list of (p1, p2): tangent line touches circle1 at p1, circle2 at p2. r2=0 -> point tangents.
    d = sub(c2, c1)
    dr, d2 = r1 - r2, norm2(d)
    h2 = d2 - dr * dr
    if d2 == 0 or h2 < 0: return []
    out = []
    for s in (-1, 1):
        k = math.sqrt(h2) * s
        v = ((d[0] * dr - d[1] * k) / d2, (d[1] * dr + d[0] * k) / d2)
        out.append(((c1[0] + v[0] * r1, c1[1] + v[1] * r1), (c2[0] + v[0] * r2, c2[1] + v[1] * r2)))
    if h2 == 0: out.pop()
    return out

def circumcircle(a, b, c):                          # center, radius of circle through 3 pts (None if collinear)
    bx, by = b[0] - a[0], b[1] - a[1]
    cx, cy = c[0] - a[0], c[1] - a[1]
    d = 2 * (bx * cy - by * cx)
    if d == 0: return None
    B, C = bx * bx + by * by, cx * cx + cy * cy
    o = (a[0] + (cy * B - by * C) / d, a[1] + (bx * C - cx * B) / d)
    return o, dist(o, a)

def min_enclosing_circle(pts):                      # Welzl (iterative, shuffled) expected O(n) -> (center, r)
    P = [(float(x), float(y)) for x, y in pts]
    random.shuffle(P)
    c, r = P[0], 0.0
    def outside(p): return dist(p, c) > r + 1e-9 * max(1.0, r)
    for i in range(len(P)):
        if not outside(P[i]): continue
        c, r = P[i], 0.0
        for j in range(i):
            if not outside(P[j]): continue
            c = ((P[i][0] + P[j][0]) / 2, (P[i][1] + P[j][1]) / 2); r = dist(c, P[i])
            for k in range(j):
                if not outside(P[k]): continue
                cc = circumcircle(P[i], P[j], P[k])
                if cc: c, r = cc
    return c, r

def circle_union_area(cs):                          # cs: list of (cx, cy, r). Green's theorem O(n^2 log n)
    cs = list(set(cs))
    keep = []
    for i, (x, y, r) in enumerate(cs):              # drop circles inside another
        if r <= 0: continue
        if not any(j != i and rj >= r and math.hypot(x - xj, y - yj) + r <= rj + EPS
                   for j, (xj, yj, rj) in enumerate(cs)): keep.append((x, y, r))
    total = 0.0
    for i, (x, y, r) in enumerate(keep):
        iv = []                                     # covered angle intervals
        for j, (xj, yj, rj) in enumerate(keep):
            if i == j: continue
            d = math.hypot(xj - x, yj - y)
            if d >= r + rj: continue
            a = math.atan2(yj - y, xj - x)
            h = math.acos(max(-1.0, min(1.0, (r * r + d * d - rj * rj) / (2 * r * d))))
            lo, hi = a - h, a + h
            if lo < -math.pi: iv += [(lo + 2 * math.pi, math.pi), (-math.pi, hi)]
            elif hi > math.pi: iv += [(lo, math.pi), (-math.pi, hi - 2 * math.pi)]
            else: iv.append((lo, hi))
        iv.sort()
        F = lambda t: r * r * t + x * r * math.sin(t) - y * r * math.cos(t)
        cur = -math.pi
        for lo, hi in iv + [(math.pi, math.pi)]:
            if lo > cur: total += F(lo) - F(cur)    # uncovered arc [cur, lo]
            cur = max(cur, hi)
    return total / 2
# Alternative: Simpson integration over x of union length of y-intervals, splitting at all
# critical x (cx +- r, pairwise intersection x's) -- see brute in tests.


# ==== 10.6 Advanced Geometry ====
# ---- closest pair (3.7 / 10.6) ----
def closest_pair(pts):                              # int pts, n>=2 -> (d2, p, q). Sweep by x + y-sorted strip. ~O(n log n)
    P = sorted(pts)
    best, bp = float('inf'), None
    strip, left = [], 0                             # strip: sorted (y, x) of pts with x within sqrt(best)
    for p in P:
        while left < len(P) and (p[0] - P[left][0]) ** 2 > best:
            q = P[left]; strip.pop(bisect_left(strip, (q[1], q[0]))); left += 1
        d = math.isqrt(best) if best != float('inf') else 1 << 62   # floats: use math.sqrt
        lo = bisect_left(strip, (p[1] - d, -float('inf')))
        hi = bisect_right(strip, (p[1] + d, float('inf')))
        for k in range(lo, hi):
            y, x = strip[k]
            dd = (x - p[0]) ** 2 + (y - p[1]) ** 2
            if dd < best:
                best, bp = dd, ((x, y), p)
                if best == 0: return 0, bp[0], bp[1]
        insort(strip, (p[1], p[0]))
    return best, bp[0], bp[1]
# D&C alternative: split by median x, recurse, merge by y, check strip neighbours (<=7 each).

# Sweep line techniques (notes): sort events by x (ties: insert before query before delete, or
# as problem demands) and maintain active set in a balanced structure / segment tree:
#   union area of rectangles: sweep x, segment tree over compressed y with cover count + covered length.
#   segment intersection: see Bentley-Ottmann / Shamos-Hoey above. closest pair: above.
#   max overlapping intervals / points covered by most circles: angular or linear events.

# Voronoi / Delaunay (notes): Delaunay = dual of Voronoi; no point inside any triangle's circumcircle.
#   naive O(n^4): triangle ijk is Delaunay iff no l inside circumcircle (exact in_circle det below).
#   O(n^2): 3D lift (x, y, x^2+y^2) -> lower hull of convex_hull_3d = Delaunay triangles
#     (faces with outward normal z < 0). Bowyer-Watson O(n^2) worst: insert point, delete triangles
#   whose circumcircle contains it, re-triangulate cavity. Fortune sweep O(n log n) - rarely coded.
#   EMST is a subset of Delaunay edges. Voronoi cell of i = half_plane_intersection of bisectors.
#   in_circle(a,b,c,d) (a,b,c ccw) > 0 iff d inside: det of rows (ax-dx, ay-dy, (ax-dx)^2+(ay-dy)^2), ...
# Point location (notes): offline -> sweep x with sorted-by-y active edges; online -> slab
#   decomposition (O(n^2) memory) or persistent BST; convex polygon -> point_in_convex.

# ---- Manhattan MST ----
def manhattan_mst(pts):                             # O(n log n): 4 octant sweeps + Fenwick(min) + Kruskal -> (weight, edges)
    n = len(pts)
    cand = []
    for tf in (lambda x, y: (x, y), lambda x, y: (y, x), lambda x, y: (-x, y), lambda x, y: (y, -x)):
        Q = [tf(x, y) for x, y in pts]
        # for each i find j with xj >= xi and yj - xj >= yi - xi minimizing xj + yj
        keys = sorted(set(y - x for x, y in Q))
        m = len(keys)
        bit = [(float('inf'), -1)] * (m + 1)        # prefix-min over reversed key index
        for i in sorted(range(n), key=lambda t: (-Q[t][0], -Q[t][1])):
            x, y = Q[i]
            k = m - bisect_left(keys, y - x)        # reversed 1-based: keys >= y-x  <->  idx <= k
            j, best = k, (float('inf'), -1)
            while j > 0:
                if bit[j] < best: best = bit[j]
                j -= j & -j
            if best[1] >= 0: cand.append((best[0] - x - y, i, best[1]))
            j, v = k, (x + y, i)
            while j <= m:
                if v < bit[j]: bit[j] = v
                j += j & -j
    cand.sort()
    par = list(range(n))
    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    total, edges = 0, []
    for w, a, b in cand:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb; total += w; edges.append((a, b, w))
    return total, edges

# ---- k-d tree (2D nearest neighbour) ----
class KDTree:                                       # build O(n log^2 n); query ~O(log n) avg (bad worst case)
    def __init__(self, pts):
        self.a = list(pts)
        st = [(0, len(self.a), 0)]
        while st:
            lo, hi, d = st.pop()
            if hi - lo <= 1: continue
            self.a[lo:hi] = sorted(self.a[lo:hi], key=lambda p: p[d])
            m = (lo + hi) // 2
            st.append((lo, m, d ^ 1)); st.append((m + 1, hi, d ^ 1))
    def nearest(self, q, skip_equal=False):         # -> (d2, point); skip_equal: ignore points == q
        best, bp = float('inf'), None
        st = [(0, len(self.a), 0, 0)]               # (lo, hi, axis, lower bound d2)
        while st:
            lo, hi, d, lb = st.pop()
            if lo >= hi or lb >= best: continue
            m = (lo + hi) // 2
            p = self.a[m]
            dd = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2
            if dd < best and not (skip_equal and p == q): best, bp = dd, p
            diff = q[d] - p[d]
            near, far = ((lo, m), (m + 1, hi)) if diff < 0 else ((m + 1, hi), (lo, m))
            st.append((far[0], far[1], d ^ 1, diff * diff))
            st.append((near[0], near[1], d ^ 1, 0))
        return best, bp

# Lattice point counting (notes): points on segment = gcd(|dx|,|dy|)+1; inside polygon = Pick
#   (pick_interior). Under a line y = (a*x+b)/c for x in [0,n): floor_sum (08_math_nt_comb.py).
#   Points in circle x^2+y^2<=r^2: sum over x of 2*isqrt(r^2-x^2)+1.
# Rotational sweep (notes): around a pivot, polar_sort other points (or events = angles where
#   something enters/leaves), sweep angle maintaining count; e.g. max points covered by a
#   half-plane/line through pivot, or for all pairs sort O(n^2) directions and keep points sorted
#   by projection (swap adjacent pair at each direction) -> O(n^2 log n) for min-area triangle etc.
# Geometric duality (notes): point (a, b) <-> line y = a*x - b. Preserves incidence & above/below
#   (point above line <=> dual line below dual point). Upper hull of points <-> lower envelope of
#   lines (CHT). Problems about lines through points become arrangements of lines.


# ==== tests ====
if __name__ == "__main__":
    from fractions import Fraction as Fr
    from itertools import combinations
    random.seed(12345)
    R = random.randint
    def rp(c=6): return (R(-c, c), R(-c, c))

    # basics
    assert cross((1, 0), (0, 1)) == 1 and orient((0, 0), (1, 0), (2, 0)) == 0
    assert abs(dist_point_line((0, 5), (-1, 0), (1, 0)) - 5) < 1e-9
    assert abs(dist_point_seg((3, 4), (-1, 0), (0, 0)) - 5) < 1e-9 and abs(dist_point_seg((0, 2), (-1, 0), (1, 0)) - 2) < 1e-9
    assert all(abs(a - b) < 1e-9 for a, b in zip(rotate((1, 0), math.pi / 2), (0, 1)))
    for _ in range(200):                            # dist_point_seg vs sampling
        p, a, b = rp(), rp(), rp()
        bf = min(dist(p, (a[0] + (b[0] - a[0]) * t / 2000, a[1] + (b[1] - a[1]) * t / 2000)) for t in range(2001))
        assert abs(dist_point_seg(p, a, b) - bf) < 1e-2
    for _ in range(100):                            # polar sort vs atan2
        pts = [p for p in (rp() for _ in range(30)) if p != (0, 0)]
        ang = [math.atan2(y, x) % (2 * math.pi) for x, y in polar_sort(pts)]
        assert all(ang[i] <= ang[i + 1] + 1e-12 for i in range(len(ang) - 1))

    # segments: vs Fraction-parameter brute
    def seg_brute(a, b, c, d):
        r, s = sub(b, a), sub(d, c)
        den = cross(r, s)
        if den != 0:
            t, u = Fr(cross(sub(c, a), s), den), Fr(cross(sub(c, a), r), den)
            return [(a[0] + r[0] * t, a[1] + r[1] * t)] if 0 <= t <= 1 and 0 <= u <= 1 else []
        if cross(r, sub(c, a)) != 0: return []
        cand = [p for p in (a, b) if on_segment(p, c, d)] + [p for p in (c, d) if on_segment(p, a, b)]
        cand = sorted(set(cand))
        return [] if not cand else ([cand[0]] if len(cand) == 1 else [cand[0], cand[-1]])
    for _ in range(5000):
        a, b, c, d = rp(3), rp(3), rp(3), rp(3)
        got = seg_inter_points(a, b, c, d)
        if len(got) == 1 and cross(sub(b, a), sub(d, c)) != 0:
            got = seg_inter_points(*[(Fr(x), Fr(y)) for x, y in (a, b, c, d)])
        assert got == seg_brute(a, b, c, d), (a, b, c, d, got)
        assert seg_intersect(a, b, c, d) == bool(got)

    # polygons: star-shaped random polygons around origin
    def star_poly(c=8, k=10):
        while True:
            pts = list({p for p in (rp(c) for _ in range(k)) if p != (0, 0)})
            P, seen = [], set()
            for p in polar_sort(pts):
                g = math.gcd(abs(p[0]), abs(p[1])); dr = (p[0] // g, p[1] // g)
                if dr not in seen: seen.add(dr); P.append(p)
            if len(P) >= 3 and all(cross(P[i - 1], P[i]) > 0 for i in range(len(P))): return P
    def tri_state(a, b, c, q):                      # ccw triangle: 1 in, 0 on, -1 out
        s = [orient(a, b, q), orient(b, c, q), orient(c, a, q)]
        return -1 if min(s) < 0 else (1 if min(s) > 0 else 0)
    for _ in range(150):
        P = star_poly()
        assert area2(P) == sum(cross(P[i - 1], P[i]) for i in range(len(P))) > 0
        cnt = {1: 0, 0: 0, -1: 0}
        for x in range(-9, 10):
            for y in range(-9, 10):
                q = (x, y)
                on = any(on_segment(q, P[i - 1], P[i]) for i in range(len(P)))
                ins = any(tri_state((0, 0), P[i - 1], P[i], q) >= 0 for i in range(len(P)))
                exp = 0 if on else (1 if ins else -1)
                got = point_in_polygon(P, q)
                assert got == exp
                assert point_in_polygon(P[::-1], q) == exp   # cw too
                cnt[got] += 1
        assert cnt[0] == boundary_points(P) and cnt[1] == pick_interior(P)

    # hull + point in convex + calipers + width
    for it in range(300):
        pts = [rp(5) for _ in range(R(1, 12))]
        H = convex_hull(pts)
        S = set(pts)
        brute = set()                               # strict vertices: not in any triangle / segment of others
        for p in S:
            o = [q for q in S if q != p]
            if any(on_segment(p, a, b) for a, b in combinations(o, 2)): continue
            if any(tri_state(*((a, b, c) if orient(a, b, c) > 0 else (a, c, b)), p) >= 0
                   for a, b, c in combinations(o, 3) if orient(a, b, c) != 0): continue
            brute.add(p)
        assert set(H) == brute and len(H) == len(brute)
        if len(H) >= 3:
            assert all(cross3(H[i - 2], H[i - 1], H[i]) > 0 for i in range(len(H)))
            Hc = convex_hull(pts, keep_collinear=True)
            assert set(Hc) == {p for p in S if point_in_polygon(H, p) == 0} and len(Hc) == len(set(Hc))
            assert all(cross3(Hc[i - 2], Hc[i - 1], Hc[i]) >= 0 for i in range(len(Hc))) and area2(Hc) == area2(H)
            for _ in range(30):
                q = rp(6)
                assert point_in_convex(H, q) == point_in_polygon(H, q)
            w = min(max(cross3(H[i - 1], H[i], p) for p in H) / dist(H[i - 1], H[i]) for i in range(len(H)))
            assert abs(hull_width(H) - w) < 1e-9
        elif len(S) >= 2 and len({orient(min(S), max(S), p) for p in S}) == 1:
            assert sorted(convex_hull(pts, True)) == sorted(S)   # all collinear
        assert hull_diameter2(H) == max([dist2(a, b) for a in S for b in S])

    # minkowski
    for _ in range(300):
        P = convex_hull([rp(5) for _ in range(R(1, 8))])
        Q = convex_hull([rp(5) for _ in range(R(1, 8))])
        M = minkowski(P, Q)
        assert convex_hull(M) == convex_hull([add(p, q) for p in P for q in Q])

    # half-plane intersection vs clipping a box
    for _ in range(300):
        B = 20
        box = [(-B, -B), (B, -B), (B, B), (-B, B)]
        Hs = [(box[i], box[(i + 1) % 4]) for i in range(4)]
        for _ in range(R(1, 7)):
            a, b = rp(10), rp(10)
            if a != b: Hs.append((a, b))
        poly = [(Fr(x), Fr(y)) for x, y in box]
        for a, b in Hs:
            poly = convex_cut(poly, a, b)
        ex = abs(area2(poly)) / 2 if len(poly) >= 3 else 0
        res = half_plane_intersection(Hs)
        got = abs(area2(res)) / 2 if res else 0
        assert abs(float(ex) - got) < 1e-6, (ex, got)
        if res: assert area2(res) >= -1e-9          # ccw

    # 3D hull
    def check3d(P):
        F = convex_hull_3d(P)
        def s3(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
        def cr(u, v): return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        def d3(u, v): return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
        E = set()
        for a, b, c in F:
            nrm = cr(s3(P[b], P[a]), s3(P[c], P[a]))
            assert nrm != (0, 0, 0)
            assert all(d3(nrm, s3(q, P[a])) <= 0 for q in P)   # supporting, outward
            for e in ((a, b), (b, c), (c, a)):
                assert e not in E; E.add(e)
        assert all((b, a) in E for a, b in E)                    # closed surface
        V = {i for f in F for i in f}
        assert len(V) - len(E) // 2 + len(F) == 2                # Euler
        return F
    for _ in range(100):                            # degenerate small coords
        P = [(R(0, 3), R(0, 3), R(0, 3)) for _ in range(R(4, 25))]
        if len({(x, y, z) for x, y, z in P}) >= 4:
            F = convex_hull_3d(P)
            if F: check3d(P)
    for _ in range(40):                             # general position: compare with brute faces
        P = list({(R(-1000, 1000), R(-1000, 1000), R(-1000, 1000)) for _ in range(R(4, 12))})
        F = check3d(P)
        def v6(a, b, c, d):
            u = [P[b][t] - P[a][t] for t in range(3)]; v = [P[c][t] - P[a][t] for t in range(3)]; w = [P[d][t] - P[a][t] for t in range(3)]
            return (u[1] * v[2] - u[2] * v[1]) * w[0] + (u[2] * v[0] - u[0] * v[2]) * w[1] + (u[0] * v[1] - u[1] * v[0]) * w[2]
        if any(v6(*q) == 0 for q in combinations(range(len(P)), 4)): continue
        bf = set()
        for t in combinations(range(len(P)), 3):
            s = {v6(*t, d) > 0 for d in range(len(P)) if d not in t}
            if len(s) == 1: bf.add(t)
        assert {tuple(sorted(f)) for f in F} == bf

    # circles
    for _ in range(500):
        c, r = rp(), R(1, 6)
        a, b = rp(), rp()
        if a == b: continue
        I = circle_line(c, r, a, b)
        dl = dist_point_line(c, a, b)
        assert len(I) == (2 if dl < r - 1e-9 else (1 if dl < r + 1e-9 else 0))
        for p in I: assert abs(dist(p, c) - r) < 1e-6 and abs(cross3(a, b, p)) < 1e-6
        c2, r2 = rp(), R(1, 6)
        if c2 == c: continue
        I = circle_circle(c, r, c2, r2)
        d = dist(c, c2)
        if abs(d - (r + r2)) > 1e-9 and abs(d - abs(r - r2)) > 1e-9:
            assert len(I) == (2 if abs(r - r2) < d < r + r2 else 0)
        for p in I: assert abs(dist(p, c) - r) < 1e-6 and abs(dist(p, c2) - r2) < 1e-6
        for p in tangents_point(c, r, a):
            assert abs(dist(p, c) - r) < 1e-6 and (p == a or abs(dot(sub(p, c), sub(p, a))) < 1e-6)
        ext, inn = tangents_circles(c, r, c2, r2), tangents_circles(c, r, c2, -r2)
        if d > r + r2 + 1e-9: assert len(ext) == 2 and len(inn) == 2
        if abs(r - r2) + 1e-9 < d < r + r2 - 1e-9: assert len(ext) == 2 and len(inn) == 0
        for p1, p2 in ext + inn:
            assert abs(dist(p1, c) - r) < 1e-6 and abs(dist(p2, c2) - r2) < 1e-6
            if dist(p1, p2) > 1e-6:
                assert abs(dist_point_line(c, p1, p2) - r) < 1e-6 and abs(dist_point_line(c2, p1, p2) - r2) < 1e-6
        for p1, p2 in ext:
            if dist(p1, p2) > 1e-6: assert orient(p1, p2, c) * orient(p1, p2, c2) >= 0
    # MEC vs brute over pairs/triples
    for _ in range(200):
        pts = [rp(10) for _ in range(R(1, 9))]
        c, r = min_enclosing_circle(pts)
        assert all(dist(p, c) <= r + 1e-7 for p in pts)
        cands = [(p, 0.0) for p in pts]
        for a, b in combinations(pts, 2):
            m = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2); cands.append((m, dist(m, a)))
        for a, b, cc in combinations(pts, 3):
            t = circumcircle(a, b, cc)
            if t: cands.append(t)
        br = min(rr for cc, rr in cands if all(dist(p, cc) <= rr + 1e-7 for p in pts))
        assert abs(br - r) < 1e-7
    # union of circles vs Simpson
    def union_simpson(cs):
        xs = set()
        for x, y, r in cs: xs.update((x - r, x + r))
        for (x1, y1, r1), (x2, y2, r2) in combinations(cs, 2):
            for p in circle_circle((x1, y1), r1, (x2, y2), r2): xs.add(p[0])
        xs = sorted(xs)
        def f(x):
            iv = []
            for cx, cy, r in cs:
                if abs(x - cx) < r:
                    h = math.sqrt(r * r - (x - cx) ** 2); iv.append((cy - h, cy + h))
            iv.sort(); tot, cur = 0.0, -1e18
            for lo, hi in iv:
                lo = max(lo, cur)
                if hi > lo: tot += hi - lo
                cur = max(cur, hi)
            return tot
        s = 0.0
        for a, b in zip(xs, xs[1:]):
            N = 400; hh = (b - a) / N
            s += hh / 3 * (f(a) + f(b) + sum((4 if k % 2 else 2) * f(a + k * hh) for k in range(1, N)))
        return s
    for _ in range(30):
        cs = [(R(-5, 5), R(-5, 5), R(1, 4)) for _ in range(R(1, 5))]
        if random.random() < 0.3: cs.append(cs[0])
        A, Bv = circle_union_area(cs), union_simpson(list(set(cs)))
        assert abs(A - Bv) < 1e-2 * max(1, Bv), (cs, A, Bv)
    assert abs(circle_union_area([(0, 0, 1)]) - math.pi) < 1e-9

    # closest pair
    for _ in range(300):
        pts = [rp(R(1, 50)) for _ in range(R(2, 40))]
        d2, p, q = closest_pair(pts)
        bf = min(dist2(pts[i], pts[j]) for i in range(len(pts)) for j in range(i))
        assert d2 == bf == dist2(p, q)
    big = [(R(0, 10 ** 9), R(0, 10 ** 9)) for _ in range(20000)]
    closest_pair(big)

    # Manhattan MST vs Prim
    for _ in range(200):
        pts = [rp(R(1, 20)) for _ in range(R(1, 30))]
        n = len(pts)
        w, e = manhattan_mst(pts)
        D = lambda i, j: abs(pts[i][0] - pts[j][0]) + abs(pts[i][1] - pts[j][1])
        inT, dd, tot = [False] * n, [float('inf')] * n, 0
        dd[0] = 0
        for _ in range(n):
            u = min((i for i in range(n) if not inT[i]), key=lambda i: dd[i])
            inT[u] = True; tot += dd[u]
            for v in range(n):
                if not inT[v]: dd[v] = min(dd[v], D(u, v))
        assert w == tot and len(e) == n - 1

    # k-d tree
    for _ in range(100):
        pts = [rp(R(1, 100)) for _ in range(R(1, 60))]
        T = KDTree(pts)
        for _ in range(20):
            q = rp(110)
            assert T.nearest(q)[0] == min(dist2(p, q) for p in pts)
            others = [p for p in pts if p != q]
            if others: assert T.nearest(q, True)[0] == min(dist2(p, q) for p in others)
    print("ALL OK")
