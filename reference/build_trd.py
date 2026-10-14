#!/usr/bin/env python3
"""Build the ICPC team reference document: trd.qmd -> trd.pdf, plus trd.md (plain Markdown for AI/reading).

Usage: python3 build_trd.py [--university "Name"] [--team "Team"] [--code-pt 6.2] [--cols 4] [--no-render]
       python3 build_trd.py --check   # trd.md is current and its code passes all self-tests

The .py files stay the source of truth. This script strips each file's
self-test block (`if __name__ == "__main__":` onward), splits the code at its
`# ==== X ====` banners, and appends the contest-relevant README sections.
After rendering it checks the ICPC TRD rules: <= 25 pages, A4, one-sided,
university name top-left and page number top-right on every page.
"""
import argparse
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MAX_PAGES = 25

# Glyphs the monospace font lacks -> ASCII.
ASCII = {
    "·": "*", "–": "-", "—": "-", "−": "-", "≤": "<=", "≥": ">=", "→": "->",
    "≡": "==", "√": "sqrt", "⇔": "<=>", "⇒": "=>", "φ": "phi", "…": "...",
    "π": "pi", "³": "^3", "²": "^2", "ö": "o", "ő": "o", "⊆": "subset of",
    "Σ": "Sum", "Π": "Prod", "⌈": "ceil(", "⌉": ")", "⌊": "floor(", "⌋": ")",
    "Θ": "Theta", "Ω": "Omega", "×": "x", "≠": "!=", "∑": "Sum",
}

TITLES = {
    "00": "Core: template, I/O, recursion, testing",
    "01": "Techniques: search, greedy, offline, heuristics",
    "02": "Data structures",
    "03": "Dynamic programming",
    "04": "Graphs",
    "05": "Flows and matching",
    "06": "Trees",
    "07": "Strings",
    "08": "Number theory, combinatorics, games",
    "09": "Linear algebra and polynomials",
    "10": "Geometry",
}

# README sections that are useful during a contest (the rest stays in README.md).
README_KEEP = ("## 2.", "## 3.")  # 1, 4, 5 duplicate notes in 00_core.py / 01_techniques.py
README_SUBS_ONLY = ("## 2.",)  # checklist duplicates 00_core.py; keep only its ### parts
README_DROP_SUB = ("### Formats & verdicts",)

BANNER = re.compile(r"^#\s*====\s*(.*?)\s*====\s*$")


def asciify(s):
    for k, v in ASCII.items():
        s = s.replace(k, v)
    return s


def tidy(lines):
    """Compact code for print: 2-space indents; blank lines only before top-level defs."""
    out = []
    for i, ln in enumerate(lines):
        ln = ln.rstrip()
        if not ln:
            nxt = next((x for x in lines[i + 1:] if x.strip()), "")
            if not out or not out[-1] or not nxt.startswith(("def ", "class ", "@")):
                continue
        n = len(ln) - len(ln.lstrip(" "))
        body = re.sub(r"(\S) {3,}#", r"\1  #", ln[n:])  # drop right-aligned comment padding
        out.append(" " * (n // 2) + body)
    while out and not out[-1]:
        out.pop()
    return out


def plain(lines):
    """Readable code for the Markdown edition: original indents, single blank lines."""
    out = []
    for ln in lines:
        ln = ln.rstrip()
        if ln or (out and out[-1]):
            out.append(ln)
    while out and not out[-1]:
        out.pop()
    return out


def code_sections(path, compact=True):
    src = open(path, encoding="utf-8").read()
    if compact:
        src = asciify(src)
    src = re.split(r"\n(?:if __name__ ==|# =+ *[Tt]ests?\b)", src)[0]
    lines = src.splitlines()
    if lines and lines[0].startswith("#!"):
        lines = lines[1:]
    head, sections, cur = [], [], None
    for ln in lines:
        m = BANNER.match(ln)
        if m:
            cur = [m.group(1), []]
            sections.append(cur)
        elif cur is None:
            head.append(ln)
        else:
            cur[1].append(ln)
    # Keep imports, constants and convention comments; drop the "# NN_name.py -- title" line.
    head_src = "\n".join(head)
    head_src = re.sub(r'^\s*(""".*?"""|\'\'\'.*?\'\'\')', "", head_src, flags=re.S)
    name = os.path.basename(path)
    head = [ln for ln in head_src.splitlines() if not ln.startswith("# " + name)]
    head = (tidy if compact else plain)(head)
    fmt = tidy if compact else plain
    return head, [(t, fmt(b)) for t, b in sections if fmt(b)]


def table_to_bullets(lines):
    rows = [ln for ln in lines if ln.startswith("|") and not re.match(r"^\|[-|: ]+\|$", ln)]
    cells = [[c.strip() for c in r.strip("|").split("|")] for r in rows]
    hdr, body = cells[0], cells[1:]
    return ["- " + "; ".join(f"{h}: {c}" if i else c for i, (h, c) in enumerate(zip(hdr, r)) if c)
            for r in body]


def readme_sections():
    text = asciify(open(os.path.join(HERE, "README.md"), encoding="utf-8").read())
    out, keep, skip_sub, buf = [], False, False, []
    for ln in text.splitlines() + ["## END"]:
        if ln.startswith("|"):
            buf.append(ln)
            continue
        if buf:
            if keep and not skip_sub:
                out.extend(table_to_bullets(buf))
            buf = []
        if ln.startswith("## "):
            keep = ln.startswith(README_KEEP)
            skip_sub = ln.startswith(README_SUBS_ONLY)
            if keep and not skip_sub:
                out.append(ln)  # level 2 under the cheat-sheet heading
            continue
        if ln.startswith("### "):
            skip_sub = ln.startswith(README_DROP_SUB)
            if keep and not skip_sub:
                out.append(ln)
            continue
        if keep and not skip_sub and ln.strip() != "---":
            out.append(ln)
    return out


def build_qmd(university, team, fs=6.4, cols=4):
    lh = round(fs * 1.05, 2)
    header_left = university + (f" -- {team}" if team else "")
    preamble = rf"""
\usepackage{{multicol}}
\setlength{{\columnsep}}{{3mm}}
\setlength{{\columnseprule}}{{0.2pt}}
\usepackage{{fancyhdr}}
\pagestyle{{fancy}}
\fancyhf{{}}
\fancyhead[L]{{\small\bfseries {header_left}}}
\fancyhead[C]{{\small Team Reference Document (Python)}}
\fancyhead[R]{{\small\bfseries \thepage}}
\renewcommand{{\headrulewidth}}{{0.3pt}}
\fancypagestyle{{plain}}{{\fancyhf{{}}\fancyhead[L]{{\small\bfseries {header_left}}}\fancyhead[C]{{\small Team Reference Document (Python)}}\fancyhead[R]{{\small\bfseries \thepage}}}}
\setmonofont[Scale=1.0,Extension=.otf,BoldFont=Inconsolatazi4-Bold]{{Inconsolatazi4-Regular}}
\usepackage{{titlesec}}
\newcommand{{\secbox}}[1]{{\colorbox{{black}}{{\parbox{{\dimexpr\columnwidth-2\fboxsep}}{{\color{{white}}#1}}}}}}
\titleformat{{\section}}[block]{{\normalfont\small\bfseries}}{{}}{{0pt}}{{\secbox}}
\titlespacing*{{\section}}{{0pt}}{{4pt}}{{2pt}}
\titleformat{{\subsection}}{{\normalfont\scriptsize\bfseries}}{{}}{{0pt}}{{\rule[0.5ex]{{3pt}}{{0.4pt}}\ }}
\titlespacing*{{\subsection}}{{0pt}}{{1pt}}{{0pt}}
\titleformat{{\subsubsection}}{{\normalfont\scriptsize\bfseries}}{{}}{{0pt}}{{}}
\titlespacing*{{\subsubsection}}{{0pt}}{{1pt}}{{0pt}}
\usepackage{{enumitem}}
\setlist{{nosep,leftmargin=*}}
\setlength{{\parskip}}{{1pt}}
\setlength{{\parindent}}{{0pt}}
\usepackage{{fvextra}}
\newcommand{{\codesize}}{{\fontsize{{{fs}}}{{{lh}}}\selectfont}}
\renewenvironment{{Shaded}}{{}}{{}}
\RecustomVerbatimEnvironment{{Highlighting}}{{Verbatim}}{{commandchars=\\\{{\}},breaklines,breakanywhere,breaksymbolleft={{\tiny$\hookrightarrow$}},fontsize=\codesize}}
\raggedcolumns
\AtBeginDocument{{\renewcommand{{\CommentTok}}[1]{{\textcolor[HTML]{{2E6B30}}{{#1}}}}}}
"""
    parts = [
        "---",
        "format:",
        "  pdf:",
        "    pdf-engine: xelatex",
        "    documentclass: article",
        "    classoption: [oneside]",
        "    papersize: a4",
        "    geometry: [landscape, margin=6mm, top=10mm, headsep=2mm, headheight=10pt]",
        "    mainfont: Noto Sans SemiCondensed",
        "    mainfontoptions: [Scale=0.95]",
        "    fontsize: 8pt",
        "    highlight-style: github",
        "    code-block-bg: false",
        "    code-block-border-left: false",
        "    number-sections: false",
        "    toc: false",
        "    colorlinks: false",
        "    include-in-header:",
        "      text: |",
        *["        " + ln for ln in preamble.strip().splitlines()],
        "---",
        "",
        r"```{=latex}",
        r"\begin{multicols*}{%d}" % cols,
        r"\thispagestyle{fancy}",
        r"{\bfseries\large Team Reference Document (Python 3 / PyPy)}\par",
        r"\scriptsize Built from tested sources (\texttt{python3 run\_tests.sh}). Sections: "
        + "; ".join(f"{k} {v}" for k, v in TITLES.items()) + r"; Cheat sheet.\par",
        r"```",
        "",
    ]
    for path in sorted(glob.glob(os.path.join(HERE, "[0-9][0-9]_*.py"))):
        key = os.path.basename(path)[:2]
        head, sections = code_sections(path)
        parts.append(f"# {key}  {TITLES.get(key, os.path.basename(path))}")
        parts.append("")
        if head:
            parts += ["```python", *head, "```", ""]
        for title, body in sections:
            parts.append("## " + title.replace("*", r"\*").replace("_", r"\_"))
            parts.append("")
            parts += ["```python", *body, "```", ""]
    parts += ["# Cheat sheet", ""] + readme_sections()
    parts += ["", r"```{=latex}", r"\end{multicols*}", r"```", ""]
    return "\n".join(parts)


def build_md():
    """Plain Markdown edition (no page limit): full code without tests, plus the whole README."""
    parts = [
        "# Python CP Reference (Team Reference Document, Markdown edition)",
        "",
        "Generated by `build_trd.py` from the self-tested files in `reference/` (tests omitted). "
        "Python 3 / PyPy, stdlib only. Conventions: indices 0-based, ranges half-open [l, r) unless noted; "
        "section numbers in `###` headings match `competitive-programming-topics.md`; "
        "`# see <file>: <name>` points to the single implementation of a technique.",
        "",
        "## Contents",
        "",
    ]
    files = sorted(glob.glob(os.path.join(HERE, "[0-9][0-9]_*.py")))
    for path in files:
        name = os.path.basename(path)
        title = TITLES.get(name[:2], name)
        anchor = re.sub(r"[^a-z0-9 _-]", "", f"{name} {title}".lower()).replace(" ", "-")
        parts.append(f"- [{name}: {title}](#{anchor})")
    parts += ["- [Topic map and cheat sheet](#topic-map-and-cheat-sheet)", ""]
    for path in files:
        name = os.path.basename(path)
        head, sections = code_sections(path, compact=False)
        parts += [f"## {name}: {TITLES.get(name[:2], name)}", ""]
        if head:
            parts += ["```python", *head, "```", ""]
        for title, body in sections:
            parts += [f"### {title}", "", "```python", *body, "```", ""]
    parts += ["## Topic map and cheat sheet", ""]
    text = open(os.path.join(HERE, "README.md"), encoding="utf-8").read()
    for ln in text.splitlines()[1:]:  # drop the README title
        if ln.startswith("#"):
            ln = "#" + ln
        parts.append(ln)
    return "\n".join(parts).rstrip() + "\n"


def check_md(md_text):
    """Run each file's self-tests against the code in the Markdown edition (proves it is complete)."""
    import contextlib
    import io
    chunks = re.split(r"\n## (\d\d_\w+\.py): ", md_text)
    ok = True
    for i in range(1, len(chunks), 2):
        name, body = chunks[i], chunks[i + 1].split("\n## Topic map")[0]
        code = "\n".join(re.findall(r"```python\n(.*?)```", body, re.S))
        src = open(os.path.join(HERE, name), encoding="utf-8").read()
        test = src[src.index("\nif __name__"):]
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                exec(compile(code + test, name, "exec"), {"__name__": "__main__"})
            passed = "ALL OK" in buf.getvalue()
        except Exception as e:  # noqa: BLE001 - report any failure
            passed = False
            print(name, "error:", e)
        print(("PASS" if passed else "FAIL"), name, "(from trd.md)")
        ok &= passed
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--university", default="University Name")
    ap.add_argument("--team", default="")
    ap.add_argument("--code-pt", type=float, default=6.2, help="code font size (pt)")
    ap.add_argument("--cols", type=int, default=4)
    ap.add_argument("--no-render", action="store_true")
    ap.add_argument("--check", action="store_true",
                    help="verify trd.md is up to date and passes every file's tests; writes nothing")
    a = ap.parse_args()
    if a.check:
        md_text = build_md()
        with open(os.path.join(HERE, "trd.md"), encoding="utf-8") as f:
            fresh = f.read() == md_text
        print("trd.md up to date" if fresh else "FAIL: trd.md is stale; run build_trd.py")
        sys.exit(0 if check_md(md_text) and fresh else 1)
    qmd = os.path.join(HERE, "trd.qmd")
    with open(qmd, "w", encoding="utf-8") as f:
        f.write(build_qmd(a.university, a.team, a.code_pt, a.cols))
    print("wrote", qmd)
    md = os.path.join(HERE, "trd.md")
    with open(md, "w", encoding="utf-8") as f:
        f.write(build_md())
    print("wrote", md)
    if a.no_render:
        return
    subprocess.run(["quarto", "render", qmd, "--to", "pdf"], check=True, cwd=HERE)
    pdf = os.path.join(HERE, "trd.pdf")
    info = subprocess.run(["pdfinfo", pdf], capture_output=True, text=True).stdout
    m_pages, m_size = re.search(r"Pages:\s+(\d+)", info), re.search(r"Page size:\s+(.*)", info)
    if not m_pages or not m_size:
        sys.exit("FAIL: could not read pdfinfo output (is poppler-utils installed?)")
    pages, size = int(m_pages.group(1)), m_size.group(1)
    print(f"pages: {pages} (limit {MAX_PAGES}); page size: {size}")
    if pages > MAX_PAGES:
        sys.exit(f"FAIL: {pages} pages > {MAX_PAGES}")
    print("OK: within ICPC TRD limits")


if __name__ == "__main__":
    main()
