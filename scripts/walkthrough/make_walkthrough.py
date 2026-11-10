"""Generate WALKTHROUGH.html (repo root) from the committed data.

Run from anywhere: python scripts/walkthrough/make_walkthrough.py
The inline charts are drawn from singlet_tower_deep.json, thooft_limit.json
and reference_spectrum_polished.json, and the errata count from ERRATA.md;
the prose lives in template.html.  scripts/test_walkthrough.py checks that
the committed WALKTHROUGH.html is what this produces.
"""
import json
import math
import re
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.special import loggamma

ROOT = Path(__file__).resolve().parents[2]


def load(p):
    with open(ROOT / p) as fh:
        return json.load(fh)


deep = load("duality/singlet_tower_deep.json")["levels"]
ref = load("thooft-target/reference_spectrum_polished.json")["rows"]
tl = load("syk-self-averaging/thooft_limit.json")["coefficients"]

PAL = ["c1", "c2", "c3", "c4", "c5", "c6", "c7"]   # O2..O8


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ---------------------------------------------------------------------------
# SVG helpers
# ---------------------------------------------------------------------------
class Plot:
    def __init__(self, w, h, x0, x1, y0, y1, logx=False, logy=False,
                 m=(52, 16, 18, 40)):
        self.w, self.h = w, h
        self.l, self.r, self.t, self.b = m
        self.x0, self.x1, self.y0, self.y1 = x0, x1, y0, y1
        self.logx, self.logy = logx, logy
        self.out = []

    def X(self, x):
        if self.logx:
            x, a, b = math.log10(x), math.log10(self.x0), math.log10(self.x1)
        else:
            a, b = self.x0, self.x1
        return self.l + (x - a) / (b - a) * (self.w - self.l - self.r)

    def Y(self, y):
        if self.logy:
            y, a, b = math.log10(y), math.log10(self.y0), math.log10(self.y1)
        else:
            a, b = self.y0, self.y1
        return self.h - self.b - (y - a) / (b - a) * (self.h - self.t - self.b)

    def axes(self, xt, yt, xl="", yl="", xfmt=str, yfmt=str):
        o = self.out
        L, R, T, B = self.l, self.w - self.r, self.t, self.h - self.b
        for v in yt:
            y = self.Y(v)
            o.append(f'<line class="grid" x1="{L}" x2="{R}" y1="{y:.1f}" y2="{y:.1f}"/>')
            o.append(f'<text class="tick" x="{L - 6}" y="{y + 4:.1f}" text-anchor="end">{yfmt(v)}</text>')
        for v in xt:
            x = self.X(v)
            o.append(f'<text class="tick" x="{x:.1f}" y="{B + 16}" text-anchor="middle">{xfmt(v)}</text>')
        o.append(f'<line class="axis" x1="{L}" x2="{R}" y1="{B}" y2="{B}"/>')
        o.append(f'<line class="axis" x1="{L}" x2="{L}" y1="{T}" y2="{B}"/>')
        if xl:
            o.append(f'<text class="lab" x="{(L + R) / 2:.1f}" y="{self.h - 6}" text-anchor="middle">{xl}</text>')
        if yl:
            o.append(f'<text class="lab" transform="translate(13,{(T + B) / 2:.1f}) rotate(-90)" text-anchor="middle">{yl}</text>')

    def line(self, xs, ys, cls, dash=False):
        pts = " ".join(f"{self.X(x):.1f},{self.Y(y):.1f}" for x, y in zip(xs, ys))
        d = ' stroke-dasharray="5 4"' if dash else ""
        self.out.append(f'<polyline class="ln {cls}" points="{pts}"{d}/>')

    def dots(self, xs, ys, cls, r=3.6):
        for x, y in zip(xs, ys):
            self.out.append(f'<circle class="dot {cls}" cx="{self.X(x):.1f}" cy="{self.Y(y):.1f}" r="{r}"/>')

    def text(self, x, y, s, cls="ann", anchor="start"):
        self.out.append(f'<text class="{cls}" x="{self.X(x):.1f}" y="{self.Y(y):.1f}" text-anchor="{anchor}">{s}</text>')

    def svg(self, label):
        return (f'<svg viewBox="0 0 {self.w} {self.h}" role="img" aria-label="{esc(label)}">'
                + "".join(self.out) + "</svg>")


# ---------------------------------------------------------------------------
# 1. the two limits
# ---------------------------------------------------------------------------
def limits_svg():
    w, h = 600, 300
    L, B, R, T = 70, 245, 570, 30
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="The two limits">']
    o.append('<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerUnits="userSpaceOnUse" markerWidth="10" markerHeight="10" orient="auto-start-reverse"><path d="M0,1 L10,5 L0,9 z" class="arrowhead"/></marker></defs>')
    o.append(f'<line class="axis" x1="{L}" y1="{B}" x2="{R}" y2="{B}" marker-end="url(#arr)"/>')
    o.append(f'<line class="axis" x1="{L}" y1="{B}" x2="{L}" y2="{T}" marker-end="url(#arr)"/>')
    o.append(f'<text class="lab" x="{R}" y="{B + 24}" text-anchor="end">1/p  →</text>')
    o.append(f'<text class="lab" x="{L - 12}" y="{T + 10}" text-anchor="end">λ</text>')
    o.append(f'<text class="tick" x="{L - 12}" y="{T + 25}" text-anchor="end">= 2p²/N</text>')
    for y in (85, 135, 185):
        o.append(f'<line class="ln c1" x1="{R - 20}" y1="{y}" x2="{L + 8}" y2="{y}" marker-end="url(#arr)"/>')
    o.append(f'<text class="ann c1t" x="{L + 16}" y="{72}">double scaling: λ fixed, N and p → ∞</text>')
    o.append(f'<text class="ann" x="{L + 16}" y="{112}">the paper\'s flat-space limit;</text>')
    o.append(f'<text class="ann" x="{L + 16}" y="{127}">chord diagrams are exact here</text>')
    for x in (440, 485, 530):
        o.append(f'<line class="ln c2" x1="{x}" y1="{T + 20}" x2="{x}" y2="{B - 8}" stroke-dasharray="6 4" marker-end="url(#arr)"/>')
    o.append(f'<text class="ann c2t" x="{428}" y="{208}" text-anchor="end">fixed p, N → ∞ (so λ → 0)</text>')
    o.append(f'<text class="ann" x="{428}" y="{224}" text-anchor="end">the paper\'s \'t Hooft limit, ᾱ = p²</text>')
    o.append(f'<circle class="dot c3" cx="{L}" cy="{B}" r="7"/>')
    o.append(f'<text class="ann c3t" x="{L + 14}" y="{B + 24}">λ → 0, p → ∞: the semiclassical corner</text>')
    o.append("</svg>")
    return "".join(o)


# ---------------------------------------------------------------------------
# 2. chord diagram
# ---------------------------------------------------------------------------
def chord_svg():
    w, h = 560, 300
    cx, cy, r = 150, 150, 110
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="A chord diagram">']
    o.append(f'<circle class="ring" cx="{cx}" cy="{cy}" r="{r}"/>')

    def pt(a):
        a = math.radians(a)
        return cx + r * math.cos(a), cy - r * math.sin(a)

    def chord(a, b, cls, dash=False):
        (x1, y1), (x2, y2) = pt(a), pt(b)
        d = ' stroke-dasharray="6 4"' if dash else ""
        o.append(f'<path class="ln {cls}" d="M{x1:.1f},{y1:.1f} Q{cx},{cy} {x2:.1f},{y2:.1f}"{d}/>')
        for x, y in ((x1, y1), (x2, y2)):
            o.append(f'<circle class="dot {cls}" cx="{x:.1f}" cy="{y:.1f}" r="4"/>')

    chord(160, 20, "c1")
    chord(110, 250, "c1")
    chord(60, 300, "c1")
    chord(200, 340, "c2", dash=True)
    o.append(f'<text class="ann" x="{cx - 120}" y="{cy + 140}">H-chords (solid), one matter chord (dashed)</text>')
    tx = 296
    rows = [("two H-chords crossing", 'weight × q,  q = e<tspan dy="-7" font-size="11">−λ</tspan>'),
            ("an H-chord crossing the matter chord", 'weight × q<tspan dy="-7" font-size="11">Δ</tspan>'),
            ("sum over every diagram", "= the exact N = ∞ moment"),
            ("Δ of a matter chord", "is its \"dimension\"")]
    y = 60
    for a, b in rows:
        o.append(f'<text class="ann" x="{tx}" y="{y}">{a}</text>')
        o.append(f'<text class="big" x="{tx}" y="{y + 22}">{b}</text>')
        y += 60
    o.append("</svg>")
    return "".join(o)


# ---------------------------------------------------------------------------
# 3. the tower: shares at lambda = 2, A_2..A_8
# ---------------------------------------------------------------------------
def shares_svg():
    ch = deep["8"]["scan_lam2"]["channels"]
    w, rowh = 620, 34
    L = 46
    h = 30 + rowh * 7 + 40
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="Sector shares at lambda = 2">']
    W = 420
    for i, n in enumerate(range(2, 9)):
        y = 20 + i * rowh
        sh = ch[str(n)]["sector_share"]
        x = L
        o.append(f'<text class="tick" x="{L - 8}" y="{y + 17}" text-anchor="end">A{chr(0x2080 + n)}</text>')
        for m in range(2, n + 1):
            v = sh[str(m)]
            ww = v * W
            o.append(f'<rect class="bar {PAL[m - 2]}" x="{x:.2f}" y="{y}" width="{max(ww, 0.0):.2f}" height="24"><title>O{m}: {v * 100:.4g}%</title></rect>')
            if ww > 34:
                o.append(f'<text class="barlab" x="{x + ww / 2:.1f}" y="{y + 16}" text-anchor="middle">{v * 100:.0f}%</text>')
            x += ww
        new = sh[str(n)]
        o.append(f'<text class="ann" x="{L + W + 10}" y="{y + 16}">new O{chr(0x2080 + n)}: {new * 100:.2g}%</text>' if n > 2 else
                 f'<text class="ann" x="{L + W + 10}" y="{y + 16}">pure O₂</text>')
    y = 20 + 7 * rowh + 12
    x = L
    for m in range(2, 9):
        o.append(f'<rect class="bar {PAL[m - 2]}" x="{x}" y="{y}" width="12" height="12"/>')
        o.append(f'<text class="tick" x="{x + 16}" y="{y + 10}">O{chr(0x2080 + m)}</text>')
        x += 60
    o.append("</svg>")
    return "".join(o)


def survival_svg():
    rows = sorted(deep["8"]["share_scaling"], key=lambda r: r["lam"])
    lams = [r["lam"] for r in rows]
    p = Plot(620, 340, 5e-4, 0.3, 1e-8, 1.5, logx=True, logy=True, m=(62, 120, 18, 44))
    p.axes([1e-3, 1e-2, 1e-1], [1e-8, 1e-6, 1e-4, 1e-2, 1],
           xl="λ (chord coupling)", yl="share of A₈'s dynamical weight",
           xfmt=lambda v: f"{v:g}", yfmt=lambda v: {1e-8: "10⁻⁸", 1e-6: "10⁻⁶", 1e-4: "10⁻⁴", 1e-2: "1%", 1: "100%"}[v])
    labels = []
    for m in range(2, 9):
        ys = [r["shares"]["8"][str(m)] for r in rows]
        cls = PAL[m - 2]
        p.line(lams, ys, cls, dash=(m % 2 == 0))
        p.dots(lams, ys, cls, r=3)
        labels.append([p.Y(ys[-1]) + 4, cls, f'O{chr(0x2080 + m)}{" fades ∝ λ" if m % 2 == 0 else " stays"}'])
    labels.sort()
    for i in range(1, len(labels)):           # keep labels 15 px apart
        labels[i][0] = max(labels[i][0], labels[i - 1][0] + 15)
    for y, cls, txt in labels:
        p.out.append(f'<text class="ann {cls}t" x="{p.X(lams[-1]) + 8:.1f}" y="{y:.1f}">{txt}</text>')
    return p.svg("Odd chords survive, even chords fade")


# ---------------------------------------------------------------------------
# 4. horizon: ladder densities, and mass spectra
# ---------------------------------------------------------------------------
def ladder_svg():
    p = Plot(420, 280, 0, 8, 0, 1.05, m=(50, 14, 14, 40))
    p.axes([0, 2, 4, 6, 8], [0, 0.25, 0.5, 0.75, 1.0],
           xl="boost frequency ω / a", yl="ρ(ω) / ρ(0)",
           yfmt=lambda v: f"{v:g}")
    x = np.linspace(0, 8, 300)
    for d, cls in ((2, "c1"), (3, "c2"), (5, "c4"), (7, "c6")):
        y = np.exp(2 * (loggamma(d + 1j * x).real - loggamma(d).real))
        p.line(x, y, cls)
    for i, (d, cls) in enumerate(((2, "c1"), (3, "c2"), (5, "c4"), (7, "c6"))):
        yy = 0.92 - 0.11 * i
        p.out.append(f'<line class="ln {cls}" x1="{p.X(5.6):.1f}" x2="{p.X(6.4):.1f}" y1="{p.Y(yy) - 4:.1f}" y2="{p.Y(yy) - 4:.1f}"/>')
        p.out.append(f'<text class="ann" x="{p.X(6.6):.1f}" y="{p.Y(yy):.1f}">Δ = {d}</text>')
    return p.svg("Quasinormal ladder densities")


def masses_svg():
    vals = [math.sqrt(float(r["value"])) for r in ref][:20]
    w, h = 620, 230
    p = Plot(w, h, 0, 4.6, 0, 1.08, m=(52, 16, 26, 42))
    p.axes([0, 1, 2, 3, 4], [], xl="mass M (units: √(πg²) for 't Hooft; schematic for the density)")
    M = np.linspace(0, 4.6, 300)
    for d, cls in ((2, "c1"), (3, "c2")):
        y = (M / 4.6) ** (2 * d - 1)
        p.line(M, y, cls)
    p.out.append(f'<text class="ann c1t" x="{p.X(0.05):.1f}" y="{p.Y(0.92):.1f}">N = ∞ singlet channels: gapless density ∝ M<tspan dy="-6" font-size="10">2Δ−1</tspan><tspan dy="6"> (Δ = 2 blue, Δ = 3 orange)</tspan></text>')
    for v in vals:
        if v <= 4.6:
            p.out.append(f'<line class="stick" x1="{p.X(v):.1f}" x2="{p.X(v):.1f}" y1="{p.Y(0)}" y2="{p.Y(0.42):.1f}"/>')
    p.out.append(f'<text class="ann c5t" x="{p.X(0.05):.1f}" y="{p.Y(0.62):.1f}">\'t Hooft model: discrete lines, the first at M₀ = 0.859</text>')
    return p.svg("Gapless density against the discrete 't Hooft tower")


# ---------------------------------------------------------------------------
# 5. fixed p
# ---------------------------------------------------------------------------
def A(p, n, k):
    return Fraction(tl[f"{p},{n},{k}"]["A"])


def fixedp_svg():
    ps = [4, 6, 8, 10, 12]
    w2 = [float(A(p, 2, 2) / A(p, 2, 0)) / (8 * p * p) for p in ps]
    ku = [float(A(p, 2, 4) * A(p, 2, 0) / A(p, 2, 2) ** 2) / 3.5 for p in ps]
    a3 = [float(A(p, 3, 4) / A(p, 3, 2)) / (40 * p * p) for p in ps]
    pl = Plot(560, 300, 3, 40, 0.4, 1.05, m=(52, 150, 18, 42))
    pl.axes([4, 8, 12, 20, 30, 40], [0.4, 0.6, 0.8, 1.0], xl="interaction order p (N → ∞ first)",
            yl="exact value ÷ ladder value", yfmt=lambda v: f"{v:g}")
    xs = np.linspace(3.6, 40, 200)
    curves = [((1 - 2 / xs), w2, "c1", "A₂ ⟨ω²⟩: 1 − 2/p"),
              ((3.5 - 1 / (2 * (xs - 2))) / 3.5, ku, "c2", "A₂ kurtosis"),
              ((10 * xs - 21) / (10 * xs), a3, "c4", "A₃ ⟨ω⁴⟩/⟨ω²⟩")]
    labels = []
    for cv, pts, cls, lab in curves:
        pl.line(xs, cv, cls, dash=True)
        pl.dots(ps, pts, cls, r=4)
        labels.append([pl.Y(cv[-1]) + 4, cls, lab])
    labels.sort()
    for i in range(1, len(labels)):
        labels[i][0] = max(labels[i][0], labels[i - 1][0] + 16)
    for y, cls, lab in labels:
        pl.out.append(f'<text class="ann {cls}t" x="{pl.X(40) + 8:.1f}" y="{y:.1f}">{lab}</text>')
    pl.out.append(f'<line class="axis" x1="{pl.X(3)}" x2="{pl.X(40)}" y1="{pl.Y(1):.1f}" y2="{pl.Y(1):.1f}" stroke-dasharray="2 3"/>')
    pl.out.append(f'<text class="ann" x="{pl.X(14):.1f}" y="{pl.Y(1) - 6:.1f}">λ → 0 ladder value</text>')
    return pl.svg("Fixed-p laws approach the ladders")


# ---------------------------------------------------------------------------
# 6. architecture
# ---------------------------------------------------------------------------
def stack_svg():
    layers = [
        ("7", "Where it stands", "verdict, open problems", "c7"),
        ("6", "The dS / OTOC axis", "scrambling sign, fold, antipode", "c6"),
        ("5", "The 't Hooft limit at fixed p", "exact shape laws", "c5"),
        ("4", "Calibration: boost vs mass", "Rindler horizon, gapless spectrum", "c4"),
        ("3", "The singlet tower, n ≤ 8", "conserved part + matter chords", "c3"),
        ("2", "Chords at N = ∞", "exact moments, two theorems", "c2"),
        ("1", "SYK at finite N", "identities, exact Wick sums, ED", "c1"),
        ("0", "The yardstick", "the 't Hooft meson spectrum", "c0"),
    ]
    w, rh = 560, 40
    h = rh * len(layers) + 10
    o = [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="The layers of the project">']
    for i, (k, a, b, cls) in enumerate(layers):
        y = 4 + i * rh
        inset = (len(layers) - 1 - i) * 6
        o.append(f'<a href="#layer{k}"><rect class="layer {cls}" x="{inset}" y="{y}" width="{w - 2 * inset}" height="{rh - 6}" rx="7"/>'
                 f'<text class="laytxt" x="{inset + 14}" y="{y + 22}">{k} · {a}</text>'
                 f'<text class="laysub" x="{w - inset - 14}" y="{y + 22}" text-anchor="end">{b}</text></a>')
    o.append("</svg>")
    return "".join(o)


SVG = dict(limits=limits_svg(), chord=chord_svg(), shares=shares_svg(),
           survival=survival_svg(), ladder=ladder_svg(), masses=masses_svg(),
           fixedp=fixedp_svg(), stack=stack_svg())

def errata_count() -> int:
    """The highest entry number in ERRATA.md (entries are E1, E2, ...)."""
    text = (ROOT / "ERRATA.md").read_text()
    return max(int(n) for n in re.findall(r"^\*\*E(\d+)\.", text, re.M))


def build_html() -> str:
    html = (Path(__file__).parent / "template.html").read_text()
    for k, v in SVG.items():
        html = html.replace("{{" + k + "}}", v)
    html = html.replace("{{errata}}", str(errata_count()))
    assert "{{" not in html, "unfilled placeholder in template.html"
    return html


def main(out=ROOT / "WALKTHROUGH.html") -> None:
    html = build_html()
    Path(out).write_text(html)
    print("wrote", len(html), "bytes; first meson", math.sqrt(float(ref[0]["value"])))


if __name__ == "__main__":
    main()
