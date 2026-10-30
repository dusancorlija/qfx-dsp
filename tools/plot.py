#!/usr/bin/env python3
"""Render a CSV of (x, series...) as a clean SVG line chart. No dependencies.

    python3 tools/plot.py docs/noise_shaping.csv docs/noise_shaping.svg
"""
import csv
import math
import sys

W, H = 760, 380
L, R, T, B = 64, 24, 48, 52
COLORS = ["#9aa0a6", "#e2531a"]
LABELS = {"round_db": "Rounded", "shaped_db": "Error feedback"}


def nice_ticks(lo, hi, n=6):
    span = hi - lo
    step = 10 ** math.floor(math.log10(span / n))
    for m in (1, 2, 5, 10):
        if span / (step * m) <= n:
            step *= m
            break
    start = math.ceil(lo / step) * step
    return [start + i * step for i in range(int((hi - start) / step) + 1)]


def main(src, dst, title="Output error spectrum, 30 Hz low-pass, Q31 (fs = 48 kHz)"):
    with open(src) as f:
        rows = list(csv.reader(f))
    head, data = rows[0], [[float(v) for v in r] for r in rows[1:]]
    xs = [r[0] for r in data]
    series = [[r[i] for r in data] for i in range(1, len(head))]

    x0, x1 = math.log10(min(xs)), math.log10(max(xs))
    flat = [v for s in series for v in s]
    y0, y1 = math.floor(min(flat) / 10) * 10, math.ceil(max(flat) / 10) * 10

    def px(x):
        return L + (math.log10(x) - x0) / (x1 - x0) * (W - L - R)

    def py(y):
        return T + (y1 - y) / (y1 - y0) * (H - T - B)

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'font-family="Inter, -apple-system, Segoe UI, sans-serif" font-size="12">',
        f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
        f'<text x="{L}" y="26" font-size="14" font-weight="600" fill="#15171a">{title}</text>',
    ]
    for y in nice_ticks(y0, y1):
        out.append(f'<line x1="{L}" x2="{W - R}" y1="{py(y):.1f}" y2="{py(y):.1f}" stroke="#eceae6"/>')
        out.append(f'<text x="{L - 8}" y="{py(y) + 4:.1f}" text-anchor="end" fill="#7a7e85">{y:g}</text>')
    decade = 10 ** math.floor(x0)
    while decade <= 10 ** x1:
        for m in (1, 2, 5):
            x = decade * m
            if 10 ** x0 <= x <= 10 ** x1:
                out.append(f'<line x1="{px(x):.1f}" x2="{px(x):.1f}" y1="{T}" y2="{H - B}" stroke="#f3f1ed"/>')
                label = f"{x / 1000:g}k" if x >= 1000 else f"{x:g}"
                out.append(f'<text x="{px(x):.1f}" y="{H - B + 18}" text-anchor="middle" fill="#7a7e85">{label}</text>')
        decade *= 10
    out.append(f'<text x="{(L + W - R) / 2}" y="{H - 12}" text-anchor="middle" fill="#3d4148">Frequency (Hz)</text>')
    out.append(f'<text transform="translate(16 {(T + H - B) / 2}) rotate(-90)" text-anchor="middle" fill="#3d4148">Error power (dB)</text>')

    for i, s in enumerate(series):
        pts = " ".join(f"{px(x):.1f},{py(y):.1f}" for x, y in zip(xs, s))
        out.append(f'<polyline points="{pts}" fill="none" stroke="{COLORS[i % 2]}" stroke-width="1.8" stroke-linejoin="round"/>')
        lx = W - R - 150
        ly = T + 8 + i * 18
        out.append(f'<line x1="{lx}" x2="{lx + 18}" y1="{ly}" y2="{ly}" stroke="{COLORS[i % 2]}" stroke-width="2.4"/>')
        out.append(f'<text x="{lx + 26}" y="{ly + 4}" fill="#15171a">{LABELS.get(head[i + 1], head[i + 1])}</text>')

    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
