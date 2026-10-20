#!/usr/bin/env python3
"""Render the Clinker HTML plots (06_clinker_windows.py --run) to PNG with headless Chromium.

Usage: 07_clinker_png.py <clinker_dir> [--scale 2] [--browser /path/to/chrome]
Writes <variant>.png next to each <variant>.html. Needs `pip install playwright` and a Chromium
(PLAYWRIGHT_BROWSERS_PATH is honoured; otherwise pass --browser). The image is cropped to the figure (no control panel).
"""
import argparse, sys
from pathlib import Path

from playwright.sync_api import sync_playwright


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("clinker_dir"); ap.add_argument("--scale", type=float, default=2)
    ap.add_argument("--browser", default=None, help="Chromium executable (default: Playwright's own)")
    a = ap.parse_args()
    fit = """(pad) => {                      // crop the figure: drop the control panel, fit the viewBox to the drawing
        const svg = document.querySelector('#root-svg'), g = svg.querySelector('g.clusterMapG');
        const r = g.getBoundingClientRect(), o = svg.getBoundingClientRect();
        const x = r.x - o.x - pad, y = r.y - o.y - pad, w = r.width + 2 * pad, h = r.height + 2 * pad;
        svg.setAttribute('viewBox', `${x} ${y} ${w} ${h}`); svg.setAttribute('width', w); svg.setAttribute('height', h);
        svg.style.background = 'white'; document.body.replaceChildren(svg);
        return [w, h];
    }"""
    with sync_playwright() as p:
        b = p.chromium.launch(executable_path=a.browser) if a.browser else p.chromium.launch()
        for html in sorted(Path(a.clinker_dir).glob("*.html")):
            pg = b.new_page(viewport={"width": 1800, "height": 1200}, device_scale_factor=a.scale)
            pg.goto(html.resolve().as_uri()); pg.wait_for_selector("g.clusterMapG", timeout=30000); pg.wait_for_timeout(1500)
            w, h = pg.evaluate(fit, 12)
            pg.set_viewport_size({"width": int(w) + 1, "height": int(h) + 1})
            pg.locator("#root-svg").screenshot(path=str(html.with_suffix(".png")))
            pg.close(); print(html.with_suffix(".png"), f"{int(w)}x{int(h)}", file=sys.stderr)
        b.close()


if __name__ == "__main__":
    main()
