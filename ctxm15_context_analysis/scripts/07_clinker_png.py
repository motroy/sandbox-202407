#!/usr/bin/env python3
"""Render the Clinker HTML plots (06_clinker_windows.py --run) to PNG/SVG with headless Chromium.

Usage: 07_clinker_png.py <clinker_dir> [--style [--grey-links]] [--scale 2] [--browser /path/to/chrome]
Plain mode writes <variant>.png (cropped to the figure, no control panel).
--style adjusts the plot first, then also writes <variant>.svg and <variant>.styled.html:
  * gene groups and their links are coloured by class (--grey-links keeps identity-shaded ribbons + identity bar) (pastel: blaCTX-M red, ISEcp1 blue, other IS purple, other ARGs orange, rest light grey)
  * legend entries are renamed from the member genes' products instead of "Group N"
  * gene labels only for blaCTX-M, ISEcp1, IS26 and ARGs (a custom label type "short"); white label halo, more row spacing
Needs `pip install playwright` and a Chromium (PLAYWRIGHT_BROWSERS_PATH is honoured; otherwise pass --browser).
"""
import argparse, json, re, sys
from collections import Counter
from pathlib import Path

from playwright.sync_api import sync_playwright

COL = {"blaCTX-M": "#f28e8e", "ISEcp1": "#8fbbe8", "IS": "#c9a9e3", "ARG": "#ffc58a", "other": "#e6e6e6"}   # pastel
ARG_RX = [(r"CTX-M-?(\d+)?", None), (r"beta-lactamase|OXA", "bla"), (r"aminoglycoside|acetyltransferase.*(AAC|aac)|nucleotidyltransferase.*aad", "aac/aad"),
          (r"dihydropteroate", "sul"), (r"dihydrofolate", "dfr"), (r"tetracycline", "tet"), (r"quinolone resistance|Qnr", "qnr"),
          (r"chloramphenicol", "cat"), (r"macrolide", "mph")]


def classify(product):
    """(class, short label) for a gene product string."""
    p = product or ""
    m = re.search(r"CTX-M-?(\d+)", p)
    if m or re.search(r"extended-spectrum.*beta-lactamase", p, re.I):
        return "blaCTX-M", f"blaCTX-M-{m.group(1)}" if m else "blaCTX-M"
    if re.search(r"ISEcp1", p): return "ISEcp1", "ISEcp1"
    m = re.findall(r"\b(IS\d+[A-Za-z]*|ISE[a-z]\w*)\b", p)       # "IS6-like element IS26 family" -> IS26
    if m or re.search(r"transposase|integrase|recombinase", p, re.I): return "IS", m[-1] if m else "tnp"
    for rx, name in ARG_RX[1:]:
        if re.search(rx, p, re.I): return "ARG", name
    return "other", " "                      # a space keeps Clinker from falling back to the locus tag


def style_html(text):
    """Return the Clinker HTML with group colours/labels and a 'short' gene label type injected."""
    i = text.index("const data=") + len("const data=")
    data, n = json.JSONDecoder().raw_decode(text[i:])
    gene = {}
    for c in data["clusters"]:
        for l in c["loci"]:
            for g in l["genes"]:
                cls, short = classify(g["names"].get("product", ""))
                if cls == "IS" and short != "IS26": short = " "          # label only the mediators: ISEcp1, IS26
                g["names"]["short"] = short; gene[g["uid"]] = (cls, g["names"].get("product", ""))
    for grp in data["groups"]:
        members = [gene[u] for u in grp["genes"] if u in gene]
        cls = Counter(c for c, _ in members).most_common(1)[0][0]
        prod = Counter(p for _, p in members if p).most_common(1)
        name = {"blaCTX-M": "blaCTX-M", "ISEcp1": "ISEcp1 transposase"}.get(cls) or (prod[0][0] if prod else grp["label"])
        if cls == "IS" and re.search(r"IS26", name): name = "IS26 (IS6 family) transposase"
        if len(name) > 38: name = name[:38].rsplit(" ", 1)[0] + "..."
        grp["label"] = f"{name} (n={len(members)})"; grp["colour"] = COL[cls]
    return text[:i] + json.dumps(data) + text[i + n:]


FIT = """(pad) => {                           // crop to the drawing, drop the control panel
    const svg = document.querySelector('#root-svg'), g = svg.querySelector('g.clusterMapG');
    const r = g.getBoundingClientRect(), o = svg.getBoundingClientRect();
    const x = r.x - o.x - pad, y = r.y - o.y - pad, w = r.width + 2 * pad, h = r.height + 2 * pad;
    svg.setAttribute('viewBox', `${x} ${y} ${w} ${h}`); svg.setAttribute('width', w); svg.setAttribute('height', h);
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg'); svg.style.background = 'white';
    document.body.replaceChildren(svg);
    return [w, h];
}"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("clinker_dir"); ap.add_argument("--scale", type=float, default=2)
    ap.add_argument("--style", action="store_true", help="colour groups, rename legend, short gene labels; also write SVG")
    ap.add_argument("--grey-links", action="store_true", help="with --style: keep identity-shaded grey ribbons (and the identity bar)")
    ap.add_argument("--browser", default=None, help="Chromium executable (default: Playwright's own)")
    a = ap.parse_args()
    with sync_playwright() as p:
        b = p.chromium.launch(executable_path=a.browser) if a.browser else p.chromium.launch()
        for html in sorted(f for f in Path(a.clinker_dir).glob("*.html") if not f.name.endswith(".styled.html")):
            src = html
            if a.style:
                src = html.with_name(html.stem + ".styled.html")
                src.write_text(style_html(html.read_text()))
            pg = b.new_page(viewport={"width": 1800, "height": 1200}, device_scale_factor=a.scale)
            pg.goto(src.resolve().as_uri()); pg.wait_for_selector("g.clusterMapG", timeout=30000); pg.wait_for_timeout(1000)
            if a.style:
                pg.check("#input-gene-labels"); pg.select_option("#select-label-type", "short")
                if not a.grey_links: pg.check("#input-link-group-colour")     # ribbons take their group's pastel colour
                for sel, val in (("#input-label-rotation", "45"), ("#input-cluster-spacing", "58"), ("#input-label-size", "10")):
                    pg.fill(sel, val); pg.dispatch_event(sel, "change")
                # white halo so labels stay readable over the black identity ribbons
                pg.add_style_tag(content="text.geneLabel{paint-order:stroke;stroke:#fff;stroke-width:3.5px;stroke-linejoin:round}")
                pg.wait_for_timeout(800)
            w, h = pg.evaluate(FIT, 14)
            if a.style: html.with_suffix(".svg").write_text(pg.evaluate("document.querySelector('#root-svg').outerHTML"))
            pg.set_viewport_size({"width": int(w) + 1, "height": int(h) + 1})
            pg.locator("#root-svg").screenshot(path=str(html.with_suffix(".png")))
            pg.close(); print(html.with_suffix(".png"), f"{int(w)}x{int(h)}", file=sys.stderr)
        b.close()


if __name__ == "__main__":
    main()
