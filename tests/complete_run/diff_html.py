"""Render a browsable HTML index of complete-run image differences.

The page is self-contained and written beside the JSON report. It does not use
``output_viewer``, which renders a fixed table with no severity ordering,
filtering, or incremental loading.
"""

from __future__ import annotations

import html
import json
import os
import re
from pathlib import Path

_FRACTION = re.compile(r"fraction: ([0-9.eE+-]+)")

_STYLE = """
:root{--bg:#f7f7f8;--panel:#fff;--ink:#1b1b1f;--muted:#5f6169;--line:#e2e3e8;
--accent:#8b2f2f;--ok:#1f7a4d;--chip:#eef0f4}
@media (prefers-color-scheme:dark){:root{--bg:#15161a;--panel:#1d1f24;--ink:#e9e9ee;
--muted:#a0a3ad;--line:#2c2f36;--accent:#e08a8a;--ok:#6fd39b;--chip:#262931}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
header{padding:24px 28px 18px;border-bottom:1px solid var(--line);background:var(--panel)}
h1{margin:0 0 4px;font-size:20px;letter-spacing:-.01em}
.sub{color:var(--muted);font-size:13px}
.stats{display:flex;gap:28px;flex-wrap:wrap;margin:16px 0 4px}
.stat b{display:block;font-size:22px;font-variant-numeric:tabular-nums}
.stat span{color:var(--muted);font-size:12px;text-transform:uppercase;letter-spacing:.04em}
.stat.ok b{color:var(--ok)}.stat.warn b{color:var(--accent)}
details.envbox{margin-top:14px;font-size:13px}
details.envbox ul{margin:8px 0 0 18px;padding:0}
.controls{position:sticky;top:0;z-index:5;background:var(--panel);
border-bottom:1px solid var(--line);padding:12px 28px;display:flex;gap:10px;
flex-wrap:wrap;align-items:center}
input[type=search]{padding:7px 10px;border:1px solid var(--line);border-radius:7px;
background:var(--bg);color:var(--ink);min-width:230px;font-size:13px}
.chip{border:1px solid var(--line);background:var(--chip);color:var(--ink);
border-radius:999px;padding:5px 11px;font-size:12px;cursor:pointer}
.chip[aria-pressed=true]{background:var(--ink);color:var(--panel);border-color:var(--ink)}
main{padding:20px 28px 60px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;
margin-bottom:14px;overflow:hidden}
.head{padding:10px 14px;display:flex;gap:14px;align-items:center;
border-bottom:1px solid var(--line)}
.frac{font-variant-numeric:tabular-nums;font-weight:650;color:var(--accent);min-width:74px}
.name{flex:1;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;
word-break:break-all}
.tag{background:var(--chip);border-radius:5px;padding:2px 7px;font-size:11px;
color:var(--muted);white-space:nowrap}
.trip{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;padding:12px 14px}
@media (max-width:900px){.trip{grid-template-columns:1fr}}
figure{margin:0}
figcaption{font-size:11px;color:var(--muted);text-transform:uppercase;
letter-spacing:.05em;margin-bottom:5px}
.trip a{display:block}
.trip img{width:100%;height:300px;object-fit:contain;border:1px solid var(--line);
border-radius:6px;background:#fff;display:block}
body.fit .trip img{height:auto}
.empty{color:var(--muted);padding:30px 0}
section.extra{margin-top:34px}
section.extra h2{font-size:15px;margin:0 0 6px}
section.extra p{color:var(--muted);margin:0 0 10px}
section.extra ul{columns:2;font-family:ui-monospace,monospace;font-size:12px;color:var(--muted)}
@media (max-width:900px){section.extra ul{columns:1}}
#sentinel{height:1px}
"""

_SCRIPT = """
const list=document.getElementById('list'),empty=document.getElementById('empty'),
shown=document.getElementById('shown'),sentinel=document.getElementById('sentinel');
const PAGE=25;let activeSet='*',query='',items=[],drawn=0;
function card(r){return `<div class="card"><div class="head">
<span class="frac">${(r.frac*100).toFixed(2)}%</span>
<span class="name">${r.path}</span><span class="tag">${r.set}</span></div>
<div class="trip">
<figure><figcaption>Baseline</figcaption><a href="${r.expected}" target="_blank" rel="noopener">
<img loading="lazy" src="${r.expected}" alt="baseline"></a></figure>
<figure><figcaption>This run</figcaption><a href="${r.actual}" target="_blank" rel="noopener">
<img loading="lazy" src="${r.actual}" alt="this run"></a></figure>
<figure><figcaption>Diff</figcaption><a href="${r.diff}" target="_blank" rel="noopener">
<img loading="lazy" src="${r.diff}" alt="diff"></a></figure></div></div>`;}
function more(){if(drawn>=items.length)return;
const next=items.slice(drawn,drawn+PAGE);
list.insertAdjacentHTML('beforeend',next.map(card).join(''));drawn+=next.length;
shown.textContent=`${drawn} of ${items.length} shown`;}
function render(){items=ROWS.filter(r=>(activeSet==='*'||r.set===activeSet)&&
(query===''||r.path.toLowerCase().includes(query)));
list.innerHTML='';drawn=0;empty.hidden=items.length>0;more();}
new IntersectionObserver(e=>{if(e[0].isIntersecting)more();},
{rootMargin:'900px'}).observe(sentinel);
document.getElementById('q').addEventListener('input',e=>{
query=e.target.value.trim().toLowerCase();render();});
document.querySelectorAll('.chip[data-set]').forEach(b=>b.addEventListener('click',()=>{
document.querySelectorAll('.chip[data-set]').forEach(o=>o.setAttribute('aria-pressed','false'));
b.setAttribute('aria-pressed','true');activeSet=b.dataset.set;render();window.scrollTo(0,0);}));
document.getElementById('size').addEventListener('click',e=>{
const on=document.body.classList.toggle('fit');
e.target.setAttribute('aria-pressed',String(on));
e.target.textContent=on?'shorter':'taller';});
render();
"""


def _mismatch_fraction(detail: str) -> float:
    """Return the mismatched-pixel fraction recorded in an issue detail."""
    match = _FRACTION.search(detail or "")

    return float(match.group(1)) if match else 0.0


def _build_rows(report: dict, root: Path) -> list[dict[str, object]]:
    """Return image-mismatch rows sorted with the largest difference first."""
    entries = sorted(
        (e for e in report["summary"]["image_mismatches"] if e.get("artifact_path")),
        key=lambda e: _mismatch_fraction(e.get("detail", "")),
        reverse=True,
    )
    rows: list[dict[str, object]] = []
    for entry in entries:
        diff = os.path.relpath(entry["artifact_path"], root)
        rows.append(
            {
                "path": entry["relative_path"],
                "set": entry["relative_path"].split("/")[0],
                "frac": _mismatch_fraction(entry.get("detail", "")),
                "diff": diff,
                "actual": diff.replace("_diff.png", "_actual.png"),
                "expected": diff.replace("_diff.png", "_expected.png"),
            }
        )

    return rows


def write_diff_html(report: dict, report_path: str | Path) -> Path | None:
    """Write an HTML index of image differences beside the JSON report.

    Parameters
    ----------
    report : dict
        The comparison report, as written to ``comparison-report.json``.
    report_path : str | Path
        Path of that JSON report; the page is written to its directory.

    Returns
    -------
    Path | None
        Path of the written page, or ``None`` when there are no image
        mismatches with diff artifacts to link.
    """
    root = Path(report_path).parent
    rows = _build_rows(report, root)
    if not rows:
        return None

    summary = report["summary"]
    sets = sorted({str(row["set"]) for row in rows})
    counts = {name: sum(1 for row in rows if row["set"] == name) for name in sets}
    environment = report.get("environment") or {}
    differences = environment.get("differences") or []
    missing_files = summary.get("missing_baseline_files", [])
    missing_images = summary.get("missing_baseline_images", [])
    extra_sets = sorted({path.split("/")[0] for path in missing_files + missing_images})

    chips = "".join(
        f'<button class="chip" data-set="{html.escape(name)}" aria-pressed="false">'
        f"{html.escape(name)} ({counts[name]})</button>"
        for name in sets
    )
    env_items = (
        "".join(f"<li><code>{html.escape(str(d))}</code></li>" for d in differences)
        or "<li>none</li>"
    )
    extra_items = "".join(
        f"<li>{html.escape(path)}</li>" for path in sorted(missing_images)
    )

    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Complete-run image diffs</title>
<style>{_STYLE}</style>
</head>
<body>
<header>
  <h1>Complete-run image diffs</h1>
  <div class="sub"><code>{html.escape(Path(report["paths"]["dev_dir"]).name)}</code>
    vs baseline <code>{html.escape(Path(report["paths"]["baseline_dir"]).name)}</code></div>
  <div class="stats">
    <div class="stat ok"><b>{len(summary["matching_files"])} / {summary["compared_file_count"]}</b>
      <span>netCDF files match</span></div>
    <div class="stat warn"><b>{len(rows)}</b><span>image mismatches</span></div>
    <div class="stat"><b>{len(summary["matching_images"])}</b><span>images matched</span></div>
    <div class="stat"><b>{summary["failure_count"]}</b><span>total findings</span></div>
  </div>
  <details class="envbox">
    <summary>Environment differences vs baseline</summary>
    <ul>{env_items}</ul>
  </details>
</header>
<div class="controls">
  <input type="search" id="q" placeholder="Filter by filename&hellip;" aria-label="Filter by filename">
  <button class="chip" data-set="*" aria-pressed="true">all ({len(rows)})</button>
  {chips}
  <button class="chip" id="size" aria-pressed="false">taller</button>
  <span class="tag" id="shown"></span>
</div>
<main>
  <div id="list"></div>
  <div id="sentinel"></div>
  <div class="empty" id="empty" hidden>Nothing matches that filter.</div>
  <section class="extra">
    <h2>Present in this run, absent from the baseline</h2>
    <p>{len(missing_files)} files and {len(missing_images)} images, under
      <code>{html.escape(", ".join(extra_sets)) or "&mdash;"}</code>. These have no diff
      images, since the baseline has nothing to compare against.</p>
    <ul>{extra_items}</ul>
  </section>
</main>
<script>
const ROWS = {json.dumps(rows)};
{_SCRIPT}
</script>
</body>
</html>
"""

    output_path = root / "index.html"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(page, encoding="utf-8")

    return output_path
