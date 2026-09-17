"""Render a browsable HTML index of complete-run image differences.

The page is self-contained and written beside the JSON report. It does not use
``output_viewer``, which renders a fixed table with no severity ordering,
filtering, or incremental loading.
"""

from __future__ import annotations

import html
import json
import os
from pathlib import Path
from typing import Sequence

from tests.complete_run.image_severity import REVIEWABLE_SEVERITIES, SEVERITY_ORDER

_SEVERITY_RANK: dict[str, int] = {
    severity: index for index, severity in enumerate(SEVERITY_ORDER)
}
_SEVERITY_DESCRIPTIONS = {
    "STRUCTURAL": "figure geometry changed substantially",
    "MAJOR": "a large part of the plot changed",
    "MODERATE": "a visible part of the plot changed",
    "MINOR": "a small or isolated change needs review",
    "NEGLIGIBLE": "cosmetic difference; counted as passed",
    "IDENTICAL": "exact image match; counted as passed",
}

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
.severity-guide{margin-top:16px;padding:12px;border:1px solid var(--line);border-radius:8px;
 background:var(--bg)}
.severity-guide h2{margin:0 0 8px;font-size:13px}.severity-guide h2 span{font-weight:400;color:var(--muted)}
.severity-table{width:100%;border-collapse:collapse;font-size:12px}.severity-table th{text-align:left;
 color:var(--muted);font-size:10px;letter-spacing:.05em;text-transform:uppercase}.severity-table th,
.severity-table td{padding:5px 8px;border-top:1px solid var(--line)}.severity-table th:first-child,
.severity-table td:first-child{padding-left:0}.severity-table th:last-child,.severity-table td:last-child{padding-right:0}
.severity-count{font-variant-numeric:tabular-nums;font-weight:700}.severity-status{color:var(--muted)}
.controls{position:sticky;top:0;z-index:5;background:var(--panel);
border-bottom:1px solid var(--line);padding:12px 28px;display:flex;gap:10px;
flex-wrap:wrap;align-items:center}
.control-group{display:flex;gap:7px;align-items:center;flex-wrap:wrap}
.control-label{color:var(--muted);font-size:11px;font-weight:700;letter-spacing:.05em;
text-transform:uppercase;white-space:nowrap}
.control-divider{height:26px;border-left:1px solid var(--line)}
input[type=search]{padding:7px 10px;border:1px solid var(--line);border-radius:7px;
background:var(--bg);color:var(--ink);min-width:230px;font-size:13px}
.chip{border:1px solid var(--line);background:var(--chip);color:var(--ink);
border-radius:999px;padding:5px 11px;font-size:12px;cursor:pointer}
.chip[aria-pressed=true]{background:var(--ink);color:var(--panel);border-color:var(--ink)}
.chip:disabled{cursor:default;opacity:.5}
main{padding:20px 28px 60px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;
margin-bottom:14px;overflow:hidden}
.head{padding:10px 14px;display:flex;gap:14px;align-items:center;
border-bottom:1px solid var(--line)}
.frac{font-variant-numeric:tabular-nums;font-weight:650;color:var(--accent);min-width:118px}
.score-label{display:block;color:var(--muted);font-size:10px;font-weight:400;text-transform:uppercase;
letter-spacing:.04em}
.name{flex:1;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;
word-break:break-all}
.tag{background:var(--chip);border-radius:5px;padding:2px 7px;font-size:11px;
color:var(--muted);white-space:nowrap}
.severity{font-weight:700}.severity-STRUCTURAL{color:#8b2f2f}.severity-MAJOR{color:#ad5a12}
.severity-MODERATE{color:#886d08}.severity-MINOR{color:#336d9a}.severity-NEGLIGIBLE{color:var(--muted)}
.severity-IDENTICAL{color:var(--ok)}
.cause{color:var(--muted);font-size:12px}
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
const PAGE=25;let activeSet='*',activeSeverity='REVIEW',query='',items=[],drawn=0;
function card(r){const cosmetic=r.severity==='NEGLIGIBLE';return `<div class="card"><div class="head">
  <span class="frac">${((cosmetic?r.raw_frac:r.frac)*100).toFixed(2)}%<span class="score-label">${cosmetic?'pixels differ':'unmatched content'}</span></span>
 <span class="name">${r.path}</span><span class="tag severity severity-${r.severity}">${r.level}. ${r.severity}</span>
<span class="tag">${r.set}</span><span class="cause">${r.cause}</span></div>
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
 (activeSeverity==='*'||(activeSeverity==='REVIEW'?r.reviewable:r.severity===activeSeverity))&&
(query===''||r.path.toLowerCase().includes(query)));
list.innerHTML='';drawn=0;empty.hidden=items.length>0;more();}
new IntersectionObserver(e=>{if(e[0].isIntersecting)more();},
{rootMargin:'900px'}).observe(sentinel);
document.getElementById('q').addEventListener('input',e=>{
query=e.target.value.trim().toLowerCase();render();});
document.querySelectorAll('.chip[data-set]').forEach(b=>b.addEventListener('click',()=>{
document.querySelectorAll('.chip[data-set]').forEach(o=>o.setAttribute('aria-pressed','false'));
b.setAttribute('aria-pressed','true');activeSet=b.dataset.set;render();window.scrollTo(0,0);}));
document.querySelectorAll('.chip[data-severity]').forEach(b=>b.addEventListener('click',()=>{
document.querySelectorAll('.chip[data-severity]').forEach(o=>o.setAttribute('aria-pressed','false'));
b.setAttribute('aria-pressed','true');activeSeverity=b.dataset.severity;render();window.scrollTo(0,0);}));
document.getElementById('size').addEventListener('click',e=>{
const on=document.body.classList.toggle('fit');
e.target.setAttribute('aria-pressed',String(on));
e.target.textContent=on?'shorter':'taller';});
render();
"""


def _build_rows(report: dict, root: Path) -> list[dict[str, object]]:
    """Return review rows and bounded cosmetic samples, worst-first."""
    summary = report["summary"]
    entries = sorted(
        (
            entry
            for entry in [
                *summary["image_mismatches"],
                *summary.get("cosmetic_samples", []),
            ]
            if entry.get("artifact_path")
        ),
        key=lambda entry: (
            _SEVERITY_RANK.get(str(entry.get("severity")), -1),
            float(entry.get("content_fraction") or 0.0),
        ),
        reverse=True,
    )
    rows: list[dict[str, object]] = []
    for entry in entries:
        diff = os.path.relpath(entry["artifact_path"], root)
        rows.append(
            {
                "path": entry["relative_path"],
                "set": entry["relative_path"].split("/")[0],
                "severity": entry.get("severity", "UNKNOWN"),
                "level": _SEVERITY_RANK.get(str(entry.get("severity")), -1) + 1,
                "cause": entry.get("cause") or "cause unavailable",
                "frac": float(entry.get("content_fraction") or 0.0),
                "raw_frac": float(entry.get("raw_fraction") or 0.0),
                "reviewable": entry.get("severity") in REVIEWABLE_SEVERITIES,
                "diff": diff,
                "actual": diff.replace("_diff.png", "_actual.png"),
                "expected": diff.replace("_diff.png", "_expected.png"),
            }
        )

    return rows


def _severity_counts(
    summary: dict, reviewable_severities: Sequence[str]
) -> dict[str, int]:
    """Return total image counts for every severity, including passing images."""
    counts: dict[str, int] = {severity: 0 for severity in SEVERITY_ORDER}
    counts["IDENTICAL"] = len(summary.get("identical_images", []))
    counts["NEGLIGIBLE"] = len(summary.get("cosmetic_images", []))
    for entry in summary["image_mismatches"]:
        severity = entry.get("severity")
        if severity in reviewable_severities:
            counts[severity] += 1

    return counts


def _severity_status(severity: str, cosmetic_sample_count: int) -> str:
    """Describe whether a severity is reviewable or represented by a sample."""
    if severity in REVIEWABLE_SEVERITIES:
        return "Needs review"
    if severity == "NEGLIGIBLE" and cosmetic_sample_count:
        return "Passed; sample available"

    return "Passed"


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
    all_severities = list(SEVERITY_ORDER)
    reviewable_severities = list(reversed(REVIEWABLE_SEVERITIES))
    severity_counts = _severity_counts(summary, reviewable_severities)
    cosmetic_sample_count = len(summary.get("cosmetic_samples", []))
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
    severity_chips = "".join(
        f'<button class="chip" data-severity="{severity}" aria-pressed="false"'
        f"{' disabled' if severity_counts[severity] == 0 else ''}>"
        f"{_SEVERITY_RANK[severity] + 1}. {severity.lower()} ({severity_counts[severity]})</button>"
        for severity in reviewable_severities
    )
    negligible_chip = (
        f'<button class="chip" data-severity="NEGLIGIBLE" aria-pressed="false"'
        f"{' disabled' if cosmetic_sample_count == 0 else ''}>"
        f"{_SEVERITY_RANK['NEGLIGIBLE'] + 1}. negligible "
        f"({severity_counts['NEGLIGIBLE']} total, {cosmetic_sample_count} sampled)</button>"
    )
    severity_table = "".join(
        "<tr>"
        f'<td class="severity severity-{severity}">{index}. {severity.title()}</td>'
        f"<td>{_SEVERITY_DESCRIPTIONS[severity]}</td>"
        f'<td class="severity-count">{severity_counts[severity]}</td>'
        f'<td class="severity-status">{_severity_status(severity, cosmetic_sample_count)}</td>'
        "</tr>"
        for index, severity in enumerate(all_severities, start=1)
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
      <span>netCDF files passing</span></div>
    <div class="stat warn"><b>{len(summary["image_mismatches"])}</b><span>images needing review</span></div>
    <div class="stat"><b>{len(summary["matching_images"])}</b>
      <span>images passing ({len(summary.get("identical_images", []))} identical, {len(summary.get("cosmetic_images", []))} cosmetic)</span></div>
    <div class="stat"><b>{summary["failure_count"] - len(summary["image_mismatches"])}</b><span>other comparison findings</span></div>
  </div>
   <section class="severity-guide" aria-label="Severity guide">
     <h2>Severity guide <span>low to high</span></h2>
     <table class="severity-table">
       <thead><tr><th>Severity</th><th>Definition</th><th>Images</th><th>Status</th></tr></thead>
       <tbody>{severity_table}</tbody>
     </table>
  </section>
  <details class="envbox">
    <summary>Environment differences vs baseline</summary>
    <ul>{env_items}</ul>
  </details>
</header>
<div class="controls">
  <input type="search" id="q" placeholder="Filter by filename&hellip;" aria-label="Filter by filename">
  <div class="control-group" aria-label="Filter by diagnostic set">
    <span class="control-label">Diagnostic set</span>
     <button class="chip" data-set="*" aria-pressed="true">available ({len(rows)})</button>
    {chips}
  </div>
  <span class="control-divider" aria-hidden="true"></span>
  <div class="control-group" aria-label="Filter by severity">
     <span class="control-label">Severity, highest first</span>
     <button class="chip" data-severity="REVIEW" aria-pressed="true">review levels ({len(summary["image_mismatches"])})</button>
     <button class="chip" data-severity="*" aria-pressed="false">all available ({len(rows)})</button>
     {severity_chips}
     {negligible_chip}
  </div>
  <span class="control-divider" aria-hidden="true"></span>
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
