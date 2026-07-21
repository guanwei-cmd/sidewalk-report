#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unified final pipeline for SIDEWALK_202412 (dataset 58791).

One consistent method:
  - Forensics on ALL segments with SW_LENG > 5000 m (geometry bbox diagonal,
    area/width equivalent length). Attribute error iff SW_LENG > 5 x bbox_diag.
  - Drop attribute-error segments; keep geometrically plausible ones.
  - Recompute per-county + national stats on kept rows:
      w bands  : SW_WTH length-weighted, denom = SW_WTH valid (>0)
      net15    : share of length with SWW_WTH < 1.5, denom = SWW_WTH valid (>0)
      CDF      : 51 values t=0.0..5.0 step 0.1, share of SWW_WTH >= t, same denom
  - Cross-check net15 == 100 - CDF[15] for every county.
  - Old-vs-new net15: old = pre-exclusion, denom = SWW_WTH not None (legacy
    analyze.py/national.py method); new = post-exclusion, denom SWW_WTH > 0.
"""
import json, math
from collections import defaultdict

BASE = '/private/tmp/claude-501/-Users-dingguanwei/1dfd1787-952e-49d3-8cba-fa3721893f8d/scratchpad/swdata'
FILES = [('p1', f'{BASE}/p1/SIDEWALK_1_202412_WGS84.geojson'),
         ('p2', f'{BASE}/p2/SIDEWALK_2_202412_WGS84.geojson')]

LATM = 110574.0
def lonm(lat): return 111320.0 * math.cos(math.radians(lat))

def geom_stats(geom):
    """bbox diagonal (m), polygon area (m^2), total ring perimeter (m)."""
    xs, ys = [], []
    def collect(c):
        if not c: return
        if isinstance(c[0], (int, float)):
            xs.append(c[0]); ys.append(c[1])
        else:
            for e in c: collect(e)
    collect(geom['coordinates'])
    if not xs: return None, None, None
    y0m, y1m = min(ys), max(ys)
    latc = (y0m + y1m) / 2.0
    lm = lonm(latc)
    diag = math.hypot((max(xs) - min(xs)) * lm, (y1m - y0m) * LATM)
    polys = geom['coordinates'] if geom['type'] == 'MultiPolygon' else [geom['coordinates']]
    area = 0.0; perim = 0.0
    for poly in polys:
        for ri, ring in enumerate(poly):
            pts = [(x * lm, y * LATM) for x, y in ring]
            a = 0.0; pp = 0.0
            for i in range(len(pts) - 1):
                x1, y1 = pts[i]; x2, y2 = pts[i + 1]
                a += x1 * y2 - x2 * y1
                pp += math.hypot(x2 - x1, y2 - y1)
            a = abs(a) / 2.0
            area += a if ri == 0 else -a   # exterior minus holes
            perim += pp
    return diag, max(area, 0.0), perim

# ---------- pass over both files ----------
rows = []            # (cty, leng, wth, sww)
big = []             # forensic candidates: dict incl. row index
file_counts = {}
parse_fail = 0
for tag, path in FILES:
    cnt = 0
    with open(path, encoding='utf-8-sig') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('{ "type": "Feature"'):
                continue
            if line.endswith(','):
                line = line[:-1]
            try:
                feat = json.loads(line)
            except json.JSONDecodeError:
                parse_fail += 1
                continue
            p = feat['properties']
            cty = p.get('COUNTY_NA') or ''
            leng = p.get('SW_LENG'); wth = p.get('SW_WTH'); sww = p.get('SWW_WTH')
            idx = len(rows)
            rows.append((cty, leng, wth, sww))
            cnt += 1
            if (leng or 0) > 5000:
                diag, area, perim = geom_stats(feat['geometry']) if feat.get('geometry') else (None, None, None)
                w_eq = wth if (wth or 0) > 0 else 3.0
                eq_len = (area / w_eq) if area else None
                big.append({'idx': idx, 'file': tag, 'cty': cty,
                            'NAME': p.get('NAME'), 'PSTART': p.get('PSTART'), 'PEND': p.get('PEND'),
                            'VILL': p.get('VILL_NAME'),
                            'SW_LENG': leng, 'SW_WTH': wth, 'SWW_WTH': sww,
                            'diag': diag, 'area': area, 'perim': perim,
                            'eq_len': eq_len, 'eq_w_used': w_eq})
    file_counts[tag] = cnt

print(f'FILE COUNTS: p1={file_counts["p1"]}  p2={file_counts["p2"]}  '
      f'sum={file_counts["p1"]+file_counts["p2"]}  parse_fail={parse_fail}')
print()

# ---------- forensics ----------
print('=== FORENSICS: all SW_LENG > 5000 m ===')
print('cty | NAME | SW_LENG | bbox_diag | eq_len(area/w) | perim/2 | leng/diag | leng/eq | leng/(perim/2) | verdict')
drop_idx = set()
def fmt(v): return f'{v:.1f}' if v is not None else 'NA'
for r in sorted(big, key=lambda r: -(r['SW_LENG'] or 0)):
    ratio_d = (r['SW_LENG'] / r['diag']) if r['diag'] else float('inf')
    ratio_e = (r['SW_LENG'] / r['eq_len']) if r['eq_len'] else float('inf')
    ratio_p = (r['SW_LENG'] / (r['perim'] / 2)) if r['perim'] else float('inf')
    # rule: attribute error iff SW_LENG > 5 x bbox diagonal,
    # EXCEPT geometry corroborates the length: for a ribbon polygon the
    # centerline length ~ perimeter/2, so keep if SW_LENG within 2x of perim/2
    # (serpentine e.g. mountain path). Empty geometry -> unverifiable -> drop.
    geom_ok = (r['perim'] or 0) > 0 and (0.5 <= ratio_p <= 2.0)
    bad = (ratio_d > 5.0) and not geom_ok
    if bad: drop_idx.add(r['idx'])
    r['ratio_d'] = ratio_d; r['ratio_e'] = ratio_e; r['ratio_p'] = ratio_p
    r['bad'] = bad; r['geom_ok'] = geom_ok
    if bad:
        verdict = 'ATTR-ERROR -> DROP' if r['diag'] is not None else 'NO-GEOMETRY -> DROP'
    elif ratio_d > 5.0:
        verdict = 'KEEP (over 5x bbox but perimeter-consistent serpentine)'
    else:
        verdict = 'PLAUSIBLE -> KEEP'
    print(f"{r['cty']} | {r['NAME']} ({r['PSTART']}->{r['PEND']}, {r['VILL']}) | "
          f"{r['SW_LENG']:.1f} | {fmt(r['diag'])} | {fmt(r['eq_len'])} (w={r['eq_w_used']}) | "
          f"{fmt(r['perim']/2 if r['perim'] else None)} | "
          f"{ratio_d:.1f}x | {ratio_e:.1f}x | {ratio_p:.1f}x | {verdict}")
print(f'\ncandidates={len(big)}  dropped={len(drop_idx)}  kept_over5000={len(big)-len(drop_idx)}')
dropped_by_cty = defaultdict(list)
for r in big:
    if r['bad']: dropped_by_cty[r['cty']].append(r)
affected = sorted(dropped_by_cty)
print('affected counties:', ', '.join(f'{c}(-{len(v)})' for c, v in sorted(dropped_by_cty.items())))
print()

# ---------- stats machinery ----------
def make_stats(rs):
    n = len(rs); tot = 0.0
    wb = [0.0, 0.0, 0.0]; wden = 0.0
    nden = 0.0; nlt = 0.0
    swws = []   # (sww, L) for CDF
    for _, l, w, s in rs:
        L = l or 0.0
        tot += L
        if w is not None and w > 0:
            wden += L
            if w < 1.5: wb[0] += L
            elif w < 2.5: wb[1] += L
            else: wb[2] += L
        if s is not None and s > 0:
            nden += L; swws.append((s, L))
            if s < 1.5: nlt += L
    cdf = []
    for i in range(51):
        t = i / 10.0
        num = sum(L for s, L in swws if s >= t - 1e-9)
        cdf.append(round(100.0 * num / nden, 1) if nden else None)
    return {'n': n, 'tot_km': tot / 1000.0,
            'w': [round(100.0 * b / wden, 1) for b in wb] if wden else None,
            'w_den_km': wden / 1000.0,
            'net15': round(100.0 * nlt / nden, 1) if nden else None,
            'net_den_km': nden / 1000.0, 'cdf': cdf}

def legacy_net15(rs):
    """old method: denom = SWW_WTH not None (0 included), no exclusion"""
    den = 0.0; lt = 0.0
    for _, l, w, s in rs:
        if s is None: continue
        L = l or 0.0
        den += L
        if s < 1.5: lt += L
    return round(100.0 * lt / den, 1) if den else None

by_cty_all = defaultdict(list)
by_cty_kept = defaultdict(list)
kept = []
for i, r in enumerate(rows):
    by_cty_all[r[0]].append(r)
    if i not in drop_idx:
        by_cty_kept[r[0]].append(r)
        kept.append(r)

CTYS = sorted(by_cty_all, key=lambda c: -sum((r[1] or 0) for r in by_cty_kept[c]))
print(f'counties present: {len(CTYS)}')

nat_all = make_stats(rows)
nat_kept = make_stats(kept)
stats_kept = {c: make_stats(by_cty_kept[c]) for c in CTYS}
old15 = {c: legacy_net15(by_cty_all[c]) for c in CTYS}
old15_nat = legacy_net15(rows)

# also count SW_WTH / SWW_WTH zero & None prevalence for method note
zw = sum(1 for _, l, w, s in rows if w == 0); nw = sum(1 for _, l, w, s in rows if w is None)
zs = sum(1 for _, l, w, s in rows if s == 0); ns = sum(1 for _, l, w, s in rows if s is None)
print(f'SW_WTH: None={nw} zero={zw} | SWW_WTH: None={ns} zero={zs}  (of {len(rows)})')
print()

# ---------- 3a WIDTH rows (affected counties) ----------
print('=== 3a WIDTH rows (affected counties, after exclusion) ===')
for c in affected:
    s = stats_kept[c]
    print(json.dumps({'city': c, 'total_km': round(s['tot_km'], 1), 'n': s['n'],
                      'w': s['w'], 'net15': s['net15']}, ensure_ascii=False))
print()

# ---------- 3b net15 all 22 ----------
print('=== 3b net15 all counties: old(legacy) -> new(converged) ===')
for c in CTYS:
    print(f'{c}\told={old15[c]}\tnew={stats_kept[c]["net15"]}')
print(f'全國\told={old15_nat}\tnew={nat_kept["net15"]}')
print()

# ---------- 3c CDF ----------
print('=== 3c CDF (affected counties + national, after exclusion) ===')
for c in affected:
    print(f'"{c}": {json.dumps(stats_kept[c]["cdf"])},')
print(f'"全國": {json.dumps(nat_kept["cdf"])}')
print()

# ---------- 3d national ----------
print('=== 3d national ===')
print(f'original n={nat_all["n"]}  dropped={len(drop_idx)}  kept n={nat_kept["n"]}')
print(f'original total={nat_all["tot_km"]:.1f} km -> kept total={nat_kept["tot_km"]:.1f} km')
print(f'kept w={nat_kept["w"]} (denom {nat_kept["w_den_km"]:.1f} km)  '
      f'net15={nat_kept["net15"]} (denom {nat_kept["net_den_km"]:.1f} km)')
print(f'[ref] pre-exclusion w={nat_all["w"]}  net15={nat_all["net15"]}')
print()

# ---------- 4 cross-check ----------
print('=== 4 cross-check net15 vs 100-CDF[15] (all counties + national) ===')
bad_ck = 0
for c in CTYS + ['全國']:
    s = nat_kept if c == '全國' else stats_kept[c]
    if s['net15'] is None: continue
    d = abs(s['net15'] - round(100.0 - s['cdf'][15], 1))
    flag = 'OK' if d <= 0.1 + 1e-9 else 'MISMATCH'
    if flag != 'OK': bad_ck += 1
    print(f'{c}: net15={s["net15"]}  100-CDF[15]={round(100.0 - s["cdf"][15], 1)}  diff={d:.2f}  {flag}')
print(f'mismatches: {bad_ck}')
print()

# ---------- per-affected-county drop detail ----------
print('=== dropped segments detail per county ===')
for c in affected:
    tot_drop = sum(r['SW_LENG'] for r in dropped_by_cty[c])
    print(f'{c}: n_dropped={len(dropped_by_cty[c])}, dropped_len={tot_drop/1000:.2f} km')
    for r in dropped_by_cty[c]:
        print(f"   {r['NAME']} {r['SW_LENG']:.0f}m (diag {fmt(r['diag'])}m, {r['ratio_d']:.0f}x, perim/2 {fmt(r['perim']/2 if r['perim'] else None)}m)")

with open(f'{BASE}/final_out.json', 'w', encoding='utf-8') as f:
    json.dump({'file_counts': file_counts, 'big': [{k: v for k, v in r.items() if k != 'idx'} for r in big],
               'affected': affected,
               'stats_kept': stats_kept, 'nat_kept': nat_kept, 'nat_all': {k: v for k, v in nat_all.items()},
               'old15': old15, 'old15_nat': old15_nat}, f, ensure_ascii=False)
print('\nsaved final_out.json')
