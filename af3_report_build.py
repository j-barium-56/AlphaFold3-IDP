#!/usr/bin/env python3
"""
Build the AF3 IDR-IDR stage-1 results report (A4, vector) from af3_report_data.json.

Colour is used by job, from the documented reference palette:
  * three categorical slots for identity (validated all-pairs on the light surface);
    a fourth hue would fail the normal-vision floor, so the memorised controls are
    encoded by marker shape (hollow ring), not colour
  * a blue<->red diverging ramp with a neutral grey midpoint for AUC / delta magnitude
    (monotone lightness per arm, asserted at build time)
  * reserved status tokens for PASS / FAIL, always with the word, never colour alone
Every heat-map cell prints its own number and every plotted value also appears in an
appendix table - the relief required for the sub-3:1 hues. Marks are thin, rules are
hairlines, dots carry a 2pt surface ring, labels are measured before they are drawn.
"""
import json, os, math, csv, statistics, sys
sys.dont_write_bytecode = True

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4
from reportlab.lib.colors import HexColor
from reportlab.pdfgen import canvas

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'af3_report_data.json')))
OUT = os.path.join(HERE, 'af3_idr_results_report.pdf')
CSV_OUT = os.path.join(HERE, 'af3_all_statistics.csv')

FD = '/System/Library/Fonts/Supplemental'
for tag, fn in [('S', 'Arial.ttf'), ('SB', 'Arial Bold.ttf'), ('SI', 'Arial Italic.ttf'),
                ('N', 'Arial Narrow.ttf'), ('NB', 'Arial Narrow Bold.ttf')]:
    pdfmetrics.registerFont(TTFont(tag, f'{FD}/{fn}'))
pdfmetrics.registerFontFamily('S', normal='S', bold='SB', italic='SI')

# ---------------------------------------------------------------- palette (documented)
SURFACE = HexColor('#fcfcfb'); INK = HexColor('#0b0b0b'); INK2 = HexColor('#52514e')
MUTED = HexColor('#898781'); GRID = HexColor('#e1e0d9'); AXIS = HexColor('#c3c2b7')
ACCENT = HexColor('#4a3aa7')                       # slot 7 violet - chrome only
POS = HexColor('#2a78d6'); NEG = HexColor('#eb6834'); SCR = HexColor('#1baf7a')
GOOD = HexColor('#0ca30c'); CRIT = HexColor('#d03b3b'); WARNC = HexColor('#fab219')
MIDG = HexColor('#f0efec'); WASH = HexColor('#f4f3f8')
BLUE_ARM = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95']
RED_ARM = ['#fbdcdc', '#f5b5b5', '#ee8b8b', '#e66767', '#e34948', '#c53433']

def _lum(h):
    """Relative luminance; accepts a hex string or a reportlab Color."""
    if hasattr(h, 'red'):
        r, g, b = h.red, h.green, h.blue
    else:
        r, g, b = [int(str(h).lstrip('#')[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    f = lambda c: c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)

for _arm in (BLUE_ARM, RED_ARM):                   # ramps must be monotone in lightness
    _ls = [_lum(c) for c in _arm]
    assert all(a > b for a, b in zip(_ls, _ls[1:])), f'ramp not monotone: {_arm}'

W, H = A4
ML, MR, MT, MB = 42, 42, 46, 40
CW = W - ML - MR

MEM = set(D['meta']['memorised'])
METRICS = [m['name'] for m in D['meta']['metrics']]
MLAB = {'iptm': 'ipTM', 'chain_pair_iptm': 'chain-pair ipTM', 'ranking_score': 'ranking score',
        'ptm': 'pTM', 'fraction_disordered': 'fraction disordered',
        'chain_pair_pae_min': 'chain-pair PAE min', 'pae_inter_mean': 'interchain PAE mean',
        'pae_inter_p05': 'interchain PAE p05', 'pae_inter_frac_lt10': 'PAE frac < 10 A',
        'pae_inter_frac_lt5': 'PAE frac < 5 A', 'contact_prob_sum': 'contact prob. sum',
        'contact_prob_density': 'contact prob. density', 'contact_prob_max': 'contact prob. max',
        'n_contacts_p50': 'contacts > 0.5', 'contact_participation_ratio': 'contact particip. ratio',
        'plddt_mean': 'mean pLDDT', 'ipsae': 'ipSAE', 'lis': 'LIS'}

CELL = {(c['group'], c['metric']): c for c in D['cells']}
REF = {r['group']: r for r in D['refs']}
TBL = {t['job_id']: t for t in D['table']}
FAMS = ['Pre-registered', 'By pair type', 'By negative type', 'By negative evidence',
        'By positive subgroup', 'Length-adjusted']
GBF = {f: [r['group'] for r in D['refs'] if r['family'] == f] for f in FAMS}
PRE, ARMS = GBF['Pre-registered'], GBF['By pair type']
SHORT = {'1. positives vs real negatives': 'Test 1  17 positives vs 31 real negatives',
         '2a. mutual folding vs hard coiled-coil': 'Test 2a  mutual folding vs hard zipper',
         '2b. fuzzy positives vs real negatives': 'Test 2b  fuzzy positives  (the result)',
         '3. positives vs gold negatives': 'Test 3  vs gold-standard negatives',
         '4. positives vs scrambles': 'Test 4  vs scrambles  (separate arm)',
         'Homotypic only': 'Homotypic arm', 'Heterotypic only': 'Heterotypic arm'}

def cell(g, m): return CELL.get((g, m), {})
def auc_of(g, m): return cell(g, m).get('auc')
def q_of(g, m): return cell(g, m).get('q')
def delta_of(g, m): return cell(g, m).get('delta')
def verd(g, m='iptm'):
    c = cell(g, m)
    return None if 'delta' not in c else ('PASS' if c.get('beats_baseline') else 'FAIL')

def f3(v, sign=False): return '--' if v is None else (f'{v:+.3f}' if sign else f'{v:.3f}')
def fci(ci, sign=False):
    if not ci or ci[0] is None: return '--'
    return f'[{ci[0]:+.2f}, {ci[1]:+.2f}]' if sign else f'[{ci[0]:.2f}, {ci[1]:.2f}]'
def fp(v): return 'n/a' if v is None else ('<0.0001' if v < 0.0001 else f'{v:.4f}')

def div_color(v, center, half):
    if v is None: return MIDG
    t = max(-1.0, min(1.0, (v - center) / half))
    if abs(t) < 0.08: return MIDG
    arm = BLUE_ARM if t >= 0 else RED_ARM
    return HexColor(arm[min(len(arm) - 1, int(abs(t) * len(arm)))])

def ink_on(fill): return HexColor('#ffffff') if _lum(fill) < 0.42 else INK


class Doc:
    def __init__(self):
        self.c = canvas.Canvas(OUT, pagesize=A4)
        self.c.setTitle('AlphaFold 3 IDR-IDR co-folding benchmark - stage 1 results')
        self.c.setSubject('Discrimination of IDR-IDR pairs')
        self.page = 0; self.section = ''; self.y = MT; self.warns = []

    # geometry ---------------------------------------------------
    def Y(self, top): return H - top
    def check(self):
        if self.y > H - MB + 2:
            self.warns.append(f'page {self.page} ({self.section}) overflows by {self.y-(H-MB):.0f}pt')

    def newpage(self, section=None, cover=False):
        if self.page: self.check(); self.c.showPage()
        self.page += 1
        if section is not None: self.section = section
        self.c.setFillColor(SURFACE); self.c.rect(0, 0, W, H, stroke=0, fill=1)
        if cover:
            self.y = MT
        else:
            self.c.setFillColor(ACCENT); self.c.rect(ML, H - 30, 20, 2.2, stroke=0, fill=1)
            self.text(ML + 26, 30, self.section.upper(), 'SB', 7, MUTED, tracking=0.5)
            self.text(W - MR, 30, 'AF3 IDR-IDR benchmark / stage 1', 'S', 7, MUTED, align='r')
            self.y = MT + 14
        self.c.setStrokeColor(GRID); self.c.setLineWidth(0.4)
        self.c.line(ML, MB - 14, W - MR, MB - 14)
        self.text(W - MR, H - MB + 22, str(self.page), 'SB', 7.5, INK2, align='r')
        self.text(ML, H - MB + 22, '60 of 62 jobs scored / 14 September 2026', 'S', 7, MUTED)

    # primitives -------------------------------------------------
    def text(self, x, top, s, font='S', size=8.5, color=INK, align='l', tracking=0):
        s = str(s)
        y = self.Y(top)
        # instrument every drawn string: clipping is caught by measurement, not by eye
        _w = pdfmetrics.stringWidth(s, font, size) + tracking * max(0, len(s) - 1)
        _left = x - (_w if align == 'r' else _w / 2 if align == 'c' else 0)
        if _left + _w > W - MR + 1.0:
            self.warns.append(f'p{self.page} text past right margin by '
                              f'{_left + _w - (W - MR):.1f}pt: {s[:56]!r}')
        if _left < ML - 1.0:
            self.warns.append(f'p{self.page} text past left margin by '
                              f'{ML - _left:.1f}pt: {s[:56]!r}')
        self.c.setFont(font, size); self.c.setFillColor(color)
        if tracking:
            # Letter-spacing is emitted as PDF Tc, which is graphics state: without the
            # save/restore it leaks into every later string on the page and silently
            # widens it past the measured width.
            self.c.saveState()
            to = self.c.beginText(_left, y)
            to.setFont(font, size)
            to.setCharSpace(tracking)
            to.textOut(s)
            self.c.drawText(to)
            self.c.restoreState()
            return
        if align == 'r': self.c.drawRightString(x, y, s)
        elif align == 'c': self.c.drawCentredString(x, y, s)
        else: self.c.drawString(x, y, s)

    def fit(self, s, font, size, maxw):
        s = str(s)
        if pdfmetrics.stringWidth(s, font, size) <= maxw: return s
        while s and pdfmetrics.stringWidth(s + '..', font, size) > maxw: s = s[:-1]
        return s + '..'

    def wrap(self, s, font, size, width):
        out, line = [], ''
        for wd in str(s).split():
            t = (line + ' ' + wd).strip()
            if pdfmetrics.stringWidth(t, font, size) > width and line:
                out.append(line); line = wd
            else: line = t
        if line: out.append(line)
        return out

    def rule(self, top, x0=ML, x1=None, color=GRID, w=0.4):
        self.c.setStrokeColor(color); self.c.setLineWidth(w)
        self.c.line(x0, self.Y(top), (W - MR) if x1 is None else x1, self.Y(top))

    def box(self, x, top, w, h, fill=None, stroke=None, r=0, lw=0.5):
        if fill is not None: self.c.setFillColor(fill)
        if stroke is not None: self.c.setStrokeColor(stroke); self.c.setLineWidth(lw)
        kw = dict(stroke=1 if stroke is not None else 0, fill=1 if fill is not None else 0)
        (self.c.roundRect if r else self.c.rect)(
            *( (x, self.Y(top) - h, w, h, r) if r else (x, self.Y(top) - h, w, h) ), **kw)

    def dot(self, x, top, r, fill, hollow=False, lw=1.6):
        y = self.Y(top)
        if hollow:
            self.c.setFillColor(SURFACE); self.c.setStrokeColor(fill); self.c.setLineWidth(1.4)
        else:
            self.c.setFillColor(fill); self.c.setStrokeColor(SURFACE); self.c.setLineWidth(lw)
        self.c.circle(x, y, r, stroke=1, fill=1)

    # composites -------------------------------------------------
    def h1(self, s, sub=None):
        self.text(ML, self.y, s, 'SB', 15, INK); self.y += 6
        self.rule(self.y, x1=ML + 30, color=ACCENT, w=1.6); self.y += 13
        if sub: self.para(sub, size=8.3, color=INK2)
        self.y += 7

    def h2(self, s, note=None):
        self.y += 3
        self.text(ML, self.y, s, 'SB', 9.4, ACCENT)
        if note:
            nx = ML + pdfmetrics.stringWidth(s, 'SB', 9.4) + 16
            self.text(nx, self.y, self.fit(note, 'SI', 7.4, W - MR - nx - 2), 'SI', 7.4, MUTED)
        self.y += 5; self.rule(self.y); self.y += 11

    def para(self, s, size=8.3, color=INK, font='S', width=None, x=ML, lead=None):
        lead = lead or size * 1.42
        for i, ln in enumerate(self.wrap(s, font, size, width or CW)):
            self.text(x, self.y + i * lead, ln, font, size, color)
        n = len(self.wrap(s, font, size, width or CW))
        self.y += n * lead
        return n * lead

    def takeaway(self, s):
        lines = self.wrap(s, 'S', 8.2, CW - 66)
        h = len(lines) * 11.6 + 9
        self.box(ML, self.y, CW, h, fill=WASH, r=2.5)
        self.text(ML + 8, self.y + 11, 'READS AS', 'SB', 6.3, ACCENT, tracking=0.5)
        for i, ln in enumerate(lines):
            self.text(ML + 58, self.y + 11 + i * 11.6, ln, 'S', 8.2, INK2)
        self.y += h + 9

    def badge(self, x, top, kind, label=None):
        col = {'PASS': GOOD, 'FAIL': CRIT, 'WARN': WARNC, 'NOTE': MUTED}[kind]
        lab = label or kind
        w = pdfmetrics.stringWidth(lab, 'SB', 6.6) + 15
        self.box(x, top, w, 10.5, stroke=col, r=5, lw=0.8)
        self.c.setFillColor(col); self.c.circle(x + 5.6, self.Y(top) - 5.2, 2, stroke=0, fill=1)
        self.text(x + 10.4, top + 7.4, lab, 'SB', 6.6, col)   # baseline inside the pill
        return w

    def tile(self, x, top, w, label, value, sub=None, color=INK, vsize=16):
        self.box(x, top, w, 60, fill=HexColor('#f5f5f2'), r=3)
        self.text(x + 9, top + 14, label.upper(), 'SB', 6.2, MUTED, tracking=0.45)
        self.text(x + 9, top + 33, self.fit(value, 'SB', vsize, w - 18), 'SB', vsize, color)
        if sub:
            for i, ln in enumerate(self.wrap(sub, 'S', 6.6, w - 18)[:2]):
                self.text(x + 9, top + 44 + i * 8, ln, 'S', 6.6, INK2)

    # tables -----------------------------------------------------
    def table(self, cols, rows, size=7.3, rowh=11.4, font='S', zebra=False, head=True, paginate=None):
        tw = sum(c['w'] for c in cols)

        def draw_head():
            self.box(ML, self.y, tw, rowh + 1, fill=HexColor('#f2f1ee'))
            x = ML
            for c in cols:
                al = c.get('align', 'l')
                tx = x + 4 if al == 'l' else (x + c['w'] - 4 if al == 'r' else x + c['w'] / 2)
                self.text(tx, self.y + rowh - 3.2, self.fit(c['title'], 'SB', size - 0.4, c['w'] - 7),
                          'SB', size - 0.4, INK2, align=al)
                x += c['w']
            self.y += rowh + 1
            self.rule(self.y, ML, ML + tw, color=AXIS, w=0.6)

        if head:
            draw_head()
        for i, r in enumerate(rows):
            if paginate and self.y + rowh + 1 > H - MB:      # row (plus a sub-header's extra 2.5pt) won't fit
                self.newpage(paginate['section'])
                if paginate.get('title'):
                    self.text(ML, self.y, paginate['title'], 'SB', 10, INK); self.y += 16
                if head:
                    draw_head()
            if r.get('_sub') is not None:
                self.y += 2.5
                self.text(ML + 1, self.y + rowh - 3.4, r['_sub'], 'SB', size - 0.1, ACCENT)
                self.y += rowh; self.rule(self.y, ML, ML + tw)
                continue
            if zebra and i % 2:
                self.box(ML, self.y, tw, rowh, fill=HexColor('#f8f8f5'))
            x = ML
            for c in cols:
                v = r.get(c['key'])
                fill = c['fill'](r) if callable(c.get('fill')) else None
                if fill is not None:
                    self.box(x + 1, self.y + 1, c['w'] - 2, rowh - 2, fill=fill)
                if c.get('badge'):
                    if v: self.badge(x + 3, self.y + 1, v)
                else:
                    s = c['fmt'](v, r) if c.get('fmt') else ('' if v is None else str(v))
                    col = c['color'](r) if callable(c.get('color')) else c.get('color', INK)
                    if fill is not None: col = ink_on(fill)
                    fnt = c.get('font', font); al = c.get('align', 'l')
                    tx = x + 4 if al == 'l' else (x + c['w'] - 4 if al == 'r' else x + c['w'] / 2)
                    self.text(tx, self.y + rowh - 3.4, self.fit(s, fnt, size, c['w'] - 8),
                              fnt, size, col, align=al)
                x += c['w']
            self.y += rowh
            self.rule(self.y, ML, ML + tw)

    # charts -----------------------------------------------------
    def forest(self, items, center=0.0, span=None, label_w=176, val_w=126, title=None,
               xlabel=None, rowh=14.5):
        if title: self.text(ML, self.y, title, 'SB', 8.2, INK); self.y += 12
        plot_x, plot_w = ML + label_w, CW - label_w - val_w
        vals = [v for it in items if it[1] is not None for v in it[1:4] if v is not None]
        span = span or (max(0.08, max(abs(v - center) for v in vals) * 1.06) if vals else 0.5)
        lo_b, hi_b = center - span, center + span
        sx = lambda v: plot_x + (max(lo_b, min(hi_b, v)) - lo_b) / (hi_b - lo_b) * plot_w
        y0 = self.y
        rows_h = sum(rowh if it[0] is not None else 12 for it in items)
        self.c.setStrokeColor(GRID); self.c.setLineWidth(0.4)
        for t in (lo_b, center - span / 2, center + span / 2, hi_b):
            self.c.line(sx(t), self.Y(y0), sx(t), self.Y(y0 + rows_h))
        self.c.setStrokeColor(AXIS); self.c.setLineWidth(0.9)
        self.c.line(sx(center), self.Y(y0), sx(center), self.Y(y0 + rows_h))
        for it in items:
            lab, est, lo, hi, kind = it
            if lab is None:                                    # family sub-heading
                self.text(ML, self.y + 8, it[4], 'SB', 6.6, ACCENT, tracking=0.4)
                self.y += 12; continue
            mid = self.y + rowh / 2 - 1
            self.text(ML, mid + 2.4, self.fit(lab, 'S', 7.4, label_w - 8), 'S', 7.4, INK)
            if est is None:
                self.text(plot_x + plot_w / 2, mid + 2.4, 'not computable', 'SI', 6.8, MUTED, align='c')
            else:
                pt = POS if est >= center else HexColor('#e34948')
                if lo is not None:
                    self.c.setStrokeColor(pt); self.c.setLineWidth(1.6); self.c.setLineCap(1)
                    self.c.line(sx(lo), self.Y(mid), sx(hi), self.Y(mid))
                self.dot(sx(est), mid, 3.3, pt)
                self.text(W - MR - 4, mid + 2.4,
                          f'{est:+.3f}  {fci([lo, hi], True)}' if lo is not None else f'{est:+.3f}',
                          'S', 6.9, INK2, align='r')
                if kind: self.badge(W - MR - val_w + 2, mid - 3.6, kind)
            self.y += rowh
        self.c.setStrokeColor(AXIS); self.c.setLineWidth(0.5)
        self.c.line(plot_x, self.Y(self.y), plot_x + plot_w, self.Y(self.y))
        for t in (lo_b, center, hi_b):
            self.text(sx(t), self.y + 8, f'{t:+.2f}', 'S', 6.3, MUTED, align='c')
        self.y += 10
        if xlabel:
            self.text(plot_x + plot_w / 2, self.y + 4, xlabel, 'SI', 6.8, MUTED, align='c')
            self.y += 10
        self.y += 6

    def strip(self, metric, h=142, title=None, logscale=False, label_jobs=(), fmt='{:.2f}',
              lo=None, hi=None):
        classes = [('Positive|fuzzy', lambda t: t['label'] == 'positive' and t['job_id'] not in MEM, POS, False),
                   ('Positive|mutual fold', lambda t: t['job_id'] in MEM, POS, True),
                   ('Negative|noncognate', lambda t: t['neg_type'] == 'noncognate', NEG, False),
                   ('Negative|cross-kingdom', lambda t: t['neg_type'] == 'cross_kingdom', NEG, False),
                   ('Negative|charge-repulsive', lambda t: t['neg_type'] == 'charge_repulsive', NEG, False),
                   ('Negative|hard zipper', lambda t: t['neg_type'] == 'hard_coiledcoil', NEG, False),
                   ('Negative|homotypic', lambda t: t['neg_type'] == 'homo_negative', NEG, False),
                   ('Scramble|(chain B shuffled)', lambda t: t['label'] == 'negative_scramble', SCR, False)]
        if title: self.text(ML, self.y, title, 'SB', 8.2, INK); self.y += 12
        vals = [t[metric] for t in D['table'] if t.get(metric) is not None]
        vmin, vmax = min(vals), max(vals)
        tf = (lambda v: math.log10(max(v, max(vmin, 1e-2)))) if logscale else (lambda v: v)
        if lo is not None and hi is not None:          # clean, round axis bounds
            lo_b, hi_b = tf(lo), tf(hi)
        else:
            lo_b, hi_b = tf(vmin), tf(vmax)
            pad = (hi_b - lo_b) * 0.09 or 0.1
            lo_b, hi_b = lo_b - pad, hi_b + pad
        y0 = self.y
        pl, prr = ML + 34, W - MR - 4
        sy = lambda v: y0 + h - (tf(v) - lo_b) / (hi_b - lo_b) * h
        self.c.setStrokeColor(GRID); self.c.setLineWidth(0.4)
        for i in range(6):
            vv = lo_b + (hi_b - lo_b) * i / 5
            yy = y0 + h - i / 5 * h
            self.c.line(pl, self.Y(yy), prr, self.Y(yy))
            raw = (10 ** vv) if logscale else vv
            self.text(pl - 5, yy + 2.3, f'{raw:,.0f}' if abs(raw) >= 100 else f'{raw:.2f}',
                      'S', 6.2, MUTED, align='r')
        bw = (prr - pl) / len(classes)
        for ci, (cl, pred, col, hollow) in enumerate(classes):
            cx = pl + bw * (ci + 0.5)
            sel = sorted([t for t in D['table'] if pred(t) and t.get(metric) is not None],
                         key=lambda t: t[metric])
            for k, t in enumerate(sel):
                self.dot(cx + ((k % 5) - 2) * 4.4, sy(t[metric]), 2.9, col, hollow=hollow)
                if t['job_id'] in label_jobs:
                    # above-right of the dot: to the side it would cross the median tick
                    self.text(cx + ((k % 5) - 2) * 4.4 + 5, sy(t[metric]) - 5.5,
                              t['job_id'], 'SB', 6.1, INK2)
            if sel:
                med = statistics.median(t[metric] for t in sel)
                self.c.setStrokeColor(INK); self.c.setLineWidth(1.3)
                self.c.line(cx - 14, self.Y(sy(med)), cx + 14, self.Y(sy(med)))
                # median value sits outside the dot cloud, left of the tick
                self.text(cx - 17, sy(med) + 2.2, fmt.format(med), 'SB', 6.2, INK, align='r')
            for li, ln in enumerate(cl.split('|')):
                self.text(cx, y0 + h + 10 + li * 7.4, ln, 'S', 6.2, INK2, align='c')
            self.text(cx, y0 + h + 10 + len(cl.split('|')) * 7.4, f'n={len(sel)}', 'S', 6.1, MUTED, align='c')
        self.c.setStrokeColor(AXIS); self.c.setLineWidth(0.5)
        self.c.line(pl, self.Y(y0 + h), prr, self.Y(y0 + h))
        self.y = y0 + h + 33

    def legend(self, items, note=None):
        x = ML
        for (lab, col, hollow) in items:
            self.dot(x + 4, self.y, 2.9, col, hollow=hollow)
            self.text(x + 11, self.y + 2.4, lab, 'S', 6.9, INK2)
            x += pdfmetrics.stringWidth(lab, 'S', 6.9) + 26
        if note: self.text(x + 4, self.y + 2.4, note, 'SI', 6.6, MUTED)
        self.y += 12

    def heatmap(self, groups, glabels, metrics, value_of, q_of_, title=None, center=0.5,
                half=0.35, label_w=150, cellh=13.5, note=None):
        if title: self.text(ML, self.y, title, 'SB', 8.2, INK); self.y += 11
        if note: self.para(note, size=7, color=MUTED); self.y += 4
        cellw = (CW - label_w) / len(metrics)
        self.y += 52                                            # room for rotated headers
        for j, m in enumerate(metrics):
            self.c.saveState()
            self.c.translate(ML + label_w + cellw * (j + 0.66), self.Y(self.y) + 4)
            self.c.rotate(52)
            self.c.setFont('S', 6.1); self.c.setFillColor(INK2)
            self.c.drawString(0, 0, self.fit(MLAB[m], 'S', 6.1, 68))
            self.c.restoreState()
        for g, gl in zip(groups, glabels):
            top = self.y
            self.text(ML, top + cellh - 4.2, self.fit(gl, 'S', 6.7, label_w - 16), 'S', 6.7, INK)
            for j, m in enumerate(metrics):
                v = value_of(g, m)
                fill = div_color(v, center, half)
                x = ML + label_w + cellw * j
                self.box(x + 1, top + 1, cellw - 2, cellh - 2, fill=fill)
                if v is not None:
                    s = f'{v:.2f}' if abs(v) < 10 else f'{v:.1f}'
                    self.text(x + cellw / 2, top + cellh - 4.3, s.replace('0.', '.'),
                              'S', 6.0, ink_on(fill), align='c')
                    qq = q_of_(g, m) if q_of_ else None
                    if qq is not None and qq < 0.05:
                        self.c.setFillColor(ink_on(fill))
                        self.c.circle(x + cellw - 4, self.Y(top + 3.8), 1.2, stroke=0, fill=1)
            self.y += cellh
        self.y += 10

    def ramp_key(self, lo_lab, hi_lab, note=None, w=150):
        steps = list(reversed(RED_ARM)) + ['#f0efec'] + BLUE_ARM
        sw = w / len(steps)
        for i, s in enumerate(steps):
            self.c.setFillColor(HexColor(s))
            self.c.rect(ML + i * sw, self.Y(self.y + 6.5), sw, 6.5, stroke=0, fill=1)
        self.text(ML, self.y + 15, lo_lab, 'S', 6.2, MUTED)
        self.text(ML + w, self.y + 15, hi_lab, 'S', 6.2, MUTED, align='r')
        if note: self.text(ML + w + 14, self.y + 5, note, 'S', 6.4, MUTED)
        self.y += 20

    def bars(self, items, lo=0.0, hi=1.0, center=None, rowh=12.6, title=None, label_w=186,
             fmt='{:.3f}'):
        if title: self.text(ML, self.y, title, 'SB', 8.2, INK); self.y += 11
        px, pw = ML + label_w, CW - label_w - 46
        sx = lambda v: px + (max(lo, min(hi, v)) - lo) / (hi - lo) * pw
        base = center if center is not None else lo
        y0 = self.y
        self.c.setStrokeColor(GRID); self.c.setLineWidth(0.4)
        for t in (lo, (lo + hi) / 2, hi):
            self.c.line(sx(t), self.Y(y0), sx(t), self.Y(y0 + len(items) * rowh))
        for (lab, v, col) in items:
            mid = self.y + rowh / 2
            self.text(ML, mid + 2.4, self.fit(lab, 'S', 7.1, label_w - 8), 'S', 7.1, INK)
            if v is not None:
                x0, x1 = sorted([sx(base), sx(v)])
                self.c.setFillColor(col)
                self.c.roundRect(x0, self.Y(mid + 3.4), max(x1 - x0, 1.2), 6.8, 2.2, stroke=0, fill=1)
                self.text(px + pw + 5, mid + 2.4, fmt.format(v), 'S', 6.9, INK2)
            self.y += rowh
        self.c.setStrokeColor(AXIS); self.c.setLineWidth(0.7)
        self.c.line(sx(base), self.Y(self.y), sx(base), self.Y(y0))
        self.c.setLineWidth(0.5)
        self.c.line(px, self.Y(self.y), px + pw, self.Y(self.y))
        for t in (lo, (lo + hi) / 2, hi):
            self.text(sx(t), self.y + 8, f'{t:.2f}', 'S', 6.3, MUTED, align='c')
        self.y += 15

    def dumbbell(self, items, lo, hi, rowh=12.6, title=None, label_w=150, a_lab='a', b_lab='b',
                 fmt='{:.2f}'):
        if title: self.text(ML, self.y, title, 'SB', 8.2, INK); self.y += 11
        px, pw = ML + label_w, CW - label_w - 104
        sx = lambda v: px + (max(lo, min(hi, v)) - lo) / (hi - lo) * pw
        light, dark = HexColor(BLUE_ARM[1]), HexColor(BLUE_ARM[4])
        y0 = self.y
        self.c.setStrokeColor(GRID); self.c.setLineWidth(0.4)
        for t in (lo, (lo + hi) / 2, hi):
            self.c.line(sx(t), self.Y(y0), sx(t), self.Y(y0 + len(items) * rowh))
        for (lab, a, b) in items:
            mid = self.y + rowh / 2
            self.text(ML, mid + 2.4, self.fit(lab, 'S', 7.1, label_w - 8), 'S', 7.1, INK)
            if a is not None and b is not None:
                self.c.setStrokeColor(AXIS); self.c.setLineWidth(1.2)
                self.c.line(sx(a), self.Y(mid), sx(b), self.Y(mid))
                self.dot(sx(a), mid, 2.9, light); self.dot(sx(b), mid, 2.9, dark)
                self.text(px + pw + 7, mid + 2.4, f'{fmt.format(a)} -> {fmt.format(b)}', 'S', 6.7, INK2)
            self.y += rowh
        self.c.setStrokeColor(AXIS); self.c.setLineWidth(0.5)
        self.c.line(px, self.Y(self.y), px + pw, self.Y(self.y))
        for t in (lo, (lo + hi) / 2, hi):
            self.text(sx(t), self.y + 8, f'{t:.2f}', 'S', 6.3, MUTED, align='c')
        self.y += 14
        self.legend([(a_lab, light, False), (b_lab, dark, False)])


# ================================================================ pages
def p_cover(d):
    d.newpage('summary', cover=True)
    d.text(ML, d.y, 'ALPHAFOLD 3  /  INTRINSICALLY DISORDERED REGIONS', 'SB', 7.4, ACCENT, tracking=1.1)
    d.y += 26
    d.text(ML, d.y, 'Can co-folding confidence tell', 'SB', 25, INK); d.y += 27
    d.text(ML, d.y, 'interacting IDR pairs apart?', 'SB', 25, INK); d.y += 20
    d.rule(d.y, x1=ML + 56, color=ACCENT, w=2.2); d.y += 18
    d.para('Stage 1 of the IDR-IDR benchmark: 60 of 62 AlphaFold Server jobs scored, five samples '
           'each. The pre-registered primary metric is ipTM of the top-ranked model and the '
           'pre-registered comparator is a leave-one-out logistic model on sequence composition '
           'alone - not chance. A test passes only if the paired-bootstrap CI on '
           'delta-AUC excludes zero.', size=8.6, color=INK2)
    d.y += 14

    # hero answer
    hero = ('Every pre-registered test fails its own criterion, and in the homotypic arm AF3 is '
            'significantly worse than the trivial baseline.')
    hero_lines = d.wrap(hero, 'S', 8, CW - 30)
    h = 46 + len(hero_lines) * 11
    d.box(ML, d.y, CW, h, fill=HexColor('#f7f2f2'), r=3)
    d.c.setFillColor(CRIT); d.c.rect(ML, d.Y(d.y + h), 3, h, stroke=0, fill=1)
    d.text(ML + 14, d.y + 17, 'THE ANSWER, ON THESE 60 JOBS', 'SB', 6.4, CRIT, tracking=0.5)
    d.text(ML + 14, d.y + 38, 'No - ipTM does not beat sequence composition.', 'SB', 14.5, INK)
    for i, ln in enumerate(hero_lines):
        d.text(ML + 14, d.y + 53 + i * 11, ln, 'S', 8, INK2)
    d.y += h + 16

    tiles = [('jobs scored', '60 / 62', 'S20, S21 pending', INK),
             ('test 1  delta AUC', f3(delta_of(PRE[0], 'iptm'), True),
              fci(cell(PRE[0], 'iptm')['delta_ci'], True) + '  FAIL', CRIT),
             ('test 2b  "the result"', f3(delta_of(PRE[2], 'iptm'), True),
              fci(cell(PRE[2], 'iptm')['delta_ci'], True) + '  FAIL', CRIT),
             ('homotypic arm', f3(delta_of('Homotypic only', 'iptm'), True),
              'CI excludes zero: worse', CRIT),
             ('GCN4 canary P21', f'{TBL["P21"]["iptm"]:.2f}', 'pipeline behaves', GOOD)]
    tw = (CW - 4 * 8) / 5
    for i, (lab, val, sub, col) in enumerate(tiles):
        d.tile(ML + i * (tw + 8), d.y, tw, lab, val, sub, col, vsize=15)
    d.y += 60 + 18

    d.h2('The pre-registered tests', 'primary metric ipTM; comparator = composition baseline')
    rows = []
    for g in PRE + ARMS:
        c, rf = cell(g, 'iptm'), REF[g]
        rows.append({'g': SHORT[g], 'n': f"{c['n_pos']}/{c['n_neg']}", 'auc': c.get('auc'),
                     'ci': c.get('ci'), 'base': c.get('base_auc'), 'd': c.get('delta'),
                     'dci': c.get('delta_ci'), 'v': verd(g)})
    d.table([
        {'title': 'test / group', 'key': 'g', 'w': 176},
        {'title': 'n pos/neg', 'key': 'n', 'w': 46, 'align': 'c'},
        {'title': 'AF3 AUC', 'key': 'auc', 'w': 42, 'align': 'r', 'fmt': lambda v, r: f3(v),
         'font': 'SB'},
        {'title': '95% CI', 'key': 'ci', 'w': 62, 'align': 'c', 'fmt': lambda v, r: fci(v)},
        {'title': 'baseline', 'key': 'base', 'w': 44, 'align': 'r', 'fmt': lambda v, r: f3(v)},
        {'title': 'delta', 'key': 'd', 'w': 40, 'align': 'r', 'fmt': lambda v, r: f3(v, True),
         'font': 'SB', 'color': lambda r: CRIT if (r['d'] or 0) < 0 else INK},
        {'title': 'delta 95% CI', 'key': 'dci', 'w': 66, 'align': 'c',
         'fmt': lambda v, r: fci(v, True)},
        {'title': 'verdict', 'key': 'v', 'w': 35, 'badge': True},
    ], rows, rowh=12.4)
    d.y += 10
    d.para('Test 4 has no baseline column by design: scrambles are composition-matched, so the '
           'composition model cannot separate them from their parents and is not a meaningful '
           'comparator in that arm.', size=7.2, color=MUTED)
    d.y += 12
    d.takeaway('The plan\'s own decision table says: delta-AUC_fuzzy near zero with AUC_fold high is the '
               'clean expected split; both low means debug first. Here delta-AUC_fuzzy is -0.11 with a CI '
               'crossing zero, and AUC_fold is 0.42 - not high. Pages 3 to 4 show why that is not a '
               'broken pipeline but a real failure mode: AF3 dimerises any zipper-competent pair, '
               'cognate or not.')


def p_forest(d):
    d.newpage('every grouping')
    d.h1('Does AF3 beat composition? Every grouping.',
         'Delta-AUC = AUC(AF3 ipTM) - AUC(composition baseline), paired bootstrap, 5,000 resamples. '
         'Right of the zero line, AF3 adds information over composition; left, it subtracts. A group '
         'passes only if the whole interval sits right of zero.')
    items = []
    for fam in ['Pre-registered', 'By pair type', 'By negative type', 'By negative evidence',
                'By positive subgroup']:
        items.append((None, None, None, None, fam))
        for g in GBF[fam]:
            c = cell(g, 'iptm')
            if 'delta' not in c:
                continue
            lab = SHORT.get(g, g) + f"  ({c['n_pos']}/{c['n_neg']})"
            items.append((lab, c['delta'], c['delta_ci'][0], c['delta_ci'][1], verd(g)))
    d.forest(items, center=0.0, label_w=196, val_w=118,
             xlabel='delta-AUC  (AF3 ipTM minus composition baseline)', rowh=13.2)
    d.takeaway('Exactly one grouping clears zero on ipTM, and it is the pipeline control: the four '
               'memorised complexes, at +0.52. Of the other 23, two exclude zero on the wrong side - '
               'the homotypic arm (-0.34) and the gold-standard negatives (-0.52), where AF3 ranks '
               'the negatives above the positives - and the rest are indistinguishable from '
               'composition.')


def p_pipeline(d):
    d.newpage('pipeline integrity')
    d.h1('Is the pipeline sound?',
         'Before reading a null result, the controls have to hold. They do: the canary folds, the '
         'memorised complexes score highest in the set, and every job ran the sequence it was '
         'supposed to.')
    d.h2('Memorised controls and the canary', 'pipeline controls, excluded from every headline number')
    rows = []
    for j in ['P21', 'P20', 'P07', 'P06']:
        t = TBL[j]
        rows.append({'j': f"{j}  {t['chain_A']} x {t['chain_B']}", 'pdb': t['in_pdb'] or '--',
                     'len': int(float(t['total_residues'])), 'iptm': t['iptm'],
                     'plddt': t['plddt_mean'], 'pae': t['pae_inter_mean'],
                     'rank': t['ranking_score'], 'sd': t['iptm_sample_sd']})
    d.table([
        {'title': 'control', 'key': 'j', 'w': 168},
        {'title': 'PDB', 'key': 'pdb', 'w': 40, 'align': 'c'},
        {'title': 'residues', 'key': 'len', 'w': 44, 'align': 'r'},
        {'title': 'ipTM', 'key': 'iptm', 'w': 42, 'align': 'r', 'fmt': lambda v, r: f'{v:.2f}',
         'font': 'SB'},
        {'title': 'mean pLDDT', 'key': 'plddt', 'w': 54, 'align': 'r', 'fmt': lambda v, r: f'{v:.1f}'},
        {'title': 'interchain PAE', 'key': 'pae', 'w': 62, 'align': 'r', 'fmt': lambda v, r: f'{v:.1f}'},
        {'title': 'ranking', 'key': 'rank', 'w': 44, 'align': 'r', 'fmt': lambda v, r: f'{v:.2f}'},
        {'title': 'sample sd', 'key': 'sd', 'w': 47, 'align': 'r', 'fmt': lambda v, r: f'{v:.3f}'},
    ], rows, rowh=12.2)
    d.y += 8
    d.para('GCN4-p1 (P21, 66 residues) is the pre-registered canary: if AF3 could not dimerise it '
           'confidently the instruction was to stop and debug. It reaches ipTM 0.73 with mean pLDDT 88 '
           'and a 4.1 A interchain PAE, on an MSA of only 23 paired sequences. All four memorised '
           'complexes land between 0.70 and 0.80, the top of the whole set.', size=7.8, color=INK2)
    d.y += 12

    d.h2('Data integrity', 'checked against the job requests the server returned')
    checks = [
        ('Submitted sequences match the pair table', 'all 60 jobs, both chains exact', 'PASS'),
        ('Samples per job', '5 of 5 for every job (300 summary files, 300 full-data files)', 'PASS'),
        ('Job-to-row matching', '60 matched by name, 0 ambiguous', 'PASS'),
        ('Model seeds', '59 jobs on seed 1; N24 on a server-assigned seed', 'WARN'),
        ('Steric clashes', 'a clash in some sample of P13 and P17 only (2 of 60)', 'WARN'),
        ('Structure templates', 'on for all 95 chain entries - server default, not a choice', 'WARN'),
        ('Jobs outstanding', 'S20 and S21 (batch 4) not yet run', 'NOTE'),
    ]
    rows = [{'c': a, 'd': b, 'v': c} for a, b, c in checks]
    d.table([{'title': 'check', 'key': 'c', 'w': 210},
             {'title': 'result', 'key': 'd', 'w': 256},
             {'title': '', 'key': 'v', 'w': 45, 'badge': True}], rows, rowh=12.2)
    d.y += 12

    d.h2('Sample spread', 'five samples from one seed, so this is sample spread, not seed replication')
    sds = {k: [TBL[j]['iptm_sample_sd'] for j in TBL
               if (TBL[j]['binding_mode'] == k and TBL[j]['label'] == 'positive')]
           for k in ('fuzzy', 'mutual_folding')}
    d.para(f"Fuzzy positives: mean within-job ipTM sd {statistics.mean(sds['fuzzy']):.4f} "
           f"(n={len(sds['fuzzy'])}).  Mutual-folding positives: {statistics.mean(sds['mutual_folding']):.4f} "
           f"(n={len(sds['mutual_folding'])}). The plan expected fuzzy pairs to be the more variable; "
           'they are not, because all five samples come from one seed and share a trunk embedding. '
           'Seed replicates are still needed to measure run-to-run variance.', size=7.8, color=INK2)
    d.y += 12
    d.takeaway('Nothing here explains the null result away. The controls behave, so the failures on '
               'pages 1 and 2 are about AF3\'s discrimination, not about the harness - with one '
               'caveat that cuts the other way: templates were on for every chain, which can only '
               'help the controls look good.')


def p_primary(d):
    d.newpage('primary metric')
    d.h1('Where every job actually lands',
         'One dot per job, ipTM of the top-ranked sample. Black tick is the class median. The dots are '
         'the whole dataset: no aggregation, nothing hidden.')
    # label sparingly: N15/N16 sit within 0.02 of N17 and their labels would collide
    d.strip('iptm', h=150, lo=0.0, hi=1.0, label_jobs={'P21', 'N17', 'N26'})
    d.legend([('interacting (positive)', POS, False), ('memorised control (hollow)', POS, True),
              ('non-interacting (negative)', NEG, False), ('scramble', SCR, False)],
             note='labelled: canary, top negatives')
    d.y += 4
    _fz = [t['iptm'] for t in D['table'] if t['label'] == 'positive' and t['job_id'] not in MEM]
    _hz = [TBL[j]['iptm'] for j in ('N15', 'N16', 'N17')]
    d.takeaway(f'The three hard-zipper negatives (N15, N16, N17) sit at {min(_hz):.2f}-{max(_hz):.2f} - '
               f'level with the memorised true complexes - while the 17 genuine fuzzy positives have a '
               f'median ipTM of {statistics.median(_fz):.2f} (range {min(_fz):.2f}-{max(_fz):.2f}). AF3 '
               'confidence tracks "can these two helices form a coiled coil", not "do these two '
               'regions interact".')
    d.h2('Primary metric, full statistics', 'permutation p, BH-corrected across the 16-metric panel')
    rows = []
    for g in PRE + ARMS:
        c = cell(g, 'iptm')
        rows.append({'g': SHORT[g], 'n': f"{c['n_pos']}/{c['n_neg']}", 'auc': c.get('auc'),
                     'ci': c.get('ci'), 'p': c.get('p'), 'q': c.get('q'),
                     'len': REF[g].get('length_auc'), 'base': c.get('base_auc')})
    d.table([
        {'title': 'group', 'key': 'g', 'w': 166},
        {'title': 'n', 'key': 'n', 'w': 44, 'align': 'c'},
        {'title': 'ipTM AUC', 'key': 'auc', 'w': 48, 'align': 'r', 'fmt': lambda v, r: f3(v),
         'font': 'SB', 'fill': lambda r: div_color(r['auc'], 0.5, 0.35)},
        {'title': '95% CI', 'key': 'ci', 'w': 56, 'align': 'c', 'fmt': lambda v, r: fci(v)},
        {'title': 'perm p', 'key': 'p', 'w': 46, 'align': 'r', 'fmt': lambda v, r: fp(v)},
        {'title': 'BH q', 'key': 'q', 'w': 46, 'align': 'r', 'fmt': lambda v, r: fp(v)},
        {'title': 'composition', 'key': 'base', 'w': 50, 'align': 'r', 'fmt': lambda v, r: f3(v)},
        {'title': 'length alone', 'key': 'len', 'w': 49, 'align': 'r', 'fmt': lambda v, r: f3(v),
         'color': lambda r: CRIT if (r['len'] or 0) > 0.8 else INK2},
    ], rows, rowh=12.2)
    d.y += 8
    d.para('The last two columns are the reference predictors on the same jobs: a composition-only '
           'model, and total chain length on its own. Where they beat the ipTM column, AF3 has added '
           'nothing a sequence length could not.', size=7.2, color=MUTED)


def p_arms(d):
    d.newpage('stratified')
    d.h1('Homotypic and heterotypic, analysed apart',
         'Plan criterion 3. The two arms have different baselines and different MSA situations, so '
         'pooling them lets chain identity leak into the answer. Each arm is scored against a '
         'baseline refitted inside that arm.')
    rows = []
    for g in ARMS:
        c = cell(g, 'iptm'); cs = cell(g, 'contact_prob_sum')
        rows.append({'g': SHORT[g], 'n': f"{c['n_pos']}/{c['n_neg']}", 'auc': c['auc'],
                     'ci': c['ci'], 'base': c['base_auc'], 'd': c['delta'], 'dci': c['delta_ci'],
                     'v': verd(g), 'cs': cs['auc'], 'csd': cs['delta']})
    d.table([
        {'title': 'arm', 'key': 'g', 'w': 104},
        {'title': 'n', 'key': 'n', 'w': 40, 'align': 'c'},
        {'title': 'ipTM AUC', 'key': 'auc', 'w': 48, 'align': 'r', 'fmt': lambda v, r: f3(v), 'font': 'SB'},
        {'title': '95% CI', 'key': 'ci', 'w': 58, 'align': 'c', 'fmt': lambda v, r: fci(v)},
        {'title': 'arm baseline', 'key': 'base', 'w': 56, 'align': 'r', 'fmt': lambda v, r: f3(v)},
        {'title': 'delta', 'key': 'd', 'w': 40, 'align': 'r', 'fmt': lambda v, r: f3(v, True),
         'font': 'SB', 'color': lambda r: CRIT if r['d'] < 0 else INK},
        {'title': 'delta 95% CI', 'key': 'dci', 'w': 64, 'align': 'c', 'fmt': lambda v, r: fci(v, True)},
        {'title': 'verdict', 'key': 'v', 'w': 35, 'badge': True},
        {'title': 'contact sum AUC', 'key': 'cs', 'w': 66, 'align': 'r', 'fmt': lambda v, r: f3(v)},
    ], rows, rowh=12.4)
    d.y += 12
    d.takeaway('The homotypic arm is where the benchmark is decided, and it is where AF3 loses '
               'outright: ipTM 0.56 against a composition baseline of 0.90, delta -0.34 with the '
               'interval clear of zero. In the heterotypic arm composition is useless (0.50) and AF3 '
               'is no better (0.49) - on only 5 positives.')

    d.h2('Why composition is so strong in the homotypic arm',
         'single-feature AUC, two-sided, from baseline_leakage.py')
    for title, key in [('Homotypic pairs  (14 positive / 11 negative)', 'homotypic'),
                       ('Heterotypic pairs  (7 positive / 20 negative)', 'heterotypic')]:
        d.text(ML, d.y, title, 'SB', 7.6, INK); d.y += 11
        items = [(x['feature'], x['auc_two_sided'],
                  HexColor('#e34948') if x['auc_two_sided'] > 0.75 else POS)
                 for x in D['single_features'][key][:6]]
        d.bars(items, lo=0.5, hi=1.0, rowh=11.4, label_w=150, fmt='{:.3f}')
        d.y += 2
    d.para('A single feature above 0.75, or a leave-one-out model above 0.80, is the documented '
           'leakage threshold. In the homotypic arm five features clear it on their own: LC domains '
           'that self-associate really are aromatic-rich and weakly charged. That is biology, not a '
           'curation error - but it means a homotypic claim needs to beat 0.85, not 0.5.',
           size=7.6, color=INK2)


def _group_table(d, fam, extra_note=None):
    rows = []
    for g in GBF[fam]:
        c, rf = cell(g, 'iptm'), REF[g]
        cs = cell(g, 'contact_prob_sum')
        rows.append({'g': g.replace('vs ', ''), 'n': f"{c['n_pos']}/{c['n_neg']}",
                     'auc': c.get('auc'), 'ci': c.get('ci'), 'q': c.get('q'),
                     'base': c.get('base_auc'), 'd': c.get('delta'), 'dci': c.get('delta_ci'),
                     'len': rf.get('length_auc'), 'cs': cs.get('auc'), 'v': verd(g),
                     'small': c.get('small_n')})
    d.table([
        {'title': 'group', 'key': 'g', 'w': 118},
        {'title': 'n', 'key': 'n', 'w': 34, 'align': 'c'},
        {'title': 'ipTM AUC', 'key': 'auc', 'w': 42, 'align': 'r', 'fmt': lambda v, r: f3(v),
         'font': 'SB', 'fill': lambda r: div_color(r['auc'], 0.5, 0.35)},
        {'title': '95% CI', 'key': 'ci', 'w': 52, 'align': 'c', 'fmt': lambda v, r: fci(v)},
        {'title': 'BH q', 'key': 'q', 'w': 40, 'align': 'r',
         'fmt': lambda v, r: 'n<3' if r.get('small') else fp(v)},
        {'title': 'compos.', 'key': 'base', 'w': 38, 'align': 'r', 'fmt': lambda v, r: f3(v)},
        {'title': 'delta', 'key': 'd', 'w': 36, 'align': 'r', 'fmt': lambda v, r: f3(v, True),
         'color': lambda r: CRIT if (r['d'] or 0) < 0 else INK},
        {'title': 'delta CI', 'key': 'dci', 'w': 58, 'align': 'c', 'fmt': lambda v, r: fci(v, True)},
        {'title': 'length', 'key': 'len', 'w': 36, 'align': 'r', 'fmt': lambda v, r: f3(v),
         'color': lambda r: CRIT if (r['len'] or 0) > 0.8 else INK2},
        {'title': 'contact sum', 'key': 'cs', 'w': 42, 'align': 'r', 'fmt': lambda v, r: f3(v)},
    ], rows, rowh=12.2)
    if extra_note:
        d.y += 7; d.para(extra_note, size=7.2, color=MUTED)


def p_negtypes(d):
    d.newpage('negative subtypes')
    d.h1('By type of negative',
         'The plan asks for negative subtypes separately, because they are not equally clean: '
         'non-cognate human pairs are the fair test, cross-kingdom pairs share the scramble arm\'s '
         'missing-MSA problem in milder form, and the hard zippers are adversarial by design.')
    _group_table(d, 'By negative type',
                 'Charge-repulsive has only 2 negatives, so no permutation test is reported for it. '
                 'Scrambles have no composition comparator by construction.')
    d.y += 10
    _cr, _ck = auc_of('vs charge_repulsive', 'iptm'), auc_of('vs cross_kingdom', 'iptm')
    _nc, _hc = auc_of('vs noncognate', 'iptm'), auc_of('vs hard_coiledcoil', 'iptm')
    _hn = auc_of('vs homo_negative', 'iptm')
    d.takeaway(f'The subtype decides the answer. AF3 scores highest against the negatives whose second '
               f'chain has the weakest evolutionary signal - charge-repulsive ({_cr:.2f}, only 2 pairs) '
               f'and cross-kingdom ({_ck:.2f}) - holds up against non-cognate human pairs ({_nc:.2f}), '
               f'then inverts completely against the hard zippers ({_hc:.2f}, where length alone scores '
               f'a perfect {REF["vs hard_coiledcoil"]["length_auc"]:.2f} because every zipper pair is '
               f'short) and against the homotypic negatives ({_hn:.2f}). The cleanest negatives and the '
               'nastiest ones disagree, which is the whole reason the subtypes are reported apart.')
    d.h2('The same groups, by what the model could be reading instead')
    items = []
    for g in GBF['By negative type']:
        rf = REF[g]
        items.append((g.replace('vs ', '') + f"  (n={rf['n_pos']}/{rf['n_neg']})",
                      rf.get('length_auc'), HexColor('#e34948') if (rf.get('length_auc') or 0) > 0.8 else POS))
    d.bars(items, lo=0.0, hi=1.0, center=0.5, rowh=12.4, label_w=196,
           title='Total chain length alone, as a classifier (0.5 = chance)')
    d.para('Bars to the right of the centre line mean the positives are longer; to the left, the '
           'negatives are. Two of these groups are separated almost perfectly by length before AF3 '
           'is consulted at all.', size=7.2, color=MUTED)


def p_evidence(d):
    d.newpage('evidence grading')
    d.h1('By how well we know the negative is negative',
         'Every negative carries an evidence grade. 21 of 31 rest on inferred or absence-only '
         'evidence; 8 are demonstrated. The plan says that if the headline and the gold-standard '
         'subset disagree, believe the gold standard.')
    _group_table(d, 'By negative evidence',
                 'Uncertain has 2 rows (H1.0 alone, ERD10 alone) and is descriptive only.')
    d.y += 10
    d.takeaway('They disagree violently - and the gold standard is the confounded one. Against the 8 '
               'demonstrated negatives ipTM inverts to 0.12, but those 8 negatives are 94-236 '
               'residues while the 17 positives are 205-882, so length alone scores 0.985 on that '
               'exact split. The inversion is a length and foldability artefact, not evidence that '
               'AF3 ranks real complexes below non-complexes.')
    d.h2('What the gold-standard subset is made of', 'the 8 demonstrated negatives, ranked by ipTM')
    rows = []
    for t in sorted([TBL[j] for j in TBL if TBL[j]['label'] == 'negative'
                     and 'demonstrated' in TBL[j]['evidence_strength']],
                    key=lambda t: -t['iptm']):
        rows.append({'j': f"{t['job_id']}  {t['chain_A']} x {t['chain_B']}",
                     'nt': t['neg_type'], 'len': int(float(t['total_residues'])),
                     'iptm': t['iptm'], 'pae': t['pae_inter_mean'], 'plddt': t['plddt_mean'],
                     'cs': t['contact_prob_sum'], 'ev': t['evidence_strength'].replace('demonstrated', 'dem.')})
    d.table([{'title': 'gold-standard negative', 'key': 'j', 'w': 150},
             {'title': 'subtype', 'key': 'nt', 'w': 74},
             {'title': 'evidence', 'key': 'ev', 'w': 72},
             {'title': 'residues', 'key': 'len', 'w': 40, 'align': 'r'},
             {'title': 'ipTM', 'key': 'iptm', 'w': 38, 'align': 'r', 'fmt': lambda v, r: f'{v:.2f}',
              'font': 'SB', 'fill': lambda r: div_color(r['iptm'], 0.45, 0.4)},
             {'title': 'PAE', 'key': 'pae', 'w': 36, 'align': 'r', 'fmt': lambda v, r: f'{v:.1f}'},
             {'title': 'pLDDT', 'key': 'plddt', 'w': 40, 'align': 'r', 'fmt': lambda v, r: f'{v:.0f}'},
             {'title': 'contact sum', 'key': 'cs', 'w': 44, 'align': 'r', 'fmt': lambda v, r: f'{v:.0f}'},
             ], rows, rowh=12.2)
    d.y += 7
    d.para('Six of the eight are short helix- or zipper-competent constructs. Three of them (N15, '
           'N16, N17) are exactly the adversarial coiled-coil negatives, and three more (N26, N27, '
           'N28) are single chains of a real complex run against themselves.', size=7.4, color=INK2)


def p_positives(d):
    d.newpage('positive subgroups')
    d.h1('By type of positive', 'All rows below are scored against the same 31 real-sequence '
         'negatives, so the groups are directly comparable.')
    _group_table(d, 'By positive subgroup',
                 'The LLPS rows are a stratification, never an outcome: a 1:1 co-fold cannot report '
                 'on phase separation, which is many-body and concentration-dependent.')
    d.y += 10
    _mem = auc_of('Mutual-folding (memorised) positives', 'iptm')
    _fhet = auc_of('Fuzzy heterotypic positives', 'iptm')
    _fhom = auc_of('Fuzzy homotypic positives', 'iptm')
    _llps = auc_of('LLPS-reported positives', 'iptm')
    _cps = cell('Fuzzy homotypic positives', 'contact_prob_sum')
    d.takeaway(
        f'The four memorised complexes are the positives AF3 separates cleanly (ipTM AUC {_mem:.2f}); '
        f'fuzzy positives sit at {_fhet:.2f} heterotypic and {_fhom:.2f} homotypic, and reported LLPS '
        f'makes no difference ({_llps:.2f}). Memorisation, not disorder, is what ipTM tracks here. One '
        f'exception is worth flagging honestly: for the 12 fuzzy homotypic positives, contact '
        f'probability sum beats composition by {_cps["delta"]:+.3f} '
        f'[{_cps["delta_ci"][0]:+.3f}, {_cps["delta_ci"][1]:+.3f}] - the only cell outside the '
        'memorised controls to clear zero in all 464. It is exploratory, uncorrected for the 464 '
        'comparisons, and the pre-registered homotypic arm (which uses homotypic negatives and the '
        'arm baseline) does not pass. Treat it as the hypothesis to test on a fresh set, not a result.')
    d.h2('The sharpest tests in the set', 'matched internal contrasts - same chains, same MSAs, only the pairing differs')
    d.para('These three comparisons do not rest on a literature absence at all. ACTR and NCBD each '
           'fold only with the other, and c-Fos forms no stable homodimer, so ACTR x ACTR, NCBD x '
           'NCBD and Fos x Fos are negatives whose chains also appear in a positive. The plan calls '
           'the outcome unambiguous either way.', size=7.8, color=INK2)
    d.y += 10
    rows = []
    for m in D['matched_internal']:
        pos = TBL[m['positive']]
        rows.append({'_sub': f"{m['positive']}  {pos['chain_A']} x {pos['chain_B']}   -   {m['why']}"})
        for j in [m['positive']] + m['negatives']:
            t = TBL[j]
            rows.append({'j': f"{j}  {t['chain_A']} x {t['chain_B']}",
                         'lab': 'positive' if t['label'] == 'positive' else 'negative',
                         'iptm': t['iptm'], 'pae': t['pae_inter_mean'], 'cs': t['contact_prob_sum'],
                         'plddt': t['plddt_mean'], 'len': int(float(t['total_residues']))})
    d.table([{'title': 'job', 'key': 'j', 'w': 196},
             {'title': 'label', 'key': 'lab', 'w': 56},
             {'title': 'residues', 'key': 'len', 'w': 46, 'align': 'r'},
             {'title': 'ipTM', 'key': 'iptm', 'w': 44, 'align': 'r', 'fmt': lambda v, r: f'{v:.2f}',
              'font': 'SB', 'fill': lambda r: div_color(r['iptm'], 0.45, 0.4)},
             {'title': 'interchain PAE', 'key': 'pae', 'w': 62, 'align': 'r', 'fmt': lambda v, r: f'{v:.1f}'},
             {'title': 'mean pLDDT', 'key': 'plddt', 'w': 54, 'align': 'r', 'fmt': lambda v, r: f'{v:.1f}'},
             {'title': 'contact sum', 'key': 'cs', 'w': 53, 'align': 'r', 'fmt': lambda v, r: f'{v:.0f}'},
             ], rows, rowh=12.2)
    d.y += 8
    n06 = sum(D['matched_internal'][0]['correct'].values())
    n07 = sum(D['matched_internal'][1]['correct'].values())
    d.takeaway(f'AF3 gets all three orderings right on ipTM, and on {n06} of 16 metrics for P06 and '
               f'{n07} of 16 for P07. This is the one genuinely encouraging result in stage 1 - with '
               'the caveat on page 11 that both positives were handed per-chain templates taken from '
               'their own cognate complex.')


def p_secondary(d):
    d.newpage('secondary panel')
    d.h1('The whole 16-metric panel',
         'ipTM may simply be floored for disordered pairs, so the plan pre-registered a 16-metric '
         'panel with permutation p-values and Benjamini-Hochberg correction. Every metric is '
         'oriented so that a higher AUC means "more interaction".')
    groups = PRE + ARMS
    HM = {PRE[0]: 'Test 1  (17/31)', PRE[1]: 'Test 2a  (4/3)', PRE[2]: 'Test 2b  (17/28)',
          PRE[3]: 'Test 3  (17/8)', PRE[4]: 'Test 4  (17/8)',
          'Homotypic only': 'Homotypic arm  (12/11)', 'Heterotypic only': 'Heterotypic arm  (5/20)'}
    d.heatmap(groups, [HM[g] for g in groups], METRICS, auc_of, q_of,
              title='AUC by metric and group', center=0.5, half=0.35, cellh=13.6,
              note='A dot in the corner marks BH q < 0.05 within that group. Cell values are AUCs; '
                   '.50 is chance.')
    d.ramp_key('0.15  inverted', 'discriminating  0.85',
               note='AUC, diverging around chance (0.50)')
    d.y += 6
    d.heatmap(groups, [HM[g] for g in groups], METRICS, delta_of, None,
              title='Delta-AUC against the composition baseline, same cells', center=0.0, half=0.3,
              cellh=13.6, note='Blank row for test 4: scrambles have no composition comparator.')
    d.ramp_key('-0.30  worse than composition', 'better  +0.30', note='delta-AUC')
    d.y += 8
    best = max(((g, m, delta_of(g, m)) for g in groups for m in METRICS
                if delta_of(g, m) is not None), key=lambda t: t[2])
    d.takeaway('Read the two grids together. Several metrics beat chance in the top grid - contact '
               'probability sum reaches 0.79 at q = 0.008, and 0.97 in the homotypic arm - yet the '
               'bottom grid is essentially all grey or red: across these seven groups not one metric '
               f'beats composition with a CI clear of zero. The best single cell is {MLAB[best[1]]} '
               f'at delta {best[2]:+.3f} (test 2a, n=4/3), and its interval still crosses zero. The '
               'one non-control exception anywhere in the 464 cells is on page 8.')


def p_length(d):
    d.newpage('confound: length')
    d.h1('Confound 1: chain length',
         'The plan flags any metric whose Spearman rho against total residue count exceeds 0.5, and '
         'requires it to be reported partialled as well as raw. Four metrics are flagged - and they '
         'are exactly the ones that carry the gold-standard inversion.')
    items = [(MLAB[c['metric']], c['rho_length'],
              HexColor('#e34948') if abs(c['rho_length']) > 0.5 else POS)
             for c in sorted(D['confounds'], key=lambda c: -abs(c['rho_length']))]
    d.bars(items, lo=-1.0, hi=1.0, center=0.0, rowh=11.2, label_w=150,
           title='Spearman rho between metric and total residue count, all 60 jobs', fmt='{:+.3f}')
    d.y += 4
    d.h2('Length alone, as a classifier, in every group')
    items = []
    for fam in ['Pre-registered', 'By pair type', 'By negative evidence']:
        for g in GBF[fam]:
            rf = REF[g]
            items.append((SHORT.get(g, g.replace('vs ', '')) + f"  ({rf['n_pos']}/{rf['n_neg']})",
                          rf.get('length_auc'),
                          HexColor('#e34948') if (rf.get('length_auc') or 0) > 0.8 else POS))
    d.bars(items, lo=0.0, hi=1.0, center=0.5, rowh=11.6, label_w=210, fmt='{:.3f}')
    d.y += 4
    d.h2('Raw versus length-adjusted', 'metrics residualised on log total residues, then rescored')
    items = []
    for g in GBF['Length-adjusted']:
        base_g = g.replace(' (adj.)', '')
        items.append((SHORT.get(base_g, base_g), auc_of(base_g, 'iptm'), auc_of(g, 'iptm')))
    d.dumbbell(items, lo=0.0, hi=1.0, rowh=12.6, label_w=182, a_lab='raw ipTM AUC',
               b_lab='length-adjusted', fmt='{:.3f}')
    d.takeaway('Adjusting for length moves the gold-standard test from 0.12 to '
               f"{auc_of('3. positives vs gold negatives (adj.)', 'iptm'):.2f} and barely touches the "
               'headline test. So the dramatic inversion was length; the null result was not. This '
               'is a sensitivity analysis, not a pre-registered test.')


def p_leakage(d):
    d.newpage('confound: information')
    d.h1('Confound 2: what the model was given',
         'Two inputs differ systematically between the arms of this benchmark, and neither was a '
         'design decision.')
    d.h2('Structure templates were on for every chain', 'server default; the batch files never set it')
    d.para('The uploaded batch files are dialect version 1 and contain no template switch. The '
           'server upgraded them to version 3 and set useStructureTemplate = true on all 95 chain '
           'entries, then returned 2 to 4 template hits per chain for all 60 jobs. For the matched '
           'internal negatives that matters directly: 6es7 - a structure containing both ACTR and '
           'NCBD, since it templates both chains of P06 - is also a template for N26 (ACTR x ACTR) '
           'and N27 (NCBD x NCBD), and 1fos, which templates both chains of P07, is a template in '
           'N15, N16, N17 and N28. AF3 applies templates within a chain, so this is the bound '
           'conformation of each partner rather than the interface itself - but it is exactly the '
           'help a zipper needs to dimerise.', size=7.9, color=INK2)
    d.y += 12
    rows = []
    for j in ['P06', 'N26', 'N27', 'P07', 'N28', 'N15', 'N17', 'P21']:
        mt = D['msa_templates'][j]; t = TBL[j]
        tmpl = '; '.join(f"{k}: {', '.join(v)}" for k, v in sorted(mt['templates'].items()))
        msa = '  '.join(f"{k.split('/')[0]} {v:,}" for k, v in sorted(mt['msa'].items())
                        if k.endswith('unpaired'))
        rows.append({'j': f"{j}  {t['chain_A']} x {t['chain_B']}", 'lab': t['label'][:8],
                     'tmpl': tmpl, 'msa': msa, 'iptm': t['iptm']})
    d.table([{'title': 'job', 'key': 'j', 'w': 158},
             {'title': 'label', 'key': 'lab', 'w': 44},
             {'title': 'template hits per chain (PDB entry)', 'key': 'tmpl', 'w': 200, 'font': 'N'},
             {'title': 'unpaired MSA depth', 'key': 'msa', 'w': 74, 'font': 'N'},
             {'title': 'ipTM', 'key': 'iptm', 'w': 35, 'align': 'r', 'fmt': lambda v, r: f'{v:.2f}',
              'font': 'SB'}], rows, rowh=11.6)
    d.y += 12
    d.h2('The scrambles really do lose their MSA', 'confirming the plan\'s reason for keeping them a separate arm')
    depths = []
    for j, mt in D['msa_templates'].items():
        for k, v in mt['msa'].items():
            if not k.endswith('unpaired'):
                continue
            is_scr = TBL[j]['label'] == 'negative_scramble' and k.startswith('b')
            depths.append((is_scr, v))
    real_d = sorted(v for s, v in depths if not s)
    scr_d = sorted(v for s, v in depths if s)
    d.para(f'Unpaired MSA depth, shuffled chains: {min(scr_d)} to {max(scr_d)} sequences '
           f'(median {int(statistics.median(scr_d))}, n={len(scr_d)}). Real chains: '
           f'{min(real_d):,} to {max(real_d):,} (median {int(statistics.median(real_d)):,}, '
           f'n={len(real_d)}). A shuffled sequence has no homologs, so low confidence there follows '
           'from missing evolutionary information rather than from a missing interface - which is '
           'why test 4 is reported apart and cannot carry the headline.', size=7.9, color=INK2)
    d.y += 12
    d.takeaway('Both confounds push the same way: they make the memorised, foldable, short, '
               'template-rich pairs look good and the long fuzzy ones look bad. The single cheapest '
               'next experiment is to re-run a handful of jobs with useStructureTemplate = false and '
               'see how much of the control signal survives.')


def p_scramble(d):
    d.newpage('scramble arm')
    d.h1('The scramble arm, paired',
         'Each scramble shares chain A with its parent positive and differs only in a shuffled chain '
         'B, so the natural analysis is paired. Eight of the ten have run; S20 and S21 are pending.')
    items = []
    for p in D['scramble_pairs']:
        lab = f"{p['parent']} -> {p['scramble']}" + ('  (memorised parent)' if p['memorised_parent'] else '')
        items.append((lab, p['scramble_vals']['iptm'], p['parent_vals']['iptm']))
    d.dumbbell(items, lo=0.0, hi=0.9, rowh=13, label_w=176, a_lab='scramble', b_lab='real parent',
               title='ipTM: parent positive versus its composition-matched scramble', fmt='{:.2f}')
    d.y += 4
    ss = D['scramble_summary']
    rows = [{'m': MLAB[m], 'w': f"{ss[m]['wins']} of {ss[m]['n']}", 'p': ss[m]['sign_p'],
             'md': ss[m]['median_delta'], 'auc': auc_of(PRE[4], m), 'q': q_of(PRE[4], m)}
            for m in METRICS]
    d.table([{'title': 'metric', 'key': 'm', 'w': 132},
             {'title': 'parent beats scramble', 'key': 'w', 'w': 92, 'align': 'c'},
             {'title': 'sign-test p', 'key': 'p', 'w': 58, 'align': 'r', 'fmt': lambda v, r: f'{v:.3f}'},
             {'title': 'median difference', 'key': 'md', 'w': 78, 'align': 'r',
              'fmt': lambda v, r: f'{v:+.3f}' if abs(v) < 10 else f'{v:+.1f}'},
             {'title': 'unpaired AUC', 'key': 'auc', 'w': 60, 'align': 'r', 'fmt': lambda v, r: f3(v),
              'fill': lambda r: div_color(r['auc'], 0.5, 0.35)},
             {'title': 'BH q', 'key': 'q', 'w': 46, 'align': 'r', 'fmt': lambda v, r: fp(v)},
             ], rows, rowh=11.6)
    d.y += 10
    d.takeaway('With 8 pairs the paired test has almost no power: ipTM favours the real parent in 6 '
               'of 8 (sign-test p = 0.29). There is no evidence either way yet on whether the little '
               'signal AF3 has depends on evolutionary information - finishing S20 and S21 makes it '
               '10 pairs, which is still thin.')


def p_local_pae(d):
    L = D['local_pae']
    LC = {(c['group'], c['metric']): c for c in L['cells']}
    def lcell(g, m): return LC.get((g, m), {})
    def lverd(g, m):
        c = lcell(g, m)
        return None if 'delta' not in c else ('PASS' if c.get('beats_baseline') else 'FAIL')
    def any_auc(g, m): return (lcell(g, m) if m in ('ipsae', 'lis') else cell(g, m)).get('auc')
    cut = L['cutoffs_angstrom']
    d.newpage('exploratory: local PAE')
    d.local_page = d.page
    d.h1('Local PAE: ipSAE and LIS',
         'Added after the stage-1 numbers had been seen, so exploratory and outside the pre-registered '
         'panel. Both use only the confident interchain pairs - ipSAE those under '
         f"{cut['ipsae']:.0f} A (Dunbrack 2025), LIS those under {cut['lis']:.0f} A (Kim et al. 2024) - so "
         'an unconfident flank cannot dilute a real interface the way it dilutes ipTM. Same top-ranked '
         'sample, same composition comparator.')

    d.h2('The pre-registered tests, rescored', 'delta-AUC against the same composition baseline')
    SH = {PRE[0]: 'Test 1', PRE[1]: 'Test 2a', PRE[2]: 'Test 2b', PRE[3]: 'Test 3', PRE[4]: 'Test 4',
          'Homotypic only': 'Homotypic arm', 'Heterotypic only': 'Heterotypic arm'}
    rows = []
    for g in PRE + ARMS:
        ci, cl = lcell(g, 'ipsae'), lcell(g, 'lis')
        rows.append({'g': SH[g], 'n': f"{ci['n_pos']}/{ci['n_neg']}", 'iptm': auc_of(g, 'iptm'),
                     'ia': ci.get('auc'), 'id': ci.get('delta'), 'idci': ci.get('delta_ci'),
                     'iv': lverd(g, 'ipsae'), 'la': cl.get('auc'), 'ld': cl.get('delta'),
                     'ldci': cl.get('delta_ci'), 'lv': lverd(g, 'lis')})
    red = lambda k: (lambda r: CRIT if (r[k] or 0) < 0 else INK)
    d.table([
        {'title': 'group', 'key': 'g', 'w': 76},
        {'title': 'n', 'key': 'n', 'w': 34, 'align': 'c'},
        {'title': 'ipTM', 'key': 'iptm', 'w': 36, 'align': 'r', 'fmt': lambda v, r: f3(v),
         'color': lambda r: INK2},
        {'title': 'ipSAE', 'key': 'ia', 'w': 38, 'align': 'r', 'fmt': lambda v, r: f3(v), 'font': 'SB'},
        {'title': 'delta', 'key': 'id', 'w': 38, 'align': 'r', 'fmt': lambda v, r: f3(v, True),
         'color': red('id')},
        {'title': 'delta 95% CI', 'key': 'idci', 'w': 62, 'align': 'c', 'fmt': lambda v, r: fci(v, True)},
        {'title': '', 'key': 'iv', 'w': 36, 'badge': True},
        {'title': 'LIS', 'key': 'la', 'w': 38, 'align': 'r', 'fmt': lambda v, r: f3(v), 'font': 'SB'},
        {'title': 'delta', 'key': 'ld', 'w': 38, 'align': 'r', 'fmt': lambda v, r: f3(v, True),
         'color': red('ld')},
        {'title': 'delta 95% CI', 'key': 'ldci', 'w': 62, 'align': 'c', 'fmt': lambda v, r: fci(v, True)},
        {'title': '', 'key': 'lv', 'w': 36, 'badge': True},
    ], rows, rowh=12.4)
    note = 'The ipTM column repeats page 1 for reference. Test 4 has no composition comparator.'
    t2a = lcell(PRE[1], 'lis').get('delta_ci')
    if t2a and abs(t2a[0]) < 1e-9:
        note += (' The LIS interval for test 2a starts at exactly zero, on 4 jobs against 3, which '
                 'does not clear it.')
    d.y += 7
    d.para(note + ' Every grouping, with q-values, is in appendix C.', size=7.2, color=MUTED)

    d.y += 10
    d.h2('Which PAE summary carries the signal',
         'AUC by group; rho is Spearman against total residues, all 60 jobs')
    RHO = {c['metric']: c['rho_length'] for c in D['confounds'] + L['confounds']}
    G4 = [(PRE[0], 'Test 1'), (PRE[2], 'Test 2b'), ('Homotypic only', 'Homotypic arm'),
          ('Heterotypic only', 'Heterotypic arm')]
    SUMS = [('iptm', 'every pair'), ('pae_inter_mean', 'every pair'),
            ('pae_inter_frac_lt10', 'every pair'), ('chain_pair_pae_min', 'single best pair'),
            ('ipsae', f"pairs < {cut['ipsae']:.0f} A"), ('lis', f"pairs < {cut['lis']:.0f} A")]
    rows = []
    for m, uses in SUMS:
        r = {'m': MLAB[m], 'u': uses, 'rho': RHO[m]}
        r.update({f'a{i}': any_auc(g, m) for i, (g, _) in enumerate(G4)})
        rows.append(r)
    cols = [{'title': 'metric', 'key': 'm', 'w': 110},
            {'title': 'scores', 'key': 'u', 'w': 78, 'color': lambda r: INK2}]
    for i, (_, lab) in enumerate(G4):
        cols.append({'title': lab, 'key': f'a{i}', 'w': 64, 'align': 'r', 'fmt': lambda v, r: f3(v),
                     'fill': (lambda k: lambda r: div_color(r[k], 0.5, 0.35))(f'a{i}')})
    cols.append({'title': 'rho vs length', 'key': 'rho', 'w': 60, 'align': 'r',
                 'fmt': lambda v, r: f'{v:+.2f}',
                 'color': lambda r: CRIT if abs(r['rho']) > 0.5 else INK2})
    d.table(cols, rows, rowh=12.2)

    d.y += 10
    d.h2('Where the confident pairs are', 'jobs with any interchain pair under each cutoff, by class')
    CLS = [('Positive, fuzzy heterotypic', lambda t: t['label'] == 'positive'
            and t['job_id'] not in MEM and t['pair_type'] == 'heterotypic'),
           ('Positive, fuzzy homotypic', lambda t: t['label'] == 'positive'
            and t['job_id'] not in MEM and t['pair_type'] == 'homotypic'),
           ('Positive, memorised control', lambda t: t['job_id'] in MEM),
           ('Negative, noncognate', lambda t: t['neg_type'] == 'noncognate'),
           ('Negative, cross-kingdom', lambda t: t['neg_type'] == 'cross_kingdom'),
           ('Negative, charge-repulsive', lambda t: t['neg_type'] == 'charge_repulsive'),
           ('Negative, hard zipper', lambda t: t['neg_type'] == 'hard_coiledcoil'),
           ('Negative, homotypic', lambda t: t['neg_type'] == 'homo_negative'),
           ('Scramble', lambda t: t['label'] == 'negative_scramble')]
    rows = []
    for lab, f in CLS:
        sel = [t for t in D['table'] if f(t)]
        top = max(sel, key=lambda t: t['ipsae'])
        rows.append({'c': lab, 'i0': f"{sum(t['ipsae'] > 0 for t in sel)} of {len(sel)}",
                     'im': statistics.median(t['ipsae'] for t in sel),
                     'top': f"{top['ipsae']:.3f}  {top['job_id']}" if top['ipsae'] > 0 else '--',
                     'l0': f"{sum(t['lis'] > 0 for t in sel)} of {len(sel)}",
                     'lm': statistics.median(t['lis'] for t in sel)})
    d.table([{'title': 'class', 'key': 'c', 'w': 150},
             {'title': f"any pair < {cut['ipsae']:.0f} A", 'key': 'i0', 'w': 76, 'align': 'c'},
             {'title': 'median ipSAE', 'key': 'im', 'w': 64, 'align': 'r',
              'fmt': lambda v, r: f'{v:.3f}', 'font': 'SB'},
             {'title': 'highest ipSAE', 'key': 'top', 'w': 76, 'align': 'r'},
             {'title': f"any pair < {cut['lis']:.0f} A", 'key': 'l0', 'w': 76, 'align': 'c'},
             {'title': 'median LIS', 'key': 'lm', 'w': 64, 'align': 'r',
              'fmt': lambda v, r: f'{v:.3f}', 'font': 'SB'}], rows, rowh=12.2)

    d.y += 10
    d.h2('The sharpest contrasts', 'matched internal negatives and the scramble arm')
    mi = {m['positive']: m for m in L['matched_internal']}
    iv = lambda p, j: mi[p]['vals'][j]['ipsae']
    ss = L['scramble_summary']['ipsae']
    hom_scr = [p for p in L['scramble_pairs']
               if TBL[p['parent']]['pair_type'] == 'homotypic' and p['parent_vals']['ipsae'] > 0]
    order = ('Both matched contrasts keep the right order on ipSAE and LIS'
             if all(all(m['correct'].values()) for m in L['matched_internal'])
             else 'The matched contrasts do not all keep the right order')
    d.para(f"{order}. Fos x Jun sharpens: ipSAE {iv('P07', 'P07'):.3f} against Fos x Fos "
           f"{iv('P07', 'N28'):.3f}, where ipTM gave {TBL['P07']['iptm']:.2f} against "
           f"{TBL['N28']['iptm']:.2f}. ACTR x NCBD scores {iv('P06', 'P06'):.3f} against "
           f"{iv('P06', 'N26'):.3f} and {iv('P06', 'N27'):.3f} for the two homodimers. In the scramble "
           f"arm the parent wins {ss['wins']} of {ss['n']} untied pairs (sign-test p = {ss['sign_p']:.3f}); "
           f"the other {ss['ties']} score exactly zero on both sides. The homotypic confident pairs vanish "
           'with the scramble ('
           + ', '.join(f"{p['parent']} {p['parent_vals']['ipsae']:.3f} -> {p['scramble_vals']['ipsae']:.3f}"
                       for p in hom_scr)
           + '), which still cannot tell sequence-specific self-pairing from an identical-chain effect.',
           size=7.6, color=INK2)

    hom = lcell('Homotypic only', 'ipsae')
    others = [c for c in L['cells'] if c.get('beats_baseline')
              and c['group'] != 'Mutual-folding (memorised) positives']
    under_min = all(any_auc(g, 'chain_pair_pae_min') >= max(any_auc(g, 'ipsae'), any_auc(g, 'lis'))
                    for g, _ in G4)
    het_pos = [t for t in D['table'] if t['label'] == 'positive' and t['job_id'] not in MEM
               and t['pair_type'] == 'heterotypic']
    hom_pos = [t for t in D['table'] if t['label'] == 'positive' and t['job_id'] not in MEM
               and t['pair_type'] == 'homotypic']
    d.y += 10
    d.takeaway(
        f"Scoring only the confident pairs changes no verdict. It lifts test 2b from "
        f"{auc_of(PRE[2], 'iptm'):.2f} to {lcell(PRE[2], 'ipsae')['auc']:.2f} and the homotypic arm from "
        f"{auc_of('Homotypic only', 'iptm'):.2f} to {hom['auc']:.2f} - enough that its delta interval now "
        f"includes zero ({hom['delta']:+.2f} [{hom['delta_ci'][0]:+.2f}, {hom['delta_ci'][1]:+.2f}]) where "
        f"ipTM's excluded it, not enough to beat composition - and drops the heterotypic arm from "
        f"{auc_of('Heterotypic only', 'iptm'):.2f} to {lcell('Heterotypic only', 'ipsae')['auc']:.2f}. "
        + ('Outside the memorised controls no cell beats composition' if not others
           else f'Outside the memorised controls {len(others)} cells beat composition')
        + (', and neither score out-ranks chain-pair PAE min, the crudest local summary already in the '
           'panel, in any of the four groups above. ' if under_min else '. ')
        + f"The class table is why: {sum(t['ipsae'] == 0 for t in het_pos)} of {len(het_pos)} "
        f"heterotypic fuzzy positives have no interchain pair under {cut['ipsae']:.0f} A, and the "
        f"homotypic ones top out at ipSAE {max(t['ipsae'] for t in hom_pos):.3f}. A local score can only "
        'sharpen an interface AF3 has actually built.')


def p_methods(d):
    d.newpage('methods')
    d.h1('How every number here was produced',
         'Nothing in this report is computed by hand. The statistics come from the project\'s own '
         'functions so the report cannot drift from the pipeline.')
    blocks = [
        ('Data', '60 of 62 AlphaFold Server jobs, five samples each, downloaded 14 September 2026. '
         'Per job the reported value is the sample with the highest ranking_score, matching '
         'parse_af3_results.py. S20 and S21 (scrambles of the two memorised homodimers) are not yet run.'),
        ('AUC', 'Mann-Whitney AUC with tie correction, on the metric oriented so higher means more '
         'interaction. 0.5 is chance; below 0.5 means the negatives outrank the positives.'),
        ('Confidence intervals', 'Percentile bootstrap, 10,000 resamples, resampling positives and '
         'negatives independently.'),
        ('p and q', 'Two-sided permutation test on the label assignment, 10,000 permutations, then '
         'Benjamini-Hochberg across the 16-metric panel within each group. Groups with fewer than 3 '
         'on either side are reported as descriptive only. ipSAE and LIS are corrected against the '
         'same 16 plus their own 2, so the panel q-values are untouched.'),
        ('Composition baseline', 'Leave-one-out logistic regression on 24 cheap pair features '
         '(NCPR, FCR, aromatic and G/S content, hydropathy, length, charge complementarity, '
         'homotypic flag - means and absolute differences per pair). Fitted on the 52 real-sequence '
         'jobs for pooled groups, and refitted inside each arm for the homotypic and heterotypic '
         'rows. Raw LOO AUC: pooled 0.700, homotypic 0.851, heterotypic 0.471 (0.529 two-sided).'),
        ('Delta-AUC', 'Paired bootstrap on AUC(AF3) - AUC(baseline), 5,000 resamples, resampling the '
         'same jobs for both score vectors. A test passes only if the 95% interval excludes zero. '
         'This is the pre-registered primary criterion; beating 0.5 is not the bar.'),
        ('Local PAE scores', 'Both come from the full PAE matrix of the same top-ranked sample. '
         f"ipSAE follows ipsae.py (Dunbrack 2025), d0res variant: for each residue, partners with PAE "
         f"< {D['local_pae']['cutoffs_angstrom']['ipsae']:.0f} A are kept, d0 is set from their number "
         '(floored at 26 residues and 1 A), 1/(1+(PAE/d0)^2) is averaged over them, and the score is '
         'the best residue in the better direction. LIS follows Kim et al. 2024: the mean of '
         f"1 - PAE/{D['local_pae']['cutoffs_angstrom']['lis']:.0f} over interchain pairs with PAE below "
         f"{D['local_pae']['cutoffs_angstrom']['lis']:.0f} A, per direction, averaged over both. Both "
         'match a residue-by-residue transcription of the reference code on all 60 jobs to 1e-16.'),
        ('Pre-registered versus exploratory', 'Pre-registered: the five subsets on page 1, the '
         'homotypic/heterotypic split, the 16-metric panel with BH correction, and the length '
         'confound flag. Exploratory, and labelled as such: the negative-subtype, evidence-grade and '
         'positive-subgroup groupings, the length-adjusted sensitivity analysis, the paired '
         f'scramble test, and the two local PAE scores on page {d.local_page} and in appendix C, '
         'which were added after the stage-1 numbers had been seen.'),
        ('Known caveats carried from the plan', 'Structure templates were on for all 95 chain '
         'entries (server default). 59 jobs ran on seed 1 and N24 on a server-assigned seed. Five '
         'samples share one seed, so spread is sample spread, not seed replication. 21 of 31 '
         'negatives rest on inferred or absence-only evidence. Length is not matched by design. Two '
         'negatives (H1.0 alone, ERD10 alone) are graded uncertain and should be confirmed or dropped. '
         'has_clash appears in the plan\'s metric list but not in the parser\'s panel; it is true in '
         'some sample of P13 and P17 only.'),
        ('Files', 'af3_report_data.json holds every number in this report; af3_all_statistics.csv is '
         'the full table - 496 panel cells (31 groups x 16 metrics), then 62 exploratory rows with '
         'metric ipsae or lis; af3_scored.csv is the per-job scored table written by '
         'parse_af3_results.py, with ipsae and lis as its last columns.'),
    ]
    for title, body in blocks:
        d.text(ML, d.y, title, 'SB', 8.2, ACCENT); d.y += 11
        d.para(body, size=7.9, color=INK2, width=CW - 8, x=ML)
        d.y += 9


def p_appendix_jobs(d):
    per_page = 44
    rows_all = sorted(D['table'], key=lambda t: t['job_id'])
    cols = [{'title': 'ID', 'key': 'job_id', 'w': 24, 'font': 'NB'},
            {'title': 'chains', 'key': 'ch', 'w': 106, 'font': 'N'},
            {'title': 'label', 'key': 'lab', 'w': 48, 'font': 'N'},
            {'title': 'class', 'key': 'cls', 'w': 58, 'font': 'N'},
            {'title': 'res', 'key': 'len', 'w': 24, 'align': 'r', 'font': 'N'},
            {'title': 'ipTM', 'key': 'iptm', 'w': 28, 'align': 'r', 'font': 'NB',
             'fmt': lambda v, r: f'{v:.2f}'},
            {'title': 'sd', 'key': 'sd', 'w': 24, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: f'{v:.3f}'},
            {'title': 'pTM', 'key': 'ptm', 'w': 26, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: f'{v:.2f}'},
            {'title': 'PAE', 'key': 'pae', 'w': 28, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: f'{v:.1f}'},
            {'title': 'pLDDT', 'key': 'plddt', 'w': 30, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: f'{v:.1f}'},
            {'title': 'C sum', 'key': 'cs', 'w': 30, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: f'{v:.0f}'},
            {'title': 'C>.5', 'key': 'nc', 'w': 26, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: f'{v:.0f}'},
            {'title': 'disord.', 'key': 'fd', 'w': 28, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: f'{v:.2f}'},
            {'title': 'compos.', 'key': 'base', 'w': 30, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: '--' if v is None else f'{v:+.2f}'}]
    for k in range(0, len(rows_all), per_page):
        d.newpage('appendix A / every job')
        if k == 0:
            d.h1('Appendix A - every job, every headline metric',
                 'All 60 scored jobs. This is the table view for every figure in the report: each '
                 'plotted dot appears here as a number.')
        else:
            d.text(ML, d.y, 'Appendix A continued', 'SB', 10, INK); d.y += 16
        rows = []
        for t in rows_all[k:k + per_page]:
            rows.append({'job_id': t['job_id'], 'ch': f"{t['chain_A']} x {t['chain_B']}",
                         'lab': t['label'].replace('negative_scramble', 'scramble'),
                         'cls': t['neg_type'] or t['binding_mode'],
                         'len': int(float(t['total_residues'])), 'iptm': t['iptm'],
                         'sd': t['iptm_sample_sd'], 'ptm': t['ptm'], 'pae': t['pae_inter_mean'],
                         'plddt': t['plddt_mean'], 'cs': t['contact_prob_sum'],
                         'nc': t['n_contacts_p50'], 'fd': t['fraction_disordered'],
                         'base': t['baseline_pooled']})
        d.table(cols, rows, size=6.8, rowh=10.4, font='N', zebra=True)
        d.y += 6
        d.para('PAE = mean interchain PAE (A). C sum = interchain contact-probability sum. C>.5 = '
               'interchain token pairs above 0.5. compos. = composition-baseline score (leave-one-out, '
               'pooled fit). sd = spread across the five samples.', size=6.8, color=MUTED)


def p_appendix_stats(d):
    cols = [{'title': 'metric', 'key': 'm', 'w': 116, 'font': 'N'},
            {'title': 'n', 'key': 'n', 'w': 40, 'align': 'c', 'font': 'N'},
            {'title': 'AUC', 'key': 'auc', 'w': 34, 'align': 'r', 'font': 'NB',
             'fmt': lambda v, r: f3(v)},
            {'title': '95% CI', 'key': 'ci', 'w': 62, 'align': 'c', 'font': 'N',
             'fmt': lambda v, r: fci(v)},
            {'title': 'perm p', 'key': 'p', 'w': 44, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: fp(v)},
            {'title': 'BH q', 'key': 'q', 'w': 44, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: fp(v)},
            {'title': 'delta', 'key': 'd', 'w': 48, 'align': 'r', 'font': 'N',
             'fmt': lambda v, r: f3(v, True)},
            {'title': 'delta 95% CI', 'key': 'dci', 'w': 64, 'align': 'c', 'font': 'N',
             'fmt': lambda v, r: fci(v, True)},
            {'title': '', 'key': 'v', 'w': 34, 'badge': True}]
    blocks = []
    for fam in FAMS:
        for g in GBF[fam]:
            rows = []
            for m in METRICS:
                c = cell(g, m)
                rows.append({'m': MLAB[m], 'n': f"{c['n_pos']}/{c['n_neg']}", 'auc': c.get('auc'),
                             'ci': c.get('ci'), 'p': c.get('p'), 'q': c.get('q'),
                             'd': c.get('delta'), 'dci': c.get('delta_ci'),
                             'v': ('PASS' if c.get('beats_baseline') else None)})
            rf = REF[g]
            head = (f"{fam}  /  {g}   ({rf['n_pos']} positive / {rf['n_neg']} negative"
                    + (f", length alone {rf['length_auc']:.3f}" if 'length_auc' in rf else '')
                    + (f", composition {rf['baseline_auc']:.3f}" if 'baseline_auc' in rf else '') + ')')
            blocks.append((head, rows))
    first = True
    i = 0
    while i < len(blocks):
        d.newpage('appendix B / all 496 statistics')
        if first:
            d.h1('Appendix B - every group, every metric',
                 '31 groupings x 16 metrics = 496 cells, each with its own bootstrap interval, '
                 'permutation p, BH-corrected q within the group, and delta against the composition '
                 'baseline (464 cells have one). A PASS badge marks a delta interval clear of zero: '
                 '15 cells qualify, 14 of them in the memorised-control block and one outside it '
                 '(fuzzy homotypic positives, contact probability sum). No pre-registered test passes.')
            first = False
        else:
            d.text(ML, d.y, 'Appendix B continued', 'SB', 10, INK); d.y += 16
        while i < len(blocks):
            head, rows = blocks[i]
            need = 14 + 10.2 * (len(rows) + 1) + 8
            if d.y + need > H - MB:
                break
            d.text(ML, d.y, d.fit(head, 'SB', 7.4, CW), 'SB', 7.4, ACCENT); d.y += 11
            d.table(cols, rows, size=6.7, rowh=10.2, font='N')
            d.y += 9
            i += 1


def p_appendix_local(d):
    L = D['local_pae']
    LC = {(c['group'], c['metric']): c for c in L['cells']}
    names = [m['name'] for m in L['metrics']]
    passed = [c for c in L['cells'] if c.get('beats_baseline')]
    n_mem = sum(1 for c in passed if c['group'] == 'Mutual-folding (memorised) positives')
    where = ((('both' if len(passed) == 2 else 'all') + ' in the memorised-control block')
             if passed and n_mem == len(passed) else f'{n_mem} of them in the memorised-control block')
    d.newpage('appendix C / local PAE')
    d.h1('Appendix C - local PAE, every group',
         f"31 groupings x ipSAE and LIS = {len(L['cells'])} exploratory cells, in the layout of appendix "
         'B. q is Benjamini-Hochberg over the group\'s 16 panel metrics plus these two. A PASS badge '
         f'marks a delta interval clear of zero: {len(passed)} cells qualify, {where}.')
    cols = [{'title': 'group', 'key': 'g', 'w': 104, 'font': 'N'},
            {'title': 'n', 'key': 'n', 'w': 30, 'align': 'c', 'font': 'N'}]
    for m in names:
        cols += [{'title': MLAB[m], 'key': f'{m}_auc', 'w': 34, 'align': 'r', 'font': 'NB',
                  'fmt': lambda v, r: f3(v)},
                 {'title': 'q', 'key': f'{m}_q', 'w': 34, 'align': 'r', 'font': 'N',
                  'fmt': lambda v, r: 'n<3' if r.get('small') else fp(v)},
                 {'title': 'delta', 'key': f'{m}_d', 'w': 34, 'align': 'r', 'font': 'N',
                  'fmt': lambda v, r: f3(v, True)},
                 {'title': 'delta 95% CI', 'key': f'{m}_dci', 'w': 54, 'align': 'c', 'font': 'N',
                  'fmt': lambda v, r: fci(v, True)},
                 {'title': '', 'key': f'{m}_v', 'w': 30, 'badge': True}]
    rows = []
    for fam in FAMS:
        rows.append({'_sub': fam})
        for g in GBF[fam]:
            r = {'g': g[3:] if g.startswith('vs ') else g}
            for m in names:
                c = LC.get((g, m), {})
                r['n'], r['small'] = f"{c.get('n_pos')}/{c.get('n_neg')}", c.get('small_n')
                r.update({f'{m}_auc': c.get('auc'), f'{m}_q': c.get('q'), f'{m}_d': c.get('delta'),
                          f'{m}_dci': c.get('delta_ci'),
                          f'{m}_v': 'PASS' if c.get('beats_baseline') else None})
            rows.append(r)
    d.table(cols, rows, size=6.7, rowh=10.2, font='N')

    d.newpage('appendix C / local PAE')
    d.text(ML, d.y, 'Appendix C continued - local PAE, every job', 'SB', 10, INK); d.y += 16
    rows = [{'job_id': t['job_id'], 'ch': f"{t['chain_A']} x {t['chain_B']}",
             'lab': t['label'].replace('negative_scramble', 'scramble'),
             'cls': t['neg_type'] or t['binding_mode'], 'len': int(float(t['total_residues'])),
             'iptm': t['iptm'], 'pmin': t['chain_pair_pae_min'], 'f10': t['pae_inter_frac_lt10'],
             'ipsae': t['ipsae'], 'lis': t['lis']}
            for t in sorted(D['table'], key=lambda t: t['job_id'])]
    d.table([{'title': 'ID', 'key': 'job_id', 'w': 24, 'font': 'NB'},
             {'title': 'chains', 'key': 'ch', 'w': 118, 'font': 'N'},
             {'title': 'label', 'key': 'lab', 'w': 50, 'font': 'N'},
             {'title': 'class', 'key': 'cls', 'w': 70, 'font': 'N'},
             {'title': 'res', 'key': 'len', 'w': 26, 'align': 'r', 'font': 'N'},
             {'title': 'ipTM', 'key': 'iptm', 'w': 32, 'align': 'r', 'font': 'N',
              'fmt': lambda v, r: f'{v:.2f}'},
             {'title': 'PAE min', 'key': 'pmin', 'w': 46, 'align': 'r', 'font': 'N',
              'fmt': lambda v, r: f'{v:.1f}'},
             {'title': 'PAE < 10 A', 'key': 'f10', 'w': 50, 'align': 'r', 'font': 'N',
              'fmt': lambda v, r: f'{v:.3f}'},
             {'title': 'ipSAE', 'key': 'ipsae', 'w': 44, 'align': 'r', 'font': 'NB',
              'fmt': lambda v, r: f'{v:.3f}'},
             {'title': 'LIS', 'key': 'lis', 'w': 44, 'align': 'r', 'font': 'NB',
              'fmt': lambda v, r: f'{v:.3f}'}], rows, size=6.8, rowh=10.4, font='N', zebra=True,
            paginate={'section': 'appendix C / local PAE',
                      'title': 'Appendix C continued - local PAE, every job'})
    d.y += 6
    d.para('PAE min = chain-pair minimum interchain PAE (A). PAE < 10 A = fraction of interchain pairs '
           'under 10 A. ipSAE and LIS as defined on the methods page; 0 means no interchain pair '
           'cleared the cutoff.', size=6.8, color=MUTED)


def write_csv():
    with open(CSV_OUT, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['family', 'group', 'metric', 'n_pos', 'n_neg', 'auc', 'ci_lo', 'ci_hi',
                    'perm_p', 'bh_q', 'baseline_auc', 'delta_auc', 'delta_ci_lo', 'delta_ci_hi',
                    'beats_baseline', 'small_n', 'length_only_auc'])
        # the 496 panel cells, then the 62 exploratory local-PAE cells (metric ipsae / lis)
        for c in D['cells'] + D['local_pae']['cells']:
            rf = REF[c['group']]
            ci = c.get('ci') or [None, None]
            dci = c.get('delta_ci') or [None, None]
            w.writerow([c['family'], c['group'], c['metric'], c['n_pos'], c['n_neg'],
                        c.get('auc'), ci[0], ci[1], c.get('p'), c.get('q'), c.get('base_auc'),
                        c.get('delta'), dci[0], dci[1], c.get('beats_baseline'),
                        c.get('small_n'), rf.get('length_auc')])


def main():
    d = Doc()
    p_cover(d); p_forest(d); p_pipeline(d); p_primary(d); p_arms(d)
    p_negtypes(d); p_evidence(d); p_positives(d); p_secondary(d)
    p_length(d); p_leakage(d); p_scramble(d); p_local_pae(d); p_methods(d)
    p_appendix_jobs(d); p_appendix_stats(d); p_appendix_local(d)
    d.check()
    d.c.save()
    write_csv()
    print(f'wrote {OUT}  ({d.page} pages)')
    print(f'wrote {CSV_OUT}')
    if d.warns:
        print('LAYOUT WARNINGS:')
        for w_ in d.warns:
            print('  ', w_)
    else:
        print('no layout overflow detected')


if __name__ == '__main__':
    main()
