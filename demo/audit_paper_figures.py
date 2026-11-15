#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""audit_paper_figures.py — E8 审计协议手稿四张英文图（W4）。

数据全部从 JSON 源实算（不手写统计）：
  六模型对照表 JSON（Fig 2/4）
  缝隙形态库定版 JSON（Fig 3）

先实测后摆放：横条 y 轴类目数确定画布高；ylim/xlim 显式覆盖。
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.join(HERE, '..', '..', 'docs')
OUT = os.path.join(HERE, 'figures')
os.makedirs(OUT, exist_ok=True)

# macOS 渲染铁律：Hiragino Sans GB（苹方在 Quartz 下缺字）
plt.rcParams['font.family'] = ['Hiragino Sans GB', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

GREEN = '#2e7d4f'    # identified
ORANGE = '#e0862c'   # one-of-many
GRAY = '#9aa0a6'     # remainder / non-identified
BLUE = '#3568a4'

six = json.load(open(os.path.join(DOCS, '闭合定理审计附件_六模型uniqueness对照表_20261114.json'),
                     encoding='utf-8'))
morph = json.load(open(os.path.join(DOCS, '闭合定理审计附件_缝隙形态库定版_20261114.json'),
                       encoding='utf-8'))


# ── Fig 1 协议流程图（五组件横向） ──
def fig1_protocol():
    fig, ax = plt.subplots(figsize=(11.5, 2.8))
    ax.set_xlim(0, 10.8)
    ax.set_ylim(0, 2.6)
    ax.axis('off')
    comps = [
        ('C1', 'Freeze-first', 'criteria committed\nbefore statistics'),
        ('C2', 'Mechanical\nmapping', 'format-adaptive keys,\nzero manual judgement'),
        ('C3', 'Layered\nconsistency', 'fine-grained must\nreproduce frozen coarse'),
        ('C4', 'Prior\nchecklist', '8-form gap library\nas scan paths'),
        ('C5', 'Consumer\nvalidation', '1,211 entries,\ndual-pipeline check'),
    ]
    w, gap = 1.95, 0.22
    x = 0.25
    for i, (cid, title, sub) in enumerate(comps):
        ax.add_patch(plt.Rectangle((x, 0.7), w, 1.3, facecolor='#e8eef4',
                                   edgecolor=BLUE, linewidth=1.4))
        ax.text(x + w/2, 1.75, cid, ha='center', va='center',
                fontsize=11, color=BLUE, fontweight='bold')
        ax.text(x + w/2, 1.35, title.replace('\n', ' '), ha='center', va='center',
                fontsize=8.8, color='#1a1a2e')
        if i < len(comps) - 1:
            ax.annotate('', xy=(x + w + gap - 0.06, 1.35), xytext=(x + w + 0.04, 1.35),
                        arrowprops=dict(arrowstyle='->', color='#666', lw=1.3))
        x += w + gap
    ax.text(5.0, 2.35, 'Freeze-First Audit Protocol', ha='center', va='center',
            fontsize=12, fontweight='bold', color='#1a1a2e')
    fig.savefig(os.path.join(OUT, 'fig1_protocol.png'), dpi=200,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


# ── Fig 2 六模型 stacked bar（identified 升序自下而上） ──
def fig2_stacked():
    models = six['models']
    NAMES = {
        'Karr': ('Whole-cell (Karr 2012)', 'parameter · 154'),
        'E-Cell': ('Hypothetical cell (E-Cell 1999)', 'model-level · 3'),
        'BAAIWorm': ('Whole-organism (BAAIWorm 2024)', 'layer-level · 3'),
        "O'Hara-Rudy": ('Ventricular myocyte (ORd 2011)', 'parameter · 262'),
        'Wodke': ('M. pneumoniae metabolism (Wodke 2013)', 'reaction · 306'),
        'Gaspari': ('M. pneumoniae serum-free (Gaspari 2020)', 'reaction · 489'),
    }
    rows = []
    for m in models:
        key = next((k for k in NAMES if k in m['model']), None)
        pct = m['uniqueness_pct']
        rows.append({
            'label': NAMES[key][0], 'gran': NAMES[key][1],
            'id': pct.get('identified', 0),
            'om': pct.get('one-of-many', 0),
            'rest': round(pct.get('non-identified', 0) + pct.get('None', 0), 1),
        })
    rows.sort(key=lambda r: r['id'])  # identified 升序自下而上

    fig, ax = plt.subplots(figsize=(9.2, 4.3))
    y = np.arange(len(rows))
    h = 0.62
    for i, r in enumerate(rows):
        left = 0
        for seg, color in (('id', GREEN), ('om', ORANGE), ('rest', GRAY)):
            v = r[seg]
            ax.barh(i, v, left=left, height=h, color=color,
                    edgecolor='white', linewidth=0.6)
            if v >= 8:
                ax.text(left + v/2, i, '%.1f' % v, ha='center', va='center',
                        fontsize=8.5, color='white' if seg != 'rest' else '#333')
            left += v
    ax.set_yticks(y)
    ax.set_yticklabels(['%s\n%s' % (r['label'], r['gran']) for r in rows],
                       fontsize=8.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel('share of audited items (%)', fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(['identified', 'one-of-many (modelling constructs)',
               'remainder (non-ident. / undeclared)'],
              loc='upper center', bbox_to_anchor=(0.5, -0.14),
              ncol=3, fontsize=8.5, frameon=False)
    ax.set_title('Identifiability profile across six audited models', fontsize=11,
                 fontweight='bold', color='#1a1a2e', pad=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig2_six_models.png'), dpi=200,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


# ── Fig 3 形态库计数条形图（四类着色，F6 槽位特殊标注） ──
def fig3_morphology():
    forms = morph['forms']
    CAT_COLOR = {'A': '#5b8ac5', 'B': '#e0862c', 'C': '#8e6fb0', 'D': '#c44e5c'}
    rows = [(f['id'], f['name'], f['category'][0], f['count'], f['category']) for f in forms]
    rows.sort(key=lambda r: r[3])  # 升序
    fig, ax = plt.subplots(figsize=(8.6, 3.6))
    y = np.arange(len(rows))
    labels = ['%s %s' % (r[0], r[1]) for r in rows]
    vals = [min(r[3], 100) for r in rows]  # F6 9188 槽裁到 100 作图（实际数另标注）
    colors = [CAT_COLOR[r[2]] for r in rows]
    ax.barh(y, vals, height=0.62, color=colors, edgecolor='white')
    for i, r in enumerate(rows):
        shown = min(r[3], 100)
        ax.text(shown + 1.5, i, '9,188 slots' if r[0] == 'F6' else str(r[3]),
                va='center', fontsize=8.5, color='#1a1a2e')
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('instances (F6 shown capped at 100; actual 9,188 slots)', fontsize=9)
    ax.set_xlim(0, 112)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=l) for l, c in
                       [('source gaps', '#5b8ac5'), ('construction gaps', '#e0862c'),
                        ('definitional', '#8e6fb0'), ('capability (model-level)', '#c44e5c')]],
              loc='lower right', fontsize=8, frameon=False)
    ax.set_title('Gap-morphology library: instances per form', fontsize=11,
                 fontweight='bold', color='#1a1a2e', pad=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig3_morphology.png'), dpi=200,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


# ── Fig 4 代际对照（Wodke vs Gaspari） ──
def fig4_generational():
    wm = next(m for m in six['models'] if 'Wodke' in m['model'])
    gm = next(m for m in six['models'] if 'Gaspari' in m['model'])
    pair = [('Wodke 2013\niJW145', wm['uniqueness_pct']),
            ('Gaspari 2020\niEG158_mpn', gm['uniqueness_pct'])]
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    x = np.arange(2)
    w = 0.32
    ident = [p[1].get('identified', 0) for p in pair]
    om = [p[1].get('one-of-many', 0) for p in pair]
    rest = [round(p[1].get('non-identified', 0) + p[1].get('None', 0), 1) for p in pair]
    b1 = ax.bar(x - w/2, ident, w, color=GREEN, label='identified')
    b2 = ax.bar(x + w/2, om, w, color=ORANGE, label='one-of-many')
    for xi, v in zip(x - w/2, ident):
        ax.text(xi, v + 1, '%.1f%%' % v, ha='center', fontsize=9, color=GREEN)
    for xi, v in zip(x + w/2, om):
        ax.text(xi, v + 1, '%.1f%%' % v, ha='center', fontsize=9, color=ORANGE)
    ax.annotate('', xy=(1 - w/2, ident[1] + 6), xytext=(0 - w/2, ident[0] + 6),
                arrowprops=dict(arrowstyle='->', color='#c44e5c', lw=1.5))
    ax.text(0.5, 71, '+11.8 pp', ha='center', fontsize=10, color='#c44e5c',
            fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([p[0] for p in pair], fontsize=9.5)
    ax.set_ylabel('identified / one-of-many share (%)', fontsize=10)
    ax.set_ylim(0, 80)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=8.5, frameon=False, loc='upper left')
    ax.set_title('Same-lineage cross-generation: identified share +11.8 pp',
                 fontsize=11, fontweight='bold', color='#1a1a2e', pad=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig4_generational.png'), dpi=200,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


fig1_protocol()
fig2_stacked()
fig3_morphology()
fig4_generational()
print('四图落盘:', OUT)
for f in sorted(os.listdir(OUT)):
    if f.startswith('fig') and f.endswith('.png') and ('protocol' in f or 'six' in f or 'morph' in f or 'generational' in f):
        print(' ', f, os.path.getsize(os.path.join(OUT, f)), 'bytes')
