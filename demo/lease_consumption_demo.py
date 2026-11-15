#!/usr/bin/env python3
"""lease_consumption_demo.py — 审计线逐条产物接入消费端演示
（消费端深化：Karr 154/ORd 262/Wodke 306/Gaspari 489 = 1,211 条
逐条判定转 census lease 条目，接入 filter_by_lease/lease_summary。

数据源（六模型审计线产物）：
  docs/闭合定理审计附件_三模型uniqueness打标_20261113.json（Karr tagged）
  docs/闭合定理审计附件_ORd2011参数冻结约束清单_20261113.json（param_level）
  docs/闭合定理审计附件_Wodke2013参数冻结约束清单_20261113.json（reaction_level）
  docs/闭合定理审计附件_Gaspari2020参数冻结约束清单_20261114.json（tagged）

来源类型映射（species_attribution 原文 → LEASE_SOURCE_TYPES）：
  mgen-explicit / own_species-explicit → own_species
  ecoli/gallisepticum/bsubtilis/mcapricolum-explicit → cross_species
  literature-unattributed → unattributed
  fit/set / fit_set → fit_set
  universal-constant → unattributed（物理常数不属物种维度，保守）

lease_axis 统一 modeling（审计租约的标准轴：模型改版即失效）；
原判定规则保留在 note（R-* 规则可追溯）。
"""
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import constraint_census_schema as S

DOCS = os.path.join(_HERE, '..', '..', '..', 'docs')

SPECIES_MAP = {
    'mgen-explicit': 'own_species',
    'own_species-explicit': 'own_species',
    'own_species': 'own_species',
    'ecoli-explicit': 'cross_species',
    'ecoli-inferred': 'cross_species',
    'gallisepticum-explicit': 'cross_species',
    'bsubtilis-explicit': 'cross_species',
    'mcapricolum-explicit': 'cross_species',
    'literature-unattributed': 'unattributed',
    'fit/set': 'fit_set',
    'fit_set': 'fit_set',
    'universal-constant': 'unattributed',
}


def _load(name):
    p = os.path.join(DOCS, name)
    with open(p, encoding='utf-8') as f:
        return json.load(f)


def import_entries():
    """四份逐条产物 → census lease 条目。"""
    entries = []

    # 1) Karr 154（三模型打标 JSON）
    tag = _load('闭合定理审计附件_三模型uniqueness打标_20261113.json')
    for t in tag['karr_2012']['tagged']:
        st = SPECIES_MAP.get(t['species_attribution'], 'unattributed')
        entries.append(S.entry(
            'karr:%s' % t['id'], 'cell', 'other',
            writer=t['name'] or t['id'],
            reader='Karr 2012 全细胞模型',
            source='§9.12 Karr 逐参数审计',
            lease=S.lease_field(st, 'modeling', 'explicit',
                                uniqueness=t['uniqueness'],
                                note='rule=%s' % t['rule'])))

    # 2) ORd 262（param_level）
    ord_d = _load('闭合定理审计附件_ORd2011参数冻结约束清单_20261113.json')
    for t in ord_d['param_level']['tagged']:
        st = SPECIES_MAP.get(t.get('species_attribution', ''), 'unattributed')
        entries.append(S.entry(
            'ord:%s' % t['param'], 'cell', 'other',
            writer=t['param'],
            reader="O'Hara-Rudy 2011（CiPA 变体）",
            source='§9.14 ORd 逐参数审计',
            lease=S.lease_field(st, 'modeling', 'explicit',
                                uniqueness=t['uniqueness'],
                                note='rule=%s' % t['rule'])))

    # 3) Wodke 306（reaction_level）
    w = _load('闭合定理审计附件_Wodke2013参数冻结约束清单_20261113.json')
    for t in w['reaction_level']['tagged']:
        st = SPECIES_MAP.get(t.get('source_type', ''), 'unattributed')
        entries.append(S.entry(
            'wodke:%s' % t['rid'], 'cell', 'other',
            writer=t.get('name') or t['rid'],
            reader='Wodke 2013 iJW145',
            source='§9.15 Wodke 逐反应审计',
            lease=S.lease_field(st, 'modeling', 'explicit',
                                uniqueness=t['uniqueness'],
                                note='rule=%s' % t['rule'])))

    # 4) Gaspari 489（tagged）
    g = _load('闭合定理审计附件_Gaspari2020参数冻结约束清单_20261114.json')
    for t in g['tagged']:
        st = SPECIES_MAP.get(t.get('source_type', ''), 'unattributed')
        entries.append(S.entry(
            'gaspari:%s' % t['rid'], 'cell', 'other',
            writer=t.get('name') or t['rid'],
            reader='Gaspari 2020 iEG158_mpn',
            source='§9.17 Gaspari 逐反应审计',
            lease=S.lease_field(st, 'modeling', 'explicit',
                                uniqueness=t['uniqueness'],
                                note='rule=%s' % t['rule'])))

    return entries


def model_of(entry_id):
    return entry_id.split(':', 1)[0]


def main():
    entries = import_entries()
    print('导入条目: %d（预期 1211）' % len(entries))
    assert len(entries) == 1211, '导入数不符'

    # 全量校验
    bad = []
    for e in entries:
        ok, err = S.validate_lease(e['lease'])
        if not ok:
            bad.append((e['entry_id'], err))
    print('validate_lease 全过: %s（bad=%d）' % (not bad, len(bad)))
    assert not bad

    # 演示 1：千条级 summary
    sm = S.lease_summary(entries)
    print('\n[演示 1] lease_summary（1211 条）:')
    print('  by_uniqueness:', sm['by_uniqueness'])
    print('  by_source_type:', sm['by_source_type'])
    print('  uniqueness_unset:', sm['uniqueness_unset'])

    # 演示 2：跨模型 identified 占比（与对照表一致性核对）
    print('\n[演示 2] 跨模型 identified 占比（filter+分组）:')
    pct = {}
    for pre in ('karr', 'ord', 'wodke', 'gaspari'):
        sub = [e for e in entries if model_of(e['entry_id']) == pre]
        ident = len(S.filter_by_lease(sub, uniqueness='identified'))
        pct[pre] = round(ident * 100.0 / len(sub), 1)
        print('  %-8s n=%-4d identified=%5.1f%%' % (pre, len(sub), pct[pre]))

    # 演示 3：one-of-many 全景（建模构造/拟合依赖全景查询）
    om = S.filter_by_lease(entries, uniqueness='one-of-many')
    om_by_model = {}
    for e in om:
        m = model_of(e['entry_id'])
        om_by_model[m] = om_by_model.get(m, 0) + 1
    print('\n[演示 3] uniqueness=one-of-many 全景 %d 条，按模型: %s'
          % (len(om), om_by_model))

    # 演示 4：AND 组合（本物种实验锚定——六模型质量核心带）
    core = S.filter_by_lease(entries, source_type='own_species',
                             uniqueness='identified')
    core_by_model = {}
    for e in core:
        m = model_of(e['entry_id'])
        core_by_model[m] = core_by_model.get(m, 0) + 1
    print('\n[演示 4] own_species ∧ identified %d 条，按模型: %s'
          % (len(core), core_by_model))

    # 演示 5：不可辨识全景
    ni = S.filter_by_lease(entries, uniqueness='non-identified')
    print('\n[演示 5] uniqueness=non-identified %d 条: %s'
          % (len(ni), [e['entry_id'] for e in ni[:6]]))

    # 报告落盘
    report = {
        'produced_by': 'lease_consumption_demo.py',
        'entries': len(entries),
        'validate_lease_fail': len(bad),
        'summary': {k: sm[k] for k in ('by_uniqueness', 'by_source_type',
                                       'uniqueness_unset')},
        'identified_pct_by_model': pct,
        'one_of_many_by_model': om_by_model,
        'own_species_identified_by_model': core_by_model,
        'non_identified': [e['entry_id'] for e in ni],
    }
    out = os.path.join(_HERE, 'lease_consumption_report.json')
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=1)
    print('\n报告落盘:', out)


if __name__ == '__main__':
    main()
