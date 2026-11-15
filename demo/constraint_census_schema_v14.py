#!/usr/bin/env python3
"""constraint_census_schema.py — 约束条目语法 schema v1.4
（v1.0 §14.7 判据一落地；v1.2 分子层几何字段+taxon；v1.3 租约
字段 lease——§9.8.4 提案级；v1.4 lease 新增 uniqueness 维度——
§9.13 三模型租约类型学落地的字段化）
+ 普查编译器 + 孤儿/强读取查询原型。

字段表（v1.0，17 项）——每条约束条目：
  identity  : entry_id, layer, constraint_type, schema_version
  roles     : writer（书写者）, reader（阅读者, None=孤儿条目）,
              directed（分子层对称互补矩阵产出=False）
  scope     : boundary（隔室归属, None=待录）
  dynamics  : setpoint（设定点值）, drift_rate（漂移率）
  accounting: cost（维护代价账目 {energy_J_per_event, note}）
  readable  : epsilon（可读性; 分子层操作定义 = pKd = -log10(kd_nm)）
  molecular : complement_score, interface_patch, kd_estimate_nm,
              kd_extrapolated（互补分数超出校准域 [0.60,0.95] 时 True）
  lease     : 租约条款（v1.3 新增，提案级——参数槽模型的租约
              条款字段化，判据反向提炼自 Karr 2012 逐参数审计；
              v1.4 新增 uniqueness 维度——§9.13 三模型租约类型学
              落地的字段化）：
              source_type（own_species/cross_species/fit_set/
              unattributed 三态+保守态）+ lease_axis（失效轴
              evolutionary/environmental/modeling）+ evidence
              （判据强度 explicit/inferred/unattributed）+
              uniqueness（可辨识性 identified/non-identified/
              one-of-many——§9.13.6 类型学三条新增提案）+
              v_drift/B_basin（ε-窗口租约参数，未定标 None）；
              经 lease_field() 构造、validate_lease() 校验
  provenance: source（产出出处）

查询原型（§14.7 判据二起步）：
  orphan_query        谁没有阅读者（"无读者条目"试查）
  strongest_readings  按读取强度排序（kd 升序 = 阅读信道最强）
  summary             层/类型分布

普查编译（§十五接线）：compile_molecular_layer() 读 surfaceome_matrix.json
（91 边，毫秒级扫描产物），逐边生成分子层 binding 条目——几何三字段
（complement_score/interface_patch/kd_estimate_nm）由扫描在编译时填充。
kd 映射锚点复制自 aging_surface_bridge.complement_to_kd（校准域
[0.60,0.95]；矩阵分数域 [0.04,0.20] 低于下锚 → kd_extrapolated=True）。

状态：判据一落地（schema 形式化）；判据二起步（22 条件目录+分子矩阵
已接入，GO/KEGG 待查）。
"""
import json
import math
import os
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA_VERSION = '1.4'
LAYERS = ('molecular', 'cell', 'organ', 'clinical', 'ecosystem')
CONSTRAINT_TYPES = ('binding', 'modification', 'boundary', 'setpoint',
                    'transport', 'degradation', 'other')
REQUIRED = ('entry_id', 'layer', 'constraint_type', 'writer', 'reader',
            'source')
GEOMETRY_FIELDS = ('complement_score', 'interface_patch', 'kd_estimate_nm')

# kd 映射锚点（复制自 aging_surface_bridge.complement_to_kd，校准域 [0.60, 0.95]）
_KD_ANCHORS = ((0.95, 0.044), (0.85, 0.20), (0.60, 1000.0))


def complement_to_kd(comp):
    """互补分数 → Kd(nM) 对数线性插值（锚点同 aging_surface_bridge）。"""
    (c_hi, k_hi), (c_lo, k_lo) = _KD_ANCHORS[0], _KD_ANCHORS[-1]
    slope = (math.log10(k_hi) - math.log10(k_lo)) / (c_hi - c_lo)
    return 10 ** (math.log10(k_lo) + slope * (comp - c_lo))


def entry(entry_id, layer, constraint_type, writer, reader, source,
          directed=False, boundary=None, setpoint=None, drift_rate=None,
          cost=None, epsilon=None, complement_score=None,
          interface_patch=None, kd_estimate_nm=None,
          kd_extrapolated=False, subtype=None,
          writer_unknown=False, reader_unknown=False,
          taxon=None, lease=None):
    """构造一条约束条目（schema 字段表全量出现——None=该字段本条目缺省）。"""
    return {
        'entry_id': entry_id, 'layer': layer,
        'constraint_type': constraint_type,
        'writer': writer, 'reader': reader, 'directed': directed,
        'boundary': boundary, 'setpoint': setpoint,
        'drift_rate': drift_rate, 'cost': cost, 'epsilon': epsilon,
        'complement_score': complement_score,
        'interface_patch': interface_patch,
        'kd_estimate_nm': kd_estimate_nm,
        'kd_extrapolated': kd_extrapolated, 'subtype': subtype,
        'writer_unknown': writer_unknown,
        'reader_unknown': reader_unknown, 'taxon': taxon,
        'lease': lease,
        'source': source, 'schema_version': SCHEMA_VERSION,
    }


def lease_field(source_type=None, lease_axis=None, evidence=None,
                uniqueness=None, v_drift=None, B_basin=None, note=None):
    """租约字段构造器（v1.3 新增，提案级；v1.4 新增 uniqueness
    维度——§9.13 三模型租约类型学落地的字段化）——参数槽模型的租约条款。

    判据反向提炼自 Karr 2012 逐参数审计（§九/思想站点 20261112）：
    三态判定 + 失效轴 + 判据强度 + (v, B) 条目化钩子。

    source_type 三态+保守态（约束参数的来源类型）：
      own_species    本物种测量——判据：注释/正文明示测量物种=
                     条目所属谱系；
      cross_species  跨谱系外推——判据：明示他物种测量，或同源
                     映射表（如 S3Y 的 E. coli Homolog 列）；
      fit_set        拟合/初始化/建模选择——无实验约束
                     （Karr 表 EC=N 或注释 "fit by sensitivity"）；
      unattributed   引文参数无物种声明——保守不判（审计纪律：
                     不猜；待引文原文核查后转前三态）。
    lease_axis 失效轴（§18.7.2 (v, B) 条目化的轴）：
      evolutionary   谱系改变即失效（源物种的选择漂移改写参数）；
      environmental  环境改变即失效（培养基/生理条件，如 LB/M9）；
      modeling       模型结构改变即失效（拟合常数随模型而改）。
    evidence 判据强度（审计可追溯性三档）：
      explicit     注释/正文显式声明；
      inferred     文档内互证（同引文在文档他处被明示，如 Feist
                   2007 经 Parameter_0092）；
      unattributed 无声明（保守态）。
    uniqueness 可辨识性（v1.4 新增——§9.13.6 类型学三条的字段
      化；三模型审计发现三种缝隙不是同一种缝隙）：
      identified   参数唯一可辨识（有充分实验约束辨识出唯一解）；
      one-of-many  解空间中的一个解（BAAIWorm Par14 作者自承
                   "one of the possible solutions"——优化有效≠
                   参数唯一可辨识）；
      non-identified 参数不可辨识（来源不明或搜索空间自由度太大）；
      None         未声明（缺口可见）。
    v_drift / B_basin：ε-窗口租约参数（漂移率/盆地宽）——未定标
      一律 None（反虚假确定性：不填数字）；二者齐备时租约期可
      计算（T_lease = B/v 下缘穿越，§18.5 保质期可计算化）。
    note：来源说明（引文锚/审计出处）。

    返回 dict（None 字段显式保留——缺口可见而非隐没）。
    """
    return {
        'source_type': source_type, 'lease_axis': lease_axis,
        'evidence': evidence, 'uniqueness': uniqueness,
        'v_drift': v_drift, 'B_basin': B_basin, 'note': note,
    }


LEASE_SOURCE_TYPES = ('own_species', 'cross_species', 'fit_set',
                      'unattributed')
LEASE_AXES = ('evolutionary', 'environmental', 'modeling')
LEASE_EVIDENCE = ('explicit', 'inferred', 'unattributed')
LEASE_UNIQUENESS = ('identified', 'non-identified', 'one-of-many')


def validate_lease(lease):
    """lease 字段校验（字段可选，填了则枚举合法）。返回 (ok, errors)。"""
    if lease is None:
        return True, []
    errors = []
    st = lease.get('source_type')
    ax = lease.get('lease_axis')
    ev = lease.get('evidence')
    if st is not None and st not in LEASE_SOURCE_TYPES:
        errors.append(f'lease.source_type 非法: {st}')
    if ax is not None and ax not in LEASE_AXES:
        errors.append(f'lease.lease_axis 非法: {ax}')
    if ev is not None and ev not in LEASE_EVIDENCE:
        errors.append(f'lease.evidence 非法: {ev}')
    uq = lease.get('uniqueness')
    if uq is not None and uq not in LEASE_UNIQUENESS:
        errors.append(f'lease.uniqueness 非法: {uq}')
    return (len(errors) == 0), errors


def validate(entries):
    """schema 校验：必填字段/层与类型枚举/entry_id 唯一。返回 (ok, errors)。"""
    errors, seen = [], set()
    for e in entries:
        eid = e.get('entry_id', '<missing>')
        for f in REQUIRED:
            if e.get(f) in (None, ''):
                errors.append(f'{eid}: 必填字段缺失 {f}')
        if e.get('layer') not in LAYERS:
            errors.append(f'{eid}: layer 非法 {e.get("layer")}')
        if e.get('constraint_type') not in CONSTRAINT_TYPES:
            errors.append(f'{eid}: constraint_type 非法')
        if eid in seen:
            errors.append(f'{eid}: entry_id 重复')
        seen.add(eid)
    return (len(errors) == 0), errors


def compile_molecular_layer(matrix_path=None):
    """§十五接线：毫秒级扫描(surfaceome 91 边) → 分子层 binding 条目。"""
    mp = matrix_path or os.path.join(HERE, 'surfaceome_matrix.json')
    m = json.load(open(mp, encoding='utf-8'))
    entries = []
    for e in m['edges']:
        a, b, comp = e['a'], e['b'], e['score']
        kd = complement_to_kd(comp)
        entries.append(entry(
            entry_id=f'mol_bind:{a}~{b}', layer='molecular',
            constraint_type='binding', writer=a, reader=b,
            directed=False, boundary=None,
            epsilon=round(-math.log10(kd), 3),
            complement_score=comp, interface_patch=f'{a}×{b}',
            kd_estimate_nm=round(kd, 3),
            kd_extrapolated=bool(comp < 0.60),
            source='surfaceome_matrix.json(毫秒级补丁互补扫描)'))
    return entries


# ── 编目条目（真实资产：生成–清除 52 条件目录 + §3.4 边界概念）──
CURATED = [
    entry('organ_setpoint:bmd', 'organ', 'setpoint',
          writer='成骨群体活性 k_form·A',
          reader='obs_dxa_bmd 观测算子（WHO 分期判据）',
          source='curated（生成-清除 52 条件目录）',
          setpoint=1.054, epsilon=0.92,
          cost={'note': '维持机构=骨重建耦合，速率由文献锚定'}),
    entry('organ_setpoint:hba1c', 'organ', 'setpoint',
          writer='升糖驱动 k_form·A（肝/肌胰岛素抵抗轴）',
          reader='obs_hba1c 观测算子（ADA 分期判据）',
          source='curated（生成-清除 52 条件目录）',
          setpoint=8.49,
          cost={'note': '未治 T2DM 稳态（k_form 修正后校准输出）'}),
    entry('organ_setpoint:uric_acid', 'organ', 'setpoint',
          writer='尿酸生成 k_form·A',
          reader='血尿酸观测算子（达标治疗 <6 mg/dL 判据）',
          source='curated（生成-清除 52 条件目录）',
          setpoint=8.41,
          cost={'note': '别嘌醇臂稳态 2.07（treat-to-target 达标）'}),
    entry('cell_boundary:plasma_membrane', 'cell', 'boundary',
          writer='封闭网络组分（磷脂合成）',
          reader='转运与膜修复过程',
          source='curated（§3.4 第一个约束）',
          boundary='细胞膜',
          cost={'note': '形态不变量：组分持续周转，边界图案持续'}),
]


def filter_by_lease(entries, source_type=None, lease_axis=None,
                    evidence=None, uniqueness=None):
    """按租约条款过滤条目（v1.3 消费端——§9.8.4 判据 (iii)；
    v1.4 扩展 uniqueness 维度）。

    任一给定维度即过滤（AND 语义）；全 None 等价于"带 lease 的
    条目"过滤。用法示例：
      filter_by_lease(entries, source_type='cross_species')  # 全部
      # 外推条目——换谱系后哪些参数失效
      filter_by_lease(entries, lease_axis='modeling')        # 全部
      # 拟合常数——模型结构改变即失效
      filter_by_lease(entries, uniqueness='one-of-many')     # v1.4
      # 解空间之一——哪些参数只是优化产物而非唯一可辨识
      filter_by_lease(entries, source_type='fit_set',
                      uniqueness='one-of-many')              # AND 组合
    """
    out = []
    for e in entries:
        lease = e.get('lease')
        if lease is None:
            continue
        if source_type is not None and lease.get('source_type') != source_type:
            continue
        if lease_axis is not None and lease.get('lease_axis') != lease_axis:
            continue
        if evidence is not None and lease.get('evidence') != evidence:
            continue
        if uniqueness is not None and lease.get('uniqueness') != uniqueness:
            continue
        out.append(e)
    return out


def lease_summary(entries):
    """账本租约概况（消费端统计——缺口可见而非隐没；
    v1.4 加 by_uniqueness 分布与 uniqueness_unset 缺口计数）。"""
    with_lease = [e for e in entries if e.get('lease')]
    by_st, by_ax, by_ev, by_uq = {}, {}, {}, {}
    for e in with_lease:
        l = e['lease']
        st, ax, ev = l.get('source_type'), l.get('lease_axis'), l.get('evidence')
        uq = l.get('uniqueness')
        if st: by_st[st] = by_st.get(st, 0) + 1
        if ax: by_ax[ax] = by_ax.get(ax, 0) + 1
        if ev: by_ev[ev] = by_ev.get(ev, 0) + 1
        if uq: by_uq[uq] = by_uq.get(uq, 0) + 1
    return {'entries': len(entries), 'with_lease': len(with_lease),
            'without_lease': len(entries) - len(with_lease),
            'by_source_type': by_st, 'by_lease_axis': by_ax,
            'by_evidence': by_ev, 'by_uniqueness': by_uq,
            'uniqueness_unset': sum(1 for e in with_lease
                                    if e['lease'].get('uniqueness') is None),
            'v_drift_set': sum(1 for e in with_lease
                               if e['lease'].get('v_drift') is not None)}


# ── E-Cell 1999 模型级租约条目（§9.8.4 判据 (ii) 实证——判读反向
# 映射：判据全部来自 §9.11 正文实抓，validate_lease 全 PASS）──
ECELL_LEASE_ENTRIES = [
    entry('ecell:model_genes_120', 'cell', 'other',
          writer='M. genitalium 基因组注释（Fraser 1995 TIGR 测序）',
          reader='E-Cell 1999 模型的转录/翻译子模型',
          source='§9.11 E-Cell 同模板审计（Bioinformatics 15:72-84）',
          lease=lease_field('own_species', 'evolutionary', 'explicit',
                            uniqueness='identified',
                            v_drift=1.678e-06, B_basin=0.05,
                            note='127 基因中 120 个来自 M. genitalium '
                                 '基因组注释（摘要显式 "Most of the genes '
                                 'are taken from Mycoplasma genitalium"）。'
                                 '(v,B) 定标（§9.8.4 判据 i，2026-11-13）：'
                                 'v=突变率/基因组/世代 1e-3（Lee 2012 '
                                 'PMID 22991466 实测）×973 世代/年（MG '
                                 '世代 9h，Karr cellCycleLength 在案）÷'
                                 '580,070bp；B=ANI 95% 物种边界 0.05'
                                 '（Richter & Rosselló-Móra 2009 DOI '
                                 '10.1073/pnas.0906412106）；'
                                 'T=B/v≈29,800 年（中性漂移上限——寄生'
                                 '瓶颈与基因组退化选择会缩短实际租约，'
                                 '量级判定用）')),
    entry('ecell:fillin_genes_7', 'cell', 'other',
          writer='建模者功能补位（磷脂途径×4+核苷酸激酶×2+'
                 'Gln-tRNA 连接酶×1）',
          reader='E-Cell 1999 模型的膜合成与核酸回收子模型',
          source='§9.11 E-Cell 同模板审计',
          lease=lease_field('fit_set', 'modeling', 'explicit',
                            uniqueness='one-of-many',
                            note='7 个基因在 M. genitalium 中未鉴定——'
                                 '建模者为补齐功能通路而加入（建模选择，'
                                 '功能归属援引其它物种同源）')),
    entry('ecell:kinetic_params', 'cell', 'other',
          writer='KEGG→EcoCyc→文献链（部分）+建模者估计（主导）',
          reader='E-Cell 1999 全部反应规则',
          source='§9.11 E-Cell 同模板审计（正文原句 "many are '
                 'unknown. We have assigned values...by estimations"）',
          lease=lease_field('fit_set', 'modeling', 'explicit',
                            uniqueness='one-of-many',
                            note='估计主导型参数槽（正文自认）；'
                                 'Barkai-Leibler 1997 鲁棒性论证为参数'
                                 '不准辩护')),
]


# ── BAAIWorm 2024 模型级租约条目（按层打标——§9.13.6 类型学
# 第三条"多尺度模型的租约必须按层打标"的字段化兑现；v1.4
# uniqueness 从提案级转实际数据）──
BAAIWORM_LEASE_ENTRIES = [
    entry('baaiworm:structure_layer', 'cell', 'other',
          writer='C. elegans 连接组/形态实验（White 1986 系）',
          reader='BAAIWorm 2024 全模型结构骨架',
          source='§9.13 BAAIWorm 审计（Nat Comput Sci 2024, '
                 'DOI 10.1038/s43588-024-00738-w）',
          lease=lease_field('own_species', 'evolutionary', 'explicit',
                            uniqueness='identified',
                            note='结构层 own_species explicit：形态/'
                                 '连接组实验真实（302 神经元+96 肌肉）'
                                 '——§9.13.6 第三条结构/动力学租约分层')),
    entry('baaiworm:conductance_layer', 'cell', 'other',
          writer='优化算法产出（data-driven tuning）',
          reader='BAAIWorm 399 神经元×双隔室×16 通道电导'
                 '（9,188 参数槽）',
          source='§9.13 BAAIWorm 审计',
          lease=lease_field('fit_set', 'modeling', 'explicit',
                            uniqueness='one-of-many',
                            note='正文自承全部电导经 optimization '
                                 'algorithm tuned——Par14 作者自承 '
                                 '"one of the possible solutions"（§9.13.3 '
                                 '声称锚）；16 家族初值文献溯源 6/16 可溯'
                                 '（Goodman 1998 PMID 9581767 经 Nicoletti '
                                 '2019 间接传承），终值均为优化产物')),
    entry('baaiworm:environment_layer', 'cell', 'other',
          writer='编译二进制内固化参数（C++ 求解器）',
          reader='BAAIWorm 环境交互子模型',
          source='§9.13.5 谬误二清点：环境层参数来源',
          lease=lease_field('unattributed', 'modeling', 'unattributed',
                            uniqueness='non-identified',
                            note='环境物理层参数封装在编译二进制内'
                                 '不可审计——不是作者隐瞒而是 C++ '
                                 '求解器常态（§9.13.5）；不可审计='
                                 '不可辨识')),
]


def orphan_query(entries):
    """§14.2 原型：谁没有阅读者（无读者条目试查）。"""
    return [e for e in entries if not e.get('reader')]


def strongest_readings(entries, k=5):
    """分子层读取信道强度排序（kd 升序 = 阅读最强）。"""
    mol = [e for e in entries if e.get('kd_estimate_nm') is not None]
    return sorted(mol, key=lambda e: e['kd_estimate_nm'])[:k]


def summary(entries):
    by_layer, by_type = {}, {}
    for e in entries:
        by_layer[e['layer']] = by_layer.get(e['layer'], 0) + 1
        by_type[e['constraint_type']] = by_type.get(
            e['constraint_type'], 0) + 1
    return {'total': len(entries), 'by_layer': by_layer,
            'by_type': by_type,
            'orphans': len(orphan_query(entries)),
            'kd_extrapolated': sum(1 for e in entries
                                   if e.get('kd_extrapolated'))}


def main():
    print('── 普查编译：分子层（surfaceome 91 边 → 条目）──')
    molecular = compile_molecular_layer()
    print(f'分子条目: {len(molecular)}（kd 外推标记 '
          f'{sum(1 for e in molecular if e["kd_extrapolated"]) }）')
    print('── 编目条目（52 条件目录 + §3.4）──')
    print(f'编目条目: {len(CURATED)}')
    entries = (molecular + CURATED + ECELL_LEASE_ENTRIES
               + BAAIWORM_LEASE_ENTRIES)
    ok, errors = validate(entries)
    print(f'── schema 校验: {"PASS" if ok else "FAIL"} ──')
    for er in errors[:5]:
        print('  !', er)
    orph = orphan_query(entries)
    print(f'── 孤儿查询（无读者条目）: {len(orph)} 条 ──')
    for e in orph[:5]:
        print('  -', e['entry_id'])
    print('── 最强读取 top5（kd 升序）──')
    for e in strongest_readings(entries, 5):
        print(f"  {e['entry_id']}: comp={e['complement_score']} "
              f"kd={e['kd_estimate_nm']}nM (外推={e['kd_extrapolated']})")
    print('── 账本摘要 ──')
    print(json.dumps(summary(entries), ensure_ascii=False))
    print('── 租约消费端（§9.8.4 判据 iii）──')
    ls = lease_summary(entries)
    print(json.dumps(ls, ensure_ascii=False))
    cross = filter_by_lease(entries, source_type='cross_species')
    print(f'filter_by_lease(source_type=cross_species): {len(cross)} 条'
          '（换谱系后失效的约束）')
    out = {'schema_version': SCHEMA_VERSION,
           'field_table': {'required': list(REQUIRED),
                           'geometry_fields': list(GEOMETRY_FIELDS),
                           'all_fields': list(entries[0].keys())},
           'compiled': datetime.now().strftime('%Y-%m-%d %H:%M'),
           'summary': summary(entries),
           'validation': {'ok': ok, 'errors': errors[:20]},
           'entries': entries}
    outp = os.path.join(HERE, 'constraint_census_v1.json')
    with open(outp, 'w', encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    book_out = os.path.normpath(os.path.join(
        HERE, '..', '..', '..', 'docs', 'book_aicb', 'code', 'output',
        'constraint_census_v1.json'))
    if os.path.isdir(os.path.dirname(book_out)):
        with open(book_out, 'w', encoding='utf-8') as f:
            json.dump(out, f, ensure_ascii=False, indent=1)
        print('副本:', book_out)
    print('输出:', outp)
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
