#!/usr/bin/env python3
"""constraint_census_schema_v15.py — schema v1.5：三态枚举扩展为四态.

v1.5 演进（2026-11-14，用户级理论决策）：LEASE_UNIQUENESS 新增第四态
'paradigmatic'——范式性不可辨识（NN 权重等深度学习范式产物：随机初始化
+非凸训练有意设计，解空间连续巨大但功能等效可复现）。

语义边界（与 non-identified 的区分——E10 预审边界发现）：
  non-identified   意外不可辨识（unintended）——参数来源不明/数据缺失，
                   影响下游约束读取可信度（来源债务，贬义）；
  paradigmatic     范式性不可辨识（paradigmatic）——建模范式有意设计
                   （NN 权重/生成式模型参数），模型功能正常可复现，
                   确定性输入-输出映射（中性）。

判读例：DeepBind weight-set → paradigmatic（训练范式使然，非来源债务）；
        ODE unattributed 参数 → non-identified（来源缺失，非有意设计）。

零破坏纪律：v1.4 三态历史数据不动（本文件为 v1.4 复制+第四态追加）；
v1.5 校验向后兼容（v1.4 三态仍合法）。
"""
import math

# ── v1.5 枚举（第四态追加）──
LEASE_UNIQUENESS = ('identified', 'non-identified', 'one-of-many',
                    'paradigmatic')

# 三态语义注释（v1.4 保留 + v1.5 第四态定义）
UNIQUENESS_SEMANTICS = {
    'identified': '参数唯一可辨识（有充分实验约束辨识出唯一解）',
    'one-of-many': '解空间中的一个解（BAAIWorm Par14 作者自承 '
                   '"one of the possible solutions"）',
    'non-identified': '意外不可辨识（unintended）——来源不明/数据缺失，'
                      '影响下游可信度（来源债务）',
    'paradigmatic': '范式性不可辨识（paradigmatic）——建模范式有意设计'
                    '（NN 权重/生成式参数），功能等效可复现（中性）',
}


def make_lease(source_type=None, lease_axis=None, evidence=None,
               uniqueness=None, note=None):
    """租约字段构造器（v1.5 四态口径，接口与 v1.4 兼容）。"""
    return {'source_type': source_type, 'lease_axis': lease_axis,
            'evidence': evidence, 'uniqueness': uniqueness, 'note': note}


def validate_lease(lease):
    """v1.5 校验（向后兼容：v1.4 三态仍合法）。"""
    if lease is None:
        return True, []
    errors = []
    uq = lease.get('uniqueness')
    if uq is not None and uq not in LEASE_UNIQUENESS:
        errors.append(f'lease.uniqueness 非法: {uq}')
    return (len(errors) == 0), errors


if __name__ == '__main__':
    # 自检：四态校验+向后兼容
    ok1, e1 = validate_lease({'uniqueness': 'identified'})
    ok2, e2 = validate_lease({'uniqueness': 'paradigmatic'})
    ok3, e3 = validate_lease({'uniqueness': '非法值'})
    assert ok1 and ok2 and not ok3, 'v1.5 自检失败'
    print('v1.5 自检：三态向后兼容+paradigmatic 合法+非法值拒绝 PASS')
