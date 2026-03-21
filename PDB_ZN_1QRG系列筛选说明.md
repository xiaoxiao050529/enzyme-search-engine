# PDB-ZN 中 1QRG 系列筛选说明（可复现）

本文档总结我们在 PDB-ZN 工作流中，如何从全量库逐步筛到约 70 个候选蛋白，并说明 1QRG 在流程中的处理方式。

## 1. 目标与口径

- 目标：在 Zn/HIS 相关蛋白中筛出与 1QRG 关注位点特征一致的一批候选。
- 口径：使用 PDB-ZN 工作流 Step1~Step5，按固定参数复现。
- 本次复现实测结果：最终 72 个候选（“70 个左右”）。

## 2. 起始数据规模

当前数据库状态（复现时）：

- total_proteins: 7688
- with_master: 72

说明：`total_proteins` 是导入到 `pdbzn.sqlite` 的总条目数。

## 3. 1QRG 在流程中的定位

- 1QRG 在主表中可查到（历史记录存在）。
- 但在 Step2 聚类前会做“本地 PDB 可用性检查”，1QRG 属于无本地 TM-align PDB 的条目，被统一剔除。
- 本次 Step2 统计中 `dropped_no_pdb_count=36`，且包含 `1QRG`。

这意味着：1QRG 是“规则锚点/参考对象”，但不直接进入当前聚类输入集合。

## 4. 本次复现参数（关键）

以下参数用于得到“约 70 个”的结果：

- Step1
  - metal_ion: `ZN`
  - residue_requirements: `HIS:3`
- Step2
  - mode: `preset_clusters`
  - keep_only_cluster_members: `true`
  - cluster_threshold: `0.7`
- Step3
  - ligand_distance: `5.0`
  - metal_distance: `8.0`
  - max_residue_per_oligomer: `180`
- Step4
  - require_his3: `true`
- Step5
  - run_fpocket: `false`（本次复现为固定统计，关闭补跑）
  - write_to_master: `false`（仅生成输出文件，不覆盖主表）

## 5. 每一步筛选依据与数量

### Step1：初筛（金属与残基规则）

筛选依据：

- 金属离子为 ZN
- 满足 HIS:3 规则

结果数量：

- Step1 输出：1265

### Step2：聚类（含无 PDB 剔除）

筛选依据：

- 使用 `preset_clusters` 聚类映射
- 在进入聚类前，统一剔除无本地 PDB 的条目

结果数量：

- Step1 输入：1265
- 剔除无 PDB：36（包含 1QRG）
- Step2 实际聚类输入：1229
- 输出簇数：148

### Step3：第二次结构筛选

筛选依据：

- 排除 ZN 5A 内大配体干扰（ligand_distance=5.0）
- 排除 ZN 8A 内其他金属干扰（metal_distance=8.0）
- 限制残基数/几聚体（max_residue_per_oligomer=180）

结果数量：

- Step3 输入簇：148
- Step3 剔除：70
- Step3 保留：78

### Step4：3HIS 几何复核（无二次聚类）

筛选依据：

- 检查 3 个 HIS 与 ZN 的几何关系
- 输出 ZN 周围 5A 全部残基、HIS-ZN 键长、键角
- 剔除对称/空间构象导致的误配

结果数量：

- Step4 复核输入：78
- Step4 剔除：6
- Step4 通过：72

### Step5：主表收口

处理依据：

- 合并 Step4 通过结果与主表字段
- 同步写入 Step4 几何字段（TriHis/键长键角/5A 残基/空间构象通过）

结果数量：

- Step5 appended_count：72
- 最终候选：72（约 70）

输出文件（本次复现）：

- `/root/PDB_ZN/zn_his_master_table_step5_1qrg72.csv`

## 6. 复现操作（推荐）

### 6.1 前端页面操作

在 `pdbzn_workflow.html`：

1. Step1 保持 `ZN + HIS:3`
2. Step2 选择 `preset_clusters`，阈值 `0.7`
3. Step3 将 `max_residue_per_oligomer` 调为 `180`
4. Step4 勾选“必须通过3HIS复核”
5. Step5 关闭 `run_fpocket`，关闭 `write_to_master`，执行收口

### 6.2 后端函数复现（同口径）

可直接调用：

- `_pdbzn_step3_filter_workflow(...)`
- `_pdbzn_step4_validate_workflow(...)`
- `_pdbzn_step5_finalize_workflow(...)`

并使用本文件第 4 节参数。

## 7. 为什么最终是“约 70”

在当前数据版本下，主要由 Step3 的 `max_residue_per_oligomer` 控制收敛速度：

- 200 时，Step4 约 74
- 180 时，Step4 约 72
- 160 时，Step4 约 69

因此把 `max_residue_per_oligomer` 设为 180，可以稳定落在“70 左右”的区间。
