# Table1 Columns

当前 `table1` 是第一批 `72` 个蛋白，对应文件：

- `backend/data/table1_master_full.json`
- `backend/data/table1_master_table.csv`

`Step5_*` 列已从 `table1` 与默认 `master table` 产物中移除。当前列含义如下。

## 标识与基础信息

| 列名 | 含义 |
| --- | --- |
| `_id` | 表内主键；通常等于 PDB 编号。 |
| `Representative` | 代表蛋白/代表结构的 PDB 编号。 |
| `Protein_Name` | 蛋白名称。 |
| `Protein_Category` | 蛋白类别或功能归类。 |
| `Species` | 物种来源。 |
| `MonomerSeq` | 单体氨基酸序列。 |
| `Length;Oligomer;Monomer` | 复合字段，含总残基长度、聚体数、单体长度。 |
| `Cluster_ID` | 当前聚类编号。 |

## 相似性与对接结果

| 列名 | 含义 |
| --- | --- |
| `Similarity_Score` | 蛋白相似性评分。 |
| `Best_Rank` | DiffDock 最优 pose 的排名。 |
| `Best_Confidence` | DiffDock 最优 pose 的置信分数。 |
| `Best_SDF` | 最优配体 pose 对应的 SDF 文件路径。 |
| `Receptor_PDB` | 受体结构文件路径。 |
| `NearestDistanceTo3HisZn` | 最优 pose 到 3HIS-Zn 活性中心的最近距离。 |
| `AllZnDistancesSorted` | 配体到 Zn 位点的距离列表，通常已排序。 |

## Zn 位点与邻域统计

| 列名 | 含义 |
| --- | --- |
| `ZN_Depth` | Zn 位点埋藏深度。 |
| `ZN_SASA` | Zn 位点溶剂可及表面积。 |
| `ZN_Depth_Rounded` | `ZN_Depth` 的近似/分箱版本。 |
| `Neighbor_Count` | Zn 周围邻近残基数量。 |
| `Neighbor_List` | Zn 周围邻近原子/残基列表及其距离等信息。 |
| `Neighbor_Profile_Similarity` | 邻域残基组成的相似性评分。 |

## 3HIS / Zn 几何相关

| 列名 | 含义 |
| --- | --- |
| `TriHisSatisfied` | 是否满足 3 个组氨酸与 Zn 协同配位的条件。 |
| `TriHisCountMax` | 检测到的最大配位组氨酸数。 |

## 最佳 pocket 统计

| 列名 | 含义 |
| --- | --- |
| `BestPocket_ID` | 选中的最佳 pocket 编号。 |
| `BestPocket_Score` | 最佳 pocket 的综合评分。 |
| `BestPocket_Druggability` | 最佳 pocket 的药物可成药性评分。 |
| `BestPocket_Volume` | 最佳 pocket 体积。 |
| `BestPocket_TotalSASA` | 最佳 pocket 总溶剂可及表面积。 |
| `BestPocket_PolarSASA` | 最佳 pocket 极性表面积。 |
| `BestPocket_ApolarSASA` | 最佳 pocket 非极性表面积。 |
| `BestPocket_AlphaSpheres` | 构成该 pocket 的 alpha spheres 数量。 |
| `BestPocket_HisCount` | 最佳 pocket 附近组氨酸数量。 |
| `BestPocket_ZnCount` | 最佳 pocket 附近 Zn 数量。 |
| `BestPocket_MinDistToZn` | 最佳 pocket 到最近 Zn 的距离。 |
| `BestPocket_ZnMatch` | 最佳 pocket 是否命中/覆盖 Zn 活性位点。 |
| `BestPocket_SelectRule` | 该最佳 pocket 被选中的规则。 |

## 苯乙酮通过性复核列

| 列名 | 含义 |
| --- | --- |
| `苯乙酮_综合判断` | 基于活性 pocket 半径与现有苯乙酮 docking pose 的综合结论。 |
| `苯乙酮_尺寸可进入` | 只看活性 pocket 几何尺寸时，苯乙酮是否有望进入 Zn 邻域。 |
| `苯乙酮_活性口袋命中Zn` | 当前选中的 best pocket 是否实际覆盖 Zn 活性位点。 |
| `苯乙酮_通道瓶颈半径_A` | 从 Zn 附近到口袋外缘估算得到的最窄通道半径，单位 Å。 |
| `苯乙酮_口袋最小半径_A` | 活性 pocket 内 alpha-sphere 的最小半径，单位 Å。 |
| `苯乙酮_O到Zn最短距离_A` | 现有苯乙酮 docking pose 中，羰基 O 到 Zn 的最短距离，单位 Å。 |
| `苯乙酮_任意原子到Zn最短距离_A` | 现有苯乙酮 docking pose 中，任意原子到 Zn 的最短距离，单位 Å。 |
| `苯乙酮_判断说明` | 对综合判断的简短中文解释。 |

## 新增 Zn 配位统计列

| 列名 | 含义 |
| --- | --- |
| `Zn_CoordResidueCount` | 与 Zn 原子相邻且具有配位/价键作用的残基总数。 |
| `Zn_CoordHisCount` | 其中属于组氨酸的配位残基数量。 |
| `Zn_CoordNonHisCount` | `Zn_CoordResidueCount - Zn_CoordHisCount`，即非组氨酸配位残基数量。 |
| `Zn_CoordResidues` | 参与 Zn 配位的残基明细，包含残基标签和距离。 |
| `Zn_CoordSite` | 采用的 Zn 位点标识，通常为 `chain:seq`。 |
| `Zn_CoordNote` | 若结构中存在多个 Zn 位点，这里记录未被主列采用的其他 Zn 位点明细。 |

## 补充说明

- `Zn_Coord*` 这组列是基于结构文件重新计算得到的，不依赖 `Step5_*` 列。
- 当前配位残基识别规则使用蛋白常见供体原子，并按 `Zn` 周围 `3.2 Å` 以内筛选。
- 如果一个结构中有多个 `Zn`，主列优先写入 `Zn_CoordNonHisCount` 最少的那个位点；其余 Zn 位点写入 `Zn_CoordNote`。
- `苯乙酮_尺寸可进入` 的阈值基于当前仓库中苯乙酮三维构象的保守估计，通过半径约为 `2.9 Å`。
- `苯乙酮_综合判断` 不是纯 pocket 体积判断，而是同时参考了 `best pocket` 是否命中 Zn、通道瓶颈半径，以及现有苯乙酮 pose 到 Zn 的最近距离。
- 若后续加入 `table2`、`table3`，建议继续沿用同一套列定义，这样前端表格和导出逻辑不用再改。
