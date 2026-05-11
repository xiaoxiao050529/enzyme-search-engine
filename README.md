# 项目完整架构与运行说明

本文件用于统一替代和整合以下几份分散说明：

- `README.md`
- `架构说明_论文素材.md`
- `PDB_ZN_1QRG系列筛选说明.md`
- `backend-port-conflict.md`

这份文档重点解决 3 个问题：

1. 代码到底是怎么组织的
2. 后端到底读哪些数据、写哪些运行时文件
3. 系统从输入到输出到底怎么跑

## 0. 推荐启动命令

如果你要正常把这个项目整体打开起来，推荐直接在项目根目录执行：

```bash
bash start.sh
```

或者：

```bash
./start.sh
```

如果你要直接以公网模式启动：

```bash
PUBLIC_MODE=1 ./start.sh
```

如需指定公网 IP 或域名用于终端提示输出：

```bash
PUBLIC_MODE=1 PUBLIC_HOST=your.domain.or.ip ./start.sh
```

这个脚本会自动：

- 启动或复用后端 API
- 启动或复用前端静态服务
- 自动处理常见端口冲突
- 输出可直接访问的页面地址

默认会打开这些地址之一：

- 工作流页：`http://127.0.0.1:8020/frontend/pdbzn_workflow.html`
- 首页：`http://127.0.0.1:8020/frontend/index.html`
- 主表页：`http://127.0.0.1:8020/frontend/master_table.html`

## 1. 先说结论

这个仓库本质上是一个“静态前端 + Python 单文件后端 + 本地数据目录 + 外部结构生信工具”的平台，面向 Zn 相关蛋白筛选、结构浏览、DiffDock 对接、fpocket 口袋分析和 TM-align 结构比对。

当前代码现实和旧文档有 4 个关键差异，需要先校正：

- 当前后端不是 Flask，而是 Python 标准库里的 `ThreadingHTTPServer`。
- 当前真正工作的数据库是 `backend/runtime/pdbzn.sqlite`，不是 `backend/pdbzn_workflow.db`。
- `backend/pdbzn_workflow.db` 当前是一个 0 字节空文件，属于历史遗留占位，不参与当前运行。
- 仓库里的静态数据和 runtime 数据并不总是完全同步。
  例如当前 `backend/data/master_table.csv` 有 121 行，而 `backend/data/data.json` / `backend/data/master_full.json` 是 72 条，`backend/runtime/pdbzn.sqlite` 中 `with_master` 也是 72。
- 因此前端静态浏览结果、workflow 动态结果和 runtime 导出结果，必须先确认“你现在看的是哪一层数据”，不能默认它们天然一致。

因此理解这个项目时，最好把数据分成两层：

- `backend/data/`：前端直接浏览的静态快照
- `backend/runtime/`：后端 API 和任务执行产生的动态状态

另外，我已经把这两层数据整理出一个便于人工查看的目录视图：

- `backend/data_catalog/`

注意：

- 这是“整理视图”，用于按类别浏览数据。
- 程序实际运行仍然读取 `backend/data/` 和 `backend/runtime/`，所以不会因为整理目录而破坏现有页面和 API。

## 1.1 当前已有数据分类

如果只看当前仓库里的“已有数据”，可以把它们按用途整理成 7 类。最容易混淆的是 `master table`，它其实不是一个单文件，而是一组彼此配套的视图文件。

### A. 数据表注册表

文件：

- `backend/data/table_registry.json`

作用：

- 这是所有“可选数据表”的目录。
- 前端 `master_table.html`、`index.html`、`pdbzn_workflow.html` 会先读这里，决定当前有哪些表可以切换。
- 当前注册表里可见 9 个表，包括 `table1`、`table2`、`table3`、`6j4c`、`1jdi` 以及几个 workflow / 手工保存生成的新表。

### B. Master Table 数据组

典型文件：

- `backend/data/table1_master_full.json`
- `backend/data/table1_data.json`
- `backend/data/table1_master_table.csv`
- `backend/data/table1_ids.json`

作用：

- 这四个文件合起来，才是一张表的完整静态快照。
- 可以把它们理解成同一批数据的 4 种视图：
- `_master_full.json`：最完整的宽表，包含 `header + rows`，给 `master_table.html` 和详情联动使用。
- `_data.json`：较轻量的 `items` 列表，给首页、DiffDock 对比页等快速浏览页面使用。
- `_master_table.csv`：便于直接打开或导出的表格版本。
- `_ids.json`：只保留该表的 PDB / Representative 标识列表，方便筛选或复用。

补充说明：

- 根目录级的 `backend/data/master_full.json`、`backend/data/data.json`、`backend/data/master_table.csv` 是默认快照。
- `table1_*` 是当前主表快照；`table2_*`、`table3_*`、`6j4c_*`、`1jdi_*`、`my_new_table_*` 这些是 workflow 或手工保存后产生的派生表。

### C. 蛋白结构文件

文件：

- `backend/data/structures/*.pdb`
- `backend/data/structures/*.cif`

当前规模：

- `72` 个 `.pdb`
- `4` 个 `.cif`

作用：

- 这是蛋白三维结构的原始坐标文件。
- 首页结构浏览、TM-align、fpocket、workflow 导入都会依赖这里的结构数据。

补充说明：

- `backend/data/structures/*_out/` 是部分结构跑 pocket / 外部工具后留下的输出目录，不是主数据表本身，而是分析产物缓存。

### D. 配体文件

文件：

- `backend/data/ligands/*.sdf`

当前规模：

- `71` 个 `.sdf`

作用：

- 存放与蛋白对应的配体构象。
- 首页加载最佳配体、DiffDock 对比页显示 pose 时都要读这些文件。

### E. DiffDock 结果数据

文件：

- `backend/data/diffdock_index.json`
- `backend/data/diffdock/<PDB>/rank*.sdf`

当前规模：

- `diffdock_index.json` 当前索引 `71` 个蛋白
- `backend/data/diffdock/` 下当前有 `71` 个 PDB 子目录

作用：

- `diffdock_index.json` 是 DiffDock 结果的目录索引，告诉前端每个蛋白有哪些 pose。
- 每个 `diffdock/<PDB>/` 目录里存具体的 `rank*.sdf` 对接结果。
- `diffdock_compare.html` 主要就是读这一组数据。

### F. fpocket / pocket 结果数据

文件：

- `backend/data/pockets/<PDB>/best_pocket_atm.pdb`
- `backend/data/pockets/<PDB>/best_pocket_vert.pqr`
- `backend/data/pocket.csv`

当前规模：

- `backend/data/pockets/` 下当前有 `72` 个 PDB 子目录

作用：

- `best_pocket_atm.pdb` 和 `best_pocket_vert.pqr` 是每个蛋白 pocket 结果的可视化输入。
- `fpocket.html` 和首页的 pocket 叠加显示会直接用这些文件。
- `pocket.csv` 可以视为 pocket 结果的汇总表。

### G. 运行时数据库与动态产物

文件 / 目录：

- `backend/runtime/`
- `backend/runtime/pdbzn.sqlite`
- `backend/pdbzn_workflow.db`

作用：

- `backend/runtime/pdbzn.sqlite` 是当前后端 workflow 真正在维护的动态数据库。
- `backend/runtime/` 里还会放任务目录、日志、临时结果和导出表。
- `backend/pdbzn_workflow.db` 是历史遗留文件，不应再当作当前主数据库理解。

### H. 数量不一致时该怎么理解

当前这些类别的数量并不完全相等，例如：

- 结构文件是 `72` 个 `.pdb` + `4` 个 `.cif`
- 配体 `.sdf` 是 `71` 个
- DiffDock 结果目录是 `71` 个
- pocket 结果目录是 `72` 个

这通常不表示仓库出错，而表示：

- 不是每个蛋白都已经完成了全部下游分析
- 某些蛋白只有结构，没有配体或 DiffDock 结果
- 某些页面读取的是“主表快照”，某些页面读取的是“下游分析产物目录”

### 一句话记忆

如果只想快速记：

- “表数据”看 `table_registry.json + *_master_full.json / *_data.json / *_csv / *_ids.json`
- “蛋白结构”看 `structures/`
- “配体”看 `ligands/`
- “DiffDock pose” 看 `diffdock_index.json + diffdock/`
- “fpocket 结果”看 `pockets/ + pocket.csv`
- “运行时状态”看 `backend/runtime/`

## 2. 仓库结构总览

建议把仓库理解成下面 6 层：

```text
.
├── frontend/
│   ├── index.html
│   ├── master_table.html
│   ├── diffdock_compare.html
│   ├── fpocket.html
│   ├── tmalign.html
│   └── pdbzn_workflow.html
├── backend/
│   ├── diffdock_api_server.py
│   ├── generate_data.py
│   ├── ingest_ranking.py
│   ├── data/
│   ├── runtime/
│   ├── pdbzn_workflow.db
│   └── __pycache__/
├── scripts/
│   └── download_pdb.sh
├── data/
│   └── pdb/                   # 额外下载的 PDB，非后端默认主路径
├── 各类说明文档 *.md
└── ids*.txt / ids72.tet 等辅助清单
```

各层职责如下：

- `frontend/`：页面和浏览器端逻辑，负责 UI、参数输入、结果展示、3D 交互、轮询任务状态。
- `backend/diffdock_api_server.py`：统一后端入口，负责 API、工作流、结构解析、任务调度、文件输出。
- `backend/data/`：仓库自带静态数据，前端即使不跑重计算，也可以直接浏览这些数据。
- `backend/runtime/`：运行时产物，包含 SQLite、任务目录、日志、导出表。
- `scripts/`：辅助脚本，目前主要是批量下载 PDB。
- 根目录说明文档：历史分散说明，现已被本文档整合。

## 3. 代码架构

## 3.1 前端架构

前端是 6 个独立 HTML 页面，没有前端框架，核心是原生 JavaScript + 3Dmol.js + `fetch()`。

页面职责如下：

- `frontend/index.html`
  - 主入口页
  - 加载 `backend/data/data.json` 或 `backend/data/master_full.json`
  - 展示结构、配体、Zn 邻域、口袋高亮、序列联动
- `frontend/master_table.html`
  - 主表浏览页
  - 读取 `backend/data/master_full.json`
  - 支持字段过滤、排序、跳转
- `frontend/diffdock_compare.html`
  - DiffDock pose 对比页
  - 读取 `backend/data/diffdock_index.json` 和 `backend/data/data.json`
  - 支持本地 pose 对比，也支持在线提交 DiffDock
- `frontend/fpocket.html`
  - 在线提交 fpocket
  - 轮询任务状态
  - 加载 pocket 参数、pocket 点云和结构文件
- `frontend/tmalign.html`
  - 在线提交 TM-align
  - 支持 1 对 1 和 1 对多
- `frontend/pdbzn_workflow.html`
  - PDB-ZN 工作流控制台
  - 负责 Step1 到 Step5 参数输入、调用 API、表格渲染

前端有两个设计特点：

- 所有页面都尽量支持自动探测后端地址。
- 页面既能直接读静态文件，也能通过后端 `/api/data/file` 来取数据。

当前常见端口约定：

- 前端静态页：`8020`
- 后端默认 API：`8015`
- 某些页面会继续尝试 `8017`、`8016`、`8021`、`8019`、`8018`

## 3.2 后端架构

当前后端几乎全部集中在 `backend/diffdock_api_server.py` 这一份文件里。它不是“薄路由 + 多模块服务”的结构，而是“单文件集成式服务”。

可以把这个文件拆成 7 个逻辑块来理解。

### A. 基础路径与环境探测

负责：

- 项目目录定位
- `backend/data`、`backend/runtime` 定位
- 外部 `PDB_ZN` 数据目录探测
- DiffDock / fpocket / TM-align 二进制探测

关键变量：

- `BASE_DIR`
- `PROJECT_DIR`
- `DATA_DIR`
- `RUNTIME_DIR`
- `PDBZN_BASE_DIR`
- `PDBZN_PDB_DIR`
- `PDBZN_DB_PATH`
- `PDBZN_EXPORT_DIR`
- `HOST`
- `PORT`

外部依赖定位函数：

- `discover_infer_py()`
- `discover_fpocket_bin()`
- `discover_tmalign_bin()`

### B. 结构解析与工具输出解析

负责读取和理解结构文件与外部工具输出。

包含：

- mmCIF / PDB 原子解析
  - `_pdbzn_mmcif_atom_rows()`
  - `_pdbzn_pdb_atom_rows()`
  - `_pdbzn_structure_atom_rows()`
- 几何计算
  - `_pdbzn_dist()`
  - `_pdbzn_angle_deg()`
  - `_pdbzn_step4_geometry_report()`
- fpocket 输出解析
  - `parse_fpocket_info_file()`
  - `parse_fpocket_vert_file()`
  - `build_fpocket_pocket_payload()`
- TM-align 日志解析
  - `parse_tmalign_log()`

### C. PDB-ZN 工作流引擎

负责 Step1 到 Step5 的主流程。

核心函数如下：

- `_pdbzn_run_workflow()`
- `_pdbzn_cluster_workflow()`
- `_pdbzn_step3_filter_workflow()`
- `_pdbzn_step4_validate_workflow()`
- `_pdbzn_step5_finalize_workflow()`

辅助能力：

- `_pdbzn_workflow_defaults()`
- `_pdbzn_cluster_rows()`
- `_pdbzn_cluster_rows_from_step1_clusters()`
- `_pdbzn_cluster_rows_tmalign()`
- `_pdbzn_step3_structure_reasons()`
- `_pdbzn_step5_final_score()`

### D. runtime 数据库管理

负责 `backend/runtime/pdbzn.sqlite` 的读写。

关键函数：

- `_pdbzn_db_connect()`
- `_pdbzn_db_ensure()`
- `_pdbzn_db_stats()`
- `_pdbzn_import_database()`

当前 `proteins` 表结构如下：

```text
pdb_id
file_name
file_path
has_zn_his_cluster
metal_ions
his_count_max
neighbor_residue_counts
similarity_score
protein_name
protein_category
details
imported_at
```

注意：

- 这个库是工作流的查询底座。
- 它保存的是“结构文件 + 主表关键字段”的导入结果。
- 当前库中有 `7688` 条记录，其中 `72` 条带 `protein_name`，也就是当前“with_master=72”的来源。

### E. 任务运行器

负责启动外部计算任务。

包括：

- `run_job()`：DiffDock
- `run_fpocket_job()`：fpocket
- `run_tmalign_job()`：TM-align
- `recover_fpocket_job()`：后端重启后恢复 fpocket 任务状态

运行方式都是：

- 创建 job 目录
- 写入输入文件
- 启动子进程
- 将 stdout/stderr 写入 `run.log`
- 轮询状态时返回 job 元信息

### F. HTTP API 路由层

通过 `class Handler(BaseHTTPRequestHandler)` 提供所有接口。

GET 类接口负责：

- 健康检查
- 读取任务状态
- 读取日志
- 读取静态数据文件
- 读取任务输出文件
- 获取工作流默认配置

POST 类接口负责：

- 提交 DiffDock / fpocket / TM-align
- 运行 PDB-ZN Step1 到 Step5
- 手动触发数据库导入

### G. 服务启动层

最后通过：

```python
server = ThreadingHTTPServer((HOST, PORT), Handler)
```

启动服务。

所以这是一个：

- 单进程
- 多线程请求处理
- 文件系统重依赖
- 本地工具强依赖

的后端。

## 3.3 其他后端脚本

除了主服务文件，后端还有两个数据构建脚本。

### `backend/generate_data.py`

作用：

- 从主表 CSV 生成前端更容易直接使用的 JSON

输出：

- `backend/data/data.json`
- `backend/data/master_full.json`

特点：

- `data.json` 是更轻的列表数据
- `master_full.json` 是完整表格视图数据

### `backend/ingest_ranking.py`

作用：

- 从 ranking CSV 提取 PDB ID
- 尝试复制或下载结构文件
- 生成 `backend/data/ranking.json`

这个脚本更像历史数据准备脚本，不是当前主工作流的核心入口。

## 3.4 辅助脚本

### `scripts/download_pdb.sh`

作用：

- 根据 ID 清单批量从 RCSB 下载 `.pdb`

现在支持：

- 自定义输入 ID 文件
- 自定义输出目录
- 自定义失败清单文件

例如：

```bash
bash scripts/download_pdb.sh ids72.txt data/pdb data/failed_ids72.txt
```

这个脚本下载到的是根目录下的 `data/pdb/`，它不是后端默认读取结构的主路径，但可作为额外结构仓库使用。

## 4. 后端数据集说明

这一节只讲“后端会读到哪些数据”。

## 4.1 `backend/data/`：静态数据层

`backend/data/` 是前端静态浏览和后端兜底查找的第一层数据源。

当前主要文件如下。

### 1. `backend/data/master_table.csv`

作用：

- 当前仓库内的主表 CSV 快照
- 是 Step5 写回的默认主表候选文件之一
- 也是 `generate_data.py` 理论上的源表

当前状态：

- 121 行
- 57 列

典型字段：

- `Representative`
- `Similarity_Score`（相对参考结构 1QRG 的分数，不是 Step2 聚类相似性）
- `Best_Confidence`
- `NearestDistanceTo3HisZn`
- `ZN_Depth`
- `Length;Oligomer;Monomer`
- `Protein_Name`
- `Protein_Category`
- `Cluster_ID`
- `Neighbor_List`
- `TriHisSatisfied`
- `TriHisCountMax`
- `Receptor_PDB`
- `Best_SDF`

### 2. `backend/data/pocket.csv`

作用：

- 预先整理好的 pocket 汇总表

当前状态：

- 1806 行
- 22 列

典型字段：

- `Representative`
- `Protein_Name`
- `Protein_Category`
- `Receptor_PDB`
- `pocket_id`
- `score`
- `druggability_score`
- `volume`
- `total_sasa`
- `pocket_his_count`
- `pocket_zn_count`

### 3. `backend/data/data.json`

作用：

- 主页面常用的轻量主表 JSON

当前状态：

- `items` 长度为 72

单条记录字段：

- `id`
- `name`
- `cluster`
- `receptor_pdb`
- `receptor_rel`
- `best_sdf`
- `species`
- `monomer_seq`
- `residue_length`
- `oligomer`
- `monomer_length`

### 4. `backend/data/master_full.json`

作用：

- Master Table 页面使用的完整行数据

当前状态：

- `rows` 长度为 72
- `header` 长度为 36

它本质上是面向前端表格展示的 JSON 化主表。

### 5. `backend/data/diffdock_index.json`

作用：

- DiffDock 对比页的索引数据
- 用于列出每个蛋白有哪些 pose 文件

当前状态：

- `items` 长度为 71
- `total_proteins = 71`

单条记录字段：

- `id`
- `pose_count`
- `poses`

每个 `pose` 包含：

- `file`
- `filename`
- `rank`
- `confidence`

### 6. `backend/data/structures/`

作用：

- 仓库内置结构文件目录
- 是后端默认的结构查找目录之一

当前状态：

- 顶层 76 个结构文件
- 其中大多数是 72 个主表蛋白的 `.pdb`
- 还包含少量额外 `.cif`
- 另有 20 个 `*_out` 子目录，存放某些结构相关衍生产物

### 7. `backend/data/ligands/`

作用：

- 仓库内置配体 SDF

当前状态：

- 71 个 `.sdf`

注意：

- 当前 72 个主表条目里，有 1 个不在这里面。
- 缺失的是 `2IMR`。

### 8. `backend/data/diffdock/`

作用：

- 每个蛋白的本地 DiffDock pose 文件目录

当前状态：

- 71 个蛋白目录
- 总大小约 1.12 MB

目录结构示例：

```text
backend/data/diffdock/1AF0/
├── rank1_confidence-0.82.sdf
├── rank2_confidence-1.34.sdf
├── ...
└── rank10_confidence-7.29.sdf
```

### 9. `backend/data/pockets/`

作用：

- 每个蛋白的 pocket 相关静态文件目录

当前状态：

- 72 个蛋白目录
- 总大小约 0.74 MB

通常每个目录包含 2 个核心文件：

- pocket 信息文件
- pocket 点云/顶点文件

## 4.2 `backend/runtime/`：运行时状态层

`backend/runtime/` 是后端执行后的动态状态目录。

当前包含 4 类内容。

### 1. `backend/runtime/pdbzn.sqlite`

作用：

- 当前工作流真实使用的 SQLite 数据库

当前状态：

- 大小约 648 KB
- `proteins` 表 7688 行
- `with_master = 72`

它是 Step1 的直接输入底座，也是 `/api/pdbzn/workflow/config` 返回的数据库状态来源。

### 2. `backend/runtime/diffdock_jobs/`

作用：

- 在线 DiffDock 提交任务目录

当前状态：

- 3 个 job 目录

典型结构：

```text
backend/runtime/diffdock_jobs/<job_id>/
├── input/
│   ├── receptor.pdb
│   └── ligand.sdf
└── run.log
```

如果推理成功，代码设计上还会有 `output/` 目录和 pose SDF。

当前样本日志显示，历史运行里曾出现环境依赖缺失：

- `ModuleNotFoundError: No module named 'yaml'`

这说明 DiffDock 在线推理依赖独立 Python 环境，不能只靠仓库静态文件。

### 3. `backend/runtime/fpocket_jobs/`

作用：

- 在线 fpocket 任务目录

当前状态：

- 27 个 job 目录
- 总大小约 24.63 MB

典型结构：

```text
backend/runtime/fpocket_jobs/<job_id>/
├── input/
│   ├── <pdb>.pdb
│   └── <pdb>_out/
│       ├── <pdb>_info.txt
│       ├── <pdb>_out.pdb
│       ├── <pdb>_pockets.pqr
│       └── pockets/
│           ├── pocket1_atm.pdb
│           ├── pocket1_vert.pqr
│           └── ...
└── run.log
```

这是 runtime 中最“重”的目录，因为 fpocket 会落很多 pocket 文件。

### 4. `backend/runtime/tmalign_jobs/`

作用：

- 在线 TM-align 任务目录

当前状态：

- 38 个 job 目录
- 总大小约 18.28 MB

典型结构：

```text
backend/runtime/tmalign_jobs/<job_id>/
├── input/
│   ├── pdb1.pdb
│   └── pdb2.pdb
└── run.log
```

日志里保留原始 TM-align 输出，后端再从日志中提取：

- `tm_score_1`
- `tm_score_2`
- `tm_score_max`
- `aligned_length`
- `rmsd`
- `seq_id`

### 5. `backend/runtime/exports/`

作用：

- 工作流导出结果

当前状态：

- 1 个导出文件
- `zn_his_master_table_step5.csv`

这个目录是 Step5 默认写导出文件的位置。

## 4.3 旧库与遗留文件

### `backend/pdbzn_workflow.db`

当前状态：

- 0 字节
- 无表结构

结论：

- 当前代码不会把它当成主数据库使用
- 它只是历史遗留文件名

### 运行状态和静态快照不一致的原因

仓库内同时存在：

- 静态主表快照
- runtime 导入库
- 历史导出 CSV
- 外部 `PDB_ZN` 路径可选接入

所以看到 `121`、`72`、`7688` 这几组数字同时存在是正常的，它们分别代表：

- 121：当前仓库里的主表 CSV 行数
- 72：当前前端静态主表 JSON 和 runtime 中有效主表条目数
- 7688：runtime SQLite 中的全部导入蛋白条目数

## 5. API 架构

当前 API 可以分成 5 组。

### 1. 健康检查

- `GET /api/diffdock/ping`
- `GET /api/fpocket/ping`
- `GET /api/tmalign/ping`
- `GET /api/pdbzn/workflow/config`

### 2. 静态数据访问

- `GET /api/data/file?path=...`

这个接口用于前端通过 API 读取 `backend/data/` 中的 JSON、结构、pose、pocket 文件。

### 3. 任务型接口

DiffDock：

- `POST /api/diffdock/submit`
- `GET /api/diffdock/status/<job_id>`
- `GET /api/diffdock/log/<job_id>`
- `GET /api/diffdock/file/<job_id>/<file>`

fpocket：

- `POST /api/fpocket/submit`
- `GET /api/fpocket/status/<job_id>`
- `GET /api/fpocket/log/<job_id>`
- `GET /api/fpocket/pockets/<job_id>`
- `GET /api/fpocket/file/<job_id>/<rel_path>`

TM-align：

- `POST /api/tmalign/submit`
- `GET /api/tmalign/status/<job_id>`
- `GET /api/tmalign/log/<job_id>`

### 4. PDB-ZN 工作流接口

- `POST /api/pdbzn/workflow/import`
- `POST /api/pdbzn/workflow/run`
- `POST /api/pdbzn/workflow/cluster`
- `POST /api/pdbzn/workflow/filter`
- `POST /api/pdbzn/workflow/validate`
- `POST /api/pdbzn/workflow/finalize`

### 5. 任务列表接口

- `GET /api/diffdock/jobs`
- `GET /api/fpocket/jobs`
- `GET /api/tmalign/jobs`

## 6. 运行方法

## 6.1 最小可运行方式

如果你只是想找“这个项目到底该敲哪条命令”，优先用：

```bash
bash start.sh
```

它会同时处理前端和后端，比手动分别启动更稳妥。

如果你只想打开页面和浏览仓库自带数据，不跑重计算，最小启动方式如下。

后端：

```bash
python3 backend/diffdock_api_server.py
```

前端：

```bash
python3 -m http.server 8020
```

然后访问：

- `http://127.0.0.1:8020/frontend/index.html`
- `http://127.0.0.1:8020/frontend/master_table.html`
- `http://127.0.0.1:8020/frontend/diffdock_compare.html`
- `http://127.0.0.1:8020/frontend/fpocket.html`
- `http://127.0.0.1:8020/frontend/tmalign.html`
- `http://127.0.0.1:8020/frontend/pdbzn_workflow.html`

## 6.2 后端默认地址

默认绑定：

- Host: `127.0.0.1`
- Port: `8015`

对应环境变量：

- `DIFFDOCK_API_HOST`
- `DIFFDOCK_API_PORT`

改端口示例：

```bash
DIFFDOCK_API_PORT=8017 python3 backend/diffdock_api_server.py
```

## 6.3 常用环境变量

### 后端服务

- `DIFFDOCK_API_HOST`
- `DIFFDOCK_API_PORT`

### 外部数据目录

- `PDBZN_BASE_DIR`
- `PDBZN_PDB_DIR`

### DiffDock

- `DIFFDOCK_INFER_PY`
- `DIFFDOCK_PYTHON_BIN`

### fpocket / TM-align

- `FPOCKET_BIN`
- `TMALIGN_BIN`

### 数据构建脚本

- `PDBZN_MASTER_TABLE`
- `PDBZN_RANKING_CSV`

## 6.4 启动后验证

先验证后端：

```bash
curl http://127.0.0.1:8015/api/pdbzn/workflow/config
```

验证 DiffDock API：

```bash
curl http://127.0.0.1:8015/api/diffdock/ping
```

如果要换端口，比如 `8017`，把上面的端口一起换掉。

## 6.5 端口冲突处理

如果出现：

```text
OSError: [Errno 98] Address already in use
```

优先按下面顺序处理：

1. 先 `curl` 检查 `8015` 是否已经有可用后端
2. 如果已有服务正常，就不要重复启动
3. 如果确实要重启，再查占用 PID
4. 如果不想影响旧服务，就换 `8017`

常用命令：

```bash
ss -ltnp '( sport = :8015 )'
ps -fp <PID>
kill <PID>
DIFFDOCK_API_PORT=8017 python3 backend/diffdock_api_server.py
```

## 6.6 外部工具依赖

### 只浏览静态数据

不强依赖：

- DiffDock
- fpocket
- TM-align

### 跑在线任务或完整工作流

需要至少满足：

- DiffDock 推理脚本可用
- fpocket 可执行文件可用
- TMalign 可执行文件可用

否则会出现：

- DiffDock 没有 `inference.py`
- fpocket / TMalign not found
- job 只有 `run.log`，没有结果文件

## 7. 运行逻辑：从输入到输出

这里按 3 条主链路讲。

## 7.1 链路一：静态浏览链路

适用页面：

- `index.html`
- `master_table.html`
- `diffdock_compare.html`

流程如下：

1. 浏览器打开 HTML 页面
2. 页面读取 `backend/data/*.json`
3. 页面根据 JSON 中的相对路径加载 `.pdb`、`.sdf`、pocket 文件
4. 3Dmol.js 在前端本地渲染
5. 用户得到结构、配体、口袋和表格展示

这条链路的特点：

- 不一定依赖 runtime
- 主要依赖 `backend/data/`
- 页面能用静态文件就优先用静态文件

## 7.2 链路二：在线任务链路

适用页面：

- `diffdock_compare.html`
- `fpocket.html`
- `tmalign.html`

以 fpocket 为例，完整路径如下：

1. 用户上传或粘贴 PDB 文本
2. 前端 `POST /api/fpocket/submit`
3. 后端创建 `backend/runtime/fpocket_jobs/<job_id>/input/<pdb>.pdb`
4. 后端后台线程启动 `fpocket`
5. `fpocket` 产出 `<pdb>_out/` 全套结果
6. 日志写入 `run.log`
7. 前端轮询 `/api/fpocket/status/<job_id>`
8. 前端再请求 `/api/fpocket/pockets/<job_id>`
9. 后端解析 `*_info.txt` 和 `pocket*_vert.pqr`
10. 前端渲染口袋参数和 pocket 点云

DiffDock 和 TM-align 也是同样的模式：

- 提交文本
- 写 job 目录
- 启动外部程序
- 写日志
- 查询状态
- 拉回结果

## 7.3 链路三：PDB-ZN 工作流链路

这是整个项目最核心的一条链。

### Step0：数据库准备

入口：

- `POST /api/pdbzn/workflow/import`
- 或 Step1 自动触发导入

逻辑：

1. 后端扫描结构文件来源
2. 读取主表 CSV
3. 按 `pdb_id` 建立 `proteins` 表
4. 从 PDB 坐标重新计算金属、ZN 邻域残基和最大配位 HIS 数
5. 合并主表中的名称、分类、相似度等辅助信息到 `pdbzn.sqlite`

输出：

- `backend/runtime/pdbzn.sqlite`

### Step1：金属和残基初筛

入口：

- `POST /api/pdbzn/workflow/run`

输入：

- 金属离子
- 残基需求，例如 `HIS:3`
- 名称过滤
- PDB ID 过滤

逻辑：

1. 从 `pdbzn.sqlite` 读 `proteins`
2. 按 `metal_ions` 过滤
3. 按 `neighbor_residue_counts` / `his_count_max` 过滤
4. 输出匹配行

输出：

- Step1 表格
- `summary.step1_count`

### Step2：聚类

入口：

- `POST /api/pdbzn/workflow/cluster`

输入：

- Step1 结果
- 聚类模式
- 阈值

逻辑：

1. 先把没有本地 PDB 的条目剔除
2. 再按模式聚类

三种模式：

- `preset_clusters`
- `recompute_tmalign`
- `recompute_neighbor`

输出：

- Cluster 表
- `dropped_no_pdb_count`
- `dropped_outside_cluster_count`

### Step3：第二次结构筛选

入口：

- `POST /api/pdbzn/workflow/filter`

输入：

- Step2 代表结构
- `ligand_distance`
- `metal_distance`
- `max_residue_per_oligomer`

逻辑：

1. 读取结构文件
2. 检查 ZN 附近是否有大配体
3. 检查 ZN 附近是否有其他金属
4. 检查单体长度/聚体限制
5. 检查 TriHis 是否至少满足基本要求

输出：

- 保留列表
- 剔除列表
- 每条剔除原因

### Step4：3HIS 几何复核

入口：

- `POST /api/pdbzn/workflow/validate`

输入：

- Step3 保留结果
- 是否必须通过 `require_his3`

逻辑：

1. 重新读取结构
2. 找 ZN
3. 找可配位 HIS 原子
4. 计算 3 个 HIS 到 ZN 的键长
5. 计算配位角
6. 统计 5A 内残基
7. 排除对称构象和几何异常

输出：

- `TriHis_Recheck_Pass`
- `Symmetry_Spatial_Filter_Pass`
- `TriHis_Residues`
- `ZN_Site`
- `His_ZN_Bond_Lengths_A`
- `His_ZN_Bond_Angles_Deg`
- `Residues_Within_5A`

### Step5：主表收口

入口：

- `POST /api/pdbzn/workflow/finalize`

输入：

- Step4 通过结果
- 是否运行 `fpocket`
- 是否写回主表
- 导出文件名

逻辑：

1. 以 Step4 通过结果作为候选集输入
2. 对每个候选重新解析本地结构并重算 3HIS / Zn 配位字段
3. 从本地 DiffDock 原始产物重新确定最佳 pose 并重算相关距离
4. 可选重新运行 `fpocket`，不再回填旧 pocket 汇总字段
5. 计算 `Step5_FinalScore`
6. 写出到 `backend/runtime/exports/<output>.csv`
7. 如果允许，再把新的 Step5 结果追加写回主表 CSV

输出：

- Step5 汇总表
- `backend/runtime/exports/zn_his_master_table_step5.csv`
- 可选写回 `backend/data/master_table.csv` 或外部主表

## 8. 1QRG 系列筛选在当前系统中的位置

这部分整合自 `PDB_ZN_1QRG系列筛选说明.md`。

当前这条“约 70 个候选”的复现口径是：

- Step1：`ZN + HIS:3`
- Step2：`recompute_tmalign`, `cluster_threshold=0.7`
- Step3：`ligand_distance=5.0`, `metal_distance=8.0`, `max_residue_per_oligomer=180`
- Step4：`require_his3=true`
- Step5：`run_fpocket=false`, `write_to_master=false`

按这个口径，文档记录的复现实测结果为：

- Step1：1265
- Step2 输入：1229
- Step2 输出簇：148
- Step3 保留：78
- Step4 通过：72
- Step5 appended：72

同时需要注意：

- `1QRG` 作为参考对象在主表中存在
- 但在 Step2 之前会因为“无本地 TM-align PDB”被剔除

也就是说：

- `1QRG` 在这个流程里更像“规则锚点”
- 不是当前聚类输入里一定参与计算的成员

## 9. 当前代码的现实约束与风险

这一节很重要，因为它解释了为什么有些旧文档会让人误判。

### 1. 后端是单文件架构

优点：

- 查逻辑快
- 部署简单

缺点：

- 代码边界弱
- 修改容易互相影响
- 工作流、任务调度、文件解析、HTTP 路由全部耦合在一起

### 2. 数据有静态快照和 runtime 双轨

优点：

- 前端可离线式浏览
- runtime 可保存历史任务

缺点：

- `master_table.csv`
- `master_full.json`
- `data.json`
- `pdbzn.sqlite`

这几层很容易不同步。

### 3. 对文件系统依赖非常重

当前系统强依赖：

- 目录命名
- 文件名模式
- `rank*_confidence-*.sdf`
- `*_info.txt`
- `pocket*_vert.pqr`
- `*_out/`

只要文件命名或目录结构变了，很多解析逻辑就会失效。

### 4. 外部工具依赖不是“软依赖”

静态浏览可以没有外部工具。

但一旦你要跑：

- DiffDock 在线推理
- fpocket 在线任务
- TMalign 在线任务
- Step5 的 `run_fpocket`

环境必须完整，否则 runtime 只会留下失败日志。

### 5. `backend/data` 里静态结果不一定代表当前 Step5 最新结果

当前仓库里就存在这个现象：

- `backend/data/master_table.csv` 是 121 行
- `backend/data/data.json` / `master_full.json` 是 72 条

这意味着：

- 前端静态浏览看到的“主表”
- 后端 workflow 看到的“主表”

可能不是同一个时间点的版本。

## 10. 推荐的理解方式

如果你以后要继续维护这个仓库，建议用下面这个脑图来理解。

### 第一层：前端页面

就是 6 个 HTML 工具页。

### 第二层：统一后端

就是 `backend/diffdock_api_server.py`。

### 第三层：两类数据

- `backend/data/` 静态快照
- `backend/runtime/` 动态状态

### 第四层：三类外部工具

- DiffDock
- fpocket
- TM-align

### 第五层：一条核心工作流

- Step1 初筛
- Step2 聚类
- Step3 结构规则筛选
- Step4 3HIS 几何复核
- Step5 主表收口

## 11. 最短使用建议

如果你的目标是“先跑通，再细化”，建议按这个顺序：

1. 先启动后端 `8015`
2. 再启动前端 `8020`
3. 先打开 `frontend/pdbzn_workflow.html`
4. 用 `/api/pdbzn/workflow/config` 看数据库状态和数据源
5. 跑 Step1 到 Step5
6. 再去 `master_table.html`、`index.html`、`diffdock_compare.html` 验证结果

如果你的目标是“只看已有结果，不跑重计算”，建议按这个顺序：

1. 启动前端静态服务
2. 启动后端，只用于 `api/data/file` 和健康检查
3. 先看 `master_table.html`
4. 再看 `index.html`
5. 最后看 `diffdock_compare.html`

## 12. 本文档和旧文档的关系

本文档已经吸收了原来 4 份文档中的核心内容：

- 从 `README.md` 吸收了页面、API、运行方式
- 从 `架构说明_论文素材.md` 吸收了架构表述，但修正了“Flask”和数据库描述
- 从 `PDB_ZN_1QRG系列筛选说明.md` 吸收了 1QRG 系列复现口径和 72 候选逻辑
- 从 `backend-port-conflict.md` 吸收了端口冲突排查和推荐端口策略

如果后面继续维护，建议优先维护本文档，其他历史说明文档只作为补充记录保留。
