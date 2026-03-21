# Data Platform 框架说明（前端 / 后端）

这个项目是一个围绕金属酶（以 Zn 相关蛋白为重点）的分析与筛选平台，包含：

- 可视化与交互页面（前端，多个 HTML 工具页）
- 统一 API 服务（后端，单进程 HTTP 服务）
- 工作流式筛选链路（Step1~Step5）
- 外部工具集成（DiffDock / fpocket / TM-align）

---

## 1. 总体架构

项目采用“静态前端 + 本地后端 API + 本地数据/工具”的结构：

- 前端：`/root/data-platform/frontend/*.html`
  - 负责参数输入、任务提交、状态轮询、结果展示、3D 可视化
- 后端：`/root/data-platform/backend/diffdock_api_server.py`
  - 负责业务逻辑、文件解析、任务编排、外部工具调用、API 输出
- 数据目录：
  - `backend/data/`：表格与中间数据（如 `data.json`、`pocket.csv` 等）
  - `/root/PDB_ZN`：PDB-ZN 原始结构数据（CIF/PDB）、聚类参考文件等
  - `backend/runtime/`：运行时任务目录与 SQLite 数据库

核心后端入口与基础路径定义见：
- [diffdock_api_server.py](file:///root/data-platform/backend/diffdock_api_server.py#L21-L33)

---

## 2. 前端在做什么

前端是多页面工具集，页面之间通过顶部导航互跳，统一调用后端 API。

### 2.1 主要页面职责

- `index.html`：主页面，加载 ranking/master 数据并进行 3D 结构交互展示  
  [index.html](file:///root/data-platform/frontend/index.html#L36-L43)
- `diffdock_compare.html`：DiffDock 对比与在线推理提交  
  [diffdock_compare.html](file:///root/data-platform/frontend/diffdock_compare.html#L39-L67)
- `fpocket.html`：提交 fpocket 任务、查看 pocket 参数与可视化  
  [fpocket.html](file:///root/data-platform/frontend/fpocket.html#L45-L61)
- `tmalign.html`：单任务与 1 对多 TM-align 评估  
  [tmalign.html](file:///root/data-platform/frontend/tmalign.html#L39-L76)
- `master_table.html`：主表筛选、排序、跳转  
  [master_table.html](file:///root/data-platform/frontend/master_table.html#L32-L44)
- `pdbzn_workflow.html`：PDB-ZN 工作流主控（Step1~Step5）  
  [pdbzn_workflow.html](file:///root/data-platform/frontend/pdbzn_workflow.html#L69-L116)

### 2.2 前端的“工作流编排”职责

`pdbzn_workflow.html` 主要做三件事：

- 采集参数并构造 payload（Step1/2/3/4/5）
- 调用后端 `/run`、`/cluster`、`/filter`、`/validate`、`/finalize`
- 渲染每一步的 `columns + rows + summary`

关键逻辑位置：
- API 地址候选与自动探测：  
  [pdbzn_workflow.html](file:///root/data-platform/frontend/pdbzn_workflow.html#L166-L171)
- Step2 payload（含聚类模式）：  
  [pdbzn_workflow.html](file:///root/data-platform/frontend/pdbzn_workflow.html#L372-L381)
- Step2 执行（聚类）：  
  [pdbzn_workflow.html](file:///root/data-platform/frontend/pdbzn_workflow.html#L536-L593)
- Step3 执行（二次筛选）：  
  [pdbzn_workflow.html](file:///root/data-platform/frontend/pdbzn_workflow.html#L594-L622)

---

## 3. 后端在做什么

后端是一个基于 `ThreadingHTTPServer` 的统一 API 服务，负责：

- 路由分发（GET/POST）
- 本地任务管理（DiffDock/fpocket/TM-align）
- 工作流计算（Step1~Step5）
- 文件读写与结果输出（CSV/SQLite）

路由总入口：
- [diffdock_api_server.py](file:///root/data-platform/backend/diffdock_api_server.py#L2124-L2460)

服务启动入口：
- [diffdock_api_server.py](file:///root/data-platform/backend/diffdock_api_server.py#L2462-L2469)

---

## 4. PDB-ZN 工作流（Step1~Step5）

默认参数定义：
- [_pdbzn_workflow_defaults](file:///root/data-platform/backend/diffdock_api_server.py#L414-L444)

### Step1：初筛（`/api/pdbzn/workflow/run`）
- 从数据库中按金属离子、残基要求、名称关键字等过滤
- 若数据库未准备好，会自动触发导入
- 实现入口：  
  [_pdbzn_run_workflow](file:///root/data-platform/backend/diffdock_api_server.py#L1101-L1218)

### Step2：聚类（`/api/pdbzn/workflow/cluster`）
- 支持三种模式：
  - `preset_clusters`：按 `optimized_tmalign_clusters1.csv` 映射
  - `recompute_tmalign`：用 TM-align 重新聚类
  - `recompute_neighbor`：按邻域特征相似度聚类
- 入口：  
  [_pdbzn_cluster_workflow](file:///root/data-platform/backend/diffdock_api_server.py#L1396-L1496)

与模式相关的关键实现：
- 预设聚类映射：  
  [_pdbzn_cluster_rows_from_step1_clusters](file:///root/data-platform/backend/diffdock_api_server.py#L582-L677)
- 邻域相似度聚类：  
  [_pdbzn_cluster_rows](file:///root/data-platform/backend/diffdock_api_server.py#L517-L579)
- TM-align 聚类：  
  [_pdbzn_cluster_rows_tmalign](file:///root/data-platform/backend/diffdock_api_server.py#L872-L944)

当前逻辑里，Step2 会在三种模式下统一剔除没有本地 PDB 结构的条目，并在 summary 输出 `dropped_no_pdb_count` / `dropped_no_pdb_members`。

### Step3：第二次筛选（`/api/pdbzn/workflow/filter`）
- 在 Step2 的簇结果上进行结构规则筛选
- 典型规则包括：大配体距离、其他金属干扰、长度/聚体限制
- 入口：  
  [_pdbzn_step3_filter_workflow](file:///root/data-platform/backend/diffdock_api_server.py#L1314-L1419)

### Step4：综合验证（`/api/pdbzn/workflow/validate`）
- 仅做 3HIS 几何复核（不再做二次聚类）
- 输出 5A 范围内全部残基、3 个 HIS 与 ZN 的键长和键角
- 对对称/空间构象误配进行剔除并记录原因
- 入口：  
  [_pdbzn_step4_validate_workflow](file:///root/data-platform/backend/diffdock_api_server.py#L1607-L1745)

### Step5：收口与主表写回（`/api/pdbzn/workflow/finalize`）
- 把 Step4 结果与主表/DiffDock/fpocket 信息整合
- 同步写入 Step4 的几何字段（TriHis、ZN位点、键长/键角、5A残基、空间构象过滤结果）
- 可选触发 fpocket 快速补全，并可写回主表
- 入口：  
  [_pdbzn_step5_finalize_workflow](file:///root/data-platform/backend/diffdock_api_server.py#L1905-L2085)

---

## 5. 任务型 API（DiffDock / fpocket / TM-align）

后端还提供独立任务接口（提交、查状态、看日志、取文件）：

- DiffDock：
  - 提交：`POST /api/diffdock/submit`
  - 状态：`GET /api/diffdock/status/{job_id}`
  - 日志：`GET /api/diffdock/log/{job_id}`
- fpocket：
  - 提交：`POST /api/fpocket/submit`
  - 状态：`GET /api/fpocket/status/{job_id}`
  - 口袋详情：`GET /api/fpocket/pockets/{job_id}`
- TM-align：
  - 提交：`POST /api/tmalign/submit`
  - 状态：`GET /api/tmalign/status/{job_id}`
  - 日志：`GET /api/tmalign/log/{job_id}`

路由集中在：
- [diffdock_api_server.py](file:///root/data-platform/backend/diffdock_api_server.py#L2124-L2460)

---

## 6. 前后端协作方式（一次完整链路）

以 `pdbzn_workflow.html` 为例：

1. 前端读取参数并调用 `workflow/config` 获取默认配置
2. 用户触发 Step1/2/3/4/5，前端逐步发起对应 POST
3. 后端返回标准结构：`steps[] + summary`
4. 前端按表格渲染结果并显示摘要（总数、代表条目等）

这个设计的好处是：页面可视化逻辑与计算逻辑解耦，后续你调整阈值、模式、筛选规则，只需改后端或 payload，不需要重写可视化框架。

---

## 7. 本地运行（最小方式）

后端：

```bash
python /root/data-platform/backend/diffdock_api_server.py
```

前端（任意静态文件服务器）：

```bash
cd /root/data-platform/frontend
python -m http.server 8020
```

浏览器打开：

- `http://127.0.0.1:8020/pdbzn_workflow.html`
- 或 `http://127.0.0.1:8020/index.html`

注意：前端页面端口与后端 API 端口不同，API 地址需要填后端端口（例如 8130 / 8015）。

---

## 8. 页面功能说明书（主页面 / fpocket / DiffDock / TM-align / 工作流 / Master Table）

这一节专门说明每个页面的定位和区别。

### 8.1 主页面（`index.html`）

定位：总览 + 快速可视化入口。

- 展示 ranking / master 数据并支持按 PDB 查询
- 进行 3D 结构可视化（蛋白、配体、金属位点）
- 适合做“先看结果，再决定去哪个专项页面深挖”

入口参考：  
[index.html](file:///root/data-platform/frontend/index.html)

### 8.2 DiffDock 对比页面（`diffdock_compare.html`）

定位：对接结果比对 + 在线推理。

- 左右双视图对比不同 PDB 或不同 pose
- 支持在线提交 DiffDock 推理并轮询任务状态
- 支持按 ZN 距离自动挑选 top pose 进行快速比较

入口参考：  
[diffdock_compare.html](file:///root/data-platform/frontend/diffdock_compare.html)

### 8.3 fpocket 页面（`fpocket.html`）

定位：口袋识别与口袋参数分析。

- 上传 PDB 提交 fpocket 任务
- 查看 pocket 参数表（如 druggability、score、volume）
- 在 3D 视图中可视化 pocket 位置

入口参考：  
[fpocket.html](file:///root/data-platform/frontend/fpocket.html)

### 8.4 TM-align 页面（`tmalign.html`）

定位：结构相似性评估。

- 单任务：1 对 1 结构比对（TM-score、RMSD 等）
- 批任务：1 对多比对并排序
- 支持不同评分归一方式（short_chain / tm1 / tm2 / tm_max）

入口参考：  
[tmalign.html](file:///root/data-platform/frontend/tmalign.html)

### 8.5 工作流页面（`pdbzn_workflow.html`）

定位：全链路筛选主控台（Step1~Step5）。

- Step1：按金属离子和残基规则初筛
- Step2：聚类（并统一剔除无 PDB 条目）
- Step3：二次结构规则筛选
- Step4：3HIS 几何复核（含 5A 残基、键长、键角、构象剔除）
- Step5：结果收口，写回主表并补全 DiffDock/fpocket

入口参考：  
[pdbzn_workflow.html](file:///root/data-platform/frontend/pdbzn_workflow.html)

### 8.6 Master Table 页面（`master_table.html`）

定位：主结果表检索、排序、人工复核。

- 全局搜索 + 指定列筛选
- 任意列排序（升/降序）
- 作为最终结果导出前核查页面

入口参考：  
[master_table.html](file:///root/data-platform/frontend/master_table.html)

### 8.7 页面怎么选（实操建议）

- 想快速看整体候选：主页面 + Master Table
- 想看口袋质量：fpocket
- 想看对接 pose 细节：DiffDock 对比
- 想做结构相似性分群/比对：TM-align
- 想跑完整筛选链路并沉淀主表：工作流页面
