# AlphaFold 3 集成说明

## 当前结论

这个仓库已经接入了 AlphaFold 页面和后端 API：

- 页面：`/frontend/alphafold.html`
- API：`/api/alphafold/ping`
- 提交：`/api/alphafold/submit`
- 导入已有结果：`/api/alphafold/import`

但是“真实 AlphaFold 3 推理”是否能在当前主机运行，取决于底层硬件和官方环境。

在本次检查中，当前机器特征是：

- 架构：`aarch64`（Kunpeng 920）
- 系统：`Kylin Linux Advanced Server V10`
- GPU：未检测到 `nvidia-smi`
- 容器：无 `docker`，只有较老的 `podman`

这意味着：

- 当前机器适合继续作为网站/API 节点
- 不适合继续作为 AlphaFold 3 真实推理节点

## 推荐部署拓扑

推荐拆成两台机器：

1. 当前主机
   - 负责前端页面、后端 API、结果管理、导入和展示

2. 独立 GPU 推理主机
   - 负责真实 AlphaFold 3 推理
   - 推荐环境：
     - `x86_64`
     - `Ubuntu 22.04`
     - `NVIDIA A100 80GB` 或 `H100 80GB`
     - `>= 64 GB RAM`
     - `~1 TB SSD` 数据库空间

## 仓库当前支持的两种接法

### 方案 A：本地 runner 模式

如果 GPU 主机上直接运行这个仓库的后端，配置以下环境变量即可：

```bash
export ALPHAFOLD3_RUNNER=/path/to/alphafold3/run_alphafold.py
export ALPHAFOLD3_PYTHON_BIN=/path/to/python
export ALPHAFOLD3_MODEL_DIR=/path/to/models
export ALPHAFOLD3_DB_DIR=/path/to/databases
```

如果官方 runner 的参数或包装脚本和仓库默认探测逻辑不一致，可以直接覆盖命令：

```bash
export ALPHAFOLD3_COMMAND='/path/to/python /path/to/run_alphafold.py --json_path={input_json} --output_dir={output_dir} --model_dir={model_dir} --db_dir={db_dir}'
```

后端启动后，页面会自动走：

- `POST /api/alphafold/submit`
- `GET /api/alphafold/status/:job_id`
- `GET /api/alphafold/result/:job_id`

### 方案 B：结果导入模式

如果 AlphaFold 3 在另一台机器或人工流程里运行，当前网站仍然可用。

在 `AlphaFold` 页面里填写：

- 结果目录路径

然后调用：

```text
POST /api/alphafold/import
```

后端会复制结果目录，并在站内：

- 解析 `.cif/.mmcif/.pdb`
- 展示结构
- 提供结果文件下载
- 尝试读取置信度 JSON

## 预检脚本

仓库提供了一个快速预检：

```bash
bash scripts/alphafold3_preflight.sh
```

它会检查：

- 架构
- OS
- NVIDIA GPU
- docker/podman
- Python
- AlphaFold 3 相关环境变量

## 当前建议

如果目标是“尽快让网站具备真实 AF3 预测能力”，不要继续在当前这台 `aarch64 + 无 NVIDIA GPU` 主机上硬装。

更合理的下一步是：

1. 准备一台符合要求的 `x86_64 + NVIDIA GPU` 主机
2. 在那台机器安装官方 AlphaFold 3 环境和数据库
3. 把当前仓库后端也部署到那台机器，或把当前站点改造成调用那台机器的 AlphaFold runner
4. 当前主机继续保留网页、数据管理和结果展示能力
