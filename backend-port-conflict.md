# Backend Port Conflict Solution

## 1. 问题现象

启动后端时出现：

```text
OSError: [Errno 98] Address already in use
```

这表示后端要绑定的监听地址已经被其他进程占用了。

在当前项目里，后端默认绑定：

- Host: `127.0.0.1`
- Port: `8015`

对应代码：

- `backend/diffdock_api_server.py`
- `HOST = os.environ.get("DIFFDOCK_API_HOST", "127.0.0.1")`
- `PORT = int(os.environ.get("DIFFDOCK_API_PORT", "8015"))`

服务启动时通过 `ThreadingHTTPServer((HOST, PORT), Handler)` 直接绑定这个地址。

## 2. 根因判断

最常见的几种情况：

- 你之前已经启动过一次后端，现在旧进程还在跑。
- 另一个终端里有人已经启动了同一个项目的后端。
- 共享机器上还有别的程序占用了 `8015`。
- 你用 `nohup`、`tmux`、`screen` 或后台 `&` 启动过服务，但忘了它还在运行。

## 3. 推荐处理策略

### 方案 A：已有服务正常，就不要重复启动

先检查 `8015` 上的现有服务是否可用：

```bash
curl http://127.0.0.1:8015/api/pdbzn/workflow/config
```

如果返回 JSON，说明后端已经在正常运行。此时不需要再执行一次：

```bash
python3 backend/diffdock_api_server.py
```

直接继续用现有服务即可。

这是最推荐的处理方式，尤其是在你刚刚已经打开过网页的情况下。

### 方案 B：需要重启同一个端口

先找到占用 `8015` 的进程：

```bash
ss -ltnp '( sport = :8015 )'
```

如果系统有 `lsof`，也可以：

```bash
lsof -iTCP:8015 -sTCP:LISTEN -n -P
```

拿到 PID 后，先确认进程是什么：

```bash
ps -fp <PID>
```

如果确认就是旧的后端进程，先优雅停止：

```bash
kill <PID>
```

等 1 到 3 秒，再确认端口是否已经释放：

```bash
ss -ltnp '( sport = :8015 )'
```

如果进程仍未退出，再强制结束：

```bash
kill -9 <PID>
```

然后重新启动：

```bash
cd /tmp/enzyme-search-engine
python3 backend/diffdock_api_server.py
```

### 方案 C：不影响现有服务，改用新端口启动

如果你不想停掉已有的 `8015` 服务，最稳妥的方法是启动到另一个端口。

这个仓库里更推荐用：

- `8017`
- `8016`

原因是多个前端页面都对这两个备用端口有兼容处理；相比之下，像 `8135` 这类端口并不是所有页面都会自动探测。

启动命令：

```bash
cd /tmp/enzyme-search-engine
DIFFDOCK_API_PORT=8017 python3 backend/diffdock_api_server.py
```

如果需要改 host：

```bash
cd /tmp/enzyme-search-engine
DIFFDOCK_API_HOST=127.0.0.1 DIFFDOCK_API_PORT=8017 python3 backend/diffdock_api_server.py
```

如果你希望其它机器也能访问，再考虑：

```bash
cd /tmp/enzyme-search-engine
DIFFDOCK_API_HOST=0.0.0.0 DIFFDOCK_API_PORT=8017 python3 backend/diffdock_api_server.py
```

只有在你明确需要远程访问时才建议使用 `0.0.0.0`。

## 4. 前端联动要求

前端页面默认把后端当成 `http://127.0.0.1:8015`。

如果你改了后端端口，需要注意：

- `pdbzn_workflow.html` 默认是 `8015`，但也会尝试一些候选端口。
- `fpocket.html`、`tmalign.html`、`diffdock_compare.html` 也默认优先找 `8015`。
- 对整个仓库来说，`8017` 和 `8016` 是更通用的备用端口。

如果页面没有自动切过去，就手动把页面里的 API 地址改成：

```text
http://127.0.0.1:8017
```

不要把前端静态页面端口 `8020` 填进 API 输入框。`8020` 是前端页面，不是后端 API。

## 5. 标准操作流程

建议按下面顺序处理：

1. 先测 `8015` 是否已经可用。
2. 如果可用，直接复用，不要重复启动。
3. 如果不可用，再查是谁占了 `8015`。
4. 如果是旧后端，就停掉后重启。
5. 如果不能停现有进程，就改用 `8017` 启动。
6. 启动成功后，再检查前端 API 地址是否一致。

## 6. 启动后的验证方法

后端正常启动时，终端会打印类似：

```text
DiffDock API listening on http://127.0.0.1:8015
```

也可以手动验证：

```bash
curl http://127.0.0.1:8015/api/pdbzn/workflow/config
```

如果改成了 `8017`，就用：

```bash
curl http://127.0.0.1:8017/api/pdbzn/workflow/config
```

浏览器侧确认：

- 前端页面正常打开
- 页面里的 API 地址和后端真实端口一致
- 点击页面按钮后不再出现连接失败

## 7. 不建议的做法

- 不要连续多次重复执行 `python3 backend/diffdock_api_server.py`
- 不要不看 PID 就直接 `killall python3`
- 不要把 `8020` 当前端端口误当成后端 API 端口
- 不要为了抢端口随便杀掉别人的进程

## 8. 最短可执行方案

如果你只想快速解决，直接按这个顺序执行：

```bash
curl http://127.0.0.1:8015/api/pdbzn/workflow/config
```

如果有返回，停止操作，直接使用现有服务。

如果没有返回：

```bash
ss -ltnp '( sport = :8015 )'
ps -fp <PID>
kill <PID>
cd /tmp/enzyme-search-engine
python3 backend/diffdock_api_server.py
```

如果你不想停旧服务：

```bash
cd /tmp/enzyme-search-engine
DIFFDOCK_API_PORT=8017 python3 backend/diffdock_api_server.py
```

然后把前端 API 地址改成：

```text
http://127.0.0.1:8017



## 启动网站的方式
```
 按这个顺序启动就行。

  1. 先看后端是不是已经在跑

  curl http://127.0.0.1:8015/api/pdbzn/workflow/config

  如果返回一段 JSON，说明后端已经启动了，不要再执行 python3 backend/diffdock_api_server.py。

  2. 如果后端没跑，再启动后端

  cd /tmp/enzyme-search-engine
  python3 backend/diffdock_api_server.py

  3. 如果又报 Address already in use，直接换端口启动

  cd /tmp/enzyme-search-engine
  DIFFDOCK_API_PORT=8017 python3 backend/diffdock_api_server.py

  4. 再启动前端静态页面

  cd /tmp/enzyme-search-engine/frontend
  python3 -m http.server 8020 --bind 127.0.0.1

  5. 浏览器打开页面

  http://127.0.0.1:8020/pdbzn_workflow.html

  如果你后端是用 8017 启动的，就把页面左上角的 API 地址改成：

  http://127.0.0.1:8017

  最稳的两套命令是：

  默认端口启动：

  cd /tmp/enzyme-search-engine
  python3 backend/diffdock_api_server.py
  cd /tmp/enzyme-search-engine/frontend
  python3 -m http.server 8020 --bind 127.0.0.1

  冲突时改端口启动：

  cd /tmp/enzyme-search-engine
  DIFFDOCK_API_PORT=8017 python3 backend/diffdock_api_server.py
  cd /tmp/enzyme-search-engine/frontend
  python3 -m http.server 8020 --bind 127.0.0.1

  如果你愿意，我可以直接再给你写一个 start.sh，以后你只需要运行一次脚本。
