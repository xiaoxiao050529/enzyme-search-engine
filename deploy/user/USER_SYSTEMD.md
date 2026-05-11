# 用户级 systemd 方案

这套方案不需要 `sudo`。

适用场景：

- 你只想把网页常驻在当前账号下
- 你可以接受服务依赖当前用户的 `systemd --user`
- 你不要求 `80/443`

默认端口：

- 前端：`8040`
- 后端：`8017`

访问地址：

- `http://你的服务器IP:8040/frontend/master_table.html`
- `http://你的服务器IP:8040/frontend/pdbzn_workflow.html`

## 安装

```bash
bash deploy/user/install-user-services.sh
```

## 启动

```bash
bash deploy/user/start-user-services.sh
```

## 查看状态

```bash
bash deploy/user/status-user-services.sh
```

## 停止

```bash
bash deploy/user/stop-user-services.sh
```

## 卸载

```bash
bash deploy/user/uninstall-user-services.sh
```

## 限制

1. 这是用户级服务，不是系统级服务。
2. 很多机器上它只能在你登录后运行。
3. 如果机器不允许用户服务跨登出保活，就不能做到真正“重启后自动对外提供服务”。
4. 不使用 `sudo` 时，通常也不能配置 `nginx`、`80/443`、系统防火墙。
