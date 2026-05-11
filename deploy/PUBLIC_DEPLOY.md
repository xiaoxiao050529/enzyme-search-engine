# 公网部署

推荐方式是：

1. 后端仅监听 `127.0.0.1:8015`
2. 前端仅监听 `127.0.0.1:8020`
3. 用 `nginx` 对外暴露 `80/443`
4. 将 `/api/` 反代到后端，其余路径反代到前端

这样比直接暴露 `8015/8020` 更稳，也更容易后续接 HTTPS。

## 一次性快速启动

```bash
PUBLIC_MODE=1 ./start.sh
```

这会让前后端直接绑定 `0.0.0.0`，适合临时测试公网访问。

## 无 sudo 用户级常驻

如果你明确不想使用 `sudo`，请直接看：

- `deploy/user/USER_SYSTEMD.md`

这套方案使用 `systemctl --user`，把服务安装到当前用户目录下。

## 直接公网常驻

如果你暂时不准备上 `nginx`，可以直接用下面这组 `systemd` 文件：

- `deploy/enzyme-backend-public.service`
- `deploy/enzyme-frontend-public.service`

默认端口是：

- 前端：`8040`
- 后端：`8017`

对应访问地址：

- `http://你的IP:8040/frontend/master_table.html`
- `http://你的IP:8040/frontend/pdbzn_workflow.html`

## 推荐长期部署

将以下文件复制到系统目录：

- `deploy/enzyme-backend.service` -> `/etc/systemd/system/enzyme-backend.service`
- `deploy/enzyme-frontend.service` -> `/etc/systemd/system/enzyme-frontend.service`
- `deploy/nginx.enzyme-search-engine.conf` -> `/etc/nginx/conf.d/enzyme-search-engine.conf`

然后执行：

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now enzyme-backend
sudo systemctl enable --now enzyme-frontend
sudo nginx -t
sudo systemctl reload nginx
```

访问地址：

- `http://你的域名/frontend/master_table.html`
- `http://你的域名/frontend/pdbzn_workflow.html`
- `http://你的域名/frontend/index.html`

## HTTPS

如果有域名，建议继续加证书。常见方式是 `certbot` 或 `acme.sh`。

## 防火墙

长期部署时只需要对外开放：

- `80`
- `443`

不建议直接开放：

- `8015`
- `8020`
