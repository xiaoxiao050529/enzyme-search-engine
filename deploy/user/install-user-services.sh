#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
USER_SYSTEMD_DIR="${XDG_CONFIG_HOME:-$HOME/.config}/systemd/user"
BACKEND_TEMPLATE="$ROOT_DIR/deploy/user/enzyme-backend.service"
FRONTEND_TEMPLATE="$ROOT_DIR/deploy/user/enzyme-frontend.service"

mkdir -p "$USER_SYSTEMD_DIR"
sed "s|__ROOT_DIR__|$ROOT_DIR|g" "$BACKEND_TEMPLATE" > "$USER_SYSTEMD_DIR/enzyme-backend.service"
sed "s|__ROOT_DIR__|$ROOT_DIR|g" "$FRONTEND_TEMPLATE" > "$USER_SYSTEMD_DIR/enzyme-frontend.service"

systemctl --user daemon-reload
systemctl --user enable enzyme-backend.service
systemctl --user enable enzyme-frontend.service

printf 'Installed user services into %s\n' "$USER_SYSTEMD_DIR"
printf 'Next: %s\n' "$ROOT_DIR/deploy/user/start-user-services.sh"
