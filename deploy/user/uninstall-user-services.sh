#!/usr/bin/env bash

set -euo pipefail

USER_SYSTEMD_DIR="${XDG_CONFIG_HOME:-$HOME/.config}/systemd/user"

systemctl --user disable --now enzyme-backend.service || true
systemctl --user disable --now enzyme-frontend.service || true
rm -f "$USER_SYSTEMD_DIR/enzyme-backend.service"
rm -f "$USER_SYSTEMD_DIR/enzyme-frontend.service"
systemctl --user daemon-reload

printf 'Removed user services from %s\n' "$USER_SYSTEMD_DIR"
