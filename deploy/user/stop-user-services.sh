#!/usr/bin/env bash

set -euo pipefail

systemctl --user stop enzyme-backend.service || true
systemctl --user stop enzyme-frontend.service || true

printf 'Stopped user services.\n'
