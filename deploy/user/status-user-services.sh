#!/usr/bin/env bash

set -euo pipefail

systemctl --user --no-pager --full status enzyme-backend.service
systemctl --user --no-pager --full status enzyme-frontend.service
