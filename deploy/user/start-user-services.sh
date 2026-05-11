#!/usr/bin/env bash

set -euo pipefail

systemctl --user daemon-reload
systemctl --user restart enzyme-backend.service
systemctl --user restart enzyme-frontend.service

printf 'Started user services.\n'
printf 'Frontend: http://%s:8040/frontend/master_table.html\n' "$(hostname -I 2>/dev/null | awk '{print $1}')"
printf 'Workflow: http://%s:8040/frontend/pdbzn_workflow.html\n' "$(hostname -I 2>/dev/null | awk '{print $1}')"
printf 'API: http://%s:8017/api/pdbzn/workflow/config\n' "$(hostname -I 2>/dev/null | awk '{print $1}')"
