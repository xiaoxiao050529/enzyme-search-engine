#!/usr/bin/env bash

set -u
set -o pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR" || exit 1

PYTHON_BIN="${PYTHON_BIN:-python3}"
CURL_BIN="${CURL_BIN:-curl}"

PUBLIC_MODE="${PUBLIC_MODE:-0}"
LOCAL_HOST="127.0.0.1"
if [[ "$PUBLIC_MODE" == "1" ]]; then
  DEFAULT_BIND_HOST="0.0.0.0"
else
  DEFAULT_BIND_HOST="127.0.0.1"
fi
BACKEND_BIND_HOST="${DIFFDOCK_API_HOST:-$DEFAULT_BIND_HOST}"
FRONTEND_BIND_HOST="${FRONTEND_BIND_HOST:-$DEFAULT_BIND_HOST}"
PUBLIC_HOST="${PUBLIC_HOST:-}"
REQUESTED_BACKEND_PORT="${DIFFDOCK_API_PORT:-8015}"
REQUESTED_FRONTEND_PORT="${FRONTEND_PORT:-8020}"

LOG_DIR="$ROOT_DIR/backend/runtime/logs"
PID_DIR="$ROOT_DIR/backend/runtime/pids"
mkdir -p "$LOG_DIR" "$PID_DIR"

info() {
  printf '[INFO] %s\n' "$*"
}

warn() {
  printf '[WARN] %s\n' "$*" >&2
}

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "missing command: $1"
}

http_ok() {
  "$CURL_BIN" -fsS --max-time 2 "$1" >/dev/null 2>&1
}

port_is_listening() {
  local port="$1"
  if command -v ss >/dev/null 2>&1; then
    ss -ltnH "( sport = :$port )" 2>/dev/null | grep -q .
    return $?
  fi
  if command -v lsof >/dev/null 2>&1; then
    lsof -nP -iTCP:"$port" -sTCP:LISTEN >/dev/null 2>&1
    return $?
  fi
  return 1
}

wait_for_url() {
  local url="$1"
  local tries="${2:-40}"
  local delay="${3:-0.25}"
  local i
  for ((i = 0; i < tries; i++)); do
    if http_ok "$url"; then
      return 0
    fi
    sleep "$delay"
  done
  return 1
}

unique_ports() {
  local -n out_ref="$1"
  shift
  local port
  declare -A seen=()
  out_ref=()
  for port in "$@"; do
    [[ -n "$port" ]] || continue
    [[ -z "${seen[$port]+x}" ]] || continue
    seen["$port"]=1
    out_ref+=("$port")
  done
}

show_recent_log() {
  local log_path="$1"
  [[ -f "$log_path" ]] || return 0
  warn "last log lines from $log_path"
  tail -n 20 "$log_path" >&2 || true
}

start_backend() {
  local port endpoint log_path pid_path pid
  local candidates=()

  unique_ports candidates "$REQUESTED_BACKEND_PORT" 8015 8017 8016

  for port in "${candidates[@]}"; do
    endpoint="http://${LOCAL_HOST}:${port}/api/pdbzn/workflow/config"
    if [[ "$PUBLIC_MODE" != "1" ]] && http_ok "$endpoint"; then
      BACKEND_PORT="$port"
      BACKEND_URL="http://${LOCAL_HOST}:${port}"
      BACKEND_STATE="reused"
      BACKEND_LOG_PATH=""
      return 0
    fi

    if port_is_listening "$port"; then
      warn "backend port $port is occupied and not responding as this API, trying next port"
      continue
    fi

    log_path="$LOG_DIR/backend-${port}.log"
    pid_path="$PID_DIR/backend-${port}.pid"
    info "starting backend on ${BACKEND_BIND_HOST}:${port}"
    nohup env DIFFDOCK_API_HOST="$BACKEND_BIND_HOST" DIFFDOCK_API_PORT="$port" \
      "$PYTHON_BIN" backend/diffdock_api_server.py >"$log_path" 2>&1 &
    pid=$!
    echo "$pid" >"$pid_path"
    disown "$pid" >/dev/null 2>&1 || true

    if wait_for_url "$endpoint" 60 0.25; then
      BACKEND_PORT="$port"
      BACKEND_URL="http://${LOCAL_HOST}:${port}"
      BACKEND_STATE="started"
      BACKEND_LOG_PATH="$log_path"
      BACKEND_PID="$pid"
      return 0
    fi

    warn "backend failed to become ready on port $port"
    show_recent_log "$log_path"
  done

  return 1
}

start_frontend() {
  local port page_url log_path pid_path pid
  local candidates=()

  unique_ports candidates "$REQUESTED_FRONTEND_PORT" 8020 8030 8031

  for port in "${candidates[@]}"; do
    page_url="http://${LOCAL_HOST}:${port}/frontend/master_table.html"
    if [[ "$PUBLIC_MODE" != "1" ]] && http_ok "$page_url"; then
      FRONTEND_PORT="$port"
      FRONTEND_URL="http://${LOCAL_HOST}:${port}"
      FRONTEND_STATE="reused"
      FRONTEND_LOG_PATH=""
      return 0
    fi

    if port_is_listening "$port"; then
      warn "frontend port $port is occupied and not serving this project, trying next port"
      continue
    fi

    log_path="$LOG_DIR/frontend-${port}.log"
    pid_path="$PID_DIR/frontend-${port}.pid"
    info "starting frontend static server on ${FRONTEND_BIND_HOST}:${port}"
    nohup "$PYTHON_BIN" -m http.server "$port" --bind "$FRONTEND_BIND_HOST" \
      >"$log_path" 2>&1 &
    pid=$!
    echo "$pid" >"$pid_path"
    disown "$pid" >/dev/null 2>&1 || true

    if wait_for_url "$page_url" 60 0.25; then
      FRONTEND_PORT="$port"
      FRONTEND_URL="http://${LOCAL_HOST}:${port}"
      FRONTEND_STATE="started"
      FRONTEND_LOG_PATH="$log_path"
      FRONTEND_PID="$pid"
      return 0
    fi

    warn "frontend failed to become ready on port $port"
    show_recent_log "$log_path"
  done

  return 1
}

require_cmd "$PYTHON_BIN"
require_cmd "$CURL_BIN"

if [[ -z "$PUBLIC_HOST" && "$BACKEND_BIND_HOST" != "127.0.0.1" ]]; then
  if command -v hostname >/dev/null 2>&1; then
    PUBLIC_HOST="$(hostname -I 2>/dev/null | awk '{print $1}')"
  fi
fi

BACKEND_PORT=""
BACKEND_URL=""
BACKEND_STATE=""
BACKEND_LOG_PATH=""
BACKEND_PID=""

FRONTEND_PORT=""
FRONTEND_URL=""
FRONTEND_STATE=""
FRONTEND_LOG_PATH=""
FRONTEND_PID=""

start_backend || die "unable to start or reuse backend; see README.md"
start_frontend || die "unable to start or reuse frontend static server"

printf '\n'
printf 'Services ready.\n'
printf 'Backend API: %s\n' "$BACKEND_URL"
printf 'Frontend:    %s/frontend/master_table.html\n' "$FRONTEND_URL"
printf 'Workflow:    %s/frontend/pdbzn_workflow.html\n' "$FRONTEND_URL"
printf 'CAVER:       %s/frontend/caver.html\n' "$FRONTEND_URL"
printf 'AlphaFold:   %s/frontend/alphafold.html\n' "$FRONTEND_URL"
printf 'Home:        %s/frontend/index.html\n' "$FRONTEND_URL"
printf 'Master:      %s/frontend/master_table.html\n' "$FRONTEND_URL"
printf '\n'

if [[ "$BACKEND_STATE" == "started" ]]; then
  printf 'Backend started by script. pid=%s log=%s\n' "$BACKEND_PID" "$BACKEND_LOG_PATH"
else
  printf 'Backend reused existing service on port %s\n' "$BACKEND_PORT"
fi

if [[ "$FRONTEND_STATE" == "started" ]]; then
  printf 'Frontend started by script. pid=%s log=%s\n' "$FRONTEND_PID" "$FRONTEND_LOG_PATH"
else
  printf 'Frontend reused existing service on port %s\n' "$FRONTEND_PORT"
fi

if [[ "$BACKEND_PORT" != "8015" ]]; then
  printf '\n'
  printf 'Note: backend is running on %s instead of 8015.\n' "$BACKEND_PORT"
  printf 'If a page does not auto-detect it, set API base to %s\n' "$BACKEND_URL"
fi

if [[ "$BACKEND_BIND_HOST" != "127.0.0.1" || "$FRONTEND_BIND_HOST" != "127.0.0.1" ]]; then
  public_name="${PUBLIC_HOST:-<server-ip-or-domain>}"
  printf '\n'
  printf 'Public mode hints:\n'
  printf '  Frontend: http://%s:%s/frontend/master_table.html\n' "$public_name" "$FRONTEND_PORT"
  printf '  Workflow: http://%s:%s/frontend/pdbzn_workflow.html\n' "$public_name" "$FRONTEND_PORT"
  printf '  API:      http://%s:%s/api/pdbzn/workflow/config\n' "$public_name" "$BACKEND_PORT"
fi
