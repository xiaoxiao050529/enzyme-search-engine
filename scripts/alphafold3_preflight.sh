#!/usr/bin/env bash

set -u
set -o pipefail

pass() { printf '[PASS] %s\n' "$*"; }
warn() { printf '[WARN] %s\n' "$*"; }
fail() { printf '[FAIL] %s\n' "$*"; }
info() { printf '[INFO] %s\n' "$*"; }

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

ARCH="$(uname -m 2>/dev/null || echo unknown)"
OS_PRETTY="$(. /etc/os-release 2>/dev/null; echo "${PRETTY_NAME:-unknown}")"
KERNEL="$(uname -r 2>/dev/null || echo unknown)"

info "AlphaFold 3 preflight"
info "Host OS: ${OS_PRETTY}"
info "Kernel: ${KERNEL}"
info "Arch: ${ARCH}"

if [[ "${ARCH}" == "x86_64" ]]; then
  pass "CPU architecture is x86_64"
else
  fail "CPU architecture is ${ARCH}; official AlphaFold 3 install docs are written for Linux with NVIDIA GPU tooling, and this repo recommends using an x86_64 GPU node"
fi

if grep -qi 'ubuntu 22.04' /etc/os-release 2>/dev/null; then
  pass "OS matches the verified Ubuntu 22.04 reference environment"
else
  warn "OS is not the verified Ubuntu 22.04 reference environment"
fi

if have_cmd nvidia-smi; then
  GPU_LINE="$(nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv,noheader 2>/dev/null | head -1)"
  if [[ -n "${GPU_LINE}" ]]; then
    pass "Detected NVIDIA GPU: ${GPU_LINE}"
  else
    warn "nvidia-smi exists but GPU query returned no rows"
  fi
else
  fail "nvidia-smi not found; no usable NVIDIA GPU stack detected"
fi

if have_cmd docker; then
  pass "docker found: $(docker --version 2>/dev/null)"
else
  warn "docker not found"
fi

if have_cmd podman; then
  warn "podman found: $(podman --version 2>/dev/null)"
fi

if have_cmd python3; then
  info "python3: $(python3 --version 2>/dev/null) ($(command -v python3))"
else
  fail "python3 not found"
fi

if have_cmd free; then
  info "Memory summary:"
  free -h
fi

if have_cmd df; then
  info "Disk summary:"
  df -h / "${HOME:-/home}" 2>/dev/null || true
fi

for var in ALPHAFOLD3_RUNNER ALPHAFOLD3_COMMAND ALPHAFOLD3_MODEL_DIR ALPHAFOLD3_DB_DIR ALPHAFOLD3_DB_DIRS; do
  val="${!var-}"
  if [[ -n "${val}" ]]; then
    pass "${var} is set: ${val}"
  else
    warn "${var} is not set"
  fi
done

cat <<'EOF'

Recommended production target for this repo:
  - Linux x86_64
  - Ubuntu 22.04
  - NVIDIA A100 80 GB or H100 80 GB
  - >= 64 GB RAM
  - ~1 TB SSD for AlphaFold 3 databases

If this host fails the GPU / architecture checks, use it as the web/API node only
and run AlphaFold 3 on a separate GPU machine, then point the site at that runner
or import completed result directories through the AlphaFold page.
EOF
