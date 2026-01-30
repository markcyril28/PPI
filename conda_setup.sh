#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ENV_NAME:-PPI}"
PY_VERSION="${PY_VERSION:-3.9}"
FOLDX_DIR="${FOLDX_DIR:-$HOME/software/foldx}"
FOLDX_BIN="$FOLDX_DIR/foldx"
ENV_HINT_FILE="${ENV_HINT_FILE:-$PWD/mutatex.env}"

usage() {
	cat <<EOF
MutateX Conda bootstrapper
Usage: ./conda_setup.sh [--help]

Environment variables:
  ENV_NAME      Target conda environment name (default: $ENV_NAME)
  PY_VERSION    Python version to install (default: $PY_VERSION)
  FOLDX_DIR     Directory for FoldX binaries (default: $FOLDX_DIR)
  ENV_HINT_FILE File where exports for mutatex.sh are written (default: $ENV_HINT_FILE)
  SOLVER        conda/mamba selector (auto-detects mamba if available)
EOF
}

log_step() {
	printf '\n[%s] %s\n' "$(date +'%H:%M:%S')" "$1"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
	usage
	exit 0
fi

SOLVER="${SOLVER:-conda}"
if command -v mamba >/dev/null 2>&1; then
	SOLVER="mamba"
fi
log_step "Using $SOLVER to create env '$ENV_NAME' (Python $PY_VERSION)"

if ! "$SOLVER" create -y -n "$ENV_NAME" "python=$PY_VERSION"; then
	log_step "$SOLVER failed, falling back to conda"
	conda create -y -n "$ENV_NAME" "python=$PY_VERSION"
fi

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

log_step "Installing MutateX and Python dependencies"
pip install --upgrade pip
pip install \
	git+https://github.com/ELELAB/mutatex.git \
	pyyaml biopython pandas numpy scipy matplotlib

log_step "Configuring FoldX"
mkdir -p "$FOLDX_DIR"
chmod +x "$FOLDX_BIN" 2>/dev/null || true
export FOLDX_BIN="$FOLDX_BIN"
export PATH="$PATH:$FOLDX_DIR"

log_step "Sanity-checking MutateX import"
python - <<'EOF'
import mutatex
print("MutateX OK:", mutatex.__file__)
EOF

ENV_PREFIX="$CONDA_BASE/envs/$ENV_NAME"
MUTATEX_BIN_PATH="$ENV_PREFIX/bin/mutatex"
PYTHON_BIN_PATH="$ENV_PREFIX/bin/python"

log_step "Writing MutateX exports to $ENV_HINT_FILE"
cat >"$ENV_HINT_FILE" <<EOF
export MUTATEX_BIN="$MUTATEX_BIN_PATH"
export PYTHON_BIN="$PYTHON_BIN_PATH"
EOF

log_step "Setup complete. Activate with: conda activate $ENV_NAME"
log_step "mutatex.sh will automatically read $ENV_HINT_FILE to reuse this environment"
