#!/usr/bin/env bash
# Public operational entrypoint; base files are read-only, output is a sidecar.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python="${PHOTONID_PYTHON:-python3}"
if [[ -z "${PHOTONID_PYTHON:-}" && -x "$script_dir/.venv/bin/python" ]]; then
  python="$script_dir/.venv/bin/python"
fi
exec "$python" "$script_dir/frontend.py" augment "$@"
