#!/usr/bin/env bash
# Public operational entrypoint; all training semantics live in the Python path.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python="${PHOTONID_PYTHON:-python3}"
if [[ -z "${PHOTONID_PYTHON:-}" && -x "$script_dir/.venv/bin/python" ]]; then
  python="$script_dir/.venv/bin/python"
fi
exec "$python" "$script_dir/frontend.py" train "$@"
