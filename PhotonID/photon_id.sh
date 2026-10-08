#!/usr/bin/env bash
# PhotonID public entrance: select an existing Python environment and forward
# the requested training/scoring command to the shared command coordinator.
# No environment activation, dependency installation or scheduler submission.
set -euo pipefail
PHOTONID_HOME="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${_CONDOR_SCRATCH_DIR:-}" ]]; then
    # Keep Python/library temporaries off the shared analysis filesystem.
    export TMPDIR="$_CONDOR_SCRATCH_DIR"
    export XDG_CACHE_HOME="$_CONDOR_SCRATCH_DIR/photonid-cache"
fi
exec "${PHOTONID_PYTHON:-python3}" "$PHOTONID_HOME/scripts/photon_id.py" "$@"
