#!/bin/sh
# Public entry point for measurement, systematic combinations and saved-result rendering.
set -eu
package=$(CDPATH= cd -- "${0%/*}" && pwd)
exec "${PHOTONJET_PYTHON:-python3}" -B "$package/scripts/final_analysis.py" "$@"
