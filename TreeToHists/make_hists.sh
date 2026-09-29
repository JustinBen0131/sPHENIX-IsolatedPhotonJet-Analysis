#!/bin/sh
# Public entry point; all physics calculations run in the C++ ROOT package.
set -eu
package=$(CDPATH= cd -- "${0%/*}" && pwd)
exec "${PHOTONJET_PYTHON:-python3}" -B "$package/scripts/make_hists.py" "$@"
