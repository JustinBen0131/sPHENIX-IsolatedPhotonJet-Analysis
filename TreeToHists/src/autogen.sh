#!/bin/sh
# Generate the package build files. Run from an empty build directory.
set -eu
srcdir=$(CDPATH= cd -- "${0%/*}" && pwd)
autoreconf --install "$srcdir"
exec "$srcdir/configure" "$@"
