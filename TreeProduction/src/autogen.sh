#!/bin/sh
# Bootstrap the standard out-of-source TreeProduction build.
# Generates Autotools machinery beside these maintained sources, then invokes
# configure in the caller's build directory with the supplied installation prefix.
# Dependencies must already be provided by the selected sPHENIX environment.
set -eu

srcdir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
for tool in autoreconf autoconf automake aclocal; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        printf 'ERROR: %s is required; use a development environment with Autotools.\n' "$tool" >&2
        exit 1
    fi
done

if [ -z "${LIBTOOLIZE:-}" ] && command -v glibtoolize >/dev/null 2>&1; then
    LIBTOOLIZE=glibtoolize
    export LIBTOOLIZE
fi

(cd "$srcdir" && autoreconf --install --force)
exec "$srcdir/configure" "$@"
