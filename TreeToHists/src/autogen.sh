#!/bin/sh
# Generate the conventional configure script and Makefile templates.
set -eu
cd "$(dirname "$0")"
autoreconf --install --force
