#!/bin/sh
set -eu
srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.
mkdir -p "$srcdir/m4"

(cd $srcdir; aclocal -I ${OFFLINE_MAIN}/share;\
libtoolize --force; automake -a --add-missing; autoconf)

$srcdir/configure "$@"
