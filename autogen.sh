#!/usr/bin/env sh

test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.	

cd "$srcdir"
git submodule update --init --recursive


autoreconf --force --install --verbose "$srcdir"
test -n "$NOCONFIGURE" || "$srcdir/configure" "$@"
