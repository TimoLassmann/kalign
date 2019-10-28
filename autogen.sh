#!/usr/bin/env sh

exists()
{
  command -v "$1" >/dev/null 2>&1
}

if ! exists libtoolize; then
    echo 'libtool not found!'
    exit 1
fi

if ! exists aclocal; then
    echo 'aclocal not found!'
    exit 1
fi

if ! exists autoheader; then
    echo 'autoheader not found!'
    exit 1
fi

if ! exists automake; then
    echo 'automake not found!'
    exit 1
fi

if ! exists autoconf; then
    echo 'autoconf not found!'
    exit 1
fi


test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

cd "$srcdir"

libtoolize --force --copy
aclocal -I m4 $AL_OPTS
autoheader $AH_OPTS
automake --add-missing --copy --gnu $AM_OPTS
autoconf $AC_OPTS


echo
echo "Now run '$srcdir/configure' and 'make' to compile."
echo

