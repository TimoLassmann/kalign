#!/usr/bin/env sh

exists()
{
  command -v "$1" >/dev/null 2>&1
}

if ! exists libtoolize; then
    if ! exists glibtoolize; then
        echo 'libtool/glibtool not found!'
        echo 'to install using homebrew:'
        echo 'brew install libtool'
        exit 1
    fi
fi

if ! exists aclocal; then
    echo 'aclocal not found!'
    echo 'to install using homebrew:'
    echo 'brew install automake'
    exit 1
fi

if ! exists autoheader; then
    echo 'autoheader not found!'
    echo 'to install using homebrew:'
    echo 'brew install automake'
    exit 1
fi

if ! exists automake; then
    echo 'automake not found!'
    echo 'to install using homebrew:'
    echo 'brew install automake'
    exit 1
fi

if ! exists autoconf; then
    echo 'autoconf not found!'
    echo 'to install using homebrew:'
    echo 'brew install automake'
    exit 1
fi



test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

cd "$srcdir"

case `uname` in Darwin*) glibtoolize --force --copy ;;
                 *) libtoolize --force  --copy ;;
esac

aclocal -I m4 $AL_OPTS
autoheader $AH_OPTS
automake --add-missing --copy --gnu $AM_OPTS
autoconf $AC_OPTS


echo
echo "Now run '$srcdir/configure' and 'make' to compile."
echo

