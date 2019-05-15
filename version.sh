#!/usr/bin/env bash
VERSION=
inside_git_repo=

inside_git_repo="$(git rev-parse --is-inside-work-tree 2>/dev/null)"

if [ "$inside_git_repo" ]; then

    VERSION=`git describe --always`
else

    if [ -f VERSION ]; then
        if [[ ! $VERSION ]]; then
            VERSION=`cat VERSION`
        fi
    else
        VERSION="[na]"
    fi
fi

echo $VERSION

