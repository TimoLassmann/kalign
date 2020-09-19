AC_DEFUN([AX_ENABLE_DEBUG],
[AC_ARG_ENABLE(debugging,
[AS_HELP_STRING([--enable-debugging],[include debugging code])
AS_HELP_STRING([--enable-debugging=x],[also set diagnostics verbosity level to <x> (1-3)])],
enable_debugging=$enableval, enable_debugging="no")

case $enable_debugging in
yes)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging. ]);;
1)  AC_DEFINE(DEBUGLEVEL, 1,[Defines debugging level 1.]);;
2)  AC_DEFINE(DEBUGLEVEL, 2,[Defines debugging level 2.]);;
3)  AC_DEFINE(DEBUGLEVEL, 3,[Defines debugging level 3.]);;
no)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging.]);;
*)  AC_MSG_ERROR([Unknown argument to --enable-debugging: $enable_debugging]);;
esac

changequote({,})
  CFLAGS=`echo "$CFLAGS" | $SED -e 's/-O[0-9s]*//g'`
  CFLAGS=`echo "$CFLAGS" | $SED -e 's/-g[0-9]*//g'`
changequote([,])


if test "$enable_debugging" != "no"; then
AC_DEFINE(DEBUG,1,[Defines debugging .])
ADD_DEBUG_COMPILE_WARNINGS
CFLAGS="-ggdb -std=gnu11 ${CFLAGS}"
CFLAGS="${CFLAGS}"
else
ADD_PRODUCTION_COMPILE_WARNINGS
CFLAGS="-O3 -std=gnu11 ${CFLAGS}"
CFLAGS="${CFLAGS}"
DEBUG=0
fi

])

AC_DEFUN([ADD_DEBUG_COMPILE_WARNINGS],
  [TLDEVEL_CFLAGS=""
  AX_CHECK_COMPILE_FLAG([-pedantic],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -pedantic"],,)
  AX_CHECK_COMPILE_FLAG([-fstack-protector],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -fstack-protector"],,)
  dnl AX_CHECK_COMPILE_FLAG([-fsanitize=safe-stack],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -fsanitize=safe-stack"],,)
  AX_CHECK_COMPILE_FLAG([-fstack-clash-protection],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -fstack-clash-protection"],,)
  AX_CHECK_COMPILE_FLAG([-D_GLIBCXX_ASSERTIONS],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -D_GLIBCXX_ASSERTIONS"],,)
  AX_CHECK_COMPILE_FLAG([-fexceptions],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -fexceptions"],,)
  AX_CHECK_COMPILE_FLAG([-fasynchronous-unwind-tables],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -fasynchronous-unwind-tables"],,)

  AX_CHECK_COMPILE_FLAG([-Waddress],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Waddress"],,)
  AX_CHECK_COMPILE_FLAG([-Wall],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wall"],,)
  AX_CHECK_COMPILE_FLAG([-Wformat],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wformat"],,)
  AX_CHECK_COMPILE_FLAG([-Wformat-security],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wformat-security"],,)
  AX_CHECK_COMPILE_FLAG([-Werror=format-security],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Werror=format-security"],,)
  AX_CHECK_COMPILE_FLAG([-Werror=implicit-function-declaration],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Werror=implicit-function-declaration"],,)
  AX_CHECK_COMPILE_FLAG([-Wswitch-default],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wswitch-default"],,)
  AX_CHECK_COMPILE_FLAG([-Wswitch-enum],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wswitch-enum"],,)
  AX_CHECK_COMPILE_FLAG([-Wswitch],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wswitch"],,)
  AX_CHECK_COMPILE_FLAG([-Wtrigraphs],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wtrigraphs"],,)
  AX_CHECK_COMPILE_FLAG([-Wtype-limits],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wtype-limits"],,)
  AX_CHECK_COMPILE_FLAG([-Wundef],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wundef"],,)
  AX_CHECK_COMPILE_FLAG([-Wuninitialized],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wuninitialized"],,)
  AX_CHECK_COMPILE_FLAG([-Wunsafe-loop-optimizations],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wunsafe-loop-optimizations"],,)
  AX_CHECK_COMPILE_FLAG([-Wextra],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wextra"],,)
  AX_CHECK_COMPILE_FLAG([-Wfloat-equal],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -Wfloat-equal"],,)

  AX_CHECK_COMPILE_FLAG([-ffunction-sections],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -ffunction-sections"],,)
  AX_CHECK_COMPILE_FLAG([-fdata-sections],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -fdata-sections"],,)
  AX_CHECK_COMPILE_FLAG([--gc-sections],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} --gc-sections"],,)

  AC_SUBST([TLDEVEL_CFLAGS], [${TLDEVEL_CFLAGS}])
])

AC_DEFUN([ADD_PRODUCTION_COMPILE_WARNINGS],
    [TLDEVEL_CFLAGS=""

    AX_CHECK_COMPILE_FLAG([-D_FORTIFY_SOURCE=2],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -D_FORTIFY_SOURCE=2"],,)
    AX_CHECK_COMPILE_FLAG([-ffunction-sections],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -ffunction-sections"],,)
    AX_CHECK_COMPILE_FLAG([-fdata-sections],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} -fdata-sections"],,)
    AX_CHECK_COMPILE_FLAG([--gc-sections],[TLDEVEL_CFLAGS="${TLDEVEL_CFLAGS} --gc-sections"],,)
    AC_SUBST([TLDEVEL_CFLAGS], [${TLDEVEL_CFLAGS}])
])
