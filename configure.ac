
AC_INIT([kalign], [m4_esyscmd_s([./version.sh])])

#AC_CONFIG_AUX_DIR([.])

AC_USE_SYSTEM_EXTENSIONS
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB
AM_SILENT_RULES([no])



#AC_CONFIG_MACRO_DIR([m4])



# M4 macros for checking of CPU features.
m4_include([m4/ax_gcc_x86_avx_xgetbv.m4])
m4_include([m4/ax_gcc_x86_cpuid.m4])
m4_include([m4/ax_check_compile_flag.m4])
m4_include([m4/ax_ext.m4])
m4_include([m4/ax_openmp.m4])
AX_EXT



AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp],[use OpenMP directives for parallelism])], enable_openmp=$enableval, enable_openmp=no)

if test "$enable_openmp" = "yes"; then
   AC_DEFINE(HAVE_OPENMP,1,[Define to enable OpenMP])
   AX_OPENMP([], [AC_MSG_ERROR([don't know how to enable OpenMP])])
fi
AC_SUBST(SIMD_FLAGS)


AC_ARG_WITH([openmp],
    AS_HELP_STRING([--without-openmp], [Do not use OpenMP]))
AS_IF([test "x$with_openmp" = "xno"], 
            [], 
            [AX_OPENMP(
                AC_DEFINE([HAVE_OPENMP],[1],[Defined if OpenMP should and can be used])
             )]
     )
AC_SUBST(OPENMP_CFLAGS)
# according to http://beagle-lib.googlecode.com/svn-history/r824/trunk/configure.ac:
# "OpenMP checker only defines for C when compiling both C and C++"
OPENMP_CXXFLAGS=$OPENMP_CFLAGS
AC_SUBST(OPENMP_CXXFLAGS)


# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL

# test for working C11 _Generic support
AN_IDENTIFIER([_Generic], [AC_C__GENERIC])
AC_DEFUN([AC_C__GENERIC],
[AC_CACHE_CHECK([for _Generic], ac_cv_c__Generic,
[AC_COMPILE_IFELSE(
   [AC_LANG_SOURCE(
      [[int
        main (int argc, char **argv)
        {
          int a = _Generic (argc, int: argc = 1);
          int *b = &_Generic (argc, default: argc);
          char ***c = _Generic (argv, int: argc, default: argv ? &argv : 0);
          _Generic (1 ? 0 : b, int: a, default: b) = &argc;
          _Generic (a = 1, default: a) = 3;
          return a + !b + !c;
        }
      ]])],
   [ac_cv_c__Generic=yes],
   [ac_cv_c__Generic=no])])
if test $ac_cv_c__Generic = no; then
AC_MSG_ERROR([We could not compile a simple C11 example. Install a new(er) compiler version (e.g.  GCC 4.9).])
fi
])# AC_C__GENERIC



AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([pthread], [main])


# Checks for header files.
AC_CHECK_HEADERS([math.h float.h stdlib.h unistd.h])
AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_HEADERS([limits.h])


AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRNLEN
AC_CHECK_FUNCS([gettimeofday pow realpath sqrt strstr])


tldevel_check_f=$ac_abs_confdir/tldevel/README

AC_CHECK_FILE([$tldevel_check_f], [],[have_tldevel=no])


if test "x${have_tldevel}" = xno; then
echo "------------------------------------------"
echo " Unable to find pthreads on this system.  "
echo " Building a single-threaded version.      "
echo "------------------------------------------"
AC_MSG_ERROR([ASFSAF])
fi

LIB_TLDEVELDIR="tldevel"
AC_SUBST(LIB_TLDEVELDIR)



AC_CANONICAL_HOST
case $host_os in
darwin* )
DEAD_STRIP="-Wl,-dead_strip"
;;
*)
DEAD_STRIP="-Wl,--gc-sections -Wl,--as-needed"
;;
esac
AC_SUBST(DEAD_STRIP)


AC_ARG_ENABLE(debugging,[AS_HELP_STRING([--enable-debugging],[include debugging code])
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

if test "$enable_debugging" != "no"; then
AC_DEFINE(DEBUG,1,[Defines debugging .])
CFLAGS="-ggdb -pedantic -Wall -std=gnu11 -fdata-sections -ffunction-sections"
else
CFLAGS="-O3 -pedantic -funroll-loops -Wall -std=gnu11 -fdata-sections -ffunction-sections"
DEBUG=0
fi

ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"

AC_ARG_ENABLE(valgrind-tests,
[--enable-valgrind-tests run self tests under valgrind],
opt_valgrind_tests="$enableval", opt_valgrind_tests="no")


OPTS="-q --trace-children=no --error-exitcode=1 --leak-check=yes"
# Run self-tests under valgrind?
if test "$opt_valgrind_tests" = "yes"; then
AC_DEFINE(DEBUG,1,[Defines debugging .])
CFLAGS="-ggdb -pedantic -Wall -m64 -std=gnu11"
AC_CHECK_PROGS(VALGRIND, valgrind)
VALGRIND="$VALGRIND $OPTS"
else
VALGRIND=
fi
AC_MSG_CHECKING([whether self tests are run under valgrind])
AC_MSG_RESULT($opt_valgrind_tests)

ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"

AC_CONFIG_FILES([Makefile
src/Makefile dev/Makefile])

AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_SUBDIRS([tldevel])

AC_OUTPUT

AC_MSG_NOTICE([---------------------------------------------------------])
AC_MSG_NOTICE([Configuration complete - $PACKAGE_NAME-$PACKAGE_VERSION])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([General Information:])
AC_MSG_NOTICE([               Version: $VERSION])
AC_MSG_NOTICE([         Configured on: $DATE])
AC_MSG_NOTICE([         Configured by: $USER])
AC_MSG_NOTICE([           Host system: $host_cpu-$host_vendor-$host_os])
AC_MSG_NOTICE([     Uname information: $UNAME_INFO])
AC_MSG_NOTICE([    Installation point: $prefix])
AC_MSG_NOTICE([            SIMD flags: $SIMD_FLAGS])
AC_MSG_NOTICE([---------------------------------------------------------])
