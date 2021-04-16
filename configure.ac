AC_INIT(kalign, 3.3.1)

#AC_CONFIG_AUX_DIR([.])

AC_USE_SYSTEM_EXTENSIONS
AM_INIT_AUTOMAKE
# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_RANLIB
AM_SILENT_RULES([no])


LT_INIT

AC_CONFIG_MACRO_DIR([m4])

m4_include([m4/ax_ext.m4])
m4_include([m4/ax_valgrind_check.m4])
m4_include([m4/ax_check_compile_flag.m4])
m4_include([m4/ax_set_compile_flags.m4])

AX_EXT

AX_ENABLE_DEBUG

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
echo "-------------------------------------------"
echo " Unable to compile C11 code on this system!"
echo "-------------------------------------------"
echo ""
echo "Install a new(er) version of GCC (e.g. GCC 4.9+)"
echo "or clang!"
echo ""
AC_MSG_ERROR([ERROR: C compiler too old!])
AC_MSG_ERROR([We could not compile a simple C11 example. Install a new(er) compiler version (e.g.  GCC 4.9).])
fi
])# AC_C__GENERIC


# Check for C11 _Generic support
AC_C__GENERIC



AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([pthread], [main])
AC_SEARCH_LIBS([clock_gettime],[rt posix])

# Checks for header files.
AC_CHECK_HEADERS([math.h float.h stdlib.h unistd.h])
AC_CHECK_HEADERS([sys/time.h])

AC_CHECK_HEADERS([mach/mach.h])

AC_CHECK_HEADERS([limits.h])


AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

AC_TYPE_INT8_T
AC_TYPE_INT32_T
AC_TYPE_INT16_T
AC_TYPE_INT64_T

AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRNLEN
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([gettimeofday pow realpath sqrt strstr memset])
AC_CHECK_FUNCS([clock_gettime])
AC_CHECK_FUNCS([gethrtime])
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([localtime_r])
AC_CHECK_FUNCS([strcspn])

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


AX_VALGRIND_DFLT([memcheck], [on])
AX_VALGRIND_CHECK


AC_CONFIG_FILES([Makefile
src/Makefile dev/Makefile])

AC_CONFIG_HEADER([src/config.h])


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
AC_MSG_NOTICE([                CFLAGS: $CFLAGS])
AC_MSG_NOTICE([               CPFLAGS: $CPPFLAGS   $TLDEVEL_CPPFLAGS])
AC_MSG_NOTICE([               LDFLAGS: $LDFLAGS])
AC_MSG_NOTICE([            SIMD flags: $SIMD_FLAGS])
AC_MSG_NOTICE([---------------------------------------------------------])
