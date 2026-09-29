AC_INIT([PhotonJetHists], [1.0], [])
AC_CONFIG_SRCDIR([PhotonJetHists.cc])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_PROG_CXX
LT_INIT
AC_PATH_PROG([ROOT_CONFIG], [root-config], [no])
AS_IF([test "x$ROOT_CONFIG" = xno], [AC_MSG_ERROR([ROOT and root-config are required])])
ROOT_CFLAGS=`$ROOT_CONFIG --cflags`
ROOT_LIBS=`$ROOT_CONFIG --libs`
AC_SUBST([ROOT_CFLAGS])
AC_SUBST([ROOT_LIBS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
