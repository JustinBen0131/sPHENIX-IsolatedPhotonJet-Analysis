AC_INIT([TreeToHists], [1.0])
AC_CONFIG_SRCDIR([TreeToHists.cc])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_PATH_PROG([ROOT_CONFIG], [root-config])
AS_IF([test -z "$ROOT_CONFIG"], [AC_MSG_ERROR([load the ROOT or sPHENIX environment first])])
ROOT_LIBDIR=`"$ROOT_CONFIG" --libdir` || AC_MSG_ERROR([root-config --libdir failed])
LDFLAGS="-Wl,-rpath,$ROOT_LIBDIR $LDFLAGS"
AC_PROG_CXX
AC_LANG([C++])
ROOT_CFLAGS=`"$ROOT_CONFIG" --cflags` || AC_MSG_ERROR([root-config --cflags failed])
ROOT_LIBS=`"$ROOT_CONFIG" --libs` || AC_MSG_ERROR([root-config --libs failed])
CPPFLAGS="$ROOT_CFLAGS $CPPFLAGS"
AS_CASE([" $ROOT_CFLAGS $CXXFLAGS "], [*" -std="*], [], [CXXFLAGS="$CXXFLAGS -std=c++17"])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if __cplusplus < 201703L
#error C++17 or newer is required
#endif
]], [[]])], [], [AC_MSG_ERROR([select C++17 or newer, consistent with ROOT])])
AC_CHECK_HEADERS([TFile.h TTreeReader.h yaml.h nlohmann/json.hpp], [],
                 [AC_MSG_ERROR([ROOT, libyaml and nlohmann/json headers are required])])
AC_SEARCH_LIBS([yaml_parser_initialize], [yaml], [], [AC_MSG_ERROR([libyaml is required])])
AC_SUBST([ROOT_LIBS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
