dnl Configure TreeProduction against the selected installation. libEventTree
dnl needs ROOT only; the producer (default) also needs the sPHENIX release.
dnl External reconstruction modules remain supplied by the installed release;
dnl only photonjet::PhotonClusterBuilder is bundled with this package.
AC_PREREQ([2.69])
AC_INIT([TreeProduction], [0.1.0])
AC_CONFIG_SRCDIR([EventTree.cc])
AC_CONFIG_AUX_DIR([.build-aux])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_PROG_CXX
AC_PATH_PROG([PYTHON], [python3])
AS_IF([test -z "$PYTHON"], [AC_MSG_ERROR([python3 is required to bind the build sources])])
AS_IF(["$PYTHON" -c 'import yaml' >/dev/null 2>&1], [],
      [AC_MSG_ERROR([$PYTHON lacks PyYAML, needed to expand the eventTree schema; set PYTHON=/path/to/python3])])
LT_INIT([disable-static])
AC_LANG([C++])

AC_ARG_ENABLE([producer],
  [AS_HELP_STRING([--disable-producer],
                  [build only libEventTree (reading, checking and extending saved trees; no sPHENIX release)])],
  [], [enable_producer=yes])
AM_CONDITIONAL([BUILD_PRODUCER], [test "x$enable_producer" != xno])

AC_ARG_VAR([OFFLINE_MAIN], [sPHENIX installation prefix])
AC_ARG_VAR([OPT_SPHENIX], [sPHENIX external dependency prefix, if separate])
AC_ARG_VAR([ROOT_CONFIG], [path to root-config for the selected release])
AC_ARG_VAR([FASTJET_CONFIG], [path to fastjet-config for the selected release])
AC_ARG_VAR([SPHENIX_LIBS], [override the sPHENIX library list for a compatible release])

AC_PATH_PROG([ROOT_CONFIG], [root-config])
AS_IF([test -z "$ROOT_CONFIG"], [AC_MSG_ERROR([root-config is required])])
ROOT_CFLAGS=`"$ROOT_CONFIG" --cflags` || AC_MSG_ERROR([root-config --cflags failed])
ROOT_LIBS=`"$ROOT_CONFIG" --libs` || AC_MSG_ERROR([root-config --libs failed])

dependency_flags=""
dependency_ldflags=""
AS_IF([test "x$enable_producer" != xno], [
  AS_IF([test -z "$OFFLINE_MAIN" || test ! -d "$OFFLINE_MAIN/include"],
        [AC_MSG_ERROR([OFFLINE_MAIN must name the configured sPHENIX installation (or use --disable-producer)])])
  AC_PATH_PROG([FASTJET_CONFIG], [fastjet-config])
  AS_IF([test -z "$FASTJET_CONFIG"], [AC_MSG_ERROR([fastjet-config is required])])
  FASTJET_CFLAGS=`"$FASTJET_CONFIG" --cxxflags` || AC_MSG_ERROR([fastjet-config --cxxflags failed])
  FASTJET_LIBS=`"$FASTJET_CONFIG" --libs` || AC_MSG_ERROR([fastjet-config --libs failed])
  dependency_flags="-I$OFFLINE_MAIN/include $FASTJET_CFLAGS"
  dependency_ldflags="-L$OFFLINE_MAIN/lib -L$OFFLINE_MAIN/lib64"
  AS_IF([test -n "$OPT_SPHENIX"], [
    dependency_flags="$dependency_flags -I$OPT_SPHENIX/include"
    dependency_ldflags="$dependency_ldflags -L$OPT_SPHENIX/lib -L$OPT_SPHENIX/lib64"
  ])
])
AS_IF([test "x$prefix" != xNONE], [
  dependency_flags="-I$prefix/include $dependency_flags"
  dependency_ldflags="-L$prefix/lib $dependency_ldflags"
])
CPPFLAGS="$dependency_flags $ROOT_CFLAGS $CPPFLAGS"
LDFLAGS="$dependency_ldflags $LDFLAGS"
dnl Respect the standard selected by ROOT or the caller (some releases use C++20).
AS_CASE([" $ROOT_CFLAGS $CXXFLAGS "], [*" -std="*], [],
        [CXXFLAGS="$CXXFLAGS -std=c++17"])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if __cplusplus < 201703L
#error TreeProduction requires C++17 or newer
#endif
]], [[]])], [], [AC_MSG_ERROR([select C++17 or newer, consistent with ROOT])])
AC_CHECK_HEADERS([TFile.h TTree.h nlohmann/json.hpp], [],
                 [AC_MSG_ERROR([ROOT or nlohmann/json headers are missing from the selected installation])])
AC_SEARCH_LIBS([dladdr], [dl])

PRODUCER_LIBS=""
AS_IF([test "x$enable_producer" != xno], [
  AS_IF([test -z "$SPHENIX_LIBS"], [
    SPHENIX_LIBS="-lcalo_reco -lCaloStatusSkimmer -lcalotrigger -lcalotrigger_io -lcentrality -lcentrality_io -lglobalvertex -lglobalvertex_io -lmbd -lmbd_io -lzdcinfo -lzdcinfo_io -ljetbackground -ljetbackground_io -ljetbase -ljetbase_io -lg4jets -lg4eval -lg4mbd -lCaloWaveformSim -lg4detectors -lphg4hit -lphhepmc -lphhepmc_io -lcalo_io -lffamodules -lffaobjects -lffarawobjects -lcdbobjects -lfun4all -lSubsysReco -lphool"
  ])
  PRODUCER_LIBS="$SPHENIX_LIBS $FASTJET_LIBS -lHepMC -lCLHEP -lTMVA -lTMVAUtils $ROOT_LIBS"
  AC_CHECK_HEADERS([fun4all/Fun4AllServer.h calobase/PhotonClusterv1.h caloreco/RawClusterBuilderTopo.h g4mbd/MbdDigitization.h calowaveformsim/CaloWaveformSim.h TMVA/RBDT.hxx], [],
                   [AC_MSG_ERROR([required production header missing or incompatible; check the selected release])])
  saved_libs="$LIBS"
  LIBS="$PRODUCER_LIBS $LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <fun4all/SubsysReco.h>
#include <nlohmann/json.hpp>
struct BuildCheck : SubsysReco { BuildCheck() : SubsysReco("build_check") {} };
]], [[BuildCheck module;
auto config = nlohmann::json::parse("{\"schema_version\":1}");
return config.is_object() ? 0 : 1;]])], [],
                [AC_MSG_ERROR([production dependencies do not link; inspect config.log])])
  LIBS="$saved_libs"
])
AC_SUBST([ROOT_LIBS])
AC_SUBST([PRODUCER_LIBS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
