AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

# --------------------------------------------------------------------
# Compiler / linker flags (prefix-resolved by ./configure)
# --------------------------------------------------------------------
# NOTE:
#  - Some sPHENIX external headers (CLHEP/HepMC2/CLHEP) trigger warnings under C++20.
#  - In many environments, -Werror may be injected via CXXFLAGS, causing those
#    external warnings to become fatal.
#  - Putting -Wno-error in AM_CXXFLAGS ensures it appears AFTER CXXFLAGS and
#    disables "warnings as errors" without removing your warning flags.
AM_CXXFLAGS  = -std=c++20 -Wall -Wextra -Wshadow -Wno-error

# Keep this as a belt-and-suspenders override in case upstream injects -Werror
CXXFLAGS += -Wno-error

AM_CPPFLAGS  = \
  -I$(top_srcdir)/../src \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS   = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

# --------------------------------------------------------------------
# Library to build and install (AuAu-only; different name to avoid any collision)
# --------------------------------------------------------------------
lib_LTLIBRARIES = libRecoilJetsAuAu.la

libRecoilJetsAuAu_la_SOURCES = \
  RecoilJets_AuAu.cc

# IMPORTANT:
#  - RecoilJets_AuAu.cc uses CaloRawClusterEval (from g4eval), so we must link -lg4eval
libRecoilJetsAuAu_la_LIBADD = \
  -lcalo_io \
  -lphool \
  -lcalotrigger \
  -lmbd_io \
  -lSubsysReco \
  -lcentrality_io \
  -lfun4all \
  -lg4dst \
  -lg4detectors_io \
  -lphg4hit \
  -lg4eval

# Avoid undefined symbols with libtool; give a basic version-info
libRecoilJetsAuAu_la_LDFLAGS = $(AM_LDFLAGS) $(LDFLAGS) -no-undefined -version-info 0:0:0

# Install public headers (AuAu-only)
pkginclude_HEADERS = \
  RecoilJets_AuAu.h \
  $(top_srcdir)/../src/PhotonClusterBuilder.h \
  $(top_srcdir)/../src/RJFastJetAlgoSubArea.h \
  $(top_srcdir)/../src/RJReplayFoundationV1.h \
  $(top_srcdir)/../src/RJDominantTruthWitnessV1.h \
  $(top_srcdir)/../src/RJSchema10DataRetentionContractV1.h \
  $(top_srcdir)/../src/RJCentralityReplayV1.h \
  $(top_srcdir)/../src/RJCollaboratorMbdV1.h \
  $(top_srcdir)/../src/RJPhotonTrainingViewV1.h \
  $(top_srcdir)/../src/RJReplayRuntimeV1.h \
  $(top_srcdir)/../src/RJGl1ScalerContract.h \
  $(top_srcdir)/../src/RJGl1TriggerContract.h \
  $(top_srcdir)/../src/RJTriggerRunInfoV1.h \
  $(top_srcdir)/../src/RJTriggerScalersV1.h \
  $(top_srcdir)/../src/RJShowerFactorialV1.h

# --------------------------------------------------------------------
# Small link test against the lib we build here
# --------------------------------------------------------------------
noinst_PROGRAMS = testexternals
testexternals_SOURCES = testexternals.C
testexternals_LDADD   = libRecoilJetsAuAu.la

# Generate a tiny main if it doesn't exist
testexternals.C:
	@echo "//*** this is a generated file. Do not commit, do not edit" > $@
	@echo "int main(){ return 0; }" >> $@

# Optional: ROOT CINT dictionary rule (kept from your file, if you need it)
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

# Clean extras
clean-local:
	rm -f testexternals.C
