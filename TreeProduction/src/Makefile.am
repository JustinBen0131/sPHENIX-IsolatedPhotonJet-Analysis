AUTOMAKE_OPTIONS = foreign

# --------------------------------------------------------------------
# Compiler / linker flags (prefix-resolved by ./configure)
# --------------------------------------------------------------------
# NOTE:
#  - Some sPHENIX external headers (CLHEP/HepMC2) trigger warnings under C++20.
#  - In many environments, -Werror may be injected via CXXFLAGS, causing those
#    external warnings to become fatal.
#  - Putting -Wno-error in AM_CXXFLAGS ensures it appears AFTER CXXFLAGS and
#    disables "warnings as errors" without removing your warning flags.
AM_CXXFLAGS  = -std=c++20 -Wall -Wextra -Wshadow -Wno-error

CXXFLAGS += -Wno-error

AM_CPPFLAGS  = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS   = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

# --------------------------------------------------------------------
# Library to build and install
# --------------------------------------------------------------------
lib_LTLIBRARIES = libRecoilJets.la

libRecoilJets_la_SOURCES = \
  RecoilJets.cc

libRecoilJets_la_LIBADD = \
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
libRecoilJets_la_LDFLAGS = $(AM_LDFLAGS) $(LDFLAGS) -no-undefined -version-info 0:0:0

# Install public headers
pkginclude_HEADERS = \
  RecoilJets.h \
  PhotonClusterBuilder.h \
  RJFastJetAlgoSubArea.h \
  RJReplayFoundationV1.h \
  RJDominantTruthWitnessV1.h \
  RJPhotonTrainingViewV1.h \
  RJReplayRuntimeV1.h \
  RJGl1ScalerContract.h \
  RJGl1TriggerContract.h \
  RJTriggerRunInfoV1.h \
  RJTriggerScalersV1.h \
  RJCentralityReplayV1.h \
  RJCollaboratorMbdV1.h \
  RJSchema10DataRetentionContractV1.h \
  RJShowerFactorialV1.h

# --------------------------------------------------------------------
# Small link test against the lib we build here
# --------------------------------------------------------------------
noinst_PROGRAMS = testexternals
testexternals_SOURCES = testexternals.C
testexternals_LDADD   = libRecoilJets.la

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
