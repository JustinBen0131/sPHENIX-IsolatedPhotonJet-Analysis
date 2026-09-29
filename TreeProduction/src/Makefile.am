# Build one library with an explicit source list and no local-test dependency.
# Package-qualified headers keep the bundled builder separate from CaloReco.
lib_LTLIBRARIES = libPhotonJetTree.la

libPhotonJetTree_la_SOURCES = \
  PhotonJetTree.cc \
  PhotonClusterBuilder.cc \
  internal/Production.cc \
  internal/Setup.cc \
  internal/Event.cc \
  internal/Objects.cc \
  internal/Truth.cc \
  internal/Output.cc

AM_CPPFLAGS = -I$(builddir) -I$(srcdir)
# Preserve the staged binary32 arithmetic used by the reference production.
# Fused multiply-add can move a saved shower moment by one float ULP.
AM_CXXFLAGS = -ffp-contract=off
libPhotonJetTree_la_LDFLAGS = -no-undefined -version-info 0:0:0

productionincludedir = $(includedir)/TreeProduction
productioninclude_HEADERS = PhotonJetTree.h PhotonClusterBuilder.h
internalincludedir = $(productionincludedir)/internal
internalinclude_HEADERS = internal/Production.h internal/Types.h

EXTRA_DIST = autogen.sh ../scripts/eventtree/writer.hxx

# Regenerate when a compiled source, installed header or build recipe changes.
# The digest describes the actual compilation inputs, including local edits.
build_identity_sources = $(libPhotonJetTree_la_SOURCES) $(productioninclude_HEADERS) $(internalinclude_HEADERS) Makefile.am configure.ac autogen.sh ../scripts/eventtree/writer.hxx
BUILT_SOURCES = BuildIdentity.h
CLEANFILES = BuildIdentity.h BuildIdentity.h.tmp
BuildIdentity.h: $(build_identity_sources)
	$(AM_V_GEN)$(PYTHON) -c 'import hashlib,pathlib,sys; r=pathlib.Path(sys.argv[1]); h=hashlib.sha256(); [(h.update(n.encode()+b"\0"),h.update((r/n).read_bytes()),h.update(b"\0")) for n in sorted(sys.argv[2:])]; print("#define PHOTONJET_SOURCE_SHA256 \""+h.hexdigest()+"\"")' "$(srcdir)" $(build_identity_sources) > $@.tmp
	mv $@.tmp $@
