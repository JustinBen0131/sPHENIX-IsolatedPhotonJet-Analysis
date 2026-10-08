bin_PROGRAMS = tree_to_hists
BUILT_SOURCES = BuildInfo.h
nodist_tree_to_hists_SOURCES = BuildInfo.h
tree_to_hists_SOURCES = TreeToHists.cc Campaign.cc Package.cc Isolation.cc Corrections.cc
tree_to_hists_LDADD = $(ROOT_LIBS) -lTreePlayer -lTMVA -lTMVAUtils
AM_CPPFLAGS = -I$(srcdir) -DTREETOHISTS_SOURCE_DIR=\"$(abs_srcdir)\"
EXTRA_DIST = TreeToHists.h PhotonJet.h Configuration.h FillPhotonJet.C FillPhotonJetResponse.C FillPhotonQA.C autogen.sh

# Embed the source hashes in the executable. Reading changed source files at
# runtime must never make an older binary appear to implement a newer plan.
code_files = TreeToHists.cc TreeToHists.h Campaign.cc Package.cc Isolation.cc Corrections.cc Configuration.h PhotonJet.h FillPhotonJet.C FillPhotonQA.C FillPhotonJetResponse.C
BuildInfo.h: $(addprefix $(srcdir)/,$(code_files)) $(srcdir)/../../TreeProduction/src/internal/Configuration.h
	@set -e; { printf '#define TREETOHISTS_COMPILED_CODE R"tth({'; separator=''; \
	for name in $(code_files) TreeProduction/Configuration.h; do \
	  file='$(srcdir)/'$$name; test "$$name" != TreeProduction/Configuration.h || file='$(srcdir)/../../TreeProduction/src/internal/Configuration.h'; \
	  hash=$$( (sha256sum "$$file" 2>/dev/null || shasum -a 256 "$$file") | cut -d' ' -f1); \
	  printf '%s"%s":"%s"' "$$separator" "$$name" "$$hash"; separator=','; \
	done; printf '})tth"\n'; } > $@.next; mv $@.next $@
CLEANFILES = BuildInfo.h
