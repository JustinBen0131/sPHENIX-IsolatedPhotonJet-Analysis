bin_PROGRAMS = tree_to_hists
BUILT_SOURCES = BuildInfo.h
nodist_tree_to_hists_SOURCES = BuildInfo.h
tree_to_hists_SOURCES = TreeToHists.cc Campaign.cc Package.cc Isolation.cc Corrections.cc
tree_to_hists_LDADD = $(ROOT_LIBS) -lTreePlayer -lTMVA -lTMVAUtils
AM_CPPFLAGS = -I$(srcdir) -DTREETOHISTS_SOURCE_DIR=\"$(abs_srcdir)\"
EXTRA_DIST = CandidateCache.h CandidateCachePackage.h MeasurementAccumulator.h CandidateMeasurementReplay.h CandidateResponseReplay.h ResponseMatching.h TreeToHists.h PhotonJet.h Preselection.h MbdTimeReference.h LeadingAssociation.h Configuration.h FillPhotonJet.C FillPhotonJetResponse.C FillPhotonQA.C autogen.sh EnergyVariation.h CentralityWeightSource.h CentralityReplay.h PhotonEnergySelection.h PhotonEnergyInputs.h NativeEnergyVariation.h PurityInputs.h SelectionReplay.h ResponseOnlyResolution.h CollectionContract.h ReplayMeasurementPackage.cc EmsComposition.h EmsNativeSelection.h NativeResponseGrid.h SelectionViewContract.h SelectionEfficiency.h SelectionFanout.cc ScientificTrace.h HistogramScientificDefinitions.h QAScientificContext.h MeasurementScientificDefinitions.h CandidateScientificContext.h ReplayScientificContext.h NativeModelDefinition.h ScoreReviewScientific.h ScoreReview.C ScientificProducerEvidence.h

# Embed the source hashes in the executable. Reading changed source files at
# runtime must never make an older binary appear to implement a newer plan.
code_files = CandidateCache.h CandidateCachePackage.h MeasurementAccumulator.h CandidateMeasurementReplay.h CandidateResponseReplay.h ResponseMatching.h TreeToHists.cc TreeToHists.h Campaign.cc Package.cc Isolation.cc Corrections.cc Configuration.h PhotonJet.h Preselection.h MbdTimeReference.h LeadingAssociation.h FillPhotonJet.C FillPhotonQA.C FillPhotonJetResponse.C EnergyVariation.h CentralityWeightSource.h CentralityReplay.h PhotonEnergySelection.h PhotonEnergyInputs.h NativeEnergyVariation.h PurityInputs.h SelectionReplay.h ResponseOnlyResolution.h CollectionContract.h EmsComposition.h EmsNativeSelection.h NativeResponseGrid.h SelectionViewContract.h SelectionEfficiency.h SelectionFanout.cc ScientificTrace.h HistogramScientificDefinitions.h QAScientificContext.h MeasurementScientificDefinitions.h CandidateScientificContext.h ReplayScientificContext.h NativeModelDefinition.h ScoreReviewScientific.h ScoreReview.C ScientificProducerEvidence.h
PYTHON_FOR_SCIENCE ?= python3
.PHONY: native_producer_witness_force
native_producer_witness_force:
BuildInfo.h: native_producer_witness_force $(srcdir)/Makefile.am $(addprefix $(srcdir)/,$(code_files)) $(srcdir)/../../TreeProduction/src/internal/Configuration.h
	@set -e; witness=$$($(PYTHON_FOR_SCIENCE) '$(srcdir)/../../FinalAnalysis/scripts/write_native_build_witness.py' --source-root '$(srcdir)/../..' --output-dir "$$(pwd)/.physics/native-source-headers"); { printf '#define TREETOHISTS_COMPILED_CODE R"tth({'; separator=''; \
	for name in $(code_files) TreeProduction/Configuration.h; do \
	  file='$(srcdir)/'$$name; test "$$name" != TreeProduction/Configuration.h || file='$(srcdir)/../../TreeProduction/src/internal/Configuration.h'; \
	  hash=$$( (sha256sum "$$file" 2>/dev/null || shasum -a 256 "$$file") | cut -d' ' -f1); \
	  printf '%s"%s":"%s"' "$$separator" "$$name" "$$hash"; separator=','; \
	done; printf '})tth"\n'; printf '#include "%s"\n' "$$witness"; } > $@.next; if ! cmp -s $@.next $@; then mv $@.next $@; fi
CLEANFILES = BuildInfo.h
