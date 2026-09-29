# The public library reuses the maintained ROOT fillers through Histograms.cc.
AM_CPPFLAGS = $(ROOT_CFLAGS) -I$(srcdir)
AM_CXXFLAGS = -std=c++17
lib_LTLIBRARIES = libPhotonJetHists.la
libPhotonJetHists_la_SOURCES = PhotonJetHists.cc internal/Selection.cc internal/Histograms.cc internal/Package.cc internal/Centrality.cc
libPhotonJetHists_la_LIBADD = $(ROOT_LIBS) -lTreePlayer
include_HEADERS = PhotonJetHists.h
noinst_HEADERS = internal/Types.h EventTreeInput.h SelectionInput.h PhotonJetSelection.h PhotonJetTruth.h
EXTRA_DIST = FillPhotonJet.C FillPhotonQA.C FillPhotonJetResponse.C FillIsolationPopulation.C SumCentrality.C
