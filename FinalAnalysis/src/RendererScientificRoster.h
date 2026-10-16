/* Explicit typed-annotation roster for the maintained drawing entrypoints.
   This declaration is not admission: qualify_plot_physics must execute every
   entrypoint and the consumer must match the exact tested producer pinset. */
#pragma once
#include <set>
#include <string>
namespace photonjet_science {
inline bool annotationRosterDeclared(const std::string& family) {
  static const std::set<std::string> families={
    "ian.atlasOverlay",
    "ian.auauCentralities",
    "ian.beforeAfter",
    "ian.centralityOverlay",
    "ian.combinatoric",
    "ian.detectorLevel",
    "ian.iterationOverlay",
    "ian.iterationRule",
    "ian.jointRefold",
    "ian.jointResponse",
    "ian.pairClosure",
    "ian.photonClosure",
    "ian.photonJetDphi",
    "ian.photonResponse",
    "ian.photonYield",
    "ian.ppVsCentral",
    "ian.priorReweighting",
    "ian.purityCorrection",
    "ian.radiusOverlay",
    "ian.recoTruthAssociation",
    "ian.statisticalBudget",
    "ian.systematicSources",
    "ian.systematicVariations",
    "qa.abcdCounts",
    "qa.bdtPerClass",
    "qa.depthDiagnostics",
    "qa.drawEtaTimeMap",
    "qa.e11Flow",
    "qa.efficiency",
    "qa.etaPhiMap",
    "qa.eventDistribution",
    "qa.featureDiagnostics",
    "qa.isolationOverlays",
    "qa.isolationStacks",
    "qa.minimumBias",
    "qa.ncbScore",
    "qa.ncbTimeSelections",
    "qa.ppg12ModelComparison",
    "qa.purityPanels",
    "qa.reweightingClosure",
    "qa.routingDiagnostics",
    "qa.savedScoreReview",
    "qa.savedWorkingPointFits",
    "qa.scoreSeparation",
    "qa.showerShapes",
    "qa.stitching",
    "qa.trainingReweighting",
    "qa.triggerStudy",
    "qa.truthIsolation",
    "qa.truthPtAfterIsolation",
    "qa.widthClassFigure",
    "qa.widthDataPanels",
    "qa.widthSummary",
    "qa.workingPointCentralityComparison",
    "qa.workingPoints"
  };
  return families.count(family)!=0;
}
}
