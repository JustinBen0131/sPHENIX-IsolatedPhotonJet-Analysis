/* IANFigures: every FinalAnalysis figure, in sPHENIX style and the IAN's layout, under the IAN's own file names. */

/* ./FinalAnalysis/final_analysis.sh draw --results results.root --qa qa.root --fig FinalAnalysis/fig */

/* The IAN is the internal analysis note whose figures this file reproduces. The measurement figures draw results.root;
   the photon and event QA figures and the saved-model figures (QAFigures.C, included before the entry point) draw the
   QA file. The output folder mirrors the IAN's fig/ folder.

   A figure draws saved quantities. What it works out itself is for display: projections, rebinning, normalisations
   and ratios, and a few derived quantities whose functions say "computed here" (simulatedCombinatoric, detectorLevel,
   systematicSources, recoTruthAssociation, iterationRule, combinatoric, jointRefold, photonYield). The results, their
   labels, their photon bins and their reported xJ ranges come from results.root. Every xJ figure shows only the
   reported range, where the distribution is fully populated (results.root derives it from the jet threshold and the
   lowest photon pT); its axis starts at 0 and ends at the range's upper edge rounded up to 0.1.

   A few places name a collision system or class and need a change for another set of classes: centralityOverlay (p+p,
   50-80% and 0-20%), the ATLAS comparison of the entry point (p+p and the class that starts at 0%), classStyle (four
   Au+Au styles), closureLines (the simulation of each package), radiusOverlay (the file names of the other radii) and
   the entry point (recoTruthAssociation for p+p, the radius figure for p+p and the class that starts at 0%, the file
   names of the IAN).

   Layout. A single-frame figure has the sPHENIX mark and its physics lines at the top left inside the frame and the
   legend at the top right, or under the header when the two do not fit side by side, and a y range opened until no
   drawn point reaches either block (ianFrame). The figures of the results set the legend right-aligned (resultFrame).
   The response maps carry the mark above the frame, the two grids of frames a boxed legend in their second frame, and
   the linear iteration panel its legend at the lower right. The physics lines are compiled from the saved metadata
   and histograms where these hold them; the photon identification, the names of the simulations, the kinds of closure
   and the names of the systematic sources are typed in the functions that write them.

   Where each figure of the note is drawn (the entry point at the end of the file holds the file names):

     purityAndCombSub/  purity correction ....................... purityCorrection
                        combinatoric subtraction ................ combinatoric
     Unfolding/         photon response, joint response ......... photonResponse, jointResponse
                        iteration stability (photons, joint) .... photonIterationRule, jointIterationRule
                        photon closure, photon yield ............ photonClosure, photonYield
                        prior reweighting ....................... priorReweighting
                        xJ before and after unfolding ........... beforeAfter
                        xJ closure, refolding ................... pairClosure, jointRefold
                        p+p against the most central Au+Au ...... ppVsCentral
     LHCComparison/     ATLAS overlay ........................... atlasOverlay
     Systematics/       uncertainty by source ................... systematicSources
     xJ_generalQA/      DATA against simulation, detector level . detectorLevel
                        reconstructed against truth-tagged ...... recoTruthAssociation
                        photon-jet Delta-phi grid ............... photonJetDphi
                        unfolded xJ for three jet radii ......... radiusOverlay
     FinalResults/      p+p with two Au+Au classes .............. centralityOverlay
                        the Au+Au classes ....................... auauCentralities
     the QA folders     ......................................... QAFigures.C (qaFigures, modelQAFigures, bdtPerClass)

   The sections of this file, in order: the IAN look (palette, sizes); the state of a run; inputs (metadata, settings,
   results, histograms); physics lines; text, header and legend; the frame; ratio pads and axis text; saving; series;
   the xJ frame and its marks; the corrections, xJ before and after unfolding and the checks of the simulation; the
   unfolded results; responses; unfolding stability and closure; the Delta-phi grid; the entry point. */
#include "sPhenixStyle.C"
#include "../src/RenderMetadata.h"
#include "../src/PlotPhysics.h"
#include "../src/RendererScientificRoster.h"
#include "../src/NativeProducerEvidence.h"
#include "../src/RootScientificTrace.h"
#include "../src/RenderScientificParameters.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TFrame.h>
#include <TPave.h>
#include <TObjString.h>
#include <TNamed.h>
#include <TPad.h>
#include <TPave.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <functional>
#include <set>
#include <algorithm>
#include <array>
#include <cctype>

namespace {

/* One result of results.root.
   key: its directory in the file. package: PP or AuAu. label: the collision system and centrality as a header writes
   them. energyLine: the collision energy as a header line (empty without metadata). fileTag, longFileTag: the two forms
   of the class in the IAN's file names (pp; auau020 and auau_0to20). centrality: the saved name of the centrality class
   (cent0_20 for Au+Au 0-20%). centralityLow: the lower edge of the Au+Au class in percent, -1 for p+p and for a class the file does not name
   by its edges. centralRank:
   the number of Au+Au classes of the file that are more central (classStyle). xjLow, xjHigh: the reported xJ range.
   jetThreshold: the recoil-jet pT threshold in GeV. photonIterations: the iterations of the photon unfolding. */
struct Result {
  std::string key, package, fileTag, longFileTag, label, energyLine, centrality;
  double xjLow = 0, xjHigh = 0, jetThreshold = 0, centralityLow = -1;
  int centralRank = 0, photonIterations = 0;
  std::string combinatoric;  /* none, simulation, or data_dphi_<shape>_<scope> */
  double sidebandLow = 0, sidebandHigh = 0;  /* the Delta-phi sideband of a data-driven combinatoric estimate, rad */
};

/* ---- the IAN look ------------------------------------------------------------------------------------------------
   The palette. Each figure gives the colours the roles its IAN original gives them; classStyle fixes the colour of a
   collision class for every figure of the results. ianPalette() sets the seven colours before any figure is drawn
   (QAFigures.C names them cBlue, cOrange, ... by reference), so nothing at file scope may copy one. A few figures add
   colours of their own where they use them. Text is sized in pixels, so every figure and every pad reads alike. */
int gBlue = 0, gVermilion = 0, gRed = 0, gGreen = 0, gOrange = 0, gPurple = 0, gGray = 0;
void ianPalette() {
  gBlue = TColor::GetColor("#0072B2"), gVermilion = TColor::GetColor("#D55E00"), gRed = TColor::GetColor("#C62828");
  gGreen = TColor::GetColor("#2CA02C"), gOrange = TColor::GetColor("#FF7F0E"), gPurple = TColor::GetColor("#6A3D9A");
  gGray = TColor::GetColor("#6C7686");
}
/* A colour for series number i of a figure with many series (the saved-model figures of QAFigures.C). */
int seriesColor(int i) {
  const int c[] = {kBlack, gBlue, gVermilion, gRed, gGreen, gOrange, gPurple, gGray};
  return c[i % 8];
}
/* The font of every text: Helvetica, sized in pixels. */
const int kFont = 43;
/* The sPHENIX mark, the header and legend lines, the small lines and the pitch of the lines, in pixels. While an
   image is written, save() draws each text at the exact size of its font step (exactTextSizes). */
const double kMarkPx = 28, kLinePx = 21.5, kSmallPx = 19.5, kStepPx = 32;
const double kInsetPx = 28;  /* from the frame to the header: clear of the axis ticks */
/* The layout numbers that the frames share, in pixels: from the frame's upper edge to the top of the sPHENIX mark; from
   the frame's right and upper edge to the legend; the width of a legend symbol; the length of an axis tick; from the
   header's lower edge to the input note. */
const double kHeaderTopPx = 14, kLegendRightPx = 14, kLegendTopPx = 10, kSymbolPx = 46, kTickPx = 14;
const double kNoteDropPx = 0.75 * kStepPx;
/* The tick labels and the titles of an axis set by pixelAxes. The xJ figures label their axes at kLinePx and title
   them at kXjTitlePx instead (xjAxes), as the IAN does: two sizes on purpose. */
const double kLabelPx = 24, kTitlePx = 26;
/* The canvases, in pixels: a single frame, a frame with a ratio pad under it, and a response map. */
const int kCanvasW = 800, kCanvasH = 600, kRatioCanvasH = 760, kMapCanvasH = 700;

/* ---- the state of a run and of the figure being drawn ------------------------------------------------------------
   gFig: the output folder. gWritten, gNotDrawn: the lines of FIGURES.txt. gLegacy: the inputs are reference inputs (the
   entry point says when). gPopulation: the status of the input coverage that the files record (PARTIAL,
   REFERENCE_INPUTS, ...; UNKNOWN when they record none). gResultMetadata, gQAMetadata: the verified metadata of the
   two input files, null without it.
   gRenderInputs, gStylePins: the hashes of the inputs and of the drawing code, for the receipts. gSemanticPackages: the
   packages the figure being drawn has read (qa/PKG, results/PKG). gVisibleText: text that is drawn in pieces, recorded
   whole for the receipt (filled by QAFigures.C). gProjectionSteps: the projections of the figure being drawn (filled
   by QAFigures.C). save() clears these three. gUnique: the counter of unique object names. */
std::string gFig;
/* Comparisons use the same drawing function and style, with explicit saved-input recipes. */
nlohmann::json gCombinatoricRecipe;
nlohmann::json gCombinatoricFrame;
std::vector<std::string> gWritten, gNotDrawn;
bool gLegacy = false;
std::string gPopulation;
std::unique_ptr<TFile> gSystematics;  /* systematics.root of final_analysis.sh systematics, when given */
using PlotJson = nlohmann::json;
PlotJson gResultMetadata, gQAMetadata, gRenderInputs, gStylePins;
PlotJson gSystematicCoverage;
std::vector<std::string> gVisibleText;
std::vector<std::string> gInputNotes;  /* what the receipt says about the state of the figure's inputs */
std::set<std::string> gSemanticPackages;
PlotJson gProjectionSteps = PlotJson::array();
#ifdef PHOTONJET_NATIVE_CAPTURE_TEST
/* Test-only sink executes maintained drawing entrypoints without writing any
   image or receipt. Production builds contain no alternate save route. */
std::function<void(TCanvas&,const std::string&,const std::string&)> gNativeCaptureSink;
#endif
/* Input observations explain changes; they deliberately do not certify the
   final plotted quantities or complete PlotPhysics coverage. Source identity
   lives on the histogram, so Clone preserves it and pointer reuse cannot lie. */
PlotJson gInputObservations = PlotJson::array();
PlotJson histogramSources(const TH1 *histogram) {
  const auto *tag = dynamic_cast<TNamed *>(histogram->GetListOfFunctions()->FindObject("ian_input_sources"));
  return tag ? PlotJson::parse(tag->GetTitle()) : PlotJson::array();
}
void bindHistogramSource(TH1 *histogram, const std::string &file, const std::string &object) {
  if (!gRenderInputs.contains(file)) return;
  const PlotJson sources = PlotJson::array({{{"input_sha256", gRenderInputs.at(file)}, {"object", object}}});
  histogram->GetListOfFunctions()->Add(new TNamed("ian_input_sources", sources.dump().c_str()));
}

/* Per-output scope: an early return or exception cannot carry scientific
   quantities into the next output. Uninstrumented functions emit UNKNOWN. */
struct PlotPhysicsScope;
PlotPhysicsScope *gPlotPhysics = nullptr;
struct PlotPhysicsScope {
  std::string output;
  photonjet_render::PlotPhysics record;
  PlotPhysicsScope(const std::string &file, const std::string &family) : output(file), record(family) {
    if (gPlotPhysics) throw std::runtime_error("nested plot physics scope");
    gPlotPhysics = this;
  }
  ~PlotPhysicsScope() { gPlotPhysics = nullptr; }
  PlotPhysicsScope(const PlotPhysicsScope &) = delete;
  PlotPhysicsScope &operator=(const PlotPhysicsScope &) = delete;
};
PlotJson physicsForOutput(const std::string &output) {
  if (gPlotPhysics && gPlotPhysics->output == output) return gPlotPhysics->record.finish(output);
  return photonjet_render::PlotPhysics().finish(output);
}

/* The last frame drawn by ianFrame and its y range. headroom, openUnderHeader, endAtAxis and whiteCentres act on this
   frame, not on the current pad; ratioFrame, the grids of frames and the response maps do not set it. */
TH1F *gFrame = nullptr;
double gFrameLow = 0, gFrameHigh = 1;
int gUnique = 0;

/* A unique name for a ROOT object. */
std::string fresh() { return Form("ian%d", ++gUnique); }

/* Forgets what the figure being drawn has read, written in pieces and projected. */
void clearFigureState() {
  gVisibleText.clear();
  gInputNotes.clear();
  gSemanticPackages.clear();
  gProjectionSteps = PlotJson::array();
  gInputObservations = PlotJson::array();
}

/* A symbol with its index as a TLatex fragment of its own, for a symbol late in a line. TLatex places an index by the
   measured width of the whole run of text before it, which the written image does not draw exactly that wide, so the
   index of a late symbol stands off from it; set apart, the symbol carries its index tightly and the difference falls
   into the word space before it. */
std::string apart(const std::string &text) { return "#scale[1]{" + text + "}"; }
/* |x| in running text. TLatex sets the opening bar against a Greek letter that follows it, so a thin space stands
   between the two. */
std::string absolute(const std::string &text) { return "|#scale[0.4]{ }" + text + "|"; }
/* "k_Bayes = n". The sign stands a thin space after the index: a word space after an index reads as two. */
std::string bayesIterations(int k) { return Form("%s#scale[0.6]{ }= %d", apart("k_{Bayes}").c_str(), k); }

/* Records a figure that is not drawn, with the reason, for FIGURES.txt. The state of the figure stays as it is. */
void missing(const std::string &name, const std::string &why) {
  gNotDrawn.push_back(name + ": " + why);
  printf("IANFigures: not drawn %s (%s)\n", name.c_str(), why.c_str());
}

/* ---- inputs: metadata, settings, results and histograms ---------------------------------------------------------- */
/* The verified metadata saved in `file`, or null when the file has none: a file without saved metadata is a reference
   input. Records the hash of the file among the inputs and throws when the file or its metadata differ from the
   receipt beside it (<path>.render.json). */
PlotJson loadRenderMetadata(TFile &file, const std::string &path) {
  fa::observeNativeInput(file);
  auto *saved = dynamic_cast<TObjString *>(file.Get("render_metadata"));
  gRenderInputs[path] = photonjet_render::sha256(path, true);
  if (!saved) return PlotJson();
  const PlotJson sealed = PlotJson::parse(saved->GetString().Data());
  const PlotJson receipt = photonjet_render::read(path + ".render.json");
  if (receipt.at("product_sha256") != gRenderInputs.at(path) || receipt.at("metadata_sha256") != sealed.at("sha256"))
    throw std::runtime_error("render product or metadata changed: " + path);
  return photonjet_render::verify(sealed);
}
/* The saved record of a package, from the metadata of the QA file (qa) or of the results file; null when the file has
   no metadata or does not name the package. key: a package as gSemanticPackages names it, qa/PKG or results/PKG. */
const PlotJson *packageRecord(const std::string &pkg, bool qa) {
  const auto &metadata = qa ? gQAMetadata : gResultMetadata;
  if (metadata.is_null() || !metadata.at("packages").contains(pkg)) return nullptr;
  return &metadata.at("packages").at(pkg);
}
const PlotJson *packageRecord(const std::string &key) {
  return packageRecord(key.substr(key.find('/') + 1), key.rfind("qa/", 0) == 0);
}
/* The selection of a package from the saved metadata, null without it. Binds the figure being drawn to the package:
   save() checks the weights of the package and records its selection in the receipt. */
const PlotJson &renderSelection(const std::string &pkg, bool qa = false) {
  static const PlotJson missingSelection;
  gSemanticPackages.insert((qa ? "qa/" : "results/") + pkg);
  const PlotJson *record = packageRecord(pkg, qa);
  return record ? record->at("selection") : missingSelection;
}
/* The collision energy of a package as a header line, from its selection; empty without it. */
std::string energyLabel(const std::string &pkg, bool qa = false) {
  const auto &m = renderSelection(pkg, qa);
  if (m.is_null()) return "";
  return Form(pkg == "PP" ? "#sqrt{s} = %g GeV" : "#sqrt{s_{NN}} = %g GeV",
              m.at("collision").at("sqrt_s_nn_gev").get<double>());
}
/* The collision system of a package and, for a centrality class, its range in percent (written with a hyphen;
   rangeDashes sets the IAN's dash where a figure has it). */
std::string collisionLabel(const std::string &pkg, const PlotJson &range) {
  const std::string system = pkg == "PP" ? "p+p" : "Au+Au";
  return range.is_null() ? system : system + Form(" %g-%g%%", range.at(0).get<double>(), range.at(1).get<double>());
}
/* A saved cut operator in TLatex grammar. */
std::string cutOperator(const PlotJson &value) {
  const auto text = value.get<std::string>();
  if (text == ">=") return "#geq";
  if (text == ">") return ">";
  throw std::runtime_error("unsupported cut operator " + text);
}

/* The settings of a result, saved as the string <key>/settings (name=value;name=value), as a map. */
std::map<std::string, std::string> settingsOf(TFile &f, const std::string &key) {
  std::map<std::string, std::string> out;
  auto *o = dynamic_cast<TObjString *>(f.Get((key + "/settings").c_str()));
  if (!o) return out;
  std::istringstream in(o->GetString().Data());
  for (std::string item; std::getline(in, item, ';');) {
    const auto eq = item.find('=');
    if (eq != std::string::npos) out[item.substr(0, eq)] = item.substr(eq + 1);
  }
  return out;
}

/* The results in file order. IAN names: p+p is "pp"; Au+Au class centA_B is "auauAB" and "auau_AtoB". With metadata,
   the label of a result is compiled from it, and settings that differ from the metadata are an error. */
std::vector<Result> resultsOf(TFile &f) {
  std::vector<Result> out;
  const size_t prefix = std::strlen("cent");  /* a class is saved under the name cent<low>_<high> */
  TIter next(f.GetListOfKeys());
  while (auto *k = static_cast<TKey *>(next())) {
    if (std::string(k->GetClassName()) != "TDirectoryFile") continue;
    auto settings = settingsOf(f, k->GetName());
    if (settings.empty()) continue;
    Result r;
    r.key = k->GetName();
    r.package = settings["package"];
    r.label = settings["label"];
    r.energyLine = energyLabel(r.package);
    const std::string c = settings["centrality"];
    r.centrality = c;
    const auto u = c.find('_');
    if (r.package == "PP") {
      r.fileTag = r.longFileTag = "pp";
    } else {
      const bool named = c.rfind("cent", 0) == 0 && u != std::string::npos;
      r.fileTag = named ? "auau" + c.substr(prefix, u - prefix) + c.substr(u + 1) : r.key;
      r.longFileTag = named ? "auau_" + c.substr(prefix, u - prefix) + "to" + c.substr(u + 1) : r.key;
      if (named) r.centralityLow = std::atof(c.substr(prefix, u - prefix).c_str());
    }
    r.jetThreshold = std::atof(settings["jet_threshold_gev"].c_str());
    std::istringstream range(settings["reported_xj"]);
    char comma;
    range >> r.xjLow >> comma >> r.xjHigh;
    r.combinatoric = settings["combinatoric"];
    std::istringstream sideband(settings["dphi_sideband"]);
    sideband >> r.sidebandLow >> comma >> r.sidebandHigh;
    r.photonIterations = std::atoi(settings["photon_iterations"].c_str());
    if (!gResultMetadata.is_null()) {
      const auto &metadata = gResultMetadata.at("results").at(r.key);
      if (metadata.at("package").get<std::string>() != r.package || metadata.at("centrality").get<std::string>() != c ||
          metadata.at("recoil_min_gev").get<double>() != r.jetThreshold ||
          metadata.at("reported_xj").at(0).get<double>() != r.xjLow ||
          metadata.at("reported_xj").at(1).get<double>() != r.xjHigh)
        throw std::runtime_error("result settings differ from the saved annotation metadata: " + r.key);
      r.label = collisionLabel(r.package, metadata.at("centrality_percent"));
      if (!metadata.at("centrality_percent").is_null()) r.centralityLow = metadata.at("centrality_percent").at(0);
    }
    out.push_back(r);
  }
  /* the rank of a class among the classes of the file, from the most central on: a class is more central when the
     lower edge in its saved name lies under the lower edge of this one */
  for (auto &r : out)
    for (const auto &other : out)
      if (other.centrality.rfind("cent", 0) == 0 &&
          std::atof(other.centrality.c_str() + prefix) < r.centralityLow - 1e-9)
        ++r.centralRank;
  return out;
}

/* A detached copy of the histogram <result>/<name> of the results file, or null. Binds the figure being drawn to the
   package of the result (see renderSelection). */
TH1 *get(TFile &f, const Result &r, const std::string &name) {
  gSemanticPackages.insert("results/" + r.package);
  auto *h = dynamic_cast<TH1 *>(f.Get((r.key + "/" + name).c_str()));
  if (!h) return nullptr;
  h = static_cast<TH1 *>(h->Clone((r.key + name + fresh()).c_str()));
  h->SetDirectory(nullptr);
  bindHistogramSource(h, f.GetName(), r.key + "/" + name);
  return h;
}

/* The bin edges of an axis. */
std::vector<double> edgesOf(const TAxis *a) {
  std::vector<double> e;
  for (int i = 1; i <= a->GetNbins() + 1; ++i) e.push_back(a->GetBinLowEdge(i));
  return e;
}

/* A bin from a to b lies inside the range [low, high] when both its edges do; the tolerance absorbs the rounding of
   saved edges. Every figure applies a reported range with this test. */
const double kEdgeTolerance = 1e-9;
bool binInside(double a, double b, double low, double high) {
  return !(a < low - kEdgeTolerance || b > high + kEdgeTolerance);
}

/* ---- the IAN's physics lines, compiled from the saved metadata and histograms ---------------------------------- */
/* The collision system, the centrality class and the collision energy of a result. */
std::string systemLine(const Result &r) { return r.label + (r.energyLine.empty() ? "" : ", " + r.energyLine); }
/* A header line or legend label with its ranges set as the IAN sets them: a dash of the length of a digit between
   two numbers (0-20% and 17-33 GeV), not a hyphen. */
std::string rangeDashes(std::string text) {
  for (size_t at = text.find('-'); at != std::string::npos; at = text.find('-', at + 1))
    if (at > 0 && at + 1 < text.size() && std::isdigit(static_cast<unsigned char>(text[at - 1])) &&
        std::isdigit(static_cast<unsigned char>(text[at + 1])))
      text.replace(at, 1, "#minus");
  return text;
}
/* The photon bins of a result: the metadata's photon_edges_gev, or for reference inputs without metadata the saved
   photon axis of the result. */
std::vector<double> photonEdges(TFile &f, const Result &r) {
  if (!gResultMetadata.is_null() && gResultMetadata.at("results").contains(r.key))
    return gResultMetadata.at("results").at(r.key).at("photon_edges_gev").get<std::vector<double>>();
  TH1 *h = get(f, r, "photons_corrected");
  return h ? edgesOf(h->GetXaxis()) : std::vector<double>{};
}
/* A photon pT interval as the IAN writes it, and the photon window of a result (all of its photon bins). */
std::string photonInterval(double low, double high) { return Form("%g #leq p_{T}^{#gamma} < %g GeV", low, high); }
std::string photonLine(TFile &f, const Result &r) {
  const auto e = photonEdges(f, r);
  return e.size() < 2 ? "" : photonInterval(e.front(), e.back());
}
/* The photon identification of the measurement, as the IAN's figures name it. */
const char *kPhotonId = "Tight, isolated #gamma";
/* The two halves of the recoil-jet definition, from the jet record of a selection: the anti-kT radius, and the pT
   threshold with its saved operator. */
std::string jetRadiusText(const PlotJson &jet) { return Form("anti-k_{T} R = %g", jet.at("radius").get<double>()); }
std::string jetThresholdText(const PlotJson &jet, double threshold) {
  return Form("p_{T}^{jet} %s %g GeV", cutOperator(jet.at("pt_operator")).c_str(), threshold);
}
/* The recoil-jet definition: the anti-kT radius from the metadata (absent for reference inputs) and the threshold. */
std::string jetLine(const Result &r) {
  const auto &selection = renderSelection(r.package);
  if (selection.is_null()) return Form("p_{T}^{jet} #geq %g GeV", r.jetThreshold);
  const auto &jet = selection.at("jet");
  return jetRadiusText(jet) + ", " + jetThresholdText(jet, r.jetThreshold);
}
/* The recoil-angle requirement of a result as the IAN writes it (>= 7pi/8), from the metadata; empty without it. */
std::string recoilCut(const Result &r) {
  const auto &selection = renderSelection(r.package);
  if (selection.is_null()) return "";
  const auto &jet = selection.at("jet");
  const auto &fraction = jet.at("delta_phi_pi_fraction");
  const std::string angle = fraction.is_null() ? Form("%.3g", jet.at("delta_phi_min_rad").get<double>())
                                               : Form("%d#pi/%d", fraction.at(0).get<int>(), fraction.at(1).get<int>());
  return cutOperator(jet.at("delta_phi_operator")) + " " + angle;
}
/* The photon and jet acceptance and the recoil angle of a result in one header line, from the metadata; empty without
   it. */
std::string recoilLine(const Result &r) {
  const auto &selection = renderSelection(r.package);
  if (selection.is_null()) return "";
  const double photon = selection.at("photon").at("abs_eta_max").get<double>();
  const double jet = selection.at("jet").at("abs_eta_max").get<double>();
  const std::string both = absolute("#eta^{#gamma,jet}"), one = absolute("#eta^{#gamma}"), other = absolute("#eta^{jet}");
  const std::string eta = photon == jet ? Form("%s < %g", both.c_str(), photon)
                                        : Form("%s < %g, %s < %g", one.c_str(), photon, other.c_str(), jet);
  return eta + ", #Delta#phi_{#gamma,jet} " + recoilCut(r);
}
/* The header lines of an xJ figure: the system, the photon window and the recoil jets. */
std::vector<std::string> xjLines(TFile &f, const Result &r) { return {systemLine(r), photonLine(f, r), jetLine(r)}; }
/* The same lines with the centrality range written with a dash, as in the IAN (0–20%). */
std::vector<std::string> xjHeader(TFile &f, const Result &r) {
  auto lines = xjLines(f, r);
  lines.front() = rangeDashes(lines.front());
  return lines;
}

/* ---- text, header and legend ------------------------------------------------------------------------------------- */
/* The width and the height of the current pad in canvas pixels; the layout helpers call them w and h. */
double padWidthPx() { return gPad->GetWw() * gPad->GetAbsWNDC(); }
double padHeightPx() { return gPad->GetWh() * gPad->GetAbsHNDC(); }
/* Axis ticks `px` canvas pixels long on both axes of a frame in the current pad (minor ticks half as long), whatever
   the size and the margins of the pad. ROOT takes the length of an x tick as a share of the pad's height times the
   share of the pad's width that the frame spans, and the length of a y tick the other way round. */
void tickLengths(TH1 *frame, double px) {
  const double across = 1 - gPad->GetLeftMargin() - gPad->GetRightMargin();
  const double up = 1 - gPad->GetTopMargin() - gPad->GetBottomMargin();
  frame->GetXaxis()->SetTickLength(px / (padHeightPx() * across));
  frame->GetYaxis()->SetTickLength(px / (padWidthPx() * up));
}
/* The width of text in the current pad, in NDC. */
double textWidth(const std::string &text, double px, int font = kFont) {
  TLatex probe(0.5, 0.5, text.c_str());
  probe.SetNDC();
  probe.SetTextFont(font);
  probe.SetTextSize(px);
  UInt_t w = 0, h = 0;
  probe.GetBoundingBox(w, h);
  return w / padWidthPx();
}
/* `text` with every `from` written as `to`. */
std::string rewritten(std::string text, const std::string &from, const std::string &to) {
  for (size_t at = text.find(from); at != std::string::npos; at = text.find(from, at + to.size()))
    text.replace(at, from.size(), to);
  return text;
}
/* Producer text in TLatex grammar: relations and the vertex symbol as the IAN writes them. */
std::string latexCuts(std::string s) {
  const std::pair<std::string, std::string> words[] = {{">=", "#geq"}, {"<=", "#leq"}, {"vtx", "#left|v_{z}#right|"}};
  for (const auto &[from, to] : words) s = rewritten(s, from, to);
  return s;
}
/* The lines that are not empty: a header leaves out a line that the inputs cannot fill. */
std::vector<std::string> nonEmpty(const std::vector<std::string> &lines) {
  std::vector<std::string> out;
  for (const auto &l : lines) if (!l.empty()) out.push_back(l);
  return out;
}
/* The state of a figure's inputs (reference inputs, a partial or an unverified population) belongs in the receipt
   of the figure and never on its canvas: records `note` for the receipt and returns the empty text the canvas gets. */
std::string forReceipt(const std::string &note) {
  if (!note.empty() && std::find(gInputNotes.begin(), gInputNotes.end(), note) == gInputNotes.end())
    gInputNotes.push_back(note);
  return "";
}
/* Records for the receipt that the inputs are not a complete physics population; the canvas carries no such line. */
std::string inputNote() {
  return forReceipt(gLegacy ? "Reference inputs (development)" : gPopulation == "PARTIAL" ? "Partial input population" :
                    gPopulation == "UNKNOWN" ? "Input coverage unverified" : "");
}
/* The header: the sPHENIX mark at (x, y), the lines below it and (note) the input note. Returns the lowest baseline
   (NDC); draw = false only measures. */
double header(double x, double y, const std::vector<std::string> &lines, bool draw = true, bool note = true) {
  const double h = padHeightPx();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(kFont);
  t.SetTextSize(kMarkPx);
  if (draw) t.DrawLatex(x, y, sPhenixInternalLabel());
  t.SetTextSize(kLinePx);
  double last = y - 0.15 * kStepPx / h;  /* the first line stands 0.15 of a pitch lower, under the larger mark */
  for (const auto &line : nonEmpty(lines)) {
    last -= kStepPx / h;
    if (draw) t.DrawLatex(x, last, line.c_str());
  }
  if (const std::string text = inputNote(); note && !text.empty()) {
    last -= kStepPx / h;
    t.SetTextColor(gGray);
    if (draw) t.DrawLatex(x, last, text.c_str());
  }
  return last;
}
/* The width of the header in the current pad (NDC): the sPHENIX mark or its widest line. */
double headerWidth(const std::vector<std::string> &lines) {
  double w = textWidth(sPhenixInternalLabel(), kMarkPx);
  for (const auto &l : nonEmpty(lines)) w = std::max(w, textWidth(l, kLinePx));
  return w;
}
/* An empty legend without a box, in the text size of the header lines. symbolShare: the share of its width that the
   symbols take. */
TLegend *ianLegend(double x1, double y1, double x2, double y2, double symbolShare) {
  auto *l = new TLegend(x1, y1, x2, y2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(kFont);
  l->SetTextSize(kLinePx);
  l->SetMargin(symbolShare);
  return l;
}
/* The same legend in a white box with a border, for a legend that stands over the frames of a grid. */
TLegend *boxedLegend(double x1, double y1, double x2, double y2, double symbolShare) {
  TLegend *l = ianLegend(x1, y1, x2, y2, symbolShare);
  l->SetBorderSize(1);
  l->SetFillStyle(1001);
  l->SetFillColor(kWhite);
  return l;
}

/* ---- the frame: header top left, legend top right, a y range opened for both ----------------------------------- */
struct Extent { double x0, x1, y; };  /* a drawn element: its horizontal extent and its top, in user coordinates */
/* Appends the extents of the points of a graph, with their bars, or of the bins of a histogram (errors: with their
   error bars). */
void extents(std::vector<Extent> &out, const TGraph *g) {
  const auto *e = dynamic_cast<const TGraphErrors *>(g);
  const auto *a = dynamic_cast<const TGraphAsymmErrors *>(g);
  for (int i = 0; i < g->GetN(); ++i) {
    const double xl = e ? e->GetEX()[i] : a ? a->GetEXlow()[i] : 0, xh = e ? e->GetEX()[i] : a ? a->GetEXhigh()[i] : 0;
    const double yh = e ? e->GetEY()[i] : a ? a->GetEYhigh()[i] : 0;
    out.push_back({g->GetX()[i] - xl, g->GetX()[i] + xh, g->GetY()[i] + std::max(0.0, yh)});
  }
}
void extents(std::vector<Extent> &out, const TH1 *h, bool errors = false) {
  for (int i = 1; i <= h->GetNbinsX(); ++i)
    out.push_back({h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(i),
                   h->GetBinContent(i) + (errors ? h->GetBinError(i) : 0)});
}

/* A single-pad frame in the IAN's layout and its empty legend: the header at the top left, the legend at the top right,
   or under the header when the two do not fit side by side, and a y range opened until no drawn element reaches either
   block. x0, x1: the x range. low: the frame's lower edge (on a log scale, a third of the smallest positive top when
   low <= 0). drawn: the extents of everything the caller will draw. lines: the header. labels: the legend texts, which
   size the legend. populationNote: the input note of this figure when it is not the general one (inputNote). Sets
   gFrame, gFrameLow and gFrameHigh. The caller adds the legend entries in `labels` order and draws the legend. */
TLegend *ianFrame(double x0, double x1, double low, const std::vector<Extent> &drawn, bool logy, const char *xtitle,
                  const char *ytitle, const std::vector<std::string> &lines, const std::vector<std::string> &labels,
                  const std::string *populationNote = nullptr) {
  gPad->SetLogy(logy);
  const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin(), T = gPad->GetTopMargin();
  const double B = gPad->GetBottomMargin();
  const double plot = 1 - T - B, h = padHeightPx(), w = padWidthPx();
  const double hx = L + kInsetPx / w, top = 1 - T - (kHeaderTopPx + kMarkPx) / h;
  const double hw = headerWidth(lines), hBottom = header(L, top, lines, false, false) - 10 / h;
  const std::string note = populationNote ? *populationNote : inputNote();
  double lw = 0;
  for (const auto &l : labels) lw = std::max(lw, textWidth(l, kLinePx));
  const double symbol = kSymbolPx / w, columnWidth = symbol + lw + 10 / w;
  const double noteWidth = note.empty() ? 0.0 : textWidth(note, kLinePx);

  double highest = -std::numeric_limits<double>::infinity(), positive = std::numeric_limits<double>::infinity();
  for (const auto &e : drawn)
    if (std::isfinite(e.y)) { highest = std::max(highest, e.y); if (e.y > 0) positive = std::min(positive, e.y); }
  if (!std::isfinite(highest)) highest = 1;
  if (logy && !(low > 0)) low = std::isfinite(positive) ? positive / 3 : 1e-3;
  const auto ndc = [&](double x) { return L + (x - x0) / (x1 - x0) * (1 - L - R); };
  const auto frameTop = [&](double y, double fraction) {  /* the top that puts y at `fraction` of the frame height */
    fraction = std::max(fraction, 0.1);
    return logy ? std::exp(std::log(low) + (std::log(std::max(y, low)) - std::log(low)) / fraction)
                : low + (y - low) / fraction;
  };

  /* The legend at the top right; a large overlay of samples takes two columns, so the plot keeps usable space. */
  int columns = labels.size() > 7 && 2 * columnWidth < 1 - L - R ? 2 : 1;
  double width = columns * columnWidth, height = ((labels.size() + columns - 1) / columns) * kStepPx / h;
  double lx2 = 1 - R - kLegendRightPx / w, lx1 = lx2 - width, ly2 = 1 - T - kLegendTopPx / h, ly1 = ly2 - height;
  if (lx1 < hx + hw + 24 / w) {  /* no room side by side: the legend goes under the header */
    /* stacked under the header it shares the data's side: two columns keep the block, and so the y range, short */
    if (columns == 1 && labels.size() >= 4 && hx + 2 * columnWidth < 1 - R - kLegendRightPx / w)
      columns = 2, width = 2 * columnWidth, height = ((labels.size() + 1) / 2) * kStepPx / h;
    lx1 = hx, lx2 = hx + width, ly2 = hBottom - 4 / h, ly1 = ly2 - height;
  }
  /* the input note: under the header, and under the legend too when it would run into it */
  double noteY = hBottom - kNoteDropPx / h;
  if (!note.empty() && hx + noteWidth > lx1 - 6 / w && noteY > ly1 - kNoteDropPx / h) noteY = ly1 - kNoteDropPx / h;
  if (!(std::min(ly1, note.empty() ? 1.0 : noteY) > B + 0.05))
    throw std::runtime_error("figure annotations do not fit the canvas");
  /* the frame top that keeps every drawn element clear of the header, the legend and the note */
  double high = frameTop(highest, 0.95);
  std::vector<std::array<double, 3>> blocks = {{hx - 6 / w, hx + hw + 6 / w, hBottom - 6 / h},
                                               {lx1 - 6 / w, lx2 + 6 / w, ly1 - 8 / h}};
  if (!note.empty()) blocks.push_back({hx - 6 / w, hx + noteWidth + 6 / w, noteY - 8 / h});
  for (const auto &e : drawn) {
    if (!std::isfinite(e.y)) continue;
    const double xa = ndc(e.x0), xb = ndc(e.x1);
    for (const auto &[ka, kb, bottom] : blocks)
      if (xb > ka && xa < kb) high = std::max(high, frameTop(e.y, (bottom - B) / plot));
  }
  if (!(high > low)) high = low + 1;
  gFrame = photonjet_science::drawFrame(*gPad,x0, low, x1, high);
  gFrameLow = low, gFrameHigh = high;
  gFrame->GetXaxis()->SetTitle(xtitle);
  gFrame->GetYaxis()->SetTitle(ytitle);
  /* counts from a thousand take a x10^n exponent instead of crowding the title; small fractions are written out */
  if (high >= 1000) gFrame->GetYaxis()->SetMaxDigits(3);
  else if (!logy) gFrame->GetYaxis()->SetNoExponent();
  /* ticks of one length on every canvas size, so they end short of the header */
  tickLengths(gFrame, kTickPx);
  header(hx, top, lines, true, false);
  if (!note.empty()) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(kFont);
    t.SetTextSize(kLinePx);
    t.SetTextColor(gGray);
    t.DrawLatex(hx, noteY, note.c_str());
  }
  auto *legend = ianLegend(lx1, ly1, lx2, ly2, symbol / columnWidth);
  legend->SetNColumns(columns);
  return legend;
}

/* ---- ratio pads and axis text ------------------------------------------------------------------------------------ */
/* A figure with a ratio panel below: the upper (70%) and lower (30%) pads. left, right: the side margins of both pads;
   0 keeps the style's. */
std::pair<TPad *, TPad *> ratioPads(TCanvas &cv, double left = 0, double right = 0) {
  cv.cd();
  auto *top = new TPad(fresh().c_str(), "", 0, 0.3, 1, 1), *bottom = new TPad(fresh().c_str(), "", 0, 0, 1, 0.3);
  top->SetBottomMargin(0.025);
  bottom->SetTopMargin(0.035);
  bottom->SetBottomMargin(0.44);
  photonjet_science::draw(top);
  photonjet_science::draw(bottom);
  for (TPad *pad : {top, bottom}) {
    if (left > 0) pad->SetLeftMargin(left);
    if (right > 0) pad->SetRightMargin(right);
  }
  return {top, bottom};
}
/* A see-through pad over the whole canvas with the given margins (fractions of the canvas), drawn. The grids of frames
   give every frame a pad of this kind, so all frames are measured in the same pixels and no label of a frame is cut at
   its neighbour. */
TPad *seeThroughPad(TCanvas &cv, double left, double right, double bottom, double top) {
  cv.cd();
  auto *pad = new TPad(fresh().c_str(), "", 0, 0, 1, 1);
  pad->SetFillStyle(4000);
  pad->SetMargin(left, right, bottom, top);
  photonjet_science::draw(pad);
  return pad;
}
/* Axis text in pixels, so the upper and the lower pad of a ratio figure read alike. */
void pixelAxes(TH1 *frame, double label = kLabelPx, double title = kTitlePx) {
  for (TAxis *a : {frame->GetXaxis(), frame->GetYaxis()})
    a->SetLabelFont(kFont), a->SetLabelSize(label), a->SetTitleFont(kFont), a->SetTitleSize(title);
  frame->GetYaxis()->SetTitleOffset(2.0);
  frame->GetXaxis()->SetTitleOffset(1.6);
}
/* No x labels and no x title: for the pad above a ratio pad, whose lower pad carries them. */
void hideXText(TH1 *frame) {
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
}
/* The axis text of the upper pad of a ratio figure. */
void upperAxes(TH1 *frame) {
  pixelAxes(frame);
  hideXText(frame);
}
/* The lower frame of a ratio figure over [y0, y1]. It does not become gFrame. */
TH1F *ratioFrame(double x0, double x1, double y0, double y1, const char *xtitle, const char *ytitle) {
  TH1F *f = photonjet_science::drawFrame(*gPad,x0, y0, x1, y1);
  f->GetXaxis()->SetTitle(xtitle);
  f->GetYaxis()->SetTitle(ytitle);
  pixelAxes(f);
  f->GetYaxis()->SetNdivisions(503);
  return f;
}
/* A dashed straight line between two points of the current pad, given in its axes or, with ndc, as fractions of the
   pad, in the one dash pattern the entry point sets for the line style. Returns the line drawn, for a legend. */
TLine *dashedLine(double x0, double y0, double x1, double y1, int color = kBlack, int width = 1, bool ndc = false) {
  TLine l;
  l.SetLineStyle(kDashed);
  l.SetLineColor(color);
  l.SetLineWidth(width);
  return ndc ? l.DrawLineNDC(x0, y0, x1, y1) : l.DrawLine(x0, y0, x1, y1);
}
/* A horizontal reference across the visible frame, in data coordinates. Read the painted pad after an axis
   zoom, so a stored histogram range, changed margins or logarithmic x axis cannot extend it into the labels. */
TLine *frameReferenceLine(double value, int color = kBlack, int width = 1, int style = kDashed) {
  gPad->Modified();
  gPad->Update();
  const double left = gPad->PadtoX(gPad->GetUxmin()), right = gPad->PadtoX(gPad->GetUxmax());
  if (!std::isfinite(left) || !std::isfinite(right) || !(left < right))
    throw std::runtime_error("reference line requires a nonempty painted frame");
  TLine line;
  line.SetLineColor(color);
  line.SetLineWidth(width);
  line.SetLineStyle(style);
  return line.DrawLine(left, value, right, value);
}
/* A photon axis with few bins is labelled at its bin edges, as in the IAN (17, 23, 33): a tick at every edge on the
   lower and the upper side of the frame and, with `labels`, the value of the edge under the axis. */
void edgeTicks(TH1 *frame, const std::vector<double> &edges, bool labels) {
  if (edges.size() < 2 || edges.size() > 7) return;
  TAxis *a = frame->GetXaxis();
  a->SetTickLength(0);  /* ROOT's ticks and labels give way to the edges; the axis title stays */
  a->SetLabelSize(0);
  const double y0 = frame->GetMinimum(), y1 = frame->GetMaximum(), h = padHeightPx();
  const bool logy = gPad->GetLogy();
  const auto at = [&](double fraction) {
    return logy ? std::exp(std::log(y0) + fraction * (std::log(y1) - std::log(y0))) : y0 + fraction * (y1 - y0);
  };
  const double plotPx = (1 - gPad->GetTopMargin() - gPad->GetBottomMargin()) * h, tick = kTickPx / plotPx;
  TLine l;
  l.SetLineWidth(1);
  TLatex t;
  t.SetTextFont(kFont);
  t.SetTextSize(kLabelPx);
  t.SetTextAlign(23);
  for (double e : edges) {
    if (e < a->GetXmin() - 1e-9 || e > a->GetXmax() + 1e-9) continue;
    l.DrawLine(e, at(0), e, at(tick));
    l.DrawLine(e, at(1 - tick), e, at(1));
    if (labels) {
      t.SetNDC(false);
      t.DrawLatex(e,
                  logy ? std::exp(std::log(y0) - 10 / plotPx * (std::log(y1) - std::log(y0)))
                       : y0 - 10 / plotPx * (y1 - y0),
                  Form("%g", e));
    }
  }
}

/* ---- saving: candidate legends, the text of a canvas, the image and its receipt ---------------------------------- */
/* Adds `hist` to `legend` under the label of its candidate population (a truth class of photon candidates) and tags the
   histogram with the population and its sample family, so that canvasText can check the label when the figure is
   saved. family: the sample the candidates come from; it is provenance and never stands for a class. style: the legend
   option. sourceSlices: the generator slices of the sample, in GeV. */
void candidateEntry(TLegend *legend, TH1 *hist, photonjet_render::CandidatePopulation population,
                    const std::string &family, const char *style = "l", const std::vector<int> &sourceSlices = {}) {
  auto identity = photonjet_render::candidatePopulation(population);
  if (family.empty()) throw std::runtime_error("candidate source family is required");
  if (population != photonjet_render::CandidatePopulation::TruthTaggedBackground && family != "inclusive+jet")
    throw std::runtime_error("Inclusive+jet requires the inclusive+jet sample family");
  identity["sample_family"] = family;
  if (!sourceSlices.empty()) identity["source_slices_gev"] = sourceSlices;
  auto *tag = dynamic_cast<TNamed *>(hist->GetListOfFunctions()->FindObject("candidate_population"));
  if (tag && PlotJson::parse(tag->GetTitle()) != identity)
    throw std::runtime_error("histogram already has a different candidate population");
  if (!tag) hist->GetListOfFunctions()->Add(new TNamed("candidate_population", identity.dump().c_str()));
  legend->AddEntry(hist, photonjet_render::candidateLegendLabel(identity).c_str(), style);
}
/* Every text of a canvas, pad by pad, appended to `text`: the TLatex texts, the legend labels and the axis titles. A
   legend entry bound to a candidate population (candidateEntry) is checked against it and the population appended to
   `populations`; an unbound entry that carries one of the candidate labels is an error. */
void canvasText(TPad &pad, std::vector<std::string> &text, PlotJson *populations = nullptr) {
  TIter next(pad.GetListOfPrimitives());
  while (auto *object = next()) {
    if (auto *child = dynamic_cast<TPad *>(object)) canvasText(*child, text, populations);
    else if (auto *latex = dynamic_cast<TLatex *>(object)) text.push_back(latex->GetTitle());
    else if (auto *legend = dynamic_cast<TLegend *>(object)) {
      TIter item(legend->GetListOfPrimitives());
      while (auto *entry = dynamic_cast<TLegendEntry *>(item())) {
        text.push_back(entry->GetLabel());
        bool bound = false;
        if (auto *hist = dynamic_cast<TH1 *>(entry->GetObject()))
          if (auto *tag = dynamic_cast<TNamed *>(hist->GetListOfFunctions()->FindObject("candidate_population"))) {
            const auto identity = PlotJson::parse(tag->GetTitle());
            photonjet_render::requireCandidateLabel(identity, entry->GetLabel());
            if (populations) populations->push_back(identity);
            bound = true;
          }
        const std::string label = entry->GetLabel();
        if (!bound && (label == "Background" || label == "Inclusive+jet" || label == "Other inclusive+jet candidates"))
          throw std::runtime_error("candidate legend requires an explicit population binding: " + label);
      }
    } else if (auto *hist = dynamic_cast<TH1 *>(object)) {
      for (auto *axis : {hist->GetXaxis(), hist->GetYaxis(), hist->GetZaxis()})
        if (*axis->GetTitle()) text.push_back(axis->GetTitle());
    }
  }
}

/* Fits the labels of every legend of a pad, nested pads included, into the legend's box as ROOT sets them. The
   configured font and size are kept as the largest size: a label that is too long or too tall makes the text of its
   legend smaller. The box of a label may be kLegendRowFill times as tall as its row. */
const double kLegendRowFill = 1.25;
void fitLegends(TPad &pad) {
  pad.cd();
  TIter next(pad.GetListOfPrimitives());
  while (auto *object = next()) {
    if (auto *child = dynamic_cast<TPad *>(object)) {
      fitLegends(*child);
      pad.cd();
    } else if (auto *legend = dynamic_cast<TLegend *>(object)) {
      const double size = legend->GetTextSize();
      if (!(size > 0)) continue;  /* ROOT's automatic sizing already fits its box. */
      const int columns = std::max(1, legend->GetNColumns());
      const int rows = (legend->GetListOfPrimitives()->GetSize() + columns - 1) / columns;
      if (!rows) continue;
      const double width = pad.GetWw() * pad.GetAbsWNDC() * (legend->GetX2NDC() - legend->GetX1NDC()) *
                           (1 - legend->GetMargin()) / columns - 6;
      const double height = pad.GetWh() * pad.GetAbsHNDC() * (legend->GetY2NDC() - legend->GetY1NDC()) / rows *
                            kLegendRowFill;
      if (!(width > 0 && height > 0)) throw std::runtime_error("legend has no usable text area");
      double scale = 1;
      TIter item(legend->GetListOfPrimitives());
      while (auto *entry = dynamic_cast<TLegendEntry *>(item())) {
        TLatex probe(0, 0, entry->GetLabel());
        probe.SetNDC();
        probe.SetTextFont(entry->GetTextFont() > 0 ? entry->GetTextFont() : legend->GetTextFont());
        probe.SetTextSize(entry->GetTextSize() > 0 ? entry->GetTextSize() : size);
        UInt_t w = 0, h = 0;
        probe.GetBoundingBox(w, h);
        if (w) scale = std::min(scale, width / w);
        if (h) scale = std::min(scale, height / h);
      }
      if (scale < 1) {
        legend->SetTextSize(size * scale);
        TIter entries(legend->GetListOfPrimitives());
        while (auto *entry = dynamic_cast<TLegendEntry *>(entries()))
          if (entry->GetTextSize() > 0) entry->SetTextSize(entry->GetTextSize() * scale);
      }
    }
  }
  pad.Modified();
}

/* A PNG file is written at kImageScale times the canvas size in pixels, so a figure is sharp at the note's print size.
   ROOT scales the text and the markers of an image by that factor; line widths and dash lengths are in pixels, so
   every line of the canvas is widened and every dash pattern lengthened by the same factor while the file is written,
   and every text is set to the exact size of its font step (exactTextSizes). widenLines returns the lines with their
   widths, to set back. */
const int kImageScale = 3;
void widenLines(TPad &pad, std::vector<std::pair<TAttLine *, Width_t>> &widened, std::set<TObject *> &seen) {
  const auto widen = [&](TObject *object) {
    auto *line = dynamic_cast<TAttLine *>(object);
    if (!line || !seen.insert(object).second) return;
    widened.push_back({line, line->GetLineWidth()});
    line->SetLineWidth(line->GetLineWidth() * kImageScale);
  };
  TIter next(pad.GetListOfPrimitives());
  while (auto *object = next()) {
    if (auto *child = dynamic_cast<TPad *>(object)) { widenLines(*child, widened, seen); continue; }
    widen(object);
    TList *attached = nullptr;  /* fitted curves of a histogram or a graph; the entries of a legend */
    if (auto *hist = dynamic_cast<TH1 *>(object)) attached = hist->GetListOfFunctions();
    else if (auto *graph = dynamic_cast<TGraph *>(object)) attached = graph->GetListOfFunctions();
    else if (auto *legend = dynamic_cast<TLegend *>(object)) attached = legend->GetListOfPrimitives();
    if (!attached) continue;
    TIter item(attached);
    while (auto *part = item()) {
      widen(part);
      if (auto *entry = dynamic_cast<TLegendEntry *>(part)) widen(entry->GetObject());
    }
  }
}

/* The text size nearest to `px` pixels that the written image shows exactly as wide as TLatex measures it. A font is
   set in whole units, of which a pixel of text size holds 0.934: on the canvas, where TLatex measures a text and
   places its parts, at the nearest whole number to the text size; in the written image, kImageScale times as large,
   at the nearest one to 0.985 of the enlarged text size. The size returned sets the image at kImageScale times the
   units of the canvas. At any other size the image shows a text up to 3.5% narrower than measured, which moves an
   index off its symbol and the end of a centred or right-aligned text off its place. From 47 px on no size is exact;
   the one returned is 1% narrow at most. */
double exactTextPx(double px) {
  const double unitsPerPx = 0.93376068, imageShare = 0.985;
  const double units = std::round(px * unitsPerPx), half = .5 / kImageScale;
  const double low = (units - half) / imageShare, high = std::min((units + half) / imageShare, units + .5);
  return (low < high ? (low + high) / 2 : units + .45) / unitsPerPx;
}

/* Sets every text of a pad that is sized in pixels (header, legend and note lines, axis titles) to the
   exact size of its font step while the image is written, so that the image shows each text as wide as it was
   measured and placed on the canvas. The canvas keeps its font steps, so nothing moves. Adds to `restore` what sets
   each size back. */
void exactTextSizes(TPad &pad, std::vector<std::function<void()>> &restore, std::set<TObject *> &seen) {
  const auto inPixels = [](int font) { return font % 10 == 3; };
  const auto text = [&](TObject *object) {
    auto *t = dynamic_cast<TAttText *>(object);
    if (!t || !inPixels(t->GetTextFont()) || !(t->GetTextSize() > 0) || !seen.insert(object).second) return;
    const float size = t->GetTextSize();
    t->SetTextSize(exactTextPx(size));
    restore.push_back([t, size] { t->SetTextSize(size); });
  };
  /* an axis title; the tick labels are single numbers, which keep their size */
  const auto axis = [&](TAxis *a) {
    const float title = a->GetTitleSize();
    if (!inPixels(a->GetTitleFont()) || !(title > 0) || !seen.insert(a).second) return;
    a->SetTitleSize(exactTextPx(title));
    restore.push_back([a, title] { a->SetTitleSize(title); });
  };
  TIter next(pad.GetListOfPrimitives());
  while (auto *object = next()) {
    if (auto *child = dynamic_cast<TPad *>(object)) { exactTextSizes(*child, restore, seen); continue; }
    text(object);
    if (auto *hist = dynamic_cast<TH1 *>(object))
      for (TAxis *a : {hist->GetXaxis(), hist->GetYaxis(), hist->GetZaxis()}) axis(a);
    if (auto *legend = dynamic_cast<TLegend *>(object)) {
      TIter entry(legend->GetListOfPrimitives());
      while (auto *line = entry()) text(line);
    }
  }
}

/* The check of the weights behind a figure, made before its image or its folder exists. `weights` given: the record of
   released cut parameters, of a model-training diagnostic or of a physics figure. None given: the saved weight record
   of every package the figure has read (gSemanticPackages); these are returned per package, for the receipt. */
PlotJson checkedWeights(const PlotJson &weights, const std::string &relative) {
  PlotJson physicsWeights = PlotJson::object();
  if (!weights.is_null()) {
    const std::string role = weights.value("role", "");
    if (role == "released_cut_parameters") {
      const auto hash = weights.value("working_points_sha256", "");
      if (weights.value("parameter_kind", "") != "released_working_points" || hash.size() != 64 ||
          hash.find_first_not_of("0123456789abcdef") != std::string::npos)
        throw std::runtime_error("released parameter authority is missing: " + relative);
    } else if (role == "model_training_diagnostic") {
      photonjet_render::requireModelTrainingDiagnostic(weights, relative);
    } else {
      photonjet_render::requirePhysicsPlotWeights(weights, relative);
    }
    return physicsWeights;
  }
  if (gSemanticPackages.empty()) photonjet_render::requirePhysicsPlotWeights(PlotJson(), relative);
  for (const auto &key : gSemanticPackages) {
    const PlotJson *record = packageRecord(key);
    const PlotJson proof = record ? record->value("weight_provenance", PlotJson()) : PlotJson();
    photonjet_render::requirePhysicsPlotWeights(proof, relative);
    physicsWeights[key] = proof;
  }
  return physicsWeights;
}

/* Writes the canvas to `path`. For a PNG the lines are widened, the dash patterns of the line styles 2 to 10 lengthened
   and the pixel-sized texts set to their exact sizes while the file is written (see kImageScale), and all of it is set
   back afterwards. */
void writeImage(TCanvas &c, const std::string &path) {
  std::vector<std::pair<TAttLine *, Width_t>> widened;
  std::vector<std::string> patterns;  /* the dash patterns of line styles 2 to 10, to set back */
  std::vector<std::function<void()>> textSizes;  /* what sets the text sizes back */
  if (path.size() > 4 && path.compare(path.size() - 4, 4, ".png") == 0) {
    std::set<TObject *> seen, sized;
    widenLines(c, widened, seen);
    exactTextSizes(c, textSizes, sized);
    gStyle->SetLineWidth(kImageScale);  /* the axes and their ticks are drawn with the style's width */
    for (int lineStyle = 2; lineStyle <= 10; ++lineStyle) {
      patterns.push_back(gStyle->GetLineStyleString(lineStyle));
      std::istringstream lengths(patterns.back());
      std::string scaled;
      for (int length; lengths >> length;) scaled += std::to_string(length * kImageScale) + " ";
      gStyle->SetLineStyleString(lineStyle, scaled.c_str());
    }
  }
  c.SaveAs(path.c_str());
  gStyle->SetLineWidth(1);
  for (size_t i = 0; i < patterns.size(); ++i) gStyle->SetLineStyleString(i + 2, patterns[i].c_str());
  for (const auto &[line, width] : widened) line->SetLineWidth(width);
  for (const auto &setBack : textSizes) setBack();
}

/* Every figure ends here: the canvas is written to gFig/relative with its receipt (the same name with .json) and listed
   for FIGURES.txt. Before the image is written the weights are checked (checkedWeights), the legends are fitted to
   their boxes, and the text of the canvas is collected, with every candidate legend entry checked against its
   population (canvasText). The receipt records the hash of the image, the inputs, the drawing code, the visible text,
   the selections of the packages the figure has read and the calling function. function, source: filled in at the call
   by their defaults; a wrapper around save() passes them on, or the receipt would name the wrapper. The state of the
   figure being drawn is cleared here. missing() does not clear it: a figure that ends there after it has read a
   package leaves its state to the next figure. */
/* What keeps a drawing inside its frames (patch 48, 2026-10-15, after the xJ axes moved to the final result's left
   edge and some pieces stayed at x = 0): the stacked full-width pads of one figure (the upper and lower pad of a ratio
   figure) span one x range, and no line or box drawn in a pad's axes reaches past the left or right edge of that pad's
   frame. Returns what is wrong, or "" for a drawing that keeps to its frames. */
std::string frameDefect(TCanvas &canvas) {
  const auto framed = [](TPad *pad) {
    for (TObject *o : *pad->GetListOfPrimitives())
      if (o->InheritsFrom(TH1::Class())) return true;
    return false;
  };
  std::vector<TPad *> pads = {&canvas}, stacked;
  for (TObject *o : *canvas.GetListOfPrimitives())
    if (auto *pad = dynamic_cast<TPad *>(o)) pads.push_back(pad);
  for (TPad *pad : pads) {
    if (!framed(pad)) continue;
    pad->Update();
    const double x0 = pad->GetUxmin(), x1 = pad->GetUxmax(), tolerance = 1e-6 * std::fabs(x1 - x0);
    for (TObject *o : *pad->GetListOfPrimitives()) {
      if (o->InheritsFrom(TPave::Class()) || o->InheritsFrom(TFrame::Class())) continue;  /* legends, text, the frame */
      double a, b;
      if (auto *line = dynamic_cast<TLine *>(o)) {
        if (line->TestBit(TLine::kLineNDC)) continue;  /* drawn in fractions of the pad */
        a = line->GetX1(), b = line->GetX2();
      } else if (auto *box = dynamic_cast<TBox *>(o)) {
        a = box->GetX1(), b = box->GetX2();
      } else {
        continue;
      }
      if (pad->GetLogx()) {
        if (!(a > 0 && b > 0)) continue;
        a = std::log10(a), b = std::log10(b);
      }
      if (std::min(a, b) < x0 - tolerance || std::max(a, b) > x1 + tolerance)
        return Form("a %s from x = %g to %g leaves its frame (%g to %g)", o->ClassName(), pad->GetLogx() ? std::pow(10, a) : a,
                    pad->GetLogx() ? std::pow(10, b) : b, pad->GetLogx() ? std::pow(10, x0) : x0,
                    pad->GetLogx() ? std::pow(10, x1) : x1);
    }
    if (pad != &canvas && pad->GetXlowNDC() < 1e-3 && pad->GetWNDC() > 0.999 && pad->GetHNDC() < 0.99)
      stacked.push_back(pad);
  }
  for (size_t i = 1; i < stacked.size(); ++i) {
    const double tolerance = 1e-6 * std::fabs(stacked[0]->GetUxmax() - stacked[0]->GetUxmin());
    if (std::fabs(stacked[i]->GetUxmin() - stacked[0]->GetUxmin()) > tolerance ||
        std::fabs(stacked[i]->GetUxmax() - stacked[0]->GetUxmax()) > tolerance)
      return Form("its stacked pads do not share one x range (%g to %g and %g to %g)", stacked[0]->GetUxmin(),
                  stacked[0]->GetUxmax(), stacked[i]->GetUxmin(), stacked[i]->GetUxmax());
  }
  return "";
}

void save(TCanvas &c, const std::string &relative, const PlotJson &weights = PlotJson(),
          const char *function = __builtin_FUNCTION(), const char *source = "plots/IANFigures.C", const char* family=nullptr) {
#ifdef PHOTONJET_NATIVE_CAPTURE_TEST
  if(gNativeCaptureSink) {
    gNativeCaptureSink(c,relative,(std::string(source).find("QAFigures")!=std::string::npos?"qa.":"ian.")+std::string(family?family:function));
    clearFigureState();gWritten.push_back(relative);return;
  }
#endif
  const PlotJson physicsWeights = checkedWeights(weights, relative);
  const std::string path = gFig + "/" + relative;
  gSystem->mkdir(gSystem->GetDirName(path.c_str()).Data(), true);
  c.Update();
  fitLegends(c);
  c.cd();
  c.Update();
  std::vector<std::string> visible;
  PlotJson populations = PlotJson::array();
  canvasText(c, visible, &populations);  /* a mislabelled candidate curve is refused before any image is written */
  if (!weights.is_null() && weights.value("role", "") == "model_training_diagnostic" &&
      std::none_of(visible.begin(), visible.end(), [](const std::string &text) {
        return text.find("Training diagnostic") != std::string::npos;
      }))
    throw std::runtime_error("training diagnostic must be visibly labelled: " + relative);
  if (const std::string defect = frameDefect(c); !defect.empty()) {  /* a broken drawing is never written */
    missing(relative, "refused by the frame check: " + defect);
    return;
  }
  photonjet_science::Collector physicsCollector;
  const PlotJson physics = physicsCollector.finish(c, relative,
      std::string(source).find("QAFigures") != std::string::npos ? "qa." + std::string(family?family:function) : "ian." + std::string(family?family:function),
      photonjet_science::annotationRosterDeclared((std::string(source).find("QAFigures") != std::string::npos ? "qa." : "ian.") + std::string(family?family:function)));
  writeImage(c, path);
  std::ifstream saved(path, std::ios::binary | std::ios::ate);
  if (!saved || saved.tellg() <= 0) throw std::runtime_error("figure was not saved: " + path);
  bool compiled = !gSemanticPackages.empty() && !gLegacy;
  PlotJson semantics = PlotJson::object();
  for (const auto &key : gSemanticPackages) {
    const PlotJson *record = packageRecord(key);
    const bool valid = record && !record->at("selection").is_null();
    compiled = compiled && valid;
    if (valid) semantics[key] = *record;
  }
  for (const auto &t : gVisibleText)
    if (std::find(visible.begin(), visible.end(), t) == visible.end()) visible.push_back(t);
  PlotJson receipt = {
      {"schema", "PhotonJetRenderReceiptV1"},
      {"image_sha256", photonjet_render::sha256(path, true)},
      {"physics", physics},
      {"scientific_definition_origins", physicsCollector.origins},
      {"native_source_witness", fa::nativeProducerWitness()},
      {"native_source_ancestry", fa::nativeProducerAncestry()},
      {"plot_recipe",
       {{"schema", "PhotonJetPlotRecipeV1"},
        {"role", "canonical"},
        {"id", "canonical"},
        {"output", relative},
        {"baseline_output", relative},
        {"source", source},
        {"function", function}}},
      {"inputs", gRenderInputs},
      {"style_and_renderer", gStylePins},
      {"visible_text", visible},
      {"candidate_populations", populations},
      {"selection", semantics},
      {"display_projections", gProjectionSteps},
      {"input_observations", {{"schema", "PlotInputObservationsV1"}, {"rows", gInputObservations}}},
      {"annotation_status", compiled ? "COMPILED_METADATA" : "LEGACY_UNVERIFIED_ANNOTATION"},
      {"input_coverage", gPopulation},
      {"input_notes", gInputNotes},
      {"layout_review", "REPRESENTATIVE_VISUAL_QA_REQUIRED"}};
  if (!gSystematicCoverage.is_null()) receipt["systematic_coverage"] = gSystematicCoverage;
  /* A physics figure (measurement or QA) records the weight records it was checked against, one per package. */
  if (weights.is_null()) receipt["weight_provenance"] = physicsWeights;
  if (std::string(function) == "combinatoric")
    receipt["display_axes"] = {{"x_min", gFrame->GetXaxis()->GetXmin()}, {"x_max", gFrame->GetXaxis()->GetXmax()},
                               {"y_min", gFrameLow}, {"y_max", gFrameHigh}};
  if (!gCombinatoricRecipe.is_null())
    for (const auto &[key, value] : gCombinatoricRecipe.items()) receipt["plot_recipe"][key] = value;
  photonjet_render::write(path + ".json", receipt);
  clearFigureState();
  gWritten.push_back(relative);
}

/* ---- series ------------------------------------------------------------------------------------------------------ */
/* The look of a series: its colour, its marker and marker size, and lines of width 2. The markers by ROOT's numbers:
   20 full circle, 21 full square, 22 full triangle, 33 full diamond, 34 full cross; 24, 25 and 27 the open circle,
   square and diamond. */
template <class T> void style(T *h, int color, int marker, double size = 1.2) {
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(size);
  h->SetLineWidth(2);
}

/* The pairs of a (photon pT, xJ) histogram summed over the photon bins, per unit xJ, divided by `photons` when
   given. */
TH1D *xjProjection(const TH2 *pairs, double photons = 0, const TH1* photonSource=nullptr) {
  auto *h = pairs->ProjectionY(fresh().c_str(), 1, pairs->GetNbinsX(), "e");
  h->SetDirectory(nullptr);
  h->Scale(photons > 0 ? 1.0 / photons : 1.0, "width");
  photonjet_science::attachTrace(h,photonjet_science::projectionOf(pairs,PlotJson::array()).then("project_y_sum_x",{{"x_first",1},{"x_last",pairs->GetNbinsX()},{"errors","ROOT ProjectionY e"}})
      .then("divide_by_bin_width").then("divide_by_photons",{{"applied",photons>0},{"photons",photons}}));
  if(photonSource)photonjet_science::attachTrace(h,photonjet_science::derivedTrace("photon_normalization",{h,photonSource},{{"photons",photons}}));
  photonjet_science::setTag(h,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"ProjectionY propagates diagonal errors; full covariance is retained in upstream inputs"));
  return h;
}

/* The bins of a distribution inside [low, high] as points, each with a horizontal bar over its bin and its bin error as
   the vertical bar. (The QA figures build their points, ratios and step lines from their own type Bins, with markers,
   over and line in QAFigures.C.) */
TGraphErrors *points(const TH1 *h, double low, double high, int color, int marker, double size = 1.2, const char *role = __builtin_FUNCTION()) {
  auto *g = new TGraphErrors;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double a = h->GetXaxis()->GetBinLowEdge(i), b = h->GetXaxis()->GetBinUpEdge(i);
    if (!binInside(a, b, low, high)) continue;
    const int n = g->GetN();
    g->SetPoint(n, 0.5 * (a + b), h->GetBinContent(i));
    g->SetPointError(n, 0.5 * (b - a), h->GetBinError(i));
  }
  style(g, color, marker, size);
  photonjet_science::bind(g,photonjet_science::seriesRole(h,role),photonjet_science::selectedTrace(h,low,high),photonjet_science::selectedHistogram(h,low,high),photonjet_science::covarianceOf(h));
  return g;
}

/* A distribution as a histogram of steps over the bins inside [low, high] and zero elsewhere: drawn, it rises from
   zero at `low` and returns to zero at `high`. jointRefold draws the refolded pairs with it. */
TH1D *steps(const TH1 *h, double low, double high, int color, const char *role = __builtin_FUNCTION()) {
  auto *s = static_cast<TH1D *>(h->Clone(fresh().c_str()));
  s->SetDirectory(nullptr);
  for (int i = 0; i <= s->GetNbinsX() + 1; ++i) {
    const double a = s->GetXaxis()->GetBinLowEdge(i), b = s->GetXaxis()->GetBinUpEdge(i);
    if (i == 0 || i > s->GetNbinsX() || !binInside(a, b, low, high)) s->SetBinContent(i, 0), s->SetBinError(i, 0);
  }
  s->GetXaxis()->SetRangeUser(low + 1e-6, high - 1e-6);
  s->SetLineColor(color);
  s->SetLineWidth(2);
  photonjet_science::bind(s,photonjet_science::seriesRole(h,role),photonjet_science::selectedTrace(h,low,high),photonjet_science::selectedHistogram(h,low,high),photonjet_science::covarianceOf(h));
  return s;
}

/* The ratio a / b bin by bin (errors of both), as points over the bins inside [low, high]. */
TGraphErrors *ratioPoints(const TH1 *a, const TH1 *b, double low, double high, int color, int marker) {
  auto *r = static_cast<TH1 *>(a->Clone(fresh().c_str()));
  r->SetDirectory(nullptr);
  r->Divide(b);
  photonjet_science::attachTrace(r,photonjet_science::Trace::combine("ROOT_TH1_ratio",{photonjet_science::selectedTrace(a,low,high),photonjet_science::selectedTrace(b,low,high)},{{"errors","independent diagonal propagation"},{"zero_denominator","ROOT returns zero content and error"}}));
  photonjet_science::setTag(r,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"TH1 Divide has no cross-bin covariance propagation"));
  return points(r, low, high, color, marker);
}

/* The steps of a distribution over the bins inside [low, high] as a line of width 2: two points per bin, at its edges.
   toZero: the line rises from zero at the first edge and falls back to zero at the last. */
TGraph *stepLine(const TH1 *h, double low, double high, int color, bool toZero, const char *role = __builtin_FUNCTION()) {
  auto *g = new TGraph;
  double last = low;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double a = h->GetXaxis()->GetBinLowEdge(i), b = h->GetXaxis()->GetBinUpEdge(i);
    if (!binInside(a, b, low, high)) continue;
    if (toZero && g->GetN() == 0) g->SetPoint(0, a, 0);
    g->SetPoint(g->GetN(), a, h->GetBinContent(i));
    g->SetPoint(g->GetN(), b, h->GetBinContent(i));
    last = b;
  }
  if (toZero && g->GetN() > 0) g->SetPoint(g->GetN(), last, 0);
  g->SetLineColor(color);
  g->SetLineWidth(2);
  photonjet_science::bind(g,photonjet_science::seriesRole(h,role),photonjet_science::selectedTrace(h,low,high),photonjet_science::selectedHistogram(h,low,high),photonjet_science::covarianceOf(h));
  return g;
}

/* A model as a step line over the bins inside [low, high] only (no drop to zero outside the reported range). */
TGraph *outline(const TH1 *h, double low, double high, int color, int lineStyle) {
  TGraph *g = stepLine(h, low, high, color, false);
  g->SetLineStyle(lineStyle);
  return g;
}

/* A background estimate as the IAN draws it: steps over the bins inside [low, high] that rise from zero at the first
   edge and fall back to zero at the last, with no line along the axis; shaded under the steps when `filled`. */
TGraph *estimateSteps(const TH1 *h, double low, double high, int color, bool filled) {
  TGraph *g = stepLine(h, low, high, color, true);
  if (filled) g->SetFillColorAlpha(color, 0.18);
  return g;
}

/* ---- the xJ frame, its axes and its marks ------------------------------------------------------------------------ */
/* The titles of an xJ axis and of the two y axes over it: pairs per unit xJ, and the same per photon. */
const char *kXjTitle = "x_{J#gamma}";
const char *kPairsTitle = "dN_{pair}/#scale[1]{dx_{J#gamma}}";
const char *kPerPhotonTitle = "(1/N_{#gamma}) dN_{pair}/#scale[1]{dx_{J#gamma}}";
/* The end of an xJ axis: the upper edge of the reported range rounded up to 0.1. */
double xjEnd(double high) { return std::ceil(high * 10 - 1e-9) / 10; }

/* The left edge of every xJ axis (Justin, 2026-10-15: every xJ figure starts where the final result starts): the
   lower edge of the reported range, where the kinematic turn-on ends, less 0.05, down to a multiple of 0.05, the rule
   of resultFrame. Set once per drawing (drawAll) from the common lower edge of the final-result figure's results;
   0 until then. */
double gXjShownLow = 0;
double xjShownLow(double low) { return low > 0 ? std::floor((low - 0.05) * 20 + 1e-9) / 20 : 0; }

/* xJ ticks every 0.3 when the axis end is a multiple of 0.3, as in the IAN (0, 0.3, ..., 1.5); otherwise ROOT's.
   SetNdivisions(n1 + 100 n2) asks for n1 primary divisions, each with n2 secondary ones; `false` takes the numbers as
   they are, where ROOT would otherwise round the steps. */
void xjTicks(TH1 *frame, double high) {
  const double end = xjEnd(high), divisions = (end - gXjShownLow) / 0.3;
  if (std::fabs(divisions - std::round(divisions)) < 1e-6 && divisions <= 10)
    frame->GetXaxis()->SetNdivisions(int(std::round(divisions)) + 3 * 100, false);
  else
    frame->GetXaxis()->SetNdivisions(505);
}

/* One series of an xJ figure and its row of the legend. object: a graph (points, or a line with a draw option) or a
   histogram (drawn as steps); null for a row of text only. text, legend: the legend text and the legend option. draw:
   the draw option of a line ("L", "LF"); none for points, which are drawn with "P" on top of the lines. */
struct Series { TObject *object; std::string text; const char *legend; const char *draw = nullptr; };

/* An xJ frame in the IAN's layout with its series: histograms and lines under the points, the legend in series order.
   upper: the frame of the pad above a ratio pad, without x labels and x title. Returns the legend, which is drawn. */
TLegend *xjFrame(const std::vector<Series> &series, double high, const char *ytitle,
                 const std::vector<std::string> &lines, bool upper = false, double low = 0) {
  std::vector<Extent> drawn;
  std::vector<std::string> labels;
  for (const auto &s : series) {
    if (auto *g = dynamic_cast<TGraph *>(s.object)) extents(drawn, g);
    else if (auto *h = dynamic_cast<TH1 *>(s.object)) extents(drawn, h);
    labels.push_back(s.text);
  }
  TLegend *leg = ianFrame(gXjShownLow, xjEnd(high), low, drawn, false, upper ? "" : kXjTitle, ytitle, lines, labels);
  if (!gCombinatoricFrame.is_null()) {
    gFrameLow = gCombinatoricFrame.at(0); gFrameHigh = gCombinatoricFrame.at(1);
    if (!(gFrameHigh > gFrameLow)) throw std::runtime_error("invalid comparison frame");
    gFrame->SetMinimum(gFrameLow); gFrame->SetMaximum(gFrameHigh);
  }
  xjTicks(gFrame, high);
  if (upper) upperAxes(gFrame);
  for (const auto &s : series) if (dynamic_cast<TH1 *>(s.object)) photonjet_science::draw(s.object,s.draw ? s.draw : "HIST SAME");
  for (const auto &s : series) if (dynamic_cast<TGraph *>(s.object) && s.draw) photonjet_science::draw(s.object,s.draw);
  for (const auto &s : series) if (dynamic_cast<TGraph *>(s.object) && !s.draw) photonjet_science::draw(s.object,"P");
  for (const auto &s : series) leg->AddEntry(s.object, s.text.c_str(), s.legend);
  photonjet_science::draw(leg);
  return leg;
}

/* The legend symbol of a point series: the marker, on its bin bar when the points carry one and with a vertical bar
   when they carry errors. purityCorrection, combinatoric, beforeAfter, detectorLevel and priorReweighting take their
   legend symbols from it, so their legends show which series carry bars; the other figures name the legend option
   themselves. */
const char *pointKey(const TGraphErrors *g) {
  bool widths = false, errors = false;
  for (int i = 0; i < g->GetN(); ++i) widths = widths || g->GetEX()[i] > 0, errors = errors || g->GetEY()[i] > 0;
  return widths ? (errors ? "lpe" : "lp") : (errors ? "pe" : "p");
}

/* The outline of an open marker in canvas pixels, as heavy as the bars of its points, and the size of the white marker
   that leaves this outline inside a filled marker of the style `filled` and the size `size`. A marker of size 1 is
   8 px high; the outline is half the difference of the two heights for a circle or a square and, the edges of the
   diamond being slanted, 1/3.6 of it for the diamond. */
const double kOutlinePx = 1.8;
double whiteSize(int filled, double size) { return size - (filled == kFullDiamond ? 3.6 : 2.0) * kOutlinePx / 8; }

/* The legend symbol of an open series: on every row of `legend` (inside the frame of this pad) that shows `open`, a
   filled marker of the style `filled` in the series colour and the smaller white one in it. */
void openKey(const TLegend *legend, const TGraphErrors *open, int filled) {
  gPad->Update();  /* the legend takes its place in the pad, and the pad's coordinates the frame's final range */
  const double width = (legend->GetX2NDC() - legend->GetX1NDC()) / legend->GetNColumns();
  const double height = (legend->GetY2NDC() - legend->GetY1NDC()) / legend->GetNRows();
  int at = 0;
  TIter next(legend->GetListOfPrimitives());
  while (auto *entry = static_cast<TLegendEntry *>(next())) {
    const int column = at % legend->GetNColumns(), row = at / legend->GetNColumns();
    ++at;
    if (entry->GetObject() != open) continue;
    const double nx = legend->GetX1NDC() + (column + 0.5 * legend->GetMargin()) * width;
    const double ny = legend->GetY2NDC() - (row + 0.5) * height;
    const double x = gPad->PadtoX(gPad->GetX1() + nx * (gPad->GetX2() - gPad->GetX1()));
    const double y = gPad->PadtoY(gPad->GetY1() + ny * (gPad->GetY2() - gPad->GetY1()));
    for (bool white : {false, true}) {
      auto *mark = new TGraph(1, &x, &y);
      mark->SetMarkerStyle(filled);
      mark->SetMarkerSize(white ? whiteSize(filled, open->GetMarkerSize()) : open->GetMarkerSize());
      mark->SetMarkerColor(white ? int(kWhite) : int(open->GetMarkerColor()));
      photonjet_science::decoration(mark,"legend marker overlay");
      photonjet_science::draw(mark,"P");
    }
  }
}

/* Open markers as the IAN draws them, white inside an outline as heavy as their bars, for a series that is drawn with
   open markers: at every point of `open` a filled marker of the style `filled` in the series colour and the smaller
   white one in it, over the bars drawn so far, and the same on the legend rows of the series in this pad; then the
   markers of `front`, which stay on top. (whiteCentres is the other way to an open marker: for a series drawn with
   full markers, with the size of the white marker given at each call.) */
void whiteInside(const TGraphErrors *open, int filled, const std::vector<TGraphErrors *> &front) {
  for (bool white : {false, true}) {
    auto *marks = new TGraph(open->GetN(), open->GetX(), open->GetY());
    marks->SetMarkerStyle(filled);
    marks->SetMarkerSize(white ? whiteSize(filled, open->GetMarkerSize()) : open->GetMarkerSize());
    marks->SetMarkerColor(white ? int(kWhite) : int(open->GetMarkerColor()));
    photonjet_science::decoration(marks,"white marker centre overlay");
    photonjet_science::draw(marks,"P");
  }
  std::vector<const TLegend *> legends;
  TIter next(gPad->GetListOfPrimitives());
  while (auto *object = next()) if (auto *legend = dynamic_cast<TLegend *>(object)) legends.push_back(legend);
  for (const auto *legend : legends) openKey(legend, open, filled);
  for (auto *g : front) photonjet_science::draw(g,"PX");
}

/* The upper end of a y axis from zero with `divisions` primary divisions, moved up where a tick label would sit on
   the frame's top edge: the last label stays a twentieth of the axis below the end. The tick step is the one the
   axis takes: the range over the divisions rounded up to 1, 2, 5 or 10, with more divisions tried while the steps
   would number half of those asked for. This is the rule of the xJ frames (headroom) and of the reco-to-truth row;
   the figures of the results and the Delta-phi grid end half a division above their last label (labelledTicks), the
   frames opened under a header a little above it (openUnderHeader), the closure and refolding ratios a little beyond
   their outer labels (kBeyondLabels). */
double clearTop(double high, int divisions) {
  for (int pass = 0; pass < 3; ++pass) {
    double step = high;
    for (int n = std::max(divisions, 2);; ++n) {
      const double raw = high / n;
      const int power = static_cast<int>(std::log10(raw)) - (raw <= 1 ? 1 : 0);
      const double mantissa = raw * std::pow(10, -power) - 1e-10;
      step = (mantissa <= 1 ? 1 : mantissa <= 2 ? 2 : mantissa <= 5 ? 5 : 10) * std::pow(10, power);
      if (divisions <= 5 || 2 * static_cast<int>(high / step + 1.00001) != divisions) break;
    }
    const double last = step * std::floor(high / step + 1e-6);
    if (high - last >= 0.049 * high) break;
    high = last / 0.95;
  }
  return high;
}

/* Room above the points as the IAN leaves it: the y range of the frame that ianFrame drew last (gFrame, whatever the
   current pad), a tenth more than its header and its legend need, with `divisions` primary divisions and no tick label
   on the frame's top edge. */
void headroom(int divisions) {
  gFrameHigh = clearTop(1.1 * gFrameHigh, divisions);
  gFrame->SetMaximum(gFrameHigh);
  gFrame->GetYaxis()->SetNdivisions(500 + divisions);
  gPad->Modified();
}

/* The axis text of an xJ figure as the IAN sets it: the tick labels at the size of the header lines and a little off
   the frame, so the two zeros of the origin stand apart; the titles larger and, when `centred`, in the middle of
   their axes (at the axis ends otherwise); ticks of `tickPx` on both axes. xTitlePx and yTitlePx: the distance of a
   title from its axis in pixels (xTitlePx 0: no x labels and no x title, for the pad above a ratio pad); yTitleSize:
   the size of the y title. A title offset counts in title sizes: ROOT's axis painter sets a title 1.3 sizes from an x
   axis and 1.6 sizes from a y axis per unit of offset.
   kXjTitlePx: the size of the axis titles. kRatioTitlePx: the y title of a ratio pad. kXTitlePx, kYTitlePx: the
   distances of the titles from their axes. */
const double kXjTitlePx = 29.5, kRatioTitlePx = kSmallPx, kXTitlePx = 38, kYTitlePx = 66;
void xjAxes(TH1 *frame, double xTitlePx, double yTitlePx, bool centred = true, double yTitleSize = kXjTitlePx,
            double tickPx = kTickPx) {
  const double w = padWidthPx(), h = padHeightPx();
  TAxis *x = frame->GetXaxis(), *y = frame->GetYaxis();
  pixelAxes(frame, kLinePx, kXjTitlePx);
  x->CenterTitle(centred), y->CenterTitle(centred);
  y->SetTitleSize(yTitleSize);
  tickLengths(frame, tickPx);
  x->SetLabelOffset(8 / h), y->SetLabelOffset(7 / w);
  x->SetTitleOffset(xTitlePx / (1.3 * kXjTitlePx)), y->SetTitleOffset(yTitlePx / (1.6 * yTitleSize));
  if (!(xTitlePx > 0)) hideXText(frame);
}

/* The title of a lower (ratio) pad, defined with the closure figures below (patch 49). */
void lowerTitles(TH1 *frame, double target = 0);

/* The margins of a single xJ frame on the current pad. left: the left margin. */
void xjPad(double left = 0.14) { gPad->SetMargin(left, 0.03, 0.16, 0.05); }

/* The frame of the two figures of pair counts once xjFrame has drawn it (gFrame): the IAN's headroom and axis text, the
   counts written in full, and the y title clear of five-digit labels. kCountLeft: the left margin of both figures, wide
   enough for these labels; kCountYTitlePx: the distance of the y title from its axis. kDiamondSize: the diamonds of
   the purity-corrected distribution in both figures, larger than the circles of the other points. */
const double kDiamondSize = 2.4, kCountLeft = 0.17, kCountYTitlePx = 92;
void countAxes() {
  headroom(6);
  xjAxes(gFrame, kXTitlePx, kCountYTitlePx);
  gFrame->GetYaxis()->SetMaxDigits(5);
}

/* ---- the corrections, xJ before and after unfolding, and the checks of the simulation --------------------------- */
/* Photon-purity correction of the reconstructed xJ distribution, in pairs per unit xJ: region A, the transferred
   region-C background alpha C and the leakage-corrected distribution, as in the IAN. The circles of region A stay on
   top of the larger diamonds, so both are seen where they coincide. Reads pairs_A, pairs_C_scaled and
   pairs_purity_corrected. */
void purityCorrection(TFile &f, const Result &r, const std::string &name) {
  auto *pairsA = dynamic_cast<TH2 *>(get(f, r, "pairs_A"));
  auto *pairsCScaled = dynamic_cast<TH2 *>(get(f, r, "pairs_C_scaled"));
  auto *purityCorrected = dynamic_cast<TH2 *>(get(f, r, "pairs_purity_corrected"));
  if (!pairsA || !pairsCScaled || !purityCorrected)
    return missing(name, "needs pairs_A, pairs_C_scaled and pairs_purity_corrected");
  TGraphErrors *regionA = points(xjProjection(pairsA), r.xjLow, r.xjHigh, kBlack, 20);
  TGraphErrors *corrected = points(xjProjection(purityCorrected), r.xjLow, r.xjHigh, gBlue, 33, kDiamondSize);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  xjPad(kCountLeft);
  xjFrame({{estimateSteps(xjProjection(pairsCScaled), r.xjLow, r.xjHigh, gVermilion, false), "#alphaC", "l", "L"},
           {regionA, "Region A", pointKey(regionA)},
           {corrected, "Leakage-corrected", pointKey(corrected)}},
          r.xjHigh, kPairsTitle, xjHeader(f, r));
  countAxes();
  photonjet_science::draw(regionA,"PX");
  save(cv, name);
}

/* Combinatorial-jet subtraction in pairs per unit xJ, in every collision system: the scaled combinatorial estimate,
   the purity-corrected distribution and the distribution after subtraction, as in the IAN, in the frame of the
   purity-correction figure. The open diamonds are larger than the circles, so they are seen around them where the
   two coincide. Reads pairs_purity_corrected, pairs_combinatoric and pairs_corrected. Computed here: whether the legend
   names a subtraction, from the integral of the saved estimate over all of its cells (not from the result's
   combinatoric setting). */
void combinatoric(TFile &f, const Result &r, const std::string &name) {
  auto *purityCorrected = dynamic_cast<TH2 *>(get(f, r, "pairs_purity_corrected"));
  auto *combinatoricPairs = dynamic_cast<TH2 *>(get(f, r, "pairs_combinatoric"));
  auto *correctedPairs = dynamic_cast<TH2 *>(get(f, r, "pairs_corrected"));
  if (!purityCorrected || !combinatoricPairs || !correctedPairs)
    return missing(name, "needs pairs_purity_corrected, pairs_combinatoric and pairs_corrected");
  TGraphErrors *purity = points(xjProjection(purityCorrected), r.xjLow, r.xjHigh, gBlue, 27, kDiamondSize);
  TGraphErrors *subtracted = points(xjProjection(correctedPairs), r.xjLow, r.xjHigh, kBlack, 20);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  xjPad(kCountLeft);
  auto lines = xjHeader(f, r);
  if (!gCombinatoricRecipe.is_null())
    for (const auto &line : gCombinatoricRecipe.at("header")) lines.push_back(line.get<std::string>());
  auto *estimate = xjProjection(combinatoricPairs);
  double minimum = 0;
  bool nonzero = false;
  for (int i = 1; i <= estimate->GetNbinsX(); ++i)
    if (binInside(estimate->GetXaxis()->GetBinLowEdge(i), estimate->GetXaxis()->GetBinUpEdge(i), r.xjLow, r.xjHigh)) {
      nonzero = nonzero || estimate->GetBinContent(i) != 0;
      minimum = std::min(minimum, estimate->GetBinContent(i));
    }
  /* A negative fitted background is a diagnostic failure, not an absent subtraction. Preserve its sign visibly. */
  xjFrame({{estimateSteps(estimate, r.xjLow, r.xjHigh, gVermilion, true),
            minimum < 0 ? "Signed combinatorial estimate" : nonzero ? "Scaled combinatorial" : "No combinatorial subtraction", "l", "LF"},
           {purity, "Purity-corrected", pointKey(purity)},
           {subtracted, "After subtraction", pointKey(subtracted)}},
          r.xjHigh, kPairsTitle, lines, false, minimum < 0 ? minimum * 1.15 : 0);
  countAxes();
  whiteInside(purity, 33, {subtracted});
  save(cv, name);
}

/* The two pads of an xJ figure with a ratio panel, in the IAN's proportions: the ratio frame a third of the upper one
   and tall enough for its title, the two frames a small gap apart, and the margins in pixels. */
std::pair<TPad *, TPad *> xjRatioPads(TCanvas &cv) {
  auto [top, bottom] = ratioPads(cv);
  const double split = 0.34, canvasPx = top->GetWh();
  top->SetPad(0, split, 1, 1);
  bottom->SetPad(0, 0, 1, split);
  for (TPad *pad : {top, bottom}) pad->SetLeftMargin(0.14), pad->SetRightMargin(0.025);
  top->SetBottomMargin(8 / (canvasPx * (1 - split)));
  bottom->SetTopMargin(8 / (canvasPx * split));
  bottom->SetBottomMargin(80 / (canvasPx * split));
  return {top, bottom};
}

/* Per-photon xJ before (purity- and combinatoric-corrected) and after unfolding, and their ratio, as in the IAN. The
   ratio is on one scale in every collision system (0 to 2.5, more only when a ratio is larger) and, as in the IAN,
   without vertical bars: both spectra come from the same pairs, and each carries its own bars above. Reads
   xj_per_photon_reco and xj_per_photon. */
void beforeAfter(TFile &f, const Result &r, const std::string &name) {
  auto *recoXj = get(f, r, "xj_per_photon_reco");
  auto *unfoldedXj = get(f, r, "xj_per_photon");
  if (!unfoldedXj || !recoXj) return missing(name, "needs xj_per_photon and xj_per_photon_reco");
  TGraphErrors *reco = points(recoXj, r.xjLow, r.xjHigh, kBlack, 24);
  /* diamonds of size 2.0, smaller than the kDiamondSize of the figures of pair counts */
  TGraphErrors *unfolded = points(unfoldedXj, r.xjLow, r.xjHigh, gBlue, 33, 2.0);
  TCanvas cv("cv", "", kCanvasW, kRatioCanvasH);
  auto [top, bottom] = xjRatioPads(cv);
  top->cd();
  xjFrame({{reco, "Corrected reconstructed", pointKey(reco)}, {unfolded, "Unfolded", pointKey(unfolded)}},
          r.xjHigh, kPerPhotonTitle, xjHeader(f, r), true);
  headroom(8);
  xjAxes(gFrame, 0, kYTitlePx);
  whiteInside(reco, 20, {unfolded});
  bottom->cd();
  TGraphErrors *ratio = ratioPoints(unfoldedXj, recoXj, r.xjLow, r.xjHigh, gBlue, 20);
  double ratioTop = 2.5;
  for (int i = 0; i < ratio->GetN(); ++i) {
    ratio->SetPointError(i, ratio->GetEX()[i], 0);
    ratioTop = std::max(ratioTop, 1.15 * ratio->GetY()[i]);
  }
  photonjet_science::copyBinding(ratio,ratio); // hidden ratio bars are display only
  TH1F *lower = ratioFrame(gXjShownLow, xjEnd(r.xjHigh), 0, ratioTop, kXjTitle, "Unf. / reco");
  xjTicks(lower, r.xjHigh);
  xjAxes(lower, kXTitlePx, kYTitlePx, true, kRatioTitlePx);
  lowerTitles(lower, kXjTitlePx);
  frameReferenceLine(1);
  photonjet_science::draw(ratio,"P");
  save(cv, name);
}

/* The statistical uncertainty of the per-photon xJ stage by stage, in percent of the value of each bin: the region-A
   pairs (Poisson), after the purity correction (the region-C subtraction adds its own counts), after the combinatorial
   subtraction (the sideband counts add theirs), and the unfolded result from the data toys alone and with the
   response toys, the error the result is drawn with. The first three are summed over the photon bins. Reads pairs_A,
   pairs_purity_corrected, pairs_corrected, xj_half_width_data and xj_per_photon. Computed here: error / value per bin;
   a bin without a positive value is left out of its line. */
void statisticalBudget(TFile &f, const Result &r, const std::string &name) {
  auto *pairsA = dynamic_cast<TH2 *>(get(f, r, "pairs_A"));
  auto *purityCorrected = dynamic_cast<TH2 *>(get(f, r, "pairs_purity_corrected"));
  auto *correctedPairs = dynamic_cast<TH2 *>(get(f, r, "pairs_corrected"));
  auto *dataHalf = get(f, r, "xj_half_width_data");
  auto *unfoldedXj = get(f, r, "xj_per_photon");
  if (!pairsA || !purityCorrected || !correctedPairs || !dataHalf || !unfoldedXj)
    return missing(name, "needs pairs_A, pairs_purity_corrected, pairs_corrected, xj_half_width_data and xj_per_photon");
  /* errors: the bin errors of `values`, or the contents of `half` (a half-width histogram on the same bins) */
  auto percent = [&](const TH1 *values, const TH1 *half, int color, int lineStyle) {
    auto *g = new TGraph;
    for (int i = 1; i <= values->GetNbinsX(); ++i) {
      const double a = values->GetXaxis()->GetBinLowEdge(i), b = values->GetXaxis()->GetBinUpEdge(i);
      const double v = values->GetBinContent(i), e = half ? half->GetBinContent(i) : values->GetBinError(i);
      if (!binInside(a, b, r.xjLow, r.xjHigh) || !(v > 0)) continue;
      g->SetPoint(g->GetN(), a, 100 * e / v);
      g->SetPoint(g->GetN(), b, 100 * e / v);
    }
    g->SetLineColor(color);
    g->SetLineStyle(lineStyle);
    g->SetLineWidth(2);
    std::vector<const TObject*> inputs{values};if(half)inputs.push_back(half);
    return photonjet_science::bindNative(g,half?"relative_error/data_toys":"relative_error/stored_variance",
      photonjet_science::derivedTrace("relative_statistical_uncertainty_percent",inputs,{{"low",r.xjLow},{"high",r.xjHigh},{"positive_denominator",true}}),
      photonjet_science::covariance("not_computed",{},0,"no covariance of relative uncertainty curve"));
  };
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  xjPad();
  xjFrame({{percent(xjProjection(pairsA), nullptr, gGray, 2), "Region A pairs", "l", "L"},
           {percent(xjProjection(purityCorrected), nullptr, gBlue, 1), "Purity-corrected", "l", "L"},
           {percent(xjProjection(correctedPairs), nullptr, gVermilion, 1), "Combinatorial-subtracted", "l", "L"},
           {percent(unfoldedXj, dataHalf, kBlack, 2), "Unfolded, data", "l", "L"},
           {percent(unfoldedXj, nullptr, kBlack, 1), "Unfolded, data and response", "l", "L"}},
          r.xjHigh, "Relative statistical uncertainty [%]", xjHeader(f, r));
  headroom(6);
  xjAxes(gFrame, kXTitlePx, kYTitlePx);
  save(cv, name);
}

/* The mean of a per-unit-xJ distribution over [low, high] and its error from the bin errors (bins independent).
   Computed here, for the header of detectorLevel. */
std::pair<double, double> meanXj(const TH1 *h, double low, double high) {
  double sum = 0, moment = 0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double a = h->GetXaxis()->GetBinLowEdge(i), b = h->GetXaxis()->GetBinUpEdge(i);
    if (!binInside(a, b, low, high)) continue;
    sum += h->GetBinContent(i) * (b - a), moment += h->GetBinContent(i) * (b - a) * 0.5 * (a + b);
  }
  if (!(sum > 0)) return {0.0, 0.0};
  const double mean = moment / sum;
  double variance = 0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double a = h->GetXaxis()->GetBinLowEdge(i), b = h->GetXaxis()->GetBinUpEdge(i);
    if (!binInside(a, b, low, high)) continue;
    variance += std::pow((0.5 * (a + b) - mean) * (b - a) * h->GetBinError(i) / sum, 2);
  }
  return {mean, std::sqrt(variance)};
}

/* Detector-level check of the simulation: the corrected reconstructed DATA pairs per photon against the simulation's
   expected ones (response column sums plus detector and boundary fakes: what the unfolding takes the measured cells
   to be), per photon and unit xJ, and DATA / simulation on one scale in every collision system (0 to 2, more only
   when a ratio with its bar is larger). A difference here is a difference of the jet energy relative to the photon
   energy, or of the modelled shape, that the unfolding carries into the result; in Au+Au it also contains the medium
   modification. The legend names the spectrum as the figure before and after unfolding does, and the ratio pad
   carries no label at its top edge, next to the upper pad's zero.
   Reads pairs_corrected, pairs_reco_expected, photons_corrected and photon_reco. Computed here, not read: the
   simulation's curve, in which each photon bin of pairs_reco_expected is divided by the simulated photons of the bin
   (photon_reco) and weighted by the DATA photon yield of that bin, so the two share one photon spectrum; and the ratio
   of the mean xJ of the two in the header, whose error is the DATA's only. */
void detectorLevel(TFile &f, const Result &r, const std::string &name) {
  auto *dataPairs = dynamic_cast<TH2 *>(get(f, r, "pairs_corrected"));
  auto *expectedPairs = dynamic_cast<TH2 *>(get(f, r, "pairs_reco_expected"));
  auto *dataPhotons = get(f, r, "photons_corrected");
  auto *simulatedPhotons = get(f, r, "photon_reco");
  if (!dataPairs || !expectedPairs || !dataPhotons || !simulatedPhotons)
    return missing(name, "needs pairs_corrected, pairs_reco_expected, photons_corrected and photon_reco");
  if (!(dataPhotons->Integral() > 0)) return missing(name, "no corrected photons");
  TH1D *data = xjProjection(dataPairs, dataPhotons->Integral(),dataPhotons);
  auto *simulation = static_cast<TH1D *>(data->Clone(fresh().c_str()));
  simulation->SetDirectory(nullptr);
  simulation->Reset();
  for (int x = 1; x <= expectedPairs->GetNbinsY(); ++x) {
    double sum = 0;
    for (int i = 1; i <= expectedPairs->GetNbinsX(); ++i)
      if (simulatedPhotons->GetBinContent(i) > 0)
        sum += expectedPairs->GetBinContent(i, x) / simulatedPhotons->GetBinContent(i) * dataPhotons->GetBinContent(i);
    simulation->SetBinContent(x, sum / dataPhotons->Integral() / expectedPairs->GetYaxis()->GetBinWidth(x));
    simulation->SetBinError(x, 0);  /* the simulated statistical error enters the result through the response toys */
  }
  photonjet_science::attachTrace(simulation,photonjet_science::derivedTrace("simulation_reweighted_to_data_photon_bins",
      {expectedPairs,simulatedPhotons,dataPhotons},{{"denominator",dataPhotons->Integral()},{"density_axis",1},{"errors","not calculated"}}));
  photonjet_science::setTag(simulation,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"simulation uncertainty is provided through response toys"));
  const auto [dataMean, dataError] = meanXj(data, r.xjLow, r.xjHigh);
  const double simulationMean = meanXj(simulation, r.xjLow, r.xjHigh).first;
  auto lines = xjHeader(f, r);
  if (simulationMean > 0)
    lines.push_back(Form("#LTx_{J#gamma}#GT data / sim. = %.3f #pm %.3f", dataMean / simulationMean,
                         dataError / simulationMean));
  TGraphErrors *measured = points(data, r.xjLow, r.xjHigh, kBlack, 20);
  TCanvas cv("cv", "", kCanvasW, kRatioCanvasH);
  if(simulationMean>0)photonjet_science::annotation(cv,"mean_data_over_simulation",
      photonjet_science::scalar(dataMean/simulationMean,"1","bin-centre density mean; DATA-only propagated error",{{"error",photonjet_science::real(dataError/simulationMean)}}),
      photonjet_science::derivedTrace("mean_xj_ratio",{data,simulation},{{"low",r.xjLow},{"high",r.xjHigh}}));
  auto [top, bottom] = xjRatioPads(cv);
  top->cd();
  xjFrame({{nullptr, "Corrected reconstructed", ""},
           {measured, "Data", pointKey(measured)},
           {outline(simulation, r.xjLow, r.xjHigh, gVermilion, 1), "Simulation", "l", "L"}},
          r.xjHigh, kPerPhotonTitle, lines, true);
  headroom(8);
  xjAxes(gFrame, 0, kYTitlePx);
  bottom->cd();
  TGraphErrors *ratio = ratioPoints(data, simulation, r.xjLow, r.xjHigh, kBlack, 20);
  /* the top of the ratio axis: 2, or the next step above the largest ratio with its bar (half steps up to 3, whole
     steps beyond) */
  double ratioTop = 2;
  for (int i = 0; i < ratio->GetN(); ++i) {
    const double v = ratio->GetY()[i], e = ratio->GetEY()[i];
    if (std::isfinite(v)) ratioTop = std::max(ratioTop, std::ceil(2 * (v + (std::isfinite(e) ? e : 0))) / 2);
  }
  const double step = ratioTop > 3 ? 1 : 0.5;
  ratioTop = step * std::ceil(ratioTop / step - 1e-9);
  TH1F *lower = ratioFrame(gXjShownLow, xjEnd(r.xjHigh), 0, ratioTop, kXjTitle, "Data / sim.");
  xjTicks(lower, r.xjHigh);
  xjAxes(lower, kXTitlePx, kYTitlePx, true, kRatioTitlePx);
  lowerTitles(lower, kXjTitlePx);
  /* a labelled tick at every step, each with two secondary divisions, and no label at the top edge */
  lower->GetYaxis()->SetNdivisions(static_cast<int>(std::round(ratioTop / step)) + 2 * 100, false);
  lower->GetYaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, " ");
  frameReferenceLine(1);
  photonjet_science::draw(ratio,"P");
  save(cv, name);
}

/* The photon+jet simulation of a result as the IAN names it (Photon+Jet 5+10+20 Sim): the slices of the samples that
   the saved metadata of its package lists in the photon+jet family; empty when it lists none. The slice of a sample
   is the number after "photonjet" in its name. Reads the metadata only: it does not bind the figure to the package. */
std::string simulationLine(const Result &r) {
  const PlotJson *package = packageRecord(r.package, false);
  if (!package || !package->contains("samples")) return "";
  const size_t stem = std::strlen("photonjet");
  std::set<int> slices;
  for (const auto &[sample, record] : package->at("samples").items()) {
    const auto at = sample.find("photonjet");
    if (record.value("family", "") == "photonjet" && at != std::string::npos &&
        std::atoi(sample.c_str() + at + stem) > 0)
      slices.insert(std::atoi(sample.c_str() + at + stem));
  }
  std::string line;
  for (int slice : slices) line += (line.empty() ? "Photon+Jet " : "+") + std::to_string(slice);
  return line.empty() ? line : line + " Sim";
}

/* A (photon E_T, xJ) histogram on fine bins summed onto the bins of `like`, with its errors in quadrature; the fine
   bins outside the range of `like` are left out. Null when an edge of `like` is not an edge of the fine histogram. */
TH2D *summedOnto(const TH2 *fine, const TH2 *like) {
  const auto x = edgesOf(like->GetXaxis()), y = edgesOf(like->GetYaxis());
  const auto shared = [](const TAxis *axis, const std::vector<double> &edges) {
    for (double edge : edges) {
      bool found = false;
      for (int b = 1; b <= axis->GetNbins() + 1; ++b) found = found || std::fabs(axis->GetBinLowEdge(b) - edge) < 1e-6;
      if (!found) return false;
    }
    return true;
  };
  if (!shared(fine->GetXaxis(), x) || !shared(fine->GetYaxis(), y)) return nullptr;
  auto *out = new TH2D(fresh().c_str(), "", x.size() - 1, x.data(), y.size() - 1, y.data());
  out->SetDirectory(nullptr);
  std::vector<double> variance((x.size() + 1) * (y.size() + 1), 0.0);
  for (int i = 1; i <= fine->GetNbinsX(); ++i)
    for (int j = 1; j <= fine->GetNbinsY(); ++j) {
      const int a = out->GetXaxis()->FindFixBin(fine->GetXaxis()->GetBinCenter(i));
      const int b = out->GetYaxis()->FindFixBin(fine->GetYaxis()->GetBinCenter(j));
      if (a < 1 || a > out->GetNbinsX() || b < 1 || b > out->GetNbinsY()) continue;
      out->SetBinContent(a, b, out->GetBinContent(a, b) + fine->GetBinContent(i, j));
      variance[a * (y.size() + 1) + b] += fine->GetBinError(i, j) * fine->GetBinError(i, j);
    }
  for (int a = 1; a <= out->GetNbinsX(); ++a)
    for (int b = 1; b <= out->GetNbinsY(); ++b) out->SetBinError(a, b, std::sqrt(variance[a * (y.size() + 1) + b]));
  photonjet_science::attachTrace(out,photonjet_science::projectionOf(fine,
      PlotJson::array({{{"axis",0},{"low",x.front()},{"high",x.back()}},{{"axis",1},{"low",y.front()},{"high",y.back()}}}))
      .then("sum_whole_cells_onto_edges",{{"x",x},{"y",y},{"variance","sum_diagonal"}}));
  photonjet_science::setTag(out,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"sum diagonal cell variances"));
  return out;
}

/* Reco-to-truth association, as in the IAN: per photon interval, the reconstructed xJ of the DATA region-A pairs
   before the purity correction, of every simulated reconstructed pair, and of the simulated pairs matched to the truth
   pair, each normalised to unit area over the xJ range shown. The tagged fraction is the matched share of the
   simulated pairs there. A pair is every recoil jet of the photon, or its leading recoil jet alone, as in the IAN's
   original, when the QA file holds that for DATA and for the simulation (qa/<package>/association/<class>/
   leading_data, leading_reco and leading_matched_reco on the packages' fine bins, summed here onto the bins of the
   association histograms; qa: the QA file, or null).
   With the association histograms of the file (association_pairs_A, association_pairs_reco_simulation and
   association_pairs_reco_matched_simulation: the IAN's six photon intervals, xJ steps of about 0.1) the figure is the
   IAN's grid of two rows of three frames over xJ from 0 to 2: the samples and the selection spread over the first
   row, the boxed legend in the first frame of the second. A file without them gives one row of the result's own photon
   bins over its reported xJ range (pairs_A, pairs_reco_simulation, pairs_reco_matched_simulation), with the selection
   in the first frame and the legend in the second. A row shares its y axis; the axis titles stand at the ends of
   their axes. The entry point draws the figure for p+p only, under the IAN's file name. Computed here: the unit-area
   densities and the tagged fraction, a ratio of two saved sums. */
void recoTruthAssociation(TFile &f, const Result &r, const std::string &name, TFile *qa) {
  auto *pairsA = dynamic_cast<TH2 *>(get(f, r, "association_pairs_A"));
  auto *all = dynamic_cast<TH2 *>(get(f, r, "association_pairs_reco_simulation"));
  auto *tagged = dynamic_cast<TH2 *>(get(f, r, "association_pairs_reco_matched_simulation"));
  const bool grid = pairsA && all && tagged;
  bool leadingOnly = false, globalLeading = false;
  if (grid && qa) {
    const std::string folder = "qa/" + r.package + "/association/" + r.centrality + "/";
    auto *definition = dynamic_cast<TObjString *>(qa->Get((folder + "leading_definition").c_str()));
    globalLeading = definition && std::string(definition->GetString().Data()) == "global_leading_v1";
    std::array<TH2 *, 3> fine{}, leading{};
    const char *saved[] = {"leading_data", "leading_reco", "leading_matched_reco"};
    for (int k = 0; k < 3; ++k)
      fine[k] = dynamic_cast<TH2 *>(qa->Get((folder + saved[k]).c_str()));
    if (globalLeading) {
      /* The saved global-leading xJ bins have width .05. The measurement's
         .41 and 1.78 edges cannot be formed from them. This diagnostic uses
         only whole saved bins; the measurement itself is unchanged. */
      const auto photonEdges = edgesOf(pairsA->GetXaxis());
      const std::vector<double> xjEdges{0, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.25, 1.5, 1.75, 2, 3};
      TH2D aligned(fresh().c_str(), "", photonEdges.size() - 1, photonEdges.data(), xjEdges.size() - 1, xjEdges.data());
      aligned.SetDirectory(nullptr);
      for (int k = 0; k < 3; ++k) {
        if (!fine[k] || !(leading[k] = summedOnto(fine[k], &aligned)))
          throw std::runtime_error("global-leading association needs all three QA maps on compatible saved bins");
      }
    } else {
      for (int k = 0; k < 3; ++k)
        if (fine[k]) leading[k] = summedOnto(fine[k], pairsA);
    }
    if (leading[0] && leading[1] && leading[2]) {
      pairsA = leading[0], all = leading[1], tagged = leading[2];
      leadingOnly = true;
      renderSelection(r.package, true);  /* the figure now draws from the QA package too */
    }
  }
  if (!grid) {
    pairsA = dynamic_cast<TH2 *>(get(f, r, "pairs_A"));
    all = dynamic_cast<TH2 *>(get(f, r, "pairs_reco_simulation"));
    tagged = dynamic_cast<TH2 *>(get(f, r, "pairs_reco_matched_simulation"));
  }
  if (!pairsA || !all || !tagged)
    return missing(name, "needs the association pairs, or pairs_A, pairs_reco_simulation and pairs_reco_matched_simulation");
  const int nPt = pairsA->GetNbinsX(), nX = pairsA->GetNbinsY();
  const auto photon = grid ? edgesOf(pairsA->GetXaxis()) : photonEdges(f, r);
  if (static_cast<int>(photon.size()) != nPt + 1)
    return missing(name, "the photon bins of the result are not recorded");
  const auto xj = edgesOf(pairsA->GetYaxis());
  /* the xJ range shown and normalised over: 0 to 2 in the grid, as in the IAN; the reported range of the result
     otherwise */
  const double shownLow = grid ? 0 : r.xjLow, shownHigh = grid ? 2 : r.xjHigh, end = grid ? 2 : xjEnd(r.xjHigh);
  auto shown = [&](int x) { return binInside(xj[x], xj[x + 1], shownLow, shownHigh); };
  /* row i of a (photon interval, xJ) histogram as a unit-area density over the bins shown; returns its sum there */
  auto density = [&](const std::function<double(int)> &content, const std::function<double(int)> &error, double *sum) {
    auto *h = new TH1D(fresh().c_str(), "", nX, xj.data());
    h->SetDirectory(nullptr);
    double total = 0;
    for (int x = 0; x < nX; ++x) if (shown(x)) total += content(x);
    for (int x = 0; x < nX; ++x) {
      const double width = xj[x + 1] - xj[x];
      h->SetBinContent(x + 1, total > 0 ? content(x) / total / width : 0.0);
      h->SetBinError(x + 1, total > 0 ? error(x) / total / width : 0.0);
    }
    if (sum) *sum = total;
    return h;
  };
  /* the three series of every photon interval and its tagged fraction; smaller markers in the frames of the grid. A
     bin without entries (below the kinematic reach of the interval, or beyond the tail) gets no marker. */
  const int magenta = TColor::GetColor("#D81B60");
  const double markerSize = grid ? 1.0 : 1.2;
  const auto populated = [](TGraphErrors *g) {
    for (int k = g->GetN() - 1; k >= 0; --k)
      if (g->GetY()[k] == 0 && g->GetEY()[k] == 0) g->RemovePoint(k);
    photonjet_science::copyBinding(g,g); // marker thinning changes no underlying measured bins
    return g;
  };
  std::vector<photonjet_science::Trace> fractionTrace(nPt);
  std::vector<PlotJson> fractionValue(nPt);
  std::vector<std::array<TGraphErrors *, 3>> series(nPt);
  std::vector<std::string> interval(nPt), fraction(nPt);
  for (int i = 0; i < nPt; ++i) {
    double matchedSum = 0, allSum = 0;
    TH1D *data = density([&](int x) { return pairsA->GetBinContent(i + 1, x + 1); },
                         [&](int x) { return pairsA->GetBinError(i + 1, x + 1); }, nullptr);
    TH1D *simulated = density([&](int x) { return all->GetBinContent(i + 1, x + 1); }, [](int) { return 0.0; },
                              &allSum);
    TH1D *matched = density([&](int x) { return tagged->GetBinContent(i + 1, x + 1); }, [](int) { return 0.0; },
                            &matchedSum);
    const auto restrictions=PlotJson::array({{{"axis",0},{"low",photon[i]},{"high",photon[i+1]}},{{"axis",1},{"low",shownLow},{"high",shownHigh}}});
    for(const auto& pair:std::vector<std::pair<TH1*,TH2*>>{{data,pairsA},{simulated,all},{matched,tagged}}) {
      photonjet_science::nameSeries(pair.first,"association/photon_bin/"+std::to_string(i));
      photonjet_science::attachTrace(pair.first,photonjet_science::projectionOf(pair.second,restrictions).then("unit_area_density",{{"range",{shownLow,shownHigh}},{"divide_bin_width",true}}));
      photonjet_science::setTag(pair.first,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"native density uses individual bin errors"));
    }
    fractionTrace[i]=photonjet_science::Trace::combine("tagged_over_all",{photonjet_science::projectionOf(tagged,restrictions),photonjet_science::projectionOf(all,restrictions)});
    fractionValue[i]=allSum>0?photonjet_science::scalar(matchedSum/allSum,"fraction","sum tagged / sum all in shown bins"):
      PlotJson{{"state","undefined"},{"reason","zero all-pair denominator"}};
    series[i] = {populated(points(data, shownLow, shownHigh, kBlack, 20, markerSize)),
                 populated(points(simulated, shownLow, shownHigh, gBlue, 24, markerSize)),
                 populated(points(matched, shownLow, shownHigh, magenta, 25, markerSize))};
    interval[i] = photonInterval(photon[i], photon[i + 1]);
    fraction[i] = allSum > 0 ? Form("Tagged fraction: %.2f", matchedSum / allSum) : "Tagged fraction: none";
  }

  /* Rows of frames, each row sharing its y axis: equal frames a small gap apart, the y labels on the first frame of a
     row and the x labels under the last row. Every frame is in its own transparent pad over the whole canvas, so the
     last x label of a frame is not cut at its neighbour. */
  const int columns = grid ? 3 : nPt, rows = (nPt + columns - 1) / columns;
  TCanvas cv("cv", "", 1400, rows > 1 ? 900 : 660);
  const double canvasW = cv.GetWw(), canvasH = cv.GetWh();
  const double leftPx = 96, rightPx = 22, topPx = 26, bottomPx = 82, gapPx = 14;
  const double frameW = (canvasW - leftPx - rightPx - (columns - 1) * gapPx) / columns;
  const double frameH = (canvasH - topPx - bottomPx - (rows - 1) * gapPx) / rows;
  std::vector<TPad *> pads;
  for (int i = 0; i < nPt; ++i) {
    const double frameX = leftPx + (i % columns) * (frameW + gapPx), frameY = topPx + (i / columns) * (frameH + gapPx);
    pads.push_back(seeThroughPad(cv, frameX / canvasW, 1 - (frameX + frameW) / canvasW,
                                 1 - (frameY + frameH) / canvasH, frameY / canvasH));
  }

  /* The text of the figure, as the IAN places it. First frame: the sPHENIX mark with the DATA and simulation samples.
     In the grid the photon and the jet selection stand at the top of the second and third frame and the boxed legend
     at the top left of the first frame of the second row; in a single row the whole selection stands in the first
     frame and the legend in the second. Every frame carries its photon interval with the tagged fraction in grey
     under it, right-aligned below the height of the legend. */
  const std::string simulationSample = simulationLine(r);
  const std::string dataSample = simulationSample.empty() ? systemLine(r)
                                 : r.label + " Data" + (r.energyLine.empty() ? "" : ", " + r.energyLine);
  /* the recoil angle with a space in its index, as this figure sets it */
  const std::string recoilAngle = "#Delta#phi_{#gamma, jet} " + recoilCut(r);
  const auto &selection = renderSelection(r.package);
  std::vector<std::vector<std::string>> lines(nPt);
  if (grid && nPt >= 3 && !selection.is_null()) {
    const auto &jet = selection.at("jet");
    lines[0] = {dataSample, simulationSample};
    lines[1] = {kPhotonId, Form("%s < %g", absolute("#eta^{#gamma}").c_str(),
                                selection.at("photon").at("abs_eta_max").get<double>())};
    lines[2] = {jetRadiusText(jet) + Form(", %s < %g", absolute("#eta^{jet}").c_str(), jet.at("abs_eta_max").get<double>()),
                jetThresholdText(jet, r.jetThreshold) + ", " + recoilAngle};
  } else {
    lines[0] = {dataSample, simulationSample, kPhotonId, jetLine(r),
                rewritten(recoilLine(r), "_{#gamma,jet}", "_{#gamma, jet}")};
  }
  const std::vector<std::string> labels = globalLeading
      ? std::vector<std::string>{"Data, global-leading jets", "Sim, global-leading jets",
                                  "Sim, truth-matched global-leading jets"}
      : leadingOnly ? std::vector<std::string>{"Data, leading recoil jets", "Sim, leading recoil jets",
                                                "Sim, truth-matched leading recoil jets"}
                    : std::vector<std::string>{"Data, before purity correction", "Sim, all pairs",
                                                "Sim, #gamma+jet truth-tagged"};
  const int keyed = grid && rows > 1 ? columns : nPt > 1 ? 1 : 0;
  /* the sizes of the photon interval and of the tagged fraction, the height of the legend box, the distance of the y
     title from its axis and the length of the ticks in these narrow frames, in pixels */
  const double intervalPx = grid ? 23 : 26, fractionPx = grid ? 18.5 : 20, legendPx = labels.size() * kStepPx + 10;
  const double yTitlePx = 58, tickPx = 10;
  struct TextPlace {
    double headerX, headerY, legendX1, legendX2, legendY1, legendY2, textX, intervalY, fractionY;
    std::vector<std::array<double, 3>> blocks;  /* left, right and lower edge of each block of text (NDC) */
  };
  const auto textPlace = [&](int i) {
    pads[i]->cd();
    const double w = padWidthPx(), h = padHeightPx(), top = 1 - gPad->GetTopMargin();
    TextPlace p;
    p.headerX = gPad->GetLeftMargin() + kInsetPx / w, p.headerY = top - (kHeaderTopPx + kMarkPx) / h;
    double under = top - 14 / h;
    if (i == 0) {
      const double bottom = header(p.headerX, p.headerY, lines[0], false) - 16 / h;
      p.blocks.push_back({p.headerX - 6 / w, p.headerX + headerWidth(lines[0]) + 6 / w, bottom});
      under = bottom;
    } else if (!lines[i].empty()) {  /* lines without the mark: the first on the mark's baseline */
      double widest = 0;
      for (const auto &l : lines[i]) widest = std::max(widest, textWidth(l, kLinePx));
      under = p.headerY - ((lines[i].size() - 1) * kStepPx + 16) / h;
      p.blocks.push_back({p.headerX - 6 / w, p.headerX + widest + 6 / w, under});
    }
    if (i == keyed) {
      double widest = 0;
      for (const auto &l : labels) widest = std::max(widest, textWidth(l, kLinePx));
      p.legendX1 = p.headerX - 8 / w, p.legendX2 = p.legendX1 + widest + 62 / w;
      p.legendY2 = under, p.legendY1 = under - legendPx / h;
      p.blocks.push_back({p.legendX1 - 6 / w, p.legendX2 + 6 / w, p.legendY1 - 14 / h});
    }
    p.textX = 1 - gPad->GetRightMargin() - 22 / w;
    p.intervalY = top - (14 + legendPx + 14 + intervalPx) / h, p.fractionY = p.intervalY - 1.5 * fractionPx / h;
    const double widest = std::max(textWidth(interval[i], intervalPx), textWidth(fraction[i], fractionPx));
    p.blocks.push_back({p.textX - widest - 8 / w, 1.0, p.fractionY - 14 / h});
    return p;
  };

  /* one y range for each row: opened until no point with its bars reaches a block of text in its pad, and a tenth
     more, the IAN's headroom */
  std::vector<TextPlace> places;
  std::vector<double> high(rows, 0.0);
  for (int i = 0; i < nPt; ++i) {
    places.push_back(textPlace(i));
    const double left = gPad->GetLeftMargin(), bottom = gPad->GetBottomMargin();
    const double plotX = 1 - left - gPad->GetRightMargin(), plotY = 1 - bottom - gPad->GetTopMargin();
    std::vector<Extent> drawn;
    for (auto *g : series[i]) extents(drawn, g);
    double &top = high[i / columns];
    for (const auto &e : drawn) {
      top = std::max(top, e.y / 0.95);
      const double xa = left + e.x0 / end * plotX, xb = left + e.x1 / end * plotX;
      for (const auto &[ka, kb, edge] : places[i].blocks)
        if (xb > ka && xa < kb) top = std::max(top, e.y / std::max((edge - bottom) / plotY, 0.1));
    }
  }
  for (double &top : high) top = top > 0 ? clearTop(1.1 * top, 6) : 1;

  /* in every pad the bars of the three series, then the open markers of the simulation and the data points on top */
  for (int i = 0; i < nPt; ++i) {
    pads[i]->cd();
    photonjet_science::annotation(*gPad,"association/tagged_fraction/"+std::to_string(i),fractionValue[i],fractionTrace[i]);
    const TextPlace &p = places[i];
    const bool firstOfRow = i % columns == 0, lastRow = i + columns >= nPt;
    TH1F *frame = photonjet_science::drawFrame(*gPad,0, 0, end, high[i / columns]);
    frame->GetXaxis()->SetTitle(lastRow ? kXjTitle : "");
    frame->GetYaxis()->SetTitle(firstOfRow ? "1/N dN/dx_{J#gamma}" : "");
    if (grid) {  /* a label every 0.5 and one decimal on both axes, as in the IAN */
      frame->GetXaxis()->SetNdivisions(504);
      frame->GetXaxis()->SetDecimals(), frame->GetYaxis()->SetDecimals();
    } else xjTicks(frame, r.xjHigh);
    xjAxes(frame, kXTitlePx, yTitlePx, false, kXjTitlePx, tickPx);
    frame->GetYaxis()->SetNdivisions(506);
    if (!firstOfRow) frame->GetYaxis()->SetLabelSize(0), frame->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    if (!lastRow) frame->GetXaxis()->SetLabelSize(0);
    const auto [data, simulated, matched] = series[i];
    photonjet_science::draw(simulated,"P"), photonjet_science::draw(matched,"P"), photonjet_science::draw(data,"P");
    whiteInside(simulated, 20, {});
    whiteInside(matched, 21, {data});
    TLatex text;
    text.SetNDC();
    text.SetTextFont(kFont);
    if (i == 0) header(p.headerX, p.headerY, lines[0]);
    else {
      text.SetTextSize(kLinePx);
      for (size_t l = 0; l < lines[i].size(); ++l)
        text.DrawLatex(p.headerX, p.headerY - l * kStepPx / padHeightPx(), lines[i][l].c_str());
    }
    text.SetTextAlign(31);
    text.SetTextSize(intervalPx);
    text.DrawLatex(p.textX, p.intervalY, interval[i].c_str());
    text.SetTextSize(fractionPx);
    text.SetTextColor(gGray);
    text.DrawLatex(p.textX, p.fractionY, fraction[i].c_str());
    if (i == keyed) {
      const double symbol = kSymbolPx / ((p.legendX2 - p.legendX1) * padWidthPx());
      TLegend *legend = boxedLegend(p.legendX1, p.legendY1, p.legendX2, p.legendY2, symbol);
      legend->AddEntry(data, labels[0].c_str(), "p");
      legend->AddEntry(simulated, labels[1].c_str(), "p");
      legend->AddEntry(matched, labels[2].c_str(), "p");
      photonjet_science::draw(legend);
      openKey(legend, simulated, 20);
      openKey(legend, matched, 21);
    }
    gPad->RedrawAxis();
  }
  save(cv, name);
}

/* ---- the unfolded results: classes on one frame, uncertainties by source, the prior, the ATLAS comparison -------- */
/* The colour, marker and marker size of a collision class, the same in every figure of the results: blue circles for
   p+p; for the Au+Au classes of the results file, from the most central on (Result::centralRank), red squares, green
   diamonds, orange triangles and purple crosses. A fifth class takes the first style again. */
std::tuple<int, int, double> classStyle(const Result &r) {
  if (r.package == "PP") return {gBlue, 20, 1.2};
  const std::tuple<int, int, double> styles[] = {{gRed, 21, 1.1}, {gGreen, 33, 1.7}, {gOrange, 22, 1.3},
                                                 {gPurple, 34, 1.4}};
  return styles[r.centralRank % 4];
}

/* Results that stand side by side in a bin are this far apart in xJ, about half a marker. */
const double kSideStep = 0.013;

/* The bins of a distribution inside [low, high] as points for a figure with several results: the marker stands `shift`
   bin widths from the bin centre and no farther than `limit` in xJ, so results that coincide stay visible while the
   markers of a wide bin stay together at its centre; the horizontal bar spans the bin. */
TGraphAsymmErrors *shiftedPoints(const TH1 *h, double low, double high, double shift, double limit, int color,
                                 int marker, double size, const char *role = __builtin_FUNCTION()) {
  auto *g = new TGraphAsymmErrors;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double a = h->GetXaxis()->GetBinLowEdge(i), b = h->GetXaxis()->GetBinUpEdge(i);
    if (!binInside(a, b, low, high)) continue;
    const int n = g->GetN();
    const double x = 0.5 * (a + b) + std::clamp(shift * (b - a), -limit, limit);
    g->SetPoint(n, x, h->GetBinContent(i));
    g->SetPointError(n, x - a, b - x, h->GetBinError(i), h->GetBinError(i));
  }
  style(g, color, marker, size);
  photonjet_science::bind(g,photonjet_science::seriesRole(h,role),photonjet_science::selectedTrace(h,low,high),photonjet_science::selectedHistogram(h,low,high),photonjet_science::covarianceOf(h));
  return g;
}

/* The look of a systematic box, the same for sPHENIX and for ATLAS: half its width, 0.11 of the bin and at most 0.012
   in xJ (a box is a little over a fifth of its bin wide and at most 0.024, the width of a marker); its opacity; the
   text of its legend row; and the grey box that stands for every result in the legend. */
double boxHalfWidth(double binWidth) { return std::min(0.11 * binWidth, 0.012); }
const double kBoxAlpha = 0.35;
const char *kSystematicLabel = "Systematic uncertainty";
TBox *systematicKey() {
  auto *key = new TBox;
  key->SetFillColorAlpha(gGray, kBoxAlpha);
  key->SetLineWidth(0);
  return key;
}

/* The systematic uncertainty of a result as a box around each of its points inside [low, high], moved from the bin
   centre like the point (shiftedPoints). The band is read from the systematics file: its upper side is
   xj_systematic_up and its lower side xj_systematic_down, which is saved as a negative shift. Null without a
   systematics file or without the band of this result. */
TGraphAsymmErrors *systematicBoxes(const std::string &key, double low, double high, double shift, double limit,
                                   int color) {
  if (!gSystematics) return nullptr;
  auto *c = dynamic_cast<TH1 *>(gSystematics->Get((key + "/xj_per_photon").c_str()));
  auto *up = dynamic_cast<TH1 *>(gSystematics->Get((key + "/xj_systematic_up").c_str()));
  auto *down = dynamic_cast<TH1 *>(gSystematics->Get((key + "/xj_systematic_down").c_str()));
  if (!c || !up || !down) return nullptr;
  auto *g = new TGraphAsymmErrors;
  for (int i = 1; i <= c->GetNbinsX(); ++i) {
    const double a = c->GetXaxis()->GetBinLowEdge(i), b = c->GetXaxis()->GetBinUpEdge(i);
    if (!binInside(a, b, low, high)) continue;
    const int n = g->GetN();
    const double width = b - a, half = boxHalfWidth(width);
    g->SetPoint(n, 0.5 * (a + b) + std::clamp(shift * width, -limit, limit), c->GetBinContent(i));
    g->SetPointError(n, half, half, -down->GetBinContent(i), up->GetBinContent(i));
  }
  g->SetFillColorAlpha(color, kBoxAlpha);
  g->SetLineColor(color);
  g->SetLineWidth(0);
  const auto trace=photonjet_science::Trace::combine("systematic_envelope",{photonjet_science::selectedTrace(c,low,high),photonjet_science::selectedTrace(up,low,high),photonjet_science::selectedTrace(down,low,high)});
  photonjet_science::bind(g,"systematic_envelope",trace,{{"central",photonjet_science::selectedHistogram(c,low,high)},{"up",photonjet_science::selectedHistogram(up,low,high)},{"down",photonjet_science::selectedHistogram(down,low,high)}},photonjet_science::covarianceOf(c));
  return g;
}

/* The photon window of several results as one header line when they share it; empty when they differ. */
std::string commonPhotonLine(TFile &f, const std::vector<const Result *> &shown) {
  std::string line;
  for (const Result *r : shown) {
    const std::string own = photonLine(f, *r);
    if (!line.empty() && own != line) return "";
    line = own;
  }
  return line;
}

/* The scale of a y axis from zero that ends half a division above its last labelled tick, for an axis that has to hold
   `need`: the step between the labelled ticks (1, 2, 5 or 10 times a power of ten, with at most 5.5 steps up to
   `need`) and the number of labelled ticks above zero. Such an axis ends at (labelled + 0.5) steps and takes
   labelled + 1 primary divisions. The figures of the results and the Delta-phi grid end their y axes so. */
std::pair<double, int> labelledTicks(double need) {
  double tick = std::pow(10, std::floor(std::log10(need / 5.5)));
  for (double m : {1.0, 2.0, 5.0, 10.0})
    if (need / (m * tick) <= 5.5) { tick *= m; break; }
  return {tick, std::max(1, int(std::ceil(need / tick - 0.5 - 1e-9)))};
}

/* A single-pad xJ frame for the figures of the results, in the IAN's layout: the header at the top left and the legend
   right-aligned inside the frame, at the top when it fits beside the header and under the height of the header's last
   line when it does not, as in the IAN's final-result figure, and both axis titles centred as there. The x axis is the
   xJ axis up to the reported `high`, with the xJ ticks. The y axis starts at zero and is opened until no drawn element
   reaches the legend or the header, lowered indices of its last line included; it ends half a division above its last
   labelled tick (labelledTicks), so figures of like content share one scale. close: for a figure that stands alone,
   the axis ends at the next fifth of a division above what the drawn elements need instead, so the highest of them
   stands close under the text. low: the lower edge of the reported range for a figure that shows that range alone,
   with 0.05 beside it on each side, instead of the axis from zero. behind: elements that only have to end inside the
   frame and may run behind the header and the legend (light systematic boxes; the caller raises the text over them
   with textOnTop). The ranges of the header lines are set with dashes. Returns the empty legend of `labels` rows; the
   caller fills and draws it. */
TLegend *resultFrame(double high, std::vector<Extent> drawn, const char *ytitle, std::vector<std::string> lines,
                     const std::vector<std::string> &labels, bool close = false, double low = 0,
                     const std::vector<Extent> &behind = {}) {
  for (auto &line : lines) line = rangeDashes(line);
  const double x0 = low > 0 ? xjShownLow(low) : gXjShownLow;
  const double x1 = low > 0 ? std::ceil((high + 0.05) * 20 - 1e-9) / 20 : xjEnd(high);
  const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin(), T = gPad->GetTopMargin();
  const double B = gPad->GetBottomMargin(), w = padWidthPx(), h = padHeightPx(), plot = 1 - T - B;
  const double hx = L + kInsetPx / w, baseline = header(L, 1 - T - (kHeaderTopPx + kMarkPx) / h, lines, false, false);
  double lw = 0, row = kStepPx;
  for (const auto &l : labels) {
    TLatex probe(0.5, 0.5, l.c_str());
    probe.SetNDC();
    probe.SetTextFont(kFont);
    probe.SetTextSize(kLinePx);
    UInt_t bw = 0, bh = 0;
    probe.GetBoundingBox(bw, bh);
    /* a row holds its label at full size: fitLegends lets a label be kLegendRowFill times as tall as its row */
    lw = std::max(lw, bw / w), row = std::max(row, std::ceil(bh / kLegendRowFill) + 1);
  }
  /* The room for the symbols of a legend beside the header, 10 px clear of the header and of the labels' end: the
     symbols may give up a third of their kSymbolPx before the legend moves down. */
  const double x2 = 1 - R - kLegendRightPx / w, hw = headerWidth(lines), room = x2 - hx - hw - lw - 20 / w;
  const bool beside = room >= 30 / w;
  const double symbol = beside ? std::min(room, kSymbolPx / w) : kSymbolPx / w, width = symbol + lw + 10 / w;
  const double xl = x2 - width;
  /* under the header the legend leaves room for the input note, which ianFrame draws under the header's last line */
  const double y2 = beside ? 1 - T - kLegendTopPx / h : baseline - (inputNote().empty() ? 14 : kNoteDropPx + 22) / h;
  const double y1 = y2 - labels.size() * row / h;
  /* the shares of the frame height left under the header (24 px under its last baseline) and under the legend */
  const double underHeader = std::max((baseline - 24 / h - B) / plot, 0.1);
  const double underLegend = std::max((y1 - 8 / h - B) / plot, 0.1);
  double need = 0;
  for (const auto &e : drawn) {
    if (!std::isfinite(e.y)) continue;
    const double a = L + (e.x0 - x0) / (x1 - x0) * (1 - L - R), b = L + (e.x1 - x0) / (x1 - x0) * (1 - L - R);
    need = std::max(need, e.y / 0.95);
    if (b > hx - 6 / w && a < hx + hw + 6 / w) need = std::max(need, e.y / underHeader);
    if (b > xl - 6 / w && a < x2 + 6 / w) need = std::max(need, e.y / underLegend);
  }
  for (const auto &e : behind)
    if (std::isfinite(e.y)) need = std::max(need, e.y / 0.95);
  if (!(need > 0)) need = 1;
  const auto [tick, labelled] = labelledTicks(need);
  const double top = close ? std::ceil(5 * need / tick - 1e-9) * tick / 5 : (labelled + 0.5) * tick;
  /* ianFrame draws the frame and the header; an extent beyond the axis end gives it the top chosen here */
  drawn.push_back({2 * x1 - x0, 2 * x1 - x0, 0.95 * top});
  delete ianFrame(x0, x1, 0, drawn, false, kXjTitle, ytitle, lines, {});
  gFrame->GetYaxis()->SetNdivisions(500 + static_cast<int>(std::ceil(top / tick - 1e-9)));
  gFrame->GetXaxis()->SetLabelOffset(0.015);  /* the zeros of the two axes stand apart at the origin */
  gFrame->GetXaxis()->CenterTitle();
  gFrame->GetYaxis()->CenterTitle();
  TLegend *legend = ianLegend(xl, y1, x2, y2, symbol / width);
  if (low > 0) gFrame->GetXaxis()->SetNdivisions(510);
  else xjTicks(gFrame, high);
  return legend;
}

/* Raises the text drawn so far in the current pad (the header) over everything else drawn so far. */
void textOnTop() {
  std::vector<TObject *> text;
  for (TObject *o : *gPad->GetListOfPrimitives())
    if (o->InheritsFrom(TLatex::Class())) text.push_back(o);
  for (TObject *o : text) o->Pop();
}

/* The lower edge of a results frame (gFrame) as the IAN's figures have it, where a point at zero ends at the axis: a
   white strip under the axis hides the half of such a marker that would hang below the frame. The caller redraws the
   axes after it. */
void endAtAxis() {
  const double below = 7 / (padHeightPx() * (1 - gPad->GetTopMargin() - gPad->GetBottomMargin()));
  auto *strip = new TBox(gFrame->GetXaxis()->GetXmin(), gFrameLow - below * (gFrameHigh - gFrameLow),
                         gFrame->GetXaxis()->GetXmax(), gFrameLow);
  strip->SetFillColor(kWhite);
  strip->SetLineWidth(0);
  photonjet_science::draw(strip);
}

/* Several results on one xJ frame, as the IAN draws its final results: each with the colour and marker of its class,
   its statistical uncertainty as bars and, when a systematics file is given, its systematic uncertainty as boxes. Two
   results stand at the bin centres, as in the IAN. Three or more would hide each other there and stand side by side;
   of three, the one that lies farthest from the other two stands in the middle, so that two results that coincide
   never stand next to each other. series: a result and its per-photon distribution, in legend order; windows: the
   legend names each result's photon window, as the IAN's p+p and Au+Au figure does. */
void resultOverlay(TFile &f, const std::vector<std::pair<const Result *, TH1 *>> &series, double low, double high,
                   const std::vector<std::string> &lines, bool windows) {
  const size_t n = series.size();
  std::vector<double> place(n, 0.0);  /* in steps from the bin centre */
  for (size_t i = 0; n > 2 && i < n; ++i) place[i] = i - 0.5 * (n - 1);
  if (n == 3) {
    std::vector<double> distance(n, 0.0);  /* the summed distance of a result to the other two over the bins shown */
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        for (int b = 1; b <= series[i].second->GetNbinsX(); ++b) {
          const TAxis *axis = series[i].second->GetXaxis();
          if (!binInside(axis->GetBinLowEdge(b), axis->GetBinUpEdge(b), low, high)) continue;
          const double other = series[j].second->GetBinContent(series[j].second->FindFixBin(axis->GetBinCenter(b)));
          distance[i] += std::fabs(series[i].second->GetBinContent(b) - other);
        }
    std::swap(place[std::max_element(distance.begin(), distance.end()) - distance.begin()], place[1]);
  }
  std::vector<Extent> drawn;
  std::vector<std::string> labels;
  std::vector<TGraphAsymmErrors *> marks, boxes;
  bool systematic = false;
  for (size_t i = 0; i < n; ++i) {
    const auto &[r, h] = series[i];
    const auto [color, marker, size] = classStyle(*r);
    /* side by side they stand kSideStep apart, as in the final result (a whole bin to one side, which the limit cuts
       to the step), and the bin is marked once, by the horizontal bar of the result at its centre */
    const double shift = place[i], limit = std::fabs(place[i]) * kSideStep;
    marks.push_back(shiftedPoints(h, low, high, shift, limit, color, marker, size));
    for (int point = 0; place[i] != 0 && point < marks.back()->GetN(); ++point) {
      marks.back()->SetPointEXlow(point, 0);
      marks.back()->SetPointEXhigh(point, 0);
    }
    boxes.push_back(systematicBoxes(r->key, low, high, shift, limit, color));
    extents(drawn, marks.back());
    if (boxes.back()) extents(drawn, boxes.back()), systematic = true;
    const auto e = photonEdges(f, *r);
    const std::string window = windows && e.size() > 1 ? Form(", %g-%g GeV photons", e.front(), e.back()) : "";
    labels.push_back(rangeDashes(r->label + window));
  }
  if (systematic) labels.push_back(kSystematicLabel);
  TLegend *leg = resultFrame(high, drawn, kPerPhotonTitle, lines, labels, true);
  for (auto *b : boxes) if (b) photonjet_science::draw(b,"2");
  for (auto *m : marks) photonjet_science::draw(m,"P");
  for (size_t i = 0; i < n; ++i) leg->AddEntry(marks[i], labels[i].c_str(), "lp");
  if (systematic) leg->AddEntry(systematicKey(), labels.back().c_str(), "f");
  photonjet_science::draw(leg);
  endAtAxis();
  gPad->RedrawAxis();
}

/* The Au+Au classes in the order of every figure that shows several of them: the most peripheral first. */
void peripheralFirst(std::vector<const Result *> &classes) {
  std::sort(classes.begin(), classes.end(),
            [](const Result *a, const Result *b) { return a->centralityLow > b->centralityLow; });
}

/* The unfolded per-photon xJ of p+p and of the most central Au+Au class (the one with the smallest lower edge) over
   the range both report, as in the IAN: statistical bars, and the systematic uncertainty as boxes when a systematics
   file is given. Reads xj_per_photon of both. */
void ppVsCentral(TFile &f, const std::vector<Result> &results, const std::string &name) {
  const Result *pp = nullptr, *central = nullptr;
  for (const auto &r : results) {
    if (r.package == "PP") pp = &r;
    else if (r.centralityLow >= 0 && (!central || r.centralityLow < central->centralityLow)) central = &r;
  }
  if (!pp || !central) return missing(name, "needs a p+p and an Au+Au result");
  auto *hp = get(f, *pp, "xj_per_photon"), *ha = get(f, *central, "xj_per_photon");
  if (!hp || !ha) return missing(name, "needs xj_per_photon of both systems");
  const double low = std::max(pp->xjLow, central->xjLow), high = std::min(pp->xjHigh, central->xjHigh);
  std::vector<std::string> lines{central->energyLine.empty() ? pp->energyLine : central->energyLine, jetLine(*pp)};
  if (jetLine(*central) != jetLine(*pp)) lines.push_back(jetLine(*central));
  lines.push_back(recoilLine(*pp));
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  resultOverlay(f, {{pp, hp}, {central, ha}}, low, high, lines, true);
  save(cv, name);
}

/* The final result: the unfolded per-photon xJ of p+p and of the 50-80% and 0-20% Au+Au classes on one plot over the
   range all three report, and the particle-level distribution of the p+p simulation as a dashed step line. Bars are
   statistical; the systematic uncertainty is drawn as boxes when a systematics file is given. The x axis shows the
   reported range alone. The three results of a bin stand 0.013 apart in xJ in legend order, about half a marker, and
   the bin is marked once, by the horizontal bar of the p+p point. The points and their bars stay clear of the header
   and the legend; the light boxes may run behind the header. The header gives the collision energy and the common
   photon window in one line. Reads xj_per_photon of the three results and xj_prior_simulation of p+p. */
void centralityOverlay(TFile &f, const std::vector<Result> &results, const std::string &name) {
  std::vector<const Result *> comparison;
  for (const auto &r : results) if (r.package == "PP") comparison.push_back(&r);
  for (const char *centrality : {"cent50_80", "cent0_20"})  /* the two classes of this figure, by their saved names */
    for (const auto &r : results)
      if (r.package == "AuAu" && r.centrality == centrality) comparison.push_back(&r);
  if (comparison.size() != 3 || comparison[0]->package != "PP" ||
      comparison[0]->jetThreshold != comparison[1]->jetThreshold ||
      comparison[0]->jetThreshold != comparison[2]->jetThreshold)
    return missing(name, "needs p+p, 0-20% and 50-80% AuAu results with one common recoil threshold");
  double low = 0, high = 1e9;
  for (const Result *r : comparison) low = std::max(low, r->xjLow), high = std::min(high, r->xjHigh);
  if (high <= low) return missing(name, "three-system comparison has no common reported xJ range");
  std::vector<std::pair<const Result *, TH1 *>> series;
  for (const Result *r : comparison)
    if (auto *h = get(f, *r, "xj_per_photon")) series.push_back({r, h});
  if (series.size() != 3) return missing(name, "needs xj_per_photon for all three comparison results");
  const Result &pp = *comparison[0], &nuclear = *comparison[1];
  TH1 *truth = get(f, pp, "xj_prior_simulation");
  if (!truth) return missing(name, "needs xj_prior_simulation of p+p");
  const std::string energy = nuclear.energyLine.empty() ? pp.energyLine : nuclear.energyLine;
  const std::string photons = commonPhotonLine(f, comparison);

  TGraph *generator = stepLine(truth, low, high, gGray, true);
  generator->SetLineStyle(2);
  const double step = kSideStep;
  std::vector<Extent> drawn, behind;
  std::vector<std::string> labels = {"PYTHIA 8 p+p"};
  std::vector<TGraphAsymmErrors *> marks, boxes;
  extents(drawn, generator);
  for (size_t i = 0; i < series.size(); ++i) {
    const auto &[r, h] = series[i];
    const auto [color, marker, size] = classStyle(*r);
    const double side = double(i) - 1;  /* a whole bin to one side, which the limit cuts to `step` */
    marks.push_back(shiftedPoints(h, low, high, side, step, color, marker, size));
    for (int point = 0; i > 0 && point < marks.back()->GetN(); ++point) {
      marks.back()->SetPointEXlow(point, 0);
      marks.back()->SetPointEXhigh(point, 0);
    }
    photonjet_science::copyBinding(marks.back(),marks.back()); // suppressing horizontal bars is cosmetic
    boxes.push_back(systematicBoxes(r->key, low, high, side, step, color));
    extents(drawn, marks.back());
    if (boxes.back()) extents(behind, boxes.back());
    const auto e = photonEdges(f, *r);
    const std::string window = photons.empty() && e.size() > 1 ? Form(", %g-%g GeV photons", e.front(), e.back()) : "";
    labels.push_back(rangeDashes(r->label + window));
  }
  if (!behind.empty()) labels.push_back(kSystematicLabel);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  TLegend *leg = resultFrame(high, drawn, kPerPhotonTitle,
                             {energy + (photons.empty() ? "" : ", " + photons), jetLine(pp), recoilLine(pp)}, labels,
                             true, low, behind);
  photonjet_science::draw(generator,"L");
  for (auto *box : boxes) if (box) photonjet_science::draw(box,"2");
  textOnTop();
  for (auto *m : marks) photonjet_science::draw(m,"P");
  leg->AddEntry(generator, labels[0].c_str(), "l");
  for (size_t i = 0; i < marks.size(); ++i) leg->AddEntry(marks[i], labels[i + 1].c_str(), i == 0 ? "lp" : "p");
  if (!behind.empty()) leg->AddEntry(systematicKey(), labels.back().c_str(), "f");
  photonjet_science::draw(leg);
  endAtAxis();
  gPad->RedrawAxis();
  save(cv, name);
}

/* The unfolded per-photon xJ of the Au+Au centrality classes on one plot, most peripheral first, over the range all of
   them report: statistical bars, and the systematic uncertainty as boxes when a systematics file is given. The header
   gives the collision energy and the common photon window in one line. */
void auauCentralities(TFile &f, const std::vector<Result> &results, const std::string &name) {
  std::vector<const Result *> classes;
  for (const auto &r : results) if (r.package == "AuAu" && r.centralityLow >= 0) classes.push_back(&r);
  peripheralFirst(classes);
  if (classes.size() < 2) return missing(name, "needs two or more Au+Au results");
  double low = 0, high = 1e9;
  for (const Result *r : classes) {
    if (r->jetThreshold != classes.front()->jetThreshold)
      return missing(name, "the Au+Au results differ in the recoil threshold");
    low = std::max(low, r->xjLow), high = std::min(high, r->xjHigh);
  }
  if (high <= low) return missing(name, "the Au+Au results have no common reported xJ range");
  std::vector<std::pair<const Result *, TH1 *>> series;
  for (const Result *r : classes)
    if (auto *h = get(f, *r, "xj_per_photon")) series.push_back({r, h});
  if (series.size() != classes.size()) return missing(name, "needs xj_per_photon of every Au+Au result");
  const Result &first = *classes.front();
  const std::string photons = commonPhotonLine(f, classes);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  resultOverlay(f, series, low, high,
                {first.energyLine + (photons.empty() ? "" : ", " + photons), jetLine(first), recoilLine(first)},
                photons.empty());
  save(cv, name);
}

/* The jet-radius dependence of the unfolded per-photon xJ of one result, as the IAN draws it: the same analysis with
   the recoil jets of each radius on one frame, with statistical bars. The nominal radius is the results file itself.
   Every other radius is a full result file of the same layout beside it, radius_r<two digits>.root, with its own
   verified metadata: the radius and the jet eta limit of a legend entry are read from the selection saved in the file
   the entry is drawn from, and a radius file that differs from the nominal in more than its jets is an error. A
   radius run has no systematic band, so none is drawn for any radius. The figure is drawn only when the file of every
   other radius is there and no two of the files hold the same radius. The radii stand side by side inside a bin, the
   smallest first; the nominal one is black. */
void radiusOverlay(TFile &f, const std::string &resultsPath, const Result &r, const std::string &name) {
  const std::vector<std::string> others = {"02", "03"};  /* the other radii, by the digits of their file names */
  const auto &nominal = renderSelection(r.package);
  if (nominal.is_null()) return missing(name, "needs the saved selection of the result");
  TH1 *own = get(f, r, "xj_per_photon");
  if (!own) return missing(name, "needs xj_per_photon");
  const std::string folder = gSystem->GetDirName(resultsPath.c_str()).Data();
  std::vector<std::string> paths;
  for (const auto &digits : others) {
    paths.push_back(folder + "/radius_r" + digits + ".root");
    if (gSystem->AccessPathName(paths.back().c_str()))
      return missing(name, "needs the results of the other jet radii beside the results file, radius_r02.root and "
                           "radius_r03.root");
  }
  /* the radius files are inputs of this figure alone: they leave the list of inputs when it is done */
  const auto forget = [&] { for (const auto &path : paths) gRenderInputs.erase(path); };
  /* a radius of the figure: its jet selection, its distribution, and its colour, marker and marker size */
  struct Radius { PlotJson jet; TH1 *h; int color, marker; double size; };
  const std::tuple<int, int, double> looks[] = {{TColor::GetColor("#D81B60"), 21, 1.1}, {gBlue, 22, 1.3}};
  std::vector<Radius> radii;
  PlotJson runs = PlotJson::array();  /* for the receipt: the file and the jet selection of every other radius */
  for (size_t i = 0; i < paths.size(); ++i) {
    std::unique_ptr<TFile> file(TFile::Open(paths[i].c_str(), "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error("IANFigures: cannot open " + paths[i]);
    const PlotJson metadata = loadRenderMetadata(*file, paths[i]);
    auto *h = dynamic_cast<TH1 *>(file->Get((r.key + "/xj_per_photon").c_str()));
    if (metadata.is_null() || !h || !metadata.at("results").contains(r.key) ||
        !metadata.at("packages").contains(r.package)) {
      forget();
      return missing(name, "a radius file lacks the result or its saved selection: " + paths[i]);
    }
    for (const char *field : {"package", "centrality", "recoil_min_gev", "reported_xj", "photon_edges_gev"})
      if (metadata.at("results").at(r.key).at(field) != gResultMetadata.at("results").at(r.key).at(field))
        throw std::runtime_error("a radius result differs from the nominal in more than its jets: " + paths[i]);
    h = static_cast<TH1 *>(h->Clone(fresh().c_str()));
    h->SetDirectory(nullptr);
    photonjet_science::attachTrace(h,photonjet_science::projectionOf(h,PlotJson::array()).inRole("radius/comparison_slot/"+std::to_string(i)));
    const PlotJson jet = metadata.at("packages").at(r.package).at("selection").at("jet");
    const auto [color, marker, size] = looks[i % 2];
    radii.push_back({jet, h, color, marker, size});
    runs.push_back({{"file", paths[i]}, {"jet", jet}});
  }
  own=static_cast<TH1*>(own->Clone(fresh().c_str()));own->SetDirectory(nullptr);
  photonjet_science::attachTrace(own,photonjet_science::projectionOf(own,PlotJson::array()).inRole("radius/nominal"));
  radii.push_back({nominal.at("jet"), own, int(kBlack), 20, 1.2});
  std::sort(radii.begin(), radii.end(), [](const Radius &a, const Radius &b) {
    return a.jet.at("radius").get<double>() < b.jet.at("radius").get<double>();
  });
  /* a comparison of radii needs a different radius in every file */
  for (size_t i = 1; i < radii.size(); ++i)
    if (radii[i].jet.at("radius").get<double>() - radii[i - 1].jet.at("radius").get<double>() < 1e-6) {
      forget();
      return missing(name, "the radius files and the results file do not hold three different jet radii");
    }
  const double low = r.xjLow, high = r.xjHigh;
  const size_t n = radii.size();
  std::vector<Extent> drawn;
  std::vector<std::string> labels;
  std::vector<TGraphAsymmErrors *> marks;
  for (size_t i = 0; i < n; ++i) {
    /* side by side as the results of resultOverlay: a step of 0.03 in xJ, the outermost at 0.3 of a narrow bin */
    const double place = i - 0.5 * (n - 1), shift = place * 0.6 / (n - 1), limit = std::fabs(place) * 0.03;
    marks.push_back(shiftedPoints(radii[i].h, low, high, shift, limit, radii[i].color, radii[i].marker, radii[i].size));
    extents(drawn, marks.back());
    labels.push_back(Form("R = %g, %s < %g", radii[i].jet.at("radius").get<double>(), absolute("#eta^{jet}").c_str(),
                          radii[i].jet.at("abs_eta_max").get<double>()));
  }
  /* the header names what the radii share: the system, the photons, and the jet threshold and recoil angle */
  const auto &jet = nominal.at("jet");
  const std::string photons = photonLine(f, r);
  const std::string photonEta =
      Form("%s < %g", absolute("#eta^{#gamma}").c_str(), nominal.at("photon").at("abs_eta_max").get<double>());
  const std::vector<std::string> lines = {
      systemLine(r), (photons.empty() ? std::string(kPhotonId) : photons) + ", " + photonEta,
      "anti-k_{T}, " + jetThresholdText(jet, r.jetThreshold) + ", #Delta#phi_{#gamma,jet} " + recoilCut(r)};
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  TLegend *leg = resultFrame(high, drawn, kPerPhotonTitle, lines, labels, true);
  for (auto *m : marks) photonjet_science::draw(m,"P");
  for (size_t i = 0; i < n; ++i) leg->AddEntry(marks[i], labels[i].c_str(), "lp");
  photonjet_science::draw(leg);
  endAtAxis();
  gPad->RedrawAxis();
  save(cv, name);
  forget();
#ifdef PHOTONJET_NATIVE_CAPTURE_TEST
  if(gNativeCaptureSink)return;
#endif
  const std::string receiptPath = gFig + "/" + name + ".json";
  PlotJson receipt = photonjet_render::read(receiptPath);
  receipt["radius_runs"] = runs;
  photonjet_render::write(receiptPath, receipt);
}

/* A relative uncertainty as a step line over the bins inside [low, high]. The line is interrupted over a bin where its
   value lies above the bin's ceiling (per bin number; infinite in a bin that is on scale): there the figure marks it
   with an arrow. Returns the pieces of the line; the first carries the line's look for a legend. */
std::vector<TGraph *> openSteps(const TH1 *h, double low, double high, const std::vector<double> &ceiling, int color,
                                int lineStyle, int lineWidth) {
  std::vector<TGraph *> pieces;
  const auto piece = [&]() {
    pieces.push_back(new TGraph);
    pieces.back()->SetLineColor(color);
    pieces.back()->SetLineStyle(lineStyle);
    pieces.back()->SetLineWidth(lineWidth);
  };
  piece();
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double a = h->GetXaxis()->GetBinLowEdge(i), b = h->GetXaxis()->GetBinUpEdge(i), v = h->GetBinContent(i);
    if (!binInside(a, b, low, high)) continue;
    if (v > ceiling[i]) {
      if (pieces.back()->GetN() > 0) piece();
      continue;
    }
    pieces.back()->SetPoint(pieces.back()->GetN(), a, v);
    pieces.back()->SetPoint(pieces.back()->GetN(), b, v);
  }
  for(auto* g:pieces)photonjet_science::bind(g,"relative_uncertainty",photonjet_science::selectedTrace(h,low,high),
      photonjet_science::selectedHistogram(h,low,high),photonjet_science::covarianceOf(h));
  return pieces;
}

/* The systematic uncertainty of a result by source, relative to its per-photon value in percent: per source the largest
   |shift| of its variations in each reported bin (the shift_<source>_<variation> histograms of the systematics file),
   the total systematic uncertainty and the statistical uncertainty of the result. Every reported bin contributes to
   the y range. Large relative uncertainties remain actual curves on the same axis, never unexplained floating
   overflow symbols or silently suppressed bins. Reads xj_per_photon, xj_systematic_up, xj_systematic_down and the shift histograms of the
   result in the systematics file. Computed here, not read: the total as the larger side of the band, the envelope of
   every source, and all of them relative to the central value. */
void systematicSources(TFile &f, const Result &r, const std::string &name) {
  if (!gSystematics) return missing(name, "needs a systematics file");
  TDirectory *dir = gSystematics->GetDirectory(r.key.c_str());
  auto *central = dir ? dynamic_cast<TH1 *>(dir->Get("xj_per_photon")) : nullptr;
  auto *bandUp = dir ? dynamic_cast<TH1 *>(dir->Get("xj_systematic_up")) : nullptr;
  auto *bandDown = dir ? dynamic_cast<TH1 *>(dir->Get("xj_systematic_down")) : nullptr;
  if (!central || !bandUp || !bandDown) return missing(name, "the systematics file has no " + r.key);
  /* the total: the larger side of the band drawn on the results */
  auto *total = static_cast<TH1D *>(bandUp->Clone(fresh().c_str()));
  total->SetDirectory(nullptr);
  for (int b = 1; b <= total->GetNbinsX(); ++b)
    total->SetBinContent(b, std::max(std::fabs(bandUp->GetBinContent(b)), std::fabs(bandDown->GetBinContent(b))));
  photonjet_science::attachTrace(total,photonjet_science::derivedTrace("max_abs_systematic_sides",{bandUp,bandDown}));
  auto relative = [&](const TH1 *h, bool errors) {  /* 100 |h| / central per reported bin */
    auto *out = static_cast<TH1D *>(central->Clone(fresh().c_str()));
    out->SetDirectory(nullptr);
    out->Reset();
    for (int b = 1; b <= out->GetNbinsX(); ++b) {
      const double value = errors ? h->GetBinError(b) : h->GetBinContent(b);
      if (central->GetBinContent(b) > 0) out->SetBinContent(b, 100 * std::fabs(value) / central->GetBinContent(b));
    }
    photonjet_science::attachTrace(out,photonjet_science::derivedTrace("relative_uncertainty_percent",{h,central},{{"numerator",errors?"bin_error":"bin_content"},{"zero_denominator","zero"}}));
    photonjet_science::setTag(out,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"no propagated covariance of uncertainty envelope"));
    return out;
  };
  /* the envelope of every source: the largest relative |shift| over its histograms shift_<source>_<variation> */
  const size_t stem = std::strlen("shift_");
  std::vector<std::string> order;
  std::map<std::string, TH1D *> bySource;
  TIter next(dir->GetListOfKeys());
  while (auto *key = static_cast<TKey *>(next())) {
    const std::string n = key->GetName();
    if (n.rfind("shift_", 0) != 0 || n.rfind('_') <= stem) continue;
    auto *shift = dynamic_cast<TH1 *>(key->ReadObj());
    if (!shift) continue;
    const std::string source = n.substr(stem, n.rfind('_') - stem);
    TH1D *one = relative(shift, false);
    if (!bySource.count(source)) { bySource[source] = one; order.push_back(source); continue; }
    TH1D *largest = bySource[source];
    const auto envelopeTrace=photonjet_science::derivedTrace("max_source_variation_envelope",{largest,one});
    for (int b = 1; b <= one->GetNbinsX(); ++b)
      largest->SetBinContent(b, std::max(largest->GetBinContent(b), one->GetBinContent(b)));
    photonjet_science::attachTrace(largest,envelopeTrace);
  }
  if (order.empty()) return missing(name, "the systematics file has no source shifts for " + r.key);
  /* The name, colour and line of each systematic source, the same in every result; related sources share a colour and
     differ in the line. A source is named by the category a reader of a jet measurement knows (photon purity,
     combinatorial jets, unfolding, physics model, jet energy scale), short enough for the legend to stand beside the
     header. A source that the systematics step adds is added to this table: one that is not listed is
     drawn grey under its saved key, in which TLatex reads an underscore as the start of an index. */
  const int pink = TColor::GetColor("#E7298A"), teal = TColor::GetColor("#17A2B8"), brown = TColor::GetColor("#8C564B");
  struct Curve { TH1D *h; std::string text; int color, lineStyle, lineWidth; };
  const std::map<std::string, Curve> look = {
      {"background_correlation", {nullptr, "Photon purity", gBlue, 1, 2}},
      {"combinatoric_sideband", {nullptr, "Combinatorial jets", gVermilion, 1, 2}},
      {"combinatoric_shape", {nullptr, "Combinatorial jet shape", gVermilion, 7, 2}},
      {"iterations", {nullptr, "Unfolding regularization", gPurple, 1, 2}},
      {"prior", {nullptr, "Physics model (prior)", pink, 7, 2}},
      {"prior_power", {nullptr, "Physics model (prior)", pink, 7, 2}},  /* runs before patch 38 */
      {"prior_form", {nullptr, "Physics model (form)", pink, 1, 2}},    /* runs before patch 38 */
      {"purity_nonclosure", {nullptr, "Photon purity closure", teal, 1, 2}},
      {"jet_energy_scale", {nullptr, "Jet energy scale", brown, 1, 2}},
      {"jet_energy_resolution", {nullptr, "Jet energy resolution", brown, 7, 2}},
      {"electromagnetic_scale", {nullptr, "EM energy scale", teal, 1, 2}},
      {"electromagnetic_resolution", {nullptr, "Photon energy resolution", teal, 7, 2}},
      {"photon_id_tight_pp", {nullptr, "Photon tight ID", gGreen, 1, 2}},
      {"photon_id_nontight_pp", {nullptr, "Photon non-tight ID", gGreen, 7, 2}},
      {"photon_id_tight_auau", {nullptr, "Photon tight ID", gGreen, 1, 2}},
      {"photon_id_nontightlower_auau", {nullptr, "Non-tight lower edge", gGreen, 7, 2}},
      {"photon_id_nontightupper_auau", {nullptr, "Non-tight upper edge", gGreen, 3, 2}},
      {"isolation_selection_pp", {nullptr, "Photon isolation", gOrange, 1, 2}},
      {"isolation_selection_auau", {nullptr, "Photon isolation (UE-sub.)", gOrange, 1, 2}},
      {"isolation_gap_pp", {nullptr, "Non-isolation gap", gOrange, 7, 2}},
      {"timing_selection", {nullptr, "Photon and jet timing", pink, 3, 2}}};
  std::vector<Curve> curves = {{relative(total, false), "Total systematic", kBlack, 1, 3},
                               {relative(central, true), "Statistical", gGray, 2, 2}};
  for (const auto &source : order) {
    curves.push_back(look.count(source) ? look.at(source) : Curve{nullptr, source, gGray, 3, 2});
    curves.back().h = bySource[source];
  }
  /* Every reported bin is visible, including a large relative uncertainty in a small-yield tail. */
  const TAxis *axis = central->GetXaxis();
  std::vector<int> reported;
  for (int b = 1; b <= axis->GetNbins(); ++b)
    if (binInside(axis->GetBinLowEdge(b), axis->GetBinUpEdge(b), r.xjLow, r.xjHigh)) reported.push_back(b);
  std::vector<Extent> drawn;
  std::vector<std::string> labels;
  for (const auto &c : curves) {
    labels.push_back(c.text);
    for (int b : reported)
      drawn.push_back({axis->GetBinLowEdge(b), axis->GetBinUpEdge(b), c.h->GetBinContent(b)});
  }
  auto lines = xjLines(f, r);
  for (auto &line : lines) line = rangeDashes(line);
  const double end = xjEnd(r.xjHigh);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  TLegend *leg = resultFrame(r.xjHigh, drawn, "Relative uncertainty [%]", lines, labels);
  const std::vector<double> ceiling(axis->GetNbins() + 2, std::numeric_limits<double>::infinity());
  std::vector<TGraph *> keys(curves.size());
  for (size_t i = curves.size(); i-- > 0;) {  /* the sources first, the total on top */
    const Curve &c = curves[i];
    const auto pieces = openSteps(c.h, r.xjLow, r.xjHigh, ceiling, c.color, c.lineStyle, c.lineWidth);
    for (auto *g : pieces) if (g->GetN() > 0) photonjet_science::draw(g,"L");
    keys[i] = pieces.front();
  }
  for (size_t i = 0; i < curves.size(); ++i) leg->AddEntry(keys[i], labels[i].c_str(), "l");
  photonjet_science::draw(leg);
  gPad->RedrawAxis();
  save(cv, name);
}

/* The bins of a distribution inside [low, high] as points with vertical bars only, `shift` bin widths from the bin
   centre: two results drawn over the steps of their priors stand side by side in every bin. */
TGraphErrors *besidePoints(const TH1 *h, double low, double high, double shift, int color, int marker) {
  TGraphErrors *g = points(h, low, high, color, marker);
  for (int i = 0; i < g->GetN(); ++i) {
    g->GetX()[i] += 2 * shift * g->GetEX()[i];
    g->GetEX()[i] = 0;
  }
  photonjet_science::copyBinding(g,g);
  return g;
}

/* The joint unfolding's prior: the simulated truth, the truth reweighted to the data and the half-reweighted truth of the
   prior systematic as step lines, and the unfolded data from the first two as points in the colour of their prior. The
   steps show the bins, so the points carry no bin bars, and the two results of a bin stand side by side. Reads
   xj_prior_simulation, xj_prior_reweighted, prior_weights, xj_per_photon (unfolded from the reweighted prior) and
   xj_per_photon_simulation_prior. Computed here: the half-reweighted prior, 1 + 0.5 (w - 1) of the simulated truth in
   every xJ bin (the weights depend on xJ only, so it is the per-photon simulated truth times that factor). */
void priorReweighting(TFile &f, const Result &r, const std::string &name) {
  auto *simulation = get(f, r, "xj_prior_simulation");
  auto *reweighted = get(f, r, "xj_prior_reweighted");
  auto *nominal = get(f, r, "xj_per_photon");
  auto *unfoldedSimulationPrior = get(f, r, "xj_per_photon_simulation_prior");
  auto *weights = get(f, r, "prior_weights");
  if (!simulation || !reweighted || !nominal || !unfoldedSimulationPrior || !weights)
    return missing(name, "needs the prior histograms");
  auto *half = static_cast<TH1 *>(simulation->Clone(fresh().c_str()));
  half->SetDirectory(nullptr);
  for (int i = 1; i <= half->GetNbinsX(); ++i)
    half->SetBinContent(i, simulation->GetBinContent(i) * (1 + 0.5 * (weights->GetBinContent(i) - 1)));
  photonjet_science::attachTrace(half,photonjet_science::derivedTrace("half_prior_reweight",{simulation,weights},{{"interpolation_fraction",.5}}));
  TGraphErrors *fromSimulated = besidePoints(unfoldedSimulationPrior, r.xjLow, r.xjHigh, -0.2, gVermilion, 24);
  TGraphErrors *fromReweighted = besidePoints(nominal, r.xjLow, r.xjHigh, 0.2, gBlue, 20);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  xjPad();
  xjFrame({{outline(simulation, r.xjLow, r.xjHigh, gVermilion, 2), "Simulation prior", "l", "L"},
           {fromSimulated, "Unfolded, simulation prior", pointKey(fromSimulated)},
           {outline(reweighted, r.xjLow, r.xjHigh, gBlue, 1), "Reweighted prior", "l", "L"},
           {outline(half, r.xjLow, r.xjHigh, gPurple, 7), "Half-reweighted prior (syst.)", "l", "L"},
           {fromReweighted, "Unfolded, reweighted prior", pointKey(fromReweighted)}},
          r.xjHigh, kPerPhotonTitle, xjHeader(f, r));
  headroom(8);
  xjAxes(gFrame, kXTitlePx, kYTitlePx);
  whiteInside(fromSimulated, 20, {fromReweighted});
  save(cv, name);
}

/* The per-photon xJ unfolded with k = 1 to 5 iterations of the joint unfolding, the photon yield at its own selected
   iterations (the saved scan), as step lines; the result's own k as points with its statistical bars; and every k over
   the selected one in the ratio pad, on one scale in every collision system (0.5 to 1.5, more only when a ratio is
   larger). Reads scan_xj_per_photon (k by xJ), xj_per_photon and the pair_iterations setting. */
void iterationOverlay(TFile &f, const Result &r, const std::string &name) {
  auto *scan = dynamic_cast<TH2 *>(get(f, r, "scan_xj_per_photon"));
  auto *nominal = get(f, r, "xj_per_photon");
  if (!scan || !nominal) return missing(name, "needs scan_xj_per_photon and xj_per_photon");
  const int selected = std::atoi(settingsOf(f, r.key)["pair_iterations"].c_str());
  if (selected < 1 || selected > scan->GetNbinsX()) return missing(name, "the selected iterations lie outside the saved scan");
  const int colours[] = {gGray, gBlue, gVermilion, gGreen, gPurple};
  std::vector<Series> series;
  std::vector<TH1 *> rows;
  for (int k = 1; k <= std::min(5, scan->GetNbinsX()); ++k) {
    auto *row = scan->ProjectionY(fresh().c_str(), k, k);
    row->SetDirectory(nullptr);
    photonjet_science::attachTrace(row,photonjet_science::projectionOf(scan,photonjet_science::window(0,scan->GetXaxis()->GetBinLowEdge(k),scan->GetXaxis()->GetBinUpEdge(k))).then("iteration_scan_row",{{"iteration",k}}));
    photonjet_science::setTag(row,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"iteration scan has diagonal errors only"));
    photonjet_science::nameSeries(row,"iteration_scan/"+std::to_string(k-1));
    rows.push_back(row);
    if (k != selected) series.push_back({outline(row, r.xjLow, r.xjHigh, colours[k - 1], 1), Form("k = %d", k), "l", "L"});
  }
  TGraphErrors *chosen = points(nominal, r.xjLow, r.xjHigh, kBlack, 20);
  series.push_back({chosen, Form("k = %d (selected)", selected), pointKey(chosen)});
  TCanvas cv("cv", "", kCanvasW, kRatioCanvasH);
  photonjet_science::annotation(cv,"iterations/selected",photonjet_science::scalar(selected,"iteration","native selected joint iterations"),photonjet_science::savedScalarTrace(f,r.key+"/physics_settings/pair_iterations"));
  auto [top, bottom] = xjRatioPads(cv);
  top->cd();
  xjFrame(series, r.xjHigh, kPerPhotonTitle, xjHeader(f, r), true);
  headroom(8);
  xjAxes(gFrame, 0, kYTitlePx);
  bottom->cd();
  double low = 0.5, high = 1.5;
  std::vector<TGraph *> ratios;
  for (size_t i = 0; i < rows.size(); ++i) {
    const int k = i + 1;
    if (k == selected) continue;
    auto *ratio = static_cast<TH1 *>(rows[i]->Clone(fresh().c_str()));
    ratio->SetDirectory(nullptr);
    for (int b = 1; b <= ratio->GetNbinsX(); ++b) {
      const double base = nominal->GetBinContent(b);
      ratio->SetBinContent(b, base > 0 ? rows[i]->GetBinContent(b) / base : 0), ratio->SetBinError(b, 0);
      const double a = ratio->GetXaxis()->GetBinLowEdge(b), c = ratio->GetXaxis()->GetBinUpEdge(b);
      if (base > 0 && binInside(a, c, r.xjLow, r.xjHigh))
        low = std::min(low, 0.9 * ratio->GetBinContent(b)), high = std::max(high, 1.1 * ratio->GetBinContent(b));
    }
    photonjet_science::attachTrace(ratio,photonjet_science::derivedTrace("iteration_over_nominal",{rows[i],nominal},{{"zero_denominator","zero"},{"errors","not computed"}}));
    photonjet_science::nameSeries(ratio,"iteration_ratio/"+std::to_string(i));
    ratios.push_back(outline(ratio, r.xjLow, r.xjHigh, colours[i], 1));
  }
  TH1F *lower = ratioFrame(gXjShownLow, xjEnd(r.xjHigh), low, high, kXjTitle, Form("k / k = %d", selected));
  xjTicks(lower, r.xjHigh);
  xjAxes(lower, kXTitlePx, kYTitlePx, true, kRatioTitlePx);
  lowerTitles(lower, kXjTitlePx);
  frameReferenceLine(1);
  for (auto *g : ratios) photonjet_science::draw(g,"L");
  save(cv, name);
}

/* The variations of every systematic source of a result as xJ overlays, one figure per source: the nominal per-photon
   xJ as points with its statistical bars, each variation (nominal plus its signed shift) as a step line, and each
   variation over the nominal in the ratio pad, on one scale (0.5 to 1.5, more only when a ratio is larger). The figure
   of source s is written to <stem><s>.png. Reads xj_per_photon of the results and of the systematics file and its
   shift_<source>_<variation> histograms. Computed here: nominal plus shift, and the ratios. */
void systematicVariations(TFile &f, const Result &r, const std::string &stem) {
  if (!gSystematics) return missing(stem + "*.png", "needs a systematics file");
  TDirectory *dir = gSystematics->GetDirectory(r.key.c_str());
  auto *central = dir ? dynamic_cast<TH1 *>(dir->Get("xj_per_photon")) : nullptr;
  auto *nominal = get(f, r, "xj_per_photon");
  if (!central || !nominal) return missing(stem + "*.png", "the systematics file has no " + r.key);
  /* names: the four first-pass sources and THE-482's SOURCE_LEGENDS_V1.json (e2_v7_delivery), PPG18-style */
  const std::map<std::string, std::string> sourceName = {
      {"prior", "Physics model (prior)"},
      {"prior_power", "Physics model (prior)"},
      {"iterations", "Unfolding regularization"},
      {"combinatoric_sideband", "Combinatorial jets"},
      {"background_correlation", "Photon purity"},
      {"purity_nonclosure", "Photon purity closure"},
      {"jet_energy_scale", "Jet energy scale"},
      {"jet_energy_resolution", "Jet energy resolution"},
      {"electromagnetic_scale", "Electromagnetic energy scale"},
      {"electromagnetic_resolution", "Photon energy resolution"},
      {"photon_id_tight_pp", "Photon tight ID (p+p)"},
      {"photon_id_nontight_pp", "Photon non-tight ID (p+p)"},
      {"photon_id_tight_auau", "Photon tight ID (Au+Au)"},
      {"photon_id_nontightlower_auau", "Photon non-tight lower edge (Au+Au)"},
      {"photon_id_nontightupper_auau", "Photon non-tight upper edge (Au+Au)"},
      {"isolation_selection_pp", "Photon isolation (p+p)"},
      {"isolation_selection_auau", "UE-subtracted photon isolation (Au+Au)"},
      {"isolation_gap_pp", "Photon non-isolation gap (p+p)"},
      {"timing_selection", "Photon and jet timing"}};
  const std::map<std::string, std::string> variationName = {
      {"prior/alternative", "Half of the reweighting"},
      {"iterations/minus", "k #minus 1"},
      {"iterations/plus", "k + 1"},
      {"combinatoric_sideband/low", "Sideband 0.8#minus2.06 rad"},
      {"combinatoric_sideband/mid", "Sideband 0.98#minus2.26 rad"},
      {"background_correlation/alternative", "Simulated R_{bkg}"},
      {"purity_nonclosure/alternative", "Non-closure scale"},
      {"jet_energy_scale/up", "Jet energy scale up"},
      {"jet_energy_scale/down", "Jet energy scale down"},
      {"jet_energy_resolution/low", "Less additional jet smearing"},
      {"jet_energy_resolution/high", "More additional jet smearing"},
      {"electromagnetic_scale/up", "EM energy scale +1.48%"},
      {"electromagnetic_scale/down", "EM energy scale #minus1.48%"},
      {"electromagnetic_resolution/low", "No additional photon smearing"},
      {"electromagnetic_resolution/high", "Wider photon energy resolution"},
      {"photon_id_tight_pp/alternative", "Tight ID intercept +0.05"},
      {"photon_id_nontight_pp/alternative", "Non#minustight lower intercept #minus0.10"},
      {"photon_id_tight_auau/eff60", "Tight ID: 60% efficiency"},
      {"photon_id_tight_auau/eff80", "Tight ID: 80% efficiency"},
      {"photon_id_nontightlower_auau/alternative", "Non#minustight lower score 0.20"},
      {"photon_id_nontightupper_auau/alternative", "Non#minustight upper edge: 90% efficiency"},
      {"isolation_selection_pp/eff70", "Isolation: 70% efficiency"},
      {"isolation_selection_pp/eff90", "Isolation: 90% efficiency"},
      {"isolation_selection_auau/eff70", "UE#minussubtracted isolation: 70% efficiency"},
      {"isolation_selection_auau/eff90", "UE#minussubtracted isolation: 90% efficiency"},
      {"isolation_gap_pp/gap01", "Non#minusisolation gap: 0.1 GeV"},
      {"isolation_gap_pp/gap10", "Non#minusisolation gap: 1.0 GeV"},
      {"timing_selection/alternative", "Joint photon and jet timing selection"}};
  /* the shifts of each source, in the order of the file */
  std::vector<std::string> order;
  std::map<std::string, std::vector<std::pair<std::string, TH1 *>>> shifts;
  TIter next(dir->GetListOfKeys());
  while (auto *key = static_cast<TKey *>(next())) {
    const std::string n = key->GetName();
    const size_t cut = n.rfind('_');
    if (n.rfind("shift_", 0) != 0 || cut <= 6) continue;
    auto *shift = dynamic_cast<TH1 *>(key->ReadObj());
    if (!shift) continue;
    const std::string source = n.substr(6, cut - 6), variation = n.substr(cut + 1);
    if (!shifts.count(source)) order.push_back(source);
    shifts[source].push_back({variation, shift});
  }
  if (order.empty()) return missing(stem + "*.png", "the systematics file has no source shifts for " + r.key);
  const int colours[] = {gBlue, gVermilion, gGreen, gPurple, gOrange};
  for (const auto &source : order) {
    std::vector<Series> series;
    std::vector<TH1 *> varied;
    for (size_t i = 0; i < shifts[source].size(); ++i) {
      const auto &[variation, shift] = shifts[source][i];
      auto *v = static_cast<TH1 *>(central->Clone(fresh().c_str()));
      v->SetDirectory(nullptr);
      v->Add(shift);
      photonjet_science::attachTrace(v,photonjet_science::derivedTrace("nominal_plus_systematic_shift",{central,shift}));
      varied.push_back(v);
      const std::string key = source + "/" + variation;
      series.push_back({outline(v, r.xjLow, r.xjHigh, colours[i % 5], 1), variationName.count(key) ? variationName.at(key) : variation,
                        "l", "L"});
    }
    TGraphErrors *chosen = points(nominal, r.xjLow, r.xjHigh, kBlack, 20);
    series.push_back({chosen, "Nominal", pointKey(chosen)});
    std::vector<std::string> lines = xjHeader(f, r);
    lines.push_back(sourceName.count(source) ? sourceName.at(source) : source);
    TCanvas cv("cv", "", kCanvasW, kRatioCanvasH);
    auto [top, bottom] = xjRatioPads(cv);
    top->cd();
    xjFrame(series, r.xjHigh, kPerPhotonTitle, lines, true);
    headroom(8);
    xjAxes(gFrame, 0, kYTitlePx);
    bottom->cd();
    double low = 0.5, high = 1.5;
    std::vector<TGraph *> ratios;
    for (size_t i = 0; i < varied.size(); ++i) {
      auto *ratio = static_cast<TH1 *>(varied[i]->Clone(fresh().c_str()));
      ratio->SetDirectory(nullptr);
      for (int b = 1; b <= ratio->GetNbinsX(); ++b) {
        const double base = central->GetBinContent(b);
        ratio->SetBinContent(b, base > 0 ? varied[i]->GetBinContent(b) / base : 0), ratio->SetBinError(b, 0);
        const double a = ratio->GetXaxis()->GetBinLowEdge(b), c = ratio->GetXaxis()->GetBinUpEdge(b);
        if (base > 0 && binInside(a, c, r.xjLow, r.xjHigh))
          low = std::min(low, 0.9 * ratio->GetBinContent(b)), high = std::max(high, 1.1 * ratio->GetBinContent(b));
      }
      photonjet_science::attachTrace(ratio,photonjet_science::derivedTrace("systematic_over_nominal",{varied[i],central},{{"zero_denominator","zero"},{"errors","not computed"}}));
      ratios.push_back(outline(ratio, r.xjLow, r.xjHigh, colours[i % 5], 1));
    }
    TH1F *lower = ratioFrame(gXjShownLow, xjEnd(r.xjHigh), low, high, kXjTitle, "Var. / nom.");
    xjTicks(lower, r.xjHigh);
    xjAxes(lower, kXTitlePx, kYTitlePx, true, kRatioTitlePx);
    lowerTitles(lower, kXjTitlePx);
    frameReferenceLine(1);
    for (auto *g : ratios) photonjet_science::draw(g,"L");
    save(cv, stem + source + ".png");
  }
}

/* The unfolded per-photon xJ of p+p (atlas: "pp") or of the most central Au+Au class ("pbpb_0_10") against ATLAS
   (Phys. Lett. B 789 (2019) 167, HEPData ins1694678 Table 1: p+p and 0-10% Pb+Pb at 5.02 TeV, plots/inputs/
   atlas_reference.json) and against the generator, p+p events at the result's own energy and selection. Bars are
   statistical, boxes systematic; the legend names the generator's and the ATLAS collision system and energy in their
   entries and the ATLAS selection in text rows under its entry. The open ATLAS markers are outlined over the sPHENIX
   points, so neither hides the other where they meet. An ATLAS bin that leaves the drawn axis is left out. Reads
   xj_per_photon and xj_prior_simulation, and the reference from the inputs folder beside this source file (the path the
   compiler was given for it); the receipts do not list the reference among their inputs. */
void atlasOverlay(TFile &f, const Result &r, const std::string &atlas, const std::string &name) {
  auto *h = get(f, r, "xj_per_photon");
  std::ifstream in(std::string(gSystem->GetDirName(__FILE__).Data()) + "/inputs/atlas_reference.json");
  if (!h || !in) return missing(name, "needs xj_per_photon and the ATLAS reference");
  const nlohmann::json ref = nlohmann::json::parse(in), &sel = ref.at("selection");
  const double end = xjEnd(r.xjHigh);
  auto *atlasPoints = new TGraphErrors;
  auto *atlasBoxes = new TGraphAsymmErrors;
  PlotJson referenceValues=PlotJson::array();
  for (const auto &b : ref.at("bins")) {
    const double lo = b.at("low"), hi = b.at("high"), y = b.at(atlas).at("value");
    const double stat = b.at(atlas).at("stat"), sys = b.at(atlas).at("sys");
    if (hi > end + 1e-9) continue;
    referenceValues.push_back({{"low",lo},{"high",hi},{"value",y},{"stat",stat},{"sys",sys}});
    const int n = atlasPoints->GetN();
    const double half = boxHalfWidth(hi - lo);  /* the width of the sPHENIX boxes */
    atlasPoints->SetPoint(n, 0.5 * (lo + hi), y), atlasPoints->SetPointError(n, 0.5 * (hi - lo), stat);
    atlasBoxes->SetPoint(n, 0.5 * (lo + hi), y);
    atlasBoxes->SetPointError(n, half, half, sys, sys);
  }
  PlotJson declarations=PlotJson::object();
  for(const auto& category:photonjet_science::categories())declarations[category]=photonjet_science::known("explicit published table and selection");
  PlotJson referenceSelection=sel;if(atlas=="pp")referenceSelection.erase("pbpb_centrality_percent");
  const auto referenceTrace=photonjet_science::Trace::source(photonjet_science::definition("published_xj_table",
      {{"publication","HEPData ins1694678 Table 1"},{"system",atlas}},
      PlotJson::array({photonjet_science::term("selection","published_selection",referenceSelection),
        photonjet_science::term("projection","published_values",referenceValues),
        photonjet_science::term("normalization","published_density","per-photon per-unit-xj")}),declarations));
  const auto referenceCov=photonjet_science::covariance("not_computed",{},0,"renderer consumes published stat/sys; no covariance calculated or used");
  photonjet_science::bind(atlasPoints,"reference_stat",referenceTrace,photonjet_science::portable(referenceValues),referenceCov);
  photonjet_science::bind(atlasBoxes,"reference_sys",referenceTrace,photonjet_science::portable(referenceValues),referenceCov);
  style(atlasPoints, kBlack, 25, 1.1);
  atlasBoxes->SetFillColorAlpha(gGray, kBoxAlpha);
  atlasBoxes->SetLineWidth(0);
  /* White squares under the open markers keep them open over their boxes. The ATLAS points keep ROOT's thin open
     square: whiteInside would give them the heavy outline of the sPHENIX open markers. */
  auto *white = new TGraph(atlasPoints->GetN(), atlasPoints->GetX(), atlasPoints->GetY());
  white->SetMarkerStyle(21);
  white->SetMarkerSize(1.0);
  white->SetMarkerColor(kWhite);
  photonjet_science::decoration(white,"white ATLAS marker centre overlay");
  const auto [color, marker, size] = classStyle(r);
  TGraphAsymmErrors *ours = shiftedPoints(h, r.xjLow, r.xjHigh, 0, 0, color, marker, size);
  TGraphAsymmErrors *oursBoxes = systematicBoxes(r.key, r.xjLow, r.xjHigh, 0, 0, color);
  TH1 *truth = get(f, r, "xj_prior_simulation");
  TGraph *generator = truth ? outline(truth, r.xjLow, r.xjHigh, gGray, 2) : nullptr;
  std::vector<Extent> drawn;
  extents(drawn, ours), extents(drawn, atlasPoints), extents(drawn, atlasBoxes);
  if (oursBoxes) extents(drawn, oursBoxes);
  if (generator) extents(drawn, generator);
  /* the legend: the three series with their collision systems and energies, the ATLAS selection as text rows under
     its entry, and the key of the boxes */
  const bool pp = atlas == "pp";
  const auto photons = sel.at("photon_pt_gev").get<std::vector<double>>();
  const auto centrality = sel.at("pbpb_centrality_percent").get<std::vector<double>>();
  const auto &own = renderSelection(r.package);
  const std::string simulated = own.is_null() ? "PYTHIA 8 p+p" :
      Form("PYTHIA 8 p+p, %g GeV", own.at("collision").at("sqrt_s_nn_gev").get<double>());
  /* the jet threshold after the centrality is set as an expression of its own: an index stands higher behind a line
     that starts with digits, and the row would be taller than in the p+p figure */
  const std::string jets = Form("p_{T}^{jet} #geq %g GeV", sel.at("jet_pt_min_gev").get<double>());
  const std::string central = rangeDashes(Form("%g-%g%%", centrality.at(0), centrality.at(1)));
  const double energy = sel.at("sqrt_s_nn_gev").get<double>() / 1000;
  TBox *key = systematicKey();
  std::vector<Series> rows = {{ours, "sPHENIX", "lp"}};
  if (generator) rows.push_back({generator, simulated, "l"});
  rows.push_back({atlasPoints, Form("ATLAS %s, %g TeV", pp ? "p+p" : "Pb+Pb", energy), "lp"});
  rows.push_back({nullptr, pp ? jets : central + ", " + apart(jets), ""});
  rows.push_back({nullptr, Form("%g < p_{T}^{#gamma} < %g GeV", photons.at(0), photons.at(1)), ""});
  rows.push_back({key, kSystematicLabel, "f"});
  std::vector<std::string> labels;
  for (const auto &row : rows) labels.push_back(row.text);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  TLegend *leg = resultFrame(r.xjHigh, drawn, kPerPhotonTitle, xjLines(f, r), labels);
  photonjet_science::draw(atlasBoxes,"2");
  if (oursBoxes) photonjet_science::draw(oursBoxes,"2");
  if (generator) photonjet_science::draw(generator,"L");
  photonjet_science::draw(atlasPoints,"P");
  photonjet_science::draw(white,"P");
  photonjet_science::draw(ours,"P");
  photonjet_science::draw(atlasPoints,"PX");
  for (const auto &row : rows) leg->AddEntry(row.object, row.text.c_str(), row.legend);
  photonjet_science::draw(leg);
  endAtAxis();
  gPad->RedrawAxis();
  save(cv, name);
}

/* ---- responses --------------------------------------------------------------------------------------------------- */
/* The header at the top left inside the frame of the current pad (the maps of QAFigures.C use it). */
void mapHeader(const std::vector<std::string> &lines) {
  header(gPad->GetLeftMargin() + kInsetPx / padWidthPx(),
         1 - gPad->GetTopMargin() - (kHeaderTopPx + kMarkPx) / padHeightPx(), lines);
}

/* A position on the map colour scale as a colour: ROOT's kBird (TColor palette 57, its nine stops in 255 colours), from
   dark blue at 0 through blue, teal and green to yellow at 1. The palette of the process is left as it is. */
int mapColour(double fraction) {
  static const int stops[9][3] = {{53, 42, 135}, {15, 92, 221}, {20, 129, 214}, {6, 164, 202}, {46, 183, 164},
                                  {135, 191, 119}, {209, 187, 89}, {254, 200, 50}, {249, 251, 14}};
  const double at = 8 * std::round(254 * std::min(1.0, std::max(0.0, fraction))) / 254;  /* 255 colours in all */
  const int i = std::min(7, static_cast<int>(at));
  const auto mix = [&](int c) {
    return static_cast<int>(std::lround(stops[i][c] + (at - i) * (stops[i + 1][c] - stops[i][c])));
  };
  return TColor::GetColor(mix(0), mix(1), mix(2));
}

/* The cells of a map as boxes in the IAN's map colours, on a logarithmic scale from `low` to `high`: a value below
   `low` takes the colour of `low`; a cell without entries is filled with the colour `empty`, or left white when `empty`
   is 0. */
void mapCells(const TH2 *m, double low, double high, int empty) {
  photonjet_science::observe(*gPad,m);
  for (int i = 1; i <= m->GetNbinsX(); ++i)
    for (int j = 1; j <= m->GetNbinsY(); ++j) {
      const double value = m->GetBinContent(i, j);
      if (!(value > 0) && !empty) continue;
      auto *cell = new TBox(m->GetXaxis()->GetBinLowEdge(i), m->GetYaxis()->GetBinLowEdge(j),
                            m->GetXaxis()->GetBinUpEdge(i), m->GetYaxis()->GetBinUpEdge(j));
      cell->SetFillColor(value > 0 ? mapColour(std::log(value / low) / std::log(high / low)) : empty);
      cell->SetLineWidth(0);
      photonjet_science::draw(cell);
    }
}

/* The colour bar of mapCells as the IAN sets it: a narrow bar apart from the frame and as tall as it, a tick and a
   label at every decade between `low` and `high`, and the title centred beside it. */
void mapScale(double low, double high, const char *title) {
  const double w = padWidthPx(), decades = std::log10(high / low);
  const double x1 = 1 - gPad->GetRightMargin() + 22 / w, x2 = x1 + 24 / w;
  const double y1 = gPad->GetBottomMargin(), y2 = 1 - gPad->GetTopMargin();
  const auto at = [&](double value) { return y1 + std::log10(value / low) / decades * (y2 - y1); };
  const int slices = 128;
  for (int s = 0; s < slices; ++s) {  /* each slice runs to the top and the next one covers it: no seams */
    auto *slice = new TPave(x1, y1 + s * (y2 - y1) / slices, x2, y2, 0, "NDC");
    slice->SetFillColor(mapColour((s + 0.5) / slices));
    photonjet_science::draw(slice);
  }
  TLine line;
  line.SetLineWidth(1);
  line.DrawLineNDC(x1, y1, x2, y1), line.DrawLineNDC(x1, y2, x2, y2);
  line.DrawLineNDC(x1, y1, x1, y2), line.DrawLineNDC(x2, y1, x2, y2);
  TLatex text;
  text.SetNDC();
  text.SetTextFont(kFont);
  text.SetTextSize(kLinePx);
  text.SetTextAlign(12);
  double widest = 0;
  const int first = int(std::floor(std::log10(low) + 1e-9)), last = int(std::floor(std::log10(high) + 1e-9));
  for (int d = first; d <= last; ++d)
    for (int n = 1; n <= 9; ++n) {  /* a long tick and a label at 10^d, short ticks at its multiples */
      const double value = n * std::pow(10.0, d);
      if (value < low * (1 - 1e-9) || value > high * (1 + 1e-9)) continue;
      line.DrawLineNDC(x2 - (n == 1 ? 10 : 5) / w, at(value), x2, at(value));
      if (n > 1) continue;
      const std::string label = Form("10^{%d}", d);
      text.DrawLatex(x2 + 6 / w, at(value), label.c_str());
      widest = std::max(widest, textWidth(label, kLinePx));
    }
  text.SetTextSize(kTitlePx);  /* the size of the axis titles */
  text.SetTextAngle(90);
  text.SetTextAlign(23);
  text.DrawLatex(x2 + 16 / w + widest, 0.5 * (y1 + y2), title);
}

/* The header of a response map above its frame, as in the IAN: the sPHENIX mark at the left, the collision system at
   the right. It writes no input note. */
void mapHeaderAbove(const Result &r) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(kFont);
  const double above = 1 - gPad->GetTopMargin() + 14 / padHeightPx();
  t.SetTextSize(kMarkPx);
  t.DrawLatex(gPad->GetLeftMargin(), above, sPhenixInternalLabel());
  t.SetTextSize(kLinePx);
  t.SetTextAlign(31);
  t.DrawLatex(1 - 14 / padWidthPx(), above, systemLine(r).c_str());
}

/* The margins of the canvas of a response map: a square frame, with room for the colour bar at the right. */
void mapMargins(TCanvas &cv) { cv.SetMargin(0.205, 0.19, 0.2, 0.083); }

/* Photon response as the IAN draws it: one equal cell per pair of bins, reconstructed on x and truth on y. The bins are
   the result's photon bins, with everything of the saved response grid below them as one cell and everything above as
   one cell, each labelled by the actual range of the grid it holds on its own axis; such a cell without entries is
   left out. Each truth bin is
   normalised to unit probability. The colour scale is the IAN's, 10^-3 to 1 in every panel: a smaller probability is
   drawn at the floor colour, a cell without entries is grey. Reads photon_response_grid. */
void photonResponse(TFile &f, const Result &r, const std::string &name) {
  auto *grid = dynamic_cast<TH2 *>(get(f, r, "photon_response_grid"));  /* saved as truth (x) by reco (y) */
  if (!grid) return missing(name, "needs photon_response_grid");
  const auto edges = photonEdges(f, r);
  const int bins = edges.size() - 1;
  /* the cell of a grid bin: 0 underflow, 1..bins the photon bins, bins + 1 overflow; -1 when it straddles an edge */
  const auto cellOf = [&](const TAxis *a, int i) {
    const double low = a->GetBinLowEdge(i), high = a->GetBinUpEdge(i);
    if (high <= edges.front() + 1e-9) return 0;
    if (low >= edges.back() - 1e-9) return bins + 1;
    for (int b = 0; b < bins; ++b) if (low >= edges[b] - 1e-9 && high <= edges[b + 1] + 1e-9) return b + 1;
    return -1;
  };
  /* the sums of the grid, truth cell by reconstructed cell */
  std::vector<std::vector<double>> sum(bins + 2, std::vector<double>(bins + 2, 0.0));
  for (int i = 1; i <= grid->GetNbinsX(); ++i)
    for (int j = 1; j <= grid->GetNbinsY(); ++j) {
      const int t = cellOf(grid->GetXaxis(), i), c = cellOf(grid->GetYaxis(), j);
      if (t < 0 || c < 0) return missing(name, "the saved response grid does not follow the result's photon bins");
      sum[t][c] += grid->GetBinContent(i, j);
    }
  std::vector<int> truth, reco;  /* the cells drawn on each axis */
  for (int c = 0; c <= bins + 1; ++c) {
    double inTruth = 0, inReco = 0;
    for (int o = 0; o <= bins + 1; ++o) inTruth += sum[c][o], inReco += sum[o][c];
    const bool margin = c == 0 || c == bins + 1;
    if (!margin || inTruth > 0) truth.push_back(c);
    if (!margin || inReco > 0) reco.push_back(c);
  }
  /* the range of a cell on axis a: its photon bin, or the part of the saved grid below or above the photon bins */
  const auto label = [&](int c, const TAxis *a) -> std::string {
    const double low = c == 0 ? a->GetXmin() : edges[c - 1], high = c == bins + 1 ? a->GetXmax() : edges[c];
    return Form("%g#minus%g", low, high);
  };
  auto *m = new TH2D(fresh().c_str(), "", reco.size(), 0, reco.size(), truth.size(), 0, truth.size());
  m->SetDirectory(nullptr);
  for (size_t y = 0; y < truth.size(); ++y) {
    double row = 0;
    for (int c : reco) row += sum[truth[y]][c];
    for (size_t x = 0; x < reco.size(); ++x) m->SetBinContent(x + 1, y + 1, row > 0 ? sum[truth[y]][reco[x]] / row : 0);
    m->GetYaxis()->SetBinLabel(y + 1, label(truth[y], grid->GetXaxis()).c_str());
  }
  for (size_t x = 0; x < reco.size(); ++x) m->GetXaxis()->SetBinLabel(x + 1, label(reco[x], grid->GetYaxis()).c_str());
  photonjet_science::bindMap(m,"photon_response_probability",photonjet_science::projectionOf(grid,PlotJson::array())
      .then("coarsen_transpose_and_normalize_truth_rows",{{"photon_edges",edges},{"truth_cells",truth},{"reco_cells",reco},{"zero_row","zero"}}),false);
  TCanvas cv("cv", "", kCanvasW, kMapCanvasH);
  mapMargins(cv);
  cv.SetTicks(0, 0);
  m->GetXaxis()->SetTitle("Reconstructed photon p_{T}^{#gamma} [GeV]");
  m->GetYaxis()->SetTitle("Truth photon p_{T}^{#gamma} [GeV]");
  /* the bin labels at the size of the header lines, and no axis ticks: the ticks are drawn at the cell centres below */
  for (TAxis *a : {m->GetXaxis(), m->GetYaxis()}) {
    a->SetLabelFont(kFont), a->SetLabelSize(kLinePx), a->SetTitleFont(kFont), a->SetTitleSize(kTitlePx);
    a->SetTickLength(0);
    a->CenterTitle();
  }
  m->GetXaxis()->SetLabelOffset(0.012);
  m->GetXaxis()->SetTitleOffset(2.05);
  m->GetYaxis()->SetTitleOffset(3.2);
  photonjet_science::draw(m,"AXIS");
  const double scaleFloor = 1e-3, nx = reco.size(), ny = truth.size();
  mapCells(m, scaleFloor, 1, TColor::GetColor("#E0E0E0"));
  /* the cell borders, and a 10 px tick at every cell centre on the four sides */
  TLine border;
  border.SetLineWidth(1);
  for (int x = 0; x <= nx; ++x) border.DrawLine(x, 0, x, ny);
  for (int y = 0; y <= ny; ++y) border.DrawLine(0, y, nx, y);
  const double tx = 10 / (padWidthPx() * (1 - cv.GetLeftMargin() - cv.GetRightMargin())) * nx;
  const double ty = 10 / (padHeightPx() * (1 - cv.GetTopMargin() - cv.GetBottomMargin())) * ny;
  for (int x = 0; x < nx; ++x) border.DrawLine(x + 0.5, 0, x + 0.5, ty), border.DrawLine(x + 0.5, ny - ty, x + 0.5, ny);
  for (int y = 0; y < ny; ++y) border.DrawLine(0, y + 0.5, tx, y + 0.5), border.DrawLine(nx - tx, y + 0.5, nx, y + 0.5);
  mapScale(scaleFloor, 1, "Migration probability");
  mapHeaderAbove(r);
  save(cv, name);
}

/* Joint (photon pT, xJ) response of the measurement window, drawn as the photon response beside it: reconstructed cell on
   x, truth cell on y, every truth cell normalised to unit probability over the reconstructed cells of the window, on the
   same scale (10^-3 to 1, a smaller probability at the floor colour, an empty cell white). The cells are the xJ bins of
   each photon bin in turn; lines separate the photon bins, which are named along both axes, and inside every photon bin
   x_Jgamma is marked at 0.5, 1 and 2, placed linearly inside the xJ bin that holds the value. Reads pair_response (truth
   by reco) and the xJ edges of xj_per_photon. Computed here: the probabilities and the tick positions. */
void jointResponse(TFile &f, const Result &r, const std::string &name) {
  auto *pairs = dynamic_cast<TH2 *>(get(f, r, "pair_response"));
  auto *xjAxisOf = get(f, r, "xj_per_photon");
  if (!pairs || !xjAxisOf) return missing(name, "needs pair_response and xj_per_photon");
  const auto edges = photonEdges(f, r);
  const std::vector<double> xj = edgesOf(xjAxisOf->GetXaxis());
  const int bins = edges.size() - 1, nX = xj.size() - 1, cells = pairs->GetNbinsX();
  if (bins < 1 || nX < 1 || cells != bins * nX || pairs->GetNbinsY() != cells)
    return missing(name, "the pair response is not the photon bins times the xJ bins");
  auto *m = new TH2D(fresh().c_str(), "", cells, 0, cells, cells, 0, cells);
  m->SetDirectory(nullptr);
  bool any = false;
  for (int t = 1; t <= cells; ++t) {
    double row = 0;
    for (int c = 1; c <= cells; ++c) row += pairs->GetBinContent(t, c);
    for (int c = 1; c <= cells; ++c) m->SetBinContent(c, t, row > 0 ? pairs->GetBinContent(t, c) / row : 0);
    any = any || row > 0;
  }
  if (!any) return missing(name, "the pair response has no entries");
  photonjet_science::bindMap(m,"joint_response_probability",photonjet_science::projectionOf(pairs,PlotJson::array())
      .then("transpose_normalize_truth_rows",{{"photon_edges",edges},{"xj_edges",xj},{"zero_row","zero"}}),false);
  TCanvas cv("cv", "", kCanvasW, kMapCanvasH);
  mapMargins(cv);
  cv.SetTicks(0, 0);
  m->GetXaxis()->SetTitle("Reconstructed x_{J#gamma}");
  m->GetYaxis()->SetTitle("Truth x_{J#gamma}");
  for (TAxis *a : {m->GetXaxis(), m->GetYaxis()}) {
    a->SetLabelSize(0), a->SetTickLength(0), a->SetTitleFont(kFont), a->SetTitleSize(kTitlePx);
    a->CenterTitle();
  }
  m->GetXaxis()->SetTitleOffset(2.6);
  m->GetYaxis()->SetTitleOffset(3.2);
  photonjet_science::draw(m,"AXIS");
  mapCells(m, 1e-3, 1, 0);  // Leave cells without entries white.
  const double L = cv.GetLeftMargin(), B = cv.GetBottomMargin();
  const double width = 1 - L - cv.GetRightMargin(), height = 1 - B - cv.GetTopMargin();
  const double tx = 10 / (padWidthPx() * width) * cells, ty = 10 / (padHeightPx() * height) * cells;
  TLine line;
  line.SetLineWidth(1);
  line.DrawLine(0, 0, cells, 0), line.DrawLine(0, cells, cells, cells);
  line.DrawLine(0, 0, 0, cells), line.DrawLine(cells, 0, cells, cells);
  TLatex t;
  t.SetNDC();
  t.SetTextFont(kFont);
  t.SetTextSize(kLinePx);
  /* the position of xJ value v in the cells of photon bin b: linear inside the xJ bin that holds it */
  const auto at = [&](int b, double v) {
    for (int i = 0; i < nX; ++i)
      if (v >= xj[i] && v < xj[i + 1]) return b * nX + i + (v - xj[i]) / (xj[i + 1] - xj[i]);
    return -1.0;
  };
  for (int b = 0; b < bins; ++b) {
    if (b) line.DrawLine(b * nX, 0, b * nX, cells), line.DrawLine(0, b * nX, cells, b * nX);
    for (double v : {0.5, 1.0, 2.0}) {
      const double p = at(b, v);
      if (p < 0) continue;
      line.DrawLine(p, 0, p, ty), line.DrawLine(p, cells - ty, p, cells);
      line.DrawLine(0, p, tx, p), line.DrawLine(cells - tx, p, cells, p);
      t.SetTextAngle(0), t.SetTextAlign(23);
      t.DrawLatex(L + p / cells * width, B - 10 / padHeightPx(), Form("%g", v));
      t.SetTextAlign(32);
      t.DrawLatex(L - 10 / padWidthPx(), B + p / cells * height, Form("%g", v));
    }
    const std::string range = photonInterval(edges[b], edges[b + 1]);
    t.SetTextAlign(22);
    t.SetTextAngle(0);
    t.DrawLatex(L + (b + 0.5) * width / bins, B - 56 / padHeightPx(), range.c_str());
    t.SetTextAngle(90);
    t.DrawLatex(L - 62 / padWidthPx(), B + (b + 0.5) * height / bins, range.c_str());
  }
  mapScale(1e-3, 1, "Migration probability");
  mapHeaderAbove(r);
  save(cv, name);
}

/* ---- unfolding stability and closure: their helpers, then the figures ------------------------------------------- */
/* The step between the labelled ticks of a linear axis from `low` to `high` with `divisions` primary divisions, as the
   axis painter rounds it: the width of a division up to 1, 2 or 5 times a power of ten. (clearTop holds the same rule
   for an axis from zero, with a guard of 1e-10 on the mantissa that this one does not have.) */
double labelStep(double low, double high, int divisions) {
  for (int n = std::max(divisions, 2);; ++n) {
    const double width = (high - low) / n;
    int power = static_cast<int>(std::log10(width));
    if (width <= 1) --power;
    const double mantissa = width / std::pow(10.0, power);
    const double step = (mantissa <= 1 ? 1 : mantissa <= 2 ? 2 : mantissa <= 5 ? 5 : 10) * std::pow(10.0, power);
    const int count = static_cast<int>(high / step + 1.00001) - static_cast<int>(std::floor(low / step));
    if (divisions <= 5 || 2 * count != divisions) return step;
  }
}

/* The white band the IAN leaves between its header and the data: the y range of the frame that ianFrame drew last
   (gFrame, whatever the current pad) is opened upward until an element that stood just under the header `lines` stands
   one more header line lower. The axis then ends clear above its last labelled tick (16 px above a decade, 0.4 steps on
   a linear axis), so no label stands on the frame's upper corner. The iteration, closure, yield and refolding figures
   open their frames so. */
void openUnderHeader(const std::vector<std::string> &lines) {
  const double h = padHeightPx(), T = gPad->GetTopMargin(), B = gPad->GetBottomMargin(), plot = 1 - T - B;
  const double under = (header(0.5, 1 - T - (kHeaderTopPx + kMarkPx) / h, lines, false, false) - 16 / h - B) / plot;
  const double stretch = under / (under - kStepPx / (plot * h));
  if (gPad->GetLogy()) {
    gFrameHigh = gFrameLow * std::pow(gFrameHigh / gFrameLow, stretch);
    const double top = std::log10(gFrameHigh), decade = std::floor(top + 1e-9);
    const double clear = 16 / (plot * h) * std::log10(gFrameHigh / gFrameLow);
    if (top - decade < clear) gFrameHigh = std::pow(10.0, decade + clear);
  } else {
    gFrameHigh = gFrameLow + stretch * (gFrameHigh - gFrameLow);
    for (int pass = 0; pass < 3; ++pass) {  /* a longer axis may take another step */
      const double step = labelStep(gFrameLow, gFrameHigh, gFrame->GetYaxis()->GetNdivisions() % 100);
      gFrameHigh = std::max(gFrameHigh, (std::floor(gFrameHigh / step + 1e-9) + 0.4) * step);
    }
  }
  gFrame->SetMaximum(gFrameHigh);
  gFrame->GetYaxis()->SetLimits(gFrameLow, gFrameHigh);
}

/* The labels of a logarithmic y axis as the IAN writes them: a power of ten at every decade of the frame, 10^0 and 10^1
   included, in place of the axis's own labels and at their size. */
void decadeLabels(TH1 *frame) {
  TAxis *y = frame->GetYaxis();
  const double low = std::log10(frame->GetMinimum()), high = std::log10(frame->GetMaximum());
  const double B = gPad->GetBottomMargin(), plot = 1 - gPad->GetTopMargin() - B;
  TLatex t;
  t.SetNDC();
  t.SetTextFont(y->GetLabelFont());
  t.SetTextSize(y->GetLabelSize());
  t.SetTextAlign(32);
  y->SetLabelSize(0);
  for (int d = static_cast<int>(std::ceil(low - 1e-9)); d <= static_cast<int>(std::floor(high + 1e-9)); ++d)
    t.DrawLatex(gPad->GetLeftMargin() - 7 / padWidthPx(), B + plot * (d - low) / (high - low), Form("10^{%d}", d));
}

/* The arrow of an annotation, as the IAN draws it: a thin line from (x1, y1) to (x2, y2) (NDC) with an open head. */
void openArrow(double x1, double y1, double x2, double y2) {
  const double w = padWidthPx(), h = padHeightPx(), back = std::atan2((y1 - y2) * h, (x1 - x2) * w);
  TLine stroke;
  stroke.SetLineWidth(1);
  stroke.DrawLineNDC(x1, y1, x2, y2);
  for (double side : {-0.42, 0.42})
    stroke.DrawLineNDC(x2, y2, x2 + 10 * std::cos(back + side) / w, y2 + 10 * std::sin(back + side) / h);
}

/* Open markers as the IAN draws them, white inside a thick edge and on top of their line, for a series `g` that is
   drawn with full markers: a white marker of the same shape and of size `inner` is set inside each of its points or,
   when a legend is given, inside the symbol in row `row` (from 0) of that one-column legend of the frame that ianFrame
   drew last (gFrame). The callers give `inner`: 0.75 in a circle of 1.2 is whiteSize; the larger circles (1.5 in 2.0)
   and the symbols of the linear iteration panel (0.9 in 1.4, 0.95 in 1.5) keep an edge a little heavier than that.
   (whiteInside is the other way to an open marker, for a series drawn with open markers.) */
void whiteCentres(const TGraph *g, double inner, TLegend *leg = nullptr, int row = 0) {
  auto *centres = new TGraph;
  if (leg) {
    /* the symbol stands in the middle of the legend's margin and of its row; before its first paint a legend holds its
       NDC corners as the box corners */
    const double L = gPad->GetLeftMargin(), B = gPad->GetBottomMargin();
    const double x = leg->GetX1() + 0.5 * leg->GetMargin() * (leg->GetX2() - leg->GetX1());
    const double y = leg->GetY2() - (row + 0.5) * (leg->GetY2() - leg->GetY1()) / leg->GetNRows();
    const double fx = (x - L) / (1 - L - gPad->GetRightMargin()), fy = (y - B) / (1 - B - gPad->GetTopMargin());
    const double x0 = gFrame->GetXaxis()->GetXmin(), x1 = gFrame->GetXaxis()->GetXmax();
    const double y0 = gFrame->GetMinimum(), y1 = gFrame->GetMaximum();
    centres->SetPoint(0, x0 + fx * (x1 - x0), gPad->GetLogy() ? y0 * std::pow(y1 / y0, fy) : y0 + fy * (y1 - y0));
  } else {
    for (int i = 0; i < g->GetN(); ++i) centres->SetPoint(i, g->GetX()[i], g->GetY()[i]);
  }
  centres->SetMarkerStyle(g->GetMarkerStyle());
  centres->SetMarkerSize(inner);
  centres->SetMarkerColor(kWhite);
  photonjet_science::decoration(centres,"white marker or legend centre overlay");
  photonjet_science::draw(centres,"P");
}

/* The frozen iteration rule of one unfolding, as in the IAN: the statistical term (the data and the response statistics
   of the rule table in quadrature), the change from k - 1 and their quadrature Q_k against k, with the selected k
   marked by a dashed line and named. Drawn from the saved rule table; a result that used configured iterations has none
   and is listed as not drawn. table: the name of the saved table, whose rows Q, S_data, S_response, D_k and selected
   run against k. lines: the header. Computed here from the saved table: the statistical term, for which the table has
   no row, and the name of the rule, "Minimum Q_k" when the selected k is the minimum of Q and "Selected" otherwise. */
void iterationRule(TFile &f, const Result &r, const char *table, const std::vector<std::string> &lines,
                   const std::string &name) {
  auto *h = dynamic_cast<TH2 *>(get(f, r, table));
  if (!h) return missing(name, "needs the frozen iteration rule (this result used configured iterations)");
  const auto row = [&](const char *label) { return h->GetYaxis()->FindBin(label); };
  const int q = row("Q"), data = row("S_data"), response = row("S_response"), change = row("D_k");
  const int chosen = row("selected");
  int selected = 0, k0 = 0, k1 = 0;
  for (int k = 1; k <= h->GetNbinsX(); ++k) {
    if (!(h->GetBinContent(k, q) > 0)) continue;
    if (!k0) k0 = k;
    k1 = k;
    if (h->GetBinContent(k, chosen) > 0.5) selected = k;
  }
  if (!k0) return missing(name, "the iteration rule table is empty");
  const std::function<double(int)> terms[] = {
      [&](int k) { return std::hypot(h->GetBinContent(k, data), h->GetBinContent(k, response)); },
      [&](int k) { return h->GetBinContent(k, change); }, [&](int k) { return h->GetBinContent(k, q); }};
  double low = std::numeric_limits<double>::infinity(), highest = 0;
  for (const auto &term : terms)
    for (int k = k0; k <= k1; ++k)
      if (const double v = term(k); v > 0) low = std::min(low, v), highest = std::max(highest, v);
  /* A log axis when every term stays within a factor of 300, about two and a half decades, as in the IAN's joint
     panel; otherwise a linear axis from zero, as in its photon panel, where the change from k - 1 falls to nothing. The
     symbols are the IAN's: filled circles, squares and diamonds on the log axis; on the linear axis, where Q_k lies on
     the statistical term, open circles and squares (white inside a thick edge) under filled circles. */
  const bool logy = low > highest / 300;
  const std::vector<std::string> labels = {"Data + response statistics",
                                           "Iteration change, " + apart("k vs.") + " " + apart("k") + " #minus 1",
                                           "Quadrature " + apart("Q") + "_{k}"};
  const int colors[] = {gBlue, gOrange, gGreen}, markers[] = {20, 21, logy ? 33 : 20};
  const double sizes[] = {1.4, logy ? 1.3 : 1.5, logy ? 2.2 : 1.4};
  const double centres[] = {logy ? 0 : 0.9, logy ? 0 : 0.95, 0};  /* the white inside of an open symbol */
  std::vector<TGraph *> curves;
  std::vector<Extent> drawn;
  for (int term = 0; term < 3; ++term) {
    auto *g = new TGraph;
    for (int k = k0; k <= k1; ++k)
      if (const double v = terms[term](k); std::isfinite(v) && (logy ? v > 0 : v >= 0)) g->SetPoint(g->GetN(), k, v);
    photonjet_science::bindNative(g,"iteration_rule/"+std::to_string(term),photonjet_science::projectionOf(h,PlotJson::array())
        .then("iteration_rule_terms",{{"term",term},{"rows",{q,data,response,change}},{"k_first",k0},{"k_last",k1}}),
        photonjet_science::covariance("not_computed",{},0,"no uncertainty of rule-table uncertainty terms"));
    style(g, colors[term], markers[term], sizes[term]);
    extents(drawn, g);
    curves.push_back(g);
  }
  /* the two lines that name the selected k: by the rule when it is the minimum of Q_k, as in the IAN's photon panel */
  bool minimum = selected > 0;
  for (int k = k0; k <= k1 && minimum; ++k) minimum = !(terms[2](k) < terms[2](selected));
  const std::string rule = minimum ? "Minimum " + apart("Q") + "_{k}" : "Selected";
  const std::string value = bayesIterations(selected);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  photonjet_science::annotation(cv,"selected_iteration",photonjet_science::scalar(selected,"iteration","stored selected row > 0.5"),
      photonjet_science::projectionOf(h,PlotJson::array()).then("selected_iteration_row",{{"row",chosen}}));
  cv.SetLeftMargin(0.125);  /* the IAN's wide frame: the header and the legend stand side by side in every class */
  cv.SetRightMargin(0.015);
  const double L = cv.GetLeftMargin(), R = cv.GetRightMargin(), B = cv.GetBottomMargin(), T = cv.GetTopMargin();
  const double plot = 1 - T - B, w = padWidthPx(), hpx = padHeightPx();
  const auto xNdc = [&](double k) { return L + (k - k0 + 0.5) / (k1 - k0 + 1) * (1 - L - R); };
  /* The legend box, sized as the frame sizes it. On the log axis it takes the frame's top right, as in the IAN's joint
     panel; on the linear axis the lower right, as in its photon panel: above the change from k - 1 where that has
     fallen to nothing and below the statistical term, unless the selected k is marked there. */
  double lw = 0;
  for (const auto &l : labels) lw = std::max(lw, textWidth(l, kLinePx));
  const double symbol = kSymbolPx / w, width = symbol + lw + 10 / w, height = labels.size() * kStepPx / hpx;
  const double x2 = 1 - R - kLegendRightPx / w;
  /* on the linear axis the two lines start nameGapPx beside the guide of the selected k */
  const double nameGapPx = 36;
  const double nameWidth = nameGapPx / w + std::max(textWidth(rule, kLinePx), textWidth(value, kLinePx));
  const bool lowerRight = !logy && xNdc(selected) + nameWidth + 12 / w < x2 - width;
  TLegend *leg = ianFrame(k0 - 0.5, k1 + 0.5, 0, drawn, logy, "Bayesian iterations, k", "Relative uncertainty / change",
                          lines, lowerRight ? std::vector<std::string>{} : labels);
  if (!logy) gFrame->GetYaxis()->SetDecimals(true), gFrame->GetYaxis()->SetNdivisions(506);
  openUnderHeader(lines);
  pixelAxes(gFrame);
  gFrame->GetXaxis()->SetTitleOffset(1.5);
  gFrame->GetYaxis()->SetTitleOffset(1.85);
  gFrame->GetXaxis()->CenterTitle(), gFrame->GetYaxis()->CenterTitle();
  gFrame->GetXaxis()->SetNdivisions(k1 - k0 + 1, 1, 1, true);  /* a tick at every k */
  /* short ticks, as in the IAN: a symbol on the axis covers the tick it stands on, and the name of the selected k at
     the foot of its guide stands clear of them */
  tickLengths(gFrame, 6);
  if (logy) decadeLabels(gFrame);
  const auto yNdc = [&](double y) {
    const double span = logy ? std::log(gFrameHigh / gFrameLow) : gFrameHigh - gFrameLow;
    return B + plot * (logy ? std::log(y / gFrameLow) : y - gFrameLow) / span;
  };
  if (lowerRight) {
    double under = 0, over = gFrameHigh;  /* the change and the statistical terms over the legend's columns */
    for (int k = k0; k <= k1; ++k)
      if (xNdc(k + 0.3) > x2 - width)
        under = std::max(under, terms[1](k)), over = std::min({over, terms[0](k), terms[2](k)});
    /* as fractions of the frame height */
    const double lowest = under / gFrameHigh, room = over / gFrameHigh - height / plot;
    const double y1 = B + plot * std::min(lowest + 0.12, 0.5 * (lowest + room));
    leg = ianLegend(x2 - width, y1, x2, y1 + height, symbol / width);
  }
  /* The selected k: a dashed guide from the axis, named beside it. On the log axis the guide runs to the height the
     header leaves free (the legend too, when it runs under it) and the name stands at its foot, as in the IAN's joint
     panel. On the linear axis the guide ends a little above the selected Q_k, and two lines beside an arrow to that
     point name it, as in the IAN's photon panel, between the change from k - 1 and the statistical term. */
  const double guide = xNdc(selected), point = selected ? yNdc(terms[2](selected)) : 0;
  const bool left = selected > 0.5 * (k0 + k1);  /* the text on the side with more room */
  const int step = left ? -1 : 1;
  const double side = step, pitch = (1 - L - R) / (k1 - k0 + 1);
  /* the room for the two lines and their arrow under the point: down to the change from k - 1 over their columns */
  double ground = B;
  const int columns = static_cast<int>(nameWidth / pitch) + 1;  /* the columns of k the two lines stand over */
  for (int k = selected; selected && k >= k0 && k <= k1 && std::abs(k - selected) <= columns; k += step)
    ground = std::max(ground, yNdc(terms[1](k)) + 8 / hpx);
  /* from the point down to the top of the first line */
  const double lead = std::min(62.0, (point - ground) * hpx - 65);
  const bool arrow = selected && !logy && lead >= 30;
  if (selected) {
    double block = header(0.5, 1 - T - (kHeaderTopPx + kMarkPx) / hpx, lines, false, false) - 28 / hpx;
    /* before its first paint a legend holds its NDC corners as the box corners */
    if (!lowerRight && guide > leg->GetX1() - 6 / w) block = std::min(block, leg->GetY1() - 12 / hpx);
    dashedLine(guide, B, guide, arrow ? std::min(block, point + 22 / hpx) : block, gGray, 1, true);
  }
  for (int term = 0; term < 3; ++term) {
    photonjet_science::draw(curves[term],"LP");
    if (centres[term] > 0) whiteCentres(curves[term], centres[term]);
  }
  if (selected) {
    TLatex text;
    text.SetNDC();
    text.SetTextFont(kFont);
    text.SetTextSize(kLinePx);
    if (arrow) {
      const double top = point - lead / hpx, tailX = guide + side * 30 / w, tailY = top + 6 / hpx;
      const double reach = std::hypot(30.0, (point - tailY) * hpx);  /* the arrow stops short of the symbol */
      openArrow(tailX, tailY, guide + side * 30 / w * 10 / reach, point - (point - tailY) * 10 / reach);
      text.SetTextAlign(left ? 31 : 11);
      text.DrawLatex(guide + side * nameGapPx / w, top - 16 / hpx, rule.c_str());
      text.DrawLatex(guide + side * nameGapPx / w, top - 45 / hpx, value.c_str());
    } else {
      text.SetTextColor(gGray);
      text.SetTextAlign(left ? 32 : 12);
      text.DrawLatex(guide + side * 9 / w, B + 29 / hpx, ("Selected " + value).c_str());
    }
  }
  for (int term = 0; term < 3; ++term) leg->AddEntry(curves[term], labels[term].c_str(), "lp");
  photonjet_science::draw(leg);
  for (int term = 0; term < 3; ++term)
    if (centres[term] > 0) whiteCentres(curves[term], centres[term], leg, term);
  save(cv, name);
}

/* The two iteration figures of a result: the rule of the photon unfolding, and the rule of the joint unfolding with
   the reported xJ range in its header. */
void photonIterationRule(TFile &f, const Result &r, const std::string &name) {
  iterationRule(f, r, "iteration_rule_photons", {systemLine(r), photonLine(f, r)}, name);
}
void jointIterationRule(TFile &f, const Result &r, const std::string &name) {
  iterationRule(f, r, "iteration_rule_pairs",
                {systemLine(r), photonLine(f, r), "Joint (p_{T}^{#gamma}, x_{J#gamma}) unfolding",
                 Form("%g #leq x_{J#gamma} #leq %g", r.xjLow, r.xjHigh)},
                name);
}

/* The header of a simulation-closure figure, as in the IAN: the collision system, the simulation that is unfolded (by
   the package of the result, in the wording of the IAN's figure) and the kind of closure. */
std::vector<std::string> closureLines(const Result &r, const char *kind) {
  return {systemLine(r), r.package == "PP" ? "PYTHIA8 simulation" : "Embedded simulation", kind};
}

/* The titles of the lower pad of a two-pad figure, as in the IAN: the x title centred under the frame, the y title
   centred on the frame and no longer than it, at its distance from the axis. */
void lowerTitles(TH1 *frame, double target) {
  /* target > 0 (patch 49, Justin Oct 15: "the y axis label not large enough for subpanel"): the size of the upper
     pad's y title, kept wherever the title fits; centred on the frame, the title may use the pad's top margin, never
     more (the pad clips it) */
  TAxis *y = frame->GetYaxis();
  const double h = padHeightPx();
  const double room = h * (1 - gPad->GetBottomMargin() - (target > 0 ? -1 : 1) * gPad->GetTopMargin());
  const double size = y->GetTitleSize(), want = target > 0 ? target : size;
  const double length = textWidth(y->GetTitle(), want) * padWidthPx();
  const double fitted = std::min(want, std::floor(2 * want * 0.96 * room / length) / 2);
  y->SetTitleOffset(y->GetTitleOffset() * size / fitted);
  y->SetTitleSize(fitted);
  y->CenterTitle();
  frame->GetXaxis()->CenterTitle();
}

/* The lower frame of a closure or refolding figure ends a little beyond its outer labels, so they stand clear of the
   pad corners: the labelled half range times this factor. */
const double kBeyondLabels = 1.12;

/* The range of a ratio frame around 1 for the points of `ratio`: labelled over the IAN's range, 0.6 to 1.4, and wider
   when a ratio with its error lies outside; the frame ends a little beyond its outer labels (kBeyondLabels). room: the
   half range of the labels over the largest distance of a point with its error from 1 (more for a frame that carries
   a legend row at its top). */
std::pair<double, double> rangeAroundOne(const TGraphErrors *ratio, double room = 1.15) {
  double half = 0.4;
  for (int i = 0; i < ratio->GetN(); ++i) {
    const double v = ratio->GetY()[i], e = ratio->GetEY()[i];
    if (std::isfinite(v)) half = std::max(half, room * (std::fabs(v - 1) + (std::isfinite(e) ? e : 0)));
  }
  for (double step : {0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0})
    if (half <= step) { half = step; break; }
  return {std::max(0.0, 1 - kBeyondLabels * half), 1 + kBeyondLabels * half};
}

/* The lower frame of a simulation-closure figure, unfolded / truth, for the points of `ratio`, over their range around
   1 (rangeAroundOne). A dashed line marks 1. */
TH1F *closureRatioFrame(double x0, double x1, const TGraphErrors *ratio, const char *xtitle, double room = 1.15) {
  const auto [low, high] = rangeAroundOne(ratio, room);
  TH1F *frame = ratioFrame(x0, x1, low, high, xtitle, "Unf. / truth");
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetDecimals(true);
  lowerTitles(frame, kXjTitlePx);
  frameReferenceLine(1);
  return frame;
}

/* The lower edge of a log axis for two spectra of points, as in the IAN: a little under the lowest point and its error
   bar, and under the decade below that point when the axis would otherwise label a single decade. 0 (the frame's own
   rule) when no point is positive. */
double logLowerEdge(const TGraphErrors *a, const TGraphErrors *b) {
  double low = std::numeric_limits<double>::infinity(), high = 0;
  for (const auto *g : {a, b})
    for (int i = 0; i < g->GetN(); ++i) {
      const double y = g->GetY()[i], bar = y - g->GetEY()[i];
      if (y > 0) low = std::min(low, bar > 0 ? bar : y / 2), high = std::max(high, y);
    }
  if (!std::isfinite(low)) return 0;
  /* the decades labelled between the edge and the points; the frame opens at least a factor 2 above the highest */
  const int decades = static_cast<int>(std::floor(std::log10(2 * high)) - std::ceil(std::log10(low / 1.8))) + 1;
  return decades < 2 ? std::min(low / 1.8, 0.9 * std::pow(10.0, std::floor(std::log10(low)))) : low / 1.8;
}

/* The side margins of the pads of the two photon-yield figures (closure and yield): the IAN's wide frame, in which the
   header and the legend stand side by side in every class. */
const double kYieldLeft = 0.14, kYieldRight = 0.03;

/* A closure test of the photon-yield unfolding, as in the IAN and named as in the sPHENIX isolated-photon note. Full
   closure test: the simulation that fills the response is the one unfolded, so the unfolding returns its truth
   exactly. Half closure test: one half of the simulated sources fills the response, the other half is unfolded and
   compared with its own truth; the bars of the unfolded yield and of the ratio are the spread of unfolded - truth
   over resamples of the sources. Upper frame, on a log scale: truth, the reconstructed yield that goes into the
   unfolding (reco (input)) and the unfolded yield. Lower frame: unfolded / truth. Truth is the larger open circle, so
   it stays visible around a coinciding unfolded point: its white centre is drawn before the unfolded points. Reads
   the photons_truth, photons_input, photons_unfolded and photons_sigma views of estimator_closure_full or of
   estimator_closure_A_to_B. A results file without the views of the full closure still gives that test from
   photon_truth and closure_photons, without the input. */
void photonClosure(TFile &f, const Result &r, const std::string &name, bool half) {
  const std::string saved = std::string("estimator_closure_") + (half ? "A_to_B" : "full") + "_photons_";
  TH1 *truthPhotons = get(f, r, saved + "truth"), *unfoldedPhotons = get(f, r, saved + "unfolded");
  TH1 *inputPhotons = get(f, r, saved + "input"), *spread = get(f, r, saved + "sigma");
  if (!half && (!truthPhotons || !unfoldedPhotons))
    truthPhotons = get(f, r, "photon_truth"), unfoldedPhotons = get(f, r, "closure_photons");
  if (!truthPhotons || !unfoldedPhotons || (half && !spread))
    return missing(name, half ? "needs the half closure of the photon yield (estimator_closure_A_to_B_photons)"
                              : "needs closure_photons and photon_truth");
  if (half)
    for (int b = 1; b <= unfoldedPhotons->GetNbinsX(); ++b) unfoldedPhotons->SetBinError(b, spread->GetBinContent(b));
  if(half)photonjet_science::attachTrace(unfoldedPhotons,photonjet_science::derivedTrace("closure_errors_from_resample_sigma",{unfoldedPhotons,spread}));
  const auto edges = edgesOf(truthPhotons->GetXaxis());
  const double x0 = edges.front(), x1 = edges.back();
  TCanvas cv("cv", "", kCanvasW, kRatioCanvasH);
  auto [top, bottom] = ratioPads(cv, kYieldLeft, kYieldRight);
  top->cd();
  TGraphErrors *truth = points(truthPhotons, x0, x1, gBlue, 20, 2.0);
  TGraphErrors *unfolded = points(unfoldedPhotons, x0, x1, gRed, 20);
  TGraphErrors *input = inputPhotons ? points(inputPhotons, x0, x1, gGray, 21) : nullptr;
  std::vector<Extent> drawn;
  extents(drawn, truth), extents(drawn, unfolded);
  if (input) extents(drawn, input);
  const auto lines = closureLines(r, half ? "Half closure test" : "Full closure test");
  std::vector<std::string> labels = {"Truth", "Unfolded"};
  if (input) labels.insert(labels.begin() + 1, "Reco (input)");
  TLegend *leg = ianFrame(x0, x1, logLowerEdge(input ? input : truth, unfolded), drawn, true, "", "Weighted photons",
                          lines, labels);
  openUnderHeader(lines);
  upperAxes(gFrame);
  gFrame->GetYaxis()->CenterTitle();
  decadeLabels(gFrame);
  edgeTicks(gFrame, edges, false);
  if (input) photonjet_science::draw(input,"P");
  photonjet_science::draw(truth,"P"), whiteCentres(truth, 1.5), photonjet_science::draw(unfolded,"P");
  leg->AddEntry(truth, "Truth", "lp");
  if (input) leg->AddEntry(input, "Reco (input)", "lp");
  leg->AddEntry(unfolded, "Unfolded", "lp");
  photonjet_science::draw(leg);
  whiteCentres(truth, 1.5, leg);
  bottom->cd();
  TGraphErrors *ratio = ratioPoints(unfoldedPhotons, truthPhotons, x0, x1, gRed, 20);
  edgeTicks(closureRatioFrame(x0, x1, ratio, "p_{T}^{#gamma, truth} [GeV]"), edges, true);
  photonjet_science::draw(ratio,"P");
  save(cv, name);
}

/* Purity-corrected reconstructed photon yields and the unfolded yields, and unfolded / reco, as in the IAN. The ratio
   axis is the same in every collision class of the results file. Reads photons_corrected and photons_unfolded, of this
   result for the points and of every result of the file for the range of the ratio axis. */
void photonYield(TFile &f, const Result &r, const std::string &name) {
  auto *correctedPhotons = get(f, r, "photons_corrected");
  auto *unfoldedPhotons = get(f, r, "photons_unfolded");
  if (!correctedPhotons || !unfoldedPhotons) return missing(name, "needs photons_corrected and photons_unfolded");
  const auto edges = edgesOf(correctedPhotons->GetXaxis());
  const double x0 = edges.front(), x1 = edges.back();
  TCanvas cv("cv", "", kCanvasW, kRatioCanvasH);
  auto [top, bottom] = ratioPads(cv, kYieldLeft, kYieldRight);
  top->cd();
  TGraphErrors *reco = points(correctedPhotons, x0, x1, gRed, 20);
  TGraphErrors *unfolded = points(unfoldedPhotons, x0, x1, kBlack, 20);
  std::vector<Extent> drawn;
  extents(drawn, reco), extents(drawn, unfolded);
  const std::vector<std::string> lines = {systemLine(r), bayesIterations(r.photonIterations)};
  const std::vector<std::string> labels = {"Corrected reconstructed", "Unfolded"};
  TLegend *leg = ianFrame(x0, x1, logLowerEdge(reco, unfolded), drawn, true, "", "Photon yield", lines, labels);
  openUnderHeader(lines);
  upperAxes(gFrame);
  gFrame->GetYaxis()->CenterTitle();
  decadeLabels(gFrame);
  edgeTicks(gFrame, edges, false);
  photonjet_science::draw(reco,"P"), photonjet_science::draw(unfolded,"P"), whiteCentres(unfolded, 0.75);
  leg->AddEntry(reco, labels[0].c_str(), "lp"), leg->AddEntry(unfolded, labels[1].c_str(), "lp");
  photonjet_science::draw(leg);
  whiteCentres(unfolded, 0.75, leg, 1);
  bottom->cd();
  TGraphErrors *ratio = ratioPoints(unfoldedPhotons, correctedPhotons, x0, x1, gBlue, 20);
  /* One ratio axis for every collision class of the results file, as the panels of the IAN's figure share one: the
     IAN's range, 0.75 to 2.55, opened in half steps until the ratio of every class lies inside with its error. The
     ratios of the other classes are computed here from their saved yields, read without binding the figure to their
     packages; a bin without corrected photons is left out. */
  double low = 0.75, high = 2.55;
  TIter next(f.GetListOfKeys());
  while (auto *key = static_cast<TKey *>(next())) {
    if (std::string(key->GetClassName()) != "TDirectoryFile") continue;
    const auto *above = dynamic_cast<TH1 *>(f.Get(Form("%s/photons_unfolded", key->GetName())));
    const auto *below = dynamic_cast<TH1 *>(f.Get(Form("%s/photons_corrected", key->GetName())));
    if (!above || !below || above->GetNbinsX() != below->GetNbinsX()) continue;
    for (int i = 1; i <= below->GetNbinsX(); ++i) {
      const double a = above->GetBinContent(i), b = below->GetBinContent(i);
      if (!(b > 0)) continue;
      const double v = a / b, e = std::hypot(above->GetBinError(i) * b, below->GetBinError(i) * a) / (b * b);
      low = std::min(low, std::max(0.0, std::floor(2 * (v - e - 0.2)) / 2 - 0.05));
      high = std::max(high, std::ceil(2 * (v + e + 0.2)) / 2 + 0.05);
    }
  }
  TH1F *lower = ratioFrame(x0, x1, low, high, "p_{T}^{#gamma} [GeV]", "Unfolded / reco (incl. eff.)");
  lower->GetYaxis()->SetNdivisions(505);
  lower->GetYaxis()->SetDecimals(true);
  lowerTitles(lower);
  edgeTicks(lower, edges, true);
  frameReferenceLine(1);
  photonjet_science::draw(ratio,"P"), whiteCentres(ratio, 0.75);
  save(cv, name);
}

/* A closure test of the joint unfolding, as in the IAN and named as the one of the photon yield (photonClosure): the
   per-photon xJ of the simulated truth and of the unfolded simulated reconstruction over the reported range, and
   unfolded / truth. Full closure test: the simulation that fills the response is unfolded and compared with its own
   truth, which the unfolding returns exactly. Half closure test: the simulated sources are split in two, one half
   fills the response and fits the prior, the other is unfolded and compared with its own truth; the bars are the
   spread of unfolded - truth over resamples of the sources. Truth is the larger open circle, so it stays visible
   around a coinciding unfolded point. Reads the unfolded and truth views of estimator_closure_full, or the unfolded,
   truth and sigma views of estimator_closure_A_to_B; the last row of a view is the integral over the photon bins,
   which is drawn. */
void pairClosure(TFile &f, const Result &r, const std::string &name, bool half) {
  const std::string saved = std::string("estimator_closure_") + (half ? "A_to_B" : "full") + "_";
  const auto view = [&](const char *what) { return dynamic_cast<TH2 *>(get(f, r, saved + what)); };
  TH2 *unfoldedViews = view("unfolded"), *truthViews = view("truth"), *spreadViews = half ? view("sigma") : nullptr;
  if (!unfoldedViews || !truthViews || (half && !spreadViews))
    return missing(name, "needs the closure views of the joint unfolding (" + saved + "unfolded, truth" +
                             (half ? ", sigma)" : ")"));
  /* the rows of a closure view: every photon bin of the result, then their integral, which is drawn */
  const int all = unfoldedViews->GetNbinsX();
  if (all != static_cast<int>(photonEdges(f, r).size()) || truthViews->GetNbinsX() != all)
    return missing(name, "the closure views are not the result's photon bins and their integral");
  const auto integral = [&](TH2 *views) {
    TH1D *row = views->ProjectionY(fresh().c_str(), all, all);
    row->SetDirectory(nullptr);
    for (int x = 1; x <= row->GetNbinsX(); ++x) row->SetBinError(x, 0);
    photonjet_science::attachTrace(row,photonjet_science::projectionOf(views,photonjet_science::window(0,views->GetXaxis()->GetBinLowEdge(all),views->GetXaxis()->GetBinUpEdge(all))).then("closure_integral_row_without_statistical_errors"));
    photonjet_science::setTag(row,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"closure row does not calculate covariance"));
    return row;
  };
  TH1D *truth = integral(truthViews), *unfolded = integral(unfoldedViews);
  for (int x = 1; half && x <= unfolded->GetNbinsX(); ++x) {
    const double spread = spreadViews->GetBinContent(all, x);
    unfolded->SetBinError(x, std::isfinite(spread) ? spread : 0);
  }
  if(half)photonjet_science::attachTrace(unfolded,photonjet_science::derivedTrace("closure_errors_from_resample_sigma",{unfolded,spreadViews}));
  TCanvas cv("cv", "", kCanvasW, kRatioCanvasH);
  auto [top, bottom] = ratioPads(cv);
  top->cd();
  TGraphErrors *truthPoints = points(truth, r.xjLow, r.xjHigh, gBlue, 20, 2.0);
  TGraphErrors *unfoldedPoints = points(unfolded, r.xjLow, r.xjHigh, gRed, 20);
  const auto lines = closureLines(r, half ? "Half closure test" : "Full closure test");
  TLegend *leg = xjFrame({{truthPoints, "Truth", "lp"}, {unfoldedPoints, "Unfolded", "lp"}}, r.xjHigh,
                         kPerPhotonTitle, lines, true);
  openUnderHeader(lines);
  gFrame->GetYaxis()->SetDecimals(true);
  gFrame->GetYaxis()->CenterTitle();
  /* truth white inside its edge, in the frame and in the legend, and the unfolded points on top again */
  whiteCentres(truthPoints, 1.5), photonjet_science::draw(unfoldedPoints,"P");
  whiteCentres(truthPoints, 1.5, leg);
  bottom->cd();
  TGraphErrors *ratio = ratioPoints(unfolded, truth, r.xjLow, r.xjHigh, gRed, 20);
  xjTicks(closureRatioFrame(gXjShownLow, xjEnd(r.xjHigh), ratio, kXjTitle), r.xjHigh);
  photonjet_science::draw(ratio,"P");
  save(cv, name);
}

/* Refolding check of the joint unfolding, as in the IAN: the corrected measured pairs and the unfolded result folded
   back through the response, per unit xJ in plain counts, and their pull, (refolded - measured) over the measured
   statistical uncertainty, computed here bin by bin. Reads pairs_corrected and pairs_refolded. */
void jointRefold(TFile &f, const Result &r, const std::string &name) {
  auto *measuredPairs = dynamic_cast<TH2 *>(get(f, r, "pairs_corrected"));
  auto *refoldedPairs = dynamic_cast<TH2 *>(get(f, r, "pairs_refolded"));
  if (!measuredPairs || !refoldedPairs) return missing(name, "needs pairs_corrected and pairs_refolded");
  TH1D *measured = xjProjection(measuredPairs), *refolded = xjProjection(refoldedPairs);
  TCanvas cv("cv", "", kCanvasW, kRatioCanvasH);
  auto [top, bottom] = ratioPads(cv);
  top->cd();
  const std::vector<std::string> lines = {systemLine(r), photonLine(f, r)};
  xjFrame({{steps(refolded, r.xjLow, r.xjHigh, gBlue), "Refolded", "l"},
           {points(measured, r.xjLow, r.xjHigh, kBlack, 20), "Corrected measured", "lp"}},
          r.xjHigh, kPairsTitle, lines, true);
  openUnderHeader(lines);
  gFrame->GetYaxis()->SetMaxDigits(5);  /* 5000 and 20000 as the IAN prints them, without a x10^3 exponent */
  gFrame->GetYaxis()->SetTitleOffset(2.4);  /* clear of five-digit labels */
  gFrame->GetYaxis()->CenterTitle();
  bottom->cd();
  auto *pull = new TGraphErrors;
  double half = 4;
  for (int i = 1; i <= measured->GetNbinsX(); ++i) {
    const double lo = measured->GetXaxis()->GetBinLowEdge(i), hi = measured->GetXaxis()->GetBinUpEdge(i);
    const double error = measured->GetBinError(i);
    if (!binInside(lo, hi, r.xjLow, r.xjHigh) || !(error > 0)) continue;
    const double value = (refolded->GetBinContent(i) - measured->GetBinContent(i)) / error;
    pull->SetPoint(pull->GetN(), 0.5 * (lo + hi), value);
    pull->SetPointError(pull->GetN() - 1, 0.5 * (hi - lo), 0);
    half = std::max(half, 2 * std::ceil((std::fabs(value) + 0.5) / 2));  /* even, so the ticks are whole numbers */
  }
  if (pull->GetN() == 0) return missing(name, "no measured bin with a statistical uncertainty in the reported range");
  photonjet_science::bindNative(pull,"refolding_pull",photonjet_science::derivedTrace("refold_minus_measured_over_measured_error",
      {refolded,measured},{{"low",r.xjLow},{"high",r.xjHigh},{"positive_error_required",true}}),
      photonjet_science::covariance("not_computed",{},0,"pull uncertainty not calculated"));
  style(pull, gBlue, 20);
  const double end = xjEnd(r.xjHigh);
  /* the frame ends a little beyond its outer labels (kBeyondLabels) */
  TH1F *lower = ratioFrame(gXjShownLow, end, -kBeyondLabels * half, kBeyondLabels * half, kXjTitle, "Pull");
  xjTicks(lower, r.xjHigh);
  lower->GetYaxis()->SetNdivisions(505);  /* every 2 over the IAN's range, -4 to 4 */
  lower->GetYaxis()->SetTitleOffset(2.4);  /* under the title of the upper pad */
  lowerTitles(lower, kXjTitlePx);
  TBox band;  /* two statistical standard deviations */
  band.SetFillColorAlpha(gGreen, 0.10);
  band.SetLineWidth(0);
  const double left = lower->GetXaxis()->GetXmin(), right = lower->GetXaxis()->GetXmax();  /* its own frame */
  band.DrawBox(left, -2, right, 2);
  frameReferenceLine(0, kBlack, 1, kSolid);
  photonjet_science::draw(pull,"P");
  gPad->RedrawAxis();
  save(cv, name);
}

/* ---- the Delta-phi grid ------------------------------------------------------------------------------------------ */
/* The Delta-phi figure counts in steps of pi/32 over [0, pi]. Region A (tight, isolated photons) is the first x bin of
   the saved abcd and leakage histograms. */
const int kDphiSteps = 32, kRegionA = 1;

/* A saved photon-jet Delta-phi distribution in steps of pi/32 from step `first` on. The saved axis splits once more at
   each end of the sideband of the combinatoric estimate; the two parts of such a step are added. */
TH1D *dphiSteps(const TH1 *saved, int first) {
  auto *h = new TH1D(fresh().c_str(), "", kDphiSteps - first, first * M_PI / kDphiSteps, M_PI);
  h->SetDirectory(nullptr);
  for (int i = 1; i <= saved->GetNbinsX(); ++i) {
    const int b = h->FindFixBin(saved->GetXaxis()->GetBinCenter(i));
    h->SetBinContent(b, h->GetBinContent(b) + saved->GetBinContent(i));
    h->SetBinError(b, std::hypot(h->GetBinError(b), saved->GetBinError(i)));
  }
  photonjet_science::attachTrace(h,photonjet_science::projectionOf(saved,photonjet_science::window(0,first*M_PI/kDphiSteps,M_PI))
      .then("sum_dphi_steps",{{"first",first},{"steps",kDphiSteps},{"variance","sum_diagonal"}}));
  photonjet_science::setTag(h,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"sum diagonal dphi variances"));
  return h;
}

/* A multiple of pi/32 as the IAN writes it on the Delta-phi axis (0, #pi/8, 3#pi/4, #pi). */
std::string piLabel(double x) {
  const int n = int(std::lround(x / M_PI * kDphiSteps));
  if (n <= 0) return "0";
  int d = kDphiSteps;
  while (n % d) d /= 2;
  return (n == d ? std::string() : std::to_string(n / d)) +
         (d == kDphiSteps ? "#pi" : "#pi/" + std::to_string(kDphiSteps / d));
}

/* Computed here, not read: the scale of the simulated combinatoric jets. The function gives the simulated combinatoric
   jets of a result (jets with no truth jet nearby) per bin of the saved Delta-phi axis and per DATA photon, at the
   level the result's subtraction gives them. The data-driven estimate counts the jets of the DATA in the Delta-phi
   sideband and carries them into the recoil region with the simulation: the simulated distribution of every jet is
   scaled to the DATA in the sideband of the result's settings, whose two ends are bin edges of the saved axis. The
   simulation template is counted per simulated region-A photon (row A of the saved leakage). data: the saved Delta-phi
   distribution of the DATA; photons: the region-A photons of the DATA. Reads dphi_simulation_combinatoric,
   dphi_simulation_truthjet and, for the simulation template, leakage. Null when the result subtracts nothing or does
   not hold the simulated distributions and its sideband. */
TH1 *simulatedCombinatoric(TFile &f, const Result &r, const TH1 *data, double photons) {
  TH1 *combinatoricJets = get(f, r, "dphi_simulation_combinatoric"), *truthJets = get(f, r, "dphi_simulation_truthjet");
  if (!combinatoricJets || !truthJets || combinatoricJets->GetNbinsX() != data->GetNbinsX()) return nullptr;
  double scale = 0;
  if (r.combinatoric.rfind("data_dphi_simulation", 0) == 0) {
    const TAxis *axis = data->GetXaxis();
    double measured = 0, simulated = 0;
    for (int i = 1; i <= axis->GetNbins(); ++i)
      if (axis->GetBinLowEdge(i) > r.sidebandLow - 1e-9 && axis->GetBinUpEdge(i) < r.sidebandHigh + 1e-9)
        measured += data->GetBinContent(i),
            simulated += combinatoricJets->GetBinContent(i) + truthJets->GetBinContent(i);
    if (simulated > 0) scale = measured / simulated / photons;
  } else if (r.combinatoric == "simulation") {
    auto *leakage = dynamic_cast<TH2 *>(get(f, r, "leakage"));
    double simulatedPhotons = 0;
    for (int y = 1; leakage && y <= leakage->GetNbinsY(); ++y) simulatedPhotons += leakage->GetBinContent(kRegionA, y);
    if (simulatedPhotons > 0) scale = 1 / simulatedPhotons;
  }
  if (!(scale > 0)) return nullptr;
  const auto scaledTrace=photonjet_science::derivedTrace("combinatoric_simulation_normalization",{combinatoricJets,truthJets,data},
      {{"method",r.combinatoric},{"scale",scale},{"photons",photons},{"sideband",{r.sidebandLow,r.sidebandHigh}}});
  combinatoricJets->Scale(scale);
  photonjet_science::attachTrace(combinatoricJets,scaledTrace);
  return combinatoricJets;
}

/* Reconstructed photon-jet azimuthal separation per region-A photon and unit Delta-phi before the purity and
   combinatoric corrections, as in the IAN: one grid of frames with shared axes, p+p first and then the Au+Au classes
   from peripheral to central, over the whole saved range. Every jet above the threshold around the leading tight
   isolated photon is counted, from the DATA event records. The dashed line is the recoil requirement of the saved
   selection and the tinted band behind it the region the xJ measurement uses. In the Au+Au frames the step line is
   the simulated combinatoric component at the level the subtraction of the class gives it (simulatedCombinatoric).
   Reads dphi_region_A and abcd of every result. A result without them, without selection metadata or without region-A
   photons is left out of the grid, with a line on the output; the file keeps its name whatever the number of frames.
   The frames carry no input note. */
void photonJetDphi(TFile &f, const std::vector<Result> &results, const std::string &name) {
  std::vector<const Result *> order, classes;
  for (const auto &r : results) (r.package == "PP" ? order : classes).push_back(&r);
  peripheralFirst(classes);
  order.insert(order.end(), classes.begin(), classes.end());
  /* a label: its text size and colour and its baseline in px under the frame top; a frame: its points, its simulated
     combinatoric component (Au+Au), its recoil requirement and its labels with their right and lower end in px from
     the frame's top left corner */
  struct Label { std::string text; double px; int color; double y; };
  struct Panel {
    const Result *r; TGraphErrors *data; TGraph *simulated; double cut, photons; std::vector<Label> text;
    double textRight, textBottom;
  };
  std::vector<Panel> panels;
  std::vector<std::string> selectionLines;  /* of the first frame */
  const double dphiStep = M_PI / kDphiSteps;
  const double smallPitchPx = 28;  /* the pitch of the small lines (kSmallPx) */
  const auto leftOut = [](const Result *r, const char *why) {
    printf("IANFigures: %s left out of the Delta-phi grid (%s)\n", r->key.c_str(), why);
  };
  double x0 = M_PI;
  for (const Result *r : order) {
    TH1 *saved = get(f, *r, "dphi_region_A");
    auto *abcd = dynamic_cast<TH2 *>(get(f, *r, "abcd"));
    const auto &selection = renderSelection(r->package);
    if (!saved || !abcd || selection.is_null()) {
      leftOut(r, "needs dphi_region_A, abcd and selection metadata");
      continue;
    }
    double photons = 0;
    for (int y = 1; y <= abcd->GetNbinsY(); ++y) photons += abcd->GetBinContent(kRegionA, y);
    if (!(photons > 0)) {
      leftOut(r, "no region-A photons");
      continue;
    }
    const int first = std::max(0, int(std::floor(saved->GetXaxis()->GetXmin() / dphiStep + 1e-6)));
    if (first >= kDphiSteps) {
      leftOut(r, "the saved Delta-phi axis starts at pi");
      continue;
    }
    TH1D *h = dphiSteps(saved, first);
    h->Scale(1.0 / photons, "width");
    photonjet_science::attachTrace(h,photonjet_science::derivedTrace("dphi_per_region_A_photon_and_width",{h,abcd},{{"photons",photons},{"region",kRegionA}}));
    TGraph *simulated = nullptr;
    if (TH1 *component = r->package == "PP" ? nullptr : simulatedCombinatoric(f, *r, saved, photons)) {
      TH1D *shown = dphiSteps(component, first);
      shown->Scale(1.0, "width");
      photonjet_science::attachTrace(shown,photonjet_science::projectionOf(shown,PlotJson::array()).then("divide_by_bin_width"));
      simulated = outline(shown, first * dphiStep, M_PI, gVermilion, 1);
    }
    x0 = std::min(x0, first * dphiStep);
    /* every frame: its system, photon window and photon count; the first also the sPHENIX mark and the selection */
    const auto &jet = selection.at("jet");
    std::string count = Form("%.0f", photons);  /* with a comma between the thousands: 1,278,959 */
    for (int at = int(count.size()) - 3; at > 0; at -= 3) count.insert(at, ",");
    std::vector<Label> text;
    double baselinePx = kHeaderTopPx + kMarkPx;
    if (panels.empty()) {
      text.push_back({sPhenixInternalLabel(), kMarkPx, kBlack, baselinePx});
      baselinePx += 1.15 * kStepPx;  /* the pitch of the header (header): 0.15 more under the larger mark */
    }
    text.push_back({rangeDashes(systemLine(*r)), kLinePx, kBlack, baselinePx});
    text.push_back({photonLine(f, *r), kLinePx, kBlack, baselinePx += kStepPx});
    text.push_back({"N_{#gamma} = " + count, kSmallPx, gGray, baselinePx += smallPitchPx});
    if (panels.empty()) {
      const double jetEta = jet.at("abs_eta_max").get<double>();
      std::string jets = jetThresholdText(jet, r->jetThreshold);
      if (jet.contains("overlap_delta_r"))
        jets += Form(", #DeltaR_{#gamma,jet} %s %g", cutOperator(jet.at("overlap_operator")).c_str(),
                     jet.at("overlap_delta_r").get<double>());
      const std::string photonEta =
          Form(", %s < %g", absolute("#eta^{#gamma}").c_str(), selection.at("photon").at("abs_eta_max").get<double>());
      const std::string radius = jetRadiusText(jet);
      selectionLines = {kPhotonId + photonEta, radius + Form(", %s < %g", absolute("#eta^{jet}").c_str(), jetEta), jets,
                        "Reco level, no purity correction", "No Au+Au combinatorial-jet subtraction"};
      baselinePx += 8;
      for (const auto &line : selectionLines) text.push_back({line, kSmallPx, kBlack, baselinePx += smallPitchPx});
    }
    const double cut = jet.at("delta_phi_min_rad").get<double>();
    panels.push_back({r, points(h, first * dphiStep, M_PI, kBlack, 20, 0.9), simulated, cut, photons, text, 0, baselinePx + 8});
  }
  if (panels.empty())
    return missing(name, "needs dphi_region_A and abcd of an event-record result with selection metadata");
  /* the canvas in px: its edges around the frames, the gap between two frames and the equal frames they leave */
  const int columns = panels.size() > 1 ? 2 : 1, rows = (panels.size() + columns - 1) / columns;
  TCanvas cv("cv", "", columns > 1 ? 1250 : kCanvasW, rows > 1 ? 520 * rows : kCanvasH);
  const double canvasW = cv.GetWw(), canvasH = cv.GetWh();
  const double leftPx = 100, rightPx = 22, bottomPx = 112, topPx = 22, gapPx = 20;
  const double frameW = (canvasW - leftPx - rightPx - (columns - 1) * gapPx) / columns;
  const double frameH = (canvasH - topPx - bottomPx - (rows - 1) * gapPx) / rows;
  for (auto &p : panels) {
    double width = 0;
    for (const auto &l : p.text) width = std::max(width, textWidth(l.text, l.px));
    p.textRight = kInsetPx + width * canvasW;
  }
  /* the boxed legend, once: under the labels of the second frame (of the only frame when there is one) */
  const size_t keyed = panels.size() > 1 ? 1 : 0;
  TGraph *legendSample = nullptr;  /* the simulated component that the legend shows: the first frame's that has one */
  for (const auto &p : panels) if (!legendSample) legendSample = p.simulated;
  const std::string cutEntry = "Nominal cut, " + absolute("#Delta#phi") + " " + recoilCut(*panels[keyed].r);
  std::vector<std::string> entries = {"Data", cutEntry, "x_{J#gamma} recoil region"};
  if (legendSample) entries.insert(entries.begin() + 1, "Simulated combinatorial jets");
  double entryWidth = 0;
  for (const auto &e : entries) entryWidth = std::max(entryWidth, textWidth(e, kLinePx));
  /* the box in px from the frame's top left corner: a symbol of 52 px, rows of 34 px */
  const double boxLeft = kInsetPx - 6, boxWidth = 66 + entryWidth * canvasW, boxTop = panels[keyed].textBottom + 6;
  const double boxHeight = 34.0 * entries.size() + 8;
  /* one y range per row of frames, opened until no point and no step of the simulated component reaches a label or
     the legend; it ends half a division above its last labelled tick, so the tick labels of two rows never meet */
  std::vector<double> high(rows, 0.0);
  std::vector<int> divisions(rows, 1);
  for (size_t i = 0; i < panels.size(); ++i) {
    const Panel &p = panels[i];
    const auto open = [&](double x, double y) {  /* for an element whose left end is x and whose top is y */
      const double px = (x - x0) / (M_PI - x0) * frameW;
      double room = 0.9 * frameH;
      if (px < p.textRight + 8) room = std::min(room, frameH - p.textBottom - 6);
      if (i == keyed && px < boxLeft + boxWidth + 8) room = std::min(room, frameH - boxTop - boxHeight - 8);
      high[i / columns] = std::max(high[i / columns], y * frameH / room);
    };
    for (int k = 0; k < p.data->GetN(); ++k)
      open(p.data->GetX()[k] - p.data->GetEX()[k], p.data->GetY()[k] + p.data->GetEY()[k]);
    for (int k = 0; p.simulated && k + 1 < p.simulated->GetN(); k += 2)  /* a step: its two ends */
      open(p.simulated->GetX()[k], p.simulated->GetY()[k]);
  }
  for (int row = 0; row < rows; ++row) {
    const double need = high[row] > 0 ? high[row] : 1;
    const auto [tick, labelled] = labelledTicks(need);
    high[row] = (labelled + 0.5) * tick, divisions[row] = labelled + 1;
  }
  const double titlePx = 28;  /* the size of the axis titles of this figure */
  for (size_t i = 0; i < panels.size(); ++i) {
    const Panel &p = panels[i];
    const int column = i % columns, row = i / columns;
    /* the frame's top left corner, px */
    const double frameX = leftPx + column * (frameW + gapPx), frameY = topPx + row * (frameH + gapPx);
    const bool xLabels = i + columns >= panels.size(), yLabels = column == 0;
    /* every frame in a see-through pad over the whole canvas, so all of them are measured in the same pixels */
    TPad *pad = seeThroughPad(cv, frameX / canvasW, 1 - (frameX + frameW) / canvasW, 1 - (frameY + frameH) / canvasH,
                              frameY / canvasH);
    pad->cd();
    const auto measured=photonjet_science::Trace::fromJson(photonjet_science::tag(p.data,"physics_series").at("trace"));
    photonjet_science::annotation(*pad,"dphi/photon_count/"+std::to_string(i),photonjet_science::scalar(p.photons,"weighted photons","region A denominator"),measured);
    photonjet_science::annotation(*pad,"dphi/recoil_reference/"+std::to_string(i),photonjet_science::scalar(p.cut,"radian","displayed recoil boundary"),
      measured.withParameter("renderer/dphi/recoil_boundary/"+p.r->package,"selection",p.cut));
    TH1F *frame = photonjet_science::drawFrame(*pad,x0, 0, M_PI, high[row]);
    for (TAxis *a : {frame->GetXaxis(), frame->GetYaxis()}) a->SetLabelFont(kFont), a->SetLabelSize(kLabelPx);
    tickLengths(frame, 9);
    /* a label at every eighth of the range, as in the IAN: 8 primary divisions, each halved */
    frame->GetXaxis()->SetNdivisions(208, false);
    frame->GetYaxis()->SetNdivisions(500 + divisions[row]);
    if (!yLabels) frame->GetYaxis()->SetLabelSize(0);
    if (!xLabels) frame->GetXaxis()->SetLabelSize(0);
    /* a frame to the right of another leaves out its first label, which would meet the last one of its neighbour */
    for (int eighth = 0; xLabels && eighth <= 8; ++eighth) {
      const std::string label = piLabel(x0 + eighth * (M_PI - x0) / 8);
      frame->GetXaxis()->ChangeLabel(eighth + 1, -1, column > 0 && eighth == 0 ? 0 : -1, -1, -1, -1, label.c_str());
    }
    auto *band = new TBox(p.cut, 0, M_PI, high[row]);
    band->SetFillColorAlpha(gBlue, 0.10);
    band->SetLineWidth(0);
    photonjet_science::draw(band);
    TLine *cut = dashedLine(p.cut, 0, p.cut, high[row], gGray, 2);
    if (p.simulated) photonjet_science::draw(p.simulated,"L");
    photonjet_science::draw(p.data,"P");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(kFont);
    for (const auto &l : p.text) {
      latex.SetTextSize(l.px);
      latex.SetTextColor(l.color);
      latex.DrawLatex((frameX + kInsetPx) / canvasW, 1 - (frameY + l.y) / canvasH, l.text.c_str());
    }
    /* the axis titles as in the IAN, on the outer frames only, each at the end of its axis */
    latex.SetTextSize(titlePx);
    latex.SetTextColor(kBlack);
    if (xLabels) {
      latex.SetTextAlign(33);
      latex.DrawLatex((frameX + frameW) / canvasW, 1 - (frameY + frameH + 46) / canvasH, absolute("#Delta#phi_{#gamma,jet}").c_str());
    }
    if (yLabels) {
      latex.SetTextAlign(33);
      latex.SetTextAngle(90);
      latex.DrawLatex((frameX - 64) / canvasW, 1 - frameY / canvasH,
                      ("1/N_{#gamma} dN_{#gammajet}/d" + absolute("#Delta#phi") + " [rad^{#minus1}#scale[0.3]{ }]").c_str());
    }
    if (i == keyed) {
      const double x1 = (frameX + boxLeft) / canvasW, y2 = 1 - (frameY + boxTop) / canvasH;
      TLegend *leg = boxedLegend(x1, y2 - boxHeight / canvasH, x1 + boxWidth / canvasW, y2, 52 / boxWidth);
      size_t entry = 0;
      leg->AddEntry(p.data, entries[entry++].c_str(), "p");
      if (legendSample) leg->AddEntry(legendSample, entries[entry++].c_str(), "l");
      leg->AddEntry(cut, entries[entry++].c_str(), "l");
      leg->AddEntry(band, entries[entry++].c_str(), "f");
      photonjet_science::draw(leg);
    }
    pad->RedrawAxis();
  }
  save(cv, name);
}

}  /* namespace */

/* The photon and event QA figures and the saved-model figures. They use the helpers above and are called from the
   entry point below. */
#include "QAFigures.C"

/* ---- the entry point --------------------------------------------------------------------------------------------- */
/* results: the results file (or ""), fig: the output folder, systematics and qa: optional systematics and QA files.
   Draw saved quantities with their verified metadata; FIGURES.txt records written and unavailable figures. */
/* ---- the resolution smearing of the simulation (Oct 16) ----------------------------------------------------------- */
/* Justin accepted the resolution-smeared nominal on Oct 16 ("I ACCEPT IT AS THE NEW NOMINAL"): the simulated jets and
   photons matched to truth carry an extra Gaussian width (jets 0.100 x truth pT at R = 0.4, sPH-JET-2026-01; photons
   PPG18), in the response only. This figure shows what that does: the relative xJ resolution of the simulated joint
   response, RMS(reco xJ) / truth xJ in each truth-xJ bin over both photon bins, with the smearing (the results) and
   without it (the reference results, the same selection unsmeared). The response is truth on x, reconstructed on y,
   photon bins times xJ bins, as jointResponse reads it. */
TH1D *responseResolution(TFile &f, const Result &r, const char *role) {
  auto *pairs = dynamic_cast<TH2 *>(get(f, r, "pair_response"));
  auto *axisOf = get(f, r, "xj_per_photon");
  if (!pairs || !axisOf) return nullptr;
  const auto edges = photonEdges(f, r);
  const std::vector<double> xj = edgesOf(axisOf->GetXaxis());
  const int bins = edges.size() - 1, nX = xj.size() - 1;
  if (bins < 1 || nX < 1 || pairs->GetNbinsX() != bins * nX || pairs->GetNbinsY() != bins * nX) return nullptr;
  auto *h = new TH1D(fresh().c_str(), "", nX, xj.data());
  h->SetDirectory(nullptr);
  for (int t = 0; t < nX; ++t) {
    const double truth = 0.5 * (xj[t] + xj[t + 1]);
    double w = 0, m = 0, m2 = 0;
    for (int pt = 0; pt < bins; ++pt)
      for (int pr = 0; pr < bins; ++pr)
        for (int c = 0; c < nX; ++c) {
          const double n = pairs->GetBinContent(pt * nX + t + 1, pr * nX + c + 1), reco = 0.5 * (xj[c] + xj[c + 1]);
          w += n, m += n * reco, m2 += n * reco * reco;
        }
    if (w > 0 && truth > 0) {
      const double mean = m / w;
      h->SetBinContent(t + 1, std::sqrt(std::max(0.0, m2 / w - mean * mean)) / truth);
    }
    h->SetBinError(t + 1, 0);
  }
  photonjet_science::attachTrace(h, photonjet_science::projectionOf(pairs, PlotJson::array())
      .then("relative_rms_of_reco_xj_per_truth_xj_bin", {{"photon_edges", edges}, {"xj_edges", xj}, {"series", role}}));
  photonjet_science::setTag(h, "physics_covariance", photonjet_science::covariance("not_computed", {}, 0,
      "a resolution summary of the response carries no statistical covariance"));
  return h;
}

void resolutionSmearing(TFile &f, TFile &reference, const std::vector<Result> &all, const std::string &name) {
  const Result *pp = nullptr, *central = nullptr;
  for (const auto &r : all) {
    if (r.package == "PP") pp = &r;
    if (r.centrality == "cent0_20") central = &r;
  }
  if (!pp || !central) return missing(name, "needs the p+p and Au+Au 0-20% results");
  TH1D *ppOff = responseResolution(reference, *pp, "pp_without_smearing");
  TH1D *ppOn = responseResolution(f, *pp, "pp_with_smearing");
  TH1D *aaOff = responseResolution(reference, *central, "auau020_without_smearing");
  TH1D *aaOn = responseResolution(f, *central, "auau020_with_smearing");
  if (!ppOff || !ppOn || !aaOff || !aaOn) return missing(name, "needs pair_response in the results and the reference");
  const double low = pp->xjLow, high = pp->xjHigh;
  TGraphErrors *a = points(ppOff, low, high, kGray + 2, 24, 1.6), *b = points(ppOn, low, high, gBlue, 20, 1.6);
  TGraphErrors *c = points(aaOff, low, high, kGray + 2, 25, 1.6), *d = points(aaOn, low, high, gVermilion, 21, 1.6);
  TCanvas cv("cv", "", kCanvasW, kCanvasH);
  xjFrame({{a, "p+p, without resolution smearing", pointKey(a)}, {b, "p+p, with resolution smearing", pointKey(b)},
           {c, "Au+Au 0-20%, without resolution smearing", pointKey(c)},
           {d, "Au+Au 0-20%, with resolution smearing", pointKey(d)}},
          high, "RMS(x_{J#gamma}^{reco}) / x_{J#gamma}^{truth}",
          {"Simulated detector response", "anti-k_{T} R = 0.4, 17 #leq p_{T}^{#gamma} < 33 GeV"});
  gFrame->GetXaxis()->SetTitle("x_{J#gamma}^{truth}");
  headroom(8);
  xjAxes(gFrame, kXTitlePx, kYTitlePx);
  save(cv, name);
}

void IANFigures(const char *results, const char *fig, const char *systematics = "", const char *qa = "",
                const char *only = "", const char *comparisons = "") {
  if (*only && std::string(only) != "combinatoric" && std::string(only) != "ncb_eta_time" &&
      std::string(only) != "systematic_sources" && std::string(only) != "ratio_frames" &&
      std::string(only) != "joint_response" && std::string(only) != "resolution_smearing")
    throw std::runtime_error("unknown figure selector");
  if (std::string(only) == "resolution_smearing" && (!(results && *results) || !*comparisons ||
      (qa && *qa) || (systematics && *systematics)))
    throw std::runtime_error("resolution_smearing requires results and a reference specification only");
  if (std::string(only) == "ncb_eta_time" && ((!qa || !*qa) || (results && *results) || (systematics && *systematics)))
    throw std::runtime_error("ncb_eta_time requires only a QA file");
  if (std::string(only) == "systematic_sources" && (!(results && *results) || !(systematics && *systematics) || (qa && *qa)))
    throw std::runtime_error("systematic_sources requires results and systematics only");
  if (*comparisons && std::string(only) != "combinatoric" && std::string(only) != "resolution_smearing")
    throw std::runtime_error("a comparison specification requires --only combinatoric or resolution_smearing");
  PlotJson comparisonSpec;
  if (*comparisons) {
    std::ifstream input(comparisons);
    if (!input) throw std::runtime_error("cannot open comparison specification");
    input >> comparisonSpec;
    const std::string expected = std::string(only) == "combinatoric" ? "PhotonJetCombinatoricComparisonsV1"
                                                                      : "PhotonJetResolutionReferenceV1";
    if (comparisonSpec.at("schema") != expected)
      throw std::runtime_error("unknown comparison specification schema");
  }
  /* nothing of an earlier call in the same process is kept */
  gWritten.clear();
  gCombinatoricRecipe = PlotJson();
  gCombinatoricFrame = PlotJson();
  gNotDrawn.clear();
  gSystematics.reset();
  gSystematicCoverage = PlotJson();
  kSystematicLabel = "Systematic uncertainty (coverage unverified)";
  gQA.reset();
  gLegacy = false;
  gSavedModelMetadata.clear();
  fa::resetNativeInputs();
  gResultMetadata = gQAMetadata = PlotJson();
  gRenderInputs = gStylePins = PlotJson::object();
  clearFigureState();
  /* the hashes of the drawing code, for the receipts */
  const auto source = std::filesystem::path(__FILE__).parent_path();
  for (const auto &name : {"IANFigures.C", "QAFigures.C", "sPhenixStyle.C", "sPhenixStyle.h"})
    gStylePins[name] = photonjet_render::sha256((source / name).string(), true);
  gStylePins["RenderMetadata.h"] = photonjet_render::sha256((source / "../src/RenderMetadata.h").string(), true);
  gStylePins["PlotPhysics.h"] = photonjet_render::sha256((source / "../src/PlotPhysics.h").string(), true);
  gStylePins["RendererScientificRoster.h"] = photonjet_render::sha256((source / "../src/RendererScientificRoster.h").string(), true);
  gStylePins["RootScientificTrace.h"] = photonjet_render::sha256((source / "../src/RootScientificTrace.h").string(), true);
  gStylePins["RenderScientificParameters.h"] = photonjet_render::sha256((source / "../src/RenderScientificParameters.h").string(), true);
  gStylePins["ScientificTrace.h"] = photonjet_render::sha256((source / "../../TreeToHists/src/ScientificTrace.h").string(), true);
  gStylePins["HistogramScientificDefinitions.h"] = photonjet_render::sha256((source / "../../TreeToHists/src/HistogramScientificDefinitions.h").string(), true);
  SetsPhenixStyle();
  /* the one pattern of every line drawn with the line style kDashed (dashedLine): dashes of 9 canvas pixels, 5 apart
     (a pattern is given in quarter pixels) */
  gStyle->SetLineStyleString(kDashed, "36 20");
  gStyle->SetOptStat(0);
  gStyle->SetImageScaling(kImageScale);
  ianPalette();
  gFig = fig;
  /* the folder is made first: a call that draws nothing still writes FIGURES.txt with the reasons */
  gSystem->mkdir(gFig.c_str(), true);
  /* the input files; one that is named and cannot be opened ends the run */
  const auto open = [](const char *path, std::unique_ptr<TFile> &file) {
    file.reset(TFile::Open(path, "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error(std::string("IANFigures: cannot open ") + path);
  };
  if (systematics && *systematics) open(systematics, gSystematics);
  if (qa && *qa) open(qa, gQA);
  std::unique_ptr<TFile> f;
  if (results && *results) open(results, f);
  if (f && f->Get("product_stage") && std::string(only) != "combinatoric")
    throw std::runtime_error("a pre-unfolding diagnostic may draw only combinatoric panels");
  if (f) gResultMetadata = loadRenderMetadata(*f, results);
  if (gQA) gQAMetadata = loadRenderMetadata(*gQA, qa);
  if (gSystematics) {
    fa::observeNativeInput(*gSystematics);
    gRenderInputs[systematics] = photonjet_render::sha256(systematics, true);
    if (auto *saved = dynamic_cast<TObjString *>(gSystematics->Get("systematic_coverage"))) {
      gSystematicCoverage = PlotJson::parse(saved->GetString().Data());
      if (gSystematicCoverage.at("schema") != "PhotonJetSystematicCoverageV1")
        throw std::runtime_error("unrecognized systematic coverage annotation");
      bool partial = false;
      for (const auto &[result, row] : gSystematicCoverage.at("results").items())
        for (const auto &[source, cell] : row.items())
          partial = partial || cell.at("status") == "not_yet_included";
      // Keep the coverage warning short enough for the result legend beside
      // its header; status wording alone must not open the result's y range.
      kSystematicLabel = partial ? "Partial systematics" : "Systematic uncertainty";
    }
  }
  /* The checks of the inputs against each other, and whether they are reference inputs (gLegacy): a file lacks its
     saved metadata, a package lacks its selection, the systematics file does not name its nominal result, or a file
     says so (the saved analysis name of the results file starts with "legacy", or the recorded input coverage is
     REFERENCE_INPUTS). */
  if (f && gSystematics) {
    auto *nominal = dynamic_cast<TObjString *>(gSystematics->Get("nominal_input_sha256"));
    if (nominal && std::string(nominal->GetString().Data()) != gRenderInputs.at(results).get<std::string>())
      throw std::runtime_error("systematic bands belong to a different nominal result");
    if (!nominal) gLegacy = true;
  }
  if (!gResultMetadata.is_null() && !gQAMetadata.is_null())
    for (const auto &[package, metadata] : gResultMetadata.at("packages").items())
      if (gQAMetadata.at("packages").contains(package) && metadata != gQAMetadata.at("packages").at(package))
        throw std::runtime_error("QA and measurement refer to different sample populations or selections: " + package);
  for (const auto &metadata : {gResultMetadata, gQAMetadata})
    if (!metadata.is_null())
      for (const auto &package : metadata.at("packages"))
        if (package.at("selection").is_null()) gLegacy = true;
  if ((f && gResultMetadata.is_null()) || (gQA && gQAMetadata.is_null())) gLegacy = true;
  if (auto *name = f ? dynamic_cast<TObjString *>(f->Get("analysis")) : nullptr)
    gLegacy = gLegacy || std::string(name->GetString().Data()).find("legacy") == 0;
  if (auto *coverage = f ? dynamic_cast<TObjString *>(f->Get("input_coverage"))
                         : gQA ? dynamic_cast<TObjString *>(gQA->Get("input_coverage")) : nullptr) {
    gPopulation = nlohmann::json::parse(coverage->GetString().Data()).value("status", "UNKNOWN");
    if (gPopulation == "REFERENCE_INPUTS") gLegacy = true;
  } else {
    gPopulation = "UNKNOWN";
  }
  const std::vector<Result> all = f ? resultsOf(*f) : std::vector<Result>();
  if (f && all.empty()) throw std::runtime_error("IANFigures: no results in the supplied file");
  clearFigureState();  /* resultsOf has read the selections of the packages; no figure is being drawn yet */
  /* every xJ axis starts where the final result starts: the common reported lower edge of its results */
  {
    double common = 0;
    for (const auto &r : all)
      if (r.package == "PP" || r.centrality == "cent50_80" || r.centrality == "cent0_20") common = std::max(common, r.xjLow);
    gXjShownLow = xjShownLow(common);
  }
  /* the resolution-smearing figure: this results file against its unsmeared reference, recorded among the inputs */
  std::unique_ptr<TFile> reference;
  if (std::string(only) == "resolution_smearing") {
    const std::string path = comparisonSpec.at("reference_results").get<std::string>();
    reference.reset(TFile::Open(path.c_str(), "READ"));
    if (!reference || reference->IsZombie()) throw std::runtime_error("IANFigures: cannot open " + path);
    loadRenderMetadata(*reference, path);
    resolutionSmearing(*f, *reference, all, "Unfolding/xj_response_resolution_smearing.png");
  }
  /* the figures of every result, under the IAN's folder and file names */
  for (const auto &r : all) {
    if (std::string(only) == "resolution_smearing") break;
    const bool pp = r.fileTag == "pp";
    if (std::string(only) == "joint_response") {
      if (r.package == "PP" || r.centrality == "cent0_20")
        jointResponse(*f, r, "Unfolding/" + r.longFileTag + "_joint_response.png");
      continue;
    }
    if (std::string(only) == "systematic_sources") {
      systematicSources(*f, r, "Systematics/" + r.fileTag + "_xj_systematic_sources.png");
      continue;
    }
    if (std::string(only) == "ratio_frames") {
      beforeAfter(*f, r, "Unfolding/" + r.fileTag + "_xj_before_after_unfolding.png");
      iterationOverlay(*f, r, "Unfolding/" + r.fileTag + "_xj_iteration_overlay.png");
      detectorLevel(*f, r, "xJ_generalQA/" + r.fileTag + "_xj_detector_level_data_vs_simulation.png");
      systematicVariations(*f, r, "Systematics/" + r.fileTag + "_xj_variations_");
      photonClosure(*f, r, "Unfolding/" + r.fileTag + "_photon_closure_overlay.png", false);
      photonClosure(*f, r, "Unfolding/" + r.fileTag + "_photon_half_closure_overlay.png", true);
      photonYield(*f, r, "Unfolding/" + r.fileTag + "_photon_yield_unfolding_nominal_2bin_k" + std::to_string(r.photonIterations) + ".png");
      pairClosure(*f, r, "Unfolding/" + r.fileTag + "_xj_closure_overlay.png", false);
      pairClosure(*f, r, "Unfolding/" + r.fileTag + "_xj_half_closure_overlay.png", true);
      jointRefold(*f, r, "Unfolding/" + r.fileTag + "_joint_refold.png");
      continue;
    }
    if (*only) {
      const std::string baseline = "purityAndCombSub/" + r.fileTag + "_xj_combinatoric_subtraction.png";
      gCombinatoricFrame = comparisonSpec.is_object() && comparisonSpec.contains("display_y_ranges")
                              ? comparisonSpec.at("display_y_ranges").at(r.key) : PlotJson();
      combinatoric(*f, r, baseline);
      PlotJson linked = PlotJson::array();
      if (!comparisonSpec.is_null()) for (const auto &alternative : comparisonSpec.at("comparisons")) {
        const std::string id = alternative.at("id"), path = alternative.at("results");
        if (id.empty() || id.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos)
          throw std::runtime_error("unsafe comparison id");
        if (photonjet_render::sha256(path, true) != alternative.at("sha256").get<std::string>())
          throw std::runtime_error("comparison result hash changed: " + id);
        const auto nominalMetadata = gResultMetadata, nominalInputs = gRenderInputs;
        const auto nominalAncestry=fa::nativeInputs();
        std::unique_ptr<TFile> other;
        open(path.c_str(), other);
        gResultMetadata = loadRenderMetadata(*other, path.c_str());
        if (gResultMetadata.at("packages") != nominalMetadata.at("packages"))
          throw std::runtime_error("comparison input packages differ from the nominal: " + id);
        const auto alternatives = resultsOf(*other);
        const auto found = std::find_if(alternatives.begin(), alternatives.end(), [&](const Result &a) { return a.key == r.key; });
        if (found == alternatives.end() || found->xjLow != r.xjLow || found->xjHigh != r.xjHigh)
          throw std::runtime_error("comparison result/range differs from the nominal: " + id);
        clearFigureState();
        const std::string output = "comparisons/" + id + "/" + baseline;
        gCombinatoricRecipe = {{"role", "comparison"}, {"id", id}, {"baseline_output", baseline},
                              {"output", output}, {"title", alternative.at("title")},
                              {"difference", alternative.at("difference")}, {"header", alternative.at("header")}};
        if (const auto *stage = dynamic_cast<TObjString *>(other->Get("product_stage")))
          gCombinatoricRecipe["product_stage"] = stage->GetString().Data();
        combinatoric(*other, *found, output);
        linked.push_back({{"id", id}, {"title", alternative.at("title")}, {"output", output},
                          {"difference", alternative.at("difference")}});
        gCombinatoricRecipe = PlotJson();
        fa::nativeInputs()=nominalAncestry;
        gResultMetadata = nominalMetadata;
        gRenderInputs = nominalInputs;
        clearFigureState();
      }
      if (!linked.empty()) {
        const std::string path = gFig + "/" + baseline + ".json";
        std::ifstream input(path);
        PlotJson receipt; input >> receipt;
        receipt["comparison_outputs"] = linked;
        receipt["comparison_specification"] = {{"path", comparisons}, {"sha256", photonjet_render::sha256(comparisons, true)}};
        photonjet_render::write(path, receipt);
      }
      continue;
    }
    purityCorrection(*f, r, "purityAndCombSub/" + r.fileTag + "_xj_purity_correction.png");
    combinatoric(*f, r, "purityAndCombSub/" + r.fileTag + "_xj_combinatoric_subtraction.png");
    photonResponse(*f, r, "Unfolding/" + r.longFileTag + "_photon_response.png");
    jointResponse(*f, r, "Unfolding/" + r.longFileTag + "_joint_response.png");
    photonIterationRule(*f, r, "Unfolding/" + r.fileTag +
                                   (pp ? "_photon_yield_iteration_stability_publication.png"
                                       : "_photon_yield_iteration_stability.png"));
    jointIterationRule(*f, r, "Unfolding/" + r.fileTag + "_joint_iteration_stability.png");
    photonClosure(*f, r, "Unfolding/" + r.fileTag + "_photon_closure_overlay.png", false);
    photonClosure(*f, r, "Unfolding/" + r.fileTag + "_photon_half_closure_overlay.png", true);
    photonYield(*f, r, "Unfolding/" + r.fileTag + "_photon_yield_unfolding_nominal_2bin_k" +
                           std::to_string(r.photonIterations) + ".png");
    priorReweighting(*f, r, "Unfolding/" + r.fileTag + "_xj_prior_reweighting.png");
    beforeAfter(*f, r, "Unfolding/" + r.fileTag + "_xj_before_after_unfolding.png");
    iterationOverlay(*f, r, "Unfolding/" + r.fileTag + "_xj_iteration_overlay.png");
    statisticalBudget(*f, r, "Unfolding/" + r.fileTag + "_xj_statistical_budget.png");
    pairClosure(*f, r, "Unfolding/" + r.fileTag + "_xj_closure_overlay.png", false);
    pairClosure(*f, r, "Unfolding/" + r.fileTag + "_xj_half_closure_overlay.png", true);
    jointRefold(*f, r, "Unfolding/" + r.fileTag + "_joint_refold.png");
    /* against ATLAS: p+p against p+p, and the Au+Au class that starts at 0% against 0-10% Pb+Pb */
    if (pp || r.centralityLow == 0)
      atlasOverlay(*f, r, pp ? "pp" : "pbpb_0_10", "LHCComparison/" + r.fileTag + "_xj_atlas_overlay.png");
    systematicSources(*f, r, "Systematics/" + r.fileTag + "_xj_systematic_sources.png");
    systematicVariations(*f, r, "Systematics/" + r.fileTag + "_xj_variations_");
    detectorLevel(*f, r, "xJ_generalQA/" + r.fileTag + "_xj_detector_level_data_vs_simulation.png");
    /* for p+p only; the file keeps the IAN's name, and the pads are the result's photon bins */
    if (pp) recoTruthAssociation(*f, r, "xJ_generalQA/pp_xj_data_reco_vs_truth_tagged_2x3.png", gQA.get());
  }
  /* the figures of several results; the Delta-phi grid keeps the IAN's file name whatever the number of its frames */
  if (f && !*only) {
    photonJetDphi(*f, all, "xJ_generalQA/pp_auau_gamma_jet_dphi_2x2.png");
    ppVsCentral(*f, all, "Unfolding/pp_vs_auau_xj_unfolded.png");
    centralityOverlay(*f, all, "FinalResults/pp_vs_auau_centrality_xj_unfolded.png");
    auauCentralities(*f, all, "FinalResults/auau_centrality_xj_unfolded.png");
    /* the radius dependence, for the two results of the IAN's figure: p+p and the Au+Au class that starts at 0% */
    for (const auto &r : all)
      if (r.fileTag == "pp" || r.centralityLow == 0)
        radiusOverlay(*f, results, r, "xJ_generalQA/" + r.fileTag + "_xj_unfolded_radius_overlay.png");
  }
  if (gQA && std::string(only) == "ncb_eta_time") {
    ncbEtaTime();
    ncbEtaTimeAfter("PP");
    ncbEtaTimeAfter("AuAu");
  }
  if (gQA && !*only) {
    auto *kind = dynamic_cast<TObjString *>(gQA->Get("analysis"));
    const std::string analysis = kind ? kind->GetString().Data() : "";
    /* A saved-model QA file holds model statistics and no event histograms, so it has figures of its own. */
    if (analysis == "saved-model-qa" || analysis == "saved_model_physics_review") {
      modelQAFigures();
      bdtPerClass();
    } else {
      qaFigures();
    }
  }
  /* FIGURES.txt is read by the pipeline (src/Pipeline.cc): a written figure is a line that starts with fig/, a figure
     that was not drawn a line that starts with "# not drawn: ". The paths are written relative to the output folder of
     a run, whose figure folder is named fig, whatever folder the `fig` argument names. */
  std::ofstream list((gFig + "/FIGURES.txt").c_str());
  list << "# Figures written by FinalAnalysis/plots/IANFigures.C from " << (f ? results : "no results file")
       << (gQA ? std::string(" and ") + qa : "") << "\n";
  for (const auto &w : gWritten) list << "fig/" << w << "\n";
  for (const auto &n : gNotDrawn) list << "# not drawn: fig/" << n << "\n";
  if (!list) throw std::runtime_error("IANFigures: cannot write FIGURES.txt");
  printf("IANFigures: %zu figures in %s, %zu not drawn\n", gWritten.size(), fig, gNotDrawn.size());
  /* the files this function owns are closed here, before ROOT's cleanup at the end of the process would delete them */
  gQA.reset();
  gSystematics.reset();
}

/* final_analysis.sh builds this file as an ordinary program, which starts here. Under ROOT (IANFigures.C+ or the
   interpreter) the function above is the entry point and this is left out. */
#if !defined(__CLING__) && !defined(__ACLIC__)
int main(int argc, char **argv) {
  if (argc != 5 && argc != 7) {
    fprintf(stderr, "usage: figures RESULTS FIGURES SYSTEMATICS QA   (\"\" for a file that is not given)\n");
    return 2;
  }
  gROOT->SetBatch(true);  /* canvases are drawn into files, never onto a screen */
  try {
    IANFigures(argv[1], argv[2], argv[3], argv[4], argc == 7 ? argv[5] : "", argc == 7 ? argv[6] : "");
  } catch (const std::exception &error) {
    fprintf(stderr, "IANFigures: %s\n", error.what());
    return 1;
  }
  return 0;
}
#endif
