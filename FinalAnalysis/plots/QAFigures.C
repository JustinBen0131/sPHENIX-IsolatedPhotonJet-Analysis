/* QAFigures: the photon and event QA figures of the IAN, from the QA file of final_analysis.sh qa. */

/* Included by IANFigures.C, whose helpers it uses (header, save, style), and called from its entry point with the QA
   file. Each figure draws the saved sums and numerics under qa/<package>/ of that file: the DATA and simulated QA
   histograms (data/, phosim/, jetsim/; the simulation summed over samples, each weighted by its cross section per
   generated event) and the ABCD counts, purity, purity closure, stitching fits, working points, trigger study and
   centrality target of the qa step. Drawing projects, rebins and normalises for display only. The classes, samples and
   triggers are the ones in the file, except where a figure's comment names its own ranges. A figure whose inputs are
   absent or empty is not drawn and is listed in FIGURES.txt with the reason; nothing is invented.

   The figures of the first half of the file, by output file (fig/README.md maps the file names to the figures of the
   note), and the function that draws each; qaFigures, at the end of the file, calls them all:
     simulation_stitching/auau_centrality_data_embeddedMC.png            centralityReweighting
     simulation_stitching/auau_vertex_data_embeddedMC.png                vertexReweighting
     event_obs/centrality.pdf, vertexz.pdf, cluster_pt.pdf, jet_pt.pdf   eventProperties (through eventDistribution)
     event_obs/cluster_eta_phi_<sample>.pdf, jet_eta_phi_<sample>.pdf    eventProperties (through etaPhiMap)
     photonOnly/auau_minbias_pretraining_before_after.png                minimumBias
     photonOnly/<system>_ncb_score_data_mc.png                           ncbScore
     photonOnly/<pp|auau>_ncb_eta_time_before.png                        ncbEtaTime
     photonOnly/<pp|auau>_ncb_time_selections.png                        ncbTimeSelections
     photonOnly/<pp|auau>_ncb_eta_time_after.png                         ncbEtaTimeAfter
     shower_shapes/dis_<shape>_cent<k>_pt0_cut1.pdf                      showerShapes
     photonOnly/auau_<signal|background|data_simulation>_<full_cluster|3x3>_width.png,
       photonOnly/auau_width_vs_centrality.png                            widths
     photonOnly/auau_bdt_score_perCent.png                               bdtPerClass
     shower_shapes/e11_e33_selection_flow.png                            e11Flow
     photonOnly/pp_auau_truth_isolation_direct_fragmentation_2x2.png     truthIsolation
     photonOnly/<system>_truth_photon_pt_after_isolation.png             truthPtAfterIsolation
     isolation/isoCentOverlay_<noUE|withUE>_pT_14_16.png                 isolationOverlays
     isolation/<system>_isoStack.png                                     isolationStacks
     photonOnly/pp_auau_abcd_region_counts_2x2.png                       abcdCounts
     photonOnly/data_purity_and_prompt_leakage_2x4_style.png             purityPanels ("purity")
     photonOnly/purity_qa_pp_auau_simCLOSURE.png                         purityPanels ("closure")
     photonOnly/pp_si_di_auau_photon_efficiency_2x2.png                  efficiency
   The simulation-stitching, trigger, working-point, score-separation and saved-model figures (stitching, triggerStudy,
   workingPoints, scoreSeparation, modelQAFigures) stand in the second half, from plainBars on.

   The saved histograms these figures read, under qa/<package>/, and their axes in the order x, y, z, which a Cut
   names 0, 1, 2. <family>: data, phosim (photon+jet simulation) or jetsim (inclusive-jet simulation); <class>: "all"
   for p+p, centA_B for Au+Au; <label>: the candidate population (labelsOf). TreeToHists fills them, and the head
   comment of its FillPhotonQA.C is the authority on their content:
     <family>/photons/<class>/<label>/
       shape_<var>            the variable, E_T, the selection stage: [0, 1) every candidate, [1, 2) preselected,
                              [2, 3) and tight, [3, 4) and isolated; a candidate fills every stage it reached
       ncb                    the NCB score, E_T: every candidate
       iso                    the isolation, E_T, the photon ID of a preselected candidate: [0, 1) tight, [1, 2)
                              non-tight
       clu_et, clu_etaphi     E_T and the stage; eta, phi and the stage
       e11e33_dt, etatime     DATA only: E11/E33, E_T and cluster time minus MBD time [ns]; eta and that time
       time_stage, etatime_stage   DATA only: that time, E_T and the stage; eta, that time and the stage (stage axis
                              as for shape_<var>; stage 0 holds the candidates of etatime)
       etatime_npb            DATA only, when the fill has it: eta, that time and the NCB score, the candidates of
                              etatime
     <family>/slices/<label>/bdt, isoslice, isonoueslice, width_<var>, width_tight_iso_<var>
                              the value, E_T, the 5% centrality slice by its number: slice k holds 5k to 5(k + 1) %
     phosim/truth/<class>/
       truth_iso              truth isolation, truth pT, the prompt class: [1, 2) direct, [2, 3) fragmentation
       eff                    truth pT, the stage of an analysis-signal truth photon: [0, 1) every one, [1, 2) matched
                              to a candidate, [2, 3) and preselected, [3, 4) preselected and isolated, [4, 5)
                              preselected and tight, [5, 6) preselected, tight and isolated; the E_T of the matched
                              candidate (below zero for a photon without one)
     <family>/events/<class>/evt_vz, jet_pt, jet_etaphi
                              vertex z; jet pT, or jet eta and phi, and the events: [0, 1) every analysis event, [1, 2)
                              those with a tight, isolated leading photon
     <family>/events/centrality   Au+Au: centrality, the weight: [0, 1) the event weight alone, [1, 2) with the
                              centrality weight
     <family>/events/calo_total   total calorimeter energy, centrality, the minimum-bias state (minimumBias)
   and, computed by the qa step: abcd/<class> (the region A to D as 0 to 3, E_T), purity/<class>/ (raw, corrected,
   leakage_B, _C, _D against E_T), closure/<class>/ (raw, corrected, truth, ratio, ratio_r_bkg against E_T) and
   centrality_target. */
#include <TColor.h>
#include <TBox.h>
#include <TF1.h>
#include <TLeaf.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <functional>
#include <limits>
#include <tuple>

namespace {

/* ---- the QA file and saving a figure --------------------------------------------------------------------------- */

/* The QA file, opened by the entry point of IANFigures.C. */
std::unique_ptr<TFile> gQA;

/* Defined with the working-point and saved-model figures in the second half of the file and called before them; the
   declaration of saveModelQA stands here for its default arguments. */
void saveModelQA(TCanvas &, const std::string &, const PlotJson &, const std::vector<std::string> &,
                 const PlotJson & = PlotJson(), const char * = __builtin_FUNCTION());
void saveReleasedWorkingPoints(TCanvas &, const std::string &, const PlotJson &, const std::string &, const char * = __builtin_FUNCTION());
bool savedScoreReview(bool truthBackground = false);
std::string modelPhotonSimulationLabel(const std::string &);
/* The receipt helpers, defined there too: a figure of the first half that has a supplement uses them. */
template <class Edit> void amendReceipt(const std::string &path, Edit edit);
PlotJson plotRecipe(const std::string &figure, const std::string &id, const std::string &role,
                    const std::string &output, const std::string &baseline, const std::string &function);

/* Saves a QA figure. Before any file is written, every package the figure draws from must carry physics event
   weights: a figure that fails is not drawn and is listed with that reason, and the other figures of the file are
   drawn (a reference-stitching QA file holds raw generator counts). save() makes the same check and writes the
   weight records into the figure's receipt. function: the name of the calling function, recorded in that receipt. */
void savePhysicsQA(TCanvas &canvas, const std::string &path, const char *function = __builtin_FUNCTION(), const char* family=__builtin_FUNCTION()) {
#ifdef PHOTONJET_NATIVE_CAPTURE_TEST
  if(gNativeCaptureSink)return save(canvas,path,PlotJson(),function,"plots/QAFigures.C",family);
#endif
  try {
    if (gSemanticPackages.empty()) photonjet_render::requirePhysicsPlotWeights(PlotJson(), path);
    for (const auto &key : gSemanticPackages) {
      const bool fromQA = key.rfind("qa/", 0) == 0;
      const auto pkg = key.substr(key.find('/') + 1);
      const auto &metadata = fromQA ? gQAMetadata : gResultMetadata;
      photonjet_render::requirePhysicsPlotWeights(!metadata.is_null() && metadata.at("packages").contains(pkg)
        ? metadata.at("packages").at(pkg).value("weight_provenance", PlotJson()) : PlotJson(), path);
    }
  } catch (const std::runtime_error &refused) {
    /* the record of the refused figure is cleared, as save() clears it, so that the next figure starts empty */
    gVisibleText.clear();
    gSemanticPackages.clear();
    gProjectionSteps = PlotJson::array();
    gInputObservations = PlotJson::array();
    return missing(path, refused.what());
  }
  save(canvas, path, PlotJson(), function, "plots/QAFigures.C",family);
}

/* The IAN palette of IANFigures.C under the names this file uses; references, because ianPalette() sets the colours
   when the drawing starts. cOrange is the palette's vermilion (gVermilion), not its gOrange. */
const int cBlack = kBlack;
const int &cBlue = gBlue, &cOrange = gVermilion, &cGreen = gGreen, &cRed = gRed, &cPurple = gPurple;

/* True when no bin inside the axes holds a content or an error. */
bool withoutEntries(const TH1 *h) {
  for (int b = 0; b < h->GetNcells(); ++b)
    if (!h->IsBinUnderflow(b) && !h->IsBinOverflow(b) && (h->GetBinContent(b) != 0 || h->GetBinError(b) != 0))
      return false;
  return true;
}

/* qa/<path> of the QA file, detached; null when absent or without entries. Reading a histogram binds the figure being
   drawn to its package, the first folder of the path: save() then checks that package's weights and records its
   selection in the receipt. photonSupport, selectionValue, the label helpers below and energyOf bind their package in
   the same way. */
TH1 *qa(const std::string &path) {
  auto *h = dynamic_cast<TH1 *>(gQA->Get(("qa/" + path).c_str()));
  if (!h || withoutEntries(h)) return nullptr;
  renderSelection(path.substr(0, path.find('/')), true);  /* called for the binding; the selection is not used here */
  h = static_cast<TH1 *>(h->Clone(fresh().c_str()));
  h->SetDirectory(nullptr);
  bindHistogramSource(h, gQA->GetName(), "qa/" + path);
  return h;
}

/* a + b, added into a and returned; either may be null. */
TH1 *joined(TH1 *a, TH1 *b) {
  if (!a) return b;
  if (b) {
    if (photonjet_science::tag(a,"physics_constituents").is_null())
      photonjet_science::sumDefinition(a,a,1.0,true);
    photonjet_science::sumDefinition(a,b,1.0,false);
    const PlotJson left = histogramSources(a), right = histogramSources(b);
    if (auto *tag = a->GetListOfFunctions()->FindObject("ian_input_sources"))
      a->GetListOfFunctions()->Remove(tag);  // ROOT owns retained objects; never delete evidence here
    if (!left.empty() && !right.empty()) {
      PlotJson sources = left;
      for (const auto &source : right) sources.push_back(source);
      a->GetListOfFunctions()->Add(new TNamed("ian_input_sources", sources.dump().c_str()));
    }
    a->Add(b);
  }
  return a;
}

/* ---- classes and candidates ------------------------------------------------------------------------------------ */

/* The centrality range in percent of a class id, from the saved selection of whichever package has that class; null
   for a class no package has and without metadata. */
PlotJson classRange(const std::string &c) {
  if (gQAMetadata.is_null()) return PlotJson();
  for (const auto &package : gQAMetadata.at("packages")) {
    const auto &selection = package.at("selection");
    const auto &classes = selection.is_null() ? package.at("classes") : selection.at("qa").at("classes");
    for (const auto &item : classes)
      if (item.at("id").get<std::string>() == c) return item.at("centrality_percent");
  }
  return PlotJson();
}
/* The lower and the upper edge of a class in percent of centrality; without a saved range they are read from the id
   centA_B, and a class that is not a centrality class (p+p "all") spans 0 to 100. */
double lowOf(const std::string &c) {
  const auto range = classRange(c);
  return !range.is_null() ? range.at(0).get<double>() : c.rfind("cent", 0) == 0 ? std::atof(c.c_str() + 4) : 0;
}
double highOf(const std::string &c) {
  const auto range = classRange(c);
  if (!range.is_null()) return range.at(1).get<double>();
  const auto u = c.find('_');
  return c.rfind("cent", 0) == 0 && u != std::string::npos ? std::atof(c.c_str() + u + 1) : 100;
}

/* The classes of a package, most central first: p+p "all", Au+Au centA_B. They are its DATA photon directories; a
   package without DATA takes them from its simulation, so its simulation-only figures are drawn. */
std::vector<std::string> classesOf(const std::string &pkg) {
  std::vector<std::string> out;
  for (const char *family : {"data", "phosim", "jetsim"}) {
    auto *d = gQA->GetDirectory(("qa/" + pkg + "/" + family + "/photons").c_str());
    if (!d) continue;
    TIter next(d->GetListOfKeys());
    while (auto *k = static_cast<TKey *>(next()))
      if (std::find(out.begin(), out.end(), k->GetName()) == out.end()) out.push_back(k->GetName());
    if (!out.empty()) break;
  }
  std::sort(out.begin(), out.end(), [](const std::string &a, const std::string &b) { return lowOf(a) < lowOf(b); });
  return out;
}

/* A centrality range as text, the two numbers joined by a hyphen: "0-10%". */
std::string percent(double a, double b) { return Form("%g-%g%%", a, b); }

/* A centrality range as the IAN prints it, the two numbers joined by a dash. */
std::string centralityRange(double low, double high) {
  return Form("%g#minus%g%%", low, high);
}

/* The system of a class as text: "Au+Au 0-10%" for a centrality class, "p+p" otherwise. */
std::string classLabel(const std::string &c) {
  return c.rfind("cent", 0) == 0 ? "Au+Au " + percent(lowOf(c), highOf(c)) : "p+p";
}

/* A centrality class as a file-name tag: "0to10". */
std::string classTag(const std::string &c) { return Form("%gto%g", lowOf(c), highOf(c)); }

/* All Au+Au classes together as a file-name tag: "0to80"; empty when the file has no Au+Au class. */
std::string spanTag() {
  const auto c = classesOf("AuAu");
  return c.empty() ? "" : Form("%gto%g", lowOf(c.front()), highOf(c.back()));
}

/* p+p, then the Au+Au classes from peripheral to central (the panels of the photon figures). */
std::vector<std::pair<std::string, std::string>> panelOrder() {
  std::vector<std::pair<std::string, std::string>> out;
  for (const auto &c : classesOf("PP")) out.push_back({"PP", c});
  const auto a = classesOf("AuAu");
  for (auto it = a.rbegin(); it != a.rend(); ++it) out.push_back({"AuAu", *it});
  return out;
}
/* p+p, then the Au+Au classes from central to peripheral (one figure each). */
std::vector<std::pair<std::string, std::string>> figureOrder() {
  std::vector<std::pair<std::string, std::string>> out;
  for (const auto &c : classesOf("PP")) out.push_back({"PP", c});
  for (const auto &c : classesOf("AuAu")) out.push_back({"AuAu", c});
  return out;
}

/* The candidate populations saved for a family, as their folder names under photons/<class>/ and slices/. DATA has
   one, "data". The photon+jet simulation has "sig", the candidates matched to an analysis-signal truth photon, and
   "phoother", every other candidate; the inclusive-jet simulation has "bkg", the candidates not matched to an
   analysis-signal truth photon, and "jetprompt", those matched to one. */
std::vector<std::string> labelsOf(const std::string &family) {
  if (family == "data") return {"data"};
  if (family == "phosim") return {"sig", "phoother"};
  return {"bkg", "jetprompt"};
}

/* The candidates of a family in one class: photons/<class>/<label>/<stem>, summed over the family's labels. */
TH1 *candidates(const std::string &pkg, const std::string &family, const std::string &c, const std::string &stem) {
  TH1 *total = nullptr;
  for (const auto &l : labelsOf(family))
    total = joined(total, qa(pkg + "/" + family + "/photons/" + c + "/" + l + "/" + stem));
  return total;
}

/* The sum over the classes of a package of the histogram that `fetch` returns for each class; null when it returns
   none. */
TH1 *overClasses(const std::string &pkg, const std::function<TH1 *(const std::string &)> &fetch) {
  TH1 *total = nullptr;
  for (const auto &c : classesOf(pkg)) total = joined(total, fetch(c));
  return total;
}

/* The candidates of a family summed over the classes of its package: photons/<class>/<label>/<stem>. */
TH1 *allCandidates(const std::string &pkg, const std::string &family, const std::string &stem) {
  return overClasses(pkg, [&](const std::string &c) { return candidates(pkg, family, c, stem); });
}

/* An event histogram of a family summed over the classes of its package: events/<class>/<stem>. */
TH1 *allEvents(const std::string &pkg, const std::string &family, const std::string &stem) {
  return overClasses(pkg, [&](const std::string &c) { return qa(pkg + "/" + family + "/events/" + c + "/" + stem); });
}

/* ---- numbers and label texts of the saved selection ------------------------------------------------------------ */

/* The collision energy of a package as text (sqrt(s) for p+p, sqrt(s_NN) for Au+Au); empty without selection
   metadata. */
std::string energyOf(const std::string &pkg) { return energyLabel(pkg, true); }

/* The photon E_T range of a package's selection in GeV; 15 to 35 for a file without selection metadata. */
std::pair<double, double> photonSupport(const std::string &pkg) {
  const auto &m = renderSelection(pkg, true);
  if (m.is_null()) return {15, 35};
  return {m.at("photon").at("et_min_gev"), m.at("photon").at("et_max_gev")};
}

/* A photon E_T window as text, with the same space on both sides of each relation. */
std::string photonWindowLabel(double low, double high) {
  return Form("%g #leq E_{T}^{#gamma,reco} < %g GeV", low, high);
}

/* The photon pseudorapidity acceptance of a package as text with sized bars; empty for a file without selection
   metadata, which the receipt then records. */
std::string photonEtaLabel(const std::string &pkg) {
  const auto &m = renderSelection(pkg, true);
  return m.is_null() ? forReceipt("Selection metadata unavailable")
                     : Form("#left|#eta^{#gamma}#right| < %g", m.at("photon").at("abs_eta_max").get<double>());
}

/* One number of the saved selection of a package, selection[section][key]; `fallback` for a file without selection
   metadata. */
double selectionValue(const std::string &pkg, const char *section, const char *key, double fallback) {
  const auto &m = renderSelection(pkg, true);
  return m.is_null() ? fallback : m.at(section).at(key).get<double>();
}

/* The QA drawing windows of a package, from its saved selection (qa.windows, written by the qa step from the settings
   QA.Windows): one window [low, high] in GeV, or a list of them; `fallback` for a file whose selection does not carry
   the key. The ABCD photon E_T bins of the qa step the same way (qa.abcd_edges_gev, the setting QA.AbcdEdges). */
const PlotJson *qaWindowRecord(const std::string &pkg, const char *key) {
  const auto &m = renderSelection(pkg, true);
  if (m.is_null() || !m.contains("qa") || !m.at("qa").contains("windows") || !m.at("qa").at("windows").contains(key))
    return nullptr;
  return &m.at("qa").at("windows").at(key);
}
std::pair<double, double> qaWindow(const std::string &pkg, const char *key, std::pair<double, double> fallback) {
  const PlotJson *w = qaWindowRecord(pkg, key);
  if (!w || !w->is_array() || w->size() != 2) return fallback;
  return {w->at(0).get<double>(), w->at(1).get<double>()};
}
std::vector<std::pair<double, double>> qaWindowList(const std::string &pkg, const char *key,
                                                    const std::vector<std::pair<double, double>> &fallback) {
  const PlotJson *w = qaWindowRecord(pkg, key);
  if (!w || !w->is_array() || w->empty()) return fallback;
  std::vector<std::pair<double, double>> out;
  for (const auto &pair : *w) {
    if (!pair.is_array() || pair.size() != 2) return fallback;
    out.push_back({pair.at(0).get<double>(), pair.at(1).get<double>()});
  }
  return out;
}
std::vector<double> qaAbcdEdges(const std::string &pkg, const std::vector<double> &fallback) {
  const auto &m = renderSelection(pkg, true);
  if (m.is_null() || !m.contains("qa") || !m.at("qa").contains("abcd_edges_gev")) return fallback;
  std::vector<double> out;
  for (const auto &e : m.at("qa").at("abcd_edges_gev")) out.push_back(e.get<double>());
  return out.size() >= 2 ? out : fallback;
}

/* The radius of the reconstructed isolation cone of a package; for a file without selection metadata 0.4 for p+p and
   0.3 for Au+Au. */
double isolationRadius(const std::string &pkg) {
  const auto &m = renderSelection(pkg, true);
  return m.is_null() ? (pkg == "PP" ? 0.4 : 0.3) : m.at("photon_id").at("isolation").at("radius").get<double>();
}

/* The photon pseudorapidity acceptance as text, "|eta^gamma| < 0.7", with plain bars of the text height. */
std::string photonEtaLine(double most) {
  return Form("%s < %g", absolute("#eta^{#gamma}").c_str(), most);
}

/* The same for the saved acceptance of a package; empty without selection metadata. */
std::string photonEtaLine(const std::string &pkg) {
  const auto &m = renderSelection(pkg, true);
  return m.is_null() ? "" : photonEtaLine(m.at("photon").at("abs_eta_max").get<double>());
}

/* The photon window as a line of a shared-axes figure, "15 <= E_T^{gamma,reco} < 35 GeV". It differs from
   photonWindowLabel by one space: the sign after the long superscript follows it without a word space, which the
   renderer already leaves there at the sizes of these figures. */
std::string gridWindowLine(double low, double high) {
  return Form("%g #leq E_{T}^{#gamma,reco}< %g GeV", low, high);
}

/* The legend text and the colour of a simulation family in the figures that show the two families side by side under
   their plain names (ncbScore, showerShapes, e11Flow). Other figures word or colour the two samples differently. */
struct SampleLook { const char *text; int color; };
SampleLook sampleLook(const std::string &family) {
  return family == "phosim" ? SampleLook{"Photon+jet simulation", cRed} : SampleLook{"Inclusive-jet simulation", cBlue};
}

/* A count with thousands separators: 27156 as "27,156". */
std::string grouped(double count) {
  std::string digits = Form("%.0f", count);
  for (int at = static_cast<int>(digits.size()) - 3; at > 0; at -= 3) digits.insert(at, ",");
  return digits;
}

/* ---- bins and projections -------------------------------------------------------------------------------------- */

/* The bins of one axis: values and variances. */
struct Bins {
  std::vector<double> edges, v, e2;
  photonjet_science::Trace physics = photonjet_science::Trace::unavailable("manually constructed bins need an explicit source binding");
  std::vector<std::string> measurementState;
  std::string seriesRole;
  std::vector<double> hiddenZeroValues;
  PlotJson covariance = photonjet_science::covariance("not_computed", {}, 0, "renderer propagates diagonal variances; no cross-bin covariance calculated");
  int n() const { return v.size(); }
  double center(int i) const { return 0.5 * (edges[i] + edges[i + 1]); }
  double width(int i) const { return edges[i + 1] - edges[i]; }
  bool inside(int i, double lo, double hi) const { return edges[i] >= lo - 1e-9 && edges[i + 1] <= hi + 1e-9; }
  double sum() const {
    double s = 0;
    for (double x : v) s += x;
    return s;
  }
};

/* A window on one axis of a saved histogram, axis 0 (x), 1 (y) or 2 (z): the whole bins inside [lo, hi] are kept. */
struct Cut { int axis; double lo, hi; };

/* The slices/ histograms carry the centrality on their third axis in slices of 5%, numbered from 0: a window on that
   axis is a centrality range in percent divided by this. */
const double kSlicePercent = 5;

/* The sum of a histogram over every axis but `axis` (no flow bins), the other axes restricted to the bins inside the
   windows given. Each projection is recorded, with its windows, in the receipt of the figure being drawn. */
Bins project(const TH1 *h, int axis, const std::vector<Cut> &windows = {}) {
  if(!h)throw std::runtime_error("scientific projection requires its retained input histogram");
  const TAxis *axes[3] = {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()};
  const int dim = h->GetDimension();
  PlotJson cuts = PlotJson::array();
  for (const auto &w : windows)
    cuts.push_back({{"axis", w.axis}, {"low", std::max(w.lo,axes[w.axis]->GetXmin())},
                   {"high", std::min(w.hi,axes[w.axis]->GetXmax())}, {"interval", "[low,high)"}});
  const PlotJson sources = histogramSources(h);
  gProjectionSteps.push_back(
      {{"histogram", h->GetTitle()}, {"sources", sources}, {"retained_axis", axis}, {"windows", cuts}, {"flow_bins", "excluded"}});
  Bins out;
  out.physics = photonjet_science::projectionOf(h,cuts).then("sum_histogram_projection",
      {{"retained_axis",axis},{"windows",cuts},{"flow_bins","excluded"},{"errors","sum_diagonal_variances"}});
  for (int b = 1; b <= axes[axis]->GetNbins() + 1; ++b) out.edges.push_back(axes[axis]->GetBinLowEdge(b));
  out.v.assign(out.edges.size() - 1, 0.0);
  out.e2 = out.v;
  out.measurementState.assign(out.n(),"measured");
  out.covariance=photonjet_science::covarianceOf(h);
  std::vector<std::vector<bool>> keep(dim);
  for (int a = 0; a < dim; ++a) {
    keep[a].assign(axes[a]->GetNbins() + 1, true);
    for (const auto &w : windows)
      if (w.axis == a)
        for (int b = 1; b <= axes[a]->GetNbins(); ++b)
          keep[a][b] = axes[a]->GetBinLowEdge(b) >= w.lo - 1e-9 && axes[a]->GetBinUpEdge(b) <= w.hi + 1e-9;
  }
  const int nx = axes[0]->GetNbins(), ny = dim > 1 ? axes[1]->GetNbins() : 1, nz = dim > 2 ? axes[2]->GetNbins() : 1;
  for (int x = 1; x <= nx; ++x)
    for (int y = 1; y <= ny; ++y)
      for (int z = 1; z <= nz; ++z) {
        const int at[3] = {x, y, z};
        bool in = true;
        for (int a = 0; a < dim && in; ++a) in = keep[a][at[a]];
        if (!in) continue;
        const int bin = h->GetBin(x, dim > 1 ? y : 0, dim > 2 ? z : 0);
        const double e = h->GetBinError(bin);
        if(!std::isfinite(h->GetBinContent(bin)) || !std::isfinite(e))out.measurementState[at[axis]-1]="missing";
        out.v[at[axis] - 1] += h->GetBinContent(bin);
        out.e2[at[axis] - 1] += e * e;
      }
  if (!sources.empty()) {
    PlotJson objects = PlotJson::array();
    for (const auto &source : sources) objects.push_back(source.at("object"));
    const std::string identity = PlotJson({{"objects", objects}, {"retained_axis", axis}, {"windows", cuts}}).dump();
    gInputObservations.push_back({{"key", identity}, {"label", identity}, {"group", objects.dump()},
      {"sources", sources}, {"bins", photonjet_render::physicsBins(out.edges, out.v, out.e2, std::vector<bool>(out.n(), true))}});
  }
  return out;
}

/* The bins times f, their variances times f squared. */
Bins scaled(Bins b, double f) {
  b.covariance={{"state","not_computed"},{"method","scale propagates diagonal variances only"},{"input_covariance",b.covariance},{"factor",photonjet_science::real(f)}};
  b.physics=b.physics.then("scale",{{"factor",f}});
  for (int i = 0; i < b.n(); ++i) b.v[i] *= f, b.e2[i] *= f * f;
  return b;
}
/* The bins divided by their sum: the fraction of the total in each bin (not divided by the bin width). Unchanged when
   the sum is not positive. */
Bins unit(const Bins &b) {
  Bins out=b.sum()>0 ? scaled(b,1/b.sum()) : b;
  out.physics=b.physics.then("unit_bin_sum",{{"sum",b.sum()},{"denominator_treatment","fixed; normalization covariance not calculated"}});
  return out;
}
/* The bins scaled so that their sum is `total`; scaled to zero when their sum is not positive. */
Bins scaledTo(const Bins &b, double total) {
  Bins out=scaled(b,b.sum()>0 ? total/b.sum() : 0.0);
  out.physics=b.physics.then("normalize_to_total",{{"sum",b.sum()},{"target",total}});return out;
}
/* Per unit of the axis: values and errors divided by the bin width. */
Bins perWidth(Bins b) {
  b.covariance={{"state","not_computed"},{"method","bin width normalization propagates diagonal variances only"},{"input_covariance",b.covariance}};
  b.physics=b.physics.then("density_per_bin_width",{{"edge_source","recorded bin edges"}});
  for (int i = 0; i < b.n(); ++i) b.v[i] /= b.width(i), b.e2[i] /= b.width(i) * b.width(i);
  return b;
}
/* The bins inside [lo, hi]. */
Bins within(const Bins &b, double lo, double hi) {
  Bins out;out.physics=b.physics.then("retain_whole_bins",{{"low",std::isfinite(lo)?PlotJson(lo):PlotJson()},
    {"high",std::isfinite(hi)?PlotJson(hi):PlotJson()},{"open_bounds","null means unbounded"}});
  out.covariance=b.covariance;out.seriesRole=b.seriesRole;
  for (int i = 0; i < b.n(); ++i)
    if (b.inside(i, lo, hi)) {
      if (out.edges.empty()) out.edges.push_back(b.edges[i]);
      out.edges.push_back(b.edges[i + 1]), out.v.push_back(b.v[i]), out.e2.push_back(b.e2[i]);
      if(!b.measurementState.empty())out.measurementState.push_back(b.measurementState.at(i));
    }
  return out;
}
/* Groups of k bins (a remainder is dropped). */
Bins rebin(const Bins &b, int k) {
  Bins out;out.physics=b.physics.then("merge_adjacent_bins",{{"group",k},{"remainder","dropped"},{"errors","sum_diagonal_variances"}});
  out.covariance={{"state","not_computed"},{"method","rebin propagates diagonal variances only"},{"input_covariance",b.covariance}};
  const int n = b.n() / k * k;
  for (int i = 0; i <= n; i += k) out.edges.push_back(b.edges[i]);
  for (int i = 0; i < n; i += k) {
    double v = 0, e2 = 0;
    for (int j = i; j < i + k; ++j) v += b.v[j], e2 += b.e2[j];
    out.v.push_back(v), out.e2.push_back(e2);
    std::string state="measured";for(int j=i;j<i+k;++j)if(!b.measurementState.empty() && b.measurementState.at(j)!="measured")state=b.measurementState.at(j);
    out.measurementState.push_back(state);
  }
  return out;
}
/* Whole stored bins merged to the width nearest to `width` (the first bin gives the stored width). */
Bins rebinToWidth(const Bins &b, double width) {
  return b.n() ? rebin(b, std::max(1L, std::lround(width / b.width(0)))) : b;
}
/* The bins of `b` summed into wider ones with the given edges: each is the sum of the whole stored bins inside it. */
Bins rebinnedTo(const Bins &b, const std::vector<double> &edges) {
  Bins out;out.physics=b.physics.then("sum_whole_bins_on_edges",{{"edges",edges},{"errors","sum_diagonal_variances"}});
  out.edges = edges;
  out.v.assign(edges.size() - 1, 0.0);
  out.e2 = out.v;
  for (size_t k = 0; k + 1 < edges.size(); ++k)
    for (int i = 0; i < b.n(); ++i)
      if (b.inside(i, edges[k], edges[k + 1])) out.v[k] += b.v[i], out.e2[k] += b.e2[i];
  return out;
}
/* a + b bin by bin, values and variances; the two share their binning. */
Bins sum(Bins a, const Bins &b) {
  a.physics=photonjet_science::Trace::combine("sum_independent",{a.physics,b.physics});
  for (int i = 0; i < a.n(); ++i) a.v[i] += b.v[i], a.e2[i] += b.e2[i];
  return a;
}
/* a / b per bin (NaN where b is not positive), without errors. */
Bins over(const Bins &a, const Bins &b) {
  Bins out = a;
  out.physics=photonjet_science::Trace::combine("ratio_without_errors",{a.physics,b.physics});
  out.measurementState.assign(a.n(),"undefined");
  for (int i = 0; i < a.n(); ++i) out.v[i] = b.v[i] > 0 ? a.v[i] / b.v[i] : std::nan(""), out.e2[i] = 0;
  for(int i=0;i<out.n();++i)if(std::isfinite(out.v[i]) && std::isfinite(out.e2[i]))out.measurementState[i]="measured";
  return out;
}
/* a / b per bin with the statistical errors of both (NaN where either is not positive). */
Bins ratioOf(const Bins &a, const Bins &b) {
  Bins out = a;
  out.physics=photonjet_science::Trace::combine("ratio_independent_errors",{a.physics,b.physics});
  out.measurementState.assign(a.n(),"undefined");
  for (int i = 0; i < a.n(); ++i) {
    const bool defined = a.v[i] > 0 && b.v[i] > 0;
    out.v[i] = defined ? a.v[i] / b.v[i] : std::nan("");
    out.e2[i] = defined ? out.v[i] * out.v[i] * (a.e2[i] / (a.v[i] * a.v[i]) + b.e2[i] / (b.v[i] * b.v[i])) : 0;
  }
  for(int i=0;i<out.n();++i)if(std::isfinite(out.v[i]) && std::isfinite(out.e2[i]))out.measurementState[i]="measured";
  return out;
}
/* The share a / (a + b) of two separate populations per bin, with the statistical error of their weighted counts (NaN
   where the sum is not positive). */
Bins shareOf(const Bins &a, const Bins &b) {
  Bins out = a;
  out.physics=photonjet_science::Trace::combine("share_independent_weighted_counts",{a.physics,b.physics});
  out.measurementState.assign(a.n(),"undefined");
  for (int i = 0; i < a.n(); ++i) {
    const double total = a.v[i] + b.v[i], share = total > 0 ? a.v[i] / total : std::nan("");
    out.v[i] = share;
    out.e2[i] = total > 0 ? ((1 - share) * (1 - share) * a.e2[i] + share * share * b.e2[i]) / (total * total) : 0;
  }
  for(int i=0;i<out.n();++i)if(std::isfinite(out.v[i]) && std::isfinite(out.e2[i]))out.measurementState[i]="measured";
  return out;
}
/* The bins of a saved one-dimensional result (purity, thresholds, spectra). */
Bins binsOf(const TH1 *h) {
  Bins out = project(h, 0);
  return out;
}
/* A distribution of counts as points: a bin without entries (no content and no error) gets no marker, so no row of
   zeros lies on the axis. */
Bins populated(Bins b) {
  b.hiddenZeroValues=b.v;
  for (int i = 0; i < b.n(); ++i)
    if (b.v[i] == 0 && b.e2[i] == 0) b.v[i] = std::nan("");
  return b;
}
/* True when at least one bin holds a finite value. */
bool anyFinite(const Bins &b) {
  for (double v : b.v) if (std::isfinite(v)) return true;
  return false;
}

/* ---- marks, lines and frames of a single-frame figure ---------------------------------------------------------- */

/* The range of the drawn values inside [x0, x1]: the smallest, the smallest positive and the largest (with errors). */
struct Span {
  double low = std::numeric_limits<double>::infinity(), positive = low, high = -low;
  void add(double v) {
    if (!std::isfinite(v)) return;
    low = std::min(low, v), high = std::max(high, v);
    if (v > 0) positive = std::min(positive, v);
  }
  void add(const Bins &b, bool errors, double x0, double x1) {
    for (int i = 0; i < b.n(); ++i)
      if (b.edges[i + 1] > x0 + 1e-9 && b.edges[i] < x1 - 1e-9) {
        const double error = errors ? std::sqrt(b.e2[i]) : 0, lower = b.v[i] - error;
        add(b.v[i]), add(b.v[i] + error);
        /* A lower error edge sets the log-scale floor only while it is a visible fraction of the value: a bar that
           reaches zero up to rounding (a single weighted entry) runs to the frame instead of stretching it over
           decades. */
        if (lower > 1e-3 * b.v[i]) add(lower);
        else if (std::isfinite(lower)) low = std::min(low, lower);
      }
  }
};

PlotJson scientificBins(const Bins &b) {
  std::vector<std::string> states=b.measurementState;
  std::vector<double> values=b.v;
  if(states.size()!=values.size())states.assign(values.size(),"measured");
  for(std::size_t i=0;i<values.size();++i) {
    if(!std::isfinite(values[i]) && b.hiddenZeroValues.size()==values.size() && b.hiddenZeroValues[i]==0 && b.e2[i]==0)values[i]=0;
    if((!std::isfinite(values[i]) || !std::isfinite(b.e2[i])) && states[i]!="undefined")states[i]="missing";
  }
  return photonjet_science::bins(b.edges,values,b.e2,states,b.covariance);
}
Bins namedBins(Bins b,const std::string& role) { b.seriesRole=role;return b; }
template<class T> T* scientificBinsObject(T* object,const Bins &b,const char* role) {
  return photonjet_science::bind(object,b.seriesRole.empty()?role:b.seriesRole,b.physics,scientificBins(b),b.covariance);
}

/* The values of the bins as a histogram line of the given colour, style and width, without their errors. */
TH1D *line(const Bins &b, int color, int style = 1, int width = 2, const char *role = __builtin_FUNCTION()) {
  auto *h = new TH1D(fresh().c_str(), "", b.n(), b.edges.data());
  h->SetDirectory(nullptr);
  for (int i = 0; i < b.n(); ++i) h->SetBinContent(i + 1, b.v[i]);
  h->SetLineColor(color);
  h->SetLineStyle(style);
  h->SetLineWidth(width);
  return scientificBinsObject(h,b,role);
}

/* Points at the bin centres (shifted by `shift`) with their errors; x errors: half the bin widths. Non-finite values
   are skipped. */
TGraphErrors *markers(const Bins &b, int color, int marker, bool errors = true, bool xerrors = false, double shift = 0,
                      double size = 1.1, const char *role = __builtin_FUNCTION()) {
  auto *g = new TGraphErrors;
  for (int i = 0; i < b.n(); ++i) {
    if (!std::isfinite(b.v[i])) continue;
    const int n = g->GetN();
    g->SetPoint(n, b.center(i) + shift, b.v[i]);
    g->SetPointError(n, xerrors ? 0.5 * b.width(i) : 0.0, errors && std::isfinite(b.e2[i]) ? std::sqrt(b.e2[i]) : 0.0);
  }
  style(g, color, marker, size);
  return scientificBinsObject(g,b,role);
}

/* The statistical uncertainty of a simulated line: a light band of its colour over each bin, without an outline. */
TH1D *uncertaintyBand(const Bins &b, int color, const char *role = __builtin_FUNCTION()) {
  TH1D *h = line(b, kWhite, 1, 1);
  for (int i = 0; i < b.n(); ++i) h->SetBinError(i + 1, std::sqrt(b.e2[i]));
  h->SetFillColorAlpha(color, 0.25);
  h->SetMarkerSize(0);
  return scientificBinsObject(h,b,role);
}

/* The drawn values inside [x0, x1] as extents, for the IAN frame's y range. */
void extents(std::vector<Extent> &out, const Bins &b, bool errors, double x0, double x1) {
  for (int i = 0; i < b.n(); ++i)
    if (b.edges[i + 1] > x0 + 1e-9 && b.edges[i] < x1 - 1e-9 && std::isfinite(b.v[i]))
      out.push_back(
          {b.edges[i], b.edges[i + 1], b.v[i] + (errors && std::isfinite(b.e2[i]) ? std::sqrt(b.e2[i]) : 0.0)});
}

/* A single-panel frame in the IAN's layout (ianFrame) for the drawn values; low: the smallest drawn value (on a log
   scale the smallest positive one). The frame starts at half of `low` on a log scale; on a linear one at zero, or 8%
   beyond a negative `low`. Returns the empty legend for `labels`. */
TLegend *frame(double x0, double x1, const std::vector<Extent> &drawn, double low, bool logy, const char *xtitle,
               const char *ytitle, const std::vector<std::string> &lines, const std::vector<std::string> &labels) {
  double bottom = 0;
  if (logy) bottom = std::isfinite(low) && low > 0 ? low / 2 : 0;
  else if (std::isfinite(low) && low < 0) bottom = 1.08 * low;
  return ianFrame(x0, x1, bottom, drawn, logy, xtitle, ytitle, lines, labels);
}

/* One curve of a single-panel figure and its legend text: a simulated line in its colour, or DATA points (black full
   circles with their error bars; xerrors: with bin-width bars as well). */
struct Curve {
  Bins bins;
  int color = cBlack;
  std::string text;
  bool points = false, xerrors = false;
};
Curve lineCurve(const Bins &b, int color, const std::string &text) {
  Curve c;
  c.bins = b, c.color = color, c.text = text;
  return c;
}
Curve dataCurve(const Bins &b, const std::string &text, bool binBars = false) {
  Curve c;
  c.bins = b, c.text = text, c.points = true, c.xerrors = binBars;
  return c;
}

/* The x-y plane of a three-dimensional histogram at one z bin, detached. zbin is the number of the bin, counted from
   1: bin 2 of a stage axis is the window [1, 2) of a Cut. */
TH2D *plane(const TH1 *h, int zbin) {
  const auto x = edgesOf(h->GetXaxis()), y = edgesOf(h->GetYaxis());
  auto *m = new TH2D(fresh().c_str(), "", x.size() - 1, x.data(), y.size() - 1, y.data());
  m->SetDirectory(nullptr);
  for (int i = 1; i < static_cast<int>(x.size()); ++i)
    for (int j = 1; j < static_cast<int>(y.size()); ++j) m->SetBinContent(i, j, h->GetBinContent(i, j, zbin));
  const auto* z=h->GetZaxis();
  return photonjet_science::bindMap(m,"plane",photonjet_science::projectionOf(h,
      photonjet_science::window(2,z->GetBinLowEdge(zbin),z->GetBinUpEdge(zbin))).then("fixed_z_plane",
      {{"retained_axes",{0,1}}}),false);
}

/* The input note of a figure drawn from simulation alone: the coverage of those packages' simulation samples. The
   DATA coverage, which the general note also counts, does not bear on such a figure. Reference inputs, an unverified
   coverage or a partial simulation population are recorded for the receipt; the canvas carries no such line. */
std::string simulationNote(const std::vector<std::string> &packages) {
  if (gLegacy) return forReceipt("Reference inputs (development)");
  PlotJson coverage;
  if (auto *saved = dynamic_cast<TObjString *>(gQA->Get("input_coverage")))
    coverage = PlotJson::parse(saved->GetString().Data());
  if (!coverage.contains("samples")) return forReceipt("Input coverage unverified");
  bool any = false;
  for (const auto &[key, record] : coverage.at("samples").items()) {
    const auto slash = key.find('/');
    if (slash == std::string::npos || key.substr(slash + 1) == "data") continue;
    if (std::find(packages.begin(), packages.end(), key.substr(0, slash)) == packages.end()) continue;
    any = true;
    if (record.value("status", "UNKNOWN") != "COMPLETE") return forReceipt("Partial simulation population");
  }
  return any ? "" : forReceipt("Input coverage unverified");
}

/* The same for the panels of a figure (panelOrder: package and class). */
std::string simulationNote(const std::vector<std::pair<std::string, std::string>> &order) {
  std::vector<std::string> packages;
  for (const auto &[package, c] : order) packages.push_back(package);
  return simulationNote(packages);
}

/* The top of a linear y axis from 0 that keeps the drawn values clear of the text of its frame: every value stays
   below `fill` of the frame height and below the lower edge of each block it stands under (a block: its left, right
   and lower edge in NDC). */
double topUnderBlocks(const std::vector<Extent> &drawn, double x0, double x1,
                      const std::vector<std::array<double, 3>> &blocks, double fill = 0.95) {
  const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin(), T = gPad->GetTopMargin();
  const double B = gPad->GetBottomMargin();
  const auto ndc = [&](double x) { return L + (x - x0) / (x1 - x0) * (1 - L - R); };
  double top = 0;
  for (const auto &e : drawn) {
    if (!std::isfinite(e.y)) continue;
    double fraction = fill;
    for (const auto &[left, right, bottom] : blocks)
      if (ndc(e.x1) > left && ndc(e.x0) < right) fraction = std::min(fraction, (bottom - B) / (1 - T - B));
    top = std::max(top, e.y / fraction);
  }
  return top > 0 ? top : 1;
}

/* The upper edge of a frame that starts at `low` for which every drawn value inside [xa, xb] stays below `fraction` of
   the frame's height: the room a text block over that range needs. Minus infinity when nothing is drawn there. */
double topBelow(const std::vector<Extent> &drawn, double low, bool logy, double xa, double xb, double fraction) {
  double top = -std::numeric_limits<double>::infinity();
  for (const auto &e : drawn) {
    if (!std::isfinite(e.y) || e.x1 <= xa || e.x0 >= xb || !(e.y > low)) continue;
    const double needed = logy ? std::exp(std::log(low) + (std::log(e.y) - std::log(low)) / fraction)
                               : low + (e.y - low) / fraction;
    top = std::max(top, needed);
  }
  return top;
}

/* The lower edge of a frame that ends at `high` for which every value of `lows` inside [xa, xb] stays above `fraction`
   of the frame's height: the room a text block under the data needs. Plus infinity when nothing is drawn there. */
double bottomAbove(const std::vector<Extent> &lows, double high, bool logy, double xa, double xb, double fraction) {
  double bottom = std::numeric_limits<double>::infinity();
  for (const auto &e : lows) {
    if (!std::isfinite(e.y) || e.x1 <= xa || e.x0 >= xb || !(e.y < high) || (logy && !(e.y > 0))) continue;
    const double needed = logy ? std::exp((std::log(e.y) - fraction * std::log(high)) / (1 - fraction))
                               : (e.y - fraction * high) / (1 - fraction);
    bottom = std::min(bottom, needed);
  }
  return bottom;
}

/* The sPHENIX mark and the lines under it in the current pad at a chosen text size: (x, y) is the start of the mark's
   baseline in NDC, pitch the distance between baselines in pixels. Returns the baseline of the last line. */
double markAndLines(double x, double y, const std::vector<std::string> &lines, double markPx, double linePx,
                    double pitch) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(43);
  t.SetTextSize(markPx);
  t.DrawLatex(x, y, sPhenixInternalLabel());
  t.SetTextSize(linePx);
  for (const auto &text : nonEmpty(lines)) t.DrawLatex(x, y -= pitch / padHeightPx(), text.c_str());
  return y;
}

/* The colour of the light guide lines that mark labelled values across a frame, as in the IAN. */
int guideColor() { return TColor::GetColor("#DDE3EA"); }

/* Light horizontal guide lines across the current frame at the given values, from x0 to x1, as the IAN draws them at
   the labelled fractions. The caller draws the axes again over them (RedrawAxis), so that the tick marks stay on the
   frame. */
void guideLines(const std::vector<double> &levels, double x0, double x1) {
  TLine level;
  level.SetLineColor(guideColor());
  for (double y : levels) level.DrawLine(x0, y, x1, y);
}

/* ---- the shared-axes grid: frames a few pixels apart, text placed in canvas pixels ----------------------------- */

/* The layout of the shared-axes figure being drawn (sharedGrid), all lengths in pixels: the outer margins and the
   space between the frames, the frame sizes, the font sizes, and the panel whose frame is the current one. */
struct GridLayout {
  TCanvas *canvas = nullptr;
  int nx = 1, ny = 1, column = 0, row = 0;                   /* column, row: the current panel, row 0 at the top */
  double left = 0, right = 0, top = 0, bottom = 0, gap = 0;  /* the outer margins and the space between frames */
  double width = 0;                                          /* of every frame */
  std::vector<double> heights;                               /* of the frames of each row */
  double label = 0, title = 0, text = 0;  /* fonts: tick labels; axis and panel titles; legends and further lines */
  double indent = 0, inset = 0;           /* from a frame's side and from its top or bottom to the text inside it */
  double reach = 0;                       /* of the block last written in the current frame, from that inset */
  double canvasWidth() const { return canvas->GetWw(); }
  double canvasHeight() const { return canvas->GetWh(); }
  double frameLeft() const { return left + column * (width + gap); }
  double frameTop() const {
    double y = top;
    for (int r = 0; r < row; ++r) y += heights[r] + gap;
    return y;
  }
  double height() const { return heights[row]; }
};
GridLayout &gridLayout() {
  static GridLayout layout;
  return layout;
}

/* Two numbers of that layout which the drawers use as well: the distance between the baselines of two lines of text
   in a frame, in units of their font size, and the gap between two blocks of text that follow each other in a frame,
   in units of the tick-label size. */
const double kGridPitch = 1.3, kGridBlockGap = 0.3;

/* Starts a figure of nx x ny frames that share their axes, as the multi-panel figures of the IAN: neighbouring frames
   stand a few pixels apart, the y axis is labelled in the left column and the x axis in the bottom row. label: the
   font size of the tick labels, which the other fonts and the margins follow. rows: the relative heights of the rows
   from the top (equal when not given). */
void sharedGrid(TCanvas &cv, int nx, int ny, double label, const std::vector<double> &rows = std::vector<double>()) {
  auto &g = gridLayout();
  g = {};
  g.canvas = &cv, g.nx = nx, g.ny = ny;
  g.label = label, g.title = std::round(1.12 * label), g.text = std::round(0.86 * label);
  g.indent = std::round(0.9 * label), g.inset = std::round(0.5 * label);
  g.left = 4.5 * label, g.right = 0.9 * label, g.top = 0.7 * label, g.bottom = 3.6 * label, g.gap = 0.5 * label;
  g.width = (g.canvasWidth() - g.left - g.right - (nx - 1) * g.gap) / nx;
  const auto share = [&](int r) { return r < static_cast<int>(rows.size()) ? rows[r] : 1.0; };
  const double frames = g.canvasHeight() - g.top - g.bottom - (ny - 1) * g.gap;
  double shares = 0;
  for (int r = 0; r < ny; ++r) shares += share(r);
  for (int r = 0; r < ny; ++r) g.heights.push_back(frames * share(r) / shares);
  cv.cd();
}

/* The columns and the rows of a figure with one panel per class in the IAN's two columns; fewer than three panels
   stand in one row. */
std::pair<int, int> twoColumns(size_t panels) {
  const int nx = panels > 2 ? 2 : panels, ny = (panels + nx - 1) / nx;
  return {nx, ny};
}

/* The frame of panel (column, row) over [x0, x1] x [y0, y1], in a transparent pad over the whole canvas that becomes
   the current pad: every panel then places its text in canvas pixels and no label is cut at a pad edge. Tick labels
   are written on the outer edges only: y in the left column, x in the bottom row. A label at the very left of an
   inner column's x axis, or at the very top of a lower row's y axis, is left out: it would run into the last label of
   the neighbouring frame. */
TH1F *gridFrame(int column, int row, double x0, double y0, double x1, double y1, bool logy = false) {
  auto &g = gridLayout();
  g.column = column, g.row = row, g.reach = 0;
  const double w = g.canvasWidth(), h = g.canvasHeight();
  g.canvas->cd();
  auto *pad = new TPad(fresh().c_str(), "", 0, 0, 1, 1);
  pad->SetFillStyle(4000);
  pad->SetFrameFillStyle(0);
  const double left = g.frameLeft(), top = g.frameTop();
  pad->SetMargin(left / w, 1 - (left + g.width) / w, 1 - (top + g.height()) / h, top / h);
  pad->SetLogy(logy);
  photonjet_science::draw(pad);
  pad->cd();
  TH1F *f = photonjet_science::drawFrame(*pad,x0, y0, x1, y1);
  TAxis *x = f->GetXaxis(), *y = f->GetYaxis();
  x->SetLabelFont(43), y->SetLabelFont(43);
  x->SetLabelSize(row + 1 == g.ny ? g.label : 0), y->SetLabelSize(column == 0 ? g.label : 0);
  x->SetLabelOffset(0.3 * g.label / h), y->SetLabelOffset(0.3 * g.label / w);
  tickLengths(f, 0.45 * g.label);
  x->SetNdivisions(505), y->SetNdivisions(508);
  if (!logy) y->SetDecimals();
  if (column > 0) x->ChangeLabelByValue(x0, -1, 0);
  if (row > 0) y->ChangeLabelByValue(y1, -1, 0);
  return f;
}

/* The axis titles of the current frame where it lies on the outside of the grid: x under a frame of the bottom row
   (centred; right: ending at the frame's right edge, as on a single frame), y at the left of a frame of the left
   column, centred on it. span: one title per axis for the whole figure instead, centred on the grid. A title longer
   than its frame or grid is written smaller, to fit. */
void gridAxisTitles(const char *xtitle, const char *ytitle, bool span = false, bool right = false) {
  const auto &g = gridLayout();
  const double w = g.canvasWidth(), h = g.canvasHeight();
  const double x0 = span ? g.left : g.frameLeft(), x1 = span ? w - g.right : x0 + g.width;
  const double y0 = span ? g.top : g.frameTop(), y1 = span ? h - g.bottom : y0 + g.height();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(43);
  if (*xtitle && (span || g.row + 1 == g.ny)) {
    t.SetTextSize(g.title * std::min(1.0, 0.96 * (x1 - x0) / (textWidth(xtitle, g.title) * w)));
    t.SetTextAlign(right ? 31 : 21);
    t.DrawLatex((right ? x1 : 0.5 * (x0 + x1)) / w, (g.bottom - 1.5 * g.label - g.title) / h, xtitle);
  }
  if (*ytitle && (span || g.column == 0)) {
    t.SetTextSize(g.title * std::min(1.0, 0.96 * (y1 - y0) / (textWidth(ytitle, g.title) * w)));
    t.SetTextAlign(21);
    t.SetTextAngle(90);
    t.DrawLatex(1.3 * g.title / w, 1 - 0.5 * (y0 + y1) / h, ytitle);
  }
}

/* The height of a block of lines of these font sizes written inside a frame from one of its edges (gridLines): the
   inset and the line pitches. */
double gridBlock(const std::vector<double> &sizes) {
  double height = gridLayout().inset;
  for (double px : sizes) height += kGridPitch * px;
  return height;
}

/* A block of lines (text and font size) inside the current frame at one of its corners: at the left or the right
   side (right aligned there), hanging from the top edge or standing on the bottom one; `skip` beyond the inset for a
   block that follows another one (its reach). */
void gridLines(const std::vector<std::pair<std::string, double>> &lines, bool right = false, bool bottom = false,
               double skip = 0, int color = cBlack) {
  auto &g = gridLayout();
  const double w = g.canvasWidth(), h = g.canvasHeight();
  double height = 0;
  for (const auto &written : lines) height += kGridPitch * written.second;
  double y = bottom ? g.frameTop() + g.height() - g.inset - skip - height : g.frameTop() + g.inset + skip;
  TLatex t;
  t.SetNDC();
  t.SetTextFont(43);
  t.SetTextAlign(right ? 31 : 11);
  t.SetTextColor(color);
  for (const auto &[text, px] : lines) {
    t.SetTextSize(px);
    t.DrawLatex((g.frameLeft() + (right ? g.width - g.indent : g.indent)) / w, 1 - (y + 0.97 * px) / h, text.c_str());
    y += kGridPitch * px;
  }
  g.reach = skip + height;
}

/* The width and the height of a legend for `labels` at font size px: a symbol and the widest label, and per row the
   pitch the tallest label needs. */
std::pair<double, double> gridLegendSize(const std::vector<std::string> &labels, double px) {
  double widest = 0, tallest = 0;
  for (const auto &label : labels) {
    TLatex probe(0, 0, label.c_str());
    probe.SetNDC();
    probe.SetTextFont(43);
    probe.SetTextSize(px);
    UInt_t width = 0, height = 0;
    probe.GetBoundingBox(width, height);
    widest = std::max(widest, double(width)), tallest = std::max(tallest, double(height));
  }
  return {2.6 * px + widest, labels.size() * std::max(kGridPitch * px, 1.1 * tallest)};
}

/* An empty legend for `labels` inside the current frame, placed like a block of lines (gridLines) and sized to its
   labels at font size px (0: the size of the lines); a legend wider than the frame is written smaller, to fit. pitch:
   the height of a row in pixels when the rows are to stand at given places (0: the height the labels need). The
   caller adds the entries in the order of `labels` and draws it. */
TLegend *gridLegend(const std::vector<std::string> &labels, bool right = false, bool bottom = false, double skip = 0,
                    double px = 0, double pitch = 0) {
  auto &g = gridLayout();
  const double w = g.canvasWidth(), h = g.canvasHeight(), room = g.width - 2 * g.indent;
  if (!(px > 0)) px = g.text;
  if (const double wide = gridLegendSize(labels, px).first; wide > room) px *= room / wide;
  const auto [width, natural] = gridLegendSize(labels, px);
  const double height = pitch > 0 ? pitch * labels.size() : natural;
  const double x = g.frameLeft() + (right ? g.width - g.indent - width : g.indent);
  const double y = bottom ? g.frameTop() + g.height() - g.inset - skip - height : g.frameTop() + g.inset + skip;
  TLegend *l = ianLegend(x / w, 1 - (y + height) / h, (x + width) / w, 1 - y / h, 2.0 * px / width);
  l->SetTextSize(px);
  g.reach = skip + height;
  return l;
}

/* The upper edge of a y range from `low` that keeps `value` under a block of `block` pixels at the top of the frames
   of a row: the room for the labels written there. */
double gridCeiling(int row, double low, double value, double block, bool logy = false) {
  const double below = 1 - std::min(block / gridLayout().heights[row], 0.6);
  return logy ? std::exp(std::log(low) + (std::log(value) - std::log(low)) / below) : low + (value - low) / below;
}

/* Points with their error bars over bin-width bars. The points may stand `shift` beside the bin centres, for two
   series that share their bins; the bin-width bars stay on the bin edges. Returns the points, for the legend. */
TGraphErrors *barredPoints(const Bins &b, int color, int marker, bool errors, double shift, double size) {
  photonjet_science::draw(markers(b, color, marker, false, true, 0, 0),"Z");
  TGraphErrors *points = markers(b, color, marker, errors, false, shift, size);
  photonjet_science::draw(points,"P");
  return points;
}

/* A grey dashed reference line in the current frame of a shared-axes figure (a threshold, unity), as the IAN draws
   them: between two points given as fractions of the frame's width and height from its lower left corner. */
void gridDashedLine(double x0, double y0, double x1, double y1) {
  const auto &g = gridLayout();
  const double w = g.canvasWidth(), h = g.canvasHeight();
  dashedLine((g.frameLeft() + x0 * g.width) / w, 1 - (g.frameTop() + (1 - y0) * g.height()) / h,
             (g.frameLeft() + x1 * g.width) / w, 1 - (g.frameTop() + (1 - y1) * g.height()) / h, kGray + 2, 2, true);
}

/* A bold title at the top left of the current frame of a shared-axes figure, on a white box: legible over the cells
   of a map. */
void gridBoxedTitle(const std::string &text, double px) {
  const auto &g = gridLayout();
  const double w = g.canvasWidth(), h = g.canvasHeight(), margin = 0.35 * px;
  const std::string bold = "#bf{" + text + "}";
  const double left = g.frameLeft() + g.indent - margin, right = left + textWidth(bold, px) * w + 2 * margin;
  const double top = g.frameTop() + g.inset - margin, bottom = top + kGridPitch * px + margin;
  auto *box = new TPave(left / w, 1 - bottom / h, right / w, 1 - top / h, 0, "NDC");
  box->SetFillColorAlpha(kWhite, 0.85);
  photonjet_science::draw(box);
  gridLines({{bold, px}});
}

/* The colour scale of a map drawn with its palette in the current frame of a shared-axes figure: a bar as tall as
   the frame in the room kept at its right, labelled like the axes, with its title along it at the canvas edge. */
void gridColorScale(TH2 *map, const char *title) {
  const auto &g = gridLayout();
  const double w = g.canvasWidth(), h = g.canvasHeight(), left = g.frameLeft() + g.width + 0.9 * g.label;
  TAxis *z = map->GetZaxis();
  z->SetLabelFont(43), z->SetLabelSize(g.label), z->SetLabelOffset(0.3 * g.label / w), z->SetTickLength(0.012);
  gPad->Update();  /* paints the map, which creates its palette */
  if (auto *bar = dynamic_cast<TPave *>(map->GetListOfFunctions()->FindObject("palette"))) {
    bar->SetX1NDC(left / w), bar->SetX2NDC((left + 1.3 * g.label) / w);
    bar->SetY1NDC(1 - (g.frameTop() + g.height()) / h), bar->SetY2NDC(1 - g.frameTop() / h);
  }
  TLatex t;
  t.SetNDC();
  t.SetTextFont(43);
  t.SetTextSize(g.title);
  t.SetTextAlign(21);
  t.SetTextAngle(90);
  t.DrawLatex(1 - 0.45 * g.title / w, 1 - (g.frameTop() + 0.5 * g.height()) / h, title);
}

/* ---- photon figures -------------------------------------------------------------------------------------------- */

/* Non-collision-background (NCB) classifier score of the photon candidates before the preselection, from 16 GeV to
   the photon maximum of the selection: DATA as points, the photon+jet and the inclusive-jet simulation (every candidate
   of each) as lines with their statistical uncertainty, each as fractions of its total, and the preselection
   threshold at 0.5 as a dashed line over the full height of the frame. One figure for p+p and one for all Au+Au
   classes together, drawn as a pair: a linear y axis from 0 to at least 1.1 in both (a fraction cannot exceed 1),
   opened further only where a drawn value would reach the header or the legend; the header at the top left and the
   legend between the threshold and the peak at 1. */
void ncbScore(const std::string &pkg) {
  const std::string name =
      "photonOnly/" + (pkg == "PP" ? std::string("pp") : "auau_" + spanTag()) + "_ncb_score_data_mc.png";
  TH1 *data = allCandidates(pkg, "data", "ncb");
  if (!data) return missing(name, "no DATA QA");
  /* GeV: qa.windows.ncb_score_gev (from 16, not from the photon minimum) */
  const auto ncbWindow = qaWindow(pkg, "ncb_score_gev", {16, photonSupport(pkg).second});
  const double low = ncbWindow.first, high = ncbWindow.second;
  /* the fractions are taken of the bins from a score of 0 up; the stored axis starts below 0 */
  auto spectrum = [&](const TH1 *h) {
    return unit(within(project(h, 0, {{1, low, high}}), 0, std::numeric_limits<double>::infinity()));
  };
  std::vector<Curve> simulated;
  for (const char *family : {"phosim", "jetsim"})
    if (TH1 *h = allCandidates(pkg, family, "ncb"))
      simulated.push_back(lineCurve(spectrum(h), sampleLook(family).color, sampleLook(family).text));
  const Bins measured = populated(spectrum(data));
  const std::string dataLabel = "DATA", band = "Simulation stat. uncertainty";

  std::vector<std::string> lines = {(pkg == "PP" ? "p+p, " : "Au+Au, ") + energyOf(pkg)}, labels;
  if (const auto classes = classesOf(pkg); pkg != "PP" && !classes.empty())
    lines.push_back(percent(lowOf(classes.front()), highOf(classes.back())) + " centrality");
  lines.push_back(photonWindowLabel(low, high));
  for (const auto &s : simulated) labels.push_back(s.text);
  labels.push_back(dataLabel);
  if (!simulated.empty()) labels.push_back(band);

  /* the score ends at 1 and the last bin of the frame stays empty, as in the IAN; the legend ends at 0.96, left of the
     two bins the distributions peak in. threshold: the NCB preselection threshold, typed here; it is not read from the
     saved selection. */
  TCanvas cv("cv", "", 960, 700);
  photonjet_science::renderSetting(cv,*gQA,"QA.Windows.NcbScore",PlotJson::array({low,high}));
  const double x0 = 0, x1 = 1.02, threshold = 0.5, legendEnd = 0.96;
  photonjet_science::annotation(cv,"ncb/reference_threshold",photonjet_science::scalar(threshold,"score","displayed NCB requirement"),measured.physics.withParameter("renderer/ncb/reference_threshold","selection",threshold));
  const double L = cv.GetLeftMargin(), R = cv.GetRightMargin(), T = cv.GetTopMargin();
  const double w = padWidthPx(), hpx = padHeightPx();
  const auto ndc = [&](double x) { return L + (x - x0) / (x1 - x0) * (1 - L - R); };
  double labelWidth = 0;
  for (const auto &l : labels) labelWidth = std::max(labelWidth, textWidth(l, kLinePx));
  const double symbol = 46 / w, lx2 = ndc(legendEnd) - 6 / w;
  const double lx1 = std::max(lx2 - symbol - labelWidth - 10 / w, ndc(threshold) + 8 / w);
  const double ly2 = 1 - T - 10 / hpx, ly1 = ly2 - labels.size() * kStepPx / hpx;
  const double hx = L + kInsetPx / w, hy = 1 - T - (14 + kMarkPx) / hpx;
  const double hBottom = header(hx, hy, lines, false) - 10 / hpx;

  /* the common range, opened only if a drawn value reached the header or the legend */
  std::vector<Extent> drawn;
  for (const auto &s : simulated) extents(drawn, s.bins, true, x0, x1);
  extents(drawn, measured, true, x0, x1);
  const std::vector<std::array<double, 3>> blocks = {{hx - 6 / w, hx + headerWidth(lines) + 6 / w, hBottom - 6 / hpx},
                                                     {lx1 - 6 / w, lx2 + 6 / w, ly1 - 8 / hpx}};
  const double top = std::max(1.1, topUnderBlocks(drawn, x0, x1, blocks, 1));

  TH1F *f = photonjet_science::drawFrame(*gPad,x0, 0, x1, top);
  f->GetXaxis()->SetTitle("NCB classifier score");
  f->GetYaxis()->SetTitle("Normalized entries");
  f->GetXaxis()->SetNdivisions(506);
  tickLengths(f, 14);
  TLegend *legend = ianLegend(lx1, ly1, lx2, ly2, symbol / (lx2 - lx1));
  for (const auto &s : simulated) photonjet_science::draw(uncertaintyBand(s.bins, s.color),"E2 SAME");
  for (const auto &s : simulated) {
    TH1D *l = line(s.bins, s.color, 1, 3);
    photonjet_science::draw(l,"HIST SAME");
    legend->AddEntry(l, s.text.c_str(), "l");
  }
  dashedLine(threshold, 0, threshold, top, kGray + 2, 2);
  TGraphErrors *g = markers(measured, cBlack, 20, true, false, 0, 0.8);
  photonjet_science::draw(g,"P");
  legend->AddEntry(g, dataLabel.c_str(), "pe");
  if (!simulated.empty())
    legend->AddEntry(uncertaintyBand(simulated.back().bins, simulated.back().color), band.c_str(), "f");
  header(hx, hy, lines);
  photonjet_science::draw(legend);
  savePhysicsQA(cv, name);
}

/* The cluster E_T window of the DATA candidates of a package as the eta-time figures state it: from the lower edge of
   the first populated E_T bin of time_stage at stage 0 (every candidate, the population of etatime), never below the
   saved QA window qa.et_min_gev, to qa.et_max_gev. Au+Au DATA holds no cluster below 12 GeV, p+p DATA starts at 10 GeV.
   The saved window alone for a file without time_stage. */
std::pair<double, double> clusterEtReach(const std::string &pkg) {
  const double low = selectionValue(pkg, "qa", "et_min_gev", 10), high = selectionValue(pkg, "qa", "et_max_gev", 40);
  TH1 *h = allCandidates(pkg, "data", "time_stage");
  if (!h) return {low, high};
  const Bins et = project(h, 1, {{2, 0, 1}});
  for (int i = 0; i < et.n(); ++i)
    if (et.v[i] > 0) return {std::max(low, et.edges[i]), high};
  return {low, high};
}

/* The map of the first two axes of a 3D histogram with its third axis restricted to the whole bins inside [lo, hi]
   (no flow bins); recorded in the receipt of the figure being drawn like a projection. */
TH2 *mapWithin(const TH1 *h, double lo, double hi) {
  gProjectionSteps.push_back({{"histogram", h->GetTitle()}, {"retained_axes", {0, 1}},
                              {"windows", {{{"axis", 2}, {"low", lo}, {"high", hi}, {"interval", "[low,high]"}}}},
                              {"flow_bins", "excluded"}});
  const TAxis *x = h->GetXaxis(), *y = h->GetYaxis(), *z = h->GetZaxis();
  std::vector<double> ex, ey;
  for (int b = 1; b <= x->GetNbins() + 1; ++b) ex.push_back(x->GetBinLowEdge(b));
  for (int b = 1; b <= y->GetNbins() + 1; ++b) ey.push_back(y->GetBinLowEdge(b));
  auto *map = new TH2D(fresh().c_str(), "", ex.size() - 1, ex.data(), ey.size() - 1, ey.data());
  map->SetDirectory(nullptr);
  map->Sumw2();
  for (int k = 1; k <= z->GetNbins(); ++k) {
    if (z->GetBinLowEdge(k) < lo - 1e-9 || z->GetBinUpEdge(k) > hi + 1e-9) continue;
    for (int i = 1; i <= x->GetNbins(); ++i)
      for (int j = 1; j <= y->GetNbins(); ++j) {
        const int bin = h->GetBin(i, j, k);
        const double e = h->GetBinError(bin), was = map->GetBinError(i, j);
        map->SetBinContent(i, j, map->GetBinContent(i, j) + h->GetBinContent(bin));
        map->SetBinError(i, j, std::sqrt(was * was + e * e));
      }
  }
  const double low=std::isfinite(lo)?lo:z->GetXmin(),high=std::isfinite(hi)?hi:z->GetXmax();
  return photonjet_science::bindMap(map,"map_within",photonjet_science::projectionOf(h,
      photonjet_science::window(2,low,high)).then("sum_z_map",{{"retained_axes",{0,1}},{"variance","sum_diagonal"}}),true);
}

/* The smallest positive cell of a map; infinity for an empty map. */
double lowestPositive(const TH1 *map) {
  double low = std::numeric_limits<double>::infinity();
  for (int c = 0; c < map->GetNcells(); ++c) if (map->GetBinContent(c) > 0) low = std::min(low, map->GetBinContent(c));
  return low;
}

/* The header lines of the eta-time figures: the system with its centrality range and energy, and the cluster E_T
   window the DATA reach (clusterEtReach). */
std::vector<std::string> etaTimeLines(const std::string &pkg) {
  const bool auau = pkg == "AuAu";
  const auto classes = classesOf(pkg);
  const std::string system = auau && !classes.empty()
                                 ? "Au+Au " + centralityRange(lowOf(classes.front()), highOf(classes.back()))
                                 : auau ? "Au+Au" : "p+p";
  const auto et = clusterEtReach(pkg);
  return {system + ", " + energyOf(pkg), Form("%g #leq E_{T}^{cluster} < %g GeV", et.first, et.second)};
}

/* One eta-time map on its own canvas, saved as `name`: drawn over its whole stored range, the time difference and the
   eta the QA file holds, on the canonical kBird logarithmic colour scale, from `zLow` to `zHigh`
   (zHigh <= 0: the map's own maximum); the header (the mark, the lines given) stands at the top left on two white
   boxes, one under the mark and one under the lines and the input note, each no wider than its text, so that little
   of the map lies under them. */
void drawEtaTimeMap(TH2 *map, const std::vector<std::string> &lines, const std::string &name, double zLow, double zHigh) {
  if(photonjet_science::tag(map,"physics_series").is_null())photonjet_science::bindMap(map,"eta_time",
      photonjet_science::projectionOf(map,PlotJson::array()),true);
  TCanvas cv("cv", "", 820, 640);
  cv.SetRightMargin(0.18);
  cv.SetLogz();
  map->GetXaxis()->SetTitle("Cluster #eta");
  map->GetYaxis()->SetTitle("Cluster time #minus MBD time [ns]");
  map->GetZaxis()->SetTitle("Clusters");
  map->GetZaxis()->CenterTitle();
  map->GetZaxis()->SetTitleOffset(1.2);
  if (std::isfinite(zLow)) map->SetMinimum(zLow);
  if (zHigh > 0) map->SetMaximum(zHigh);
  photonjet_science::draw(map,"COLZ");
  const double w = padWidthPx(), hpx = padHeightPx();
  const double x = cv.GetLeftMargin() + kInsetPx / w, y = 1 - cv.GetTopMargin() - (14 + kMarkPx) / hpx;
  const auto box = [&](double width, double bottom, double top) {
    auto *b = new TPave(x - 5 / w, bottom, x + width + 5 / w, top, 0, "NDC");
    b->SetFillColorAlpha(kWhite, 0.9);
    b->SetLineWidth(0);
    photonjet_science::draw(b);
  };
  std::vector<std::string> below = lines;
  below.push_back(inputNote());
  double wide = 0;
  for (const auto &l : nonEmpty(below)) wide = std::max(wide, textWidth(l, kLinePx));
  box(textWidth(sPhenixInternalLabel(), kMarkPx), y - 9 / hpx, y + (kMarkPx + 1) / hpx);
  box(wide, header(x, y, lines, false) - 9 / hpx, y - 9 / hpx);
  mapHeader(lines);
  savePhysicsQA(cv, name);
}

/* Cluster time minus MBD time against cluster eta for the DATA candidates before the NCB requirement, one map per
   system: p+p (pp_ncb_eta_time_before) and Au+Au summed over its classes (auau_ncb_eta_time_before, panel 2 of the
   figure, Justin 2026-10-14), drawn by drawEtaTimeMap on the map's own colour range. The E_T line states the window the
   DATA reach (12 GeV and up in Au+Au; Justin 2026-10-14, "go ahead and fix that"). */
void ncbEtaTime(const std::string &pkg);
void ncbEtaTime() {
  ncbEtaTime("PP");
  ncbEtaTime("AuAu");
}
void ncbEtaTime(const std::string &pkg) {
  const bool auau = pkg == "AuAu";
  const std::string name = auau ? "photonOnly/auau_ncb_eta_time_before.png" : "photonOnly/pp_ncb_eta_time_before.png";
  TH1 *h = allCandidates(pkg, "data", "etatime");
  if (!h) return missing(name, auau ? "no Au+Au DATA QA" : "no p+p DATA QA");
  auto *map = static_cast<TH2 *>(h);
  drawEtaTimeMap(map, etaTimeLines(pkg), name, lowestPositive(map), 0);
}

/* The IAN's second NCB TODO (Justin 2026-10-14): the eta-time map of the same DATA candidates after the NCB requirement,
   one per system (pp_ncb_eta_time_after, auau_ncb_eta_time_after), on the colour range of the map before it
   (ncbEtaTime), so that the two compare cell for cell; a line says the share of the candidates kept. From etatime_npb
   (the NCB score as third axis): the NCB requirement alone, s_NCB > 0.5, the threshold typed as in ncbScore. A file
   without it but with etatime_stage gets the map after the preselection (stage 1), which holds the NCB requirement and
   the three shower-shape requirements, and says so on its own line. */
void ncbEtaTimeAfter(const std::string &pkg) {
  const bool auau = pkg == "AuAu";
  const std::string name = auau ? "photonOnly/auau_ncb_eta_time_after.png" : "photonOnly/pp_ncb_eta_time_after.png";
  TH1 *before = allCandidates(pkg, "data", "etatime");
  if (!before) return missing(name, auau ? "no Au+Au DATA QA" : "no p+p DATA QA");
  const double threshold = 0.5, inf = std::numeric_limits<double>::infinity();
  TH2 *all = nullptr, *after = nullptr;
  std::string step;
  if (TH1 *byScore = allCandidates(pkg, "data", "etatime_npb")) {
    all = mapWithin(byScore, -inf, inf), after = mapWithin(byScore, threshold, inf);
    step = Form("After the NCB requirement, s_{NCB} > %g", threshold);
  } else if (TH1 *byStage = allCandidates(pkg, "data", "etatime_stage")) {
    all = mapWithin(byStage, 0, 1), after = mapWithin(byStage, 1, 2);
    step = "After the preselection";
  } else {
    return missing(name, "needs the cluster time against the NCB score (etatime_npb) in the QA fill");
  }
  auto lines = etaTimeLines(pkg);
  lines.push_back(step);
  if (all->Integral() > 0) lines.push_back(Form("%.0f%% of the clusters kept", 100 * after->Integral() / all->Integral()));
  if(all->Integral()>0)photonjet_science::attachAnnotation(after,"retained_cluster_fraction",
      photonjet_science::scalar(100*after->Integral()/all->Integral(),"percent","sum_after / sum_before"),
      photonjet_science::derivedTrace("retained_cluster_fraction",{after,all}));
  drawEtaTimeMap(after, lines, name, lowestPositive(before), before->GetMaximum());
}

/* The IAN's first NCB TODO (Justin 2026-10-14): the cluster time minus MBD time of the DATA candidates for different
   requirements on the NCB score, none, s_NCB > 0.3, the nominal 0.5 and 0.7, one figure per system
   (pp_ncb_time_selections, auau_ncb_time_selections), as clusters per ns on a logarithmic axis over the whole stored
   time range; the legend gives each requirement with the share of the candidates it keeps; the header names the
   vertex window, the cluster acceptance and that no other selection is applied. Needs etatime_npb (the NCB score as
   third axis), which the fill adds; not drawn without it. */
void ncbTimeSelections(const std::string &pkg) {
  const bool auau = pkg == "AuAu";
  const std::string name = auau ? "photonOnly/auau_ncb_time_selections.png" : "photonOnly/pp_ncb_time_selections.png";
  TH1 *h = allCandidates(pkg, "data", "etatime_npb");
  if (!h) return missing(name, "needs the cluster time against the NCB score (etatime_npb) in the QA fill");
  const double inf = std::numeric_limits<double>::infinity();
  struct Requirement { double cut; int color, width; std::string text; };
  const std::vector<Requirement> requirements = {{-inf, cBlack, 2, "No NCB requirement"}, {0.3, cBlue, 2, "s_{NCB} > 0.3"},
                                                 {0.5, cRed, 3, "s_{NCB} > 0.5 (nominal)"}, {0.7, cGreen, 2, "s_{NCB} > 0.7"}};
  std::vector<Bins> shown;
  std::vector<std::string> labels;
  std::vector<double> retained;
  double total = 0, most = 0, least = inf;
  for (const auto &r : requirements) {
    const Bins counts = project(h, 1, {{2, r.cut, inf}});
    if (shown.empty()) total = counts.sum();
    const Bins b = perWidth(counts);
    for (int i = 0; i < b.n(); ++i)
      if (b.v[i] > 0) most = std::max(most, b.v[i]), least = std::min(least, b.v[i]);
    retained.push_back(total>0?100*counts.sum()/total:std::nan(""));
    shown.push_back(b);
    shown.back().seriesRole="ncb_requirement/"+std::to_string(shown.size()-1);
    labels.push_back(shown.size() == 1 || total <= 0 ? r.text : Form("%s, %.0f%% kept", r.text.c_str(), 100 * counts.sum() / total));
  }
  if (!(most > 0)) return missing(name, "no DATA candidates in etatime_npb");
  /* the cuts the clusters pass, from the saved selection (Justin, Oct 15: name the cuts, not "DATA") */
  auto lines = etaTimeLines(pkg);
  lines.push_back(Form("%s < %g cm, %s < %g", absolute("z_{vtx}").c_str(),
                       selectionValue(pkg, "event", "vertex_abs_z_max_cm", pkg == "PP" ? 60 : 10),
                       absolute("#eta^{cluster}").c_str(), selectionValue(pkg, "photon", "abs_eta_max", 0.7)));
  lines.push_back("Before the preselection");

  TCanvas cv("cv", "", 960, 700);
  for(std::size_t i=1;i<shown.size();++i)if(std::isfinite(retained[i]))
    photonjet_science::annotation(cv,"ncb/retained/"+std::to_string(i),photonjet_science::scalar(retained[i],"percent","selected sum / unselected sum"),
      photonjet_science::Trace::combine("retained_fraction",{shown[i].physics,shown[0].physics}));
  cv.SetLogy();
  const double L = cv.GetLeftMargin(), R = cv.GetRightMargin(), T = cv.GetTopMargin();
  const double w = padWidthPx(), hpx = padHeightPx();
  double labelWidth = 0;
  for (const auto &l : labels) labelWidth = std::max(labelWidth, textWidth(l, kLinePx));
  const double symbol = 46 / w, lx2 = 1 - R - 10 / w, lx1 = lx2 - symbol - labelWidth - 10 / w;
  const double ly2 = 1 - T - 10 / hpx, ly1 = ly2 - labels.size() * kStepPx / hpx;
  const double hx = L + kInsetPx / w, hy = 1 - T - (14 + kMarkPx) / hpx;
  /* the curves fill the lower 78% of the frame in the logarithm; the header and the legend keep the two top corners,
     clear of the narrow central peak (60% left three empty decades over the p+p peak, Justin 2026-10-15) */
  const double bottom = 0.5 * least, top = bottom * std::pow(most / bottom, 1 / 0.78);
  TH1F *f = photonjet_science::drawFrame(*gPad,shown.front().edges.front(), bottom, shown.front().edges.back(), top);
  f->GetXaxis()->SetTitle("Cluster time #minus MBD time [ns]");
  f->GetYaxis()->SetTitle("Clusters / ns");
  tickLengths(f, 14);
  TLegend *legend = ianLegend(lx1, ly1, lx2, ly2, symbol / (lx2 - lx1));
  for (size_t k = 0; k < shown.size(); ++k) {
    TH1D *l = line(shown[k], requirements[k].color, 1, requirements[k].width);
    photonjet_science::draw(l,"HIST SAME");
    legend->AddEntry(l, labels[k].c_str(), "l");
  }
  header(hx, hy, lines);
  photonjet_science::draw(legend);
  savePhysicsQA(cv, name);
}

/* Shower shapes and the photon-ID score after the preselection, 15-20 GeV: DATA as points against two simulated
   lines, the candidates of the photon+jet simulation matched to a signal photon ("sig") and the candidates of the
   inclusive-jet simulation not matched to one ("bkg"), each as fractions of its total over the whole stored range of
   the variable; one figure per variable and Au+Au class (dis_<shape>_cent<k>_pt0_cut1). Each variable is shown in the
   IAN's bin width, whole stored bins merged (the fractions E_T,2 and E_T,3 are stored in bins of 0.02, so three make
   0.06 where the IAN shows 0.05), and a bin without DATA entries gets no point. The window is the IAN's, except where
   the preselection itself ends the distribution: E32/E35 starts at 0.8 and the full-cluster width in eta ends at 0.6,
   so E32/E35 is shown from 0.75 and the two full-cluster widths, drawn as a pair, from 0 to 1. Header and legend are
   sized for the note, where nine of these figures stand three to a row: the legend at the top right, the mark and two
   short lines beside it, and the photon selection on a line of its own under them. */
void showerShapes() {
  /* the file-name slot, the stored variable, the axis title, the window shown and the stored bins merged into one */
  struct Shape { const char *slot, *var, *xtitle; double low, high; int merged; };
  const std::vector<Shape> shapes = {
      {"e11_to_e33", "e11_over_e33", "E_{11}#scale[0.5]{ }/E_{33}", 0, 1, 5},
      {"e32_to_e35", "e32_over_e35", "E_{32}#scale[0.5]{ }/E_{35}", 0.75, 1, 1},
      {"bdt", "bdt", "Photon-ID score s_{BDT}", 0, 1, 2}, {"weta_cogx", "weta_cogx", "w_{#eta}^{COG}", 0, 1, 2},
      {"wphi_cogx", "wphi_cogx", "w_{#phi}^{COG}", 0, 1, 2}, {"et1", "et1", "E_{T,1} fraction", 0.5, 1, 2},
      {"et2", "et2", "E_{T,2} fraction", -0.5, 1, 3}, {"et3", "et3", "E_{T,3} fraction", -0.5, 1, 3},
      {"et4", "et4", "E_{T,4} fraction", 0, 0.25, 2}};
  const auto classes = classesOf("AuAu");
  for (size_t k = 0; k < classes.size(); ++k)
    for (const auto &s : shapes) {
      const std::string name = Form("shower_shapes/dis_%s_cent%zu_pt0_cut1.pdf", s.slot, k), c = classes[k];
      TH1 *data = qa("AuAu/data/photons/" + c + "/data/shape_" + s.var);
      if (!data) { missing(name, "no Au+Au DATA QA"); continue; }
      /* GeV: the photon window of these figures, qa.windows.shower_shapes_gev */
      const auto shapeWindow = qaWindow("AuAu", "shower_shapes_gev", {15, 20});
      const double low = shapeWindow.first, high = shapeWindow.second;
      const Cut afterPreselection = {2, 1, 2};  /* of the stage axis */
      auto spectrum = [&](const TH1 *h) {
        return within(rebin(unit(project(h, 0, {{1, low, high}, afterPreselection})), s.merged), s.low, s.high);
      };
      std::vector<Curve> simulated;
      /* the legend in the words of the photon-ID figures (Justin, Oct 15: "label it 'Background-tagged: Inclusive+jet
         12/20/30/40' ... no where does it specify embedded simulation samples"): the embedded samples by their
         generator thresholds in GeV, typed here as in the BDT-score figure; the background line is bound to its
         candidate population when the legend is drawn (candidateEntry) */
      for (const auto &[family, population] : {std::make_pair("phosim", "sig"), std::make_pair("jetsim", "bkg")})
        if (TH1 *h = qa(std::string("AuAu/") + family + "/photons/" + c + "/" + population + "/shape_" + s.var))
          simulated.push_back(lineCurve(spectrum(h), sampleLook(family).color,
                                        std::string(family) == "phosim" ? "Prompt photon: photon+jet 12/20"
                                                                        : "Background-tagged: inclusive+jet 12/20/30/40"));
      const Bins measured = populated(spectrum(data));

      /* a square canvas, as in the IAN, on which the header and legend lines stand at about three quarters of the
         size of the tick labels; its axis text is sized so that the y title stays clear of the tick labels, its ticks
         so that they end short of the header */
      TCanvas cv("cv", "", 720, 720);
      photonjet_science::renderSetting(cv,*gQA,"QA.Windows.ShowerShapes",PlotJson::array({shapeWindow.first,shapeWindow.second}));
      std::vector<Extent> drawn;
      std::vector<std::string> labels = {"DATA"};
      extents(drawn, measured, true, s.low, s.high);
      for (const auto &m : simulated) extents(drawn, m.bins, false, s.low, s.high), labels.push_back(m.text);
      const std::string eta = photonEtaLine("AuAu");
      const std::vector<std::string> lines = {
          "Au+Au, " + energyOf("AuAu"), percent(lowOf(c), highOf(c)) + ", after preselection",
          photonWindowLabel(low, high) + (eta.empty() ? "" : ", " + eta), "Embedded PYTHIA 8 simulation"};
      const double L = cv.GetLeftMargin(), R = cv.GetRightMargin(), T = cv.GetTopMargin();
      const double w = padWidthPx(), hpx = padHeightPx();
      double labelWidth = 0, wide = headerWidth(lines);
      for (const auto &l : labels) labelWidth = std::max(labelWidth, textWidth(l, kLinePx));
      if (!inputNote().empty()) wide = std::max(wide, textWidth(inputNote(), kLinePx));
      const double symbol = 46 / w, hx = L + kInsetPx / w, hy = 1 - T - (14 + kMarkPx) / hpx;
      const double hBottom = header(hx, hy, lines, false) - 10 / hpx;
      double lx2 = 1 - R - 14 / w, lx1 = lx2 - symbol - labelWidth - 10 / w, ly2 = 1 - T - 10 / hpx;
      if (lx1 < hx + wide + 20 / w)  /* the legend does not fit beside the header: it stands under it */
        lx1 = hx - 8 / w, lx2 = lx1 + symbol + labelWidth + 10 / w, ly2 = hBottom - 4 / hpx;
      const double ly1 = ly2 - labels.size() * kStepPx / hpx;
      const double top = topUnderBlocks(drawn, s.low, s.high, {{hx - 6 / w, hx + wide + 6 / w, hBottom - 10 / hpx},
                                                              {lx1 - 6 / w, lx2 + 6 / w, ly1 - 10 / hpx}});
      TH1F *f = photonjet_science::drawFrame(*gPad,s.low, 0, s.high, top);
      f->GetXaxis()->SetTitle(s.xtitle);
      f->GetYaxis()->SetTitle("Normalized counts");
      for (TAxis *a : {f->GetXaxis(), f->GetYaxis()})
        a->SetLabelSize(0.042), a->SetTitleSize(0.044), a->SetTickLength(0.022);
      f->GetYaxis()->SetTitleOffset(1.7);
      f->GetYaxis()->SetNoExponent();
      f->GetXaxis()->SetNdivisions(505);
      header(hx, hy, lines);
      TLegend *legend = ianLegend(lx1, ly1, lx2, ly2, symbol / (lx2 - lx1));
      TGraphErrors *g = markers(measured, cBlack, 20, true, false, 0, 1.4);
      legend->AddEntry(g, labels.front().c_str(), "pe");
      for (const auto &m : simulated) {
        TH1D *l = line(m.bins, m.color);
        photonjet_science::draw(l,"HIST SAME");
        if (m.text.rfind("Background-tagged", 0) == 0)
          candidateEntry(legend, l, photonjet_render::CandidatePopulation::TruthTaggedBackground, "inclusive+jet", "l",
                         {12, 20, 30, 40});
        else
          legend->AddEntry(l, m.text.c_str(), "l");
      }
      photonjet_science::draw(g,"P");
      photonjet_science::draw(legend);
      savePhysicsQA(cv, name);
    }
}

/* ---- Figure 10: shower width against centrality (widths and the helpers beside it) -----------------------------

   The figure shows that the full-cluster widths of Au+Au photon candidates grow strongly toward central collisions,
   where the underlying event spreads the cluster, while the 3x3 widths, the two widths that the 14-feature Au+Au
   photon-ID model adds, change much less. Each centrality class keeps one colour in every view (widthClassColor):
   its simulation is a line of that colour and its DATA are points of that colour. The files, all under photonOnly/:
     auau_signal_<full_cluster|3x3>_width.png        the IAN's figure: the simulated candidates of the photon+jet
                                                      samples matched to a truth photon, preselected, one line per class
     auau_background_<full_cluster|3x3>_width.png    the same for the background candidates of the inclusive-jet
                                                      samples, each line with its statistical band
     auau_data_simulation_<full_cluster|3x3>_width.png   DATA over simulation, one frame per class, tight and isolated
                                                      on both sides; a supplement of the signal figure of that width
     auau_width_vs_centrality.png                    the mean width against centrality, full cluster and 3x3 in both
                                                      directions, DATA points over simulated lines; a supplement of
                                                      the signal full-cluster figure
   Every distribution is drawn as fractions of its total. The stored bins are 0.01 wide; each view merges whole stored
   bins to the width its statistics carry as a smooth curve: the truth-matched photons hold a few hundred thousand
   effective candidates per class, the inclusive-jet background, after its sample weights, about four thousand. DATA
   are compared like for like, tight and isolated against tight and isolated simulated photons, and are not purity
   corrected. The full-cluster widths of DATA lie below those of the simulation in every class; that known difference
   is shown as it is, and no population is changed to close it. */

/* The colour of class k of the width figures, from the most central: blue, orange (the palette's vermilion), green. */
int widthClassColor(size_t k) {
  const int colors[] = {cBlue, cOrange, cGreen, cRed, cPurple};
  return colors[k % 5];
}

/* The windows of a slices/ width histogram for a centrality range in percent: the photon E_T window of the package
   and the whole 5% slices of the range. */
std::vector<Cut> widthCuts(double low, double high) {
  const auto domain = photonSupport("AuAu");
  return {{1, domain.first, domain.second}, {2, low / kSlicePercent, high / kSlicePercent}};
}

/* The embedded Au+Au samples of a simulation family as the working-point figures name them ("Photon+jet 12/20",
   "Inclusive jet 12/20/30/40"): the number in the name of each sample the Au+Au package of the QA file records in
   that family (photonjet, inclusivejet); the bare family name without sample records. */
std::string auauFamilySamples(const std::string &family) {
  const std::string name = family == "photonjet" ? "Photon+jet" : "Inclusive jet";
  if (gQAMetadata.is_null() || !gQAMetadata.at("packages").contains("AuAu")) return name;
  std::set<int> thresholds;
  const auto samples = gQAMetadata.at("packages").at("AuAu").value("samples", PlotJson::object());
  for (const auto &sample : samples.items()) {
    if (!sample.value().is_object() || sample.value().value("family", "") != family) continue;
    const auto digits = sample.key().find_first_of("0123456789");
    if (digits != std::string::npos) thresholds.insert(std::stoi(sample.key().substr(digits)));
  }
  std::string list;
  for (int t : thresholds) list += (list.empty() ? " " : "/") + std::to_string(t);
  return name + list;
}

/* The mean of a distribution of counts, each bin at its centre, and the statistical error of that mean from the
   variances of the bins; NaN for a distribution without entries. */
std::pair<double, double> meanWidth(const Bins &b) {
  double total = 0, first = 0;
  for (int i = 0; i < b.n(); ++i) total += b.v[i], first += b.v[i] * b.center(i);
  if (!(total > 0)) return {std::nan(""), std::nan("")};
  const double mean = first / total;
  double variance = 0;
  for (int i = 0; i < b.n(); ++i) variance += b.e2[i] * (b.center(i) - mean) * (b.center(i) - mean);
  return {mean, std::sqrt(variance) / total};
}

/* True when the figure at `path` was written in this run: its receipt stands in the output folder. */
bool receiptWritten(const std::string &path) { return std::filesystem::exists(gFig + "/" + path + ".json"); }

/* One of the IAN's width figures on a single frame: the simulated candidates of `path` (a slices/ width histogram
   after the preselection), one line per class as fractions of its total in bins of `bin` (whole stored bins merged)
   over [0, 0.7]; band: each line with its statistical uncertainty as a light band of its colour. The header keeps
   the IAN's lines at the top left, the centrality legend stands at the top right. */
void widthClassFigure(const std::string &name, const std::string &path, const char *xtitle, const char *text,
                      double bin, bool band) {
  TH1 *simulated = qa(path);
  const auto classes = classesOf("AuAu");
  if (!simulated || classes.empty()) return missing(name, "needs the Au+Au shower-width QA of the simulation");
  const auto domain = photonSupport("AuAu");
  const double x0 = 0, x1 = 0.7;
  double width = bin;
  std::vector<Bins> lines;
  std::vector<Extent> drawn;
  std::vector<std::string> labels = {"Centrality"};
  for (const auto &c : classes) {
    lines.push_back(rebinToWidth(unit(project(simulated, 0, widthCuts(lowOf(c), highOf(c)))), bin));
    lines.back().seriesRole="width_class/"+std::to_string(lines.size()-1);
    width = lines.back().width(0);
    extents(drawn, lines.back(), band, x0, x1);
    labels.push_back(percent(lowOf(c), highOf(c)));
  }

  TCanvas cv("cv", "", 960, 640);
  const std::vector<std::string> headerLines = {
      "Au+Au, " + energyOf("AuAu") + ", " + photonWindowLabel(domain.first, domain.second), text};
  const std::string ytitle = Form("Unit-normalized candidates / %g", width);
  const double L = cv.GetLeftMargin(), R = cv.GetRightMargin(), T = cv.GetTopMargin();
  const double w = padWidthPx(), hpx = padHeightPx();
  /* the widths of the widest legend line (its heading included) and of the header */
  double labelWidth = 0, wide = headerWidth(headerLines);
  for (const auto &l : labels) labelWidth = std::max(labelWidth, textWidth(l, kLinePx));
  if (!inputNote().empty()) wide = std::max(wide, textWidth(inputNote(), kLinePx));
  const double symbol = 46 / w, lx2 = 1 - R - 14 / w, lx1 = lx2 - symbol - labelWidth - 10 / w;
  const double ly2 = 1 - T - 10 / hpx, ly1 = ly2 - labels.size() * kStepPx / hpx;
  const double hx = L + kInsetPx / w, hy = 1 - T - (14 + kMarkPx) / hpx;
  const double hBottom = header(hx, hy, headerLines, false) - 10 / hpx;
  const std::vector<std::array<double, 3>> blocks = {{hx - 6 / w, hx + wide + 6 / w, hBottom - 6 / hpx},
                                                     {lx1 - 6 / w, lx2 + 6 / w, ly1 - 8 / hpx}};
  TH1F *f = photonjet_science::drawFrame(*gPad,x0, 0, x1, topUnderBlocks(drawn, x0, x1, blocks));
  f->GetXaxis()->SetTitle(xtitle);
  f->GetYaxis()->SetTitle(ytitle.c_str());
  f->GetYaxis()->SetNoExponent();
  f->GetYaxis()->SetNdivisions(506);  /* a handful of labelled values, as in the IAN */
  tickLengths(f, 14);
  header(hx, hy, headerLines);
  TLegend *legend = ianLegend(lx1, ly1, lx2, ly2, symbol / (lx2 - lx1));
  legend->SetHeader(labels.front().c_str(), "C");
  if (band)
    for (size_t i = 0; i < classes.size(); ++i) photonjet_science::draw(uncertaintyBand(lines[i], widthClassColor(i)),"E2 SAME");
  for (size_t i = 0; i < classes.size(); ++i) {
    TH1D *l = line(lines[i], widthClassColor(i), 1, 3);
    photonjet_science::draw(l,"HIST SAME");
    legend->AddEntry(l, labels[i + 1].c_str(), "l");
  }
  photonjet_science::draw(legend);
  gPad->RedrawAxis();
  savePhysicsQA(cv, name, "widths");
}

/* DATA over simulation for one width, one frame per class side by side on shared axes, each in the colour of its
   class: the tight and isolated simulated candidates matched to a truth photon as a line with its statistical band,
   the tight and isolated DATA candidates (not purity corrected) as points, each as fractions of its total in bins of
   `bin` (whole stored bins merged) over [0, high]. A bin without DATA entries gets no point. The class stands at the
   top right of its frame with the mean width of DATA and of the simulation under it (over the stored range, without
   flow bins; at this stage no candidate lies beyond it); the mark and the selection stand at the top left of the
   first frame, the key of the line and the points at the top left of the second. One frame per class keeps the DATA
   of a class on its own line: on one frame the narrower full-cluster DATA of a class would stand on the simulated
   line of the next class. */
void widthDataPanels(const std::string &name, const std::string &var, const char *xtitle, double bin, double high) {
  PlotPhysicsScope physics(name, "qa.width_data_panels");
  auto &record = physics.record;
  const auto real = photonjet_render::physicsDouble;
  record.definition("projection", {{"retained_axis", 0}, {"flow_bins", "excluded"},
    {"bin_window_rule", "whole stored bins within low/high, edge epsilon 1e-9"},
    {"normalization", "unit sum before rebin and visible-range restriction"},
    {"rebin", "whole bins nearest requested width"}, {"requested_width", real(bin)},
    {"visible_low", real(0)}, {"visible_high", real(high)},
    {"data_marker_rule", "omit zero content and zero variance bins"},
    {"means", "weighted bin centers over full stored range excluding flow"}});

  TH1 *simulated = qa("AuAu/phosim/slices/sig/width_tight_iso_" + var);
  TH1 *data = qa("AuAu/data/slices/data/width_tight_iso_" + var);
  const auto classes = classesOf("AuAu");
  if (!simulated || !data || classes.empty())
    return missing(name, "needs the tight and isolated Au+Au shower-width QA of DATA and of the simulation");
  const auto domain = photonSupport("AuAu");
  struct Panel { Bins line, points; double dataMean, simulatedMean; };
  std::vector<Panel> panels;
  Span s;
  double width = bin;
  for (const auto &c : classes) {
    const Bins fineSimulated = project(simulated, 0, widthCuts(lowOf(c), highOf(c)));
    const Bins fineData = project(data, 0, widthCuts(lowOf(c), highOf(c)));
    Panel p;
    p.line = within(rebinToWidth(unit(fineSimulated), bin), 0, high);
    const Bins denseData = within(rebinToWidth(unit(fineData), bin), 0, high);
    p.points = populated(denseData);
    p.line.seriesRole="width/"+std::to_string(panels.size())+"/simulation";
    p.points.seriesRole="width/"+std::to_string(panels.size())+"/data";
    p.dataMean = meanWidth(fineData).first;
    p.simulatedMean = meanWidth(fineSimulated).first;
    const std::string key = c;
    std::vector<bool> all(p.line.n(), true), dataMask;
    for (int i=0; i<denseData.n(); ++i) dataMask.push_back(denseData.v[i] != 0 || denseData.e2[i] != 0);
    record.quantity(key + "/simulation", photonjet_render::physicsBins(p.line.edges,p.line.v,p.line.e2,all));
    record.quantity(key + "/data", photonjet_render::physicsBins(denseData.edges,denseData.v,denseData.e2,dataMask));
    record.quantity(key + "/means", {{"data",real(p.dataMean)},{"simulation",real(p.simulatedMean)}});
    record.definition(key + "/window", {{"centrality_low",real(lowOf(c))},
      {"centrality_high",real(highOf(c))}, {"photon_et_low",real(domain.first)},
      {"photon_et_high",real(domain.second)}});

    if (p.line.n()) width = p.line.width(0);
    s.add(p.line, true, 0, high), s.add(p.points, true, 0, high);
    panels.push_back(p);
  }

  const int n = classes.size();
  TCanvas cv("cv", "", 560 * n, 560);
  sharedGrid(cv, n, 1, 22);
  const auto &g = gridLayout();
  using Lines = std::vector<std::pair<std::string, double>>;
  const double classSize = g.title + 5, systemSize = 21.5, lineSize = 19.5;
  const double between = kGridBlockGap * g.label, clear = 0.8 * g.label;
  const Lines system = {{sPhenixInternalLabel(), g.title + 1}, {"Au+Au, " + energyOf("AuAu"), systemSize}};
  const Lines selection = {{gridWindowLine(domain.first, domain.second), lineSize},
                           {"Tight and isolated candidates", lineSize},
                           {"Embedded " + auauFamilySamples("photonjet") + " simulation", lineSize}};
  const std::vector<std::string> labels = {"Simulation, signal photons", "DATA"};
  const int keyPanel = n > 1 ? 1 : 0;  /* with one class the key stands under the selection */
  const auto means = [&](const Panel &p) {
    return Lines{{Form("Mean: DATA %.3f", p.dataMean), lineSize}, {Form("simulation %.3f", p.simulatedMean), lineSize}};
  };
  /* the y range opens until the highest bin stands clear of the tallest block of text */
  const double head = gridBlock({g.title + 1, systemSize, lineSize, lineSize}) + between;
  const double key = gridLegendSize(labels, lineSize).second;
  double tallest = gridBlock({classSize, lineSize, lineSize}) + between;
  tallest = std::max(tallest, keyPanel == 0 ? head + between + key : std::max(head, g.inset + key));
  const double ceiling = gridCeiling(0, 0, s.high, tallest + clear);
  for (int k = 0; k < n; ++k) {
    const int color = widthClassColor(k);
    TH1F *axes = gridFrame(k, 0, 0, 0, high, ceiling);
    axes->GetXaxis()->SetNdivisions(505);
    gridAxisTitles(xtitle, Form("Unit-normalized candidates / %g", width), false, true);
    photonjet_science::draw(uncertaintyBand(panels[k].line, color),"E2 SAME");
    photonjet_science::draw(line(panels[k].line, color, 1, 3),"HIST SAME");
    photonjet_science::draw(markers(panels[k].points, color, 20, true, false, 0, 1.2),"P");
    gridLines({{"#bf{" + percent(lowOf(classes[k]), highOf(classes[k])) + "}", classSize}}, true, false, 0, color);
    gridLines(means(panels[k]), true, false, g.reach + between);
    photonjet_science::annotation(*gPad,"mean_width/data/"+std::to_string(k),photonjet_science::scalar(panels[k].dataMean,"1","weighted bin centres over full stored range"),panels[k].points.physics);
    photonjet_science::annotation(*gPad,"mean_width/simulation/"+std::to_string(k),photonjet_science::scalar(panels[k].simulatedMean,"1","weighted bin centres over full stored range"),panels[k].line.physics);
    double skip = 0;
    if (k == 0) {
      gridLines(system);
      gridLines(selection, false, false, g.reach + between);
      skip = g.reach + between;
    }
    if (k == keyPanel) {
      TLegend *l = gridLegend(labels, false, false, skip, lineSize);
      auto *lineKey = new TGraph;
      lineKey->SetLineColor(cBlack), lineKey->SetLineWidth(3);
      l->AddEntry(lineKey, labels[0].c_str(), "l");
      auto *pointKey = new TGraphErrors(1);
      style(pointKey, cBlack, 20, 1.2);
      l->AddEntry(pointKey, labels[1].c_str(), "pe");
      photonjet_science::draw(l);
    }
    gPad->RedrawAxis();
  }
  record.dependency("qa_histogram_keys", {
    {"simulation", "qa/AuAu/phosim/slices/sig/width_tight_iso_" + var},
    {"data", "qa/AuAu/data/slices/data/width_tight_iso_" + var}});
  // Inputs are provenance, not physics identity. Full per-histogram definition
  // ancestry is not available yet: do NOT mark selection/models/weights complete.
  record.dependency("input_products", gRenderInputs);
  for (const char *field : {"quantities","errors","scalar_annotations","normalization","projection"})
    record.coverage(field,"complete","widthDataPanels: final Bins arrays, validity masks, exact means and operations");
  savePhysicsQA(cv, name, "widths");
}

/* The summary of the width figures: the mean width against centrality, one frame per direction (eta, phi) on shared
   axes. Each class is cut into pieces of 10% (two stored slices; a class narrower than that is one piece). Per piece
   the full-cluster width (solid line, full circles) and the 3x3 width (dashed line, open circles): the tight and
   isolated simulated candidates matched to a truth photon as a step line, the tight and isolated DATA (not purity
   corrected) as points with their statistical errors, both in the colour of the class. The means are taken over the
   stored range without flow bins; at this stage no candidate lies beyond it. Light guide lines mark the edges between
   the classes, which are named in their colours under the points. The mark and the selection stand at the top left of
   the first frame, the direction at the top right of each, the key at the top right of the last. The receipt records
   every mean drawn (width_means). */
void widthSummary(const std::string &name) {
  struct Direction { const char *full, *small, *title; };
  const std::vector<Direction> directions = {{"weta_cogx", "weta33_cogx", "#eta direction"},
                                             {"wphi_cogx", "wphi33_cogx", "#phi direction"}};
  const auto classes = classesOf("AuAu");
  if (classes.empty()) return missing(name, "needs the Au+Au centrality classes");
  const double step = 10;  /* percent */
  struct Piece { size_t cls; double low, high; };
  std::vector<Piece> pieces;
  for (size_t k = 0; k < classes.size(); ++k)
    for (double a = lowOf(classes[k]); a < highOf(classes[k]) - 1e-9; a += step)
      pieces.push_back({k, a, std::min(a + step, highOf(classes[k]))});
  /* per direction and width (0 the full cluster, 1 the 3x3 window) and per sample (0 simulation, 1 DATA): the mean
     and its error of every piece */
  using Means = std::vector<std::pair<double, double>>;
  std::vector<std::array<std::array<Means, 2>, 2>> means(directions.size());
  std::vector<std::array<std::array<std::vector<photonjet_science::Trace>,2>,2>> meanTraces(directions.size());
  PlotJson record = PlotJson::array();
  Span s;
  for (size_t d = 0; d < directions.size(); ++d)
    for (int z = 0; z < 2; ++z) {
      const std::string var = z == 0 ? directions[d].full : directions[d].small;
      TH1 *simulated = qa("AuAu/phosim/slices/sig/width_tight_iso_" + var);
      TH1 *data = qa("AuAu/data/slices/data/width_tight_iso_" + var);
      if (!simulated || !data)
        return missing(name, "needs the tight and isolated Au+Au shower-width QA of DATA and of the simulation");
      for (const auto &p : pieces) {
        const auto cuts = widthCuts(p.low, p.high);
        const auto simBins=project(simulated,0,cuts),dataBins=project(data,0,cuts);
        const auto simulatedMean = meanWidth(simBins), dataMean = meanWidth(dataBins);
        meanTraces[d][z][0].push_back(simBins.physics);meanTraces[d][z][1].push_back(dataBins.physics);
        means[d][z][0].push_back(simulatedMean);
        means[d][z][1].push_back(dataMean);
        s.add(simulatedMean.first), s.add(dataMean.first + dataMean.second);
        record.push_back(PlotJson{{"variable", var},
                                  {"centrality_percent", {p.low, p.high}},
                                  {"simulation_mean", simulatedMean.first},
                                  {"simulation_error", simulatedMean.second},
                                  {"data_mean", dataMean.first},
                                  {"data_error", dataMean.second}});
      }
    }

  const auto domain = photonSupport("AuAu");
  const int n = directions.size();
  TCanvas cv("cv", "", 600 * n, 600);
  sharedGrid(cv, n, 1, 22);
  const auto &g = gridLayout();
  using Lines = std::vector<std::pair<std::string, double>>;
  const double titleSize = g.title + 5, systemSize = 21.5, lineSize = 19.5;
  const double between = kGridBlockGap * g.label, clear = 0.8 * g.label;
  const Lines system = {{sPhenixInternalLabel(), g.title + 1}, {"Au+Au, " + energyOf("AuAu"), systemSize}};
  const Lines selection = {{gridWindowLine(domain.first, domain.second), lineSize},
                           {"Tight and isolated candidates", lineSize}};
  const std::vector<std::string> labels = {"Full cluster", "3#times3"};
  const Lines key = {{"Points: DATA", lineSize},
                     {"Lines: embedded " + auauFamilySamples("photonjet") + " sim, signal photons", lineSize}};
  /* the y range opens until the highest mean stands clear of the tallest block of text */
  const double head = gridBlock({g.title + 1, systemSize, lineSize, lineSize}) + between;
  const double right = gridBlock({titleSize}) + between + gridLegendSize(labels, lineSize).second + between +
                       kGridPitch * lineSize * key.size();
  const double ceiling = gridCeiling(0, 0, s.high, std::max(head, right) + clear);
  const double x0 = lowOf(classes.front()), x1 = highOf(classes.back());
  for (int d = 0; d < n; ++d) {
    TH1F *axes = gridFrame(d, 0, x0, 0, x1, ceiling);
    axes->GetXaxis()->SetNdivisions(508);
    gridAxisTitles("Centrality [%]", "Mean shower width", false, true);
    TLine guide;
    guide.SetLineColor(guideColor());
    for (size_t k = 1; k < classes.size(); ++k) guide.DrawLine(lowOf(classes[k]), 0, lowOf(classes[k]), ceiling);
    for (int z = 0; z < 2; ++z)
      for (size_t k = 0; k < classes.size(); ++k) {
        auto *steps = new TGraph;
        auto *points = new TGraphErrors;
        for (size_t i = 0; i < pieces.size(); ++i) {
          if (pieces[i].cls != k) continue;
          const auto [simulatedMean, simulatedError] = means[d][z][0][i];
          const auto [dataMean, dataError] = means[d][z][1][i];
          if (std::isfinite(simulatedMean)) {
            steps->SetPoint(steps->GetN(), pieces[i].low, simulatedMean);
            steps->SetPoint(steps->GetN(), pieces[i].high, simulatedMean);
          }
          if (std::isfinite(dataMean)) {
            const int at = points->GetN();
            points->SetPoint(at, 0.5 * (pieces[i].low + pieces[i].high), dataMean);
            points->SetPointError(at, 0, std::isfinite(dataError) ? dataError : 0.0);
          }
        }
        for(int sample=0;sample<2;++sample) {
          std::vector<photonjet_science::Trace> traces;PlotJson values=PlotJson::array();
          for(std::size_t i=0;i<pieces.size();++i)if(pieces[i].cls==k) {
            traces.push_back(meanTraces[d][z][sample][i]);const auto mean=means[d][z][sample][i];
            values.push_back({{"low",pieces[i].low},{"high",pieces[i].high},{"mean",std::isfinite(mean.first)?PlotJson(mean.first):PlotJson()},
              {"error",std::isfinite(mean.second)?PlotJson(mean.second):PlotJson()},{"state",std::isfinite(mean.first)&&std::isfinite(mean.second)?"measured":"missing"}});
          }
          photonjet_science::bind(sample?static_cast<TGraph*>(points):steps,"width_mean/"+std::to_string(d)+"/"+std::to_string(z)+"/"+std::to_string(k)+"/"+std::to_string(sample),
            photonjet_science::Trace::combine("weighted_bin_center_means",traces,{{"error","sum variance*(center-mean)^2 / total^2"}}),
            photonjet_science::portable(values),photonjet_science::covariance("not_computed",{},0,"no covariance between disjoint centrality pieces"));
        }
        const int color = widthClassColor(k);
        steps->SetLineColor(color), steps->SetLineWidth(3), steps->SetLineStyle(z == 0 ? 1 : 7);
        if (steps->GetN()) photonjet_science::draw(steps,"L");
        style(points, color, z == 0 ? 20 : 71, 1.5);  /* 71: the open circle with lines three wide */
        if (points->GetN()) photonjet_science::draw(points,"P");
      }
    TLatex t;  /* the classes, named in their colours under the points */
    t.SetTextFont(43), t.SetTextSize(lineSize), t.SetTextAlign(21);
    for (size_t k = 0; k < classes.size(); ++k) {
      t.SetTextColor(widthClassColor(k));
      t.DrawLatex(0.5 * (lowOf(classes[k]) + highOf(classes[k])), 0.04 * ceiling,
                  percent(lowOf(classes[k]), highOf(classes[k])).c_str());
    }
    gridLines({{directions[d].title, titleSize}}, true);
    if (d == n - 1) {
      TLegend *l = gridLegend(labels, true, false, g.reach + between, lineSize);
      auto *fullKey = new TGraph, *smallKey = new TGraph;
      style(fullKey, cBlack, 20, 1.5), fullKey->SetLineWidth(3);
      style(smallKey, cBlack, 71, 1.5), smallKey->SetLineWidth(3), smallKey->SetLineStyle(7);
      l->AddEntry(fullKey, labels[0].c_str(), "lp");
      l->AddEntry(smallKey, labels[1].c_str(), "lp");
      photonjet_science::draw(l);
      gridLines(key, true, false, g.reach + between);
    }
    if (d == 0) {
      gridLines(system);
      gridLines(selection, false, false, g.reach + between);
    }
    gPad->RedrawAxis();
  }
  savePhysicsQA(cv, name, "widths");
  if (receiptWritten(name)) amendReceipt(name, [&](PlotJson &receipt) { receipt["width_means"] = record; });
}

/* Figure 10: the IAN's four width figures, the two DATA figures and the summary (see the head of this group). */
void widths() {
  struct Variable { const char *var, *slot, *xtitle; };
  const std::vector<Variable> variables = {{"weta_cogx", "full_cluster", "w_{#eta}^{COG}"},
                                           {"weta33_cogx", "3x3", "w_{#eta}^{3#times3, COG}"}};
  /* a view of the IAN's figure: the word for it in the file name, the family and the folder under slices/ of its
     simulated candidates, the header line that names it, the bin width for each variable and whether its lines carry
     their statistical band */
  struct View { const char *who, *family, *label; std::string text; double bins[2]; bool band; };
  const std::vector<View> views = {
      {"signal", "phosim", "sig", "Embedded " + auauFamilySamples("photonjet") + " sim, signal photons, preselected",
       {0.02, 0.01}, false},
      {"background", "jetsim", "bkg", "Embedded " + auauFamilySamples("inclusivejet") + " sim, background, preselected",
       {0.04, 0.02}, true}};
  /* the DATA figures: the bin width and the upper end of the axis for each variable */
  const double dataBins[] = {0.02, 0.01}, dataHigh[] = {0.5, 0.3};
  const auto fileOf = [](const char *who, const char *slot) {
    return std::string("photonOnly/auau_") + who + "_" + slot + "_width.png";
  };
  const std::string summary = "photonOnly/auau_width_vs_centrality.png";
  for (const auto &view : views)
    for (size_t i = 0; i < variables.size(); ++i)
      widthClassFigure(fileOf(view.who, variables[i].slot),
                       std::string("AuAu/") + view.family + "/slices/" + view.label + "/width_" + variables[i].var,
                       variables[i].xtitle, view.text.c_str(), view.bins[i], view.band);
  for (size_t i = 0; i < variables.size(); ++i)
    widthDataPanels(fileOf("data_simulation", variables[i].slot), variables[i].var, variables[i].xtitle, dataBins[i],
                    dataHigh[i]);
  widthSummary(summary);

  /* the DATA figures and the summary are supplements of the signal figures: the receipts of both sides say so */
  for (size_t i = 0; i < variables.size(); ++i) {
    const std::string signal = fileOf("signal", variables[i].slot), overlay = fileOf("data_simulation", variables[i].slot);
    const std::string var = variables[i].var;
    PlotJson supplements = PlotJson::array();
    if (receiptWritten(overlay)) {
      supplements.push_back(PlotJson{{"id", "data_simulation"},
                                     {"title", "DATA over simulated photons, tight and isolated, one frame per class"},
                                     {"output", overlay}, {"function", "widths"}, {"position", "after"}});
      amendReceipt(overlay, [&](PlotJson &receipt) {
        receipt["plot_recipe"] =
            plotRecipe("auau_width_centrality", "data_simulation", "supplemental", overlay, signal, "widths");
        receipt["data_overlay"] = {
            {"data", "AuAu/data/slices/data/width_tight_iso_" + var + ": tight and isolated Au+Au DATA candidates, "
                     "not purity corrected"},
            {"simulation", "AuAu/phosim/slices/sig/width_tight_iso_" + var + ": tight and isolated simulated "
                           "candidates matched to an analysis-signal truth photon"},
            {"normalisation", "each distribution as fractions of its total"},
            {"means", "over the stored range, flow bins excluded"}};
      });
    }
    if (i == 0 && receiptWritten(summary)) {
      supplements.push_back(PlotJson{{"id", "width_vs_centrality"},
                                     {"title", "Mean width against centrality: full cluster and 3x3, DATA and simulation"},
                                     {"output", summary}, {"function", "widths"}, {"position", "after"}});
      amendReceipt(summary, [&](PlotJson &receipt) {
        receipt["plot_recipe"] =
            plotRecipe("auau_width_centrality", "width_vs_centrality", "supplemental", summary, signal, "widths");
      });
    }
    if (!supplements.empty() && receiptWritten(signal))
      amendReceipt(signal, [&](PlotJson &receipt) { receipt["supplemental_outputs"] = supplements; });
  }
}

/* The photon-ID score of the Au+Au simulation after the preselection, without an isolation requirement, one panel per
   centrality class on shared axes as in the IAN's score figure. When the QA file holds the saved score record of its
   default model (model_qa/<model>/application/review), savedScoreReview draws this figure from that record, with a
   second figure for the truth-tagged background, and nothing below runs. Otherwise the figure is drawn from the slice
   histograms: the candidates of the photon+jet simulation matched to an analysis-signal photon ("sig") against the
   candidates of the inclusive-jet simulation not matched to one ("bkg"), each curve as fractions of its total and
   shown in bins of 0.02. The inclusive-jet curve holds every candidate that is not matched to an analysis-signal
   photon; it is not a pure truth-background tag. The labels stand at the top of the frames, which the y range is
   opened for: the class in every panel, under it the mark and the selection in the first panel, the legend in the
   second and the sample in the last.
   When the fill states its photon-ID lines (photon_id.optional_model.working_point_choice, mode "explicit"), each
   panel also shows the regions of the photon ID at the middle of its class: non-tight from non_tight_lower to the
   non-tight upper line (dashed edges), excluded up to the tight line, tight above it (solid), named and without counts.
   The tight, isolated DATA of the class are drawn as points when the package holds them
   (AuAu/data/photons/<class>/data/scoreiso_fine, score in 0.01 bins, when the fill has it, else scoreiso, in 0.05
   bins: score, isolation, E_T; Justin, Oct 14: three points in 50-80% are too few): isolated means the isolation bins whose centre
   lies below the isolation line at the middle of the class, the E_T window that of the curves; the points are scaled
   to the prompt curve over the DATA bins wholly above the tight cut and shown per the curves' bin width. */
void bdtPerClass() {
  if (savedScoreReview()) return;
  const std::string name = "photonOnly/auau_bdt_score_perCent.png";
  TH1 *sig = qa("AuAu/phosim/slices/sig/bdt"), *bkg = qa("AuAu/jetsim/slices/bkg/bdt");
  const auto classes = classesOf("AuAu");
  if (!sig || !bkg || classes.empty()) return missing(name, "needs Au+Au photon and inclusive simulation QA");
  std::vector<std::pair<Bins, Bins>> spectra;
  const auto domain = photonSupport("AuAu");
  Span s;
  double width = 0.02;
  for (const auto &c : classes) {
    const std::vector<Cut> cuts = {{1, domain.first, domain.second},
                                   {2, lowOf(c) / kSlicePercent, highOf(c) / kSlicePercent}};
    auto shape = [&](const TH1 *h) { return rebinToWidth(unit(project(h, 0, cuts)), 0.02); };
    spectra.push_back({shape(sig), shape(bkg)});
    width = spectra.back().first.width(0);
    s.add(spectra.back().first, false, 0, 1), s.add(spectra.back().second, false, 0, 1);
  }
  /* the regions of the photon ID and the DATA: read before the layout, which makes room for the DATA note */
  const auto &selected = renderSelection("AuAu", true);
  struct Lines2 { double a = std::nan(""), b = std::nan(""); double at(double c) const { return a + b * c; } };
  Lines2 tightLine, upperLine, isoLine;
  double lowerEdge = std::nan("");
  auto lineOf = [](const PlotJson &r, Lines2 &l) {
    if (r.is_object() && r.contains("intercept") && r.contains("slope"))
      l.a = r.at("intercept").get<double>(), l.b = r.at("slope").get<double>();
  };
  if (!selected.is_null() && selected.contains("photon_id")) {
    const auto &id = selected.at("photon_id");
    if (id.contains("optional_model") && id.at("optional_model").is_object() &&
        id.at("optional_model").contains("working_point_choice")) {
      const auto &choice = id.at("optional_model").at("working_point_choice");
      if (choice.is_object() && choice.value("mode", std::string()) == "explicit") {
        if (choice.contains("tight")) lineOf(choice.at("tight"), tightLine);
        if (choice.contains("non_tight_upper")) lineOf(choice.at("non_tight_upper"), upperLine);
        if (choice.contains("non_tight_lower") && choice.at("non_tight_lower").is_number())
          lowerEdge = choice.at("non_tight_lower").get<double>();
      }
    }
    if (id.contains("isolation") && id.at("isolation").contains("threshold")) lineOf(id.at("isolation").at("threshold"), isoLine);
  }
  const bool regions = std::isfinite(tightLine.a) && std::isfinite(upperLine.a) && std::isfinite(lowerEdge);
  /* per class: the DATA points (centre, value, error) in the units of the curves, or none */
  std::vector<std::vector<std::array<double, 3>>> dataPoints(classes.size());
  std::vector<photonjet_science::Trace> dataPointTrace(classes.size());
  std::string overlayStem = "scoreiso";   /* the histogram the points came from, for the receipt */
  for (size_t k = 0; regions && std::isfinite(isoLine.a) && k < classes.size(); ++k) {
    const double mid = (lowOf(classes[k]) + highOf(classes[k])) / 2, tight = tightLine.at(mid);
    TH1 *scoreiso = qa("AuAu/data/photons/" + classes[k] + "/data/scoreiso_fine");
    if (scoreiso) overlayStem = "scoreiso_fine";
    else scoreiso = qa("AuAu/data/photons/" + classes[k] + "/data/scoreiso");
    if (!scoreiso || scoreiso->GetDimension() != 3) continue;
    const double isoHigh = std::floor(isoLine.at(mid) + 0.5);  /* whole 1 GeV bins whose centre is below the line */
    const Bins d = project(scoreiso, 0, {{1, scoreiso->GetYaxis()->GetXmin(), isoHigh}, {2, domain.first, domain.second}});
    /* the prompt curve above the same edge, from the fine projection */
    const std::vector<Cut> cuts = {{1, domain.first, domain.second},
                                   {2, lowOf(classes[k]) / kSlicePercent, highOf(classes[k]) / kSlicePercent}};
    const Bins prompt = unit(project(sig, 0, cuts));
    int first = -1;
    for (int i = 0; i < d.n() && first < 0; ++i)
      if (d.edges[i] >= tight - 1e-9) first = i;
    if (first < 0) continue;
    double inData = 0, inPrompt = 0;
    for (int i = first; i < d.n(); ++i) inData += d.v[i];
    for (int i = 0; i < prompt.n(); ++i)
      if (prompt.edges[i] >= d.edges[first] - 1e-9) inPrompt += prompt.v[i];
    if (!(inData > 0 && inPrompt > 0)) continue;
    dataPointTrace[k]=photonjet_science::Trace::combine("data_scaled_to_prompt_above_tight",{d.physics,prompt.physics},
      {{"tight",tight},{"isolation_high",isoHigh},{"data_sum",inData},{"prompt_sum",inPrompt},{"curve_bin_width",width}});
    for (int i = first; i < d.n(); ++i) {
      const double per = width / d.width(i);  /* shown per the curves' bin width */
      dataPoints[k].push_back({d.center(i), d.v[i] * inPrompt / inData * per, std::sqrt(d.e2[i]) * inPrompt / inData * per});
      s.add(d.v[i] * inPrompt / inData * per + std::sqrt(d.e2[i]) * inPrompt / inData * per);
    }
  }
  bool anyData = false;
  for (const auto &p : dataPoints) anyData = anyData || !p.empty();
  const int n = classes.size();
  TCanvas cv("cv", "", 520 * n, 560);
  const std::string simulation = simulationNote({"AuAu"});
  sharedGrid(cv, n, 1, 22);
  const auto &g = gridLayout();
  using Lines = std::vector<std::pair<std::string, double>>;
  /* the text sizes: half pixels, which set a symbol and its indices tight (see kLinePx) */
  const double classSize = g.title + 5, systemSize = 21.5, lineSize = 19.5;
  const double between = kGridBlockGap * g.label, clear = 0.8 * g.label;
  const int lineWidth = 3;
  const double vertex = selectionValue("AuAu", "event", "vertex_abs_z_max_cm", 10);
  const Lines system = {{sPhenixInternalLabel(), g.title + 1}, {"Au+Au, " + energyOf("AuAu"), systemSize}};
  const Lines selection = {{gridWindowLine(domain.first, domain.second), lineSize},
                           {photonEtaLine(selectionValue("AuAu", "photon", "abs_eta_max", 0.7)) +
                                Form(", |%s| < %g cm", apart("z_{vtx}").c_str(), vertex), lineSize}};
  Lines sample = {{"Embedded simulation", lineSize},
                  {"Preselected, no isolation requirement", lineSize},
                  {"Each curve unit-normalized", lineSize}};
  if (anyData) sample.push_back({"DATA: scaled to the prompt curve above the tight cut", lineSize});
  const auto inclusive = photonjet_render::CandidatePopulation::OtherInclusiveJetCandidates;
  std::vector<std::string> labels = {
      "Prompt photons", photonjet_render::candidatePopulation(inclusive).at("label").get<std::string>()};
  if (anyData) labels.push_back("DATA: tight, isolated");
  /* the y range opens until the highest bin stands clear of the tallest block of labels */
  const int legendPanel = std::min(1, n - 1);
  double tallest = 0;
  for (int k = 0; k < n; ++k) {
    std::vector<double> sizes = {classSize};
    double further = 0;
    if (k == 0) {
      for (const Lines &block : {system, selection})
        for (const auto &written : block) sizes.push_back(written.second);
      further += between;
    }
    if (k == legendPanel) further += between + gridLegendSize(labels, lineSize).second;
    if (k == n - 1) {
      sizes.insert(sizes.end(), sample.size() + (simulation.empty() ? 0 : 1), lineSize);
      further += between;
    }
    tallest = std::max(tallest, gridBlock(sizes) + further);
  }
  const double ceiling = gridCeiling(0, 0, s.high, tallest + clear);
  for (int k = 0; k < n; ++k) {
    TH1F *axes = gridFrame(k, 0, 0, 0, 1, ceiling);
    axes->GetXaxis()->SetNdivisions(405);  /* labelled every 0.2 as in the IAN, a tick every 0.05 */
    axes->GetXaxis()->SetDecimals();
    /* the symbol is set apart, so that its index sits tight to it at the title size (apart) */
    gridAxisTitles(("BDT score, " + apart("s_{BDT}")).c_str(), Form("Weighted fraction of candidates / %g", width));
    if (regions) {
      /* the regions at the middle of the class, under the curves and up to just above the highest of them */
      const double mid = (lowOf(classes[k]) + highOf(classes[k])) / 2, top = 1.06 * s.high;
      const double low = lowerEdge, upper = upperLine.at(mid), tight = tightLine.at(mid);
      auto band = [](double x1, double x2, double y, int color, double alpha) {
        auto *b = new TBox(x1, 0, x2, y);
        b->SetFillColorAlpha(color, alpha);
        b->SetLineColorAlpha(color, 0);
        photonjet_science::draw(b);
      };
      band(low, upper, top, TColor::GetColor("#e78300"), .11);
      band(upper, tight, top, kGray, .18);
      band(tight, 1, top, TColor::GetColor("#2d71df"), .10);
      for (auto [x, style, color] : {std::tuple<double, int, int>{low, kDashed, TColor::GetColor("#925505")},
                                     {upper, kDashed, TColor::GetColor("#78609b")},
                                     {tight, 1, TColor::GetColor("#234d80")}}) {
        auto *l = new TLine(x, 0, x, top);
        l->SetLineStyle(style), l->SetLineColor(color), l->SetLineWidth(3);
        photonjet_science::draw(l);
      }
      auto name = [&](double x, double y, const char *text, int color, double angle) {
        auto *t = new TLatex(x, y, text);
        t->SetTextFont(63), t->SetTextSize(lineSize), t->SetTextColor(color), t->SetTextAlign(angle ? 12 : 21);
        t->SetTextAngle(angle);
        photonjet_science::draw(t);
      };
      name((low + upper) / 2, 0.86 * top, "Non-tight", TColor::GetColor("#925505"), 0);
      name((tight + 1) / 2, 0.86 * top, "Tight", TColor::GetColor("#234d80"), 0);
      name((upper + tight) / 2, 0.30 * top, "Excluded", kGray + 2, 90);
    }
    TH1D *prompt = line(spectra[k].first, cBlue, 1, lineWidth);
    TH1D *others = line(spectra[k].second, cOrange, 1, lineWidth);
    photonjet_science::draw(others,"HIST SAME");
    photonjet_science::draw(prompt,"HIST SAME");
    TGraphErrors *points = nullptr;
    if (!dataPoints[k].empty()) {
      points = new TGraphErrors;
      for (const auto &[x, y, e] : dataPoints[k]) {
        points->SetPoint(points->GetN(), x, y);
        points->SetPointError(points->GetN() - 1, 0, e);
      }
      photonjet_science::bindNative(points,"bdt_data/"+std::to_string(k),dataPointTrace[k],
        photonjet_science::covariance("not_computed",{},0,"diagonal DATA bars; normalization covariance not calculated"));
      points->SetMarkerStyle(20), points->SetMarkerSize(1.25), points->SetMarkerColor(cBlack);
      points->SetLineColor(cBlack), points->SetLineWidth(2);
      photonjet_science::draw(points,"P SAME");
    }
    gridLines({{"#bf{" + percent(lowOf(classes[k]), highOf(classes[k])) + "}", classSize}});
    if (k == 0) {
      gridLines(system, false, false, g.reach);
      gridLines(selection, false, false, g.reach + between);
    }
    if (k == legendPanel) {
      TLegend *l = gridLegend(labels, false, false, g.reach + between, lineSize);
      l->AddEntry(prompt, labels[0].c_str(), "l");
      candidateEntry(l, others, inclusive, "inclusive+jet");
      if (anyData) {
        auto *key = new TGraphErrors(1);
        key->SetMarkerStyle(20), key->SetMarkerSize(1.25), key->SetMarkerColor(cBlack), key->SetLineColor(cBlack);
        l->AddEntry(key, labels.back().c_str(), "pe");
      }
      photonjet_science::draw(l);
    }
    if (k == n - 1) {
      gridLines(sample, false, false, g.reach + between, kGray + 3);
      if (!simulation.empty()) gridLines({{simulation, lineSize}}, false, false, g.reach, kGray + 2);
    }
  }
  savePhysicsQA(cv, name);
  if (regions)
    amendReceipt(name, [&](PlotJson &receipt) {
      receipt["photon_id_regions"] = {
        {"source", "photon_id.optional_model.working_point_choice (mode explicit)"},
        {"tight", {{"intercept", tightLine.a}, {"slope", tightLine.b}}},
        {"non_tight_upper", {{"intercept", upperLine.a}, {"slope", upperLine.b}}},
        {"non_tight_lower", lowerEdge}, {"evaluated_at", "middle of each centrality class"}, {"counts", false}};
      if (anyData)
        receipt["data_overlay"] = {
          {"histograms", "AuAu/data/photons/<class>/data/" + overlayStem},
          {"population", "tight, isolated Au+Au DATA candidates: whole 1 GeV isolation bins with their centre below "
                         "the isolation line at the middle of the class"},
          {"scaling", "to the prompt curve over the DATA bins wholly above the tight cut, per the curves' bin width"}};
    });
}

/* E11/E33 of the photon candidates through the selection: before the preselection, after it and after the tight
   photon ID (the columns), for each Au+Au centrality class and the p+p reference at the bottom (the rows), on shared
   axes with one y range as in the IAN. DATA as points against two simulated lines, the candidates of the photon+jet
   simulation matched to a signal photon ("sig") and those of the inclusive-jet simulation not matched to one ("bkg"),
   each as fractions of its total and shown in bins of 0.04 (whole stored bins merged; a bin without DATA entries gets
   no point). The first column adds, in green, the shape of the DATA candidates whose cluster time minus MBD time lies
   in [-40, -5) ns, a sample enriched in non-collision background (NCB): as in the IAN it is scaled to one height in
   every row, its highest bin at 0.035, so that it reads beside the other distributions, and the number of tagged
   candidates is written under the row label. The timing window is typed here; it is not read from the saved
   selection. The column titles stand on one line in the top row and the row labels in the first column; the mark with
   the selection and the legend stand in the corners of the top row that the distributions leave free. No text stands
   on a guide line: a line of text that would do so starts under the guide line, and the legend entries stand between
   the guide lines. */
void e11Flow() {
  const std::string name = "shower_shapes/e11_e33_selection_flow.png";
  std::vector<std::tuple<std::string, std::string, std::string>> rows;
  for (const auto &c : classesOf("AuAu")) rows.push_back({"AuAu", c, classLabel(c)});
  for (const auto &c : classesOf("PP")) rows.push_back({"PP", c, "p+p reference"});
  bool anyData = false;
  for (const auto &[pkg, c, text] : rows)
    anyData = anyData || qa(pkg + "/data/photons/" + c + "/data/shape_e11_over_e33");
  if (!anyData) return missing(name, "no DATA QA");
  /* the IAN's bin width (about 0.04) */
  auto coarse = [](const Bins &b) { return rebinToWidth(b, 0.04); };
  /* one panel: the simulated lines with their colours, the DATA, and the timing-tagged DATA (timed: the panel has
     the timing of its candidates) with their number and their scaled shape */
  struct Panel {
    std::vector<std::pair<Bins, int>> simulated;
    Bins data, tagged;
    bool timed = false;
    double tags = 0;
  };
  std::vector<std::array<Panel, 3>> panels(rows.size());
  std::vector<std::string> selections;  /* of each row: the photon window and acceptance */
  Span s;
  double width = 0.04;
  const double tagHeight = 0.035;  /* of the highest bin of the timing-tagged shape, in every row */
  bool anyTagged = false;
  for (size_t r = 0; r < rows.size(); ++r) {
    const auto &[pkg, c, text] = rows[r];
    const auto domain = photonSupport(pkg);
    selections.push_back(gridWindowLine(domain.first, domain.second) + ", " +
                         photonEtaLine(selectionValue(pkg, "photon", "abs_eta_max", 0.7)));
    const std::string photons = "/photons/" + c;
    TH1 *data = qa(pkg + "/data" + photons + "/data/shape_e11_over_e33");
    for (int stage = 0; stage < 3; ++stage) {  /* the stage axis: every candidate, preselected, and tight */
      Panel &p = panels[r][stage];
      const std::vector<Cut> cuts = {{1, domain.first, domain.second}, {2, double(stage), stage + 1.0}};
      if (TH1 *h = qa(pkg + "/phosim" + photons + "/sig/shape_e11_over_e33"))
        p.simulated.push_back({unit(coarse(project(h, 0, cuts))), sampleLook("phosim").color});
      if (TH1 *h = qa(pkg + "/jetsim" + photons + "/bkg/shape_e11_over_e33"))
        p.simulated.push_back({unit(coarse(project(h, 0, cuts))), sampleLook("jetsim").color});
      if (data) {
        const Bins raw = coarse(project(data, 0, cuts));
        p.data = unit(raw);
        width = raw.width(0);
      }
      if (data && stage == 0)
        if (TH1 *t = qa(pkg + "/data" + photons + "/data/e11e33_dt")) {
          const Bins tagged = coarse(project(t, 0, {{1, domain.first, domain.second}, {2, -40, -5}}));
          const double highest = tagged.n() ? *std::max_element(tagged.v.begin(), tagged.v.end()) : 0;
          p.timed = true, p.tags = tagged.sum();
          if (highest > 0) p.tagged = scaled(tagged, tagHeight / highest), anyTagged = true;
        }
      const auto scope="e11/stage/"+std::to_string(stage)+"/row/"+std::to_string(r);
      p.data.seriesRole=scope+"/data";p.tagged.seriesRole=scope+"/timing_tag";
      for(auto& l:p.simulated)l.first.seriesRole=scope+"/simulation";
      for (const auto &l : p.simulated) s.add(l.first, false, 0, 1);
      s.add(p.data, true, 0, 1), s.add(p.tagged, false, 0, 1);
    }
  }
  /* the selection is written once, under the mark, when every row shares it, and under each row label otherwise */
  bool oneSelection = true;
  for (const auto &selection : selections) oneSelection = oneSelection && selection == selections.front();
  const int ny = rows.size();
  TCanvas cv("cv", "", 1600, 110 + 280 * ny);
  sharedGrid(cv, 3, ny, 20);
  const auto &g = gridLayout();
  /* the text sizes: half pixels, which set a symbol and its indices tight (see kLinePx) */
  const double titleSize = 24, rowSize = 21.5, lineSize = 17.5, markerSize = 0.9, between = kGridBlockGap * g.label;
  const int lineWidth = 3;
  /* the IAN's y range, 0 to 0.25, opened in steps of 0.05 if a bin stands higher */
  const double step = 0.05, ceiling = std::max(0.25, step * std::ceil(s.high / step - 1e-6));
  const char *titles[] = {"Before preselection", "After preselection", "After tight ID"};
  std::vector<std::string> labels = {"DATA", sampleLook("phosim").text, sampleLook("jetsim").text};
  if (anyTagged) labels.push_back("NCB-tagged DATA (timing, scaled)");
  const auto bold = [](const std::string &text) { return "#bf{" + text + "}"; };
  /* a light line at every labelled fraction; a frame's top is one of them, so they stand `band` pixels apart */
  std::vector<double> levels;
  for (double y = step; y < ceiling - 1e-6; y += step) levels.push_back(y);
  using Lines = std::vector<std::tuple<std::string, double, int>>;  /* text, font size and colour */
  /* lines hanging from the top of the current frame at its left or right side; one that would stand on a guide line
     starts under it */
  const auto stacked = [&](const Lines &lines, bool right) {
    const double band = g.height() * step / ceiling;
    double y = g.inset;
    for (const auto &[text, px, color] : lines) {
      const double guide = band * std::ceil((y + 0.1 * px) / band);  /* the first one under the top of the glyphs */
      if (guide < y + 1.25 * px) y = guide;
      gridLines({{text, px}}, right, false, y - g.inset, color);
      y += kGridPitch * px;
    }
  };
  for (int r = 0; r < ny; ++r)
    for (int stage = 0; stage < 3; ++stage) {
      const Panel &p = panels[r][stage];
      TH1F *axes = gridFrame(stage, r, 0, 0, 1, ceiling);
      axes->GetXaxis()->SetNdivisions(405);  /* labelled every 0.2 as in the IAN, a tick every 0.05 */
      axes->GetXaxis()->SetDecimals();
      guideLines(levels, 0, 1);
      gPad->RedrawAxis();
      for (const auto &l : p.simulated) photonjet_science::draw(line(l.first, l.second, 1, lineWidth),"HIST SAME");
      if (p.tagged.n()) photonjet_science::draw(line(p.tagged, cGreen, 1, lineWidth),"HIST SAME");
      if (p.data.n()) photonjet_science::draw(markers(populated(p.data), cBlack, 20, true, false, 0, markerSize),"P");
      /* at the top left: the column title in the top row; in the first column the row label and the tagged number */
      Lines lines;
      if (r == 0) lines.push_back({bold(titles[stage]), titleSize, cBlack});
      if (stage == 0) {
        lines.push_back({bold(std::get<2>(rows[r])), rowSize, cBlack});
        if (!oneSelection) lines.push_back({selections[r], lineSize, cBlack});
        if (p.timed) {
          lines.push_back({bold("NCB timing tag: " + grouped(p.tags)), lineSize, cGreen});
          photonjet_science::annotation(*gPad,"e11_flow/tagged_count/"+std::to_string(r),
            photonjet_science::scalar(p.tags,"weighted candidates","timing tag count"),p.tagged.physics);
        }
      }
      stacked(lines, false);
      const double reach = g.reach;
      if (r == 0 && stage == 0) {
        Lines mark = {{sPhenixInternalLabel(), titleSize, cBlack}};
        if (oneSelection) mark.push_back({selections[r], lineSize, cBlack});
        stacked(mark, true);
      }
      if (r == 0 && stage == 2) {
        /* the legend starts at the first guide line under the column title, a whole number of rows to each band */
        const double band = g.height() * step / ceiling;
        const double top = band * std::ceil((g.inset + reach + between) / band);
        const double pitch = band / std::max(1.0, std::round(band / (kGridPitch * lineSize)));
        const Bins key{{0, 1}, {0}, {0}};
        TLegend *l = gridLegend(labels, false, false, top - g.inset, lineSize, pitch);
        l->AddEntry(markers(key, cBlack, 20, false, false, 0, markerSize), labels[0].c_str(), "p");
        l->AddEntry(line(key, sampleLook("phosim").color, 1, lineWidth), labels[1].c_str(), "l");
        l->AddEntry(line(key, sampleLook("jetsim").color, 1, lineWidth), labels[2].c_str(), "l");
        if (anyTagged) l->AddEntry(line(key, cGreen, 1, lineWidth), labels[3].c_str(), "l");
        photonjet_science::draw(l);
      }
    }
  gridAxisTitles("E_{11}#scale[0.5]{ }/E_{33}", Form("Fraction of candidates per %g", width), true);
  savePhysicsQA(cv, name);
}

/* Truth isolation of the prompt photons of the photon simulation: the cumulative fraction below a threshold, direct
   and fragmentation photons, 15-20 and 25-30 GeV of truth pT; one panel per class on shared axes. The fractions and
   their binomial errors are computed here, at drawing time, from the saved truth_iso histogram: a fraction is taken of
   the photons inside the stored isolation range (no overflow), and a threshold counts the stored bins that end at or
   below it. Light guide lines mark the labelled fractions and a dashed line the nominal threshold, as in the IAN. The
   sample lines over the legend state the truth cone and the acceptance of the first panel's package. */
void truthIsolation() {
  const std::string name = "photonOnly/pp_auau_truth_isolation_direct_fragmentation_2x2.png";
  const auto order = panelOrder();
  std::vector<TH1 *> hs;
  bool any = false;
  for (const auto &[pkg, c] : order)
    hs.push_back(qa(pkg + "/phosim/truth/" + c + "/truth_iso")), any = any || hs.back();
  if (!any) return missing(name, "no photon-family simulation QA");
  const int nx = 2, ny = (order.size() + 1) / 2;
  TCanvas cv("cv", "", 1100, 440 * ny);
  const std::string simulation = simulationNote(order);
  sharedGrid(cv, nx, ny, 22);
  const auto &g = gridLayout();
  /* the IAN's range: the fractions lie between 0.8 and 1 (full: 1 and the marker on it); above them the frames
     open for the line of the panel titles */
  const double floor = 0.77, full = 1.006, ceiling = gridCeiling(0, floor, full, gridBlock({g.title}));
  /* the thresholds, 1 to 10 GeV in steps of 1 GeV, each at the middle of a unit of the x axis */
  const int thresholds = 10;
  const double x0 = 0.5, x1 = thresholds + 0.5;
  const std::string labelled = order.front().first;  /* the package whose truth cone and acceptance the legend states */
  /* the prompt class on the third axis of truth_iso with its colour and name, and the truth-pT windows with their
     markers */
  const std::vector<std::tuple<int, int, std::string>> kinds = {{1, cBlue, "Direct"}, {2, cOrange, "Fragmentation"}};
  std::vector<std::tuple<double, double, int>> windows;   /* qa.windows.truth_isolation_gev, markers 20, 21, ... */
  for (const auto &[lo, hi] : qaWindowList(labelled, "truth_isolation_gev", {{15.0, 20.0}, {25.0, 30.0}}))
    windows.emplace_back(lo, hi, 20 + static_cast<int>(windows.size()));
  PlotJson usedWindows=PlotJson::array();for(const auto& [lo,hi,marker]:windows)usedWindows.push_back({lo,hi});
  photonjet_science::renderSetting(cv,*gQA,"QA.Windows.TruthIsolation",usedWindows);
  /* the sample and its input note, written over the legend entries as lines without a symbol */
  const double radius = selectionValue(labelled, "truth_photon", "isolation_radius", 0.3);
  const double etaMax = selectionValue(labelled, "qa", "stored_abs_eta_max", 0.7);
  std::vector<std::string> sample = {"Photon+jet simulation", Form("Truth isolation R = %g", radius),
                                     "Truth " + photonEtaLine(etaMax)};
  if (!simulation.empty()) sample.push_back(simulation);
  std::vector<std::string> labels;
  for (const auto &kind : kinds)
    for (const auto &window : windows)
      labels.push_back(Form("%s, %g-%g GeV", std::get<2>(kind).c_str(), std::get<0>(window), std::get<1>(window)));
  for (size_t p = 0; p < order.size(); ++p) {
    const std::string &pkg = order[p].first, &c = order[p].second;
    TH1F *axes = gridFrame(p % nx, p / nx, x0, floor, x1, ceiling);
    axes->GetXaxis()->SetNdivisions(110);  /* a tick at every threshold, labelled as in the IAN: 1, 2, 4, 6, 8, 10 */
    for (int odd = 3; odd < 10; odd += 2) axes->GetXaxis()->ChangeLabelByValue(odd, -1, 0);
    axes->GetYaxis()->SetNdivisions(106);
    /* a light line at every labelled fraction and the dashed line at the nominal threshold, the axes over them */
    guideLines({0.80, 0.85, 0.90, 0.95, 1.00}, x0, x1);
    /* the nominal threshold of the panel's package as a fraction of the frame width */
    const double nominalAt = (selectionValue(pkg, "truth_photon", "isolation_max_gev", 4) - x0) / (x1 - x0);
    gridDashedLine(nominalAt, 0, nominalAt, (full - floor) / (ceiling - floor));
    gPad->RedrawAxis();
    gridLines({{classLabel(c) + ", " + energyOf(pkg), g.title}});
    /* the first panel: the mark on the title's line, and the sample and the legend at the lower right, which the
       rising fractions leave free. The rows stand between the guide lines: a whole number of them to the band
       between two guide lines, counted from the lowest one. */
    TLegend *l = nullptr;
    if (p == 0) {
      gridLines({{sPhenixInternalLabel(), g.title}}, true);
      std::vector<std::string> block = sample;
      if (hs[p]) block.insert(block.end(), labels.begin(), labels.end());
      const double perUnit = g.height() / (ceiling - floor), band = 0.05 * perUnit, lowest = (0.80 - floor) * perUnit;
      const double pitch = band / std::max(1.0, std::round(band / (kGridPitch * g.text)));
      const double foot = lowest - pitch * std::floor((lowest - g.inset) / pitch);
      l = gridLegend(block, true, true, foot - g.inset, 0, pitch);
      for (const auto &text : sample) {
        TLegendEntry *written = l->AddEntry(static_cast<TObject *>(nullptr), text.c_str(), "");
        if (text == simulation) written->SetTextColor(kGray + 2);
      }
    }
    if (l && !hs[p]) photonjet_science::draw(l);
    if (!hs[p]) continue;
    size_t entry = 0;
    for (const auto &[kind, color, who] : kinds)
      for (const auto &[lo, hi, marker] : windows) {
        const std::string &label = labels[entry++];
        const Bins v = project(hs[p], 0, {{1, lo, hi}, {2, double(kind), kind + 1.0}});
        const double total = v.sum();
        if (!(total > 0)) continue;
        /* weighted photons: the binomial error takes the effective number of entries, (sum w)^2 / sum w^2 */
        double squares = 0;
        for (double e2 : v.e2) squares += e2;
        const double entries = squares > 0 ? total * total / squares : total;
        auto *points = new TGraphErrors;
        for (int cut = 1; cut <= thresholds; ++cut) {
          double below = 0;
          for (int i = 0; i < v.n(); ++i) if (v.edges[i + 1] <= cut + 1e-9) below += v.v[i];
          const double f = below / total;
          points->SetPoint(cut - 1, cut, f);
          points->SetPointError(cut - 1, 0, std::sqrt(std::max(f * (1 - f), 0.0) / std::max(entries, 1.0)));
        }
        photonjet_science::bindNative(points,"truth_isolation/"+std::to_string(entry-1),
          v.physics.then("cumulative_fraction",{{"thresholds",thresholds},{"total",total},{"effective_entries",entries},{"error","binomial effective count"}}),
          photonjet_science::covariance("not_computed",{},0,"cross-threshold covariance not calculated"));
        style(points, color, marker, 1.2);
        photonjet_science::draw(points,"P");
        if (l) l->AddEntry(points, label.c_str(), "pe");
      }
    if (l) photonjet_science::draw(l);
  }
  gridAxisTitles("Truth isolation threshold E_{T,cut}^{iso} [GeV]",
                 "Cumulative fraction of truth photons passing E_{T}^{iso} < E_{T,cut}^{iso}", true);
  savePhysicsQA(cv, name);
}

/* The prompt photons below the truth-isolation maximum of the selection against truth pT, over the photon range of
   the selection: all, direct and fragmentation photons above, the share of direct photons below, in a lower frame 40%
   as tall as the upper one; one figure per class, all laid out alike as in the IAN. The share and its statistical
   error are computed here, at drawing time, from the two saved spectra (shareOf). The header stands at the top left
   over a band the y range is opened for, the legend at the top right and the selection at the bottom left, which a
   falling spectrum leaves free. The lowest point sits at 8% of the frame height in every figure, or a little higher
   where a decade would otherwise fall on the bottom edge. Light lines mark the decades and the labelled shares: they
   lie under the axis ticks and give way to the selection text. The lower panel is labelled every 0.04 around the drawn
   shares, with a long tick at each label and a short one every 0.01. */
void truthPtAfterIsolation() {
  for (const auto &[pkg, c] : figureOrder()) {
    const std::string name = "photonOnly/" + (pkg == "PP" ? std::string("pp") : "auau_" + classTag(c)) +
                             "_truth_photon_pt_after_isolation.png";
    TH1 *h = qa(pkg + "/phosim/truth/" + c + "/truth_iso");
    if (!h) { missing(name, "no photon-family simulation QA"); continue; }
    const double isoMax = selectionValue(pkg, "truth_photon", "isolation_max_gev", 4);
    const auto domain = photonSupport(pkg);
    /* the third axis of truth_iso is the prompt class: [1, 2) direct, [2, 3) fragmentation */
    const Bins direct = namedBins(within(project(h, 1, {{0, 0, isoMax}, {2, 1, 2}}), domain.first, domain.second),"truth/direct"),
               fragmentation = namedBins(within(project(h, 1, {{0, 0, isoMax}, {2, 2, 3}}), domain.first, domain.second),"truth/fragmentation");
    const Bins all = namedBins(sum(direct, fragmentation),"truth/prompt_sum");
    Span s;
    for (const Bins *b : {&all, &direct, &fragmentation}) s.add(perWidth(*b), true, domain.first, domain.second);
    if (!(s.positive < s.high)) { missing(name, "no isolated prompt photons in the photon range"); continue; }
    TCanvas cv("cv", "", 760, 700);
    auto [top, bottom] = ratioPads(cv);
    /* the lower pad takes 37% of the canvas; its margins keep their pixels, and so the x labels and title their room */
    const double split = 0.37, tick = 6, markerSize = 0.85;
    bottom->SetBottomMargin(bottom->GetBottomMargin() * bottom->GetHNDC() / split);
    bottom->SetTopMargin(bottom->GetTopMargin() * bottom->GetHNDC() / split);
    top->SetPad(0, split, 1, 1);
    bottom->SetPad(0, 0, 1, split);
    top->cd();
    gPad->SetLogy();
    const double w = padWidthPx(), hpx = padHeightPx();
    const double base = gPad->GetBottomMargin(), ceiling = 1 - gPad->GetTopMargin();
    const double frameLeft = gPad->GetLeftMargin(), frameWidth = 1 - frameLeft - gPad->GetRightMargin();
    /* the header starts and the legend text ends at one distance from the frame */
    const double left = frameLeft + kInsetPx / w, right = frameLeft + frameWidth - (kInsetPx - 10) / w;
    const double mark = ceiling - (14 + kMarkPx) / hpx, legendTop = ceiling - 10 / hpx;
    const std::vector<std::string> system = {classLabel(c) + ", " + energyOf(pkg)};
    const std::string simulation =
        std::string(pkg == "PP" ? "PYTHIA 8" : "Embedded PYTHIA 8") + " photon+jet simulation";
    const std::string note = simulationNote({pkg});
    /* the simulation line 27 px under the system line and the input note 25 px under that */
    const double simulationY = header(left, mark, system, false, false) - 27 / hpx, noteY = simulationY - 25 / hpx;
    /* the y range: the lowest point at 8% of the frame height, the highest 22 px under the lowest of these lines
       (clear: that height as a fraction of the frame); a decade within 6% of the bottom edge would put its label on
       the corner and its line on the axis, so the frame is then set by that decade at 6% instead. logHeight: the
       natural logarithm of the ratio of the frame's upper to its lower edge. */
    const double clear = ((note.empty() ? simulationY : noteY) - 22 / hpx - base) / (ceiling - base);
    double logHeight = std::log(s.high / s.positive) / (clear - 0.08), low = s.positive * std::exp(-0.08 * logHeight);
    if (const double decade = std::pow(10, std::ceil(std::log10(low))); std::log(decade / low) < 0.06 * logHeight)
      logHeight = std::log(s.high / decade) / (clear - 0.06), low = decade * std::exp(-0.06 * logHeight);
    TH1F *f = photonjet_science::drawFrame(*gPad,domain.first, low, domain.second, low * std::exp(logHeight));
    upperAxes(f);
    f->GetXaxis()->SetNdivisions(505);
    tickLengths(f, tick);
    /* the selection at the bottom left, and the part of the frame its two lines take; the bars of the second line are
       plain ones of one height, and the vertex symbol is set apart, so that its index sits on it (apart). The lines
       are written at a half-pixel size, at which an index sits tight to the text before it (see kLinePx), on
       baselines 52 and 20 px above the bottom of the frame. */
    const std::string isolation = Form("E_{T,truth}^{iso} < %g GeV, R = %g", isoMax,
                                       selectionValue(pkg, "truth_photon", "isolation_radius", 0.3));
    const std::string acceptance =
        Form("%s < %g, |%s| < %g cm", absolute("#eta^{#gamma}_{truth}").c_str(), selectionValue(pkg, "qa", "stored_abs_eta_max", 0.7),
             apart("z_{vtx}").c_str(), selectionValue(pkg, "event", "vertex_abs_z_max_cm", 60));
    const double selectionPx = 19.5, isolationY = base + 52 / hpx, acceptanceY = base + 20 / hpx;
    const double selectionWidth = std::max(textWidth(isolation, selectionPx), textWidth(acceptance, selectionPx));
    const double selectionEnd =
        domain.first + (domain.second - domain.first) * (left + selectionWidth + 4 / w - frameLeft) / frameWidth;
    TLine grid;  /* the IAN's light lines at the decades the points reach */
    grid.SetLineColor(guideColor());
    for (double y = std::pow(10, std::ceil(std::log10(low))); y < s.high; y *= 10) {
      const double at = base + (ceiling - base) * std::log(y / low) / logHeight;
      /* a line from 14 px under the lower baseline to 30 px over the upper one would cross the selection text */
      const bool behindText = at > acceptanceY - 14 / hpx && at < isolationY + 30 / hpx;
      grid.DrawLine(behindText ? selectionEnd : domain.first, y, domain.second, y);
    }
    gPad->RedrawAxis();
    const char *allLabel = "All prompt photons";  /* the widest of the three entries, which sizes the legend */
    const double symbol = 46 / w, entry = symbol + textWidth(allLabel, kLinePx) + 10 / w;
    TLegend *leg = ianLegend(right - entry, legendTop - 3 * kStepPx / hpx, right, legendTop, symbol / entry);
    for (const auto &[b, color, marker, text] :
         {std::make_tuple(&all, cBlack, 20, allLabel), std::make_tuple(&direct, cRed, 21, "Direct"),
          std::make_tuple(&fragmentation, cBlue, 22, "Fragmentation")}) {
      TGraphErrors *g = markers(perWidth(*b), color, marker, true, false, 0, markerSize);
      photonjet_science::draw(g,"P");
      leg->AddEntry(g, text, "ep");
    }
    photonjet_science::draw(leg);
    header(left, mark, system, true, false);
    TLatex t;
    t.SetNDC();
    t.SetTextFont(43);
    t.SetTextSize(17.5);  /* the simulation line ends well short of the legend's symbols */
    t.SetTextColor(gGray);
    t.DrawLatex(left, simulationY, simulation.c_str());
    if (!note.empty()) t.DrawLatex(left, noteY, note.c_str());
    t.SetTextSize(selectionPx);
    t.SetTextColor(cBlack);
    t.DrawLatex(left, isolationY, isolation.c_str());
    t.DrawLatex(left, acceptanceY, acceptance.c_str());
    bottom->cd();
    /* the share moves by a few percent over the range: labels every 0.04 with at least half a step around the shares,
       and a frame that ends half a step beyond the outer labels */
    const Bins fraction = shareOf(direct, fragmentation);
    Span share;
    share.add(fraction, true, domain.first, domain.second);
    const double step = 0.04, first = step * std::floor((share.low - step / 2) / step),
                 last = step * std::ceil((share.high + step / 2) / step);
    const int halves = 2 * std::lround((last - first) / step) + 2;
    TH1F *r =
        ratioFrame(domain.first, domain.second, first - step / 2, last + step / 2, "Truth photon p_{T} [GeV]", "");
    r->GetXaxis()->SetNdivisions(505);
    r->GetXaxis()->CenterTitle();
    r->GetYaxis()->SetNdivisions(halves + 200, false);
    r->GetYaxis()->SetDecimals();
    for (int i = 1; i <= halves + 1; i += 2) r->GetYaxis()->ChangeLabel(i, -1, 0);
    tickLengths(r, tick);
    r->GetYaxis()->SetTickLength(0);  /* the y ticks are drawn here: long at the labels, short every 0.01 between */
    for (double y = first; y < last + step / 4; y += step) grid.DrawLine(domain.first, y, domain.second, y);
    gPad->RedrawAxis();
    TLine shareTick;
    const double perPx = (domain.second - domain.first) / (padWidthPx() * frameWidth);
    for (int k = std::lround((first - step / 2) / 0.01) + 1; 0.01 * k < last + step / 2 - 1e-9; ++k) {
      const double y = 0.01 * k, length = (std::fabs(std::remainder(y - first, step)) < 1e-9 ? tick : tick / 2) * perPx;
      shareTick.DrawLine(domain.first, y, domain.first + length, y);
      shareTick.DrawLine(domain.second - length, y, domain.second, y);
    }
    photonjet_science::draw(markers(fraction, cRed, 21, true, false, 0, markerSize),"P");
    /* the y titles, each centred on its frame as in the IAN: on the canvas, where the lower one is not cut by its
       pad */
    cv.cd();
    TLatex title;
    title.SetNDC();
    title.SetTextFont(43);
    title.SetTextSize(26);  /* the size of the axis titles (pixelAxes) */
    title.SetTextAngle(90);
    title.SetTextAlign(22);
    /* centred 32 px from the left edge of the canvas, at the middle of each frame */
    for (const auto &[pad, text] :
         {std::make_pair(top, "Weighted photons / GeV"), std::make_pair(bottom, "Direct / total")})
      title.DrawLatex(32 / w,
                      pad->GetYlowNDC() + pad->GetHNDC() * (1 + pad->GetBottomMargin() - pad->GetTopMargin()) / 2,
                      text);
    savePhysicsQA(cv, name);
  }
}

/* A Gaussian through the core of a drawn distribution, in its colour: the bins around the highest one that reach 40%
   of its height (about 1.35 standard deviations to either side of a Gaussian peak), fitted to the bin integrals with
   their statistical errors. The fit is made here, at drawing time, on the displayed bins; its parameters are not
   saved by the qa step. Null when the core has fewer than four bins or the fit does not converge. */
TF1 *coreGaussian(const Bins &b, int color) {
  const double coreFraction = 0.4;  /* of the height of the highest bin: the bins that reach it are the core */
  const double startWidth = 0.4;    /* the width the fit starts from, as a fraction of the fitted range */
  int peak = -1;
  for (int i = 0; i < b.n(); ++i) if (std::isfinite(b.v[i]) && (peak < 0 || b.v[i] > b.v[peak])) peak = i;
  if (peak < 0) return nullptr;
  int first = peak, last = peak;
  while (first > 0 && b.v[first - 1] >= coreFraction * b.v[peak]) --first;
  while (last + 1 < b.n() && b.v[last + 1] >= coreFraction * b.v[peak]) ++last;
  if (last - first < 3) return nullptr;
  TH1D core(fresh().c_str(), "", last - first + 1, b.edges.data() + first);
  core.SetDirectory(nullptr);
  for (int i = first; i <= last; ++i)
    core.SetBinContent(i - first + 1, b.v[i]), core.SetBinError(i - first + 1, std::sqrt(b.e2[i]));
  auto *fit = new TF1(fresh().c_str(), "gaus", b.edges[first], b.edges[last + 1]);
  fit->SetParameters(b.v[peak], b.center(peak), startWidth * (b.edges[last + 1] - b.edges[first]));
  const TFitResultPtr fitResult = core.Fit(fit, "QNRIS");
  const int status = fitResult;
  if (status != 0) { delete fit; return nullptr; }
  fit->SetLineColor(color);
  fit->SetLineWidth(2);
  photonjet_science::bindFit(fit,"isolation_core_gaussian/"+b.seriesRole,b.physics,fitResult,
    {{"options","QNRI"},{"first_bin",first},{"last_bin",last},{"start_width_fraction",startWidth},
     {"initial",{b.v[peak],b.center(peak),startWidth*(b.edges[last+1]-b.edges[first])}}});
  return fit;
}

/* A fitted function as a smooth line over its own range, in its colour and width. */
TGraph *smoothCurve(const TF1 *f) {
  const int steps = 200;
  auto *g = new TGraph(steps + 1);
  for (int i = 0; i <= steps; ++i) {
    const double x = f->GetXmin() + (f->GetXmax() - f->GetXmin()) * i / steps;
    g->SetPoint(i, x, f->Eval(x));
  }
  g->SetLineColor(f->GetLineColor());
  g->SetLineWidth(f->GetLineWidth());
  const auto fit=photonjet_science::tag(f,"physics_series");
  if(!fit.is_null())photonjet_science::bind(g,fit.at("role"),photonjet_science::Trace::fromJson(fit.at("trace")),fit.at("value"),fit.at("covariance"));
  return g;
}

/* The legend of the two isolation figures (isolationOverlays, isolationStacks), smaller than the header's: the width
   of a symbol, the height of a row and the text size, in pixels. */
const double kIsolationKeySymbolPx = 40, kIsolationKeyRowPx = 28, kIsolationKeyTextPx = 20;

/* Reconstructed isolation of Au+Au DATA candidates at 14-16 GeV in five centrality ranges, each as fractions of its
   total in 1 GeV bins, without and with the underlying-event subtraction, on the IAN's wide canvas with the text at
   the right, each on its own y range. The five ranges are typed here on the 5% centrality axis of the saved slices;
   they are not the classes of the file. A Gaussian is fitted to the core of each range at drawing time and drawn
   through it (coreGaussian). With the subtraction a dashed line marks the Gaussian mean of the most peripheral range,
   whose fitted mean and width are quoted on the figure; they are not saved by the qa step, and when that fit does not
   converge neither the line nor the numbers are drawn. The first text line is the trigger by its name. */
void isolationOverlays() {
  /* GeV of E_T (qa.windows.isolation_overlay_gev; the file name keeps the IAN's _pT_14_16), GeV of isolation, px */
  const auto overlayWindow = qaWindow("AuAu", "isolation_overlay_gev", {14, 16});
  const double low = overlayWindow.first, high = overlayWindow.second, x0 = -10, x1 = 40, tick = 14;
  const std::string folder = "AuAu/data/slices/data/";
  struct Range { std::string text; int color; Bins bins; };
  /* the distributions of a saved slice histogram in the IAN's colours, central to peripheral, in 1 GeV bins: ten of
     the stored 0.1 GeV bins to one. lo, hi: the first slice of a range and the one after its last. */
  const auto rangesOf = [&](const TH1 *h) {
    std::vector<Range> out;
    for (const auto &[lo, hi, color] : {std::make_tuple(0, 2, kRed + 1), std::make_tuple(2, 4, kBlue + 1),
                                         std::make_tuple(4, 8, kGreen + 2), std::make_tuple(8, 12, kMagenta + 1),
                                         std::make_tuple(12, 16, kOrange + 1)}) {
      const Bins b = populated(unit(rebin(project(h, 0, {{1, low, high}, {2, double(lo), double(hi)}}), 10)));
      out.push_back({percent(kSlicePercent * lo, kSlicePercent * hi), color,namedBins(b,"isolation/centrality_slot/"+std::to_string(out.size()))});
    }
    return out;
  };
  for (const auto &[kind, stem, note] :
       {std::make_tuple("noUE", "isonoueslice", "without UE subtraction"),
        std::make_tuple("withUE", "isoslice", "with SUB1 UE subtraction")}) {
    const std::string name = std::string("isolation/isoCentOverlay_") + kind + "_pT_14_16.png";
    TH1 *h = qa(folder + stem);
    if (!h) { missing(name, "no Au+Au DATA QA"); continue; }
    const bool subtracted = std::string(stem) == "isoslice";
    const std::vector<Range> ranges = rangesOf(h);
    Span s;
    for (const auto &r : ranges) s.add(r.bins, true, x0, x1);
    TCanvas cv("cv", "", 880, 615);
    photonjet_science::renderSetting(cv,*gQA,"QA.Windows.IsolationOverlay",PlotJson::array({low,high}));
    cv.SetLeftMargin(0.125);
    cv.SetRightMargin(0.03);
    const double top = s.high / 0.94;
    TH1F *f = photonjet_science::drawFrame(*gPad,x0, 0, x1, top);
    f->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
    f->GetYaxis()->SetTitle("Normalized to unit area");
    f->GetXaxis()->SetNdivisions(505);
    f->GetYaxis()->SetTitleOffset(1.1);
    tickLengths(f, tick);
    std::vector<TGraphErrors *> points;
    for (const auto &r : ranges) {
      points.push_back(markers(r.bins, r.color, 20, true, true, 0, 0.9));
      photonjet_science::draw(points.back(),"P");
    }
    TF1 *peripheral = nullptr;  /* after the loop: the fit of the last, most peripheral range (null if it failed) */
    for (const auto &r : ranges) {
      peripheral = coreGaussian(r.bins, r.color);
      if (peripheral) photonjet_science::draw(smoothCurve(peripheral),"L");
    }
    if (subtracted && peripheral) {
      TLine at;
      at.SetLineColor(ranges.back().color);
      at.SetLineStyle(kDashed);
      at.SetLineWidth(2);
      at.DrawLine(peripheral->GetParameter(1), 0, peripheral->GetParameter(1), top);
    }
    /* the text at the right, as in the IAN: the selection, the legend, then the sPHENIX mark and the system, set in
       from the frame by the tick length and a margin; at half-pixel sizes (18.5, 19.5, 21.5) a symbol or an index sits
       tight to the text before it */
    const double w = padWidthPx(), hpx = padHeightPx(), right = 1 - gPad->GetRightMargin() - (tick + 24) / w;
    double y = 1 - gPad->GetTopMargin() - (tick + 4) / hpx;
    TLatex t;
    t.SetNDC();
    t.SetTextFont(43);
    t.SetTextAlign(31);
    const auto write = [&](const std::string &text, double px, int color) {
      y -= 1.35 * px / hpx;
      t.SetTextSize(px);
      t.SetTextColor(color);
      t.DrawLatex(right, y, text.c_str());
    };
    if (const auto &selection = renderSelection("AuAu", true); !selection.is_null()) {
      /* the trigger's own vertex condition stays in the trigger's words: it is not the analysis vertex window */
      std::string trigger = selection.at("event").at("trigger_name").get<std::string>();
      if (const auto at = trigger.find(">="); at != std::string::npos) trigger.replace(at, 2, "#geq");
      write(trigger, 19.5, cBlack);
    }
    write(photonWindowLabel(low, high) + Form(", #DeltaR_{cone} < %g", isolationRadius("AuAu")), 19.5, cBlack);
    write(note, 19.5, cBlack);
    if(subtracted && peripheral)for(int parameter:{1,2})photonjet_science::annotation(*gPad,"isolation_fit_parameter/"+std::to_string(parameter),
      photonjet_science::scalar(peripheral->GetParameter(parameter),"GeV","Gaussian fit parameter",{{"error",photonjet_science::real(peripheral->GetParError(parameter))}}),
      photonjet_science::Trace::fromJson(photonjet_science::tag(peripheral,"physics_series").at("trace")));
    if (subtracted && peripheral)
      for (const auto &[greek, parameter] : {std::make_pair("#mu", 1), std::make_pair("#sigma", 2)})
        write(Form("%s^{Gauss} (%s) = %.2f #pm %.2f GeV", greek, ranges.back().text.c_str(),
                   peripheral->GetParameter(parameter), peripheral->GetParError(parameter)),
              18.5, ranges.back().color);
    /* five entries in two columns make three rows; a column is as wide as a symbol and a label of the width of
       "10-20%" */
    const double symbol = kIsolationKeySymbolPx / w, rows = 3 * kIsolationKeyRowPx / hpx;
    const double column = symbol + textWidth("10-20%", kIsolationKeyTextPx) + 16 / w;
    y -= 12 / hpx;
    TLegend *leg = ianLegend(right - 2 * column, y - rows, right, y, symbol / column);
    leg->SetNColumns(2);
    leg->SetTextSize(kIsolationKeyTextPx);
    for (size_t i = 0; i < ranges.size(); ++i) leg->AddEntry(points[i], ranges[i].text.c_str(), "ep");
    photonjet_science::draw(leg);
    y -= rows + 10 / hpx;
    write(sPhenixInternalLabel(), kMarkPx, cBlack);
    write("Au+Au, " + energyOf("AuAu"), kLinePx, cBlack);
    if (const std::string input = inputNote(); !input.empty()) write(input, kLinePx, gGray);
    savePhysicsQA(cv, name);
  }
}

/* The isolation requirement of the Au+Au selection as the IAN writes it: the cone energy, with its radius and its
   underlying-event subtraction, below the saved straight line in the centrality percentile c. Empty for p+p, whose
   selection stores the decision per photon and no line. */
std::string isolationRequirement(const std::string &pkg) {
  const auto &m = renderSelection(pkg, true);
  if (pkg == "PP" || m.is_null()) return "";
  const auto &isolation = m.at("photon_id").at("isolation");
  const auto &threshold = isolation.at("threshold");
  if (!threshold.is_object()) return "";
  const double intercept = threshold.at("intercept").get<double>(), slope = threshold.at("slope").get<double>();
  const bool sub1 = isolation.at("subtraction_method").get<int>() == 2;  /* the producer's number for SUB1 */
  return Form("E_{T,R=%g}^{iso%s} < (%.2f %s %.3g c) GeV", isolation.at("radius").get<double>(), sub1 ? ",SUB1" : "",
              intercept, slope < 0 ? "#minus" : "+", std::fabs(slope));
}

/* Reconstructed isolation of the tight DATA candidates over the photon range of the selection, against the non-tight
   DATA scaled to it above 6 GeV and the prompt photons of the photon simulation scaled to the remainder, stacked; one
   figure per class, all laid out alike on the IAN's portrait canvas: the header at the top left, the legend at the
   top right and the highest point at 70% of the frame height. The two scale factors are computed here, at drawing
   time, from the saved isolation histograms (the sums run over the stored isolation range, without overflow): they
   are a normalisation for display, not the ABCD estimate of the analysis. The last header line is the isolation
   requirement where the selection carries its line (Au+Au) and the cone otherwise. p+p in 0.5 GeV bins from -2 to
   15 GeV; Au+Au in 1 GeV bins up to 20 GeV, from -20 GeV or from where the saved isolation axis starts. */
void isolationStacks() {
  const double sideband = 6.0;  /* GeV of isolation: above it the non-tight DATA are scaled to the tight DATA */
  const Cut tightId = {2, 0, 1}, nonTightId = {2, 1, 2};  /* of the photon-ID axis of iso */
  for (const auto &[pkg, c] : figureOrder()) {
    const std::string name =
        "isolation/" + (pkg == "PP" ? std::string("pp_nominal_topo") : "auau_" + classTag(c)) + "_isoStack.png";
    TH1 *data = qa(pkg + "/data/photons/" + c + "/data/iso"), *sig = qa(pkg + "/phosim/photons/" + c + "/sig/iso");
    if (!data || !sig) { missing(name, "needs DATA and photon-family simulation QA"); continue; }
    const auto domain = photonSupport(pkg);
    const Bins tight = project(data, 0, {{1, domain.first, domain.second}, tightId}),
               nonTight = project(data, 0, {{1, domain.first, domain.second}, nonTightId});
    Bins signal = project(sig, 0, {{1, domain.first, domain.second}, tightId});
    double tightSide = 0, nonTightSide = 0;
    for (int i = 0; i < tight.n(); ++i)
      if (tight.edges[i] >= sideband) tightSide += tight.v[i], nonTightSide += nonTight.v[i];
    const Bins background = scaled(nonTight, nonTightSide > 0 ? tightSide / nonTightSide : 0.0);
    signal = scaled(signal, signal.sum() > 0 ? (tight.sum() - background.sum()) / signal.sum() : 0.0);
    const bool pp = pkg == "PP";
    const double width = pp ? 0.5 : 1.0;  /* GeV per drawn bin, as in the IAN */
    const int k = std::lround(width / tight.width(0));
    const Bins t = namedBins(rebin(tight, k),"isolation/tight"), b = namedBins(rebin(background, k),"isolation/non_tight_scaled"), stacked = namedBins(sum(rebin(signal, k), b),"isolation/signal_plus_background");
    const double x0 = pp ? -2 : std::max(-20.0, t.edges.front()), x1 = pp ? 15 : 20;
    Span s;
    s.add(t, true, x0, x1), s.add(stacked, false, x0, x1);
    TCanvas cv("cv", "", 760, 1064);
    cv.SetMargin(0.16, 0.03, 0.095, 0.025);
    const double w = padWidthPx(), hpx = padHeightPx();
    const double base = gPad->GetBottomMargin(), ceiling = 1 - gPad->GetTopMargin();
    const double left = gPad->GetLeftMargin() + kInsetPx / w, right = 1 - gPad->GetRightMargin() - 14 / w;
    const double mark = ceiling - (14 + kMarkPx) / hpx, legendTop = ceiling - 16 / hpx;
    const std::string requirement = isolationRequirement(pkg);
    const std::string cone = Form("Isolation cone R = %g%s", isolationRadius(pkg), pp ? "" : ", UE subtracted");
    const std::vector<std::string> lines = {classLabel(c) + ", " + energyOf(pkg),
                                            photonWindowLabel(domain.first, domain.second), photonEtaLabel(pkg),
                                            requirement.empty() ? cone : requirement};
    const double peak = std::min(0.70, (header(left, mark, lines, false) - 16 / hpx - base) / (ceiling - base));
    TH1F *f = photonjet_science::drawFrame(*gPad,x0, 0, x1, s.high / peak);
    f->GetXaxis()->SetTitle("E_{T}^{iso,reco} [GeV]");
    f->GetYaxis()->SetTitle(Form("Counts / %g GeV", width));
    pixelAxes(f);
    f->GetXaxis()->CenterTitle();
    f->GetYaxis()->CenterTitle();
    tickLengths(f, 8);  /* short of the header, the legend and the filled stack, as in the IAN */
    const int edge = TColor::GetColor("#3454B4");
    /* both components outlined in blue, as in the IAN; the legend shows the control sample with its own red edge */
    TH1D *all = line(stacked, edge), *under = line(b, edge),
         *control = line(Bins{{0, 1}, {0}, {0}}, TColor::GetColor("#C8605A"));
    all->SetFillColor(TColor::GetColor("#C4C8F2"));
    under->SetFillColor(TColor::GetColor("#F6C6C3"));
    control->SetFillColor(under->GetFillColor());
    photonjet_science::draw(all,"HIST SAME");
    photonjet_science::draw(under,"HIST SAME");
    TGraphErrors *g = markers(t, cBlack, 20, true, false, 0, 1.0);
    photonjet_science::draw(g,"P");
    gPad->RedrawAxis();
    const char *controlLabel = "Non-tight data control";  /* the widest of the three entries, which sizes the legend */
    const double symbol = kIsolationKeySymbolPx / w, rows = 3 * kIsolationKeyRowPx / hpx;
    const double entry = symbol + textWidth(controlLabel, kIsolationKeyTextPx) + 10 / w;
    TLegend *leg = ianLegend(right - entry, legendTop - rows, right, legendTop, symbol / entry);
    leg->SetTextSize(kIsolationKeyTextPx);
    leg->AddEntry(g, "Tight-photon data", "p");
    leg->AddEntry(control, controlLabel, "f");
    leg->AddEntry(all, "Signal MC", "f");
    photonjet_science::draw(leg);
    TLatex sample;
    sample.SetNDC();
    sample.SetTextFont(43);
    sample.SetTextSize(17);
    sample.SetTextAlign(31);
    sample.DrawLatex(1 - gPad->GetRightMargin() - kInsetPx / w, legendTop - rows - 24 / hpx,
                     pp ? "Signal MC: PYTHIA 8 photon+jet" : "Signal MC: embedded photon+jet");
    header(left, mark, lines);
    savePhysicsQA(cv, name);
  }
}

/* The photon E_T range of the ABCD, purity and efficiency panels in GeV. It is typed here, one x axis for p+p and
   Au+Au, and is not read from the saved selection. */
const double kPanelEtLow = 15, kPanelEtHigh = 35;

/* DATA event-leading candidates per region (A-D) in each ABCD photon bin, as counts per bin with their statistical
   errors; one panel per class on shared axes, each with the mark and its title at the top right as in the IAN. The
   first axis of the saved abcd histogram is the region, A to D as 0 to 3. */
void abcdCounts() {
  const std::string name = "photonOnly/pp_auau_abcd_region_counts_2x2.png";
  const auto order = panelOrder();
  std::vector<TH1 *> hs;
  std::vector<std::array<Bins, 4>> yields(order.size());
  Span s;
  for (size_t p = 0; p < order.size(); ++p) {
    hs.push_back(qa(order[p].first + "/abcd/" + order[p].second));
    if (!hs.back()) continue;
    for (int r = 0; r < 4; ++r) yields[p][r] = project(hs.back(), 1, {{0, double(r), r + 1.0}});
    for (const Bins &b : yields[p]) s.add(b, true, kPanelEtLow, kPanelEtHigh);
  }
  if (!std::isfinite(s.high)) return missing(name, "no DATA analysis fills");
  const auto [nx, ny] = twoColumns(order.size());
  TCanvas cv("cv", "", 550 * nx, 450 * ny);
  sharedGrid(cv, nx, ny, 22);
  const auto &g = gridLayout();
  /* the regions: colour, marker, marker size (equal areas) and legend entry */
  const std::array<std::tuple<int, int, double, const char *>, 4> regions = {
      {{cBlack, 20, 1.6, "A: tight, isolated"}, {cRed, 21, 1.5, "B: tight, non-isolated"},
       {cBlue, 22, 1.7, "C: non-tight, isolated"}, {cGreen, 33, 2.3, "D: non-tight, non-isolated"}}};
  std::vector<std::string> labels, titles;
  for (const auto &region : regions) labels.push_back(std::get<3>(region));
  for (const auto &[pkg, c] : order) {
    const std::string centrality = pkg == "PP" ? "" : ", " + percent(lowOf(c), highOf(c));
    titles.push_back((pkg == "PP" ? "p+p" : "Au+Au") + std::string(" DATA, ") + energyOf(pkg) + centrality);
  }
  /* The labels stand at the top right, which the falling yields leave free: in every panel the mark and under it
     the panel's title, as in the IAN, and in the first panel the legend under them. The y range starts just below the
     lowest yield and opens at the top until the yields stand clear of the labels: under each line and under the
     legend, the yields of the bins it reaches over. The title is written at a half-pixel size, which sets a symbol
     and its index tight. */
  const double low = s.positive / 1.6, markHeight = g.title + 3, titleHeight = g.text + 2.5;
  const double between = kGridBlockGap * g.label, clear = 0.8 * g.label;
  const auto [legendWidth, legendHeight] = gridLegendSize(labels, g.text);
  const double markWidth = textWidth(sPhenixInternalLabel(), markHeight) * g.canvasWidth();
  double ceiling = 1.5 * s.high;
  for (size_t p = 0; p < order.size(); ++p) {
    if (!hs[p]) continue;
    /* the widths of the mark, the title and the legend, and the depths their lower edges reach from the frame's top */
    std::vector<std::pair<double, double>> pieces = {
        {markWidth, gridBlock({markHeight})},
        {textWidth(titles[p], titleHeight) * g.canvasWidth(), gridBlock({markHeight, titleHeight})}};
    if (p == 0) pieces.push_back({legendWidth, pieces.back().second + between + legendHeight});
    for (const auto &[wide, depth] : pieces) {
      /* the yields under the piece: from its left end, as a photon E_T, to the right edge of the frame, and half a
         GeV more for the marker of a bin that ends at the piece */
      Span under;
      for (const Bins &b : yields[p])
        under.add(b, true, kPanelEtHigh - (wide + g.indent) / g.width * (kPanelEtHigh - kPanelEtLow) - 0.5,
                  kPanelEtHigh);
      if (!std::isfinite(under.high)) continue;
      ceiling = std::max(ceiling, gridCeiling(p / nx, low, under.high, depth + clear, true));
    }
  }
  for (size_t p = 0; p < order.size(); ++p) {
    gridFrame(p % nx, p / nx, kPanelEtLow, low, kPanelEtHigh, ceiling, true);
    gridAxisTitles("E_{T}^{#gamma,reco} [GeV]", "Event-leading photon candidates / bin", false, true);
    gridLines({{sPhenixInternalLabel(), markHeight}, {titles[p], titleHeight}}, true);
    if (!hs[p]) continue;
    TLegend *l = p == 0 ? gridLegend(labels, true, false, g.reach + between) : nullptr;
    for (int r = 0; r < 4; ++r) {
      const auto &[color, marker, size, text] = regions[r];
      TGraphErrors *points = markers(namedBins(yields[p][r],"abcd/region/"+std::to_string(r)), color, marker, true, true, 0, size);
      photonjet_science::draw(points,"P");
      if (l) l->AddEntry(points, text, "p");
    }
    if (l) l->SetBorderSize(1), photonjet_science::draw(l);
  }
  savePhysicsQA(cv, name);
}

/* The selection lines the purity figures write into the free space of a panel, from the saved selection of a
   package: its acceptance and vertex range ("acceptance"), its tight and non-tight photon identification
   ("identification") or its isolation cone with the isolated and non-isolated regions ("isolation"; I(c) is the
   centrality-dependent threshold). When the selection lacks one of the values of a block, the whole block is left
   out. */
std::vector<std::string> puritySelectionLines(const std::string &pkg, const std::string &what) {
  const auto &m = renderSelection(pkg, true);
  std::vector<std::string> lines;
  if (m.is_null()) return lines;
  const std::string score = apart("s_{BDT}"), cone = apart("E_{T}^{iso}"), vertex = apart("z_{vtx}");
  try {
    if (what == "acceptance") {
      lines.push_back(photonEtaLine(m.at("photon").at("abs_eta_max").get<double>()) +
                      Form(", |%s| < %g cm", vertex.c_str(), m.at("event").at("vertex_abs_z_max_cm").get<double>()));
    } else if (what == "identification") {
      const auto &id = m.at("photon_id");
      if (id.contains("optional_model")) {
        const auto &choice = id.at("optional_model").at("working_point_choice");
        lines.push_back("Tight: " + score + " > " + choice.at("tight").get<std::string>() + "(c)");
        lines.push_back(Form("Non-tight: %g < ", choice.at("non_tight_lower").get<double>()) + score + " < " +
                        choice.at("non_tight_upper").get<std::string>() + "(c)");
      } else if (id.value("model_choice", "") == "ppg12_original") {
        lines.push_back("PPG12 tight / non-tight cuts");
      }
    } else if (what == "isolation") {
      const auto &isolation = m.at("photon_id").at("isolation");
      const std::string system = pkg == "PP" ? "p+p" : "Au+Au";
      lines.push_back(system + Form(" isolation: R = %g", isolation.at("radius").get<double>()));
      if (isolation.at("threshold").is_object() && isolation.at("sideband_gap_gev").is_number()) {
        const double gap = isolation.at("sideband_gap_gev").get<double>();
        lines.push_back("Isolated: " + cone + " < I(c)");
        lines.push_back("Non-isolated: " + cone + " > I(c)" + (gap > 0 ? Form(" + %g GeV", gap) : ""));
      }
    }
  } catch (const PlotJson::exception &) {
    lines.clear();
  }
  return lines;
}

/* Two rows of panels per class on shared axes, every frame on the plain white of the canvas (the purity figures). The
   DATA figure: the raw and the leakage-corrected ABCD purity above, the prompt-signal leakage fractions of the
   photon+jet simulation below. The closure figure: the same purities of the simulation with the truth purity above,
   the leakage-corrected ABCD purity over the truth purity below, around a dashed line at one. Every series is
   drawn with full markers and told from the others
   by its colour: the raw purity black, the truth purity red and the leakage-corrected purity blue, as is its ratio to
   the truth in the lower row of the closure figure. The leakage-corrected purity, the point of both figures, is drawn
   last. The leakage fractions are drawn without error bars. The selection the IAN states in the free space of the
   panels is written from the saved selection: at the top of the lower frames in the DATA figure, under the headings
   in the closure figure. The y ranges given are the least ones: they open further where the labels at the top of a
   frame or a point need the room; an error bar longer than the lower range ends at the frame. Each row's axis title
   says what the row shows. what: "purity" (the DATA figure) or "closure", the folder qa/<package>/<what>/<class>/ the
   points are read from. */
void purityPanels(const std::string &name, const char *what, const char *upperTitle, double upperHigh,
                  const char *lowerTitle, double lowerLow, double lowerHigh) {
  const auto order = panelOrder();
  const bool closure = std::string(what) == "closure";
  /* the legend entries: the purities of the upper row and the leakage fractions of the lower */
  const char *rawLabel = "Raw ABCD", *correctedLabel = "Leakage-corrected ABCD", *truthLabel = "Truth purity";
  const auto leakageLabel = [](int k) { return std::string(Form("S_{%c}/S_{A}", "BCD"[k])); };
  /* The saved points of one class. errors: the saved purity carries statistical errors, which are then drawn.
     ratioErrors: the closure ratio was saved with its errors (otherwise it is the quotient of the saved purities).
     upper, lower: the values of the upper and of the lower row with the error bars drawn there: the labels stay clear
     of them, and the lower range opens to hold the lower ones. */
  struct View {
    bool drawn = false, errors = false, ratioErrors = false;
    Bins raw, corrected, truth, ratio, leakage[3];
    Span upper, lower;
  };
  std::vector<View> views(order.size());
  size_t legendPanel = order.size();  /* the legends sit in the first panel with points */
  for (size_t p = 0; p < order.size(); ++p) {
    const std::string dir = order[p].first + "/" + what + "/" + order[p].second + "/";
    View &v = views[p];
    TH1 *raw = qa(dir + "raw");
    if (raw) v.raw = binsOf(raw);
    v.drawn = raw && anyFinite(v.raw);
    if (!v.drawn) continue;
    legendPanel = std::min(legendPanel, p);
    v.errors = raw->GetSumw2N() > 0;
    if (TH1 *h = qa(dir + "corrected")) v.corrected = binsOf(h);
    v.upper.add(v.raw, v.errors, kPanelEtLow, kPanelEtHigh);
    v.upper.add(v.corrected, v.errors, kPanelEtLow, kPanelEtHigh);
    if (closure) {
      if (TH1 *h = qa(dir + "truth")) v.truth = binsOf(h);
      v.upper.add(v.truth, v.errors, kPanelEtLow, kPanelEtHigh);
      if (TH1 *h = qa(dir + "ratio")) v.ratio = binsOf(h), v.ratioErrors = true;
      else if (v.truth.n() && v.corrected.n()) v.ratio = over(v.corrected, v.truth);
      v.lower.add(v.ratio, v.ratioErrors, kPanelEtLow, kPanelEtHigh);
    } else
      for (int k = 0; k < 3; ++k) {
        if (TH1 *h = qa(dir + "leakage_" + "BCD"[k])) v.leakage[k] = binsOf(h);
        v.lower.add(v.leakage[k], false, kPanelEtLow, kPanelEtHigh);
      }
  }
  if (legendPanel == order.size())
    return missing(name, closure ? "needs both simulation families and the responses"
                                 : "needs DATA analysis fills and the simulation responses");
  TCanvas cv("cv", "", 450 * order.size(), 900);
  const std::string note = closure ? simulationNote(order) : inputNote();
  sharedGrid(cv, order.size(), 2, 26, closure ? std::vector<double>{3, 2} : std::vector<double>{});
  const auto &g = gridLayout();
  /* the legends and the selection lines at a half-pixel size, which sets a symbol and its index tight, and at one
     pitch */
  const double markerSize = 1.5, lineSize = g.text - 2.5, pitch = kGridPitch * lineSize;
  const double between = kGridBlockGap * g.label, clear = between;
  /* per panel: the lines at the top left of its upper frame, the mark in the first panel, then the system and its
     collision energy (on two lines where one is too long for the frame) */
  std::vector<std::vector<std::pair<std::string, double>>> headings(order.size());
  for (size_t p = 0; p < order.size(); ++p) {
    const std::string system = classLabel(order[p].second), energy = energyOf(order[p].first);
    const bool oneLine = textWidth(system + ", " + energy, g.title) * g.canvasWidth() < g.width - 2 * g.indent;
    if (p == 0) headings[p].push_back({sPhenixInternalLabel(), g.title + 2});
    if (oneLine) headings[p].push_back({system + ", " + energy, g.title});
    else headings[p].push_back({system, g.title}), headings[p].push_back({energy, g.text});
  }
  /* The selection lines the IAN writes into the free space of the panels. DATA figure: the p+p cuts and acceptance
     in the p+p panel; the Au+Au acceptance, photon identification and isolation, one block in each Au+Au panel (the
     last panel takes what is left when there are fewer panels than blocks). Closure figure: the Au+Au photon
     identification, in the same panel as in the DATA figure. */
  std::vector<std::vector<std::pair<std::string, double>>> selections(order.size());
  const std::array<std::pair<const char *, const char *>, 3> blocks = {
      {{"acceptance", "Au+Au photon selection"}, {"identification", "Au+Au photon identification"}, {"isolation", ""}}};
  size_t auauPanels = 0;  /* the Au+Au panels */
  for (const auto &panel : order) auauPanels += panel.first != "PP";
  for (size_t p = 0, auau = 0; p < order.size(); ++p) {
    const std::string &pkg = order[p].first;
    std::vector<std::string> lines;
    if (pkg == "PP") {
      if (closure) continue;
      for (const auto &text : puritySelectionLines(pkg, "identification")) lines.push_back("p+p: " + text);
      for (const auto &text : puritySelectionLines(pkg, "acceptance")) lines.push_back(text);
    } else {
      const size_t k = auau++;
      for (size_t j = 0; j < blocks.size(); ++j) {
        if (std::min(j, auauPanels - 1) != k || (closure && j != 1)) continue;
        const auto block = puritySelectionLines(pkg, blocks[j].first);
        if (!block.empty() && *blocks[j].second) lines.push_back(blocks[j].second);
        lines.insert(lines.end(), block.begin(), block.end());
      }
    }
    for (const auto &text : lines) selections[p].push_back({text, lineSize});
  }
  std::vector<std::string> upperLabels = {rawLabel}, lowerLabels;
  if (views[legendPanel].corrected.n()) upperLabels.push_back(correctedLabel);
  if (views[legendPanel].truth.n()) upperLabels.push_back(truthLabel);
  for (int k = 0; k < 3; ++k) if (views[legendPanel].leakage[k].n()) lowerLabels.push_back(leakageLabel(k));
  /* The y ranges: the lower one holds every point with its error bar, and both open at the top until the points and
     their error bars stand clear of the lines and the legend above them. */
  double upperCeiling = upperHigh, lowerFloor = lowerLow, lowerCeiling = lowerHigh;
  for (const View &v : views) {
    if (!v.drawn || !std::isfinite(v.lower.high)) continue;
    if (v.lower.low - 0.03 < lowerFloor) lowerFloor = std::max(0.0, std::floor(10 * (v.lower.low - 0.03)) / 10);
    if (v.lower.high + 0.03 > lowerCeiling) lowerCeiling = std::ceil(10 * (v.lower.high + 0.03)) / 10;
  }
  std::vector<std::pair<std::string, double>> lowerLines;  /* over the legend of the lower row */
  if (!closure) lowerLines.push_back({"Photon+jet simulation", lineSize});
  for (size_t p = 0; p < order.size(); ++p) {
    const View &v = views[p];
    if (!v.drawn) continue;
    const double selection = pitch * selections[p].size();
    std::vector<double> sizes;
    for (const auto &written : headings[p]) sizes.push_back(written.second);
    if (p == 0 && !note.empty()) sizes.push_back(g.text);
    double upper = gridBlock(sizes);
    if (p == legendPanel) upper += between + pitch * upperLabels.size();
    if (closure && selection > 0) upper += between + selection;
    upperCeiling = std::max(upperCeiling, gridCeiling(0, 0, v.upper.high, upper + clear));
    double lower = 0;  /* of the lower frame's labels, from the inset */
    if (p == legendPanel && !lowerLabels.empty()) {
      lower = pitch * lowerLines.size();
      lower += (lower > 0 ? between : 0) + pitch * lowerLabels.size();
    }
    if (!closure && selection > 0) lower += (lower > 0 ? between : 0) + selection;
    if (lower > 0 && std::isfinite(v.lower.high))
      lowerCeiling = std::max(lowerCeiling, gridCeiling(1, lowerFloor, v.lower.high, g.inset + lower + clear));
  }
  /* a range ends on the next multiple of 0.05 (the one below, when that is less than a pixel away): no tick label
     then stands just under the top of a frame */
  for (double *ceiling : {&upperCeiling, &lowerCeiling}) *ceiling = 0.05 * std::ceil(*ceiling / 0.05 - 0.05);
  for (size_t p = 0; p < order.size(); ++p) {
    const View &v = views[p];
    gridFrame(p, 0, kPanelEtLow, 0, kPanelEtHigh, upperCeiling);
    gridAxisTitles("", upperTitle);
    gridLines(headings[p]);
    if (p == 0 && !note.empty()) gridLines({{note, g.text}}, false, false, g.reach, kGray + 2);
    if (v.drawn) {
      TLegend *upperLegend =
          p == legendPanel ? gridLegend(upperLabels, false, false, g.reach + between, lineSize, pitch) : nullptr;
      const char *option = v.errors ? "lep" : "lp";
      TGraphErrors *raw = markers(v.raw, cBlack, 20, v.errors, true, 0, markerSize);
      photonjet_science::draw(raw,"P");
      TGraphErrors *truth = v.truth.n() ? markers(v.truth, cRed, 20, v.errors, true, 0, markerSize) : nullptr;
      if (truth) photonjet_science::draw(truth,"P");
      TGraphErrors *corrected =
          v.corrected.n() ? markers(v.corrected, cBlue, 20, v.errors, true, 0, markerSize) : nullptr;
      if (corrected) photonjet_science::draw(corrected,"P");
      if (upperLegend) {
        upperLegend->AddEntry(raw, rawLabel, option);
        if (corrected) upperLegend->AddEntry(corrected, correctedLabel, option);
        if (truth) upperLegend->AddEntry(truth, truthLabel, option);
        photonjet_science::draw(upperLegend);
      }
    }
    if (closure && !selections[p].empty()) gridLines(selections[p], false, false, g.reach + between);
    gridFrame(p, 1, kPanelEtLow, lowerFloor, kPanelEtHigh, lowerCeiling);
    gridAxisTitles("E_{T}^{#gamma,reco} [GeV]", lowerTitle);
    if (!v.drawn) continue;
    TLegend *lowerLegend = nullptr;
    if (p == legendPanel && !lowerLabels.empty()) {
      gridLines(lowerLines);
      lowerLegend = gridLegend(lowerLabels, false, false, lowerLines.empty() ? 0 : g.reach + between, lineSize, pitch);
    }
    if (!closure && !selections[p].empty()) gridLines(selections[p], false, false, g.reach > 0 ? g.reach + between : 0);
    if (closure) {
      const double one = (1 - lowerFloor) / (lowerCeiling - lowerFloor);
      gridDashedLine(0, one, 1, one);
      if (v.ratio.n()) barredPoints(v.ratio, cBlue, 20, v.ratioErrors, 0, markerSize);
    } else {
      const int colors[] = {cBlack, cRed, cBlue};
      for (int k = 0; k < 3; ++k) {
        if (!v.leakage[k].n()) continue;
        TGraphErrors *leakage = markers(v.leakage[k], colors[k], 20, false, true, 0, markerSize - 0.2);
        photonjet_science::draw(leakage,"P");
        if (lowerLegend) lowerLegend->AddEntry(leakage, leakageLabel(k).c_str(), "lp");
      }
    }
    if (lowerLegend) photonjet_science::draw(lowerLegend);
  }
  savePhysicsQA(cv, name);
}

/* Photon efficiency of the analysis-signal truth photons of the photon simulation against truth pT; one panel per
   class on shared axes. The reconstruction efficiency (black squares), its product with the identification
   efficiency of the reconstructed isolated photons (magenta circles) and the efficiency of the full tight and
   isolated selection (green crosses), in the ABCD photon bins of the qa step as in the IAN: 2 GeV wide up to 25 GeV
   and one bin from 25 to 35 GeV. Each efficiency is a ratio of the sums of the stored 1 GeV bins a bin covers,
   computed here at drawing time; it is drawn with bin-width bars and without a statistical uncertainty. */
void efficiency() {
  const std::string name = "photonOnly/pp_si_di_auau_photon_efficiency_2x2.png";
  const auto order = panelOrder();
  std::vector<TH1 *> hs;
  bool any = false;
  for (const auto &[pkg, c] : order) hs.push_back(qa(pkg + "/phosim/truth/" + c + "/eff")), any = any || hs.back();
  if (!any) return missing(name, "no photon-family simulation QA");
  const auto [nx, ny] = twoColumns(order.size());
  TCanvas cv("cv", "", 550 * nx, 450 * ny);
  const std::string simulation = simulationNote(order);
  sharedGrid(cv, nx, ny, 22);
  const auto &g = gridLayout();
  /* the ABCD photon bins of the qa step (qa.abcd_edges_gev, the setting QA.AbcdEdges): the truth pT is drawn in the
     bins of the reconstructed E_T */
  const std::vector<double> edges = qaAbcdEdges(order.front().first, {15, 17, 19, 21, 23, 25, 35});
  photonjet_science::renderSetting(cv,*gQA,"QA.AbcdEdges",edges);
  const int magenta = TColor::GetColor("#CC2E91");
  const std::vector<std::string> labels = {"#varepsilon_{reco}", "#varepsilon_{reco} #varepsilon_{ID|iso}",
                                           "#varepsilon_{sel}"};
  for (size_t p = 0; p < order.size(); ++p) {
    const std::string &pkg = order[p].first, &c = order[p].second;
    gridFrame(p % nx, p / nx, kPanelEtLow, 0, kPanelEtHigh, 1.1);
    gridAxisTitles("Truth photon p_{T} [GeV]", "Photon efficiency", false, true);
    /* the labels stand on the bottom edge, as in the IAN: the efficiencies leave it free in every system */
    std::vector<std::pair<std::string, double>> lines;
    if (p == 0) lines.push_back({sPhenixInternalLabel(), g.title + 3});
    lines.push_back({classLabel(c) + ", " + energyOf(pkg), g.title});
    lines.push_back({pkg == "PP" ? "Photon+jet simulation" : "Embedded photon+jet simulation", g.text});
    double above = 0;
    if (p == 0 && !simulation.empty()) gridLines({{simulation, g.text}}, false, true, 0, kGray + 2), above = g.reach;
    gridLines(lines, false, true, above);
    if (!hs[p]) continue;
    /* the stages on the second axis of eff: 0 every analysis-signal truth photon, 1 matched to a candidate, 2 and
       preselected, 3 and isolated, 4 preselected and tight, 5 preselected, tight and isolated */
    std::vector<Bins> stage;
    for (int k = 0; k < 6; ++k) stage.push_back(rebinnedTo(project(hs[p], 0, {{1, double(k), k + 1.0}}), edges));
    Bins reco = over(stage[1], stage[0]), idIso = over(stage[5], stage[3]), selected = over(stage[5], stage[0]);
    reco.seriesRole="efficiency/reconstruction";selected.seriesRole="efficiency/full_selection";
    Bins both = reco;
    for (int i = 0; i < both.n(); ++i) both.v[i] = reco.v[i] * idIso.v[i];
    both.physics=photonjet_science::Trace::combine("multiply_efficiencies",{reco.physics,idIso.physics},{{"errors","not displayed; retain source errors"}});
    both.seriesRole="efficiency/reco_times_id_given_iso";
    /* the legend at the lower right of the first panel, in larger type for its symbols and indices */
    TLegend *l = p == 0 ? gridLegend(labels, true, true, 0, g.title - 2) : nullptr;
    const std::tuple<const Bins *, int, int> series[] = {{&reco, cBlack, 21}, {&both, magenta, 20},
                                                         {&selected, cGreen, 34}};
    for (size_t k = 0; k < labels.size(); ++k) {
      const auto &[b, color, marker] = series[k];
      TGraphErrors *points = markers(*b, color, marker, false, true, 0, 1.5);  /* with the bin widths as bars */
      photonjet_science::draw(points,"P");
      if (l) l->AddEntry(points, labels[k].c_str(), "p");
    }
    if (l) photonjet_science::draw(l);
  }
  savePhysicsQA(cv, name);
}

/* ---- event figures --------------------------------------------------------------------------------------------- */

/* The number of neighbouring stored bins that make one cell of a map: the divisor of the bin count that brings the
   cell closest to `size`, so that whole bins are merged and none is left over. */
int cellBins(const std::vector<double> &edges, double size) {
  const int n = edges.size() - 1;
  const double stored = (edges.back() - edges.front()) / n;
  int best = 1;
  for (int k = 2; k <= n; ++k)
    if (n % k == 0 && std::abs(k * stored - size) < std::abs(best * stored - size)) best = k;
  return best;
}

/* Total calorimeter energy against centrality in the embedded photon+jet simulation, before and after the
   minimum-bias classifier requirement: two maps side by side on one energy axis with one logarithmic colour scale,
   as in the IAN, over 0-80% centrality. The cells are as coarse as the IAN's, 5% of centrality by an eighteenth of
   the energy range to the nearest whole number of stored bins, each the sum of the stored bins it covers; the colour
   scale spans every filled cell of the two maps. Each map's title stands on a white box at its top left; the mark
   and the sample stand at the upper right of the second map, which the events passing the classifier leave free. */
void minimumBias() {
  const std::string name = "photonOnly/auau_minbias_pretraining_before_after.png";
  TH1 *h = qa("AuAu/phosim/events/calo_total");
  if (!h) return missing(name, "no Au+Au photon-family simulation QA");
  const auto energy = edgesOf(h->GetXaxis()), centrality = edgesOf(h->GetYaxis());
  const double peripheral = 80;  /* the last centrality shown */
  std::vector<double> energyCells, centralityCells;
  for (size_t i = 0, k = cellBins(energy, (energy.back() - energy.front()) / 18); i < energy.size(); i += k)
    energyCells.push_back(energy[i]);
  for (size_t i = 0, k = cellBins(centrality, 5); i < centrality.size(); i += k)
    centralityCells.push_back(centrality[i]);
  TH2D *maps[2];
  double low = std::numeric_limits<double>::infinity(), top = 0;  /* the smallest and the largest filled cell shown */
  for (int m = 0; m < 2; ++m) {
    maps[m] = new TH2D(fresh().c_str(), "", centralityCells.size() - 1, centralityCells.data(), energyCells.size() - 1,
                       energyCells.data());
    maps[m]->SetDirectory(nullptr);
    for (int e = 1; e < static_cast<int>(energy.size()); ++e)
      for (int c = 1; c < static_cast<int>(centrality.size()); ++c) {
        /* the four bins of the third axis are the minimum-bias states 0 to 3: the classifier's decision plus twice
           the embedded event's own. The second map takes bins 2 and 4, the states in which the classifier passes. */
        double content = 0;
        for (int state = 1; state <= 4; ++state)
          if (m == 0 || state == 2 || state == 4) content += h->GetBinContent(e, c, state);
        if (content != 0) maps[m]->Fill(h->GetYaxis()->GetBinCenter(c), h->GetXaxis()->GetBinCenter(e), content);
      }
    for (int c = 1; c <= maps[m]->GetNbinsX() && maps[m]->GetXaxis()->GetBinUpEdge(c) < peripheral + 1e-9; ++c)
      for (int e = 1; e <= maps[m]->GetNbinsY(); ++e) {
        const double content = maps[m]->GetBinContent(c, e);
        if (content > 0) low = std::min(low, content), top = std::max(top, content);
      }
  }
  if (!(top > 0)) return missing(name, "the calorimeter-energy map of the Au+Au photon simulation has no entries");
  TCanvas cv("cv", "", 1500, 660);
  const std::string simulation = simulationNote({"AuAu"});
  sharedGrid(cv, 2, 1, 24);
  auto &g = gridLayout();
  /* the frames leave room at the right for the colour scale: the bar, its labels and its title */
  const double scaleRoom = 5.6 * g.label;
  g.right += scaleRoom, g.width -= scaleRoom / g.nx;
  const char *titles[] = {"Before: no classifier requirement", "After: minimum-bias pass"};
  for (int m = 0; m < 2; ++m) {
    gridFrame(m, 0, 0, energy.front(), peripheral, energy.back());
    gPad->SetLogz();
    gridAxisTitles("Centrality percentile [%]", "E^{total}_{calo} [GeV]");
    photonjet_science::bindMap(maps[m],"minimum_bias/"+std::to_string(m),photonjet_science::projectionOf(h,PlotJson::array())
      .then("minimum_bias_energy_map",{{"pass_required",m==1},{"passed_state_bins",{2,4}},{"centrality_edges",centralityCells},{"energy_edges",energyCells}}),false);
    maps[m]->SetMinimum(low);
    maps[m]->SetMaximum(top);
    photonjet_science::draw(maps[m],m == 1 ? "COLZ SAME" : "COL SAME");
    gPad->RedrawAxis();
    gridBoxedTitle(titles[m], g.title);
  }
  const double lineSize = 21.5;  /* a half-pixel size, which sets a symbol and its indices tight (see kLinePx) */
  gridLines({{sPhenixInternalLabel(), g.title + 1},
             {"Au+Au embedded photon+jet simulation", lineSize},
             {energyOf("AuAu"), lineSize}},
            true, false, g.reach + kGridBlockGap * g.label);
  if (!simulation.empty()) gridLines({{simulation, lineSize}}, true, false, g.reach, kGray + 2);
  gridColorScale(maps[1], "Weighted events / bin");
  savePhysicsQA(cv, name);
}

/* A simulation family of a reweighting figure: its legend entry, its colour, the marker and the shift inside a bin of
   its ratio, its distribution after the weights and, where the file holds it, before them (no bins: not drawn), and
   the weighted family over the DATA the weights aim at (no bins: the figure has no ratio frame). */
struct ReweightedFamily {
  std::string entry;
  int color, marker;
  double shift;
  Bins before, after, ratio;
  Bins ratioBefore = {};  /* the ratio before the weights, drawn with an open marker; empty: not drawn */
};

/* a / b in the bins of a, each paired with the bin of b that has the same edges (NaN where b has none): the two may
   cover different ranges, as the simulated vertex before the weights does */
Bins ratioOnEdges(const Bins &a, const Bins &b) {
  Bins matched = a;
  for (int i = 0; i < a.n(); ++i) {
    matched.v[i] = 0, matched.e2[i] = 0;
    for (int j = 0; j < b.n(); ++j)
      if (std::fabs(b.edges[j] - a.edges[i]) < 1e-9 && std::fabs(b.edges[j + 1] - a.edges[i + 1]) < 1e-9)
        matched.v[i] = b.v[j], matched.e2[i] = b.e2[j];
  }
  return ratioOf(a, matched);
}

/* An event weight of the Au+Au embedded simulation at work, the one layout of the centrality and of the vertex
   figure. The frame of the distributions: the DATA distribution the weights make the simulation follow, as points,
   and each simulation family before the weights (dashed) and after them (solid), normalised by the caller. Under it,
   when the families carry a `ratio`, a second frame with each weighted simulation over DATA; the families stand a
   little apart inside a bin, so their error bars do not cover each other. Without ratios the figure is the one frame,
   on the same canvas. The colour names the simulation and the line style the stage; where no family has its
   distribution before the weights, the header says that the simulation carries them.
   weights: the name of the weights in the legend ("centrality weights"). peaked: a distribution that is high across
   the range has the header at the top left and the legend at the top right over it; one that falls has both at the
   top right, the legend above the header, as the IAN's centrality figure. Horizontal guide lines, as in the IAN; they
   stop at the text, which stands on white. */
void reweightingClosure(const std::string &name, const Bins &target, const std::vector<ReweightedFamily> &families,
                        double x0, double x1, const char *xtitle, const char *ytitle, const std::string &dataEntry,
                        const std::string &weights, const std::vector<std::string> &lines, bool peaked,
                        const char *ratioTitle = "Weighted / DATA", double cut = -1,
                        const char *function = __builtin_FUNCTION()) {
  bool anyBefore = false, anyRatio = false;
  for (const auto &family : families)
    anyBefore = anyBefore || family.before.n() > 0, anyRatio = anyRatio || family.ratio.n() > 0;
  /* the legend names the two stages when both are drawn; a figure of the weighted simulation alone says so in its
     header instead */
  std::vector<std::string> entries = {dataEntry}, headerLines = lines;
  for (const auto &family : families) entries.push_back(family.entry);
  if (anyBefore) entries.push_back("Before " + weights), entries.push_back("After " + weights);
  else headerLines.push_back("Simulation with its " + weights);
  const double textPx = 30;  /* the axis labels and titles of both frames */
  TCanvas cv("cv", "", 1050, 900);
  const double padLeft = 0.125, padRight = 0.03, split = 0.33;  /* split: the canvas height share of a ratio pad */
  const double topPx = 0.045 * (1 - split) * 900, bottomPx = 0.35 * split * 900;  /* the outer margins, either way */
  TPad *lower = nullptr;
  if (anyRatio) {
    auto pads = ratioPads(cv, padLeft, padRight);
    lower = pads.second;
    pads.first->SetPad(0, split, 1, 1), lower->SetPad(0, 0, 1, split);
    pads.first->SetTopMargin(0.045), pads.first->SetBottomMargin(0.025);
    lower->SetTopMargin(0.04), lower->SetBottomMargin(0.35);
    pads.first->cd();
  } else {
    cv.SetMargin(padLeft, padRight, bottomPx / 900, topPx / 900);
  }
  gPad->SetGridy();
  const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin();
  const double T = gPad->GetTopMargin(), B = gPad->GetBottomMargin();
  const double w = padWidthPx(), h = padHeightPx(), plot = 1 - T - B;
  std::vector<Extent> drawn;
  double labelWidth = 0, peak = 0;
  extents(drawn, target, true, x0, x1);
  for (const auto &family : families) {
    if (family.before.n() > 0) extents(drawn, family.before, false, x0, x1);
    extents(drawn, family.after, false, x0, x1);
  }
  for (const auto &entry : entries) labelWidth = std::max(labelWidth, textWidth(entry, kLinePx));
  for (const auto &e : drawn) peak = std::max(peak, e.y);
  /* the legend at the top right; the header under it with a gap, from the same left edge, or at the top left */
  const double symbol = 46 / w, legendWidth = symbol + labelWidth + 10 / w, right = 1 - R - kInsetPx / w;
  const double ly2 = 1 - T - 14 / h, ly1 = ly2 - entries.size() * kStepPx / h;
  const double legendLeft = peaked ? right - legendWidth : right - std::max(legendWidth, headerWidth(headerLines));
  const double headerLeft = peaked ? L + kInsetPx / w : legendLeft;
  const double mark = peaked ? 1 - T - (14 + kMarkPx) / h : ly1 - (30 + kMarkPx) / h;
  const double headerBottom = mark - (headerLines.size() * kStepPx + 12) / h;
  /* the white fields of the text: their left, right and lower edge */
  std::vector<std::array<double, 3>> blocks;
  if (peaked) {
    /* the header, the legend, and the band between them down to the higher of their lower edges, so that no guide
       line shows between the two */
    blocks.push_back({L, headerLeft + headerWidth(headerLines) + 12 / w, headerBottom});
    blocks.push_back({legendLeft - 12 / w, 1 - R, ly1 - 10 / h});
    blocks.push_back({L, 1 - R, std::max(headerBottom, ly1 - 10 / h)});
  } else {
    blocks.push_back({legendLeft - 12 / w, 1 - R, headerBottom});
  }
  const auto userX = [&](double ndc) { return x0 + (ndc - L) / (1 - L - R) * (x1 - x0); };
  /* the frame starts just below zero, so a point of a nearly empty bin is drawn whole */
  const double low = -0.03 * peak;
  double high = low + 1.1 * (peak - low);
  for (const auto &[left, edge, bottom] : blocks)
    high = std::max(high, topBelow(drawn, low, false, userX(left), userX(edge), (bottom - B) / plot));
  TH1F *f = photonjet_science::drawFrame(*gPad,x0, low, x1, high);
  f->GetYaxis()->SetTitle(ytitle);
  pixelAxes(f, textPx, textPx);
  if (anyRatio) {
    hideXText(f);
  } else {
    f->GetXaxis()->SetTitle(xtitle);
    f->GetXaxis()->SetTitleOffset(1.45);
    f->GetXaxis()->CenterTitle();
  }
  f->GetYaxis()->SetTitleOffset(2.1);
  f->GetYaxis()->CenterTitle();
  tickLengths(f, 14);
  const auto userY = [&](double ndc) { return low + (ndc - B) / plot * (high - low); };
  TBox backing;
  backing.SetFillColor(kWhite);
  backing.SetLineColor(kWhite);
  for (const auto &[left, edge, bottom] : blocks)
    backing.DrawBox(userX(left + 2 / w), userY(bottom + 2 / h), userX(edge), userY(ly2 + 6 / h));
  TLegend *leg = ianLegend(legendLeft, ly1, legendLeft + legendWidth, ly2, symbol / legendWidth);
  /* the analysis cut, dotted, over the whole frame (cut > 0 only) */
  const auto cutLines = [&](double y0, double y1) {
    if (!(cut > 0)) return;
    for (double x : {-cut, cut}) {
      auto *l = new TLine(x, y0, x, y1);
      l->SetLineStyle(kDotted), l->SetLineColor(kGray + 2), l->SetLineWidth(2);
      photonjet_science::draw(l);
    }
  };
  cutLines(low, high);
  /* the simulation as steps over the bins, with no drop to zero at the ends of the range */
  std::vector<TGraph *> afterLines;
  for (const auto &family : families) {
    if (family.before.n() > 0)
      photonjet_science::draw(outline(line(namedBins(family.before,"reweighting/before"), family.color), x0, x1, family.color, kDashed),"L");
    afterLines.push_back(outline(line(namedBins(family.after,"reweighting/after"), family.color), x0, x1, family.color, kSolid));
  }
  for (auto *afterLine : afterLines) photonjet_science::draw(afterLine,"L");
  TGraphErrors *points = markers(namedBins(populated(target),"reweighting/target"), cBlack, 20, true, false, 0, 1.3);
  photonjet_science::draw(points,"P");  /* DATA last, over the simulation */
  /* the ratios, on one range around 1 that holds every family */
  std::vector<TGraphErrors *> ratios;
  double ratioLow = 1, ratioHigh = 1;
  bool anyRatioBefore = false;
  for (const auto &family : families) {
    if (family.ratioBefore.n() > 0) {   /* open marker, a little to the other side inside the bin */
      anyRatioBefore = true;
      /* ROOT styles 72/73 retain the open square/triangle with a three-pixel outline. */
      ratios.push_back(markers(namedBins(family.ratioBefore,"reweighting/ratio_before"), family.color, family.marker + 51, true, false, -family.shift, 1.4));
      const auto [lowOne, highOne] = rangeAroundOne(ratios.back(), 1.0);
      ratioLow = std::min(ratioLow, lowOne), ratioHigh = std::max(ratioHigh, highOne);
    }
    if (family.ratio.n() == 0) continue;
    ratios.push_back(markers(namedBins(family.ratio,"reweighting/ratio_after"), family.color, family.marker, true, false, family.shift, 1.4));
    const auto [lowOne, highOne] = rangeAroundOne(ratios.back(), 1.0);
    ratioLow = std::min(ratioLow, lowOne), ratioHigh = std::max(ratioHigh, highOne);
  }
  /* a family's entry shows its line and, with a ratio frame, the marker of its ratio; the stages are named by a
     black line each */
  leg->AddEntry(points, entries[0].c_str(), "p");
  for (size_t i = 0; i < families.size(); ++i) {
    auto *key = new TGraph;
    style(key, families[i].color, families[i].marker, 1.3);
    key->SetLineWidth(afterLines[i]->GetLineWidth());
    leg->AddEntry(key, entries[1 + i].c_str(), anyRatio ? "lp" : "l");
  }
  for (int stage : anyBefore ? std::vector<int>{kDashed, kSolid} : std::vector<int>{}) {
    auto *key = new TGraph;
    key->SetLineColor(cBlack);
    key->SetLineStyle(stage);
    key->SetLineWidth(afterLines.front()->GetLineWidth());
    if (anyRatioBefore) key->SetMarkerStyle(stage == kDashed ? 71 : 20), key->SetMarkerColor(cBlack), key->SetMarkerSize(1.3);
    leg->AddEntry(key, entries[entries.size() - (stage == kDashed ? 2 : 1)].c_str(), anyRatioBefore ? "lp" : "l");
  }
  photonjet_science::draw(leg);
  markAndLines(headerLeft, mark, headerLines, kMarkPx, kLinePx, kStepPx);
  gPad->RedrawAxis();
  if (anyRatio) {
    lower->cd();
    TH1F *rf = ratioFrame(x0, x1, ratioLow, ratioHigh, xtitle, ratioTitle);
    pixelAxes(rf, textPx, textPx);
    rf->GetXaxis()->SetTitleOffset(1.45);
    rf->GetYaxis()->SetTitleOffset(2.1);
    rf->GetYaxis()->SetNdivisions(505);
    rf->GetYaxis()->SetDecimals(true);
    tickLengths(rf, 14);
    lowerTitles(rf);
    frameReferenceLine(1);
    cutLines(ratioLow, ratioHigh);
    for (auto it = ratios.rbegin(); it != ratios.rend(); ++it) photonjet_science::draw((*it),"P");
  }
  /* light solid guide lines at the y ticks of the upper frame; the style's own grid look is put back after saving */
  const int gridColor = gStyle->GetGridColor(), gridStyle = gStyle->GetGridStyle();
  gStyle->SetGridColor(kGray);
  gStyle->SetGridStyle(1);
  savePhysicsQA(cv, name, function);
  gStyle->SetGridColor(gridColor);
  gStyle->SetGridStyle(gridStyle);
}

/* The trigger and the vertex window of the Au+Au package as one header line, the trigger named as the IAN does, by
   its photon and MBD parts. */
std::string triggerAndVertexLine() {
  const auto &selection = renderSelection("AuAu", true);
  std::string trigger;
  if (!selection.is_null()) trigger = selection.value("event", PlotJson::object()).value("trigger_name", trigger);
  if (trigger.find(" NS") != std::string::npos) trigger.resize(trigger.find(" NS"));
  const std::string vertex = Form("|z_{vtx}| < %g cm", selectionValue("AuAu", "event", "vertex_abs_z_max_cm", 10));
  return trigger.empty() ? vertex : latexCuts(trigger) + ", " + vertex;
}

/* The centrality weights of the Au+Au embedded simulation, in the 5% slices they are defined in (reweightingClosure).
   DATA: the centrality of the tight, isolated photon candidates in the photon E_T range of the selection, the
   population the weights make the simulated photons follow. Simulation: the candidates of each family under the same
   selection, the prompt photons of the photon+jet samples and the background candidates of the inclusive-jet samples,
   the population each family gives the analysis. The saved simulation carries the weights. The weight is one number
   per slice, the ratio of the two weight rows of a family's simulated events' centrality, so a slice of the weighted
   candidates divided by it is that slice before the weights. */
void centralityReweighting() {
  const std::string name = "simulation_stitching/auau_centrality_data_embeddedMC.png";
  /* every tight, isolated candidate fills each saved width histogram of its slice once: one of them counts them */
  const std::string counted = "width_tight_iso_weta_cogx";
  TH1 *data = qa("AuAu/data/slices/data/" + counted);
  const double shownEnd = 80;  /* percent: the end of the centrality range shown */
  const auto et = photonSupport("AuAu");
  /* the candidates of the photon E_T range per slice; the third axis numbers the slices, which become percent */
  const auto perSlice = [&](const TH1 *h) {
    Bins b = project(h, 2, {{1, et.first, et.second}});
    for (double &edge : b.edges) edge *= kSlicePercent;
    return within(b, 0, shownEnd);
  };
  /* a family: its folder, the folder under slices/ of its candidates, and its look in the figure */
  std::vector<ReweightedFamily> families = {{"Embedded photon+jet", cBlue, 21, -0.7, {}, {}},
                                            {"Embedded inclusive jet", cOrange, 22, 0.7, {}, {}}};
  const char *folders[] = {"AuAu/phosim", "AuAu/jetsim"}, *populations[] = {"sig", "bkg"};
  for (size_t i = 0; i < families.size(); ++i) {
    const std::string folder = folders[i];
    TH1 *candidates = qa(folder + "/slices/" + populations[i] + "/" + counted), *events = qa(folder + "/events/centrality");
    if (!data || !candidates || !events)
      return missing(name, "needs the tight, isolated candidates per centrality slice of Au+Au DATA and of both "
                           "embedded simulations, and the centrality of the simulated events");
    /* the second axis of the simulated events' centrality is the weight: [0, 1) the event weight alone, [1, 2) with
       the centrality weight */
    const auto eventsPerSlice = [&](double row) {
      return within(rebinToWidth(project(events, 0, {{1, row, row + 1}}), kSlicePercent), 0, shownEnd);
    };
    const Bins weighted = perSlice(candidates), weight = over(eventsPerSlice(1), eventsPerSlice(0));
    if (weight.n() != weighted.n())
      return missing(name, "the centrality slices of the simulated events and of the candidates differ");
    Bins unweighted = weighted;
    for (int k = 0; k < weighted.n(); ++k) {
      const double w = weight.v[k] > 0 ? weight.v[k] : std::numeric_limits<double>::infinity();
      unweighted.v[k] /= w, unweighted.e2[k] /= w * w;
    }
    families[i].before = unit(unweighted), families[i].after = unit(weighted);
  }
  const Bins target = unit(perSlice(data));
  for (auto &family : families) family.ratio = ratioOf(family.after, target);
  reweightingClosure(name, target, families, 0, shownEnd, "Centrality percentile [%]",
                     "Fraction of candidates per 5% bin", "DATA, selected photon candidates", "centrality weights",
                     {"Au+Au, " + energyOf("AuAu"), triggerAndVertexLine(), "Tight, isolated photon candidates",
                      photonWindowLabel(et.first, et.second), "Weight target: the DATA candidates"},
                     false);
}

/* The vertex weights of the Au+Au embedded simulation, the centrality figure's counterpart (reweightingClosure), as
   one frame over the whole recorded vertex range in 2 cm bins: the reconstructed vertex of the DATA events of the
   classes and of the simulated events of each family with their full weight, every distribution with unit area
   inside the vertex window, where the simulation lives. DATA: the events that pass the event selection but for its
   vertex window (events/vertex_selected_wide), the population the weights make the simulation follow, where the file
   holds it; else every recorded event of the classes (their evt_vz, summed), whose shape inside the window differs
   from that of the selected events by a few percent. Before the vertex weights: the row [0, 1) of events/vertex,
   the same simulated events with the full weight without its vertex factor, where the file holds it. */
void vertexReweighting() {
  const std::string name = "simulation_stitching/auau_vertex_data_embeddedMC.png";
  TH1 *data = qa("AuAu/data/events/vertex_selected_wide");
  if (!data) data = allEvents("AuAu", "data", "evt_vz");
  if (!data) return missing(name, "needs the vertex of the Au+Au DATA events");
  /* the analysis vertex window only (Justin, Oct 15, on this panel: "This should just be -10 to 10 cm for
     simplicity") */
  const double window = selectionValue("AuAu", "event", "vertex_abs_z_max_cm", 10), shownEnd = window, binCm = 2.0;
  /* a saved vertex distribution in the bins shown, with unit area inside the vertex window */
  const auto shown = [&](const TH1 *h) {
    const Bins b = within(rebinToWidth(project(h, 0), binCm), -shownEnd, shownEnd);
    const double inside = within(b, -window, window).sum();
    return scaled(b, inside > 0 ? 1 / inside : 0.0);
  };
  std::vector<ReweightedFamily> families = {{"PYTHIA 8 #gamma+jet embedded", cBlue, 21, 0.25, {}, {}, {}},
                                            {"PYTHIA 8 jet embedded", cOrange, 22, -0.25, {}, {}, {}}};
  const char *names[] = {"phosim", "jetsim"};
  for (size_t i = 0; i < families.size(); ++i) {
    TH1 *weighted = allEvents("AuAu", names[i], "evt_vz");
    if (!weighted) return missing(name, "needs the vertex of the simulated events of every centrality class");
    families[i].after = shown(weighted);
    if (TH1 *rows = qa(std::string("AuAu/") + names[i] + "/events/vertex")) {
      const Bins b = within(rebinToWidth(project(rows, 0, {{1, 0, 1}}), binCm), -shownEnd, shownEnd);
      const double inside = within(b, -window, window).sum();   /* unit area inside the window, as DATA */
      families[i].before = scaled(b, inside > 0 ? 1 / inside : 0.0);
    }
  }
  /* DATA over each simulation, before and after the weights (the after one is the validation), inside the analysis
     window, where both are normalised; outside it the DATA (before the event selection) and the generation walls of
     the simulation make the ratio meaningless */
  const Bins target = shown(data);
  for (auto &family : families) {
    family.ratio = within(ratioOnEdges(target, family.after), -window, window);
    if (family.before.n() > 0) family.ratioBefore = within(ratioOnEdges(target, family.before), -window, window);
  }
  const auto classes = classesOf("AuAu");
  double low = lowOf(classes.front()), high = highOf(classes.front());
  for (const auto &c : classes) low = std::min(low, lowOf(c)), high = std::max(high, highOf(c));
  reweightingClosure(name, target, families, -shownEnd, shownEnd, "Vertex z [cm]",
                     Form("Fraction of events per %g cm bin", binCm), "DATA", "vertex weights",
                     {"Au+Au " + rangeDashes(percent(low, high)) + ", " + energyOf("AuAu"),
                      Form("Unit area inside |z_{vtx}| < %g cm", window)},
                     true, "DATA / Simulation", window);
}

/* One distribution of the event-property figures in the IAN's arrangement: the sPHENIX mark, the system line and the
   legend (DATA first) as one block at the top right, where a falling distribution is low; for a peaked distribution
   (split) the header at the top left and the legend at the top right, over its two wings. The y range keeps every
   drawn value clear of them. `corner` lines are written at the lower left, under a falling spectrum. Simulation as
   lines over the bins they fill, DATA as points; secondLine: the line style of the second simulation (dashed where
   the two lie on each other, so both show). */
void eventDistribution(const std::vector<Curve> &curves, double x0, double x1, bool logy, const char *xtitle,
                       const char *ytitle, const std::vector<std::string> &lines,
                       const std::vector<std::string> &corner, const std::string &name, bool split = false,
                       int secondLine = 1) {
  TCanvas cv("cv", "", 700, 600);
  gPad->SetLogy(logy);
  const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin();
  const double T = gPad->GetTopMargin(), B = gPad->GetBottomMargin();
  const double w = padWidthPx(), h = padHeightPx(), plot = 1 - T - B;
  /* a line is drawn from its first to its last filled bin; a DATA bin without entries has no point */
  const auto shownBins = [](const Curve &c) {
    if (c.points) return populated(c.bins);
    const Bins &b = c.bins;
    int first = 0, last = b.n() - 1;
    while (first < last && b.v[first] == 0) ++first;
    while (last > first && b.v[last] == 0) --last;
    return b.n() ? within(b, b.edges[first], b.edges[last + 1]) : b;
  };
  Span s;
  std::vector<Extent> drawn, lows;
  std::vector<Bins> shown;
  double labelWidth = 0;
  for (const auto &c : curves) {
    shown.push_back(shownBins(c));
    const Bins &b = shown.back();
    const bool bars = c.points;  /* points carry their error bars, lines none */
    s.add(b, bars, x0, x1), extents(drawn, b, bars, x0, x1);
    for (int i = 0; i < b.n(); ++i)
      if (std::isfinite(b.v[i]) && b.v[i] != 0)
        lows.push_back({b.edges[i], b.edges[i + 1], b.v[i] - (bars ? std::sqrt(b.e2[i]) : 0.0)});
    labelWidth = std::max(labelWidth, textWidth(c.text, kLinePx));
  }
  const std::string note = inputNote();
  const double symbol = 46 / w, legendWidth = symbol + labelWidth + 10 / w, right = 1 - R - kInsetPx / w;
  const double headWidth = std::max(headerWidth(lines), note.empty() ? 0.0 : textWidth(note, kLinePx));
  const double hx = split ? L + kInsetPx / w : right - std::max(headWidth, legendWidth);
  const double top = 1 - T - (14 + kMarkPx) / h, headBottom = header(hx, top, lines, false);
  const double lx1 = split ? right - legendWidth : hx, ly2 = split ? 1 - T - 12 / h : headBottom - 10 / h;
  const double ly1 = ly2 - curves.size() * kStepPx / h;
  const auto user = [&](double ndc) { return x0 + (ndc - L) / (1 - L - R) * (x1 - x0); };
  std::vector<std::array<double, 3>> blocks;  /* from, to (user x), the lower edge as a fraction of the frame height */
  blocks.push_back({user(lx1 - 10 / w), x1, (ly1 - 10 / h - B) / plot});
  if (split) blocks.push_back({x0, user(hx + headWidth + 10 / w), (headBottom - 12 / h - B) / plot});
  double low = 0;
  if (logy) low = std::isfinite(s.positive) ? s.positive / 2 : 1e-3;
  else if (std::isfinite(s.low)) low = std::min(0.0, 1.08 * s.low);
  /* on a linear scale the highest value stays a fifth of the frame under its top, as in the IAN */
  const auto frameTop = [&] {
    double high = topBelow(drawn, low, logy, x0, x1, logy ? 0.95 : 0.8);
    for (const auto &b : blocks) high = std::max(high, topBelow(drawn, low, logy, b[0], b[1], b[2]));
    return high > low ? high : low + 1;
  };
  double high = frameTop();
  /* The lines at the lower left start right of the edge with which a simulated line rises inside the frame. The
     spectrum stays above them: the frame is opened downwards for them where needed. */
  double cornerLeft = L + kInsetPx / w, cornerWidth = 0;
  for (size_t i = 0; i < curves.size(); ++i)
    if (!curves[i].points && shown[i].n() && shown[i].edges.front() < x0 + 0.25 * (x1 - x0))
      cornerLeft = std::max(cornerLeft, L + (shown[i].edges.front() - x0) / (x1 - x0) * (1 - L - R) + 12 / w);
  for (const auto &text : corner) cornerWidth = std::max(cornerWidth, textWidth(text, kLinePx));
  if (!corner.empty()) {
    const double needed =
        bottomAbove(lows, high, logy, user(cornerLeft - 6 / w), user(cornerLeft + cornerWidth + 10 / w),
                    (14 + corner.size() * kStepPx) / h / plot);
    if (needed < low) low = needed, high = frameTop();
  }
  TH1F *f = photonjet_science::drawFrame(*gPad,x0, low, x1, high);
  f->GetXaxis()->SetTitle(xtitle);
  f->GetYaxis()->SetTitle(ytitle);
  f->GetYaxis()->SetMaxDigits(4);  /* counts are labelled with one to three digits in units of 10^3, 10^6, ... */
  header(hx, top, lines);
  TLatex t;
  t.SetNDC();
  t.SetTextFont(43);
  t.SetTextSize(kLinePx);
  for (size_t i = 0; i < corner.size(); ++i)
    t.DrawLatex(cornerLeft, B + (14 + (corner.size() - 1 - i) * kStepPx) / h, corner[i].c_str());
  TLegend *leg = ianLegend(lx1, ly1, lx1 + legendWidth, ly2, symbol / legendWidth);
  std::vector<TObject *> objects(curves.size(), nullptr);
  int drawnLines = 0;
  for (size_t i = 0; i < curves.size(); ++i)
    if (!curves[i].points) {
      TH1D *hist = line(shown[i], curves[i].color, drawnLines++ ? secondLine : 1);
      /* a line that starts on the frame's left edge is drawn without its end strokes, so the frame keeps its colour;
         where it ends inside the frame its fall to the axis is drawn again */
      const bool onEdge = shown[i].n() && std::fabs(shown[i].edges.front() - x0) < 1e-9 * (x1 - x0);
      photonjet_science::draw(hist,onEdge ? "HIST SAME ][" : "HIST SAME");
      if (onEdge && shown[i].edges.back() < x1 - 1e-9 * (x1 - x0)) {
        TLine fall;
        fall.SetLineColor(hist->GetLineColor());
        fall.SetLineStyle(hist->GetLineStyle());
        fall.SetLineWidth(hist->GetLineWidth());
        fall.DrawLine(shown[i].edges.back(), low, shown[i].edges.back(), shown[i].v.back());
      }
      objects[i] = hist;
    }
  for (size_t i = 0; i < curves.size(); ++i)
    if (curves[i].points) {
      const Curve &c = curves[i];
      TGraphErrors *g = markers(shown[i], c.color, 20, true, c.xerrors);
      photonjet_science::draw(g,"P");
      objects[i] = g;
    }
  for (bool points : {true, false})
    for (size_t i = 0; i < curves.size(); ++i)
      if (curves[i].points == points) leg->AddEntry(objects[i], curves[i].text.c_str(), points ? "p" : "l");
  photonjet_science::draw(leg);
  savePhysicsQA(cv, name);
}

/* An eta-phi map of the event-property figures in the IAN's format. The frame is the layout of the detector on every
   map, |eta| < 1.1 by -pi < phi < pi, whatever part of it the cells cover: the cluster and the jet maps of DATA and of
   the simulation then read against each other, and against the IAN's, place by place. The colour scale stands at the
   right of the frame with its title. The header (three lines at most) stands above the frame, where it hides no cell:
   the sPHENIX mark over the frame's left edge with the first line, the sample, ending at the right edge of the canvas;
   under them the second line at the left and the third at the right. A cell without entries inside the map takes the
   colour of zero, as in the IAN. The colour scale starts at zero. It ends at the highest cell, as the IAN's does, or,
   with midScale, at twice the median of the filled cells: a map whose cells lie at one level then shows that level in
   the middle of the scale, and the DATA and the simulated map of one quantity read alike. */
void etaPhiMap(TH2 *m, const char *ztitle, const std::vector<std::string> &lines, const std::string &name,
               bool midScale) {
  const auto written = nonEmpty(lines);
  if (written.size() > 3) throw std::runtime_error("the header of an eta-phi map holds three lines at most: " + name);
  /* One canvas for every map of the family. The frame has the proportions of the IAN's; the right margin holds the
     colour scale, its labels of up to four characters and its title. The top margin holds the two rows of the header,
     placed as the header of a response map is: the lower row kHeaderTopPx above the frame; the upper row a line pitch
     higher and, as under every mark, 0.15 of a pitch more; over it the height of the mark and a clear strip. */
  TCanvas cv("cv", "", 780, 664);
  const double rowRisePx = 1.15 * kStepPx, headerRoomPx = kHeaderTopPx + rowRisePx + kMarkPx + kLegendTopPx;
  cv.SetRightMargin(0.23);
  cv.SetTopMargin(headerRoomPx / padHeightPx());
  /* ROOT leaves a cell of exactly zero unpainted: an empty cell takes a value far below one entry instead */
  if(photonjet_science::tag(m,"physics_series").is_null())photonjet_science::bindMap(m,"eta_phi",photonjet_science::projectionOf(m,PlotJson::array()),true);
  const double trace = 1e-9 * m->GetMaximum();
  std::vector<double> filled;
  for (int i = 1; i <= m->GetNbinsX(); ++i)
    for (int j = 1; j <= m->GetNbinsY(); ++j) {
      if (m->GetBinContent(i, j) > 0) filled.push_back(m->GetBinContent(i, j));
      if (m->GetBinContent(i, j) == 0) m->SetBinContent(i, j, trace);
    }
  photonjet_science::copyBinding(m,m); // zero-cell paint floor leaves the measured zero in the binding
  m->SetMinimum(0);
  if (midScale && !filled.empty()) {
    std::nth_element(filled.begin(), filled.begin() + filled.size() / 2, filled.end());
    m->SetMaximum(2 * filled[filled.size() / 2]);
  }
  const double etaEdge = 1.1;  /* the pseudorapidity the calorimeters reach */
  TH1F *f = photonjet_science::drawFrame(*gPad,-etaEdge, -M_PI, etaEdge, M_PI);
  f->GetXaxis()->SetTitle("#eta");
  f->GetYaxis()->SetTitle("#phi [rad]");
  f->GetXaxis()->SetNdivisions(505);  /* a label at every half unit of eta: closer labels would run into each other */
  tickLengths(f, kTickPx);
  m->GetZaxis()->SetTitle(ztitle);
  m->GetZaxis()->SetTitleOffset(1.65);
  photonjet_science::draw(m,"COLZ SAME");
  gPad->RedrawAxis();
  /* the two rows of the header; a line at the right ends as far from the edge of the canvas as the lower row stands
     above the frame */
  const double w = padWidthPx(), h = padHeightPx();
  const double left = gPad->GetLeftMargin(), right = 1 - kHeaderTopPx / w;
  const double lower = 1 - gPad->GetTopMargin() + kHeaderTopPx / h, upper = lower + rowRisePx / h;
  TLatex t;
  t.SetNDC();
  t.SetTextFont(kFont);
  t.SetTextSize(kMarkPx);
  t.DrawLatex(left, upper, sPhenixInternalLabel());
  t.SetTextSize(kLinePx);
  for (size_t i = 0; i < written.size(); ++i) {
    const bool atRight = i != 1;
    t.SetTextAlign(atRight ? 31 : 11);
    t.DrawLatex(atRight ? right : left, i == 0 ? upper : lower, written[i].c_str());
  }
  inputNote();  /* the state of the inputs goes into the receipt of the figure; the canvas carries no such line */
  savePhysicsQA(cv, name);
}

/* Au+Au event properties, all classes together, under event_obs/: the centrality in the analysis classes
   (centrality.pdf), the vertex (vertexz.pdf), the cluster E_T (cluster_pt.pdf) and the jet pT (jet_pt.pdf) of DATA
   against the photon and the inclusive-jet simulation, and the eta-phi maps of the clusters and the jets for DATA and
   for the photon simulation (cluster_eta_phi_<data|photon>.pdf, jet_eta_phi_<data|photon>.pdf). Each simulated
   distribution is scaled to the DATA over the range shown, the vertex inside the analysis vertex window; the maps are
   sums without a scale. The colours and the legend wording are the IAN's. The distributions are drawn by
   eventDistribution and the maps by etaPhiMap, whose names the receipts of these figures carry. */
void eventProperties() {
  const auto classes = classesOf("AuAu");
  std::string system = "Au+Au";
  if (!classes.empty()) system += " " + centralityRange(lowOf(classes.front()), highOf(classes.back()));
  const std::string energy = energyOf("AuAu"), span = system + ", " + energy;
  struct Look { std::string family; int color; std::string text; };
  const std::vector<Look> sims = {{"phosim", cBlue, "Pythia #gamma Overlay"}, {"jetsim", cRed, "Pythia jet Overlay"}};
  /* centrality, in the analysis classes as in the IAN: the DATA events of the analysis trigger
     (qa/AuAu/centrality_target; data/events/centrality holds every recorded event, before the trigger selection)
     against the simulation without its centrality weights */
  if (TH1 *data = qa("AuAu/centrality_target")) {
    std::vector<double> edges;
    for (const auto &c : classes) edges.push_back(lowOf(c));
    edges.push_back(classes.empty() ? 80 : highOf(classes.back()));
    if (classes.empty()) edges.insert(edges.begin(), 0);
    const Bins d = rebinnedTo(binsOf(data), edges);
    std::vector<Curve> curves;
    for (const auto &s : sims)
      if (TH1 *h = qa("AuAu/" + s.family + "/events/centrality")) {
        /* the second axis is the weight mode: bin 1 the event weight alone, bin 2 with the centrality weight */
        const TAxis *mode = h->GetYaxis();
        if (h->GetDimension() != 2 || mode->GetNbins() != 2 || mode->GetXmin() != 0 || mode->GetXmax() != 2)
          throw std::runtime_error("simulation centrality QA has no weight-mode axis: " + s.family);
        const Bins b = rebinnedTo(project(h, 0, {{1, 0, 1}}), edges);
        curves.push_back(lineCurve(scaledTo(b, d.sum()), s.color, s.text));
      }
    curves.push_back(dataCurve(d, "Data", true));
    eventDistribution(curves, edges.front(), edges.back(), false, "Centrality [%]", "Events", {span}, {},
                      "event_obs/centrality.pdf", false, 2);
  } else missing("event_obs/centrality.pdf", "no Au+Au DATA QA");
  /* vertex, over the whole range of the saved histograms in 2 cm bins */
  if (TH1 *data = allEvents("AuAu", "data", "evt_vz")) {
    /* The DATA fill covers the vertex before the event selection, the simulation only the selected range: the
       simulation is scaled to the DATA inside the selected range, so the shapes compare there. */
    const double zmax = selectionValue("AuAu", "event", "vertex_abs_z_max_cm", 10);
    const double binCm = 2.0;
    const std::string perBin = Form("Events / %g cm", binCm);
    const Bins d = binsOf(data);
    const double inside = within(d, -zmax, zmax).sum();
    std::vector<Curve> curves;
    for (const auto &s : sims)
      if (TH1 *h = allEvents("AuAu", s.family, "evt_vz")) {
        const Bins b = binsOf(h);
        const double own = within(b, -zmax, zmax).sum();
        curves.push_back(lineCurve(rebinToWidth(scaled(b, own > 0 ? inside / own : 0.0), binCm), s.color, s.text));
      }
    curves.push_back(dataCurve(rebinToWidth(d, binCm), "Data"));
    eventDistribution(curves, d.edges.front(), d.edges.back(), false, "Vertex z [cm]", perBin.c_str(), {span}, {},
                      "event_obs/vertexz.pdf", true, 2);
  } else missing("event_obs/vertexz.pdf", "no Au+Au DATA QA");
  /* The Au+Au DATA trees hold photon candidates from 12 GeV and the simulation from lower, so the cluster spectra are
     shown and scaled to the DATA from 12 GeV. A cluster map holds the candidates above the threshold of its fill
     (10 GeV in a package that does not record it); the DATA map has none below the threshold of its trees. */
  const double treeEtMin = 12;
  const auto &saved = renderSelection("AuAu", true);
  const bool recorded = !saved.is_null() && saved.at("qa").contains("cluster_map_et_min_gev");
  const double mapEtMin = recorded ? saved.at("qa").at("cluster_map_et_min_gev").get<double>() : 10;
  /* cluster ET of the preselected candidates, up to the end of the IAN's axis */
  const Cut preselected = {1, 1, 2};  /* of the stage axis of clu_et */
  if (TH1 *data = allCandidates("AuAu", "data", "clu_et")) {
    const double from = treeEtMin, to = 50;
    const Bins d = within(project(data, 0, {preselected}), from, to);
    std::vector<Curve> curves;
    for (const auto &s : sims)
      if (TH1 *h = allCandidates("AuAu", s.family, "clu_et")) {
        const Bins b = within(project(h, 0, {preselected}), from, to);
        curves.push_back(lineCurve(scaledTo(b, d.sum()), s.color, s.text));
      }
    curves.push_back(dataCurve(d, "Data"));
    eventDistribution(curves, from, to, true, "Cluster E_{T} [GeV]", "Clusters", {}, {span, "Preselected candidates"},
                      "event_obs/cluster_pt.pdf");
  } else missing("event_obs/cluster_pt.pdf", "no Au+Au DATA QA");
  /* cluster and jet maps, their cells merged in whole bins to the IAN's granularity; the jet map over the jet
     acceptance of the analysis, as the cluster map is over the photon acceptance, and both on the frame of the whole
     detector (etaPhiMap). The jet maps keep the IAN's colour scale, from zero to their highest cell. The cluster maps
     end theirs at twice their median cell: the simulated clusters fill every cell of the acceptance to one level,
     which then stands in the middle of the scale as the bulk of the DATA map does. */
  const double jetPtMin = selectionValue("AuAu", "qa", "jet_map_pt_min_gev", 5);
  const double jetEtaMax = selectionValue("AuAu", "jet", "abs_eta_max", 0.7);
  /* the cells of a map merged in place to about etaStep by phiStep, in whole stored bins with none left over */
  const auto merged = [](TH2D *map, double etaStep, double phiStep) {
    const auto group = [](const TAxis *a, double step) {
      int k = std::max(1, static_cast<int>(std::lround(step / a->GetBinWidth(1))));
      while (a->GetNbins() % k) --k;
      return k;
    };
    const int etaGroup=group(map->GetXaxis(),etaStep),phiGroup=group(map->GetYaxis(),phiStep);
    const auto trace=photonjet_science::projectionOf(map,PlotJson::array()).then("rebin_map",{{"eta_group",etaGroup},{"phi_group",phiGroup}});
    map->Rebin2D(etaGroup,phiGroup);
    return photonjet_science::bindMap(map,"eta_phi_map",trace,false);
  };
  /* the whole eta bins of a map inside |eta| < limit, as a new map (the map itself when no whole bin lies inside) */
  const auto central = [](TH2D *map, double limit) {
    const TAxis *eta = map->GetXaxis();
    std::vector<double> edges;
    int firstBin = 0;
    for (int i = 1; i <= eta->GetNbins(); ++i)
      if (eta->GetBinLowEdge(i) >= -limit - 1e-9 && eta->GetBinUpEdge(i) <= limit + 1e-9) {
        if (edges.empty()) edges.push_back(eta->GetBinLowEdge(i)), firstBin = i;
        edges.push_back(eta->GetBinUpEdge(i));
      }
    if (edges.size() < 2) return map;
    const auto phi = edgesOf(map->GetYaxis());
    auto *inside = new TH2D(fresh().c_str(), "", edges.size() - 1, edges.data(), phi.size() - 1, phi.data());
    inside->SetDirectory(nullptr);
    for (int i = 1; i <= inside->GetNbinsX(); ++i)
      for (int j = 1; j <= inside->GetNbinsY(); ++j)
        inside->SetBinContent(i, j, map->GetBinContent(firstBin + i - 1, j));
    return photonjet_science::bindMap(inside,"eta_phi_map",photonjet_science::projectionOf(map,
      photonjet_science::window(0,edges.front(),edges.back())).then("whole_eta_cells",{{"edges",edges}}),false);
  };
  for (const auto &m : std::vector<Look>{{"data", 0, "Data"}, {"phosim", 0, "Pythia #gamma Overlay"}}) {
    const bool measured = m.family == "data";
    const std::string tag = measured ? "data" : "photon", population = m.text + ", " + system;
    const std::string clusters = "event_obs/cluster_eta_phi_" + tag + ".pdf";
    const std::string jets = "event_obs/jet_eta_phi_" + tag + ".pdf";
    /* the planes drawn are bin 2 of the third axis: the preselected clusters, and the jets of the events with a
       tight, isolated leading photon */
    if (TH1 *h = allCandidates("AuAu", m.family, "clu_etaphi"))
      etaPhiMap(merged(plane(h, 2), 0.05, 2 * M_PI / 64), measured ? "Clusters" : "Weighted clusters",
                {population, Form("E_{T}^{cluster} > %g GeV, preselected", measured ? std::max(mapEtMin, treeEtMin) : mapEtMin)},
                clusters, true);
    else missing(clusters, "no Au+Au " + m.family + " QA");
    if (TH1 *h = allEvents("AuAu", m.family, "jet_etaphi"))
      etaPhiMap(merged(central(plane(h, 2), jetEtaMax), 0.2, 2 * M_PI / 8), measured ? "Jets" : "Weighted jets",
                {population, Form("p_{T}^{jet} #geq %g GeV, %s < %g", jetPtMin, absolute("#eta^{jet}").c_str(), jetEtaMax),
                 "Events with a tight, isolated photon"}, jets, false);
    else missing(jets, "no Au+Au " + m.family + " QA");
  }
  /* jets of the events with a tight, isolated photon, from the jet threshold of the maps in 2 GeV bins */
  if (TH1 *data = allEvents("AuAu", "data", "jet_pt")) {
    const Cut withPhoton = {1, 1, 2};  /* of the events axis of jet_pt: those with a tight, isolated leading photon */
    auto shown = [&](const Bins &b) { return rebinToWidth(within(b, jetPtMin, b.edges.back()), 2.0); };
    const Bins all = project(data, 0, {withPhoton}), d = shown(all);
    std::vector<Curve> curves;
    for (const auto &s : sims)
      if (TH1 *h = allEvents("AuAu", s.family, "jet_pt")) {
        const Bins b = shown(project(h, 0, {withPhoton}));
        curves.push_back(lineCurve(scaledTo(b, d.sum()), s.color, s.text));
      }
    curves.push_back(dataCurve(d, "Data", true));
    const std::string jetText = Form("anti-k_{T} R = %g jets", selectionValue("AuAu", "jet", "radius", 0.4));
    eventDistribution(curves, all.edges.front(), all.edges.back(), true, "Jet p_{T} [GeV]", "Jets", {},
                      {span, jetText, "Events with a tight, isolated photon"}, "event_obs/jet_pt.pdf");
  } else missing("event_obs/jet_pt.pdf", "no Au+Au DATA QA");
}

/* ---- simulation stitching -------------------------------------------------------------------------------------- */

/* Text with its absolute-value bars as plain strokes of the text height: the sized delimiters are drawn short and
   thin, with a gap, beside a symbol that carries an index. */
std::string plainBars(std::string text) {
  for (const std::string sized : {"#left|", "#right|"})
    for (size_t at = text.find(sized); at != std::string::npos; at = text.find(sized))
      text.replace(at, sized.size(), sized == "#left|" ? "|#scale[0.4]{ }" : "|");
  return text;
}

/* The fewest points a display fit (stitching) or a guide curve (turnOnCurve) is drawn through. */
const int kFewestFitPoints = 6;

/* The JSON record saved as text at `path` of the QA file; a null value when the file has none. The saved-model
   sections below read their records through it as well. */
PlotJson savedJson(const std::string &path) {
  auto *record = dynamic_cast<TObjString *>(gQA->Get(path.c_str()));
  return record ? PlotJson::parse(record->GetString().Data()) : PlotJson();
}

/* The stitched simulated spectrum of the samples of `family` (phosim: photon+jet, jetsim: inclusive jet) of the
   package `pkg`: at generator level (`generator`: against the truth photon or truth-jet p_T the samples are stitched
   on) or of the reconstructed leading tight, isolated photons. Read from qa/<pkg>/stitching/<family>_generator or
   <family>_reco: `total`, one `sample_<name>` per sample and the record `recipe`. Drawn in pb per GeV.
   Every sample where the QA file has it (open rings, one colour each, in the IAN's colour order), the stitched sum
   (filled black points inside the rings) with a descriptive fit, and sum / fit below. One layout for p+p and Au+Au,
   generator truth and reconstructed tags, the IAN's: a tall canvas, 1 GeV bins over the IAN's range, the header at
   the top left of a y range opened for it, the fit and the legend in one column under it at the right, the ratio from
   0.5 to 1.5. The fit is the IAN's display fit: a smooth curve through the stitched points shown, by least squares on
   their logarithms with equal weights (the modified power law for a generator spectrum, the log-parabola in
   E_T / 20 GeV for the reconstructed tags); the RMS is that of sum / fit - 1 over the same points.
   Computed here for display, not saved by the qa step: that fit and its RMS, the sum of the single- and the
   double-interaction part of a p+p sample and, where the sample records carry no ownership window (p+p), the
   boundaries between the samples.
   An incomplete set of samples is noted in gray under the header; the coverage of each figure is in its receipt. */
void stitching(const std::string &pkg, const std::string &family, bool generator) {
  const bool pp = pkg == "PP", photons = family == "phosim";
  const std::string name = std::string("simulation_stitching/") + (pp ? "pp" : "auau") + "_" +
                           (photons ? "photonjet" : "inclusivejet") + "_" + (generator ? "generator" : "reco_tag") +
                           "_stitching.png";
  TDirectory *dir =
    gQA->GetDirectory(("qa/" + pkg + "/stitching/" + family + (generator ? "_generator" : "_reco")).c_str());
  auto *total = dir ? dynamic_cast<TH1 *>(dir->Get("total")) : nullptr;
  if (!total) return missing(name, "no " + pkg + " " + family + " simulation QA");
  PlotJson recipe = PlotJson::object();
  if (auto *saved = dynamic_cast<TObjString *>(dir->Get("recipe")))
    recipe = PlotJson::parse(saved->GetString().Data());
  const bool reference = recipe.value("mode", "nominal") == "reference";
  gProjectionSteps.push_back({{"operation","stitching_recipe"},{"recipe",recipe}});
  /* Shown over the IAN's ranges: a generator spectrum of the photon+jet samples to 40 GeV and of the inclusive-jet
     samples to 50 GeV, p+p from 5 and 9 GeV, Au+Au from 15 GeV; the reconstructed tags over the QA step's range, or
     from 15 to 35 GeV when the recipe gives none. */
  double ptLow = 15, ptHigh = 35;
  if (generator) ptLow = pp ? (photons ? 5 : 9) : 15, ptHigh = photons ? 40 : 50;
  else if (recipe.contains("range_gev")) {
    ptLow = recipe.at("range_gev")[0].get<double>();
    ptHigh = recipe.at("range_gev")[1].get<double>();
  }
  /* Drawn in 1 GeV bins, per GeV: a finer saved spectrum (the packages' 0.5 GeV generator bins) is summed in whole
     groups starting at a GeV edge. */
  int summed = 1;
  double binWidth = 0;
  const auto display = [&](const TH1 *h) {
    Bins b = binsOf(h);
    const double width = b.n() ? b.width(0) : 1.0;
    const int k = static_cast<int>(std::lround(1.0 / width));
    bool whole = k > 1 && std::fabs(k * width - 1.0) < 1e-9 && std::fabs(b.edges[0] - std::round(b.edges[0])) < 1e-9;
    for (int i = 1; whole && i < b.n(); ++i) whole = std::fabs(b.width(i) - width) < 1e-9;
    if (whole) b = rebin(b, k), summed = k;
    binWidth = b.n() ? b.width(0) : 0;
    return within(perWidth(b), ptLow, ptHigh);
  };
  /* a sample: its name, first as saved and then as the legend shows it, and its display bins */
  struct Sample { std::string label; Bins bins; };
  std::vector<Sample> samples;
  const PlotJson coverage = savedJson("input_coverage");
  bool complete = true;
  PlotJson plottedSamples = PlotJson::array();
  /* the ownership windows of a family, from the sample records of the saved metadata where they carry them: the lower
     edges above the first are the boundaries between its samples */
  std::set<double> windowLows;
  bool filterPhoton = false;
  TIter next(dir->GetListOfKeys());
  while (auto *k = static_cast<TKey *>(next()))
    if (std::string(k->GetName()).rfind("sample_", 0) == 0) {
      const std::string sample = std::string(k->GetName()).substr(7), key = pkg + "/" + sample;
      samples.push_back({sample, display(static_cast<TH1 *>(dir->Get(k->GetName())))});
      complete = complete && coverage.contains("samples") && coverage.at("samples").contains(key)
        && coverage.at("samples").at(key).value("status", "UNKNOWN") == "COMPLETE";
      plottedSamples.push_back(key);
      if (gQAMetadata.is_null() || !gQAMetadata.at("packages").contains(pkg)) continue;
      const PlotJson &package = gQAMetadata.at("packages").at(pkg);
      if (!package.contains("samples") || !package.at("samples").contains(sample)) continue;
      const PlotJson normalization = package.at("samples").at(sample).value("normalization", PlotJson::object());
      if (!normalization.contains("ownership")) continue;
      const PlotJson windowLow = normalization.at("ownership").at("window_gev").at(0);
      if (windowLow.is_number()) windowLows.insert(windowLow.get<double>());
      filterPhoton = filterPhoton || normalization.at("ownership").value("owner", "") == "generator_filter_photon_pt";
    }
  /* An incomplete population is recorded for the receipt (as inputNote does for every other figure). */
  const std::string populationNote =
    complete ? std::string() : forReceipt(std::string("Partial ") + (pp ? "p+p" : "Au+Au") + " simulation");
  gProjectionSteps.push_back({{"operation","figure_input_coverage"},{"samples",plottedSamples},
                              {"status",complete ? "COMPLETE" : "PARTIAL"}});
  /* A p+p sample is saved in two parts, single-interaction and double-interaction (suffix _double), each with its
     physics weight: the two are added bin by bin without rescaling and drawn as one sample, named as in the IAN
     (run28_photonjet12 as Photon12, run28_jet20 as Jet20). */
  if (pp) {
    std::map<std::string, Bins> merged;
    for (const auto &sample : samples) {
      std::string key = sample.label;
      const auto suffix = key.rfind("_double");
      if (suffix != std::string::npos && suffix + 7 == key.size()) key.resize(suffix);
      if (!merged.count(key)) merged[key] = sample.bins;
      else {
        Bins &slice = merged.at(key);
        if (slice.edges != sample.bins.edges) throw std::runtime_error("stitching slice binning differs");
        for (int b = 0; b < slice.n(); ++b) slice.v[b] += sample.bins.v[b], slice.e2[b] += sample.bins.e2[b];
      }
    }
    samples.clear();
    for (const auto &[key, bins] : merged) {
      const std::string prefix = photons ? "run28_photonjet" : "run28_jet";
      const std::string label = key.rfind(prefix, 0) == 0 ? (photons ? "Photon" : "Jet") + key.substr(prefix.size())
                                                         : key;
      samples.push_back({label, bins});
    }
  } else  /* photon12, jet20: named as the p+p samples are (Photon12, Jet20) */
    for (auto &sample : samples)
      if (!sample.label.empty() && sample.label[0] >= 'a' && sample.label[0] <= 'z') sample.label[0] += 'A' - 'a';
  /* in the order of the first number in each label (Photon12 before Photon20); by name where a label has none */
  std::sort(samples.begin(), samples.end(), [](const auto &a, const auto &b) {
    const auto ai = a.label.find_first_of("0123456789"), bi = b.label.find_first_of("0123456789");
    return ai != std::string::npos && bi != std::string::npos
             ? std::atoi(a.label.c_str() + ai) < std::atoi(b.label.c_str() + bi)
             : a.label < b.label;
  });
  const Bins all = display(total);
  if (summed > 1)
    gProjectionSteps.push_back({{"operation","display_rebin"},{"bins_summed",summed},{"bin_width_gev",binWidth}});
  /* An empty stitched sum is listed, not drawn as an empty canvas. */
  int first = -1, last = -1;
  for (int i = 0; i < all.n(); ++i)
    if (all.v[i] > 0) { if (first < 0) first = i; last = i; }
  if (first < 0) return missing(name, "no " + pkg + " " + family + " entries in the displayed range");
  const double x0 = all.edges[first], x1 = ptHigh;  /* from the first filled bin: no empty strip at the left */
  /* The display fit: linear least squares on the logarithms of the stitched points shown, each with the same weight.
     In the coefficients of the fit, the generator form exp[a - (b + c ln x + d x) ln x] reads
     exp[p0 + p1 ln x + p2 (ln x)^2 + p3 x ln x]. */
  const bool parabola = !generator;
  TGraph logarithms;
  for (int i = first; i <= last; ++i)
    if (all.v[i] > 0) logarithms.SetPoint(logarithms.GetN(), all.center(i), std::log(all.v[i]));
  TF1 *model = nullptr;
  double rms = std::nan("");
  if (logarithms.GetN() >= kFewestFitPoints) {
    const char *terms = parabola ? "1 ++ log(x/20) ++ log(x/20)*log(x/20)" : "1 ++ log(x) ++ log(x)*log(x) ++ x*log(x)";
    TF1 linear(fresh().c_str(), terms, x0, all.edges[last + 1]);
    const TFitResultPtr fitResult=logarithms.Fit(&linear, "QNS");
    model = new TF1(fresh().c_str(), parabola ? "exp([0]+[1]*log(x/20)+[2]*pow(log(x/20),2))"
                                              : "exp([0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*x*log(x))", x0, x1);
    for (int p = 0; p < linear.GetNpar(); ++p) model->SetParameter(p, linear.GetParameter(p));
    model->SetLineColor(cBlack), model->SetLineStyle(1), model->SetLineWidth(2), model->SetNpx(300);
    double squares = 0;
    for (int i = 0; i < logarithms.GetN(); ++i)
      squares += std::pow(std::exp(logarithms.GetY()[i]) / model->Eval(logarithms.GetX()[i]) - 1, 2);
    rms = std::sqrt(squares / logarithms.GetN());
    photonjet_science::bindFit(model,"stitched_generator_fit",all.physics,fitResult,
        {{"options","QN"},{"fitted_formula",terms},{"input_transform","log(y)"},{"weights","equal"},
         {"range",{x0,all.edges[last+1]}},{"relative_rms",rms}});
    gProjectionSteps.push_back({{"operation","display_fit"},
                                {"function",parabola ? "log_parabola_et_over_20" : "modified_power_law"},
                                {"range_gev",{x0, all.edges[last + 1]}},{"points",logarithms.GetN()},
                                {"weights","equal, on the logarithm of the stitched sum"},{"rms",rms}});
  }
  /* Where the sample records carry no window (p+p), the boundaries are read off the spectra: the edges at which the
     largest share of the stitched sum passes to the next sample. */
  std::vector<double> boundaries;
  if (windowLows.size() > 1) boundaries.assign(std::next(windowLows.begin()), windowLows.end());
  else if (generator) {
    int owner = -1;
    bool ordered = true;
    for (int i = first; i <= last; ++i) {
      int largest = -1;
      for (size_t s = 0; s < samples.size(); ++s) {
        const double share = samples[s].bins.v[i];
        if (share > 0 && (largest < 0 || share > samples[largest].bins.v[i])) largest = s;
      }
      if (largest < 0 || largest == owner) continue;
      if (owner >= 0) boundaries.push_back(all.edges[i]), ordered = ordered && largest > owner;
      owner = largest;
    }
    if (!ordered) boundaries.clear();
  }
  /* The IAN's proportions and text sizes on its tall canvas: the upper frame two thirds of the height, the ratio
     frame a fifth, the x title directly under its labels. The small text has the size of the shared header lines
     (kLinePx, a half pixel on purpose: the exponent and the indices of the fit note are then set close to their
     symbols), and the axis ticks are short enough to end clear of the header. */
  const double markPx = 35, linePx = 27, smallPx = kLinePx, pitch = 39, labelPx = 27, titlePx = 32, tickPx = 13;
  TCanvas cv("cv", "", 900, 950);
  auto [upper, lower] = ratioPads(cv);
  upper->SetPad(0, 0.32, 1, 1), lower->SetPad(0, 0, 1, 0.32);
  upper->SetMargin(0.175, 0.045, 0.02, 0.035), lower->SetMargin(0.175, 0.045, 0.30, 0.03);
  upper->cd();
  gPad->SetLogy();
  /* On the log scale a bin without entries has no point. A sample point whose statistical error is about as large as
     its value (above 90% of it: its bar would run more than a decade down) carries no information and is not drawn
     either. */
  const auto filledBins = [](Bins b) {
    for (int i = 0; i < b.n(); ++i) if (!(b.v[i] > 0)) b.v[i] = std::nan("");
    return b;
  };
  const auto informative = [](Bins b) {
    for (int i = 0; i < b.n(); ++i) if (!(b.v[i] > 0) || !(b.e2[i] < 0.81 * b.v[i] * b.v[i])) b.v[i] = std::nan("");
    return b;
  };
  const Bins stitched = filledBins(all);
  for (auto &sample : samples) sample.bins = informative(sample.bins);
  /* The frame starts at the decade below the stitched sum for a generator spectrum, where every sample follows the
     sum inside its window, and at the decade below the lowest sample point for the reconstructed tags, whose samples
     overlap over decades. */
  double least = std::numeric_limits<double>::infinity();
  for (double v : stitched.v) if (v > 0) least = std::min(least, v);
  if (!generator)
    for (const auto &sample : samples)
      for (double v : sample.bins.v) if (v > 0) least = std::min(least, v);
  const double low = std::pow(10, std::floor(std::log10(least / 2)));
  std::vector<Extent> drawn;
  extents(drawn, stitched, true, x0, x1);
  for (const auto &sample : samples) extents(drawn, sample.bins, true, x0, x1);
  /* The most central class of the package. The header of an embedded Au+Au generator spectrum names the events it
     holds, as the IAN's does: the population the saved recipe states (the centrality range and the vertex window of
     the events counted), or for a file without that record the span of all classes and the vertex window of the
     analysis. The header of the reconstructed tags gives the range of the most central class alone. */
  const std::string mostCentral = classesOf(pkg).empty() ? "" : classesOf(pkg).front();
  std::string level = generator ? "Generator level, no reconstruction" : "Leading tight + isolated photon tag";
  if (generator && !pp && recipe.contains("population"))
    level = "Generator truth; " + rangeDashes(rewritten(recipe.at("population").get<std::string>(), "z_vtx", "z_{vtx}"));
  else if (generator && !pp && !mostCentral.empty())
    level = "Generator truth; " + centralityRange(lowOf(mostCentral), highOf(classesOf(pkg).back()))
            + Form(", |z_{vtx}| < %g cm", selectionValue(pkg,"event","vertex_abs_z_max_cm",10));
  std::vector<std::string> headerLines = {
    (pp ? "p+p PYTHIA 8,#scale[0.5]{ }" : "Au+Au embedded PYTHIA 8,#scale[0.9]{ }") + energyOf(pkg),
    photons ? "Photon+jet samples" : "Inclusive-jet samples", level};
  if (!generator) {
    const double vertex = reference ? recipe.at("vertex_abs_z_max_cm").get<double>()
                                    : selectionValue(pkg,"event","vertex_abs_z_max_cm",pp ? 60 : 10);
    headerLines.push_back(plainBars(photonEtaLabel(pkg)) + Form(", |z_{vtx}| < %g cm", vertex)
                          + (pp ? std::string() : ", " + centralityRange(lowOf(mostCentral), highOf(mostCentral))));
  }
  if (generator && !boundaries.empty()) {  /* as the IAN states them: the boundaries between the samples */
    std::string list;
    for (double edge : boundaries) list += (list.empty() ? "" : ", ") + std::string(Form("%g", edge));
    const char *ownership = boundaries.size() > 1 ? "Ownership boundaries " : "Ownership boundary ";
    headerLines.push_back((filterPhoton ? "Generator-filter boundary " : ownership) + list + " GeV");
  }
  /* the weights of the spectrum, in the IAN's smaller line under the header */
  const std::string interactions = reference ? "single-interaction samples" : "single- and double-interaction samples";
  std::string weights = (generator ? "Cross-section weighted; " : "Physics weighted; ") + interactions;
  if (!pp && generator) weights = "Cross-section weighted; no centrality reweighting";
  if (!pp && !generator) weights = "Cross-section #times vertex #times centrality weights";
  const std::vector<std::string> smallLines = nonEmpty({weights, populationNote});
  const std::string formula = parabola ? "exp[a + b ln x + c (ln x)^{2}]"
                                       : "exp[a #minus (b + c ln x + d x) ln x]";
  const std::string variable = std::string("x = ") + (generator ? "p_{T}" : "E_{T}")
                               + (parabola ? "/(20 GeV)" : "/(1 GeV)");
  /* The header at the top left. At the right, from the header's lower edge down, as in the IAN: the fit, then the
     legend in one column. The frame is opened until the points clear all three. */
  const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin(), T = gPad->GetTopMargin(), B = gPad->GetBottomMargin();
  const double w = padWidthPx(), h = padHeightPx(), plot = 1 - T - B;
  const double hx = L + 30 / w, markY = 1 - T - (14 + markPx) / h, legendPx = 24, row = 25 / h, smallStep = 31;
  const double headBottom = markY - (headerLines.size() * pitch + smallLines.size() * smallStep) / h;
  double headWidth = textWidth(sPhenixInternalLabel(), markPx), legendText = textWidth("Stitched sum", legendPx);
  for (const auto &text : headerLines) headWidth = std::max(headWidth, textWidth(text, linePx));
  for (const auto &text : smallLines) headWidth = std::max(headWidth, textWidth(text, smallPx));
  for (const auto &sample : samples) legendText = std::max(legendText, textWidth(sample.label, legendPx));
  const double symbol = 44 / w, columnWidth = symbol + legendText + 14 / w, lx2 = 1 - R - 30 / w;
  const double fitWidth = symbol + std::max(textWidth(formula, smallPx), textWidth(variable, smallPx)) + 14 / w;
  const double fitTop = headBottom - 4 / h, pointsTop = fitTop - (model ? 56 / h : 0);
  const double legendBottom = pointsTop - (samples.size() + 1) * row;
  const auto user = [&](double ndc) { return x0 + (ndc - L) / (1 - L - R) * (x1 - x0); };
  /* the frame top that keeps the points under a text block from `left` (NDC) to the right edge, down to `bottom` */
  const auto topUnder = [&](double left, double bottom) {
    return topBelow(drawn, low, true, user(left - 12 / w), x1, (bottom - 10 / h - B) / plot);
  };
  double high = topBelow(drawn, low, true, x0, x1, 0.97);
  const double headRight = user(hx + headWidth + 12 / w);
  high = std::max(high, topBelow(drawn, low, true, x0, headRight, (headBottom - 40 / h - B) / plot));
  if (model) high = std::max(high, topUnder(lx2 - fitWidth, pointsTop));
  high = std::max(high, topUnder(lx2 - columnWidth, legendBottom));
  if (!(high > low)) high = 10 * low;
  /* Every stitched spectrum is in pb: the QA step scales each sample's weighted events by its pb per generated
     event. The y title follows the weights: the cross section of the reference recipe, the analysis weights else. */
  TH1F *f = photonjet_science::drawFrame(*gPad,x0, low, x1, high);
  /* a cross section where the samples are drawn from their raw generator counts (p+p: the reference mode, and the
     nominal mode with its single- and double-interaction mix) */
  const bool crossSection = reference || (pp && recipe.contains("interactions"));
  const char *spectrum = crossSection ? "d#sigma/dp_{T} [pb/GeV]" : "Weighted generator spectrum [pb/GeV]";
  f->GetYaxis()->SetTitle(generator ? spectrum : "Weighted photon tags [pb/GeV]");
  pixelAxes(f, labelPx, titlePx);
  f->GetXaxis()->SetLabelSize(0);
  f->GetXaxis()->SetTitleSize(0);
  tickLengths(f, tickPx);
  double y = markAndLines(hx, markY, headerLines, markPx, linePx, pitch);
  TLatex note;
  note.SetNDC();
  note.SetTextFont(43);
  note.SetTextSize(smallPx);
  for (const auto &text : smallLines) {
    note.SetTextColor(text == populationNote ? gGray : cBlack);
    note.DrawLatex(hx, y -= smallStep / h, text.c_str());
  }
  /* the IAN's colour order: vermilion (cOrange), blue, green, red, purple, then three further colours for a family
     of more than five samples; the highest photon sample purple. The rings are open circles with the heaviest
     outline the marker set has, so they stay visible around the black points. */
  const int colors[] = {cOrange, cBlue, cGreen, cRed, cPurple, kOrange + 3, kCyan + 2, kSpring + 5};
  std::vector<TGraphErrors *> rings;
  for (size_t i = 0; i < samples.size(); ++i) {
    const int color = photons && i > 0 && i + 1 == samples.size() ? cPurple : colors[i % 8];
    rings.push_back(markers(samples[i].bins, color, 107, true, false, 0, 1.7));
    photonjet_science::draw(rings.back(),"P");
  }
  TGraphErrors *stitchedPoints = markers(stitched, cBlack, 20, true, false, 0, 1.15);
  photonjet_science::draw(stitchedPoints,"P");
  if (model) {
    photonjet_science::draw(model,"SAME");
    TLegend *fitted = ianLegend(lx2 - fitWidth, pointsTop, lx2, fitTop, symbol / fitWidth);
    fitted->SetTextSize(smallPx);
    fitted->AddEntry(model, ("#splitline{" + formula + "}{" + variable + "}").c_str(), "l");
    photonjet_science::draw(fitted);
  }
  TLegend *leg = ianLegend(lx2 - columnWidth, legendBottom, lx2, pointsTop, symbol / columnWidth);
  leg->SetTextSize(legendPx);
  leg->AddEntry(stitchedPoints, "Stitched sum", "p");
  for (size_t i = 0; i < samples.size(); ++i) leg->AddEntry(rings[i], samples[i].label.c_str(), "p");
  photonjet_science::draw(leg);
  lower->cd();
  const double ratioLow = 0.5, ratioHigh = 1.5;
  /* The photon of a generator spectrum. Au+Au photon+jet samples are stitched on the generator-filter photon. Their
     spectrum is of that photon unless the recipe names the maximum signal truth photon as its coordinate: then every
     sample shows its rise below the boundary, as in the IAN. */
  const std::string pt = apart("p_{T}") + " [GeV]";
  const bool truthPhoton = recipe.value("coordinate", "").find("maximum signal truth-photon") != std::string::npos;
  std::string xtitle = "Reconstructed photon " + apart("E_{T}") + " [GeV]";
  if (generator && photons)
    xtitle = (filterPhoton && !truthPhoton ? "Generator-filter photon " : "Maximum truth photon ") + pt;
  else if (generator) xtitle = "Maximum anti-" + apart("k_{T}") + " R = 0.4 truth-jet " + pt;
  TH1F *rf = ratioFrame(x0, x1, ratioLow, ratioHigh, xtitle.c_str(), "Spectrum / fit");
  pixelAxes(rf, labelPx, titlePx);
  rf->GetYaxis()->SetNdivisions(505);
  rf->GetXaxis()->SetTitleOffset(1.15);
  tickLengths(rf, tickPx);
  /* about five labelled divisions, as in the IAN: every 5 GeV over a range of up to 25 GeV, every 10 GeV beyond */
  f->GetXaxis()->SetNdivisions(505), rf->GetXaxis()->SetNdivisions(505);
  frameReferenceLine(1);
  if (model) {
    Bins r = all;
    for (int i = 0; i < r.n(); ++i) {
      const double m = model->Eval(r.center(i));
      r.v[i] = all.v[i] > 0 ? all.v[i] / m : std::nan(""), r.e2[i] = all.e2[i] / (m * m);
    }
    r.physics=photonjet_science::Trace::combine("spectrum_over_fitted_function",{all.physics,
      photonjet_science::Trace::fromJson(photonjet_science::tag(model,"physics_series").at("trace"))},{{"fit_parameters",photonjet_science::nativeValues(model)}});
    r.seriesRole="stitching/spectrum_over_fit";
    photonjet_science::annotation(*gPad,"stitching/relative_rms",photonjet_science::scalar(rms,"fraction","root mean square of spectrum/fit minus one"),r.physics);
    photonjet_science::draw(markers(r, cBlack, 20, true, false, 0, 1.15),"P");
    /* The bin width at the top left and the RMS at the top right, as in the IAN, both under the ticks of the frame;
       each at the bottom instead when a ratio point or its error bar reaches it at the top and the bottom is free. */
    const double Lr = gPad->GetLeftMargin(), Rr = gPad->GetRightMargin(), Tr = gPad->GetTopMargin();
    const double Br = gPad->GetBottomMargin(), wr = padWidthPx(), hr = padHeightPx();
    const double inset = tickPx + 4;
    const double depth = (inset + smallPx + 6) / hr / (1 - Tr - Br) * (ratioHigh - ratioLow);
    const auto reached = [&](double from, double to, bool top) {
      for (int i = 0; i < r.n(); ++i) {
        if (!std::isfinite(r.v[i])) continue;
        const double x = Lr + (r.center(i) - x0) / (x1 - x0) * (1 - Lr - Rr);
        const double e = std::isfinite(r.e2[i]) ? std::sqrt(r.e2[i]) : 0.0;
        if (x < from - 8 / wr || x > to + 8 / wr) continue;
        if (top ? r.v[i] + e > ratioHigh - depth : r.v[i] - e < ratioLow + depth) return true;
      }
      return false;
    };
    TLatex t;
    t.SetNDC();
    t.SetTextFont(43);
    t.SetTextSize(smallPx);
    const auto corner = [&](const std::string &text, bool left) {
      const double width = textWidth(text, smallPx), from = left ? Lr + 30 / wr : 1 - Rr - 40 / wr - width;
      const bool top = !reached(from, from + width, true) || reached(from, from + width, false);
      t.SetTextAlign(top ? 13 : 11);
      t.DrawLatex(from, top ? 1 - Tr - inset / hr : Br + inset / hr, text.c_str());
    };
    if (binWidth > 0) corner(Form("%g GeV bins; statistical errors", binWidth), true);
    if (std::isfinite(rms)) corner(Form("RMS %.1f%%", 100 * rms), false);
  }
  savePhysicsQA(cv, name);
}

/* ---- trigger study --------------------------------------------------------------------------------------------- */

/* The curve drawn through a photon-trigger turn-on, as in the IAN: a generalized sigmoid that rises from a floor to
   the trigger's plateau, by least squares through the measured points with equal weights. It guides the eye; the
   plateau and the turn-on point quoted on the canvas are computed from the points. Null with too few points. */
TF1 *turnOnCurve(const Bins &eff, double plateau, int color) {
  TGraph measured;
  for (int b = 0; b < eff.n(); ++b)
    if (std::isfinite(eff.v[b])) measured.SetPoint(measured.GetN(), eff.center(b), eff.v[b]);
  if (measured.GetN() < kFewestFitPoints || !std::isfinite(plateau)) return nullptr;
  const double first = measured.GetX()[0], last = measured.GetX()[measured.GetN() - 1];
  double middle = first;
  for (int i = 0; i < measured.GetN() && measured.GetY()[i] < 0.5 * plateau; ++i) middle = measured.GetX()[i];
  /* the parameters: [0] the floor, [1] the plateau, fixed to the measured one, [2] the midpoint, which starts at the
     last point before the efficiency first reaches half the plateau, [3] the width in GeV, [4] the exponent */
  auto *f = new TF1(fresh().c_str(), "[0] + ([1] - [0]) / pow(1 + exp(-(x - [2]) / [3]), [4])", first, last);
  f->SetParameters(0, plateau, middle, 0.8, 1);
  f->FixParameter(1, plateau);
  f->SetParLimits(0, 0, 0.3);
  f->SetParLimits(2, first - 5, last);
  f->SetParLimits(3, 0.05, 10);
  f->SetParLimits(4, 0.05, 50);
  const TFitResultPtr fitResult=measured.Fit(f, "QNS");
  f->SetLineColor(color);
  f->SetLineWidth(3);
  f->SetNpx(400);
  photonjet_science::bindFit(f,"trigger_turn_on_guide",eff.physics,fitResult,
      {{"options","QN"},{"weights","equal"},{"fixed_plateau",plateau},{"initial",{0.,plateau,middle,.8,1.}},
       {"limits",{{"0",{0.,.3}},{"2",{first-5,last}},{"3",{.05,10.}},{"4",{.05,50.}}}}});
  return f;
}

/* Two figures per package from the trigger study of the qa step, over the runs in which every trigger of the study is
   enabled (the qa step's selection). Read from qa/<pkg>/triggers: the record `study` ("key=value;" with runs, first,
   reference and correction), `reference`, `reference_denominator` and, per photon trigger, `spectrum_<trigger>` and
   `numerator_<trigger>`.
   First figure: the maximum-cluster-energy spectra of the reference trigger and of the photon triggers. Second
   figure: the photon-trigger turn-ons against the reference. Each is one frame in the arrangement of its IAN
   original: the wide Au+Au canvas with the legend and the header over the falling spectra and the turn-on legend at
   the lower right; the p+p canvas with the run line and the legend at the top and the sPHENIX mark in the free corner
   below.
   Computed here from the saved numerators and the reference denominator, not saved by the qa step: the efficiencies
   with their binomial errors, each plateau (the error-weighted mean above the lowest photon E_T of the analysis), the
   energy at which the efficiency reaches 90% of it, and the guide curves (turnOnCurve). */
void triggerStudy(const std::string &pkg) {
  const bool pp = pkg == "PP";
  const std::string spectra = pp ? "event_obs/pp_max_cluster_energy_overlay.png"
                                 : "event_obs/hMaxClusterEnergy_groupOverlay.png";
  const std::string turnOn = pp ? "event_obs/pp_max_cluster_energy_turnon.png"
                                : "event_obs/hMaxClusterEnergy_groupTurnOnOverlay.png";
  TDirectory *dir = gQA->GetDirectory(("qa/" + pkg + "/triggers").c_str());
  auto *study = dir ? dynamic_cast<TObjString *>(dir->Get("study")) : nullptr;
  std::map<std::string, std::string> info;
  if (study) {
    std::istringstream in(study->GetString().Data());
    for (std::string item; std::getline(in, item, ';');)
      if (item.find('=') != std::string::npos) info[item.substr(0, item.find('='))] = item.substr(item.find('=') + 1);
  }
  const int runs = std::atoi(info["runs"].c_str());
  if (runs <= 0) {
    missing(spectra, "no run with every trigger of the study enabled");
    return missing(turnOn, "no run with every trigger of the study enabled");
  }
  std::vector<std::string> photons;
  TIter next(dir->GetListOfKeys());
  while (auto *k = static_cast<TKey *>(next()))
    if (std::string(k->GetName()).rfind("spectrum_", 0) == 0) photons.push_back(std::string(k->GetName()).substr(9));
  std::sort(photons.begin(), photons.end(), [](const std::string &a, const std::string &b) {
    return std::atof(a.c_str() + a.find_first_of("0123456789")) < std::atof(b.c_str() + b.find_first_of("0123456789"));
  });
  /* corrected: the study is prescale-corrected run by run (correction=run_live_over_scaled) and the first figure
     draws the corrected spectra (spectrum_<trigger>); else it is the live-bit study, whose first figure draws the
     reference-trigger events with the photon trigger live (numerator_<trigger>). */
  const bool corrected = info["correction"] == "run_live_over_scaled";
  const std::string system = (pp ? "p+p, " : "Au+Au, ") + energyOf(pkg);
  /* the p+p study starts with the photon-trigger firmware update at run 47289, and its run line says so as the IAN's */
  const int firmwareUpdateRun = 47289;
  std::string runLine = Form("Photon-trigger study, %d runs", runs);
  if (pp) runLine = std::atoi(info["first"].c_str()) >= firmwareUpdateRun
                      ? Form("%d runs (After firmware update at run %d)", runs, firmwareUpdateRun)
                      : Form("%d runs (from run %s)", runs, info["first"].c_str());
  const std::vector<std::string> lines = nonEmpty({system, corrected ? "Prescale-corrected per run" : "", runLine});
  /* the IAN's colours: Au+Au blue and red, p+p blue, green and magenta, from the lowest photon threshold up */
  const std::vector<int> colors = pp ? std::vector<int>{cBlue, cGreen, kMagenta + 1, cRed}
                                     : std::vector<int>{cBlue, cRed, cGreen, cPurple};
  const auto saved = [&](const std::string &name) { return binsOf(static_cast<TH1 *>(dir->Get(name.c_str()))); };
  std::vector<Curve> curves = {lineCurve(saved("reference"), cBlack, plainBars(latexCuts(info["reference"])))};
  /* the live-bit study draws what its efficiency divides: the reference-trigger events and those with the photon
     trigger live */
  for (size_t i = 0; i < photons.size(); ++i)
    curves.push_back(lineCurve(saved((corrected ? "spectrum_" : "numerator_") + photons[i]), colors[i % 4],
                               corrected ? photons[i] : photons[i] + " + MBD"));
  const Bins &reference = curves.front().bins;
  /* The study requires a cluster of at least 1 GeV: the spectra are drawn from the first filled bin. The Au+Au axis
     starts there and the p+p axis at zero, as in the IAN. */
  int firstFilled = 0;
  while (firstFilled + 1 < reference.n() && !(reference.v[firstFilled] > 0)) ++firstFilled;
  const double energyLow = pp ? reference.edges.front() : reference.edges[firstFilled];
  const double energyHigh = reference.edges.back();
  const std::string energyTitle = pp ? "Maximum cluster energy [GeV]"
                                     : "Max cluster energy [GeV], " + apart("E_{clus}") + " #geq 1 GeV";
  /* The canvases are 800 high and 1067 (p+p, 4:3) or 1200 (Au+Au, 3:2) wide.
     Text and axis sizes in proportion to the canvas width: the two systems are printed side by side at one width.
     The symbol of the Au+Au title is a fragment of its own (apart), so that its index sits on it as in the IAN, and
     the axis ticks are short enough to end clear of the header. */
  const double width = pp ? 1067 : 1200, scale = width / 1200;
  const double markPx = 37 * scale, linePx = 30 * scale, pitch = 41 * scale, labelPx = 34 * scale;
  const double titlePx = 37.5 * scale, tickPx = 16 * scale;
  /* the frame of either figure on the current canvas, with its titles, text sizes and ticks */
  const auto studyFrame = [&](double low, double high, const char *ytitle) {
    TH1F *f = photonjet_science::drawFrame(*gPad,energyLow, low, energyHigh, high);
    f->GetXaxis()->SetTitle(energyTitle.c_str());
    f->GetYaxis()->SetTitle(ytitle);
    pixelAxes(f, labelPx, titlePx);
    f->GetXaxis()->SetTitleOffset(1.25);
    f->GetYaxis()->SetTitleOffset(1.35);
    tickLengths(f, tickPx);
  };
  {
    TCanvas cv("cv", "", width, 800);
    cv.SetMargin(0.13, 0.04, 0.15, 0.04);
    cv.SetLogy();
    const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin(), T = gPad->GetTopMargin(), B = gPad->GetBottomMargin();
    const double w = padWidthPx(), h = padHeightPx(), plot = 1 - T - B;
    std::vector<Extent> drawn;
    std::vector<Bins> shown;
    double labelWidth = 0, least = std::numeric_limits<double>::infinity();
    for (const auto &c : curves) {
      shown.push_back(within(c.bins, reference.edges[firstFilled], energyHigh));
      extents(drawn, shown.back(), false, energyLow, energyHigh);
      for (double v : shown.back().v) if (v > 0) least = std::min(least, v);
      labelWidth = std::max(labelWidth, textWidth(c.text, linePx));
    }
    /* At the right, over the falling spectra: Au+Au the legend with the header under it, p+p the run line with the
       legend under it, both ending at one right edge. The p+p mark and system line stand at the lower left, under
       the spectra. */
    const std::vector<std::string> under = pp ? std::vector<std::string>{system} : lines;
    const double symbol = 60 * scale / w, legendWidth = symbol + labelWidth + 12 / w, right = 1 - R - 36 / w;
    const double runY = 1 - T - (18 + linePx) / h, ly2 = pp ? runY - 14 / h : 1 - T - 16 / h;
    const double ly1 = ly2 - curves.size() * pitch / h;
    double headWidth = textWidth(sPhenixInternalLabel(), markPx);
    for (const auto &text : under) headWidth = std::max(headWidth, textWidth(text, linePx));
    const double runLeft = right - textWidth(runLine, linePx);
    const double blockLeft = right - (pp ? legendWidth : std::max(legendWidth, headWidth));
    const double markY = pp ? B + (22 + under.size() * pitch) / h : ly1 - (22 + markPx) / h;
    const double hx = pp ? L + 0.14 * (1 - L - R) : blockLeft;
    const double blockBottom = pp ? ly1 - 12 / h : markY - (under.size() * pitch + 14) / h;
    const auto user = [&](double ndc) { return energyLow + (ndc - L) / (1 - L - R) * (energyHigh - energyLow); };
    double low = std::isfinite(least) ? least / 1.6 : 1;
    const auto frameTop = [&] {
      const double underBlock = topBelow(drawn, low, true, user(blockLeft - 14 / w), energyHigh, (blockBottom - B) / plot);
      double high = std::max(topBelow(drawn, low, true, energyLow, energyHigh, 0.93), underBlock);
      const double runBottom = (runY - 14 / h - B) / plot;
      if (pp) high = std::max(high, topBelow(drawn, low, true, user(runLeft - 14 / w), energyHigh, runBottom));
      return high > low ? high : 10 * low;
    };
    double high = frameTop();
    if (pp) {  /* the frame is opened downwards until the spectra clear the mark at the lower left */
      const double needed = bottomAbove(drawn, high, true, user(hx - 10 / w), user(hx + headWidth + 10 / w),
                                        (markY + (markPx + 12) / h - B) / plot);
      if (needed < low) low = needed, high = frameTop();
    }
    studyFrame(low, high, corrected ? "Prescale-corrected events" : "Events");
    TLegend *leg = ianLegend(blockLeft, ly1, blockLeft + legendWidth, ly2, symbol / legendWidth);
    leg->SetTextSize(linePx);
    for (size_t i = 0; i < curves.size(); ++i) {
      TH1D *hist = line(shown[i], curves[i].color, 1, 3);
      photonjet_science::draw(hist,"HIST SAME ][");
      /* a spectrum that starts inside the frame rises from the axis at its first bin, as in the IAN */
      if (shown[i].n() && shown[i].edges.front() > energyLow + 1e-9 && shown[i].v.front() > low) {
        TLine rise;
        rise.SetLineColor(hist->GetLineColor());
        rise.SetLineWidth(hist->GetLineWidth());
        rise.DrawLine(shown[i].edges.front(), low, shown[i].edges.front(), shown[i].v.front());
      }
      leg->AddEntry(hist, curves[i].text.c_str(), "l");
    }
    photonjet_science::draw(leg);
    markAndLines(hx, markY, under, markPx, linePx, pitch);
    if (pp) {
      TLatex t;
      t.SetNDC();
      t.SetTextFont(43);
      t.SetTextSize(linePx);
      t.DrawLatex(runLeft, runY, runLine.c_str());
    }
    savePhysicsQA(cv, spectra);
  }
  const Bins denominator = binsOf(static_cast<TH1 *>(dir->Get("reference_denominator")));
  struct TurnOn { Bins eff; double plateau, crossing; int color; std::string trigger, text; };
  std::vector<TurnOn> turnOns;
  /* The plateau of a trigger: its error-weighted mean efficiency above the lowest photon E_T of the analysis. The
     turn-on point: where the efficiency first reaches 90% of that plateau, interpolated between bin centres (the IAN
     quotes 90% of the plateau; a prescale-corrected plateau lies a few percent below one). */
  const double plateauFrom = photonSupport(pkg).first;
  const double turnOnShare = 0.9;  /* of the plateau; the canvas texts write it as 90% */
  for (size_t i = 0; i < photons.size(); ++i) {
    Bins eff = binsOf(static_cast<TH1 *>(dir->Get(("numerator_" + photons[i]).c_str())));
    for (int b = 0; b < eff.n(); ++b) {
      const double n = denominator.v[b], e = n > 0 ? eff.v[b] / n : std::nan("");
      /* The variance of the efficiency: for prescale-corrected, that is weighted, counts the binomial variance of a
         ratio of weighted counts as TH1::Divide with option B computes it; for plain counts e(1 - e)/n. */
      const double numVariance = eff.e2[b];
      eff.v[b] = e;
      eff.e2[b] = n > 0 ? (corrected ? (e == 1 ? 0
                                               : std::fabs((1 - 2 * e) * numVariance + e * e * denominator.e2[b]) /
                                                   (n * n))
                                     : std::max(e * (1 - e), 0.0) / n)
                        : 0;
    }
    eff.physics=photonjet_science::Trace::combine("trigger_efficiency",{eff.physics,denominator.physics},
      {{"variance",corrected?"weighted_binomial_TH1_B":"binomial"},{"zero_denominator","undefined"}});
    eff.seriesRole="trigger/efficiency/"+std::to_string(i);
    eff.measurementState.assign(eff.n(),"undefined");
    for(int b=0;b<eff.n();++b)eff.measurementState[b]=std::isfinite(eff.v[b])?"measured":"undefined";
    double weights = 0, weighted = 0;
    for (int b = 0; b < eff.n(); ++b)
      if (eff.center(b) >= plateauFrom && std::isfinite(eff.v[b]) && eff.e2[b] > 0)
        weights += 1 / eff.e2[b], weighted += eff.v[b] / eff.e2[b];
    const double plateau = weights > 0 ? weighted / weights : std::nan("");
    double crossing = std::nan("");
    for (int b = 1; b < eff.n() && std::isnan(crossing) && std::isfinite(plateau); ++b)
      if (std::isfinite(eff.v[b - 1]) && std::isfinite(eff.v[b]) && eff.v[b - 1] < turnOnShare * plateau &&
          eff.v[b] >= turnOnShare * plateau)
        crossing = eff.center(b - 1) + (turnOnShare * plateau - eff.v[b - 1]) / (eff.v[b] - eff.v[b - 1]) *
                                         (eff.center(b) - eff.center(b - 1));
    /* (the narrower Au+Au corner takes the short form; its note says what the 90% is of) */
    std::string text = photons[i] + " / MBD";
    if (std::isfinite(crossing)) text += Form(pp ? ", 90%% of plateau at %.1f GeV" : ", 90%% at %.1f GeV", crossing);
    turnOns.push_back({eff, plateau, crossing, colors[i % 4], photons[i], text});
  }
  renderSelection(pkg, true);  /* a second figure of the same package: each saved figure names its inputs */
  TCanvas cv("cv", "", width, 800);
  cv.SetMargin(0.13, 0.04, 0.15, 0.04);
  const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin(), T = gPad->GetTopMargin(), B = gPad->GetBottomMargin();
  const double w = padWidthPx(), h = padHeightPx(), plot = 1 - T - B;
  double labelWidth = 0;
  for (const auto &turn : turnOns) labelWidth = std::max(labelWidth, textWidth(turn.text, linePx));
  const double symbol = 44 * scale / w, legendWidth = symbol + labelWidth + 12 / w;
  const double left = L + 36 / w, right = 1 - R - 28 / w;
  /* Above the unity line at the top left: Au+Au the header, p+p the run line with the legend under it. The y range
     ends where that block clears the unity line, and not below the IAN's 1.4 (p+p, a labelled tick on the frame) and
     1.2 (Au+Au). */
  const double markY = 1 - T - (14 + markPx) / h, runY = 1 - T - (14 + linePx) / h;
  const double ly2 = runY - 10 / h, ly1 = ly2 - turnOns.size() * pitch / h;
  const double blockBottom = pp ? ly1 : markY - (lines.size() * pitch + 10) / h;
  const double high = std::max(pp ? 1.4 : 1.2, 1.04 / ((blockBottom - B) / plot));
  studyFrame(0, high, "Trigger / MBD efficiency");
  frameReferenceLine(1);
  TLatex t;
  t.SetNDC();
  t.SetTextFont(43);
  t.SetTextSize(linePx);
  TLegend *leg = nullptr;
  if (pp) {  /* the sPHENIX mark and the system in the free lower right corner */
    t.DrawLatex(left + symbol, runY, runLine.c_str());
    leg = ianLegend(left, ly1, left + legendWidth, ly2, symbol / legendWidth);
    const double markWidth = std::max(textWidth(sPhenixInternalLabel(), markPx), textWidth(system, linePx));
    markAndLines(right - markWidth, B + (30 + pitch) / h, {system}, markPx, linePx, pitch);
  } else {  /* at the lower right, under the plateau, in one column: the plateau values, what is drawn, the legend */
    markAndLines(left, markY, lines, markPx, linePx, pitch);
    const char *curveNote = "Curves: generalized-sigmoid fits", *lineNote = "Dashed lines: 90% of the plateau";
    const double noteWidth = std::max(textWidth(curveNote, linePx), textWidth(lineNote, linePx));
    const double column = right - std::max(labelWidth + 12 / w, noteWidth);
    const double legendTop = B + (26 + turnOns.size() * pitch) / h;
    leg = ianLegend(column - symbol, B + 26 / h, column - symbol + legendWidth, legendTop, symbol / legendWidth);
    double y = legendTop + 18 / h;
    t.DrawLatex(column, y, lineNote);
    t.DrawLatex(column, y += pitch / h, curveNote);
    y += 14 / h;
    for (auto turn = turnOns.rbegin(); turn != turnOns.rend(); ++turn) {
      if (!std::isfinite(turn->plateau)) continue;
      t.SetTextColor(turn->color);
      t.DrawLatex(column, y += pitch / h, Form("%s plateau = %.2f", turn->trigger.c_str(), turn->plateau));
    }
  }
  leg->SetTextSize(linePx);
  for(std::size_t i=0;i<turnOns.size();++i) {
    const auto& turn=turnOns[i];
    const auto trace=turn.eff.physics.then("efficiency_plateau_and_crossing",{{"plateau_from",plateauFrom},{"crossing_share",turnOnShare},{"plateau_method","inverse_variance_mean"},{"crossing_method","linear_between_first_crossing_bin_centres"}});
    photonjet_science::annotation(cv,"trigger/plateau/"+std::to_string(i),photonjet_science::scalar(turn.plateau,"efficiency","inverse variance mean above threshold"),trace);
    photonjet_science::annotation(cv,"trigger/crossing/"+std::to_string(i),std::isfinite(turn.crossing)?photonjet_science::scalar(turn.crossing,"GeV","first interpolated crossing of plateau fraction"):PlotJson{{"state","undefined"},{"reason","no finite crossing"}},trace);
  }
  /* the turn-on point of each trigger, a dashed line up to 90% of its plateau; the curves; the points with their bin
     widths */
  for (const auto &turn : turnOns)
    if (std::isfinite(turn.crossing))
      dashedLine(turn.crossing, 0, turn.crossing, turnOnShare * turn.plateau, turn.color, 3);
  for (const auto &turn : turnOns)
    if (TF1 *curve = turnOnCurve(turn.eff, turn.plateau, turn.color)) photonjet_science::draw(curve,"SAME");
  for (size_t i = 0; i < turnOns.size(); ++i) {
    /* Au+Au as in the IAN: filled triangles, then open circles (with a heavy outline and a little larger, so they
       show on their curve) */
    const int marker = pp ? 20 : i == 0 ? 23 : 89;
    const double size = (marker == 89 ? 1.9 : 1.5) * scale;
    TGraphErrors *g = markers(turnOns[i].eff, turnOns[i].color, marker, true, true, 0, size);
    photonjet_science::draw(g,"P");
    leg->AddEntry(g, turnOns[i].text.c_str(), "p");
  }
  photonjet_science::draw(leg);
  savePhysicsQA(cv, turnOn);
}

/* ---- saved photon-ID models: records and receipts -------------------------------------------------------------- */

/* The training of the photon-ID models (PhotonID) computes the statistics of a model; the qa step imports its tables,
   histograms and JSON records into the QA file under model_qa/, and the figures of the saved models draw them. Where
   a figure works a number out of them (a difference between models, a mean, a maximum), its own comment says so.
   The E_T range a model is trained in is its own, whatever the photon bins of the measurement.
   `index` is, in every function that takes it, the parsed model_qa/collection record: default_model, models,
   comparisons, system and the hash of the collection.
   Inputs: a comparison whose models the file does not hold is listed through missing(); a figure whose optional
   table or record is absent is left out without a listing. A table, histogram or record that a figure needs and
   that is absent, incomplete or inconsistent mostly stops the run with an error instead (the throw statements),
   before FIGURES.txt is written.
   Each PNG has a .json receipt, written by save() and completed by the savers below. Its plot_recipe.function names
   the function that drew the figure. */

/* The diagnostics record of a saved model (model_qa/<name>/meta/diagnostics). The record is large: it is parsed once
   and kept here, and the entry point empties the store before it draws. */
std::map<std::string, PlotJson> gSavedModelMetadata;
const PlotJson &modelMetadata(const std::string &name) {
  const auto cached = gSavedModelMetadata.find(name);
  if (cached != gSavedModelMetadata.end()) return cached->second;
  PlotJson record = savedJson("model_qa/" + name + "/meta/diagnostics");
  if (record.is_null()) throw std::runtime_error("missing saved PhotonID diagnostics metadata: " + name);
  return gSavedModelMetadata.emplace(name, record).first->second;
}

/* True for the record of a p+p model or collection (model metadata or `index`); a record that names no system is
   Au+Au. */
bool isPP(const PlotJson &record) { return record.value("system", "AuAu") == "PP"; }

/* The simulation of a system as a header line names it. */
std::string simulationName(bool pp) { return pp ? "p+p simulation" : "Au+Au embedded simulation"; }

/* The thresholds of the photon+jet samples as a sample line lists them: "12/20". */
std::string thresholdList(const std::set<int> &thresholds) {
  std::string list;
  for (int value : thresholds) list += (list.empty() ? "" : "/") + std::to_string(value);
  return list;
}

/* The thresholds of the photon+jet samples in the saved normalisation record of a model (`review`: its
   application/review record, sample_normalization). A threshold is the number after "Photon" or "photonjet" in the
   key of a sample; the name of the model is not read. */
std::set<int> photonJetThresholds(const PlotJson &review) {
  std::set<int> thresholds;
  for (const auto &sample : review.at("sample_normalization").items()) {
    if (sample.value().at("definition").value("family", "") != "photonjet") continue;
    const auto key = sample.key();
    auto first = key.find("Photon");
    if (first != std::string::npos) first += 6;
    else {
      first = key.find("photonjet");
      if (first != std::string::npos) first += 9;
    }
    const auto last = first == std::string::npos ? first : key.find_first_not_of("0123456789", first);
    if (first == std::string::npos || first == last)
      throw std::runtime_error("missing photon sample threshold: " + key);
    thresholds.insert(std::stoi(key.substr(first, last - first)));
  }
  if (thresholds.empty()) throw std::runtime_error("no photon-family source for simulation label");
  return thresholds;
}

/* The sample line of a saved model, "Au+Au embedded simulation, Photon+jet 12/20": the system of the model and the
   thresholds of the photon+jet samples it was applied to (photonJetThresholds). */
std::string modelPhotonSimulationLabel(const std::string &name) {
  const auto review = savedJson("model_qa/" + name + "/application/review");
  if (review.is_null()) throw std::runtime_error("missing photon sample label provenance");
  const std::set<int> thresholds = photonJetThresholds(review);
  return simulationName(isPP(modelMetadata(name))) + ", Photon+jet " + thresholdList(thresholds);
}

/* A saved table of a model (a tree under model_qa/<name>/) as one record per row. Every column is read as a number,
   the first element of its leaf: a flag such as `valid` is 0 or 1. */
std::vector<PlotJson> modelTable(const std::string &name, const std::string &table) {
  auto *tree = dynamic_cast<TTree *>(gQA->Get(("model_qa/" + name + "/" + table).c_str()));
  if (!tree) throw std::runtime_error("missing saved PhotonID table: " + table);
  std::vector<PlotJson> rows;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    PlotJson row = PlotJson::object();
    TIter next(tree->GetListOfLeaves());
    while (auto *leaf = static_cast<TLeaf *>(next())) row[leaf->GetName()] = leaf->GetValue();
    rows.push_back(row);
  }
  return rows;
}

/* The source file that the receipts of the figures saved here name. */
const char *const kQASource = "plots/QAFigures.C";

/* Reads the receipt of the figure at `path` (the .json beside the PNG), lets `edit` change it and writes it back. */
template <class Edit> void amendReceipt(const std::string &path, Edit edit) {
#ifdef PHOTONJET_NATIVE_CAPTURE_TEST
  if(gNativeCaptureSink)return;
#endif
  const std::string file = gFig + "/" + path + ".json";
  PlotJson receipt = photonjet_render::read(file);
  edit(receipt);
  photonjet_render::write(file, receipt);
}

/* The plot_recipe record of a receipt: the figure of the note that `output` belongs to (`figure`), which version of
   that figure it is (`id` and `role`; `baseline` is the output of the canonical version) and the function that drew
   it. */
PlotJson plotRecipe(const std::string &figure, const std::string &id, const std::string &role,
                    const std::string &output, const std::string &baseline, const std::string &function) {
  return {{"schema", "PhotonJetPlotRecipeV1"},
          {"figure", figure},
          {"id", id},
          {"role", role},
          {"output", output},
          {"baseline_output", baseline},
          {"source", kQASource},
          {"function", function}};
}

/* Adds to the receipt of the figure at `path` the identity of the model collection `index` and of the models `names`
   the figure shows: the hash of the collection, the entries of the models, the import validation of the QA file
   (model_qa/validation) and, for a model that has one, its application review. The per-source records of a review
   stay in the QA file; the receipt keeps their number and their hash. The qualification is SEALED_MODEL_QA only for
   a collection whose import validation passed with the same hash and that is no format test. */
void modelQAReceipt(const std::string &path, const PlotJson &index, const std::vector<std::string> &names) {
  const auto validation = savedJson("model_qa/validation");
  const bool imported = !validation.is_null() && validation.value("status", "") == "PASS" &&
                        validation.at("collection_sha256") == index.at("collection_sha256") &&
                        index.value("qualification", "").find("FORMAT_TEST") == std::string::npos;
  PlotJson models = PlotJson::array();
  for (const auto &name : names)
    for (const auto &entry : index.at("models"))
      if (entry.at("name") == name) models.push_back(entry);
  amendReceipt(path, [&](PlotJson &receipt) {
    receipt["kind"] = "PHOTON_ID_MODEL_QA";
    receipt["annotation_status"] = "SAVED_MODEL_METADATA";
    receipt["model_qa"] = {{"collection_sha256", index.at("collection_sha256")},
                           {"models", models},
                           {"import_validation", validation},
                           {"qualification", imported ? "SEALED_MODEL_QA" : "FORMAT_TEST_ONLY"},
                           {"scientific_acceptance", false}};
    for (const auto &name : names) {
      auto record = savedJson("model_qa/" + name + "/application/review");
      if (record.is_null()) continue;
      if (record.value("schema", "") == "PhotonIDWorkingPointReviewV1" && record.contains("sources")) {
        record["source_count"] = record.at("sources").size();
        record["source_inventory_sha256"] = photonjet_render::sha256(record.at("sources").dump());
        record.erase("sources");
      }
      receipt["model_qa"]["application_review"][name] = record;
    }
    if (isPP(index)) receipt["model_qa"]["comparisons"] = index.at("comparisons");
  });
}

/* Saves a figure of physics-weighted candidates of the models `names`. save() checks the weight record `weights`
   before it writes anything; the receipt then gets the model identity and the origin of the weights (of a
   working-point review record its hash, its number of sources and its sample normalisation in place of the whole
   record). */
void saveModelQA(TCanvas &canvas, const std::string &path, const PlotJson &index, const std::vector<std::string> &names,
                 const PlotJson &weights, const char *function) {
  gSemanticPackages.clear();
  save(canvas, path, weights, function, kQASource);
  modelQAReceipt(path, index, names);
  amendReceipt(path, [&](PlotJson &receipt) {
    if (weights.value("schema", "") == "PhotonIDWorkingPointReviewV1")
      receipt["weight_provenance"] = {{"schema", "PhotonJetVerifiedNativePlotWeightsV1"},
                                      {"role", "physics_event_source"},
                                      {"full_record_sha256", photonjet_render::sha256(weights.dump())},
                                      {"source_count", weights.at("sources").size()},
                                      {"sample_normalization", weights.at("sample_normalization")},
                                      {"training_weights_used", false},
                                      {"record_in_input_root", "model_qa/" + names.at(0) + "/application/review"}};
    else receipt["weight_provenance"] = weights;
  });
}

/* Saves a figure of the released working points of the model `name` (model_qa/<name>/meta/working_points): cut
   values and calibration quantiles, which describe the learned model and are no weighted candidate yields. The
   receipt carries the hash of the working-point record, the model identity, the saved derivation where the file has
   one, and a recipe that follows from the file name: the threshold grid (pt_centrality_grid, drawn by workingPoints,
   with the centrality comparison as its supplement) or else the centrality lines (savedWorkingPointFits).
   workingPointCentralityComparison saves through this function and then writes its own recipe over that one. */
void saveReleasedWorkingPoints(TCanvas &canvas, const std::string &path, const PlotJson &index,
                               const std::string &name, const char *function) {
  const auto parameters = savedJson("model_qa/" + name + "/meta/working_points");
  if (parameters.is_null()) throw std::runtime_error("missing released working-point parameters: " + name);
  const bool grid = path.find("pt_centrality_grid") != std::string::npos;
  gSemanticPackages.clear();
  save(canvas, path, {{"role", "released_cut_parameters"},
                      {"parameter_kind", "released_working_points"},
                      {"working_points_sha256", photonjet_render::sha256(parameters.dump())}}, function, kQASource);
  modelQAReceipt(path, index, {name});
  amendReceipt(path, [&](PlotJson &receipt) {
    receipt["kind"] = "RELEASED_WORKING_POINT_PARAMETERS";
    receipt["working_points_sha256"] = photonjet_render::sha256(parameters.dump());
    receipt["plot_recipe"] = plotRecipe(grid ? "auau_bdt_threshold_grid" : "auau_bdt_working_point_fits", "canonical",
                                        "canonical", path, path, grid ? "workingPoints" : "savedWorkingPointFits");
    if (grid)
      receipt["supplemental_outputs"] = PlotJson::array(
        {{{"id", "centrality_extremes"},
          {"title", "0–10% and 75–80% centrality"},
          {"output", "photonOnly/auau_h70_bdt_centrality_comparison.png"},
          {"function", "workingPointCentralityComparison"},
          {"position", "after"}}});
    const auto derivation = savedJson("model_qa/" + name + "/meta/working_point_derivation");
    if (!derivation.is_null()) receipt["working_point_derivation"] = derivation;
    receipt["parameter_panel"] = grid ? "fine_bin_thresholds" : "released_centrality_fit";
  });
}

/* What the receipt of each training diagnostic states: how what it shows is weighted, and the sample it shows. */
struct TrainingDiagnostic { const char *weighting, *sample; };
const std::map<std::string, TrainingDiagnostic> kTrainingDiagnostics = {
  {"training_reweighting", {"raw_and_training_weighted", "prepared sample (training and held-out sets)"}},
  {"feature_separation", {"training_weighted", "prepared sample (training and held-out sets)"}},
  {"split_gain", {"learned_model_split_gain", "prepared sample (training and held-out sets)"}},
  {"model_routing", {"training_weighted", "common held-out population"}},
  {"training_capacity", {"training_weighted", "prepared sample (training and held-out sets)"}},
  {"paired_score_mapping", {"unit_counts", "paired held-out test photons"}},
  {"paired_score_difference", {"unit_counts", "paired held-out test photons"}}};

/* Saves a figure that shows how the models `names` were trained. `diagnostic` is one of the keys of
   kTrainingDiagnostics; it names what the figure shows and sets the weighting and the sample that the receipt states.
   Such a figure shows the sample and the weighting its kTrainingDiagnostics row names: the weights describe the
   learning and are no physics yields. save() checks the record made here before it writes anything, and
   refuses the figure unless a line on the canvas reads "Training diagnostic". The receipt names the function that
   calls this one as the function that drew the figure. */
void saveTrainingDiagnostic(TCanvas &canvas, const std::string &path, const PlotJson &index,
                            const std::vector<std::string> &names, const char *diagnostic,
                            const char *function = __builtin_FUNCTION()) {
  const std::string kind = diagnostic;
  const auto stated = kTrainingDiagnostics.find(kind);
  if (stated == kTrainingDiagnostics.end()) throw std::runtime_error("unknown training diagnostic: " + kind);
  const bool pairs = kind == "paired_score_mapping" || kind == "paired_score_difference";
  PlotJson models = PlotJson::array();
  for (const auto &name : names) {
    const auto &meta = modelMetadata(name);
    const auto &identity = meta.at("identity");
    models.push_back({{"name", name},
                      {"population_sha256", identity.at("population_sha256")},
                      {"model_artifacts_sha256", identity.at("model_artifacts_sha256")},
                      {"training_weights_sha256", identity.at("training_weights_sha256")},
                      {"event_split_sha256", identity.at("event_split_sha256")},
                      {"split_sha256", identity.at("split_sha256")}});
  }
  PlotJson proof = {{"schema", "PhotonIDTrainingDiagnosticV1"},
                    {"role", "model_training_diagnostic"},
                    {"diagnostic", diagnostic},
                    {"physics_distribution", false},
                    {"weighting", stated->second.weighting},
                    {"dataset_role", stated->second.sample},
                    {"models", models}};
  if (pairs) {  /* a paired comparison with the original classifier: its model and the number of pairs */
    const auto &comparison = modelMetadata(names.at(0)).at("original_comparison");
    proof["reference_model_sha256"] = comparison.at("reference_model_sha256");
    proof["paired_rows"] = comparison.at("rows");
  }
  gSemanticPackages.clear();
  save(canvas, path, proof, function, kQASource);
  modelQAReceipt(path, index, names);
  amendReceipt(path, [&](PlotJson &receipt) {
    receipt["kind"] = "MODEL_TRAINING_DIAGNOSTIC";
    receipt["weight_provenance"] = proof;
  });
}

/* ---- working points -------------------------------------------------------------------------------------------- */

/* The figures of this section and the score-region figure further down are laid out in canvas pixels counted from the
   top left of the canvas: pixelPanel places a frame, pixelLabel a text, pixelDashes a dashed line and pixelKey a
   legend entry. */

/* The 70, 80 and 90% working points as every working-point figure draws them: filled circles (kWorkingPointMarker)
   in blue, green and vermilion (cOrange). Each entry: the efficiency in percent and its colour. The green is a bluish
   green of this section, beside the green of the palette. */
const int kWorkingPointMarker = 20;
std::array<std::pair<int, int>, 3> workingPointLevels() {
  return {{{70, cBlue}, {80, TColor::GetColor("#009E73")}, {90, cOrange}}};
}

/* The sample line of the Au+Au working-point figures, "Au+Au embedded simulation, Photon+jet 12/20": the photon+jet
   samples behind the saved model `model` or, without a model, behind the Au+Au package of the QA file (from its
   sample records: the number in the name of each photon+jet sample, and "embedded" where a record says so). */
std::string photonJetSampleLine(const std::string &model) {
  if (!model.empty()) return modelPhotonSimulationLabel(model);
  std::set<int> thresholds;
  bool embedded = false;
  if (!gQAMetadata.is_null() && gQAMetadata.at("packages").contains("AuAu")) {
    const auto samples = gQAMetadata.at("packages").at("AuAu").value("samples", PlotJson::object());
    for (const auto &sample : samples.items()) {
      if (!sample.value().is_object() || sample.value().value("family", "") != "photonjet") continue;
      const auto digits = sample.key().find_first_of("0123456789");
      if (digits != std::string::npos) thresholds.insert(std::stoi(sample.key().substr(digits)));
      const auto normalization = sample.value().value("normalization", PlotJson::object());
      embedded = embedded || normalization.value("source_factor", "").find("embedded") != std::string::npos;
    }
  }
  if (thresholds.empty()) return "Au+Au simulation";
  return (embedded ? simulationName(false) : "Au+Au simulation") + ", Photon+jet " + thresholdList(thresholds);
}

/* The distance in pixels from an axis of a panel placed by pixelPanel to its tick labels. */
const double kPixelLabelGap = 9;

/* One panel of a figure laid out in canvas pixels: a transparent pad over the whole canvas whose frame is the pixel
   rectangle from (x0, y0) to (x1, y1), y counted from the top of the canvas, with tick labels of `labelPx` and ticks
   `tickPx` long. Panels placed this way have exactly equal frames whatever labels they carry, no label is cut at a
   panel's edge, and every text of the figure is placed in canvas pixels. Returns the frame over the given axis
   ranges, without titles. */
TH1F *pixelPanel(TCanvas &canvas, double x0, double y0, double x1, double y1, double xlow, double ylow, double xhigh,
                 double yhigh, double labelPx, double tickPx = 3) {
  canvas.cd();
  const double w = canvas.GetWw(), h = canvas.GetWh();
  auto *pad = new TPad(fresh().c_str(), "", 0, 0, 1, 1);
  pad->SetFillStyle(4000);
  pad->SetMargin(x0 / w, 1 - x1 / w, 1 - y1 / h, y0 / h);
  pad->SetTicks(1, 1);
  photonjet_science::draw(pad);
  pad->cd();
  TH1F *axes = photonjet_science::drawFrame(*pad,xlow, ylow, xhigh, yhigh);
  for (TAxis *a : {axes->GetXaxis(), axes->GetYaxis()})
    a->SetLabelFont(43), a->SetLabelSize(labelPx), a->SetTitleSize(0);
  tickLengths(axes, tickPx);
  axes->GetXaxis()->SetLabelOffset(kPixelLabelGap / h), axes->GetYaxis()->SetLabelOffset(kPixelLabelGap / w);
  return axes;
}

/* The division code of a linear axis that spans `range`: major ticks at the smallest of 1, 2 and 5 times a power of
   ten that gives at most `most` intervals, each divided in five, or in four when the step is 2. */
int axisDivisions(double range, double most) {
  double step = std::pow(10, std::floor(std::log10(range / most)));
  int minor = 5;
  for (double factor : {1.0, 2.0, 5.0, 10.0})
    if (range / (factor * step) <= most) {
      step *= factor;
      minor = factor == 2.0 ? 4 : 5;
      break;
    }
  return 100 * minor + int(std::ceil(range / step - 1e-9));
}

/* Text at the canvas pixel (x, y), y counted from the top, in a pad that covers the canvas: `px` high, regular or
   bold, with ROOT's alignment code and rotation. */
void pixelLabel(double x, double y, const std::string &text, double px, int align = 11, int color = cBlack,
                bool bold = false, double angle = 0) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(bold ? 63 : 43);
  t.SetTextSize(px);
  t.SetTextAlign(align);
  t.SetTextColor(color);
  t.SetTextAngle(angle);
  t.DrawLatex(x / padWidthPx(), 1 - y / padHeightPx(), text.c_str());
}

/* The upper edge of a linear y range from `low` that keeps `value` under a block of text which reaches `blockPx`
   below the top of a frame `framePx` high. */
double raisedTop(double low, double value, double blockPx, double framePx) {
  return low + (value - low) / (1 - blockPx / framePx);
}

/* The working-point thresholds against centrality with their straight lines, drawn on the current canvas as in the IAN:
   the sPHENIX mark and the sample lines `headerLines` (and the input note, when there is one) at the top left, the
   relation and the line of each working point at the top right in its colour, and a y range opened until neither
   block lies over a point or a line. thresholds: (centrality, threshold) per working point; lines: intercept and slope
   (NaN without a saved line); digits: the decimals of the intercept. applied: the lines the analysis applies, from the
   selection record of the packages, each with the name the legend gives it ("applied in the analysis"). A working
   point is its signal efficiency and its line follows from the simulation at hand, so an applied line is expected to
   be one of the lines derived here. It is that line when both parameters agree within one unit of the last digit the
   legend prints: the level is then drawn as the applied line and its entry says so. An applied line that agrees with
   none is drawn black and dashed with its own entry: a line taken from another pass is then seen, never hidden. */
void workingPointFitPanel(double low, double high,
                          const std::array<std::vector<std::pair<double, double>>, 3> &thresholds,
                          const std::array<std::pair<double, double>, 3> &derived, const std::string &ytitle,
                          const std::string &relation, int digits, const std::vector<std::string> &headerLines,
                          const std::string &note,
                          const std::vector<std::pair<std::pair<double, double>, std::string>> &applied = {},
                          const std::array<photonjet_science::Trace,3>& traces = {}) {
  const double w = padWidthPx(), h = padHeightPx();
  const double left = 142, right = 30, top = 20, bottom = 116, fw = w - left - right, fh = h - top - bottom;
  const double markPx = 40, linePx = 28, step = 42, inset = 30, symbolPx = 84;
  gPad->SetMargin(left / w, right / w, bottom / h, top / h);
  gPad->SetTicks(1, 1);
  const auto levels = workingPointLevels();
  /* a line as the legend prints it */
  const auto printed = [&](const std::pair<double, double> &line) {
    return std::string(Form("%.*f %s %.5f C", digits, line.first, line.second < 0 ? "#minus" : "+",
                            std::fabs(line.second)));
  };
  /* The applied lines: each is the line of the level it agrees with, or a line of its own (`own`: line, entry). */
  std::array<std::pair<double, double>, 3> lines = derived;  /* as drawn: a level an applied line takes is that line */
  std::array<std::string, 3> appliedAs;                       /* the name of the applied line a level stands for */
  std::vector<std::pair<std::pair<double, double>, std::string>> own;
  const double interceptStep = std::pow(10.0, -digits), slopeStep = 1e-5;
  for (const auto &[line, name] : applied) {
    if (!std::isfinite(line.first) || !std::isfinite(line.second)) continue;
    int level = -1;
    for (size_t i = 0; i < levels.size() && level < 0; ++i)
      if (appliedAs[i].empty() && std::isfinite(derived[i].first) &&
          std::fabs(derived[i].first - line.first) <= interceptStep &&
          std::fabs(derived[i].second - line.second) <= slopeStep)
        level = i;
    if (level >= 0) {
      lines[level] = line, appliedAs[level] = name;
    } else {
      std::string capital = name;
      if (!capital.empty()) capital[0] = std::toupper(capital[0]);
      own.push_back({line, capital + ": " + printed(line)});
    }
  }
  /* the two blocks: their widths and how far below the top of the frame they reach */
  std::vector<std::string> entries;
  double entryWidthPx = textWidth(relation, linePx + 4) * w;
  for (size_t i = 0; i < levels.size(); ++i) {
    const int level = std::get<0>(levels[i]);
    entries.push_back(!std::isfinite(lines[i].first) ? std::string(Form("%d%%", level))
                      : appliedAs[i].empty() ? Form("%d%%: ", level) + printed(lines[i])
                      : Form("%d%% (%s): ", level, appliedAs[i].c_str()) + printed(lines[i]));
    entryWidthPx = std::max(entryWidthPx, textWidth(entries.back(), linePx) * w);
  }
  for (const auto &[line, entry] : own) entryWidthPx = std::max(entryWidthPx, textWidth(entry, linePx) * w);
  const int ownStyles[] = {kDashed, 9, kDashDotted};  /* dashed, long dashes, dash-dot */
  const double legendWidthPx = symbolPx + entryWidthPx + 10, legendTop = 58,
               legendDepth = legendTop + (levels.size() + own.size()) * step + 4;
  double headerWidthPx = textWidth(sPhenixInternalLabel(), markPx) * w;
  for (const auto &l : headerLines) headerWidthPx = std::max(headerWidthPx, textWidth(l, linePx) * w);
  if (!note.empty()) headerWidthPx = std::max(headerWidthPx, textWidth(note, linePx) * w);
  const double headerDepth = 18 + markPx + (headerLines.size() + (note.empty() ? 0 : 1)) * step;
  /* the y range: the points and the lines, opened at the top for the two blocks */
  double lowest = std::numeric_limits<double>::infinity(), highest = -lowest;
  for (size_t i = 0; i < levels.size(); ++i) {
    for (const auto &p : thresholds[i]) lowest = std::min(lowest, p.second), highest = std::max(highest, p.second);
    if (!std::isfinite(lines[i].first)) continue;
    for (double x : {low, high}) {
      const double y = lines[i].first + lines[i].second * x;
      lowest = std::min(lowest, y), highest = std::max(highest, y);
    }
  }
  for (const auto &[line, entry] : own)
    for (double x : {low, high}) {
      const double y = line.first + line.second * x;
      lowest = std::min(lowest, y), highest = std::max(highest, y);
    }
  if (!std::isfinite(lowest)) lowest = 0, highest = 1;
  /* the highest point or line over the stretch [x0, x1] of the frame, in pixels from its left edge */
  const auto under = [&](double x0, double x1) {
    double y = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < levels.size(); ++i) {
      for (const auto &p : thresholds[i]) {
        const double x = (p.first - low) / (high - low) * fw;
        if (x > x0 - 14 && x < x1 + 14) y = std::max(y, p.second);
      }
      if (!std::isfinite(lines[i].first)) continue;
      for (double x : {std::max(x0, 0.0), std::min(x1, fw)})
        y = std::max(y, lines[i].first + lines[i].second * (low + x / fw * (high - low)));
    }
    for (const auto &[line, entry] : own)
      for (double x : {std::max(x0, 0.0), std::min(x1, fw)})
        y = std::max(y, line.first + line.second * (low + x / fw * (high - low)));
    return y;
  };
  const double span = std::max(highest - lowest, 1e-3), yLow = lowest - 0.09 * span;
  /* a block stays 18 px clear of what is drawn under it */
  const auto opened = [&](double y, double depth) { return raisedTop(yLow, y, depth + 18, fh); };
  const double yHigh = std::max({highest + 0.06 * span,
                                 opened(under(inset - 8, inset + headerWidthPx + 8), headerDepth),
                                 opened(under(fw - 24 - legendWidthPx - 8, fw), legendDepth)});
  TH1F *axes = photonjet_science::drawFrame(*gPad,low, yLow, high, yHigh);
  for (TAxis *a : {axes->GetXaxis(), axes->GetYaxis()}) a->SetLabelFont(43), a->SetLabelSize(32), a->SetTitleSize(0);
  axes->GetXaxis()->SetNdivisions(axisDivisions(high - low, 9));
  axes->GetYaxis()->SetNdivisions(axisDivisions(yHigh - yLow, 6.5));
  tickLengths(axes, 14);
  axes->GetXaxis()->SetLabelOffset(10 / h), axes->GetYaxis()->SetLabelOffset(10 / w);
  pixelLabel(left + fw / 2, h - 40, "Centrality percentile [%]", 40, 22);
  pixelLabel(38, top + fh / 2, ytitle, 40, 22, cBlack, false, 90);
  /* the lines under the points, and the legend: each line's parameters in the colour of its working point */
  for (size_t i = 0; i < levels.size(); ++i) {
    if (!std::isfinite(lines[i].first)) continue;
    photonjet_science::annotation(*gPad,"working_point/line/"+std::to_string(i),
      {{"intercept",lines[i].first},{"slope",lines[i].second},{"range",{low,high}},{"efficiency_percent",std::get<0>(levels[i])}},
      traces[i].then("display_linear_working_point"));
    TLine l;
    l.SetLineColor(std::get<1>(levels[i])), l.SetLineWidth(4);
    l.DrawLine(low, lines[i].first + lines[i].second * low, high, lines[i].first + lines[i].second * high);
  }
  for (size_t i = 0; i < own.size(); ++i) {
    const auto &line = own[i].first;
    photonjet_science::annotation(*gPad,"working_point/applied/"+std::to_string(i),
      {{"intercept",line.first},{"slope",line.second},{"range",{low,high}}},
      traces[0].withParameter("photon_id/applied_working_point/"+std::to_string(i),"selection",{line.first,line.second}));
    TLine l;
    l.SetLineColor(cBlack), l.SetLineWidth(4), l.SetLineStyle(ownStyles[i % 3]);
    l.DrawLine(low, line.first + line.second * low, high, line.first + line.second * high);
  }
  const double legendX = left + fw - 24 - legendWidthPx;
  TLegend *legend = ianLegend(legendX / w, 1 - (top + legendDepth - 4) / h, (legendX + legendWidthPx) / w,
                              1 - (top + legendTop) / h, symbolPx / legendWidthPx);
  legend->SetTextSize(linePx);
  for (size_t i = 0; i < levels.size(); ++i) {
    const int color = std::get<1>(levels[i]);
    if (thresholds[i].empty()) continue;
    auto *g = new TGraph;
    for (const auto &p : thresholds[i]) g->SetPoint(g->GetN(), p.first, p.second);
    photonjet_science::bindNative(g,"working_point/thresholds/"+std::to_string(i),traces[i],
      photonjet_science::covariance("not_computed",{},0,"saved quantile thresholds have no covariance"));
    style(g, color, kWorkingPointMarker, 1.9);
    g->SetLineWidth(4);
    photonjet_science::draw(g,"P");
    legend->AddEntry(g, entries[i].c_str(), std::isfinite(lines[i].first) ? "lp" : "p")->SetTextColor(color);
  }
  for (size_t i = 0; i < own.size(); ++i) {
    auto *key = new TLine;
    key->SetLineColor(cBlack), key->SetLineWidth(4), key->SetLineStyle(ownStyles[i % 3]);
    legend->AddEntry(key, own[i].second.c_str(), "l")->SetTextColor(cBlack);
  }
  photonjet_science::draw(legend);
  pixelLabel(legendX, top + 44, relation, linePx + 4);
  double y = top + 16 + 0.8 * markPx;
  pixelLabel(left + inset, y, sPhenixInternalLabel(), markPx);
  for (const auto &l : headerLines) pixelLabel(left + inset, y += step, l, linePx);
  if (!note.empty()) pixelLabel(left + inset, y += step, note, linePx, 11, gGray);
}

/* A dashed line of a figure laid out in canvas pixels, `px` canvas pixels long: solid pieces twice as long as the gaps
   between them, a whole piece at either end. Drawn piece by piece, the dashes keep their proportions at every size the
   figure is written at. The ends are given in the axes of the current pad or, with `pixels`, in canvas pixels counted
   from the top left of a pad that covers the canvas. `px` is given by the caller and must be the length of that
   segment on the canvas: it sets the number of pieces. */
void pixelDashes(double x0, double y0, double x1, double y1, double px, int color, bool pixels = false) {
  const int width = 3;
  const double dash = 4 * width, gap = 2 * width;
  const int pieces = std::max(1, int(std::lround((px + gap) / (dash + gap))));
  const double piece = 1 / (1.5 * pieces - 0.5), w = padWidthPx(), h = padHeightPx();
  TLine l;
  l.SetLineColor(color), l.SetLineWidth(width);
  for (int i = 0; i < pieces; ++i) {
    const double from = 1.5 * i * piece, to = from + piece;
    const double xa = x0 + from * (x1 - x0), ya = y0 + from * (y1 - y0);
    const double xb = x0 + to * (x1 - x0), yb = y0 + to * (y1 - y0);
    if (pixels) l.DrawLineNDC(xa / w, 1 - ya / h, xb / w, 1 - yb / h);
    else l.DrawLine(xa, ya, xb, yb);
  }
}

/* One legend entry of a figure whose summary lines are dashed, at the canvas pixel (x, y), y counted from the top, in
   the frame of a panel placed by pixelPanel: the marker on a piece of its dashed line, `keyPx` long from x, and the
   text after it, centred on the line. */
void pixelKey(double x, double y, double keyPx, int color, int marker, double size, const std::string &text,
              double px) {
  pixelDashes(x, y, x + keyPx, y, keyPx, color, true);
  /* the marker is a graph point: its canvas pixel in the axes of the panel */
  const double w = padWidthPx(), h = padHeightPx(), fx = gPad->GetLeftMargin() * w, fy = gPad->GetTopMargin() * h;
  const double fw = w - fx - gPad->GetRightMargin() * w, fh = h - fy - gPad->GetBottomMargin() * h;
  auto *g = new TGraph;
  g->SetPoint(0, gPad->GetUxmin() + (x + keyPx / 2 - fx) / fw * (gPad->GetUxmax() - gPad->GetUxmin()),
              gPad->GetUymax() - (y - fy) / fh * (gPad->GetUymax() - gPad->GetUymin()));
  style(g, color, marker, size);
  photonjet_science::decoration(g,"working point legend marker");
  photonjet_science::draw(g,"P");
  pixelLabel(x + keyPx + 10, y, text, px, 12);
}

/* The photon-ID score (bdt) or isolation (isolation) thresholds of the 70, 80 and 90% working points. First figure: the
   thresholds against photon E_T in one panel per saved centrality slice, four panels to a row, equal frames with
   shared axes and each slice's median as a dashed line. Second figure: the medians against centrality with the
   straight lines of the qa step and the lines the analysis applies, read from the selection record of the Au+Au
   packages (the isolation line; for the score, the explicit pair of the fill): each named in the legend entry of the
   derived line it agrees with, or drawn black and dashed when it agrees with none.
   Read from qa/AuAu/working_points/<what>: thresholds_<level> (centrality slice against photon E_T), median_<level>
   (per slice) and line_<level> (bin 1: the intercept, bin 2: the slope). A file without that directory but with a
   saved model gives the score thresholds of its default model (model_qa/<model>/application) and the first figure
   only; savedWorkingPointFits draws the released lines of the model. */
void workingPoints(const std::string &what) {
  const bool bdt = what == "bdt";
  const double radius = bdt ? 0.0 : isolationRadius("AuAu");
  const std::string radiusTag = Form("%02d", int(std::round(radius * 10)));
  const std::string gridPng = bdt ? "photonOnly/auau_h70_bdt_pt_centrality_grid.png"
                                  : "isolation/auau_r" + radiusTag + "_isolation_pt_centrality_grid.png";
  const std::string fitsPng = bdt ? "photonOnly/auau_h70_bdt_centrality_fits.png"
                                  : "isolation/auau_r" + radiusTag + "_isolation_centrality_fits.png";
  TDirectory *dir = gQA->GetDirectory(("qa/AuAu/working_points/" + what).c_str());
  PlotJson modelIndex;
  std::string modelName;
  if (!dir && bdt) {
    modelIndex = savedJson("model_qa/collection");
    if (!modelIndex.is_null()) {
      modelName = modelIndex.at("default_model");
      dir = gQA->GetDirectory(("model_qa/" + modelName + "/application").c_str());
    }
  }
  if (!dir) {
    missing(gridPng, "no Au+Au photon-family simulation QA");
    /* savedWorkingPointFits draws the lines of a saved model under the same file name: the figure is listed as not
       drawn only while no figure of that name is written */
    if (std::find(gWritten.begin(), gWritten.end(), fitsPng) == gWritten.end())
      missing(fitsPng, "no Au+Au photon-family simulation QA");
    return;
  }
  const std::string ytitle = bdt ? "BDT score threshold" : Form("R = %g isolation threshold [GeV]", radius);
  const auto levels = workingPointLevels();
  auto *thresholds70 = dynamic_cast<TH2 *>(dir->Get("thresholds_70"));
  const int slices = thresholds70 ? thresholds70->GetNbinsX() : 0;
  if (!slices) {
    missing(gridPng, "missing threshold grid");
    return missing(fitsPng, "missing threshold grid");
  }
  const double ptLow = thresholds70->GetYaxis()->GetXmin(), ptHigh = thresholds70->GetYaxis()->GetXmax();
  auto *median70 = dynamic_cast<TH1 *>(dir->Get("median_70"));
  if (!median70 || median70->GetNbinsX() != slices) throw std::runtime_error("working-point centrality axes differ");
  const TAxis *sliceAxis = median70->GetXaxis();  /* the centrality slices, for their edges */
  const std::string sample = photonJetSampleLine(modelName);
  const std::string simulation = modelName.empty() ? simulationNote({"AuAu"}) : "";
  /* per slice */
  {
    /* One layout for the score and the isolation thresholds, from a saved model or from the QA file: equal frames on
       shared axes, the centrality range above each panel, the axis titles once on the canvas. The first row carries
       the sPHENIX mark with the sample line at the top left of its first panel and the legend at the top of its last;
       its y range opens until both stand clear of the points, and every other row ends just above its points. */
    const int rows = (slices + 3) / 4, outer = std::min(slices, 4);
    TCanvas cv("cv", "", 1800, 115 + 265 * rows);
    const double w = cv.GetWw(), h = cv.GetWh();
    /* The sizes of the axis text are chosen first and the margins follow from them, so that large tick labels and
       axis titles stand clear of each other and of the frames. Left: the y title, centred on yTitleX, then a tick
       label of four characters. Bottom: the x labels, then the x title, centred xTitleUp above the lower edge of the
       canvas. Between two frames of a row: the halves of the last x label of one and of the first of the next, and
       most of a label's height between them, so that the two read as two numbers. Above every row: the centrality
       range of its panels. The titles have the largest size at which the written image sets an index exactly on its
       symbol (exactTextPx). The ticks are tickPx long, a length that shows on a frame of this size. */
    const double labelPx = 32, titlePx = 46, rangePx = 32, tickPx = 8;
    const double firstLabel = textWidth(Form("%g", ptLow), labelPx) * w,
                 lastLabel = textWidth(Form("%g", ptHigh), labelPx) * w;
    const double yTitleX = 10 + titlePx / 2, xTitleUp = 12 + 0.6 * titlePx;
    const double left = yTitleX + titlePx / 2 + 6 + textWidth("0.00", labelPx) * w + kPixelLabelGap,
                 right = lastLabel / 2 + 8, top = 16 + 1.1 * rangePx,
                 bottom = kPixelLabelGap + 0.85 * labelPx + 18 + titlePx / 2 + xTitleUp,
                 gapX = (firstLabel + lastLabel) / 2 + 0.8 * labelPx, gapY = 24 + 1.1 * rangePx;
    const double fw = (w - left - right - 3 * gapX) / 4, fh = (h - top - bottom - (rows - 1) * gapY) / rows;
    /* the sample line at sampleFullPx, or smaller where it would fill more than the width of a frame */
    const double sampleFullPx = 17;
    const double markPx = kMarkPx,
                 samplePx = std::min(sampleFullPx, sampleFullPx * (fw - 40) / (textWidth(sample, sampleFullPx) * w)),
                 notePx = 15;
    const double headerDepth = 14 + markPx + 8 + samplePx + (simulation.empty() ? 0 : 6 + notePx);
    const auto highest = [&](int slice) {  /* the highest threshold of a slice */
      double y = -std::numeric_limits<double>::infinity();
      for (const auto &level : levels)
        if (auto *t = dynamic_cast<TH2 *>(dir->Get(Form("thresholds_%d", std::get<0>(level)))))
          for (int b = 1; b <= t->GetNbinsY(); ++b)
            if (std::isfinite(t->GetBinContent(slice, b))) y = std::max(y, t->GetBinContent(slice, b));
      return y;
    };
    const int markPanel = 1, legendPanel = outer;
    const double legendTop = markPanel == legendPanel ? headerDepth + 6 : 10, legendDepth = legendTop + 32;
    const double keyPx = 42, entryPx = 22;
    for (int row = 0; row < rows; ++row) {
      Span s;
      for (const auto &level : levels)
        if (auto *t = dynamic_cast<TH2 *>(dir->Get(Form("thresholds_%d", std::get<0>(level)))))
          for (int slice = 4 * row + 1; slice <= std::min(slices, 4 * row + 4); ++slice)
            for (int etBin = 1; etBin <= t->GetNbinsY(); ++etBin) s.add(t->GetBinContent(slice, etBin));
      const bool any = std::isfinite(s.high);
      const double span = any ? std::max(s.high - s.low, 1e-3) : 1, yLow = (any ? s.low : 0) - 0.12 * span;
      double yHigh = (any ? s.high : 1) + 0.14 * span;
      if (row == 0) {
        /* a block stays 16 px clear of what is drawn under it */
        const auto opened = [&](double y, double depth) { return raisedTop(yLow, y, depth + 16, fh); };
        yHigh = std::max({yHigh, opened(highest(markPanel), headerDepth), opened(highest(legendPanel), legendDepth)});
      }
      for (int col = 0; col < 4 && 4 * row + col < slices; ++col) {
        const int slice = 4 * row + col + 1;
        const double fx = left + col * (fw + gapX), fy = top + row * (fh + gapY);
        TH1F *axes = pixelPanel(cv, fx, fy, fx + fw, fy + fh, ptLow, yLow, ptHigh, yHigh, labelPx, tickPx);
        axes->GetXaxis()->SetNdivisions(505);
        axes->GetYaxis()->SetNdivisions(axisDivisions(yHigh - yLow, 6.2));
        axes->GetYaxis()->SetDecimals();
        if (row + 1 != rows) axes->GetXaxis()->SetLabelSize(0);
        if (col != 0) axes->GetYaxis()->SetLabelSize(0);
        const std::string range = centralityRange(sliceAxis->GetBinLowEdge(slice), sliceAxis->GetBinUpEdge(slice));
        pixelLabel(fx + fw / 2, fy - 0.4 * rangePx, range, rangePx, 21, cBlack, true);
        /* the dashed medians under the points */
        for (const auto &[level, color] : levels) {
          auto *median = dynamic_cast<TH1 *>(dir->Get(Form("median_%d", level)));
          if (!median || !std::isfinite(median->GetBinContent(slice))) continue;
          photonjet_science::annotation(*gPad,"working_point/median/"+std::to_string(slice)+"/"+std::to_string(level),
            photonjet_science::scalar(median->GetBinContent(slice),bdt?"score":"GeV","saved median"),
            photonjet_science::selectedTrace(median,sliceAxis->GetBinLowEdge(slice),sliceAxis->GetBinUpEdge(slice)));
          pixelDashes(ptLow, median->GetBinContent(slice), ptHigh, median->GetBinContent(slice), fw, color);
        }
        for (const auto &[level, color] : levels) {
          auto *t = dynamic_cast<TH2 *>(dir->Get(Form("thresholds_%d", level)));
          if (!t) continue;
          auto *g = new TGraph;
          for (int etBin = 1; etBin <= t->GetNbinsY(); ++etBin) {
            const double value = t->GetBinContent(slice, etBin);
            if (std::isfinite(value)) g->SetPoint(g->GetN(), t->GetYaxis()->GetBinCenter(etBin), value);
          }
          const auto trace=photonjet_science::projectionOf(t,photonjet_science::window(0,t->GetXaxis()->GetBinLowEdge(slice),t->GetXaxis()->GetBinUpEdge(slice)));
          photonjet_science::bindNative(g,"working_point/grid/"+std::to_string(slice)+"/"+std::to_string(level),trace,
            photonjet_science::covariance("not_computed",{},0,"saved threshold grid has no covariance"));
          style(g, color, kWorkingPointMarker, 1.3);
          if (g->GetN()) photonjet_science::draw(g,"P");
        }
        if (row == 0 && slice == legendPanel) {
          /* the legend: each marker on its dashed line, the three entries spread over the width of the frame */
          double entryWidth = 0;
          for (const auto &level : levels)
            entryWidth = std::max(entryWidth, textWidth(Form("%d%%", std::get<0>(level)), entryPx) * w);
          const double each = keyPx + 10 + entryWidth, between = (fw - 40 - levels.size() * each) / (levels.size() - 1);
          for (size_t i = 0; i < levels.size(); ++i) {
            const auto [level, color] = levels[i];
            pixelKey(fx + 20 + i * (each + between), fy + legendTop + 16, keyPx, color, kWorkingPointMarker, 1.3,
                     Form("%d%%", level), entryPx);
          }
        }
        if (row == 0 && slice == markPanel) {
          double y = fy + 14 + 0.78 * markPx;
          pixelLabel(fx + 22, y, sPhenixInternalLabel(), markPx);
          pixelLabel(fx + 22, y += 8 + 0.86 * samplePx, sample, samplePx);
          if (!simulation.empty()) pixelLabel(fx + 22, y += 6 + notePx, simulation, notePx, 11, gGray);
        }
      }
    }
    cv.cd();
    /* the axis titles once, centred on the grid; a y title longer than the grid is high is written smaller, to fit */
    const double gridHeight = h - top - bottom;
    const double yTitlePx = std::min(titlePx, titlePx * 0.96 * gridHeight / (textWidth(ytitle, titlePx) * w));
    pixelLabel(left + (w - left - right) / 2, h - xTitleUp, "Photon E_{T} [GeV]", titlePx, 22);
    pixelLabel(yTitleX, top + gridHeight / 2, ytitle, yTitlePx, 22, cBlack, false, 90);
    pixelLabel(w - right, h - 30, "Y scale shared within each row", 19, 31, gGray);
    if (!modelName.empty()) saveReleasedWorkingPoints(cv, gridPng, modelIndex, modelName);
    else {
      renderSelection("AuAu", true);  /* the threshold grid reads the Au+Au QA directly: name its inputs */
      savePhysicsQA(cv, gridPng);
    }
  }
  /* A saved model has its released lines, which are not refitted from the finer display grid: savedWorkingPointFits
     draws them. */
  if (!modelName.empty()) return;
  /* the slice medians against centrality, with the saved straight lines */
  std::array<std::vector<std::pair<double, double>>, 3> thresholds;
  std::array<std::pair<double, double>, 3> lines;
  std::array<photonjet_science::Trace,3> traces;
  for (size_t i = 0; i < levels.size(); ++i) {
    const int level = std::get<0>(levels[i]);
    auto *median = dynamic_cast<TH1 *>(dir->Get(Form("median_%d", level)));
    auto *fit = dynamic_cast<TH1 *>(dir->Get(Form("line_%d", level)));
    lines[i] = {std::nan(""), std::nan("")};
    if (!median) continue;
    traces[i]=photonjet_science::derivedTrace("saved_working_point",{median,fit});
    const Bins medians = binsOf(median);
    for (int k = 0; k < medians.n(); ++k)
      if (std::isfinite(medians.v[k])) thresholds[i].push_back({medians.center(k), medians.v[k]});
    if (fit) lines[i] = {fit->GetBinContent(1), fit->GetBinContent(2)};
  }
  TCanvas cv("cv", "", 1500, 900);
  /* the medians and their lines read the Au+Au QA directly: name its inputs */
  const auto &selection = renderSelection("AuAu", true);
  /* The lines the analysis applies stand in the selection record of the packages; the lines of this figure are
     derived from the simulation at hand. The panel names the derived line an applied one agrees with, or draws the
     applied one beside them when it agrees with none: the figure and the analysis cannot differ unseen.
     Isolation: photon_id.isolation.threshold. Score: the pair of a fill that states its working points explicitly,
     photon_id.optional_model.working_point_choice with mode "explicit" (tight, non_tight_upper). A fill that takes the
     working points of a model release names no line of its own there, and none is drawn. */
  std::vector<std::pair<std::pair<double, double>, std::string>> applied;
  const auto lineOf = [](const PlotJson &record) {
    return std::pair<double, double>{record.at("intercept").get<double>(), record.at("slope").get<double>()};
  };
  if (!selection.is_null() && selection.contains("photon_id")) {
    const auto &id = selection.at("photon_id");
    if (!bdt && id.contains("isolation") && id.at("isolation").contains("threshold") &&
        id.at("isolation").at("threshold").is_object())
      applied.push_back({lineOf(id.at("isolation").at("threshold")), "applied in the analysis"});
    if (bdt && id.contains("optional_model") && id.at("optional_model").is_object() &&
        id.at("optional_model").contains("working_point_choice")) {
      const auto &choice = id.at("optional_model").at("working_point_choice");
      if (choice.is_object() && choice.value("mode", std::string()) == "explicit") {
        if (choice.contains("tight") && choice.at("tight").is_object())
          applied.push_back({lineOf(choice.at("tight")), "applied, tight"});
        if (choice.contains("non_tight_upper") && choice.at("non_tight_upper").is_object())
          applied.push_back({lineOf(choice.at("non_tight_upper")), "applied, non-tight edge"});
      }
    }
  }
  workingPointFitPanel(sliceAxis->GetXmin(), sliceAxis->GetXmax(), thresholds, lines, ytitle,
                       bdt ? "s_{BDT}^{cut} = a + bC" : "E_{T}^{iso,cut} = a + bC", bdt ? 4 : 3,
                       {sample, Form("%g #leq E_{T}^{#gamma} < %g GeV", ptLow, ptHigh)}, simulation, applied, traces);
  savePhysicsQA(cv, fitsPng);
}

/* The BDT-score thresholds of the 70, 80 and 90% working points against photon E_T in the two centrality slices of the
   saved comparison, the most central and the most peripheral, side by side on one y axis, with the median of each
   working point over the E_T bins as a dashed line. The thresholds are the saved quantiles of the signal scores pooled
   over each slice: a 0-10% threshold is never the average of two 5% thresholds.
   Read from model_qa/<name>/application/centrality_comparison: thresholds_<level> (slice against photon E_T),
   median_<level> (per slice) and the record `derivation`, which names the two slices. */
void workingPointCentralityComparison(const std::string &name, const PlotJson &index) {
  const std::string output = "photonOnly/auau_h70_bdt_centrality_comparison.png";
  const std::string folder = "model_qa/" + name + "/application/centrality_comparison";
  auto *dir = gQA->GetDirectory(folder.c_str());
  if (!dir) return missing(output, "missing saved pooled centrality-comparison calibration");
  const auto derivation = savedJson(folder + "/derivation");
  if (derivation.is_null()) throw std::runtime_error("missing centrality-comparison calibration provenance");
  const auto levels = workingPointLevels();
  double low = 1, high = 0, ptLow = 0, ptHigh = 0;  /* the E_T range: set from the saved thresholds */
  for (const auto &level : levels) {
    auto *t = dynamic_cast<TH2 *>(dir->Get(Form("thresholds_%d", std::get<0>(level))));
    auto *m = dynamic_cast<TH1 *>(dir->Get(Form("median_%d", std::get<0>(level))));
    if (!t || t->GetNbinsX() != 2 || !m || m->GetNbinsX() != 2)
      throw std::runtime_error("incomplete pooled centrality-comparison thresholds");
    ptLow = t->GetYaxis()->GetXmin();
    ptHigh = t->GetYaxis()->GetXmax();
    for (int c = 1; c <= 2; ++c)
      for (int e = 1; e <= t->GetNbinsY(); ++e) {
        const double value = t->GetBinContent(c, e);
        if (!std::isfinite(value)) throw std::runtime_error("nonfinite centrality-comparison threshold");
        low = std::min(low, value);
        high = std::max(high, value);
      }
  }
  const double span = std::max(high - low, .10);
  /* Two equal frames on one y axis. Each carries its centrality range at the top left, the first the sPHENIX mark and
     the sample line under it, the second the legend at the top right; the y range opens until every block stands clear
     of the points and of the dashed medians beneath it. */
  TCanvas canvas("wp_centrality_comparison", "", 1800, 720);
  const double w = canvas.GetWw(), h = canvas.GetWh(), left = 124, right = 26, top = 14, bottom = 100, gap = 50;
  const double fw = (w - left - right - gap) / 2, fh = h - top - bottom;
  const double inset = 32, rangePx = 40, markPx = 32, samplePx = 24, legendPx = 24, keyPx = 54, rowPx = 36,
               legendTop = 26;
  const std::string sample = photonJetSampleLine(name);
  std::vector<std::string> ranges, entries;
  for (const auto &slice : derivation.at("centrality_slices"))
    ranges.push_back(centralityRange(slice.at("low").get<double>(), slice.at("high_exclusive").get<double>()));
  if (ranges.size() != 2) throw std::runtime_error("centrality comparison needs two slices");
  double entryWidthPx = 0;
  for (const auto &level : levels) {
    entries.push_back(Form("%d%% signal efficiency", std::get<0>(level)));
    entryWidthPx = std::max(entryWidthPx, textWidth(entries.back(), legendPx) * w);
  }
  const double legendWidth = keyPx + 10 + entryWidthPx, legendDepth = legendTop + levels.size() * rowPx,
               rangeDepth = 30 + rangePx;
  const double headerDepth = rangeDepth + 12 + markPx + 10 + samplePx,
               rangeWidth = textWidth(ranges[1], rangePx, 63) * w;
  const double headerWidthPx = std::max({textWidth(ranges[0], rangePx, 63),
                                         textWidth(sPhenixInternalLabel(), markPx), textWidth(sample, samplePx)}) * w;
  /* the highest point or median of slice c over the stretch [x0, x1] of its frame, in pixels from its left edge */
  const auto under = [&](int c, double x0, double x1) {
    double y = -std::numeric_limits<double>::infinity();
    for (const auto &level : levels) {
      auto *t = dynamic_cast<TH2 *>(dir->Get(Form("thresholds_%d", std::get<0>(level))));
      y = std::max(y, dynamic_cast<TH1 *>(dir->Get(Form("median_%d", std::get<0>(level))))->GetBinContent(c));
      for (int e = 1; e <= t->GetNbinsY(); ++e) {
        const double x = (t->GetYaxis()->GetBinCenter(e) - ptLow) / (ptHigh - ptLow) * fw;
        if (x > x0 - 16 && x < x1 + 16) y = std::max(y, t->GetBinContent(c, e));
      }
    }
    return y;
  };
  const double yLow = low - .10 * span;
  /* a block stays 18 px clear of what is drawn under it */
  const auto opened = [&](double y, double depth) { return raisedTop(yLow, y, depth + 18, fh); };
  const double yHigh = std::max({high + .06 * span, opened(under(1, inset, inset + headerWidthPx), headerDepth),
                                 opened(under(2, inset, inset + rangeWidth), rangeDepth),
                                 opened(under(2, fw - 18 - legendWidth, fw), legendDepth)});
  for (int c = 1; c <= 2; ++c) {
    const double fx = left + (c - 1) * (fw + gap);
    TH1F *axes = pixelPanel(canvas, fx, top, fx + fw, top + fh, ptLow, yLow, ptHigh, yHigh, 28);
    axes->GetXaxis()->SetNdivisions(505);
    axes->GetYaxis()->SetNdivisions(axisDivisions(yHigh - yLow, 6.5));
    if (c == 2) axes->GetYaxis()->SetLabelSize(0);
    for (const auto &[level, color] : levels) {
      const double value = dynamic_cast<TH1 *>(dir->Get(Form("median_%d", level)))->GetBinContent(c);
      auto* median=dynamic_cast<TH1*>(dir->Get(Form("median_%d",level)));
      photonjet_science::annotation(*gPad,"working_point/median/"+std::to_string(c)+"/"+std::to_string(level),
        photonjet_science::scalar(value,"score","saved pooled centrality median"),
        photonjet_science::selectedTrace(median,median->GetXaxis()->GetBinLowEdge(c),median->GetXaxis()->GetBinUpEdge(c)));
      pixelDashes(ptLow, value, ptHigh, value, fw, color);
    }
    for (size_t i = 0; i < levels.size(); ++i) {
      const auto [level, color] = levels[i];
      auto *t = dynamic_cast<TH2 *>(dir->Get(Form("thresholds_%d", level)));
      auto *graph = new TGraph;
      for (int e = 1; e <= t->GetNbinsY(); ++e)
        graph->SetPoint(graph->GetN(), t->GetYaxis()->GetBinCenter(e), t->GetBinContent(c, e));
      photonjet_science::bindNative(graph,"working_point/comparison/"+std::to_string(c)+"/"+std::to_string(i),
        photonjet_science::projectionOf(t,photonjet_science::window(0,t->GetXaxis()->GetBinLowEdge(c),t->GetXaxis()->GetBinUpEdge(c))),
        photonjet_science::covariance("not_computed",{},0,"pooled thresholds have no covariance"));
      style(graph, color, kWorkingPointMarker, 1.8);
      photonjet_science::draw(graph,"P");
      /* the legend: each marker on its dashed line */
      if (c == 2)
        pixelKey(fx + fw - 18 - legendWidth, top + legendTop + (i + .5) * rowPx, keyPx, color, kWorkingPointMarker,
                 1.8, entries[i], legendPx);
    }
    pixelLabel(fx + inset, top + 30 + .74 * rangePx, ranges[c - 1], rangePx, 11, cBlack, true);
    if (c == 1) {
      pixelLabel(fx + inset, top + rangeDepth + 12 + .78 * markPx, sPhenixInternalLabel(), markPx);
      pixelLabel(fx + inset, top + headerDepth - .2 * samplePx, sample, samplePx);
    }
    pixelLabel(fx + fw / 2, h - 36, "Photon E_{T} [GeV]", 38, 22);
  }
  canvas.cd();
  pixelLabel(34, top + fh / 2, "BDT score threshold", 38, 22, cBlack, false, 90);
  saveReleasedWorkingPoints(canvas, output, index, name);
  amendReceipt(output, [&](PlotJson &receipt) {
    receipt["plot_recipe"] = plotRecipe("auau_bdt_threshold_grid", "centrality_extremes", "supplemental", output,
                                        "photonOnly/auau_h70_bdt_pt_centrality_grid.png",
                                        "workingPointCentralityComparison");
    receipt["parameter_panel"] = "centrality_comparison";
    receipt["calibration_comparison"] = derivation;
  });
}

/* The BDT-score thresholds of the 70, 80 and 90% working points of a saved model against centrality, with the straight
   lines s = a + bC that define its cuts. The points are the thresholds of the centrality bins the lines were fitted to;
   points and lines are drawn as saved, nothing is refitted.
   Read from model_qa/<name>/meta/working_points. The figure is drawn only for a record with one cell whose axis is
   the centrality (cent_value); without the record, or for any other, the function returns without a figure and
   without listing one. */
void savedWorkingPointFits(const std::string &name, const PlotJson &index) {
  const auto wp = savedJson("model_qa/" + name + "/meta/working_points");
  if (wp.is_null()) return;
  if (wp.at("cells").size() != 1) return;
  const auto &cell = wp.at("cells").begin().value();
  if (cell.value("axis", "") != "cent_value") return;
  const auto &bounds = cell.at("bounds");
  const double low = bounds.at(0), high = bounds.at(1);
  const auto levels = workingPointLevels();
  std::array<std::vector<std::pair<double, double>>, 3> thresholds;
  std::array<std::pair<double, double>, 3> lines;
  std::array<photonjet_science::Trace,3> traces;
  for (size_t i = 0; i < levels.size(); ++i) {
    const auto &cut = cell.at("cuts").at("WP" + std::to_string(std::get<0>(levels[i])));
    traces[i]=photonjet_science::recordTrace(*gQA,"model_qa/"+name+"/meta/working_points")
      .then("released_working_point",{{"efficiency_percent",std::get<0>(levels[i])},{"cut",cut}});
    lines[i] = {cut.at("intercept").get<double>(), cut.at("slope").get<double>()};
    for (const auto &row : cut.at("bins"))
      thresholds[i].push_back({row.at("center").get<double>(), row.at("threshold").get<double>()});
  }
  const auto &domain = wp.at("calibration_domain");
  const std::string window = Form("%g #leq E_{T}^{#gamma} < %g GeV", domain.at("et_min").get<double>(),
                                  domain.at("et_max_exclusive").get<double>());
  TCanvas canvas("released_wp", "", 1500, 900);
  workingPointFitPanel(low, high, thresholds, lines, "BDT score threshold", "s_{BDT}^{cut} = a + bC", 4,
                       {photonJetSampleLine(name), window}, "", {}, traces);
  saveReleasedWorkingPoints(canvas, "photonOnly/auau_h70_bdt_centrality_fits.png", index, name);
}

/* ---- saved photon-ID models: training diagnostics -------------------------------------------------------------- */

/* A training diagnostic shows the sample and the weighting that its row of kTrainingDiagnostics names: the
   weighting describes the learning and never stands for physics yields. Each of these canvases carries a line with
   the words "Training diagnostic", which save() requires.
   Their text is placed in NDC of the current pad: diagnosticText draws a line whose indices must sit on their symbols,
   at an exactTextPx size; legendRow draws a row of legend entries. */

/* The width of `text` at `px` pixels as it is drawn, in NDC of the current pad. It is measured in the font of scalable
   precision: the box TLatex reports for a font sized in pixels takes an index at the full size of the text. */
double drawnWidth(const std::string &text, double px, bool bold = false) {
  return textWidth(text, px / std::min(padWidthPx(), padHeightPx()), bold ? 62 : 42);
}

/* A line of text of a training-diagnostic figure, `px` high (an exactTextPx size), anchored at (x, y) in NDC of the
   current pad as the TLatex alignment code `align` says, its vertical centre and top being those of a capital letter;
   upwards along the y axis with vertical = true. Returns the last piece drawn, for a colour.
   A word that carries an index starts a piece of its own, so that no index is placed by the measured width of a long
   text before it; the pieces are laid one word space apart. Every piece is placed by its left end and its drawn
   width, also in centred and in right-aligned text, where TLatex would shift it by its width on the canvas. */
TLatex *diagnosticText(double x, double y, const std::string &text, double px, int align = 11, bool bold = false,
                       bool vertical = false) {
  const int horizontal = align / 10;  /* 1: from the anchor on, 2: centred on it, 3: up to it */
  const double w = padWidthPx(), h = padHeightPx(), space = .26 * px;  /* the font's word space */
  struct Piece { std::string text; double width = 0; };  /* width: as drawn, px */
  std::vector<Piece> pieces(1);
  int depth = 0;
  for (std::size_t at = 0, start = 0; at <= text.size(); ++at) {
    if (at < text.size()) depth += text[at] == '{' ? 1 : text[at] == '}' ? -1 : 0;
    if (at < text.size() && (text[at] != ' ' || depth != 0)) continue;
    const std::string word = text.substr(start, at - start);
    if (!pieces.back().text.empty() && word.find_first_of("_^") != std::string::npos) pieces.emplace_back();
    pieces.back().text += (pieces.back().text.empty() ? "" : " ") + word;
    start = at + 1;
  }
  double total = -space;
  for (auto &piece : pieces) {
    piece.width = drawnWidth(piece.text, px, bold) * w;
    total += piece.width + space;
  }
  if (pieces.size() > 1) gVisibleText.push_back(text);
  TLatex t, *last = nullptr;
  t.SetNDC();
  t.SetTextFont(bold ? 63 : 43);
  t.SetTextSize(px);
  t.SetTextAlign(11);
  t.SetTextAngle(vertical ? 90 : 0);
  /* px from the anchor: down to the baseline, and along the text to the left end of a piece */
  const double drop = (align % 10 == 2 ? .36 : align % 10 == 3 ? .75 : 0) * px;
  double along = horizontal == 2 ? -total / 2 : horizontal == 3 ? -total : 0;
  for (const Piece &piece : pieces) {
    last = vertical ? t.DrawLatex(x + drop / w, y + along / h, piece.text.c_str())
                    : t.DrawLatex(x + along / w, y - drop / h, piece.text.c_str());
    along += piece.width + space;
  }
  return last;
}

/* One row of legend entries (object, label, TLegend option) in the current pad, from x (NDC) to the right and centred
   on y. Every entry is as wide as its own text, so the gaps between the entries are equal. Returns the right edge of
   the row; draw = false only measures it. */
double legendRow(double x, double y, const std::vector<std::tuple<TObject *, std::string, const char *>> &entries,
                 double px, bool draw = true) {
  const double w = padWidthPx(), h = padHeightPx(), symbol = 1.7 * px / w, gap = 1.5 * px / w;
  for (const auto &[object, label, option] : entries) {
    if (draw) {
      /* the box holds the label as the saving code measures it, which may be wider than it is drawn */
      const double width = symbol + textWidth(label, px) + 14 / w;
      TLegend *entry = ianLegend(x, y - .65 * px / h, x + width, y + .65 * px / h, symbol / width);
      entry->SetTextSize(px);
      entry->AddEntry(object, label.c_str(), option);
      photonjet_science::draw(entry);
    }
    x += symbol + drawnWidth(label, px) + gap;
  }
  return x - gap;
}

/* Free room for a block of text in a horizontal bar chart of rows one unit high: taken[i] is the x range that the bar
   and the value of row i (counted from the bottom) cover. The block spans [xa, xb] and is `height` rows high. Returns
   the y of its upper edge in the highest (fromTop) or the lowest run of free rows that holds it, or NaN if none does.
   A run that reaches the top row (i = n - 1) may extend `above` of a row above it, one that reaches the bottom row
   (i = 0) `under` of a row below it: the room that the ticks of the frame's edges leave there. */
double barChartRoom(const std::vector<std::pair<double, double>> &taken, double xa, double xb, double height,
                    bool fromTop, double above, double under) {
  const int n = taken.size();
  const auto isFree = [&](int i) { return taken[i].second <= xa || taken[i].first >= xb; };
  std::vector<std::pair<double, double>> runs;  /* from the bottom */
  for (int i = 0; i < n; ++i) {
    if (!isFree(i)) continue;
    int j = i;
    while (j + 1 < n && isFree(j + 1)) ++j;
    runs.push_back({i == 0 ? -under : 1. * i, j == n - 1 ? n + above : j + 1.});
    i = j;
  }
  if (fromTop) std::reverse(runs.begin(), runs.end());
  for (const auto &[low, high] : runs) {
    if (high - low < height) continue;
    const double pad = std::min(.2, (high - low - height) / 2);
    return fromTop ? high - pad : low + pad + height;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

/* The spacing of the minor ticks on a linear axis of that span with `primary` primary divisions, each of five. A frame
   whose limits are multiples of it carries a tick on the frame line and none a hair inside it. */
double minorTickStep(double span, int primary) {
  const auto rounded = [](double raw) {
    const double decade = std::pow(10, std::floor(std::log10(raw))), digit = raw / decade;
    return decade * (digit <= 1 ? 1 : digit <= 2 ? 2 : digit <= 5 ? 5 : 10);
  };
  return rounded(rounded(span / primary) / 5);
}

/* The kinematic reweighting of the BDT training sample, once in cluster E_T and once in cluster eta. Left: the signal
   and the background before reweighting, each as a unit density. Right: the ratio of the two weighted densities with
   a band of +-5% about one and the largest deviation from one. Two IAN frames of equal size side by side, each with
   the sPHENIX mark, the sample and its legend inside. The E_T panels show the E_T intervals of the model's diagnostic
   grid, as the IAN does: whole saved bins are merged.
   Read from model_qa/<name>/training: <et|eta>_<signal|background>_<raw|weighted>, and the training E_T range and
   the diagnostic grid of the model's metadata. The densities, their ratio and the largest |S/B - 1| written in the
   header are computed here from those histograms; they are not saved. */
void trainingReweighting(const std::string &name, const PlotJson &index, const std::string &prefix) {
  const auto &meta = modelMetadata(name);
  const auto range = meta.at("et_ranges_gev").at("training");
  for (const std::string variable : {"et", "eta"}) {
    const bool pp = isPP(meta);
    const std::string system = pp ? "pp" : "auau";
    const std::string relative = prefix + "photonOnly/" + system + "_training_" + variable + "_reweighting.png";
    std::vector<double> grid;  /* the intervals shown; empty: the saved bins */
    if (variable == "et" && meta.contains("grid_et_edges")) grid = meta.at("grid_et_edges").get<std::vector<double>>();
    bool merged = false;
    auto hist = [&](const std::string &klass, const std::string &weight) {
      const std::string path = "model_qa/" + name + "/training/" + variable + "_" + klass + "_" + weight;
      auto *saved = dynamic_cast<TH1 *>(gQA->Get(path.c_str()));
      if (!saved || !(saved->Integral() > 0)) throw std::runtime_error("empty saved training histogram");
      /* a grid merges the saved bins only if it spans their range and every edge of it is a saved bin edge */
      const TAxis *axis = saved->GetXaxis();
      const double tolerance = 1e-6 * axis->GetBinWidth(1);
      bool merge = grid.size() > 2 && std::abs(grid.front() - axis->GetXmin()) < tolerance &&
                   std::abs(grid.back() - axis->GetXmax()) < tolerance;
      for (std::size_t i = 0; merge && i < grid.size(); ++i) {
        const double nearest = axis->GetBinLowEdge(axis->FindFixBin(grid[i] + tolerance));
        merge = std::abs(nearest - grid[i]) < tolerance && (i == 0 || grid[i] > grid[i - 1] + tolerance);
      }
      auto *out = merge ? saved->Rebin(grid.size() - 1, fresh().c_str(), grid.data())
                        : static_cast<TH1 *>(saved->Clone(fresh().c_str()));
      merged = merge;
      out->SetDirectory(nullptr);
      const auto trace=photonjet_science::projectionOf(saved,PlotJson::array()).then("training_unit_density",{{"rebin_edges",merge?PlotJson(grid):PlotJson()},{"integral",out->Integral()},{"divide_width",true}});
      out->Scale(1 / out->Integral(), "width");
      photonjet_science::attachTrace(out,trace);
      photonjet_science::bindNative(out,"training/"+variable+"/"+klass+"/"+weight,trace,
        photonjet_science::covariance("not_computed",{},0,"normalization propagates individual bin errors only"));
      return out;
    };
    auto *signal = hist("signal", "raw"), *background = hist("background", "raw");
    auto *ratio = hist("signal", "weighted"), *den = hist("background", "weighted");
    const auto ratioTrace=photonjet_science::derivedTrace("training_density_ratio",{ratio,den},{{"error_method","ROOT_TH1_Divide"}});
    ratio->Divide(den);
    photonjet_science::bindNative(ratio,"training/"+variable+"/ratio",ratioTrace,
      photonjet_science::covariance("not_computed",{},0,"TH1 ratio diagonal errors only"));
    const char *xlabel = variable == "et" ? "Cluster E_{T} [GeV]" : "Cluster #eta";
    const std::string sample = simulationName(pp);
    const std::string window = Form("%g #leq E_{T} < %g GeV", range.at(0).get<double>(), range.at(1).get<double>());
    const double x0 = signal->GetXaxis()->GetXmin(), x1 = signal->GetXaxis()->GetXmax();
    /* x labels as in the IAN: every 5 GeV, and the eta range in four steps */
    const int divisions = variable == "et" ? 500 + std::max(2, int(std::lround((x1 - x0) / 5))) : -504;
    /* The axes of both frames: the titles centred as in the IAN, every label of an axis written with the same number
       of decimals, as the IAN writes them (0.90, 0.95, 1.00, 1.05, 1.10), and the x labels lowered until the first
       of them is clear of the lowest y label, which stands at the same corner of the frame. */
    const auto axes = [&]() {
      gFrame->GetXaxis()->SetNdivisions(divisions);
      gFrame->GetXaxis()->SetLabelOffset(.02);
      for (TAxis *axis : {gFrame->GetXaxis(), gFrame->GetYaxis()}) axis->SetDecimals(), axis->CenterTitle();
    };
    TCanvas canvas(fresh().c_str(), "", 1560, 650);
    canvas.Divide(2, 1, .001, .001);

    canvas.cd(1);
    style(signal, cRed, 20);
    style(background, cBlue, 20);
    signal->SetLineWidth(3);
    background->SetLineWidth(3);
    {
      std::vector<Extent> drawn;
      extents(drawn, signal);
      extents(drawn, background);
      /* named as the classes of the training sample: the saving code takes an entry "Background" alone only for a
         candidate population of the physics samples */
      const std::vector<std::string> labels = {"Training signal", "Training background"};
      TLegend *legend = ianFrame(x0, x1, 0, drawn, false, xlabel, "Unweighted density", {sample, window}, labels);
      axes();
      /* the densities span the frame: their first and last strokes would lie on its left and right edge */
      photonjet_science::draw(signal,"HIST SAME ][");
      photonjet_science::draw(background,"HIST SAME ][");
      /* the legend a tick length further from the right axis, so that no tick reads as a dash after an entry */
      const double clear = 18 / padWidthPx();
      legend->SetX1(legend->GetX1() - clear);
      legend->SetX2(legend->GetX2() - clear);
      legend->AddEntry(signal, labels[0].c_str(), "l");
      legend->AddEntry(background, labels[1].c_str(), "l");
      photonjet_science::draw(legend);
      gPad->RedrawAxis();
    }

    canvas.cd(2);
    /* the largest |S/B - 1|, and the span of the points with their errors and of the band */
    const double bandLow = .95, bandHigh = 1.05;  /* the reference band: 5% either side of one, as its entry says */
    double maximum = .0, lowest = bandLow, highest = bandHigh;
    for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin)
      if (den->GetBinContent(bin) > 0) {
        maximum = std::max(maximum, std::abs(ratio->GetBinContent(bin) - 1));
        lowest = std::min(lowest, ratio->GetBinContent(bin) - ratio->GetBinError(bin));
        highest = std::max(highest, ratio->GetBinContent(bin) + ratio->GetBinError(bin));
      }
    photonjet_science::annotation(*gPad,"training/"+variable+"/maximum_deviation",photonjet_science::scalar(maximum,"fraction","max absolute S/B minus one over positive denominator bins"),ratioTrace);
    photonjet_science::annotation(*gPad,"training/reference_band",{{"low",bandLow},{"high",bandHigh}},ratioTrace.then("reference_band"));
    const double span = std::max(.12, 1.35 * maximum);
    highest = std::max(highest, 1 + span);  /* the frame reaches as far above one as below it, as in the IAN */
    style(ratio, cPurple, 20, 1.4);
    {
      const std::vector<std::string> lines = {sample, window,
        Form("Training diagnostic: #color[%d]{max |S/B #minus 1| = %.1f%%}", cPurple, 100 * maximum)};
      /* The band's legend entry stands at the bottom right, as in the IAN: the frame starts low enough for the lowest
         point and the band to stay above it, whatever room the header takes at the top. */
      const double w = padWidthPx(), h = padHeightPx();
      const double bottom = gPad->GetBottomMargin(), top = 1 - gPad->GetTopMargin();
      const double headerBottom = header(gPad->GetLeftMargin(), top - (14 + kMarkPx) / h, lines, false, false);
      /* the shares of the frame's height under the header and, relative to it, under the top of the legend entry */
      const double below = (headerBottom - 16 / h - bottom) / (top - bottom);
      const double entry = (24 + kStepPx) / h / (top - bottom) / below;
      const double low = std::min(1 - span, (lowest - entry * highest) / (1 - entry));
      /* the points with their errors, taken half a header line higher than they are: the frame then opens until the
         highest of them is that far under the header */
      std::vector<Extent> drawn;
      extents(drawn, ratio, true);
      const double raised = below / (below - .5 * kStepPx / h / (top - bottom));
      for (auto &e : drawn) e.y = low + (e.y - low) * raised;
      drawn.push_back({x0, x0, 1 + span});
      ianFrame(x0, x1, low, drawn, false, xlabel, "Weighted signal / background", lines, {});
      axes();
      auto *band = new TBox(x0, bandLow, x1, bandHigh);
      band->SetFillColorAlpha(cGreen, .15);
      band->SetLineWidth(0);
      photonjet_science::draw(band,"SAME");
      auto *unity = new TLine(x0, 1, x1, 1);
      unity->SetLineWidth(2);
      photonjet_science::draw(unity);
      photonjet_science::draw(ratio,"E X0 SAME");
      const std::string label = "#pm5% reference band";
      /* as the legend of ianFrame: 14 px from the right axis, a symbol of 46 px and 10 px before the text */
      const double right = 1 - gPad->GetRightMargin() - 14 / w, width = 56 / w + textWidth(label, kLinePx);
      TLegend *legend = ianLegend(right - width, bottom + 10 / h, right, bottom + (10 + kStepPx) / h, 46 / w / width);
      legend->AddEntry(band, label.c_str(), "f");
      photonjet_science::draw(legend);
      gPad->RedrawAxis();
    }
    gProjectionSteps.push_back({{"operation", "display_normalization"},
                                {"scope", "each class to unit density; ratio of saved learner-weighted densities"}});
    if (merged)
      gProjectionSteps.push_back({{"operation", "display_binning"},
                                  {"scope", "whole saved bins merged into the E_T intervals of the diagnostic grid"}});
    saveTrainingDiagnostic(canvas, relative, index, {name}, "training_reweighting");
  }
}

/* The inputs the Au+Au BDT has beyond those of the p+p baseline, by their names in the saved feature order. */
const std::set<std::string> kAuAuOnlyInputs = {"cluster_weta33_cogx", "cluster_wphi33_cogx", "centrality"};

/* The label of a BDT input beside its row of a bar chart, by its name in the saved feature order; an input without a
   label is shown by its name. A thin space sets the multiplication sign of a tower window (E_{3x2}) apart from its two
   digits, and the slash of a ratio apart from the index before it. */
std::string featureLabel(const std::string &name) {
  static const std::string times = "#scale[0.35]{ }#times#scale[0.35]{ }", slash = "#scale[0.5]{ }/";
  static const std::map<std::string, std::string> labels = {
    {"cluster_Et", "Cluster E_{T}"}, {"cluster_Eta", "Cluster #eta"}, {"vertexz", "Vertex z"},
    {"cluster_weta_cogx", "Full-cluster w_{#eta}"}, {"cluster_wphi_cogx", "Full-cluster w_{#phi}"},
    {"cluster_weta33_cogx", "3#times3 w_{#eta}"}, {"cluster_wphi33_cogx", "3#times3 w_{#phi}"},
    {"e11_over_e33", "E_{1" + times + "1}" + slash + "E_{3" + times + "3}"},
    {"e32_over_e35", "E_{3" + times + "2}" + slash + "E_{3" + times + "5}"},
    {"cluster_et1", "Leading tower E_{T}"}, {"cluster_et2", "2nd tower E_{T}"},
    {"cluster_et3", "3rd tower E_{T}"}, {"cluster_et4", "4th tower E_{T}"}, {"centrality", "Centrality"}};
  const auto found = labels.find(name);
  return found == labels.end() ? name : found->second;
}

/* The inputs of the photon-ID BDT as horizontal bars, in the order of their share of the split gain, as in the IAN.
   One figure gives the standardized difference of the signal and background means in the training-weighted sample
   (blue: signal larger; vermilion, cOrange: background larger), the other the share of the total split gain
   (vermilion: the inputs added for Au+Au, kAuAuOnlyInputs). The bars fill the frame; the sPHENIX mark, the sample and
   the legend stand inside it where no bar or value reaches.
   Read from the table `features` of the model (one row per input, in the saved feature order: the two values with
   their validity flags) and from the feature order and the training E_T range of its metadata. */
void featureDiagnostics(const std::string &name, const PlotJson &index, const std::string &prefix) {
  const auto &meta = modelMetadata(name);
  const auto &features = meta.at("feature_order");
  const auto rows = modelTable(name, "features");
  if (rows.size() != features.size()) throw std::runtime_error("saved feature diagnostic size differs");
  std::vector<std::size_t> order;
  for (std::size_t i = 0; i < rows.size(); ++i) order.push_back(i);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return rows[a].at("total_gain_fraction").get<double>() > rows[b].at("total_gain_fraction").get<double>();
  });
  const bool pp = isPP(meta);
  const std::string system = pp ? "pp" : "auau";
  const auto range = meta.at("et_ranges_gev").at("training");
  const int n = rows.size();
  const double textPx = exactTextPx(34), titlePx = exactTextPx(39), markPx = exactTextPx(45), pitch = 42;
  for (bool gain : {false, true}) {
    /* The separation is that of the sample the learner sees: in the unweighted spectra E_T and eta would seem to
       discriminate through the sampling of the two classes alone. */
    const char *validKey = gain ? "gain_valid" : "weighted_separation_valid";
    const char *valueKey = gain ? "total_gain_fraction" : "weighted_standardized_mean_difference";
    const std::string stem = gain ? "feature_split_gain" : "single_feature_separation";
    const std::string relative = prefix + "photonOnly/" + system + "_bdt_" + stem + ".png";
    double maximum = gain ? 1 : 1.1;
    for (const auto &row : rows)
      if (row.at(validKey).get<double>())
        maximum = std::max(maximum, std::abs(row.at(valueKey).get<double>() * (gain ? 100 : 1)) * (gain ? 1.12 : 1.2));
    const double step = minorTickStep((gain ? 1 : 2) * maximum, 10);  /* the axis ends on a tick */
    maximum = step * std::ceil(maximum / step - 1e-9);
    const double xa = gain ? 0 : -maximum, xb = maximum;

    TCanvas canvas(fresh().c_str(), "", 1800, 1220);
    for(std::size_t i=0;i<rows.size();++i) {
      const bool valid=rows[i].at(validKey).get<double>();
      const PlotJson measurement=valid?PlotJson{{"state","measured"},{"feature",features.at(i)},{"value",rows[i].at(valueKey)},{"scale",gain?100:1}}:
        PlotJson{{"state","missing"},{"feature",features.at(i)},{"reason","producer marked feature metric invalid"}};
      photonjet_science::annotation(canvas,"feature/"+std::to_string(i)+"/"+valueKey,measurement,
        photonjet_science::recordTrace(*gQA,"model_qa/"+name+"/features",valueKey).then("feature_metric",{{"feature_slot",i},{"metric",valueKey}}));
    }
    const double w = padWidthPx(), h = padHeightPx();
    double labelWidth = 0;
    for (const auto &feature : features)
      labelWidth = std::max(labelWidth, drawnWidth(featureLabel(feature.get<std::string>()), textPx));
    canvas.SetMargin(labelWidth + 44 / w, 30 / w, 130 / h, 22 / h);
    canvas.SetTicks(1, 0);
    const double left = canvas.GetLeftMargin(), bottom = canvas.GetBottomMargin();
    const double plotWidth = 1 - left - canvas.GetRightMargin(), plotHeight = 1 - bottom - canvas.GetTopMargin();
    const double perPx = (xb - xa) / (w * plotWidth);  /* x-axis units per pixel */

    /* each row: its bar, its value and the x range the two cover */
    struct Row { std::string label, value; double end = 0, at = 0; int align = 12, color = cBlue; bool bar = false; };
    std::vector<Row> drawn(n);
    std::vector<std::pair<double, double>> taken(n);
    bool base = false, added = false;
    for (int rank = 0; rank < n; ++rank) {
      const auto i = order[rank];
      const std::string feature = features.at(i);
      Row &row = drawn[rank];
      row.label = featureLabel(feature);
      if (!rows[i].at(validKey).get<double>()) row.value = "unavailable", row.at = 14 * perPx;
      else {
        const double value = rows[i].at(valueKey).get<double>() * (gain ? 100 : 1);
        const bool extra = kAuAuOnlyInputs.count(feature) > 0;
        (extra ? added : base) = true;
        row.bar = true, row.end = value;
        row.color = gain ? (extra ? cOrange : cBlue) : (value < 0 ? cOrange : cBlue);
        const bool zero = !gain && std::abs(value) < .005;  /* reads 0.00: written as in the IAN */
        if (gain) row.value = Form("%.1f", value);
        else row.value = zero ? "#approx0" : Form("%s%.2f", value < 0 ? "#minus" : "+", std::abs(value));
        if (value < 0 && !zero) row.align = 32;
        row.at = (zero ? 0 : value) + (row.align == 32 ? -14 : 14) * perPx;
      }
      const double reach = row.at + (row.align == 32 ? -1 : 1) * (textWidth(row.value, textPx) * w + 26) * perPx;
      taken[n - 1 - rank] = {std::min({reach, row.end, -12 * perPx}), std::max({reach, row.end, 12 * perPx})};
    }

    /* The header: beside the first rows where they leave room, for which the frame may open by 0.6 of a row more
       above them; else in the highest run of free rows, at the left or at the right; in a band above the bars if none
       holds it. `ticks` is the room the ticks of the upper and of the lower edge take, in rows. */
    const std::vector<std::string> lines = {simulationName(pp),
      Form("%g #leq E_{T}^{cluster} < %g GeV", range.at(0).get<double>(), range.at(1).get<double>()),
      gain ? "Training diagnostic: split gain" : "Training diagnostic: input separation"};
    /* rowPx: the height of a row for a frame that holds the n rows, half a row beyond each outer row and the 0.6
       of a row it may open for the header */
    const double rowPx = h * plotHeight / (n + 1.6), inset = 28 * perPx, tickPx = 20, ticks = (tickPx + 8) / rowPx;
    const double edge = .5 - ticks;  /* the room between an outer row and the ticks */
    double blockWidth = textWidth(sPhenixInternalLabel(), markPx);
    for (const auto &text : lines) blockWidth = std::max(blockWidth, textWidth(text, textPx));
    blockWidth *= w * perPx;
    const double blockRows = (markPx + lines.size() * pitch + 6) / rowPx;
    bool headerLeft = false;
    double headerTop = std::numeric_limits<double>::quiet_NaN();
    for (double above : {edge, edge + .6}) {
      const double atRight = barChartRoom(taken, xb - inset - blockWidth, xb - inset, blockRows, true, above, edge);
      const double atLeft = barChartRoom(taken, xa + inset, xa + inset + blockWidth, blockRows, true, above, edge);
      headerLeft = std::isfinite(atLeft) && !(atRight > atLeft);
      headerTop = headerLeft ? atLeft : atRight;
      if (headerTop >= n) break;
    }
    if (!std::isfinite(headerTop)) headerLeft = true, headerTop = n + .1 + blockRows;
    double top = std::max(n + .5, headerTop + ticks);
    const double headerX = headerLeft ? xa + inset : xb - inset - blockWidth;
    for (int i = 0; i < n; ++i)
      if (i + 1 > headerTop - blockRows && i < headerTop)
        taken[i] = {std::min(taken[i].first, headerX), std::max(taken[i].second, headerX + blockWidth)};

    /* the legend: at the bottom right as in the IAN, or where the bars leave room */
    TH1D blue(fresh().c_str(), "", 1, 0, 1), orange(fresh().c_str(), "", 1, 0, 1);
    blue.SetFillColor(cBlue), blue.SetLineColor(cBlue), orange.SetFillColor(cOrange), orange.SetLineColor(cOrange);
    const std::vector<std::tuple<TObject *, std::string, const char *>> entries = {
      {&blue, gain ? "Baseline inputs" : "Signal mean larger", "f"},
      {&orange, gain ? "Additional Au+Au inputs" : "Background mean larger", "f"}};
    const bool showLegend = !gain || (base && added);  /* the split gain of a p+p model has one class of inputs only */
    const double symbol = 1.7 * textPx / w;
    double legendWidth = 0;  /* two stacked entries, or (split gain) one row as in the IAN */
    for (const auto &entry : entries)
      legendWidth = std::max(legendWidth, symbol + textWidth(std::get<1>(entry), textPx) + 14 / w);
    if (gain) legendWidth = legendRow(0, 0, entries, textPx, false);
    legendWidth *= w * perPx;
    const double legendRows = (gain ? 1.3 * textPx : entries.size() * pitch) / rowPx;
    bool legendLeft = false;
    double legendTop = barChartRoom(taken, xb - inset - legendWidth, xb - inset, legendRows, false, edge, edge);
    if (!std::isfinite(legendTop)) {
      legendLeft = true;
      legendTop = barChartRoom(taken, xa + inset, xa + inset + legendWidth, legendRows, false, edge, edge);
    }
    if (showLegend && !std::isfinite(legendTop)) {
      legendLeft = false, legendTop = n + .1 + legendRows;
      top = std::max(top, legendTop + ticks);
    }
    const double legendX = legendLeft ? xa + inset : xb - inset - legendWidth;

    TH1F *axes = photonjet_science::drawFrame(*gPad,xa, -.5, xb, top);
    axes->GetXaxis()->SetLabelFont(43);
    axes->GetXaxis()->SetLabelSize(textPx);
    axes->GetXaxis()->SetLabelOffset(.008);
    /* ticks tickPx long: a tick length is a share of the pad's height, scaled by the share of its width the axis
       spans */
    axes->GetXaxis()->SetTickLength(tickPx / (h * plotWidth));
    axes->GetYaxis()->SetLabelSize(0);
    axes->GetYaxis()->SetTickLength(0);
    const auto ndcX = [&](double x) { return left + (x - xa) / (xb - xa) * plotWidth; };
    const auto ndcY = [&](double y) { return bottom + (y + .5) / (top + .5) * plotHeight; };
    for (int rank = 0; rank < n; ++rank) {
      const Row &row = drawn[rank];
      const double y = n - rank - .5;
      if (row.bar) {
        auto *bar = new TBox(std::min(0., row.end), y - .34, std::max(0., row.end), y + .34);
        bar->SetFillColor(row.color);
        bar->SetLineWidth(0);
        photonjet_science::draw(bar);
      }
      diagnosticText(ndcX(row.at), ndcY(y), row.value, textPx, row.align);
      diagnosticText(left - 16 / w, ndcY(y), row.label, textPx, 32);
    }
    if (!gain) {
      auto *zero = new TLine(0, -.5, 0, top);
      zero->SetLineWidth(2);
      photonjet_science::draw(zero);
    }
    const char *xtitle = gain ? "Fraction of total split gain [%]"
                              : "Standardized signal #minus background mean difference";
    diagnosticText(ndcX((xa + xb) / 2), 20 / h, xtitle, titlePx, 21);
    double y = ndcY(headerTop) - .8 * markPx / h;
    diagnosticText(ndcX(headerX), y, sPhenixInternalLabel(), markPx);
    for (const auto &text : lines) diagnosticText(ndcX(headerX), y -= pitch / h, text, textPx);
    if (showLegend && gain) legendRow(ndcX(legendX), ndcY(legendTop - legendRows / 2), entries, textPx);
    else if (showLegend) {
      const double x1 = ndcX(legendX), x2 = ndcX(legendX + legendWidth);
      TLegend *stacked = ianLegend(x1, ndcY(legendTop - legendRows), x2, ndcY(legendTop), symbol / (x2 - x1));
      stacked->SetTextSize(textPx);
      for (const auto &[object, label, option] : entries) stacked->AddEntry(object, label.c_str(), option);
      photonjet_science::draw(stacked);
    }
    gPad->RedrawAxis();
    const char *diagnostic = gain ? "split_gain" : "feature_separation";
    saveTrainingDiagnostic(canvas, relative, index, {name}, diagnostic);
  }
}

/* Classifiers trained in restricted E_T and centrality regions against the global BDT, as in the IAN: one row of
   grouped bars per centrality class, 1000 x (routed AUC - global AUC) on the common held-out sample in every cluster
   E_T interval, one colour per routing scheme. The three rows share one y range and one y title. The sPHENIX mark, the
   class labels and the legend stand where the bars leave the frames free: above the bars, or below them when the
   differences are mostly negative.
   Read from the table `holdout_et_c3` (one row per E_T interval and centrality class: the AUC and its validity flag)
   of the reference model and of every routed model of comparisons.routing. The differences are computed here. */
void routingDiagnostics(const PlotJson &index) {
  const std::string output = "photonOnly/auau_bdt_routing_auc_by_et_centrality.png";
  const auto group = index.at("comparisons").at("routing");
  const std::string reference = group.at("reference");
  if (group.at("models").empty()) return missing(output, "no matched routed model releases");
  const auto baseline = modelTable(reference, "holdout_et_c3");
  std::vector<std::string> names = {reference};
  std::vector<std::vector<PlotJson>> tables;
  double low = 0, high = 0;
  for (const auto &item : group.at("models")) {
    const std::string name = item;
    const auto rows = modelTable(name, "holdout_et_c3");
    if (rows.size() != baseline.size()) throw std::runtime_error("routing QA grids differ");
    for (std::size_t i = 0; i < rows.size(); ++i) {
      for (const char *key : {"et_low", "et_high", "centrality_low", "centrality_high"})
        if (rows[i].at(key) != baseline[i].at(key)) throw std::runtime_error("routing QA boundaries differ");
      if (rows[i].at("valid").get<double>() && baseline[i].at("valid").get<double>()) {
        const double delta = 1000 * (rows[i].at("auc").get<double>() - baseline[i].at("auc").get<double>());
        if (!std::isfinite(delta)) throw std::runtime_error("nonfinite routing AUC difference");
        low = std::min(low, delta);
        high = std::max(high, delta);
      }
    }
    names.push_back(name);
    tables.push_back(rows);
  }
  /* the legend label and the colour of each routing scheme, by the saved name of its model (the second purple is a
     lighter one of this figure); a model that is not listed is shown by its name in the next series colour */
  const std::map<std::string, std::pair<std::string, int>> styles = {
    {"centCoarse_v1", {"C3", cPurple}}, {"centFine_v1", {"C7", TColor::GetColor("#A65FD4")}},
    {"ETcentCoarse_v1", {"E_{T} #times C3", cOrange}}, {"ETcentFine_v1", {"E_{T} #times C7", cGreen}}};
  /* the centrality classes of the three rows, in percent: typed here, and each must be the centrality_low and
     centrality_high of cells of the table */
  const std::array<std::pair<double, double>, 3> classes = {{{0, 20}, {20, 50}, {50, 80}}};

  /* the bars of each row on an x axis of one unit per E_T interval; a cell without a valid AUC has no bar */
  struct Bar { double left, right, value; int series; bool valid; };
  std::array<std::vector<std::size_t>, 3> bins;
  std::array<std::vector<Bar>, 3> bars;
  std::vector<bool> shown(tables.size(), false);
  for (int c = 0; c < 3; ++c) {
    for (std::size_t i = 0; i < baseline.size(); ++i) {
      const auto &cell = baseline[i];
      if (cell.at("centrality_low") == classes[c].first && cell.at("centrality_high") == classes[c].second)
        bins[c].push_back(i);
    }
    if (bins[c].empty()) throw std::runtime_error("empty routing centrality row");
    for (std::size_t m = 0; m < tables.size(); ++m)
      for (std::size_t b = 0; b < bins[c].size(); ++b) {
        const auto i = bins[c][b];
        const auto &row = tables[m][i];
        const double left = b + .12 + m * .76 / tables.size();
        const bool valid = row.at("valid").get<double>() && baseline[i].at("valid").get<double>();
        const double delta = valid ? 1000 * (row.at("auc").get<double>() - baseline[i].at("auc").get<double>()) : 0;
        bars[c].push_back({left, left + .70 / tables.size(), delta, int(m), valid});
        if (valid) shown[m] = true;
      }
  }
  std::vector<std::tuple<TObject *, std::string, const char *>> entries;  /* the legend: a box per routing scheme */
  std::vector<int> colors;
  for (std::size_t m = 0; m < tables.size(); ++m) {
    const auto look = styles.find(names[m + 1]);
    colors.push_back(look == styles.end() ? seriesColor(m + 1) : look->second.second);
    auto *sample = new TBox(0, 0, 1, 1);
    sample->SetFillColor(colors.back());
    sample->SetLineColor(colors.back());
    if (shown[m]) entries.push_back({sample, look == styles.end() ? names[m + 1] : look->second.first, "f"});
  }

  /* Three frames of equal size, laid out in pixels: the y title at the left, the E_T intervals and the x title under
     the last frame. */
  TCanvas canvas(fresh().c_str(), "", 2050, 1430);
  const double w = padWidthPx(), h = padHeightPx();
  const double tickLabelPx = exactTextPx(36), titlePx = exactTextPx(45), markPx = exactTextPx(52);
  const double linePx = exactTextPx(40), classPx = exactTextPx(45);
  const double leftPx = 190, rightPx = 28, topPx = 16, bottomPx = 150, gapPx = 34, inset = 30, edge = 22;
  const double framePx = (h - topPx - bottomPx - 2 * gapPx) / 3, widthPx = w - leftPx - rightPx;
  const std::vector<std::string> lines = {simulationName(false), "Training diagnostic: restricted-region BDTs"};
  const auto rowLabel = [&](int c) {
    return std::string(Form("Au+Au %g#scale[0.4]{ }#minus#scale[0.4]{ }%g%%", classes[c].first, classes[c].second));
  };
  const auto stacked = [](const std::vector<double> &sizes) {  /* baselines of stacked lines, px below their top */
    std::vector<double> y;
    double cursor = 0;
    for (double size : sizes) y.push_back(cursor + .8 * size), cursor += 1.3 * size;
    return y;
  };
  const std::vector<double> headerSizes = {markPx, linePx, linePx, classPx}, classSizes = {classPx};
  double blockWidth = std::max(textWidth(sPhenixInternalLabel(), markPx), textWidth(rowLabel(0), classPx, 63));
  for (const auto &text : lines) blockWidth = std::max(blockWidth, textWidth(text, linePx));
  blockWidth *= w;
  const double legendX = inset + blockWidth + 70, legendWidth = legendRow(0, 0, entries, linePx, false) * w;

  /* The y range, shared by the rows: the bars with a margin, opened at the top or at the bottom, whichever takes less,
     until the labels of every row are clear of the bars they stand over. A block: the header and the legend of the
     first row, the class label of the others; px from the frame's left edge, and px high with its margins. */
  struct Block { int row; double xa, xb, height; };
  std::vector<Block> blocks = {{0, 0, inset + blockWidth + 30, edge + stacked(headerSizes).back() + .25 * classPx + 16},
                               {0, legendX - 30, legendX + legendWidth + 30, edge + 1.3 * classPx + 16}};
  for (int c = 1; c < 3; ++c)
    blocks.push_back({c, 0, inset + textWidth(rowLabel(c), classPx, 63) * w + 30, edge + 1.05 * classPx + 16});
  const double span = std::max(high - low, 3.0), floorY = low - .07 * span, ceilingY = high + .07 * span;
  double lower = floorY, upper = ceilingY;
  for (const auto &block : blocks) {
    const double fraction = block.height / framePx, perUnit = widthPx / bins[block.row].size();
    if (!(fraction < .8)) throw std::runtime_error("routing labels do not fit the canvas");
    double deepest = 0, tallest = 0;  /* the bars this block stands over */
    for (const auto &bar : bars[block.row])
      if (bar.valid && bar.right * perUnit > block.xa && bar.left * perUnit < block.xb)
        deepest = std::min(deepest, bar.value), tallest = std::max(tallest, bar.value);
    lower = std::min(lower, (deepest - fraction * ceilingY) / (1 - fraction));
    upper = std::max(upper, (tallest - fraction * floorY) / (1 - fraction));
  }
  const bool below = ceilingY - lower < upper - floorY;
  const double loose = below ? lower : floorY, wide = below ? ceilingY : upper;
  /* both limits on a tick of the y axis */
  const double step = minorTickStep(wide - loose, 5);
  const double minimum = step * std::floor(loose / step + 1e-9), maximum = step * std::ceil(wide / step - 1e-9);

  for (int c = 0; c < 3; ++c) {
    const double above = c == 0 ? topPx : gapPx / 2, beneath = c == 2 ? bottomPx : gapPx / 2;
    const double padTop = topPx + c * (framePx + gapPx) - above, padPx = above + framePx + beneath;
    canvas.cd();
    auto *pad = new TPad(fresh().c_str(), "", 0, 1 - (padTop + padPx) / h, 1, 1 - padTop / h);
    pad->SetMargin(leftPx / w, rightPx / w, beneath / padPx, above / padPx);
    pad->SetTicks(1, 1);
    photonjet_science::draw(pad);
    pad->cd();
    for(std::size_t m=0;m<tables.size();++m)for(const auto i:bins[c]) {
      const auto& a=tables[m][i];const auto& b=baseline[i];
      const bool valid=a.at("valid").get<double>() && b.at("valid").get<double>();
      const PlotJson value=valid?PlotJson{{"state","measured"},{"delta_auc_times1000",1000*(a.at("auc").get<double>()-b.at("auc").get<double>())},{"et",{a.at("et_low"),a.at("et_high")}},{"centrality",{a.at("centrality_low"),a.at("centrality_high")}}}:
        PlotJson{{"state","missing"},{"reason","invalid AUC cell"}};
      photonjet_science::annotation(*gPad,"routing/"+std::to_string(m)+"/"+std::to_string(c)+"/"+std::to_string(i),value,
        photonjet_science::Trace::combine("routed_minus_reference_auc",{photonjet_science::recordTrace(*gQA,"model_qa/"+names[m+1]+"/holdout_et_c3"),photonjet_science::recordTrace(*gQA,"model_qa/"+reference+"/holdout_et_c3")},{{"cell",i}}));
    }
    const int groups = bins[c].size();
    TH1F *axes = photonjet_science::drawFrame(*pad,0, minimum, groups, maximum);
    axes->GetXaxis()->SetNdivisions(groups, false);
    axes->GetXaxis()->SetLabelSize(0);
    /* the two tick lengths are ROOT's own fractions of the pad, not pixels; the x ticks only part the groups of
       bars */
    axes->GetXaxis()->SetTickLength(.03);
    axes->GetYaxis()->SetLabelFont(43);
    axes->GetYaxis()->SetLabelSize(tickLabelPx);
    axes->GetYaxis()->SetNdivisions(505);
    axes->GetYaxis()->SetTickLength(.008);
    auto *zero = new TLine(0, 0, groups, 0);
    zero->SetLineWidth(2);
    photonjet_science::draw(zero);
    const auto atX = [&](double px) { return (leftPx + px) / w; };         /* px right of the frame's left edge */
    const auto atY = [&](double px) { return 1 - (above + px) / padPx; };  /* px below the frame's upper edge */
    const double zeroPx = maximum / (maximum - minimum) * framePx, groupPx = widthPx / groups;
    for (const auto &bar : bars[c]) {
      if (!bar.valid) {
        const double centre = (bar.left + bar.right) / 2 * groupPx;
        diagnosticText(atX(centre), atY(zeroPx - 8), "NA", 26, 21)->SetTextColor(kGray + 2);
        continue;
      }
      auto *box = new TBox(bar.left, std::min(0., bar.value), bar.right, std::max(0., bar.value));
      box->SetFillColor(colors[bar.series]);
      box->SetLineWidth(0);
      photonjet_science::draw(box);
    }
    /* the labels, at the top or at the bottom of the frame */
    const auto &sizes = c == 0 ? headerSizes : classSizes;
    const auto baselines = stacked(sizes);
    const double shift = below ? framePx - edge - baselines.back() - .25 * classPx : edge;
    if (c == 0) {
      diagnosticText(atX(inset), atY(shift + baselines[0]), sPhenixInternalLabel(), markPx);
      for (std::size_t i = 0; i < lines.size(); ++i)
        diagnosticText(atX(inset), atY(shift + baselines[i + 1]), lines[i], linePx);
      /* the legend on the line of the class label (below the bars) or of the mark (above them) */
      const double centre = below ? baselines.back() - .35 * classPx : baselines[0] - .35 * markPx;
      legendRow(atX(legendX), atY(shift + centre), entries, linePx);
    }
    diagnosticText(atX(inset), atY(shift + baselines.back()), rowLabel(c), classPx, 11, true);
    if (c == 2) {
      for (int b = 0; b < groups; ++b) {
        const auto &cell = baseline[bins[c][b]];
        const double from = cell.at("et_low").get<double>(), to = cell.at("et_high").get<double>();
        /* a thin space on either side of the dash keeps it off the digits */
        const std::string interval = Form("%g#scale[0.4]{ }#minus#scale[0.4]{ }%g", from, to);
        diagnosticText(atX((b + .5) * groupPx), atY(framePx + 14), interval, tickLabelPx, 23);
      }
      diagnosticText(atX(widthPx / 2), 24 / padPx, "Cluster E_{T} interval [GeV]", titlePx, 21);
    }
    pad->RedrawAxis();
  }
  canvas.cd();
  const char *ytitle = "#DeltaAUC relative to global BDT [#times10^{#minus3}]";
  diagnosticText(44 / w, 1 - (topPx + 1.5 * framePx + gapPx) / h, ytitle, titlePx, 22, false, true);
  saveTrainingDiagnostic(canvas, output, index, names, "model_routing");
  amendReceipt(output, [&](PlotJson &receipt) {
    receipt["comparison_statistic"] = {{"reference", reference},
                                       {"formula", "1000 * (routed held-out AUC - reference held-out AUC)"},
                                       {"invalid_cells", "NA; never zero"}};
  });
}

/* The capacity scan of the Au+Au photon-ID BDT against the maximum tree depth, in the three panels of the IAN: (a) the
   AUC and (b) the log loss of the training and the held-out sample, (c) the differences between the two samples as
   bars. A dashed line marks the depth of the nominal classifier (nominalDepth, typed in the function). The panel names
   stand in a thin line above the frames; the sPHENIX mark, the sample and the legends inside them.
   Read for every model of comparisons.depth: the table `performance` (per split: the AUC and the log loss with a
   validity flag) and max_depth of the metadata. The differences of panel (c) are computed here. */
void depthDiagnostics(const PlotJson &index) {
  const std::string output = "photonOnly/auau_bdt_tree_depth_train_holdout_diagnostic.png";
  const auto group = index.at("comparisons").at("depth");
  const auto names = group.at("models").get<std::vector<std::string>>();
  if (names.size() < 2) return missing(output, "need at least two matched depths for a capacity scan");
  /* per depth: the AUC of the training sample (split 0) and of the held-out sample (split 1), then the log loss of
     the two in the same order */
  std::vector<std::pair<double, std::array<double, 4>>> scan;
  for (const auto &name : names) {
    const auto rows = modelTable(name, "performance");
    const double depth = modelMetadata(name).at("max_depth");
    std::array<double, 4> values{};
    std::set<int> splits;
    for (const auto &row : rows) {
      const int split = row.at("split").get<double>();
      if (!row.at("valid").get<double>() || (split != 0 && split != 1) || !splits.insert(split).second)
        return missing(output, "undefined or duplicate train/held-out statistics");
      values[split] = row.at("auc");
      values[split + 2] = row.at("logloss");
      if (!std::isfinite(values[split]) || !std::isfinite(values[split + 2]))
        throw std::runtime_error("nonfinite capacity scan metric");
    }
    if (splits.size() != 2) return missing(output, "incomplete train/held-out statistics");
    scan.push_back({depth, values});
  }
  std::vector<photonjet_science::Trace> scanInputs;
  for(const auto& name:names)scanInputs.push_back(photonjet_science::recordTrace(*gQA,"model_qa/"+name+"/performance"));
  const auto scanTrace=photonjet_science::Trace::combine("tree_depth_comparison",scanInputs);
  std::sort(scan.begin(), scan.end());
  for (std::size_t i = 1; i < scan.size(); ++i)
    if (scan[i].first == scan[i - 1].first) throw std::runtime_error("duplicate capacity scan depth");

  const auto &meta = modelMetadata(names.front());
  const auto range = meta.at("et_ranges_gev").at("training");
  const std::vector<std::string> lines = {
    simulationName(isPP(meta)),
    Form("%g #leq E_{T}^{#gamma} < %g GeV", range.at(0).get<double>(), range.at(1).get<double>()),
    "Training diagnostic: tree-depth scan"};
  const char *titles[] = {"(a) Classifier discrimination", "(b) Predictive log loss", "(c) Generalization gap"};
  const char *ytitles[] = {"AUC", "Log loss", "Train-holdout gap"};
  const double tickLabelPx = exactTextPx(30), titlePx = exactTextPx(34), markPx = exactTextPx(40),
               linePx = exactTextPx(30);
  const double pitch = 40, x0 = scan.front().first - .45, x1 = scan.back().first + .45, nominalDepth = 4;
  const double tickLength = 16, inset = 34;  /* px: the labels start a clear gap behind the ticks of the y axis */
  /* the samples of the legends */
  auto *train = new TGraph, *holdout = new TGraph;
  style(train, cOrange, 20, 1.6);
  style(holdout, cBlue, 21, 1.6);
  train->SetLineWidth(3);
  holdout->SetLineWidth(3);
  auto *aucGap = new TBox(0, 0, 1, 1), *lossGap = new TBox(0, 0, 1, 1);
  aucGap->SetFillColor(cBlue), aucGap->SetLineColor(cBlue), lossGap->SetFillColor(cRed), lossGap->SetLineColor(cRed);
  auto *reference = new TLine;  /* the look of the dashed line at the nominal depth, for the legend */
  reference->SetLineStyle(kDashed);
  reference->SetLineColor(kGray + 2);
  reference->SetLineWidth(2);
  const std::vector<std::tuple<TObject *, std::string, const char *>> curves = {{train, "Training", "lp"},
                                                                                {holdout, "Holdout", "lp"}};
  const std::vector<std::tuple<TObject *, std::string, const char *>> gaps = {
    {aucGap, "AUC_{train} #minus AUC_{holdout}", "f"},
    {lossGap, "Loss_{holdout} #minus Loss_{train}", "f"},
    {reference, Form("Depth-%g reference", nominalDepth), "l"}};

  TCanvas canvas(fresh().c_str(), "", 1850, 1300);
  for (int panel = 0; panel < 3; ++panel) {
    canvas.cd();  /* (a) and (b) side by side above (c) */
    auto *pad = new TPad(fresh().c_str(), "", panel == 1 ? .5 : 0, panel < 2 ? .47 : 0, panel == 0 ? .5 : 1,
                         panel < 2 ? 1 : .47);
    photonjet_science::draw(pad);
    pad->cd();
    const double w = padWidthPx(), h = padHeightPx();
    const double leftPx = 150, rightPx = 24, topPx = 64, bottomPx = panel < 2 ? 72 : 128;
    const double frameW = w - leftPx - rightPx, frameH = h - topPx - bottomPx;
    pad->SetMargin(leftPx / w, rightPx / w, bottomPx / h, topPx / h);
    pad->SetTicks(1, 1);

    std::vector<Extent> drawn;  /* the points, or the ends of the bars */
    double low = 1e300, high = -1e300;
    for (const auto &[depth, v] : scan) {
      const std::array<double, 2> shown = panel < 2 ? std::array<double, 2>{v[2 * panel], v[2 * panel + 1]}
                                                    : std::array<double, 2>{v[0] - v[1], v[3] - v[2]};
      for (double value : shown) {
        low = std::min(low, value);
        high = std::max(high, value);
        drawn.push_back({depth - .4, depth + .4, value});
      }
    }
    const double span = std::max(high - low, panel == 0 ? .015 : .005);
    const double minimum = panel < 2 ? std::max(0., low - .08 * span) : std::min(0., low) - .10 * span;

    /* the labels at the top of the frame: the header (a, left), the two curves (b, right), the legend row (c, left);
       the frame opens until the points under them are clear */
    double blockWidth = 0, blockHeight = 0;  /* px, with the margins around the text */
    if (panel == 0) {
      blockWidth = textWidth(sPhenixInternalLabel(), markPx);
      for (const auto &text : lines) blockWidth = std::max(blockWidth, textWidth(text, linePx));
      blockWidth = inset + blockWidth * w + 16;
      blockHeight = 16 + .8 * markPx + lines.size() * pitch + .25 * linePx + 12;
    } else if (panel == 1) {
      for (const auto &curve : curves) blockWidth = std::max(blockWidth, textWidth(std::get<1>(curve), linePx));
      blockWidth = 14 + 1.9 * linePx + blockWidth * w + 14 + 16, blockHeight = 12 + curves.size() * pitch + 12;
    } else blockWidth = inset + legendRow(0, 0, gaps, linePx, false) * w + 16, blockHeight = 18 + 1.3 * linePx + 12;
    const double perPx = (x1 - x0) / frameW;
    const double xa = panel == 1 ? x1 - blockWidth * perPx : x0, xb = panel == 1 ? x1 : x0 + blockWidth * perPx;
    /* the upper edge: a margin above the highest value, and every value under the labels below their lower edge */
    const double maximum = std::max(high + .06 * span,
                                    topBelow(drawn, minimum, false, xa, xb, 1 - blockHeight / frameH));

    TH1F *axes = photonjet_science::drawFrame(*pad,x0, minimum, x1, maximum);
    for (TAxis *axis : {axes->GetXaxis(), axes->GetYaxis()}) axis->SetLabelFont(43), axis->SetLabelSize(tickLabelPx);
    axes->GetXaxis()->SetNdivisions(10);
    axes->GetYaxis()->SetNdivisions(505);
    axes->GetYaxis()->SetDecimals();
    /* ticks of tickLength px: a tick length is a share of the pad's other side, scaled by the share the axis spans */
    axes->GetXaxis()->SetTickLength(tickLength * w / (h * frameW));
    axes->GetYaxis()->SetTickLength(tickLength * h / (w * frameH));
    diagnosticText(36 / w, (bottomPx + frameH / 2) / h, ytitles[panel], titlePx, 22, false, true);
    if (panel == 2) diagnosticText((leftPx + frameW / 2) / w, 22 / h, "XGBoost maximum tree depth", titlePx, 21);
    diagnosticText((leftPx + 4) / w, 1 - (topPx - 16) / h, titles[panel], titlePx, 11, true);

    photonjet_science::annotation(*gPad,"depth/reference",{{"value",nominalDepth}},scanTrace.then("reference_depth"));
    if (x0 < nominalDepth && nominalDepth < x1) {  /* the nominal depth, up to the labels it would cross */
      const double under = xa < nominalDepth && nominalDepth < xb ? blockHeight : 0;
      const double x = (leftPx + (nominalDepth - x0) / (x1 - x0) * frameW) / w;
      dashedLine(x, bottomPx / h, x, (bottomPx + frameH - under) / h, reference->GetLineColor(),
                 reference->GetLineWidth(), true);
    }
    if (panel < 2) {
      auto *training = static_cast<TGraph *>(train->Clone()), *heldOut = static_cast<TGraph *>(holdout->Clone());
      for (std::size_t i = 0; i < scan.size(); ++i) {
        training->SetPoint(i, scan[i].first, scan[i].second[2 * panel]);
        heldOut->SetPoint(i, scan[i].first, scan[i].second[2 * panel + 1]);
      }
      const auto cov=photonjet_science::covariance("not_computed",{},0,"saved performance metrics have no covariance");
      photonjet_science::bindNative(training,"depth/train/"+std::to_string(panel),scanTrace,cov);
      photonjet_science::bindNative(heldOut,"depth/holdout/"+std::to_string(panel),scanTrace,cov);
      photonjet_science::draw(training,"LP SAME");
      photonjet_science::draw(heldOut,"LP SAME");
    } else {
      for (const auto &[depth, v] : scan) {
        photonjet_science::annotation(*gPad,"depth/gaps/"+std::to_string(std::find_if(scan.begin(),scan.end(),[&](const auto& row){return row.first==depth;})-scan.begin()),
          {{"depth",depth},{"auc_train_minus_holdout",v[0]-v[1]},{"loss_holdout_minus_train",v[3]-v[2]}},scanTrace.then("generalization_gaps"));
        auto *auc = new TBox(depth - .32, 0, depth - .02, v[0] - v[1]);
        auc->SetFillColor(cBlue);
        auc->SetLineWidth(0);
        photonjet_science::draw(auc);
        auto *loss = new TBox(depth + .02, 0, depth + .32, v[3] - v[2]);
        loss->SetFillColor(cRed);
        loss->SetLineWidth(0);
        photonjet_science::draw(loss);
      }
      auto *zero = new TLine(x0, 0, x1, 0);
      zero->SetLineWidth(2);
      photonjet_science::draw(zero);
    }
    if (panel == 0) {
      double y = 1 - (topPx + 16 + .8 * markPx) / h;
      diagnosticText((leftPx + inset) / w, y, sPhenixInternalLabel(), markPx);
      for (const auto &text : lines) diagnosticText((leftPx + inset) / w, y -= pitch / h, text, linePx);
    } else if (panel == 1) {
      const double right = 1 - (rightPx + 14) / w, top = 1 - (topPx + 12) / h, width = (blockWidth - 30) / w;
      TLegend *legend = ianLegend(right - width, top - curves.size() * pitch / h, right, top, 1.9 * linePx / w / width);
      legend->SetTextSize(linePx);
      for (const auto &[object, label, option] : curves) legend->AddEntry(object, label.c_str(), option);
      photonjet_science::draw(legend);
    } else legendRow((leftPx + inset) / w, 1 - (topPx + 18 + .65 * linePx) / h, gaps, linePx);
    pad->RedrawAxis();
  }
  saveTrainingDiagnostic(canvas, output, index, names, "training_capacity");
}

/* The retrained p+p photon-ID BDT against the original PPG12 classifier on the same held-out candidates, as in the
   IAN. First figure: the two scores of every candidate as a map on a log colour scale, with the diagonal of equal
   scores, the number of pairs and their correlation. Second figure: the distribution of the score difference for the
   signal and for the background, each normalised to one, with the mean absolute difference of all pairs.
   Read from model_qa/<model>/original_comparison: the tree `paired_scores` (the two scores, their difference and the
   class of every pair) and `signal_difference`, which gives the axis of the difference distributions only: its
   content is not used. Filled and computed here from the tree, not saved: the map, the two difference distributions,
   the correlation and the mean absolute difference. */
void ppg12ModelComparison(const PlotJson &index) {
  const std::string name = index.at("comparisons").at("ppg12").at("model");
  const auto &meta = modelMetadata(name);
  const std::string folder = "model_qa/" + name + "/original_comparison/";
  auto *tree = dynamic_cast<TTree *>(gQA->Get((folder + "paired_scores").c_str()));
  if (!tree || tree->GetEntries() != meta.at("original_comparison").at("rows").get<Long64_t>())
    throw std::runtime_error("saved PPG12 pair count differs");
  auto *original = tree->GetLeaf("original_score"), *retrained = tree->GetLeaf("retrained_score");
  auto *difference = tree->GetLeaf("difference"), *label = tree->GetLeaf("label");
  if (!original || !retrained || !difference || !label) throw std::runtime_error("missing saved paired score columns");
  /* the axes end just above one, so that a score of exactly 1 falls in the last bin */
  const double upper = std::nextafter(1., std::numeric_limits<double>::infinity());
  TH2D pairs(fresh().c_str(), "", 80, 0, upper, 80, 0, upper);
  pairs.SetDirectory(nullptr);
  pairs.SetStats(false);
  auto *templateHist = dynamic_cast<TH1 *>(gQA->Get((folder + "signal_difference").c_str()));
  if (!templateHist) throw std::runtime_error("missing paired difference display axes");
  auto *signal = static_cast<TH1 *>(templateHist->Clone(fresh().c_str()));
  auto *background = static_cast<TH1 *>(templateHist->Clone(fresh().c_str()));
  signal->SetDirectory(nullptr);
  background->SetDirectory(nullptr);
  signal->Reset();
  background->Reset();
  signal->Sumw2();
  background->Sumw2();
  double sumA = 0, sumB = 0, sumAA = 0, sumBB = 0, sumAB = 0, sumAbs = 0;  /* for the two numbers on the figures */
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    const double a = original->GetValue(), b = retrained->GetValue(), delta = difference->GetValue();
    const int klass = label->GetValue();
    if (!std::isfinite(a) || !std::isfinite(b) || a < 0 || a > 1 || b < 0 || b > 1 || !std::isfinite(delta) ||
        std::abs(delta - (b - a)) > 1e-12 || (klass != 0 && klass != 1))
      throw std::runtime_error("invalid saved PPG12 paired scores");
    /* Two model functions on the same held-out candidates: every pair counts once. The learner's weight column is
       not used; the counts are no yields. */
    pairs.Fill(a, b);
    (klass ? signal : background)->Fill(delta);
    sumA += a, sumB += b, sumAA += a * a, sumBB += b * b, sumAB += a * b, sumAbs += std::abs(delta);
  }
  if (pairs.GetEntries() != tree->GetEntries() ||
      signal->GetEntries() + background->GetEntries() != tree->GetEntries())
    throw std::runtime_error("paired diagnostic count differs");
  /* The two numbers the IAN writes for the pairs, each pair counted once: the linear correlation coefficient of the
     two scores and the mean absolute score difference. */
  const double count = tree->GetEntries();
  const double correlation = (count * sumAB - sumA * sumB) /
                             std::sqrt((count * sumAA - sumA * sumA) * (count * sumBB - sumB * sumB));
  const double meanDifference = sumAbs / count;
  const auto trace=photonjet_science::projectionOf(tree,PlotJson::array()).then("paired_scores_unit_counts",{{"count",count}});
  photonjet_science::bindMap(&pairs,"paired_scores/map",trace,true);
  const std::string pairCount = grouped(tree->GetEntries());
  {
    TCanvas canvas(fresh().c_str(), "", 1000, 860);
    photonjet_science::annotation(canvas,"paired_scores/count",photonjet_science::scalar(count,"candidates","tree entries"),trace);
    photonjet_science::annotation(canvas,"paired_scores/correlation",photonjet_science::scalar(correlation,"correlation","Pearson correlation of paired scores"),trace);
    const double w = padWidthPx(), h = padHeightPx();
    const double left = .135, top = .04, bottom = .14, plot = 1 - top - bottom;
    const double right = left + plot * h / w, across = right - left;  /* a square frame */
    canvas.SetMargin(left, 1 - right, bottom, top);
    canvas.SetTicks(1, 1);
    canvas.SetLogz();
    pairs.GetXaxis()->SetTitle("Original PPG12 BaseV3E score");
    pairs.GetYaxis()->SetTitle("Retrained BaseV3E score");
    pairs.GetZaxis()->SetTitle("Held-out candidates / score bin");
    for (TAxis *axis : {pairs.GetXaxis(), pairs.GetYaxis(), pairs.GetZaxis()})
      axis->SetLabelSize(30 / h), axis->SetTitleSize(34 / h), axis->CenterTitle();
    pairs.GetXaxis()->SetNdivisions(505);
    pairs.GetYaxis()->SetNdivisions(505);
    /* the labels a little off the frame: the zeros of the two axes stand apart at the origin */
    pairs.GetXaxis()->SetLabelOffset(.014);
    pairs.GetYaxis()->SetLabelOffset(.012);
    pairs.GetXaxis()->SetTitleOffset(1.5);
    pairs.GetYaxis()->SetTitleOffset(1.5);
    pairs.GetZaxis()->SetTitleOffset(1.35);
    pairs.SetMinimum(1);
    photonjet_science::draw(&pairs,"COLZ");
    /* the diagonal of equal scores, dashed black on a white line: visible on every colour */
    auto *under = new TLine(0, 0, 1, 1);
    under->SetLineColor(kWhite);
    under->SetLineWidth(4);
    photonjet_science::draw(under);
    dashedLine(left, bottom, right, 1 - top, cBlack, 2, true);
    /* The header at the top left on a light panel, as in the IAN. The panel starts in the corner of the frame and ends
       on bin edges, so that it cuts no bin; the bins under it stay visible through it. */
    const std::vector<std::string> lines = {simulationName(true) + ", " + pairCount + " pairs",
                                            Form("Training diagnostic: #rho = %.4f", correlation)};
    const double x = left + kInsetPx / w, y = 1 - top - (14 + kMarkPx) / h;
    const int binsX = pairs.GetNbinsX(), binsY = pairs.GetNbinsY();
    const double panelRight = std::ceil((x + headerWidth(lines) + 12 / w - left) / across * binsX) / binsX;
    const double panelBottom = std::floor((header(x, y, lines, false, false) - 12 / h - bottom) / plot * binsY) / binsY;
    auto *panel = new TBox(0, panelBottom, panelRight, 1);
    panel->SetFillColorAlpha(kWhite, .88);
    panel->SetLineWidth(0);
    photonjet_science::draw(panel);
    gPad->RedrawAxis();
    header(x, y, lines, true, false);
    const char *output = "photonOnly/pp_bdt_original_vs_retrained_scores.png";
    saveTrainingDiagnostic(canvas, output, index, {name}, "paired_score_mapping");
  }
  {
    if (!(signal->Integral() > 0) || !(background->Integral() > 0))
      throw std::runtime_error("empty paired diagnostic class");
    signal->Scale(1 / signal->Integral());
    background->Scale(1 / background->Integral());
    /* The window of the IAN, +-0.32, in whole bins: from the bin that holds -0.31 to the bin that holds +0.31, which
       end at +-0.32 on an axis with bin edges there; wider if a class keeps over 0.2% of its pairs in a bin beyond
       it. */
    const TAxis *axis = signal->GetXaxis();
    int first = std::max(1, axis->FindFixBin(-.31)), last = std::min(axis->GetNbins(), axis->FindFixBin(.31));
    for (int bin = 1; bin <= axis->GetNbins(); ++bin)
      if (std::max(signal->GetBinContent(bin), background->GetBinContent(bin)) > .002)
        first = std::min(first, bin), last = std::max(last, bin);
    const double xa = axis->GetBinLowEdge(first), xb = axis->GetBinUpEdge(last);
    TCanvas canvas(fresh().c_str(), "", 1000, 860);
    style(signal, cRed, 20);
    style(background, cBlue, 21);
    signal->SetLineWidth(3);
    background->SetLineWidth(3);
    /* Every bin is taken 15% higher than it is: the frame then opens until the peak stands that far under the legend,
       whose sample line would else read as one more step of the histogram. */
    std::vector<Extent> drawn;
    for (int bin = first; bin <= last; ++bin)
      for (const TH1 *hist : {signal, background})
        drawn.push_back({axis->GetBinLowEdge(bin), axis->GetBinUpEdge(bin), 1.15 * hist->GetBinContent(bin)});
    const std::vector<std::string> lines = {simulationName(true), "Training diagnostic: held-out pairs",
                                            Form("mean %s = %.4f", absolute("#Deltas").c_str(), meanDifference)};
    const std::vector<std::string> labels = {"Prompt photon+jet signal", "Nonprompt inclusive-jet background"};
    const std::string ytitle = Form("Fraction / %g score bin", axis->GetBinWidth(first));
    TLegend *legend = ianFrame(xa, xb, 0, drawn, false, "", ytitle.c_str(), lines, labels);
    const double w = padWidthPx(), h = padHeightPx();
    const double left = gPad->GetLeftMargin(), right = 1 - gPad->GetRightMargin();
    const double bottom = gPad->GetBottomMargin(), top = 1 - gPad->GetTopMargin();
    for (TAxis *a : {gFrame->GetXaxis(), gFrame->GetYaxis()})
      a->SetLabelSize(30 / h), a->SetTitleSize(34 / h), a->SetTitleOffset(1.7), a->CenterTitle();
    /* the x title under the middle of its axis, as in the IAN, with its two indices kept on their symbols */
    const char *xtitle = "Score difference, s_{retrained} #minus s_{original}";
    diagnosticText((left + right) / 2, bottom - 96 / h, xtitle, exactTextPx(34), 21);
    /* The dashed line of equal scores, up to the legend or the header if it would cross one of them. (The legend is
       not drawn yet: its corners are still the NDC ones it was made with.) */
    const double zeroX = left + (0 - xa) / (xb - xa) * (right - left);
    double zeroTop = top;
    if (zeroX > legend->GetX1() && zeroX < legend->GetX2()) zeroTop = legend->GetY1() - 8 / h;
    else if (zeroX < left + kInsetPx / w + headerWidth(lines))
      zeroTop = header(left, top - (14 + kMarkPx) / h, lines, false, false) - 8 / h;
    dashedLine(zeroX, bottom, zeroX, zeroTop, cBlack, 2, true);
    const auto cov=photonjet_science::covariance("not_computed",{},0,"unit area histograms retain diagonal errors");
    photonjet_science::bindNative(signal,"paired_scores/signal_difference",trace.then("label_selection_and_unit_area",{{"label",1}}),cov);
    photonjet_science::bindNative(background,"paired_scores/background_difference",trace.then("label_selection_and_unit_area",{{"label",0}}),cov);
    photonjet_science::annotation(canvas,"paired_scores/mean_absolute_difference",photonjet_science::scalar(meanDifference,"score","mean absolute paired score difference"),trace);
    photonjet_science::draw(signal,"HIST SAME");
    photonjet_science::draw(background,"HIST SAME");
    legend->AddEntry(signal, labels[0].c_str(), "l");
    legend->AddEntry(background, labels[1].c_str(), "l");
    photonjet_science::draw(legend);
    gPad->RedrawAxis();
    const char *output = "photonOnly/pp_bdt_retraining_score_difference.png";
    saveTrainingDiagnostic(canvas, output, index, {name}, "paired_score_difference");
  }
}

/* ---- saved photon-ID models: score regions --------------------------------------------------------------------- */

/* The photon-ID score of prompt photons and of inclusive-jet candidates in embedded Au+Au simulation, one panel per
   centrality class on one y axis, as in the IAN: the non-tight, excluded and tight score regions are shaded between
   their cuts at the middle of the class, and each region carries the saved number of candidates in it, for which the
   cuts are applied candidate by candidate.
   Read from model_qa/<default model>/application: the record `review` (the classes with their region counts, the
   working points, the calibration domain and the sample normalisation) and, per class i, score_<i>_signal and
   score_<i>_background or score_<i>_truth_background, each drawn normalised to unit area.
   The call without argument draws two figures: every inclusive-jet candidate as the background and then, by calling
   itself with truthBackground, the truth-tagged background.
   The number written in a region counts the candidates of the two populations the figure draws, the prompt photons
   and its background, each candidate sorted by the cuts at its own centrality: below a score of 0.1, non-tight up to
   the cut of the 80% working point, excluded up to the cut of the 70% working point, tight above it. The record
   holds these numbers per class, under `regions` for the prompt photons with every inclusive-jet candidate and under
   `regions_truth_background` for the prompt photons with the truth-tagged background. No histogram of the file gives
   them: the score histograms hold weighted sums in bins that the centrality-dependent cuts cross. A figure whose
   record lacks the counts of its populations carries the region names alone.
   DATA: when the folder holds data_tight_isolated_score_<i> (the tight, isolated Au+Au DATA candidates of the
   measurement, counts per score bin; FinalAnalysis has no part in making them), they are drawn as points over the
   curves with one legend entry. Tight is a cut on this score, so they lie above the tight cut; they are scaled to the
   prompt curve over the bins wholly above the highest tight cut of the class, which compares the two shapes there.
   Returns false, drawing nothing, when the file has no saved model with an application record: bdtPerClass, which
   calls this function first, then draws photonOnly/auau_bdt_score_perCent.png from the package QA in its own
   layout. */
bool savedScoreReview(bool truthBackground) {
  const auto index = savedJson("model_qa/collection");
  if (index.is_null()) return false;
  const std::string name = index.at("default_model");
  auto *dir = gQA->GetDirectory(("model_qa/" + name + "/application").c_str());
  if (!dir) return false;
  const auto review = savedJson("model_qa/" + name + "/application/review");
  if (review.is_null()) return false;
  const auto weights = review.contains("weight_provenance") ? review.at("weight_provenance") : review;
  const std::string output = truthBackground ? "photonOnly/auau_bdt_score_truth_background_perCent.png"
                                             : "photonOnly/auau_bdt_score_perCent.png";
  /* a figure without physics weights is refused here, before anything is drawn; save() makes the check again */
  photonjet_render::requirePhysicsPlotWeights(weights, output);
  /* The layout holds three classes: the mark and the selection stand in the first frame, the legend in the second
     and the notes in the third. */
  const int panelCount = 3;
  /* The non-tight region runs from a score of nonTightLow to the cut of the 80% working point and the tight region
     from the cut of the 70% working point up; the scores between the two cuts are excluded. */
  const double nonTightLow = .1;
  const int nonTightLevel = 80, tightLevel = 70;
  const auto &cell = review.at("working_points").at("cells").begin().value();
  auto cut = [&](int level, double c) {  /* the cut of a working point at the centrality c */
    const auto &v = cell.at("cuts").at("WP" + std::to_string(level));
    return v.at("intercept").get<double>() + v.at("slope").get<double>() * c;
  };
  /* the region counts of the populations this figure draws, shown when every class of the record holds them */
  const char *regionRecord = truthBackground ? "regions_truth_background" : "regions";
  bool counted = true;
  for (const auto &r : review.at("classes")) counted = counted && r.contains(regionRecord);
  const int blue = TColor::GetColor("#2d71df"), orange = TColor::GetColor("#e78300"),
            purple = TColor::GetColor("#78609b");
  const int darkBlue = TColor::GetColor("#234d80"), darkOrange = TColor::GetColor("#925505");
  const int hatch = TColor::GetColor("#dcdcdc");
  std::vector<std::pair<TH1 *, TH1 *>> histograms;
  std::vector<TH1 *> data;  /* per class, or null */
  double peak = 0, dataPeak = 0;  /* the highest bin of a curve; the highest DATA point with its error */
  for (int i = 0; i < panelCount; ++i) {
    auto *a = dynamic_cast<TH1 *>(dir->Get(Form("score_%d_signal", i)));
    auto *b = dynamic_cast<TH1 *>(
      dir->Get(Form(truthBackground ? "score_%d_truth_background" : "score_%d_background", i)));
    if (!a || !b || a->Integral() <= 0 || b->Integral() <= 0)
      throw std::runtime_error("incomplete saved score review");
    a = static_cast<TH1 *>(a->Clone(fresh().c_str()));
    b = static_cast<TH1 *>(b->Clone(fresh().c_str()));
    a->SetDirectory(nullptr);
    b->SetDirectory(nullptr);
    const auto traceA=photonjet_science::projectionOf(a,PlotJson::array()).then("unit_area",{{"integral",a->Integral()}});
    const auto traceB=photonjet_science::projectionOf(b,PlotJson::array()).then("unit_area",{{"integral",b->Integral()}});
    a->Scale(1 / a->Integral());
    b->Scale(1 / b->Integral());
    const auto cov=photonjet_science::covariance("not_computed",{},0,"unit area histograms retain diagonal errors");
    photonjet_science::attachTrace(a,traceA);photonjet_science::attachTrace(b,traceB);
    photonjet_science::bindNative(a,"score_review/signal/"+std::to_string(i),traceA,cov);
    photonjet_science::bindNative(b,"score_review/background/"+std::to_string(i),traceB,cov);
    peak = std::max({peak, a->GetMaximum(), b->GetMaximum()});
    histograms.push_back({a, b});
    /* the DATA of this class, when the file holds it: scaled to the prompt curve over the bins wholly above the tight
       cut at the far end of the class, where the cut is highest */
    TH1 *points = nullptr;
    if (auto *saved = dynamic_cast<TH1 *>(dir->Get(Form("data_tight_isolated_score_%d", i)))) {
      if (saved->GetNbinsX() != a->GetNbinsX()) throw std::runtime_error("the DATA score bins differ from the review's");
      const int from = a->GetXaxis()->FindFixBin(cut(tightLevel, review.at("classes").at(i).at("high").get<double>())) + 1;
      const double inData = saved->Integral(from, saved->GetNbinsX()), inPrompt = a->Integral(from, a->GetNbinsX());
      if (inData > 0 && inPrompt > 0) {
        points = static_cast<TH1 *>(saved->Clone(fresh().c_str()));
        points->SetDirectory(nullptr);
        const auto dataTrace=photonjet_science::derivedTrace("data_scale_to_prompt",{saved,a},{{"first_bin",from},{"data_integral",inData},{"prompt_integral",inPrompt}});
        points->Scale(inPrompt / inData);
        photonjet_science::bindNative(points,"score_review/data/"+std::to_string(i),dataTrace,cov);
        for (int bin = 1; bin <= points->GetNbinsX(); ++bin)
          dataPeak = std::max(dataPeak, points->GetBinContent(bin) + points->GetBinError(bin));
      }
    }
    data.push_back(points);
  }
  /* Three joined frames on one y axis. The top of every frame holds the header, the shaded regions end below it, and
     the region names with their counts stand in the top of the shading, above every curve. The y axis ends at 1.85
     times the highest bin of any curve and not below 0.235, and high enough for the highest DATA point with its error
     to end under the counts (at 0.57 of the axis); the shading ends at 0.70 of it. */
  const double ymax = std::max({.235, peak * 1.85, dataPeak / .57}), shadeTop = ymax * .70;
  TCanvas cv(truthBackground ? "saved_bdt_truth_background_review" : "saved_bdt_review", "", 1920, 690);
  const double w = cv.GetWw(), h = cv.GetWh(), left = 128, right = 24, top = 12, bottom = 98,
               fw = (w - left - right) / panelCount, fh = h - top - bottom;
  const double bandTop = top + fh * (1 - shadeTop / ymax), lineTop = ymax * (1 - (bandTop - top + 68) / fh),
               inset = 22;
  /* the two parts of the model's sample line: the simulation, for the header, and the thresholds of its photon+jet
     samples, for the legend */
  const std::string simulation = simulationName(isPP(modelMetadata(name))),
                    sources = thresholdList(photonJetThresholds(review));
  const auto &domain = review.at("working_points").at("calibration_domain");
  std::string selection = Form("%g #leq E_{T}^{#gamma} < %g GeV,  ", domain.at("et_min").get<double>(),
                               domain.at("et_max_exclusive").get<double>());
  selection += Form("%s < %g,  ", absolute("#eta^{#gamma}").c_str(), domain.at("eta_abs_max").get<double>());
  selection += Form("|z#lower[0.35]{#scale[0.67]{vtx}}| < %g cm", domain.at("vertex_abs_max_cm").get<double>());
  const auto population = truthBackground ? photonjet_render::CandidatePopulation::TruthTaggedBackground
                                          : photonjet_render::CandidatePopulation::AllInclusiveJetCandidates;
  const std::string binWidth = Form("%g", histograms.at(0).first->GetXaxis()->GetBinWidth(1));
  auto box = [](double x1, double y1, double x2, double y2, int color, double alpha) {
    auto *b = new TBox(x1, y1, x2, y2);
    b->SetFillColorAlpha(color, alpha);
    b->SetLineColorAlpha(color, 0);
    photonjet_science::draw(b);
  };
  bool dataShown = false;  /* set when a class draws its DATA points */
  /* the text of the counts at its largest, and the text of the excluded region */
  const double countFullPx = 19, excludedPx = 18;
  double countPx = countFullPx;
  for (int i = 0; i < panelCount; ++i) {
    const double fx = left + i * fw;
    TH1F *axes = pixelPanel(cv, fx, top, fx + fw, top + fh, 0, 0, 1, ymax, 26);
    axes->GetXaxis()->SetNdivisions(505);
    axes->GetYaxis()->SetNdivisions(506);
    axes->GetXaxis()->SetDecimals();
    axes->GetYaxis()->SetDecimals();
    /* the frames touch: y labels at the first only, and a frame's first x label gives way to its neighbour's last */
    if (i) {
      axes->GetYaxis()->SetLabelSize(0);
      axes->GetXaxis()->ChangeLabel(1, -1, 0);
    }
    const auto &row = review.at("classes").at(i);
    double low = row.at("low"), high = row.at("high"), mid = (low + high) / 2.;
    double nt = cut(nonTightLevel, mid), tight = cut(tightLevel, mid);
    const auto reviewTrace=photonjet_science::recordTrace(*gQA,"model_qa/"+name+"/application/review").then("score_review_regions",{{"class_slot",i},{"background",truthBackground?"truth_tagged":"inclusive"}});
    photonjet_science::annotation(*gPad,"score_review/cuts/"+std::to_string(i),{{"centrality",{low,high}},{"non_tight_low",nonTightLow},{"non_tight_upper",nt},{"tight",tight}},reviewTrace);
    if(counted)photonjet_science::annotation(*gPad,"score_review/counts/"+std::to_string(i),row.at(regionRecord),reviewTrace);
    box(nonTightLow, 0, nt, shadeTop, orange, .11);
    box(nt, 0, tight, shadeTop, kGray, .15);
    box(tight, 0, 1, shadeTop, blue, .10);
    /* the excluded region is hatched faintly, as in the IAN: lines at 45 degrees on the canvas, eight pixels apart */
    const double bandWidth = (tight - nt) * fw, bandHeight = shadeTop / ymax * fh;
    for (double y = 8 - bandWidth; y < bandHeight; y += 8) {
      const double from = std::max(y, 0.0), to = std::min(y + bandWidth, bandHeight);
      auto *l = new TLine(nt + (from - y) / fw, from / fh * ymax, nt + (to - y) / fw, to / fh * ymax);
      l->SetLineColor(hatch);
      l->SetLineWidth(1);
      photonjet_science::draw(l);
    }
    /* the lower edge of the non-tight region and its upper edge dashed, the tight cut solid, as in the IAN */
    pixelDashes(nonTightLow, 0, nonTightLow, lineTop, lineTop / ymax * fh, darkOrange);
    pixelDashes(nt, 0, nt, lineTop, lineTop / ymax * fh, purple);
    auto *tightCut = new TLine(tight, 0, tight, lineTop);
    tightCut->SetLineColor(darkBlue);
    tightCut->SetLineWidth(3);
    photonjet_science::draw(tightCut);
    auto [a, b] = histograms[i];
    a->SetLineColor(blue);
    b->SetLineColor(orange);
    a->SetLineWidth(3);
    b->SetLineWidth(3);
    photonjet_science::draw(a,"HIST SAME");
    photonjet_science::draw(b,"HIST SAME");
    TH1 *points = data[i];
    if (points) {
      points->SetMarkerStyle(20);
      points->SetMarkerSize(1.25);
      points->SetMarkerColor(cBlack);
      points->SetLineColor(cBlack);
      points->SetLineWidth(2);
      photonjet_science::draw(points,"E1 X0 SAME");
      dataShown = true;
    }
    pixelLabel(fx + inset, top + 52, centralityRange(low, high), 40, 11, cBlack, true);
    if (i == 0) {
      pixelLabel(fx + inset, top + 88, sPhenixInternalLabel(), 30);
      pixelLabel(fx + inset, top + 117, simulation, 22);
      pixelLabel(fx + inset, top + 145, selection, 21);
    }
    if (i == 1) {
      TLegend *legend = ianLegend((fx + inset) / w, 1 - (top + (points ? 164 : 132)) / h, (fx + fw - inset) / w,
                                  1 - (top + 68) / h, .13);
      legend->SetTextSize(21);
      legend->AddEntry(a, ("Prompt photon: photon+jet " + sources).c_str(), "l");
      /* the inclusive-jet samples behind the background curve, by their thresholds in GeV: typed here, not read from
         the file */
      candidateEntry(legend, b, population, "inclusive+jet", "l", {12, 20, 30, 40});
      if (points) legend->AddEntry(points, "DATA: tight, isolated", "pe");
      photonjet_science::draw(legend);
    }
    if (i == 2) {
      /* with DATA beside the simulated curves the notes say which of the two each describes */
      pixelLabel(fx + inset, top + 88, dataShown ? "Simulation: preselected photon candidates," : "Preselected photon candidates",
                 21, 11, kGray + 3);
      pixelLabel(fx + inset, top + 117, dataShown ? "no isolation requirement; each class unit-normalized"
                                                  : "No isolation requirement; each class unit-normalized", 21, 11, kGray + 3);
      if (dataShown)
        pixelLabel(fx + inset, top + 146, "DATA: scaled to the prompt curve above the tight cut", 21, 11, kGray + 3);
    }
    /* The region names in the top of the shading, their counts under them; the tight region is centred between its cut
       and the ticks at the right edge of the frame, and the excluded region reads upwards. The counts have one size in
       all three panels, at which the longest stands clear of its cut and of the ticks. `ticks` is the y tick length
       the frame carries, a share of the pad's width, times that width. */
    const double ticks = axes->GetYaxis()->GetTickLength() * w;
    const double nonTightX = fx + (nonTightLow + nt) / 2 * fw, excludedX = fx + (nt + tight) / 2 * fw,
                 tightX = fx + (tight + 1) / 2 * fw - ticks / 4;
    const double excludedY = top + fh * .66;
    if (i == 0 && counted)
      for (const auto &r : review.at("classes")) {
        const double c = (r.at("low").get<double>() + r.at("high").get<double>()) / 2;
        const std::string longest = "N = " + grouped(r.at(regionRecord).at("tight").get<long long>());
        countPx = std::min(countPx, countFullPx * ((1 - cut(tightLevel, c)) * fw - ticks - 14) /
                                      (textWidth(longest, countFullPx) * w));
      }
    pixelLabel(nonTightX, bandTop + 31, "Non-tight", 25, 21, darkOrange, true);
    pixelLabel(tightX, bandTop + 31, "Tight", 25, 21, darkBlue, true);
    std::string excluded = "#bf{Excluded}";
    if (counted) {
      const auto count = [&](const char *region) {
        return "N = " + grouped(row.at(regionRecord).at(region).get<long long>());
      };
      pixelLabel(nonTightX, bandTop + 57, count("non_tight"), countPx, 21, darkOrange);
      pixelLabel(tightX, bandTop + 57, count("tight"), countPx, 21, darkBlue);
      excluded += Form("   #scale[%.3f]{%s}", countPx / excludedPx, count("gap").c_str());
    }
    pixelLabel(excludedX, excludedY, excluded, excludedPx, 22, kGray + 3, false, 90);
    pixelLabel(fx + fw / 2, h - 34, "BDT score, #it{s}_{BDT}", 32, 22);
    gPad->RedrawAxis();
  }
  cv.cd();
  pixelLabel(32, top + fh / 2, "Weighted fraction of candidates / " + binWidth, 32, 22, cBlack, false, 90);
  saveModelQA(cv, output, index, {name}, weights);
  amendReceipt(output, [&](PlotJson &receipt) {
    /* Typed here, not read from the file: the selection of the saved score histograms as the PhotonID export fills
       them. */
    receipt["score_population"] = {
      {"signal", "exclusive source-role-1 prompt-photon winners"},
      {"background", truthBackground ? "truth-tagged nonprompt inclusive-jet candidates"
                                     : "all selected inclusive-jet candidates"},
      {"background_histogram", truthBackground ? "score_i_truth_background" : "score_i_background"},
      {"layout", "one row, three centrality panels; no lower cut diagram"},
      {"region_counts_shown", counted},
      {"region_counts_record", counted ? regionRecord : "not saved"},
      {"preselection", "pho_presel_active == 1"},
      {"isolation_cut", false},
      {"centrality_reweighting", false},
      {"vertex_reweighting", false}};
    if (dataShown)
      receipt["data_overlay"] = {
        {"histograms", "data_tight_isolated_score_i"},
        {"population", "tight, isolated Au+Au DATA candidates of the measurement"},
        {"scaling", "to the prompt curve over the bins wholly above the highest tight cut of the class"},
        {"record", savedJson("model_qa/" + name + "/application/data_overlay")}};
    /* Both figures of this function carry one recipe, auau_bdt_score: the truth-tagged background is its baseline and
       every inclusive-jet candidate its comparison. */
    receipt["plot_recipe"] = plotRecipe("auau_bdt_score", truthBackground ? "canonical" : "all_inclusive_jet",
                                        truthBackground ? "canonical" : "comparison", output,
                                        "photonOnly/auau_bdt_score_truth_background_perCent.png", "savedScoreReview");
    receipt["plot_recipe"]["selector"] = truthBackground ? "truth_tagged_background" : "all_inclusive_jet_candidates";
    receipt["plot_recipe"]["title"] = truthBackground ? "Truth-tagged background" : "All inclusive-jet candidates";
    if (truthBackground)
      receipt["comparison_outputs"] = PlotJson::array(
        {{{"id", "all_inclusive_jet"},
          {"title", "All inclusive-jet candidates"},
          {"output", "photonOnly/auau_bdt_score_perCent.png"},
          {"difference",
           "Orange includes all selected inclusive-jet candidates instead of only truth-tagged background."}}});
  });
  if (!truthBackground) savedScoreReview(true);
  return true;
}

/* ---- the figures of the QA file, in drawing order -------------------------------------------------------------- */

/* From a figure of the note to the function of this part of the file that draws it (the photon and event figures are
   drawn by the functions further up; qaFigures lists them all):
     simulation_stitching/<pp|auau>_<photonjet|inclusivejet>_<generator|reco_tag>_stitching.png      stitching
     event_obs/pp_max_cluster_energy_overlay.png, ..._turnon.png                                      triggerStudy
     event_obs/hMaxClusterEnergy_groupOverlay.png, ..._groupTurnOnOverlay.png                         triggerStudy
     photonOnly/auau_h70_bdt_pt_centrality_grid.png, isolation/auau_r<NN>_isolation_..._grid.png      workingPoints
     isolation/auau_r<NN>_isolation_centrality_fits.png                                               workingPoints
     photonOnly/auau_h70_bdt_centrality_fits.png                        workingPoints or savedWorkingPointFits
     photonOnly/auau_h70_bdt_centrality_comparison.png                  workingPointCentralityComparison
     photonOnly/<pp|auau>_training_<et|eta>_reweighting.png             trainingReweighting
     photonOnly/<pp|auau>_bdt_single_feature_separation.png, ..._feature_split_gain.png               featureDiagnostics
     photonOnly/auau_bdt_routing_auc_by_et_centrality.png               routingDiagnostics
     photonOnly/auau_bdt_tree_depth_train_holdout_diagnostic.png        depthDiagnostics
     photonOnly/pp_bdt_original_vs_retrained_scores.png, photonOnly/pp_bdt_retraining_score_difference.png
                                                                        ppg12ModelComparison
     photonOnly/auau_bdt_score_perCent.png             savedScoreReview, or bdtPerClass without a saved model
     photonOnly/auau_bdt_score_truth_background_perCent.png             savedScoreReview
     photonOnly/pp_auau_bdt_heldout_score_separation.png, photonOnly/pp_auau_bdt_heldout_roc.png      scoreSeparation
   The figures of a saved model other than the default one are written under models/<name>/. */

/* The figures of the saved photon-ID models (model_qa/ of the QA file). Per model with a saved feature table: the
   training reweighting and the two input diagnostics; a model other than the default one writes them under
   models/<name>/. For a p+p collection: the comparison with the PPG12 classifier where it is saved; the two figures
   of the p+p against the Au+Au score are listed as not drawn, since they need the physics-weighted score QA of the
   two packages (scoreSeparation). For an Au+Au collection: the threshold grid, the centrality comparison and the
   lines of the working points, and the routing comparison and the depth scan where their tables are saved.
   While these figures are drawn the input-coverage word of their receipts is SAVED_MODEL_QA, with which the header
   writes no note of an incomplete input population on their canvases. */
void modelQAFigures() {
  const auto index = savedJson("model_qa/collection");
  if (index.is_null()) return;
  const std::string nominal = index.at("default_model");
  const std::string population = gPopulation;
  gPopulation = "SAVED_MODEL_QA";
  const std::string reason =
    "saved learner QA has training or unnormalized candidate weights; physics event/source provenance required";
  for (const auto &entry : index.at("models")) {
    const std::string name = entry.at("name"), prefix = name == nominal ? "" : "models/" + name + "/";
    /* a model with a saved feature table has the training diagnostics */
    if (gQA->Get(("model_qa/" + name + "/features").c_str())) {
      trainingReweighting(name, index, prefix);
      featureDiagnostics(name, index, prefix);
    }
  }
  if (isPP(index)) {
    for (const std::string stem : {"pp_auau_bdt_heldout_score_separation", "pp_auau_bdt_heldout_roc"})
      missing("photonOnly/" + stem + ".png", reason);
    if (index.at("comparisons").contains("ppg12")) ppg12ModelComparison(index);
  } else {
    workingPoints("bdt");
    workingPointCentralityComparison(nominal, index);
    savedWorkingPointFits(nominal, index);
    if (gQA->Get(("model_qa/" + nominal + "/holdout_et_c3").c_str())) routingDiagnostics(index);
    if (gQA->Get(("model_qa/" + nominal + "/performance").c_str())) depthDiagnostics(index);
  }
  gPopulation = population;
}

/* The photon-ID score of the p+p model against that of the Au+Au model in the most peripheral class, as in the IAN,
   from the simulated candidates of the pass after the preselection with their physics weights: prompt photons of the
   photon+jet simulation and the other candidates of the inclusive-jet simulation. Read from
   qa/<pkg>/phosim/slices/sig/bdt and qa/<pkg>/jetsim/slices/bkg/bdt of the two packages. First figure: the four
   score distributions, each normalised to unit area, on a log scale. Second figure: for a cut on the score, the
   signal efficiency against the background acceptance, with the area under each curve.
   The curves of the second figure and their areas are computed here from the binned scores, not saved by the qa
   step: the cut runs over the bin edges of the saved score axis, and the area is the trapezoid sum over them. */
void scoreSeparation() {
  const std::string scores = "photonOnly/pp_auau_bdt_heldout_score_separation.png",
                    roc = "photonOnly/pp_auau_bdt_heldout_roc.png";
  const auto classes = classesOf("AuAu");
  struct System { std::string pkg, c, label; int style, color; Bins signal, background; };
  std::vector<System> systems = {{"PP", "", "p+p", 1, cBlue, {}, {}}};
  if (!classes.empty()) {
    const auto peripheral = *std::max_element(classes.begin(), classes.end(),
                                              [](const std::string &a, const std::string &b) { return lowOf(a) < lowOf(b); });
    systems.push_back({"AuAu", peripheral, classLabel(peripheral), kDashed, cOrange, {}, {}});
  }
  /* one photon window for both systems: the range every model is applied in */
  double low = 0, high = 1e9;
  for (const auto &s : systems)
    low = std::max(low, photonSupport(s.pkg).first), high = std::min(high, photonSupport(s.pkg).second);
  for (auto &s : systems) {
    TH1 *signal = qa(s.pkg + "/phosim/slices/sig/bdt"), *background = qa(s.pkg + "/jetsim/slices/bkg/bdt");
    if (!signal || !background) {
      missing(scores, "needs the photon+jet and inclusive-jet score QA of p+p and Au+Au");
      return missing(roc, "needs the photon+jet and inclusive-jet score QA of p+p and Au+Au");
    }
    std::vector<Cut> window = {{1, low, high}};
    /* axis 2 of the saved histogram counts the centrality in slices of 5%: the class 50-80% is its range 10 to 16 */
    if (!s.c.empty()) window.push_back({2, lowOf(s.c) / 5, highOf(s.c) / 5});
    s.signal = unit(project(signal, 0, window)), s.background = unit(project(background, 0, window));
  }
  if (systems.size() < 2) { missing(scores, "needs p+p and Au+Au"); return missing(roc, "needs p+p and Au+Au"); }
  const std::vector<std::string> lines = {"p+p, " + energyOf("PP"), systems[1].label + ", " + energyOf("AuAu"),
                                          "Simulation, " + photonEtaLabel("AuAu"), photonWindowLabel(low, high)};

  {  /* the score distributions in bins of 0.02 */
    const int merge = std::max(1, int(std::lround(0.02 / systems[0].signal.width(0))));
    std::vector<Extent> drawn;
    std::vector<std::pair<Bins, std::pair<int, int>>> curves;  /* distribution, (colour, line style) */
    for (const auto &s : systems) {
      curves.push_back({rebin(s.signal, merge), {cRed, s.style}});
      curves.push_back({rebin(s.background, merge), {cBlue, s.style}});
      extents(drawn, curves[curves.size() - 2].first, false, 0, 1), extents(drawn, curves.back().first, false, 0, 1);
    }
    TCanvas cv("cv", "", 800, 740);
    const std::vector<std::string> labels = {"Prompt photon+jet", "Nonprompt inclusive jet", systems[0].label + " model",
                                             "Au+Au model"};
    TLegend *legend = frame(0, 1, drawn, 1e-4, true, "System-specific BDT score",
                            Form("Weighted fraction / %g score bin", merge * systems[0].signal.width(0)), lines, labels);
    gFrame->GetXaxis()->SetNdivisions(505);
    for (const auto &[bins, look] : curves) photonjet_science::draw(line(bins, look.first, look.second, 3),"HIST SAME");
    const std::array<std::pair<int, int>, 4> keys = {{{cRed, 1}, {cBlue, 1}, {cBlack, 1}, {cBlack, kDashed}}};
    for (size_t i = 0; i < labels.size(); ++i) {
      auto *key = new TLine;
      key->SetLineColor(keys[i].first), key->SetLineStyle(keys[i].second), key->SetLineWidth(3);
      legend->AddEntry(key, labels[i].c_str(), "l");
    }
    photonjet_science::draw(legend);
    savePhysicsQA(cv, scores);
  }

  /* signal efficiency against background acceptance for a cut at every bin edge of the stored score axis; a second
     figure of the same packages: each saved figure names its inputs */
  for (const auto &s : systems) renderSelection(s.pkg, true);
  TCanvas cv("cv", "", 800, 740);
  TH1F *axes = photonjet_science::drawFrame(*gPad,0, 0, 1, 1.03);
  axes->GetXaxis()->SetTitle("Background acceptance");
  axes->GetYaxis()->SetTitle("Signal efficiency");
  axes->GetXaxis()->SetNdivisions(505), axes->GetYaxis()->SetNdivisions(506);
  const double w = padWidthPx(), h = padHeightPx(), right = 1 - gPad->GetRightMargin(), bottom = gPad->GetBottomMargin();
  auto *chance = new TLine(0, 0, 1, 1);
  chance->SetLineStyle(kDashed), chance->SetLineColor(kGray + 2), chance->SetLineWidth(2);
  photonjet_science::draw(chance);
  std::vector<std::string> labels;
  std::vector<TGraph *> graphs;
  for (const auto &s : systems) {
    auto *g = new TGraph;
    double passSignal = 0, passBackground = 0, area = 0;
    g->SetPoint(0, 0, 0);
    for (int i = s.signal.n() - 1; i >= 0; --i) {  /* from the highest score down: the cut loosens */
      const double backgroundBefore = passBackground, signalBefore = passSignal;
      passSignal += s.signal.v[i], passBackground += s.background.v[i];
      area += 0.5 * (passSignal + signalBefore) * (passBackground - backgroundBefore);
      g->SetPoint(g->GetN(), passBackground, passSignal);
    }
    const auto rocTrace=photonjet_science::Trace::combine("cumulative_score_roc",{s.signal.physics,s.background.physics},{{"integration","descending score edges; trapezoid"}});
    photonjet_science::bindNative(g,"score_roc",rocTrace,photonjet_science::covariance("not_computed",{},0,"ROC and AUC uncertainty not calculated"));
    photonjet_science::annotation(*gPad,"score_auc",photonjet_science::scalar(area,"1","trapezoid ROC integral"),rocTrace);
    g->SetLineColor(s.color), g->SetLineWidth(3), g->SetLineStyle(1);
    photonjet_science::draw(g,"L");
    graphs.push_back(g);
    labels.push_back(Form("%s, AUC = %.3f", s.label.c_str(), area));
  }
  labels.push_back("Random classifier");
  /* the curves run along the top left, so the legend stands at the bottom right and the header above it; the legend
     as that of ianFrame: a symbol of 46 px, 10 px before the text, 14 px from the right axis */
  double text = 0;
  for (const auto &l : labels) text = std::max(text, textWidth(l, kLinePx));
  const double symbol = 46 / w, width = symbol + text + 10 / w, x2 = right - 14 / w, y1 = bottom + 18 / h;
  TLegend *legend = ianLegend(x2 - width, y1, x2, y1 + labels.size() * kStepPx / h, symbol / width);
  for (size_t i = 0; i < graphs.size(); ++i) legend->AddEntry(graphs[i], labels[i].c_str(), "l");
  legend->AddEntry(chance, labels.back().c_str(), "l");
  photonjet_science::draw(legend);
  const double headerBottom = y1 + (labels.size() + 1) * kStepPx / h;
  header(right - 14 / w - headerWidth(lines), headerBottom + (lines.size() + 0.15) * kStepPx / h, lines);
  savePhysicsQA(cv, roc);
}

/* Draws every QA figure the file has inputs for, from the saved quantities. A figure of the package QA without its
   inputs is listed through missing(). The saved-model figures list a comparison without its models through missing()
   and leave out a figure whose optional table is absent; a table, histogram or record of a saved model that a figure
   needs and that is absent, incomplete or inconsistent mostly stops the run with an error (see the saved-model
   sections above). */
void qaFigures() {
  modelQAFigures();
  for (const std::string pkg : {"AuAu", "PP"}) triggerStudy(pkg);
  for (const std::string pkg : {"PP", "AuAu"})
    for (const std::string family : {"phosim", "jetsim"})
      for (bool generator : {true, false}) stitching(pkg, family, generator);
  centralityReweighting();
  vertexReweighting();
  eventProperties();
  minimumBias();
  ncbScore("PP");
  ncbScore("AuAu");
  ncbEtaTime();
  for (const std::string pkg : {"PP", "AuAu"}) ncbTimeSelections(pkg), ncbEtaTimeAfter(pkg);
  showerShapes();
  widths();
  workingPoints("bdt");
  bdtPerClass();
  scoreSeparation();
  e11Flow();
  truthIsolation();
  truthPtAfterIsolation();
  isolationOverlays();
  workingPoints("isolation");
  isolationStacks();
  abcdCounts();
  /* the numbers given to purityPanels: the least upper end of the upper y axis, and the least range of the lower */
  purityPanels("photonOnly/data_purity_and_prompt_leakage_2x4_style.png", "purity", "DATA photon purity", 1.2,
               "Prompt-signal leakage fraction", 0, 0.6);
  purityPanels("photonOnly/purity_qa_pp_auau_simCLOSURE.png", "closure", "Photon purity (simulation)", 1.3,
               "Corrected ABCD / truth", 0.4, 1.7);
  efficiency();
}

}  /* namespace */
