/* IANFigures: every FinalAnalysis figure, in sPHENIX style and the IAN's layout, under the IAN's own file names. */

/* ./FinalAnalysis/final_analysis.sh draw --results results.root --qa qa.root --fig FinalAnalysis/fig */

/* The measurement figures draw results.root; the photon and event QA figures (QAFigures.C) draw the QA file. Each
   figure draws saved quantities and computes nothing beyond projections, ratios and normalisations for display.
   The results, their labels and their reported xJ ranges all come from results.root, so a new binning or centrality
   class needs no change here. Every xJ figure shows only the reported range, where the distribution is fully
   populated (results.root derives it from the jet threshold and the lowest photon pT); its axis starts at 0 and ends
   at the range's upper edge rounded up to 0.1. The output folder mirrors the IAN's fig/ folder.

   Every figure follows the IAN's layout: the sPHENIX mark and its physics lines at the top left inside the frame, the
   legend at the top right inside the frame, a y range opened until no drawn point reaches either block, and one
   palette. The physics lines are compiled from the saved metadata and histograms, never typed. */
#include "sPhenixStyle.C"
#include "../src/RenderMetadata.h"

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
#include <TObjString.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <array>
#include <cctype>

namespace {

/* One result of results.root, its IAN file-name tags and its reported xJ range. */
struct Result {
  std::string key, package, tag, longTag, label, energy;
  double low = 0, high = 0, jetThreshold = 0, centralityLow = -1;  /* centralityLow: Au+Au only */
  int pairIterations = 0, photonIterations = 0;
  std::string combinatoric;  /* none, simulation, or data_dphi_<shape>_<scope> */
};

/* ---- the IAN look ------------------------------------------------------------------------------------------------
   The palette is the IAN's: black for DATA, region A and truth; blue for corrected distributions and p+p; vermilion
   for background estimates; red for Au+Au and unfolded simulation; green, orange and purple for the iteration-rule
   terms. Text of the header and the legend is sized in pixels, so every figure and every pad reads alike. */
int gBlue = kAzure + 2, gVermilion = kOrange + 7, gRed = kRed + 1, gGreen = kGreen + 2, gOrange = kOrange + 1,
    gPurple = kViolet + 1, gGray = kGray + 2;
void ianPalette() {
  gBlue = TColor::GetColor("#0072B2"), gVermilion = TColor::GetColor("#D55E00"), gRed = TColor::GetColor("#C62828");
  gGreen = TColor::GetColor("#2CA02C"), gOrange = TColor::GetColor("#FF7F0E"), gPurple = TColor::GetColor("#6A3D9A");
  gGray = TColor::GetColor("#6C7686");
}
int seriesColor(int i) {
  const int c[] = {kBlack, gBlue, gVermilion, gRed, gGreen, gOrange, gPurple, gGray};
  return c[i % 8];
}
const int kMarkers[] = {20, 21, 33, 34, 24, 25, 27, 28};
const double kMarkPx = 28, kLinePx = 22, kStepPx = 29;  /* sPHENIX mark, header and legend lines, line pitch */

std::string gFig;
std::vector<std::string> gWritten, gNotDrawn;
bool gLegacy = false;
std::string gPopulation;
std::unique_ptr<TFile> gSystematics;  /* systematics.root of final_analysis.sh systematics, when given */
using PlotJson = nlohmann::json;
PlotJson gResultMetadata, gQAMetadata, gRenderInputs, gStylePins;
std::vector<std::string> gVisibleText;
std::set<std::string> gSemanticPackages;
PlotJson gProjectionSteps = PlotJson::array();
TH1F *gFrame = nullptr;  /* the last frame drawn by ianFrame and its y range */
double gFrameLow = 0, gFrameHigh = 1;
int gUnique = 0;

std::string fresh() { return Form("ian%d", ++gUnique); }

void missing(const std::string &name, const std::string &why) {
  gNotDrawn.push_back(name + ": " + why);
  printf("IANFigures: not drawn %s (%s)\n", name.c_str(), why.c_str());
}

PlotJson loadRenderMetadata(TFile &file, const std::string &path) {
  auto *saved = dynamic_cast<TObjString *>(file.Get("render_metadata"));
  gRenderInputs[path] = photonjet_render::sha256(path, true);
  if (!saved) return PlotJson();  /* historical references are explicitly diagnostic */
  const PlotJson sealed = PlotJson::parse(saved->GetString().Data());
  const PlotJson receipt = photonjet_render::read(path+".render.json");
  if (receipt.at("product_sha256") != gRenderInputs.at(path) || receipt.at("metadata_sha256") != sealed.at("sha256"))
    throw std::runtime_error("render product or metadata changed: "+path);
  return photonjet_render::verify(sealed);
}
const PlotJson &renderSelection(const std::string &pkg, bool qa = false) {
  static const PlotJson missingSelection;
  const auto &metadata = qa ? gQAMetadata : gResultMetadata;
  gSemanticPackages.insert((qa ? "qa/" : "results/")+pkg);
  if (metadata.is_null() || !metadata.at("packages").contains(pkg)) return missingSelection;
  return metadata.at("packages").at(pkg).at("selection");
}
std::string energyLabel(const std::string &pkg, bool qa = false) {
  const auto &m = renderSelection(pkg, qa);
  if (m.is_null()) return "";
  return Form(pkg == "PP" ? "#sqrt{s} = %g GeV" : "#sqrt{s_{NN}} = %g GeV", m.at("collision").at("sqrt_s_nn_gev").get<double>());
}
std::string collisionLabel(const std::string &pkg, const PlotJson &range) {
  const std::string system = pkg == "PP" ? "p+p" : "Au+Au";
  return range.is_null() ? system : system+Form(" %g-%g%%", range.at(0).get<double>(), range.at(1).get<double>());
}
std::string cutOperator(const PlotJson &value) {
  const auto text = value.get<std::string>();
  if (text == ">=") return "#geq";
  if (text == ">") return ">";
  throw std::runtime_error("unsupported cut operator "+text);
}

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

/* The results in file order. IAN names: p+p is "pp"; Au+Au class centA_B is "auauAB" and "auau_AtoB". */
std::vector<Result> resultsOf(TFile &f) {
  std::vector<Result> out;
  TIter next(f.GetListOfKeys());
  while (auto *k = static_cast<TKey *>(next())) {
    if (std::string(k->GetClassName()) != "TDirectoryFile") continue;
    auto s = settingsOf(f, k->GetName());
    if (s.empty()) continue;
    Result r;
    r.key = k->GetName();
    r.package = s["package"];
    r.label = s["label"];
    const std::string c = s["centrality"];
    const auto u = c.find('_');
    if (s["package"] == "PP") {
      r.tag = r.longTag = "pp";
      r.energy = energyLabel("PP");
    } else {
      const bool named = c.rfind("cent", 0) == 0 && u != std::string::npos;
      r.tag = named ? "auau" + c.substr(4, u - 4) + c.substr(u + 1) : r.key;
      r.longTag = named ? "auau_" + c.substr(4, u - 4) + "to" + c.substr(u + 1) : r.key;
      r.energy = energyLabel(s["package"]);
      if (named) r.centralityLow = std::atof(c.substr(4, u - 4).c_str());
    }
    r.jetThreshold = std::atof(s["jet_threshold_gev"].c_str());
    std::istringstream range(s["reported_xj"]);
    char comma;
    range >> r.low >> comma >> r.high;
    r.pairIterations = std::atoi(s["pair_iterations"].c_str());
    r.combinatoric = s["combinatoric"];
    r.photonIterations = std::atoi(s["photon_iterations"].c_str());
    if (!gResultMetadata.is_null()) {
      const auto &metadata = gResultMetadata.at("results").at(r.key);
      if (metadata.at("package").get<std::string>() != r.package || metadata.at("centrality").get<std::string>() != c ||
          metadata.at("recoil_min_gev").get<double>() != r.jetThreshold ||
          metadata.at("reported_xj").at(0).get<double>() != r.low || metadata.at("reported_xj").at(1).get<double>() != r.high)
        throw std::runtime_error("result settings differ from the saved annotation metadata: "+r.key);
      r.label = collisionLabel(r.package, metadata.at("centrality_percent"));
      if (!metadata.at("centrality_percent").is_null()) r.centralityLow = metadata.at("centrality_percent").at(0);
    }
    out.push_back(r);
  }
  return out;
}

TH1 *get(TFile &f, const Result &r, const std::string &name) {
  auto *h = dynamic_cast<TH1 *>(f.Get((r.key + "/" + name).c_str()));
  if (!h) return nullptr;
  h = static_cast<TH1 *>(h->Clone((r.key + name + fresh()).c_str()));
  h->SetDirectory(nullptr);
  return h;
}

std::vector<double> edgesOf(const TAxis *a) {
  std::vector<double> e;
  for (int i = 1; i <= a->GetNbins() + 1; ++i) e.push_back(a->GetBinLowEdge(i));
  return e;
}

/* ---- the IAN's physics lines, compiled from the saved metadata and histograms ---------------------------------- */
std::string systemLine(const Result &r) { return r.label + (r.energy.empty() ? "" : ", " + r.energy); }
/* The photon bins of a result: the metadata's photon_edges_gev, or for reference inputs without metadata the saved
   photon axis of the result (f may be null when only metadata is wanted). */
std::vector<double> photonEdges(TFile *f, const Result &r) {
  if (!gResultMetadata.is_null() && gResultMetadata.at("results").contains(r.key))
    return gResultMetadata.at("results").at(r.key).at("photon_edges_gev").get<std::vector<double>>();
  TH1 *h = f ? get(*f, r, "photons_corrected") : nullptr;
  return h ? edgesOf(h->GetXaxis()) : std::vector<double>{};
}
std::vector<double> photonEdges(TFile &f, const Result &r) { return photonEdges(&f, r); }
std::string photonLine(TFile *f, const Result &r) {
  const auto e = photonEdges(f, r);
  return e.size() < 2 ? "" : Form("%g #leq p_{T}^{#gamma} < %g GeV", e.front(), e.back());
}
std::string photonLine(TFile &f, const Result &r) { return photonLine(&f, r); }
/* The recoil-jet definition: the anti-kT radius from the metadata (absent for reference inputs) and the threshold. */
std::string jetLine(const Result &r) {
  const auto &s = renderSelection(r.package);
  if (s.is_null()) return Form("p_{T}^{jet} #geq %g GeV", r.jetThreshold);
  const auto &j = s.at("jet");
  return Form("anti-k_{T} R = %g, p_{T}^{jet} %s %g GeV", j.at("radius").get<double>(), cutOperator(j.at("pt_operator")).c_str(),
              r.jetThreshold);
}
/* The photon and jet acceptance and the recoil angle, from the metadata. */
std::string acceptanceLine(const Result &r) {
  const auto &s = renderSelection(r.package);
  if (s.is_null()) return "";
  const auto &j = s.at("jet");
  const auto &fraction = j.at("delta_phi_pi_fraction");
  const std::string angle = fraction.is_null() ? Form("%.3g", j.at("delta_phi_min_rad").get<double>())
                                               : Form("%d#pi/%d", fraction.at(0).get<int>(), fraction.at(1).get<int>());
  const double photon = s.at("photon").at("abs_eta_max").get<double>(), jet = j.at("abs_eta_max").get<double>();
  const std::string eta = photon == jet ? Form("#left|#eta^{#gamma,jet}#right| < %g", photon)
                                        : Form("#left|#eta^{#gamma}#right| < %g, #left|#eta^{jet}#right| < %g", photon, jet);
  return eta + ", #Delta#phi_{#gamma jet} " + cutOperator(j.at("delta_phi_operator")) + " " + angle;
}
std::vector<std::string> xjLines(TFile &f, const Result &r) { return {systemLine(r), photonLine(f, r), jetLine(r)}; }

/* ---- text, header and legend ------------------------------------------------------------------------------------- */
double padWidthPx() { return gPad->GetWw() * gPad->GetAbsWNDC(); }
double padHeightPx() { return gPad->GetWh() * gPad->GetAbsHNDC(); }
/* The width of text in the current pad, in NDC. */
double textWidth(const std::string &text, double px, int font = 43) {
  TLatex probe(0.5, 0.5, text.c_str());
  probe.SetNDC();
  probe.SetTextFont(font);
  probe.SetTextSize(px);
  UInt_t w = 0, h = 0;
  probe.GetBoundingBox(w, h);
  return w / padWidthPx();
}
/* Producer text in TLatex grammar: relations and the vertex symbol as the IAN writes them. */
std::string latexCuts(std::string s) {
  const std::pair<std::string, std::string> words[] = {{">=", "#geq"}, {"<=", "#leq"}, {"vtx", "#left|v_{z}#right|"}};
  for (const auto &[from, to] : words)
    for (size_t at = s.find(from); at != std::string::npos; at = s.find(from, at + to.size())) s.replace(at, from.size(), to);
  return s;
}
std::vector<std::string> nonEmpty(const std::vector<std::string> &lines) {
  std::vector<std::string> out;
  for (const auto &l : lines) if (!l.empty()) out.push_back(l);
  return out;
}
/* The note under the header for inputs that are not a complete physics population. */
std::string inputNote() {
  return gLegacy ? "Reference inputs (development)" : gPopulation == "PARTIAL" ? "Partial input population" :
         gPopulation == "UNKNOWN" ? "Input coverage unverified" : "";
}
/* The header: the sPHENIX mark at (x, y), the lines below it and (note) the input note. Returns the lowest baseline
   (NDC); draw = false only measures. */
double header(double x, double y, const std::vector<std::string> &lines, bool draw = true, bool note = true) {
  const double h = padHeightPx();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(43);
  t.SetTextSize(kMarkPx);
  if (draw) { t.DrawLatex(x, y, sPhenixInternalLabel()); gVisibleText.push_back(sPhenixInternalLabel()); }
  t.SetTextSize(kLinePx);
  double last = y - 0.15 * kStepPx / h;
  for (const auto &line : nonEmpty(lines)) {
    last -= kStepPx / h;
    if (draw) { t.DrawLatex(x, last, line.c_str()); gVisibleText.push_back(line); }
  }
  if (const std::string text = inputNote(); note && !text.empty()) {
    last -= kStepPx / h;
    t.SetTextColor(gGray);
    if (draw) { t.DrawLatex(x, last, text.c_str()); gVisibleText.push_back(text); }
  }
  return last;
}
double headerWidth(const std::vector<std::string> &lines) {
  double w = textWidth(sPhenixInternalLabel(), kMarkPx);
  for (const auto &l : nonEmpty(lines)) w = std::max(w, textWidth(l, kLinePx));
  return w;
}
TLegend *ianLegend(double x1, double y1, double x2, double y2, double symbol) {
  auto *l = new TLegend(x1, y1, x2, y2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(43);
  l->SetTextSize(kLinePx);
  l->SetMargin(symbol);
  return l;
}

/* ---- the frame: header top left, legend top right, a y range opened for both ----------------------------------- */
struct Extent { double x0, x1, y; };  /* a drawn element: its horizontal extent and its top, in user coordinates */
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

/* A single-pad frame in the IAN's layout and its empty legend: the header at the top left, the legend at the top right
   (under the header when the two do not fit side by side), and a y range opened until no drawn element reaches either
   block. low: the frame's lower edge (on a log scale, a third of the smallest positive top when low <= 0). The caller
   adds the legend entries in `labels` order and draws the legend. */
TLegend *ianFrame(double x0, double x1, double low, const std::vector<Extent> &drawn, bool logy, const char *xtitle,
                  const char *ytitle, const std::vector<std::string> &lines, const std::vector<std::string> &labels) {
  gPad->SetLogy(logy);
  const double L = gPad->GetLeftMargin(), R = gPad->GetRightMargin(), T = gPad->GetTopMargin(), B = gPad->GetBottomMargin();
  const double plot = 1 - T - B, h = padHeightPx(), w = padWidthPx();
  const double hx = L + 18 / w, top = 1 - T - (14 + kMarkPx) / h;
  const double hw = headerWidth(lines), hBottom = header(hx, top, lines, false, false) - 10 / h;
  const std::string note = inputNote();
  double lw = 0;
  for (const auto &l : labels) lw = std::max(lw, textWidth(l, kLinePx));
  const double symbol = 46 / w, width = symbol + lw + 10 / w, height = labels.size() * kStepPx / h;
  double lx2 = 1 - R - 14 / w, lx1 = lx2 - width, ly2 = 1 - T - 10 / h, ly1 = ly2 - height;
  if (lx1 < hx + hw + 24 / w) {  /* no room side by side: the legend goes under the header */
    lx1 = hx, lx2 = hx + width, ly2 = hBottom - 4 / h, ly1 = ly2 - height;
  }
  /* the input note: under the header, and under the legend too when it would run into it */
  const double noteWidth = note.empty() ? 0.0 : textWidth(note, kLinePx);
  double noteY = hBottom - 0.75 * kStepPx / h;
  if (!note.empty() && hx + noteWidth > lx1 - 6 / w && noteY > ly1 - 0.75 * kStepPx / h) noteY = ly1 - 0.75 * kStepPx / h;
  if (std::min(ly1, note.empty() ? 1.0 : noteY) <= B + 0.05) throw std::runtime_error("figure annotations do not fit the canvas");

  double highest = -std::numeric_limits<double>::infinity(), positive = std::numeric_limits<double>::infinity();
  for (const auto &e : drawn)
    if (std::isfinite(e.y)) { highest = std::max(highest, e.y); if (e.y > 0) positive = std::min(positive, e.y); }
  if (!std::isfinite(highest)) highest = 1;
  if (logy && !(low > 0)) low = std::isfinite(positive) ? positive / 3 : 1e-3;
  const auto ndc = [&](double x) { return L + (x - x0) / (x1 - x0) * (1 - L - R); };
  const auto frameTop = [&](double y, double fraction) {  /* the top that puts y at `fraction` of the frame height */
    fraction = std::max(fraction, 0.1);
    return logy ? std::exp(std::log(low) + (std::log(std::max(y, low)) - std::log(low)) / fraction) : low + (y - low) / fraction;
  };
  double high = frameTop(highest, 0.95);
  std::vector<std::array<double, 3>> blocks = {{hx - 6 / w, hx + hw + 6 / w, hBottom - 6 / h}, {lx1 - 6 / w, lx2 + 6 / w, ly1 - 8 / h}};
  if (!note.empty()) blocks.push_back({hx - 6 / w, hx + noteWidth + 6 / w, noteY - 8 / h});
  for (const auto &e : drawn) {
    if (!std::isfinite(e.y)) continue;
    const double a = ndc(e.x0), b = ndc(e.x1);
    for (const auto &[ka, kb, bottom] : blocks)
      if (b > ka && a < kb) high = std::max(high, frameTop(e.y, (bottom - B) / plot));
  }
  if (!(high > low)) high = low + 1;
  gFrame = gPad->DrawFrame(x0, low, x1, high);
  gFrameLow = low, gFrameHigh = high;
  gFrame->GetXaxis()->SetTitle(xtitle);
  gFrame->GetYaxis()->SetTitle(ytitle);
  gFrame->GetYaxis()->SetMaxDigits(4);  /* five-digit counts take a x10^n exponent instead of crowding the title */
  header(hx, top, lines, true, false);
  if (!note.empty()) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(43);
    t.SetTextSize(kLinePx);
    t.SetTextColor(gGray);
    t.DrawLatex(hx, noteY, note.c_str());
    gVisibleText.push_back(note);
  }
  return ianLegend(lx1, ly1, lx2, ly2, symbol / width);
}

/* ---- ratio pads ---------------------------------------------------------------------------------------------------- */
/* A figure with a ratio panel below: the upper (70%) and lower (30%) pads. */
std::pair<TPad *, TPad *> ratioPads(TCanvas &cv) {
  cv.cd();
  auto *top = new TPad(fresh().c_str(), "", 0, 0.3, 1, 1), *bottom = new TPad(fresh().c_str(), "", 0, 0, 1, 0.3);
  top->SetBottomMargin(0.025);
  bottom->SetTopMargin(0.035);
  bottom->SetBottomMargin(0.44);
  top->Draw();
  bottom->Draw();
  return {top, bottom};
}
/* Axis text in pixels, so the upper and the lower pad of a ratio figure read alike. */
void pixelAxes(TH1 *frame, double label = 24, double title = 26) {
  for (TAxis *a : {frame->GetXaxis(), frame->GetYaxis()})
    a->SetLabelFont(43), a->SetLabelSize(label), a->SetTitleFont(43), a->SetTitleSize(title);
  frame->GetYaxis()->SetTitleOffset(2.0);
  frame->GetXaxis()->SetTitleOffset(1.6);
}
/* The upper pad of a ratio figure: no x labels; the lower pad carries them. */
void upperAxes(TH1 *frame) {
  pixelAxes(frame);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
}
/* The lower frame of a ratio figure over [y0, y1]. */
TH1F *ratioFrame(double x0, double x1, double y0, double y1, const char *xtitle, const char *ytitle) {
  TH1F *f = gPad->DrawFrame(x0, y0, x1, y1);
  f->GetXaxis()->SetTitle(xtitle);
  f->GetYaxis()->SetTitle(ytitle);
  pixelAxes(f);
  f->GetYaxis()->SetNdivisions(503);
  return f;
}
/* The lower frame of a ratio figure, symmetric about 1 and wide enough for every ratio with its error, rounded to a
   readable step; a dashed line marks 1. */
TH1F *ratioFrame(double x0, double x1, const std::vector<std::pair<double, double>> &ratios, const char *xtitle, const char *ytitle) {
  double half = 0.1;
  for (const auto &[v, e] : ratios) if (std::isfinite(v)) half = std::max(half, std::fabs(v - 1) + (std::isfinite(e) ? e : 0));
  const double steps[] = {0.1, 0.2, 0.25, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0};
  for (double s : steps) if (half * 1.15 <= s) { half = s; break; }
  TH1F *f = ratioFrame(x0, x1, std::max(0.0, 1 - half), 1 + half, xtitle, ytitle);
  TLine one;
  one.SetLineStyle(2);
  one.SetLineWidth(1);
  one.DrawLine(x0, 1, x1, 1);
  return f;
}
/* A photon axis with few bins is labelled at its bin edges, as in the IAN (17, 23, 33). */
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
  const double plot = (1 - gPad->GetTopMargin() - gPad->GetBottomMargin()) * h, tick = 14 / plot;
  TLine l;
  l.SetLineWidth(1);
  TLatex t;
  t.SetTextFont(43);
  t.SetTextSize(24);
  t.SetTextAlign(23);
  for (double e : edges) {
    if (e < a->GetXmin() - 1e-9 || e > a->GetXmax() + 1e-9) continue;
    l.DrawLine(e, at(0), e, at(tick));
    l.DrawLine(e, at(1 - tick), e, at(1));
    if (labels) {
      t.SetNDC(false);
      t.DrawLatex(e, logy ? std::exp(std::log(y0) - 10 / plot * (std::log(y1) - std::log(y0))) : y0 - 10 / plot * (y1 - y0), Form("%g", e));
      gVisibleText.push_back(Form("%g", e));
    }
  }
}

/* ---- saving --------------------------------------------------------------------------------------------------------- */
void canvasText(TPad &pad, std::vector<std::string> &text) {
  TIter next(pad.GetListOfPrimitives());
  while (auto *object = next()) {
    if (auto *child = dynamic_cast<TPad *>(object)) canvasText(*child, text);
    else if (auto *latex = dynamic_cast<TLatex *>(object)) text.push_back(latex->GetTitle());
    else if (auto *legend = dynamic_cast<TLegend *>(object)) {
      TIter item(legend->GetListOfPrimitives());
      while (auto *entry = dynamic_cast<TLegendEntry *>(item())) text.push_back(entry->GetLabel());
    } else if (auto *hist = dynamic_cast<TH1 *>(object)) {
      for (auto *axis : {hist->GetXaxis(), hist->GetYaxis(), hist->GetZaxis()}) if (*axis->GetTitle()) text.push_back(axis->GetTitle());
    }
  }
}

/* Fit the actual ROOT/TLatex glyphs inside each legend, including nested pads.
   Keep the configured font and maximum size; longer labels may only shrink. */
void fitLegends(TPad &pad) {
  pad.cd();
  TIter next(pad.GetListOfPrimitives());
  while (auto *object = next()) {
    if (auto *child = dynamic_cast<TPad *>(object)) { fitLegends(*child); pad.cd(); }
    else if (auto *legend = dynamic_cast<TLegend *>(object)) {
      const double size=legend->GetTextSize();
      if (!(size>0)) continue;  /* ROOT's automatic sizing already fits its box. */
      const int columns=std::max(1,legend->GetNColumns());
      const int rows=(legend->GetListOfPrimitives()->GetSize()+columns-1)/columns;
      if (!rows) continue;
      const double width=pad.GetWw()*pad.GetAbsWNDC()*(legend->GetX2NDC()-legend->GetX1NDC())*
                         (1-legend->GetMargin())/columns-6;
      const double height=pad.GetWh()*pad.GetAbsHNDC()*(legend->GetY2NDC()-legend->GetY1NDC())/rows*0.85;
      if (!(width>0 && height>0)) throw std::runtime_error("legend has no usable text area");
      double scale=1;
      TIter item(legend->GetListOfPrimitives());
      while (auto *entry=dynamic_cast<TLegendEntry *>(item())) {
        TLatex probe(0,0,entry->GetLabel()); probe.SetNDC();
        probe.SetTextFont(entry->GetTextFont()>0 ? entry->GetTextFont() : legend->GetTextFont());
        probe.SetTextSize(entry->GetTextSize()>0 ? entry->GetTextSize() : size);
        UInt_t w=0,h=0; probe.GetBoundingBox(w,h);
        if (w) scale=std::min(scale,width/w);
        if (h) scale=std::min(scale,height/h);
      }
      if (scale<1) {
        legend->SetTextSize(size*scale);
        TIter entries(legend->GetListOfPrimitives());
        while (auto *entry=dynamic_cast<TLegendEntry *>(entries()))
          if (entry->GetTextSize()>0) entry->SetTextSize(entry->GetTextSize()*scale);
      }
    }
  }
  pad.Modified();
}

void save(TCanvas &c, const std::string &relative) {
  const std::string path = gFig + "/" + relative;
  gSystem->mkdir(gSystem->GetDirName(path.c_str()).Data(), true);
  c.Update(); fitLegends(c); c.cd(); c.Update();
  c.SaveAs(path.c_str());
  std::ifstream saved(path, std::ios::binary | std::ios::ate);
  if (!saved || saved.tellg() <= 0) throw std::runtime_error("figure was not saved: " + path);
  bool compiled = !gSemanticPackages.empty() && !gLegacy;
  PlotJson semantics = PlotJson::object();
  for (const auto &key : gSemanticPackages) {
    const bool qa = key.rfind("qa/",0)==0;
    const auto pkg = key.substr(key.find('/')+1);
    const auto &metadata = qa ? gQAMetadata : gResultMetadata;
    const bool valid = !metadata.is_null() && metadata.at("packages").contains(pkg) && !metadata.at("packages").at(pkg).at("selection").is_null();
    compiled = compiled && valid;
    if (valid) semantics[key] = metadata.at("packages").at(pkg);
  }
  std::vector<std::string> visible; canvasText(c, visible);
  for (const auto &t : gVisibleText) if (std::find(visible.begin(),visible.end(),t)==visible.end()) visible.push_back(t);
  photonjet_render::write(path+".json", {{"schema", "PhotonJetRenderReceiptV1"}, {"image_sha256", photonjet_render::sha256(path, true)},
    {"inputs", gRenderInputs}, {"style_and_renderer", gStylePins}, {"visible_text", visible}, {"selection", semantics}, {"display_projections", gProjectionSteps},
    {"annotation_status", compiled ? "COMPILED_METADATA" : "LEGACY_UNVERIFIED_ANNOTATION"},
    {"input_coverage", gPopulation}, {"layout_review", "REPRESENTATIVE_VISUAL_QA_REQUIRED"}});
  gVisibleText.clear(); gSemanticPackages.clear(); gProjectionSteps = PlotJson::array();
  gWritten.push_back(relative);
}

/* ---- series ----------------------------------------------------------------------------------------------------------- */
template <class T> void style(T *h, int color, int marker, double size = 1.2) {
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(size);
  h->SetLineWidth(2);
}

/* pairs summed over photon bins per unit xJ, divided by `photons` when given */
TH1D *xjProjection(const TH2 *pairs, double photons, const char *name) {
  auto *h = pairs->ProjectionY(Form("%s%s", name, fresh().c_str()), 1, pairs->GetNbinsX(), "e");
  h->SetDirectory(nullptr);
  h->Scale(photons > 0 ? 1.0 / photons : 1.0, "width");
  return h;
}

/* The bins of a distribution inside [low, high], as points; widths: horizontal bars over each bin. */
TGraphErrors *points(const TH1 *h, double low, double high, int color, int marker, bool widths = true, double size = 1.2) {
  auto *g = new TGraphErrors;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double a = h->GetXaxis()->GetBinLowEdge(i), b = h->GetXaxis()->GetBinUpEdge(i);
    if (a < low - 1e-9 || b > high + 1e-9) continue;
    const int n = g->GetN();
    g->SetPoint(n, 0.5 * (a + b), h->GetBinContent(i));
    g->SetPointError(n, widths ? 0.5 * (b - a) : 0.0, h->GetBinError(i));
  }
  style(g, color, marker, size);
  return g;
}

/* A background estimate as the IAN draws it: a step line over the bins inside [low, high], filled when `filled`. */
TH1D *steps(const TH1 *h, double low, double high, int color, bool filled) {
  auto *s = static_cast<TH1D *>(h->Clone(fresh().c_str()));
  s->SetDirectory(nullptr);
  for (int i = 0; i <= s->GetNbinsX() + 1; ++i) {
    const double a = s->GetXaxis()->GetBinLowEdge(i), b = s->GetXaxis()->GetBinUpEdge(i);
    if (i == 0 || i > s->GetNbinsX() || a < low - 1e-9 || b > high + 1e-9) s->SetBinContent(i, 0), s->SetBinError(i, 0);
  }
  s->SetLineColor(color);
  s->SetLineWidth(2);
  if (filled) s->SetFillColorAlpha(color, 0.18);
  return s;
}

/* The ratio a / b bin by bin (errors of both), as points. */
TGraphErrors *ratioPoints(const TH1 *a, const TH1 *b, double low, double high, int color, int marker,
                          std::vector<std::pair<double, double>> &values) {
  auto *r = static_cast<TH1 *>(a->Clone(fresh().c_str()));
  r->SetDirectory(nullptr);
  r->Divide(b);
  auto *g = points(r, low, high, color, marker);
  for (int i = 0; i < g->GetN(); ++i) values.push_back({g->GetY()[i], g->GetEY()[i]});
  return g;
}

/* The systematic band of a result as boxes over the reported bins inside [low, high]; null without systematics. */
TGraphAsymmErrors *systematicBoxes(const std::string &key, double low, double high, int color) {
  if (!gSystematics) return nullptr;
  auto *c = dynamic_cast<TH1 *>(gSystematics->Get((key + "/xj_per_photon").c_str()));
  auto *up = dynamic_cast<TH1 *>(gSystematics->Get((key + "/xj_systematic_up").c_str()));
  auto *down = dynamic_cast<TH1 *>(gSystematics->Get((key + "/xj_systematic_down").c_str()));
  if (!c || !up || !down) return nullptr;
  auto *g = new TGraphAsymmErrors;
  for (int i = 1; i <= c->GetNbinsX(); ++i) {
    const double a = c->GetXaxis()->GetBinLowEdge(i), b = c->GetXaxis()->GetBinUpEdge(i);
    if (a < low - 1e-9 || b > high + 1e-9) continue;
    const int n = g->GetN();
    g->SetPoint(n, 0.5 * (a + b), c->GetBinContent(i));
    g->SetPointError(n, 0.5 * (b - a), 0.5 * (b - a), -down->GetBinContent(i), up->GetBinContent(i));
  }
  g->SetFillColorAlpha(color, 0.3);
  g->SetLineColor(color);
  return g;
}

/* ---- xJ figures ---------------------------------------------------------------------------------------------------- */
const char *kPairs = "dN_{pair}/dx_{J#gamma}", *kPerPhoton = "(1/N_{#gamma}) dN_{pair}/dx_{J#gamma}";
double xjEnd(double high) { return std::ceil(high * 10 - 1e-9) / 10; }

/* xJ ticks every 0.3 when the axis end is a multiple of 0.3, as in the IAN (0, 0.3, ..., 1.5); otherwise ROOT's. */
void xjTicks(TH1 *frame, double high) {
  const double end = xjEnd(high), steps = end / 0.3;
  if (std::fabs(steps - std::round(steps)) < 1e-6 && steps <= 10) frame->GetXaxis()->SetNdivisions(int(std::round(steps)) + 300, false);
  else frame->GetXaxis()->SetNdivisions(505);
}

/* One entry of an xJ figure: a graph (points) or a histogram (a background estimate), its legend text and style. */
struct Series { TObject *object; std::string text; const char *legend; };

/* An xJ frame in the IAN's layout with its series: histograms under the points, the legend in series order. */
void xjFrame(const std::vector<Series> &series, double high, const char *ytitle, const std::vector<std::string> &lines,
             const std::vector<TGraphAsymmErrors *> &boxes = {}, bool upper = false) {
  std::vector<Extent> drawn;
  std::vector<std::string> labels;
  for (const auto &s : series) {
    if (auto *g = dynamic_cast<TGraph *>(s.object)) extents(drawn, g);
    else if (auto *h = dynamic_cast<TH1 *>(s.object)) extents(drawn, h);
    labels.push_back(s.text);
  }
  for (auto *b : boxes) if (b) extents(drawn, b);
  TLegend *leg = ianFrame(0, xjEnd(high), 0, drawn, false, upper ? "" : "x_{J#gamma}", ytitle, lines, labels);
  xjTicks(gFrame, high);
  if (upper) upperAxes(gFrame);
  for (auto *b : boxes) if (b) b->Draw("2");
  for (const auto &s : series) if (dynamic_cast<TH1 *>(s.object)) s.object->Draw("HIST SAME");
  for (const auto &s : series) if (dynamic_cast<TGraph *>(s.object)) s.object->Draw("P");
  for (const auto &s : series) leg->AddEntry(s.object, s.text.c_str(), s.legend);
  leg->Draw();
}

/* Photon-purity correction of the reconstructed xJ distribution, in pairs per unit xJ: region A, the transferred region-C
   background alpha C and the leakage-corrected distribution, as in the IAN. */
void purityCorrection(TFile &f, const Result &r, const std::string &name) {
  auto *a = dynamic_cast<TH2 *>(get(f, r, "pairs_A"));
  auto *c = dynamic_cast<TH2 *>(get(f, r, "pairs_C_scaled"));
  auto *p = dynamic_cast<TH2 *>(get(f, r, "pairs_purity_corrected"));
  if (!a || !c || !p) return missing(name, "needs pairs_A, pairs_C_scaled and pairs_purity_corrected");
  TCanvas cv("cv", "", 700, 600);
  xjFrame({{steps(xjProjection(c, 0, "hc"), r.low, r.high, gVermilion, false), "#alphaC", "l"},
           {points(xjProjection(a, 0, "ha"), r.low, r.high, kBlack, 20), "Region A", "lp"},
           {points(xjProjection(p, 0, "hp"), r.low, r.high, gBlue, 33, true, 1.6), "Leakage-corrected", "lp"}},
          r.high, kPairs, xjLines(f, r));
  save(cv, name);
}

/* Combinatorial-jet subtraction in pairs per unit xJ: the scaled combinatorial estimate, the purity-corrected
   distribution and the distribution after subtraction, as in the IAN. */
void combinatoric(TFile &f, const Result &r, const std::string &name) {
  auto *p = dynamic_cast<TH2 *>(get(f, r, "pairs_purity_corrected"));
  auto *k = dynamic_cast<TH2 *>(get(f, r, "pairs_combinatoric"));
  auto *s = dynamic_cast<TH2 *>(get(f, r, "pairs_corrected"));
  if (!p || !k || !s) return missing(name, "needs pairs_purity_corrected, pairs_combinatoric and pairs_corrected");
  const std::string estimate = r.combinatoric.rfind("data_dphi", 0) == 0 ? "Scaled combinatorial (#Delta#phi sideband)"
                               : k->Integral() > 0 ? "Scaled combinatorial (simulation)" : "No combinatorial subtraction";
  TCanvas cv("cv", "", 700, 600);
  xjFrame({{steps(xjProjection(k, 0, "hk"), r.low, r.high, gVermilion, true), estimate, "f"},
           {points(xjProjection(p, 0, "hp"), r.low, r.high, gBlue, 27, true, 1.6), "Purity-corrected", "lp"},
           {points(xjProjection(s, 0, "hs"), r.low, r.high, kBlack, 20), "After subtraction", "lp"}},
          r.high, kPairs, xjLines(f, r));
  save(cv, name);
}

/* Per-photon xJ before (purity- and combinatoric-corrected) and after unfolding, and their ratio, as in the IAN. */
void beforeAfter(TFile &f, const Result &r, const std::string &name) {
  auto *b = get(f, r, "xj_per_photon_reco");
  auto *a = get(f, r, "xj_per_photon");
  if (!a || !b) return missing(name, "needs xj_per_photon and xj_per_photon_reco");
  TCanvas cv("cv", "", 700, 760);
  auto [top, bottom] = ratioPads(cv);
  top->cd();
  xjFrame({{points(b, r.low, r.high, kBlack, 24), "Corrected reconstructed", "lp"},
           {points(a, r.low, r.high, gBlue, 20), "Unfolded", "lp"}}, r.high, kPerPhoton, xjLines(f, r), {}, true);
  bottom->cd();
  std::vector<std::pair<double, double>> values;
  TGraphErrors *ratio = ratioPoints(a, b, r.low, r.high, gBlue, 20, values);
  xjTicks(ratioFrame(0, xjEnd(r.high), values, "x_{J#gamma}", "Unfolded / reco"), r.high);
  ratio->Draw("P");
  save(cv, name);
}

/* The unfolded per-photon xJ of p+p and of the most central Au+Au class over the range both report, as in the IAN. */
void ppVsCentral(TFile &f, const std::vector<Result> &results, const std::string &name) {
  const Result *pp = nullptr, *central = nullptr;
  for (const auto &r : results) {
    if (r.package == "PP") pp = &r;
    else if (r.centralityLow >= 0 && (!central || r.centralityLow < central->centralityLow)) central = &r;
  }
  if (!pp || !central) return missing(name, "needs a p+p and an Au+Au result");
  auto *hp = get(f, *pp, "xj_per_photon"), *ha = get(f, *central, "xj_per_photon");
  if (!hp || !ha) return missing(name, "needs xj_per_photon of both systems");
  const double low = std::max(pp->low, central->low), high = std::min(pp->high, central->high);
  const auto window = [&](const Result &r) {
    const auto e = photonEdges(f, r);
    return e.size() < 2 ? std::string() : Form(", %g-%g GeV photons", e.front(), e.back());
  };
  std::vector<std::string> lines{central->energy.empty() ? pp->energy : central->energy, jetLine(*pp)};
  if (jetLine(*central) != jetLine(*pp)) lines.push_back(jetLine(*central));
  lines.push_back(acceptanceLine(*pp));
  TCanvas cv("cv", "", 700, 600);
  xjFrame({{points(hp, low, high, gBlue, 20), pp->label + window(*pp), "lp"},
           {points(ha, low, high, gRed, 21), central->label + window(*central), "lp"}},
          high, kPerPhoton, lines,
          {systematicBoxes(pp->key, low, high, gBlue), systematicBoxes(central->key, low, high, gRed)});
  save(cv, name);
}

/* The unfolded per-photon xJ of every result on one plot, over the range all of them report, with the systematic
   boxes when a systematics file is given. */
void centralityOverlay(TFile &f, const std::vector<Result> &results, const std::string &name) {
  if (results.empty()) return;
  double low = 0, high = 1e9;
  for (const auto &r : results) { low = std::max(low, r.low); high = std::min(high, r.high); }
  std::vector<Series> series;
  std::vector<TGraphAsymmErrors *> boxes;
  int i = 0;
  for (const auto &r : results)
    if (auto *h = get(f, r, "xj_per_photon")) {
      series.push_back({points(h, low, high, seriesColor(i), kMarkers[i % 8], true, 1.2), r.label, "lp"});
      boxes.push_back(systematicBoxes(r.key, low, high, seriesColor(i)));
      ++i;
    }
  if (series.empty()) return missing(name, "needs xj_per_photon");
  const Result &first = results.front();
  TCanvas cv("cv", "", 700, 600);
  xjFrame(series, high, kPerPhoton, {first.energy, jetLine(first), acceptanceLine(first),
                                     gSystematics ? "Bars statistical, boxes systematic" : "Statistical uncertainties"}, boxes);
  save(cv, name);
}

/* The joint unfolding's prior: the simulated truth, the truth reweighted to the data, and the unfolded data from each. */
void priorReweighting(TFile &f, const Result &r, const std::string &name) {
  auto *simulation = get(f, r, "xj_prior_simulation");
  auto *reweighted = get(f, r, "xj_prior_reweighted");
  auto *nominal = get(f, r, "xj_per_photon");
  auto *fromSimulation = get(f, r, "xj_per_photon_simulation_prior");
  if (!simulation || !reweighted || !nominal || !fromSimulation) return missing(name, "needs the prior histograms");
  TCanvas cv("cv", "", 700, 600);
  xjFrame({{points(simulation, r.low, r.high, gGray, 25), "Simulation truth", "lp"},
           {points(reweighted, r.low, r.high, gVermilion, 26), "Truth reweighted to the data", "lp"},
           {points(fromSimulation, r.low, r.high, gRed, 24), "Unfolded, simulation prior", "lp"},
           {points(nominal, r.low, r.high, gBlue, 20), "Unfolded, reweighted prior", "lp"}},
          r.high, kPerPhoton, xjLines(f, r));
  save(cv, name);
}

/* The unfolded per-photon xJ of p+p (atlas: "pp") or of the most central Au+Au class ("pbpb_0_10") against ATLAS
   (Phys. Lett. B 789 (2019) 167, HEPData ins1694678 Table 1: p+p and 0-10% Pb+Pb at 5.02 TeV, plots/inputs/
   atlas_reference.json), the ATLAS systematic uncertainty shaded. */
void atlasOverlay(TFile &f, const Result &r, const std::string &atlas, const std::string &name) {
  auto *h = get(f, r, "xj_per_photon");
  std::ifstream in(std::string(gSystem->GetDirName(__FILE__).Data()) + "/inputs/atlas_reference.json");
  if (!h || !in) return missing(name, "needs xj_per_photon and the ATLAS reference");
  const nlohmann::json ref = nlohmann::json::parse(in), &sel = ref.at("selection");
  auto *atlasPoints = new TGraphErrors;
  auto *atlasBoxes = new TGraphAsymmErrors;
  for (const auto &b : ref.at("bins")) {
    const double lo = b.at("low"), hi = b.at("high"), y = b.at(atlas).at("value"), stat = b.at(atlas).at("stat"), sys = b.at(atlas).at("sys");
    const int n = atlasPoints->GetN();
    atlasPoints->SetPoint(n, 0.5 * (lo + hi), y), atlasPoints->SetPointError(n, 0.5 * (hi - lo), stat);
    atlasBoxes->SetPoint(n, 0.5 * (lo + hi), y), atlasBoxes->SetPointError(n, 0.5 * (hi - lo), 0.5 * (hi - lo), sys, sys);
  }
  style(atlasPoints, gVermilion, 25);
  atlasBoxes->SetFillColorAlpha(gVermilion, 0.2);
  atlasBoxes->SetLineWidth(0);
  const auto photons = sel.at("photon_pt_gev");
  const std::string system = atlas == "pp" ? "p+p" : Form("Pb+Pb %g-%g%%", sel.at("pbpb_centrality_percent")[0].get<double>(),
                                                           sel.at("pbpb_centrality_percent")[1].get<double>());
  std::vector<TGraphAsymmErrors *> shaded = {atlasBoxes};
  if (auto *ours = systematicBoxes(r.key, r.low, r.high, gBlue)) shaded.push_back(ours);
  TCanvas cv("cv", "", 700, 600);
  xjFrame({{points(h, r.low, r.high, gBlue, 20), systemLine(r), "lp"},
           {atlasPoints, Form("ATLAS %s, %g TeV", system.c_str(), sel.at("sqrt_s_nn_gev").get<double>() / 1000), "lp"}},
          r.high, kPerPhoton,
          {photonLine(f, r), jetLine(r), Form("ATLAS: %g < p_{T}^{#gamma} < %g GeV, p_{T}^{jet} #geq %g GeV", photons[0].get<double>(),
                                              photons[1].get<double>(), sel.at("jet_pt_min_gev").get<double>())},
          shaded);
  save(cv, name);
}

/* ---- responses ------------------------------------------------------------------------------------------------------ */
/* Dashed lines at the result's photon bins over a response map, on both axes, below yTop (the data's upper edge). */
void binLines(const TH2 *m, const std::vector<double> &edges, double yTop) {
  TLine l;
  l.SetLineStyle(2);
  l.SetLineWidth(1);
  l.SetLineColor(kGray + 2);
  const double x0 = m->GetXaxis()->GetXmin(), x1 = m->GetXaxis()->GetXmax(), y0 = m->GetYaxis()->GetXmin(), y1 = yTop;
  for (double e : edges) {
    if (e > x0 + 1e-9 && e < x1 - 1e-9) l.DrawLine(e, y0, e, y1);
    if (e > y0 + 1e-9 && e < y1 - 1e-9) l.DrawLine(x0, e, x1, e);
  }
}
/* A map with a data-free band added above its data, tall enough for the header `lines` at the top left of the current
   pad (the y axis is extended, nothing is hidden). */
TH2D *withHeaderBand(const TH2 *m, const std::vector<std::string> &lines) {
  const double h = padHeightPx(), plot = 1 - gPad->GetTopMargin() - gPad->GetBottomMargin();
  const double top = 1 - gPad->GetTopMargin() - (14 + kMarkPx) / h;
  const double fraction = std::min(0.5, (1 - gPad->GetTopMargin() - (header(0.5, top, lines, false) - 14 / h)) / plot);
  const auto x = edgesOf(m->GetXaxis());
  auto y = edgesOf(m->GetYaxis());
  y.push_back(y.front() + (y.back() - y.front()) / (1 - fraction));
  auto *b = new TH2D(fresh().c_str(), "", x.size() - 1, x.data(), y.size() - 1, y.data());
  b->SetDirectory(nullptr);
  for (int i = 1; i <= m->GetNbinsX(); ++i)
    for (int j = 1; j <= m->GetNbinsY(); ++j) b->SetBinContent(i, j, m->GetBinContent(i, j)), b->SetBinError(i, j, m->GetBinError(i, j));
  return b;
}
/* The header at the top left inside the frame of the current pad. */
void mapHeader(const std::vector<std::string> &lines) {
  header(gPad->GetLeftMargin() + 18 / padWidthPx(), 1 - gPad->GetTopMargin() - (14 + kMarkPx) / padHeightPx(), lines);
}

/* The colour floor of a log map: the decade at or below its smallest positive cell. */
double logFloor(const TH2 *m) {
  double low = std::numeric_limits<double>::infinity();
  for (int b = 0; b < m->GetNcells(); ++b) if (m->GetBinContent(b) > 0) low = std::min(low, m->GetBinContent(b));
  return std::isfinite(low) ? std::pow(10, std::floor(std::log10(low))) : 1e-3;
}

/* Photon response on its native grid, reconstructed on x and truth on y (the house style), each truth bin normalised to
   unit probability, on a log colour scale whose floor follows the data; the result's photon bins are dashed. */
void photonResponse(TFile &f, const Result &r, const std::string &name) {
  auto *grid = dynamic_cast<TH2 *>(get(f, r, "photon_response_grid"));  /* saved as truth (x) by reco (y) */
  if (!grid) return missing(name, "needs photon_response_grid");
  const auto truth = edgesOf(grid->GetXaxis()), reco = edgesOf(grid->GetYaxis());
  auto *m = new TH2D(fresh().c_str(), "", reco.size() - 1, reco.data(), truth.size() - 1, truth.data());
  m->SetDirectory(nullptr);
  for (int i = 1; i <= grid->GetNbinsX(); ++i) {
    double row = 0;
    for (int j = 1; j <= grid->GetNbinsY(); ++j) row += grid->GetBinContent(i, j);
    for (int j = 1; j <= grid->GetNbinsY(); ++j) m->SetBinContent(j, i, row > 0 ? grid->GetBinContent(i, j) / row : 0);
  }
  TCanvas cv("cv", "", 760, 680);
  cv.SetRightMargin(0.19);
  cv.SetLogz();
  const std::vector<std::string> lines{systemLine(r)};
  const double truthTop = m->GetYaxis()->GetXmax();
  m = withHeaderBand(m, lines);
  m->SetMinimum(logFloor(m));
  m->SetMaximum(1);
  m->GetXaxis()->SetTitle("Reconstructed photon p_{T}^{#gamma} [GeV]");
  m->GetYaxis()->SetTitle("Truth photon p_{T}^{#gamma} [GeV]");
  m->GetZaxis()->SetTitle("Migration probability");
  m->GetZaxis()->SetTitleOffset(1.35);
  m->Draw("COLZ");
  binLines(m, photonEdges(f, r), truthTop);
  mapHeader(lines);
  save(cv, name);
}

/* Joint (photon pT, xJ) response of the measurement window: truth cell (x) against reco cell (y), as in the IAN. */
void jointResponse(TFile &f, const Result &r, const std::string &name) {
  auto *m = dynamic_cast<TH2 *>(get(f, r, "pair_response"));
  if (!m) return missing(name, "needs pair_response");
  TCanvas cv("cv", "", 760, 680);
  cv.SetRightMargin(0.19);
  cv.SetLogz();
  m = withHeaderBand(m, {systemLine(r)});
  m->SetMinimum(logFloor(m));
  m->GetXaxis()->SetTitle("Global bin (truth: p_{T}^{#gamma}, x_{J#gamma})");
  m->GetYaxis()->SetTitle("Global bin (reco: p_{T}^{#gamma}, x_{J#gamma})");
  m->GetZaxis()->SetTitle("Weighted matched pairs");
  m->GetZaxis()->SetTitleOffset(1.35);
  m->Draw("COLZ");
  mapHeader({systemLine(r)});
  save(cv, name);
}

/* ---- unfolding stability and closure --------------------------------------------------------------------------- */
/* The frozen iteration rule of one unfolding, as in the IAN: the data and response statistical terms, the change from
   k - 1 and their quadrature Q_k against k, with the selected k marked. Drawn from the saved rule table; a result that
   used configured iterations has none and is listed as not drawn. */
void iterationRule(TFile &f, const Result &r, const char *table, const std::vector<std::string> &lines, const std::string &name) {
  auto *h = dynamic_cast<TH2 *>(get(f, r, table));
  if (!h) return missing(name, "needs the frozen iteration rule (this result used configured iterations)");
  const auto row = [&](const char *label) { return h->GetYaxis()->FindBin(label); };
  const int q = row("Q"), data = row("S_data"), response = row("S_response"), change = row("D_k"), chosen = row("selected");
  int selected = 0, k0 = 0, k1 = 0;
  for (int k = 1; k <= h->GetNbinsX(); ++k) {
    if (!(h->GetBinContent(k, q) > 0)) continue;
    if (!k0) k0 = k;
    k1 = k;
    if (h->GetBinContent(k, chosen) > 0.5) selected = k;
  }
  if (!k0) return missing(name, "the iteration rule table is empty");
  const struct { int row; const char *text; int color, marker; } terms[] = {
      {data, "Data statistics", gBlue, 24}, {response, "Response statistics", gPurple, 26},
      {change, "Iteration change, k vs. k#minus1", gOrange, 25}, {q, "Quadrature Q_{k}", gGreen, 20}};
  std::vector<Series> series;
  std::vector<Extent> drawn;
  double low = std::numeric_limits<double>::infinity(), highest = 0;
  for (const auto &t : terms) {
    auto *g = new TGraph;
    for (int k = k0; k <= k1; ++k) {
      const double v = h->GetBinContent(k, t.row);
      if (std::isfinite(v) && v >= 0) g->SetPoint(g->GetN(), k, v);
      if (v > 0) low = std::min(low, v), highest = std::max(highest, v);
    }
    style(g, t.color, t.marker, 1.3);
    extents(drawn, g);
    series.push_back({g, t.text, "lp"});
  }
  const bool logy = highest > 30 * low;
  std::vector<std::string> labels;
  for (const auto &s : series) labels.push_back(s.text);
  TCanvas cv("cv", "", 700, 600);
  TLegend *leg = ianFrame(k0 - 0.4, k1 + 0.4, 0, drawn, logy, "Bayesian iterations, k", "Relative uncertainty / change",
                          lines, labels);
  gFrame->GetXaxis()->SetNdivisions(std::min(k1 - k0 + 1, 10), true);
  if (selected) {
    TLine l;
    l.SetLineStyle(2);
    l.SetLineWidth(1);
    l.DrawLine(selected, gFrameLow, selected, gFrameHigh);
    TLatex t;
    t.SetTextFont(43);
    t.SetTextSize(kLinePx);
    t.DrawLatex(selected + 0.08, logy ? gFrameLow * std::pow(gFrameHigh / gFrameLow, 0.05) : gFrameLow + 0.05 * (gFrameHigh - gFrameLow),
                Form("Selected k_{Bayes} = %d", selected));
    gVisibleText.push_back(Form("Selected k_{Bayes} = %d", selected));
  }
  for (const auto &s : series) s.object->Draw("LP");
  for (const auto &s : series) leg->AddEntry(s.object, s.text.c_str(), s.legend);
  leg->Draw();
  save(cv, name);
}

/* Simulation closure of the photon-yield unfolding, as in the IAN: truth and the unfolded simulated reconstruction on a
   log scale, and unfolded / truth. */
void photonClosure(TFile &f, const Result &r, const std::string &name) {
  auto *u = get(f, r, "closure_photons");
  auto *t = get(f, r, "photon_truth");
  if (!u || !t) return missing(name, "needs closure_photons and photon_truth");
  const auto edges = edgesOf(t->GetXaxis());
  const double x0 = edges.front(), x1 = edges.back();
  TCanvas cv("cv", "", 700, 760);
  auto [top, bottom] = ratioPads(cv);
  top->cd();
  TGraphErrors *truth = points(t, x0, x1, gBlue, 24), *unfolded = points(u, x0, x1, gRed, 20);
  std::vector<Extent> drawn;
  extents(drawn, truth), extents(drawn, unfolded);
  TLegend *leg = ianFrame(x0, x1, 0, drawn, true, "", "Weighted photons", {systemLine(r), "Simulation, same-population closure"},
                          {"Truth", "Unfolded"});
  upperAxes(gFrame);
  edgeTicks(gFrame, edges, false);
  truth->Draw("P"), unfolded->Draw("P");
  leg->AddEntry(truth, "Truth", "lp"), leg->AddEntry(unfolded, "Unfolded", "lp");
  leg->Draw();
  bottom->cd();
  std::vector<std::pair<double, double>> values;
  TGraphErrors *ratio = ratioPoints(u, t, x0, x1, gRed, 20, values);
  edgeTicks(ratioFrame(x0, x1, values, "p_{T}^{#gamma, truth} [GeV]", "Unfolded / truth"), edges, true);
  ratio->Draw("P");
  save(cv, name);
}

/* Purity-corrected reconstructed photon yields and the unfolded yields, and unfolded / reco, as in the IAN. */
void photonYield(TFile &f, const Result &r, const std::string &name) {
  auto *c = get(f, r, "photons_corrected");
  auto *u = get(f, r, "photons_unfolded");
  if (!c || !u) return missing(name, "needs photons_corrected and photons_unfolded");
  const auto edges = edgesOf(c->GetXaxis());
  const double x0 = edges.front(), x1 = edges.back();
  TCanvas cv("cv", "", 700, 760);
  auto [top, bottom] = ratioPads(cv);
  top->cd();
  TGraphErrors *reco = points(c, x0, x1, gRed, 20), *unfolded = points(u, x0, x1, kBlack, 24);
  std::vector<Extent> drawn;
  extents(drawn, reco), extents(drawn, unfolded);
  TLegend *leg = ianFrame(x0, x1, 0, drawn, true, "", "Photon yield", {systemLine(r), Form("k_{Bayes} = %d", r.photonIterations)},
                          {"Corrected reconstructed", "Unfolded"});
  upperAxes(gFrame);
  edgeTicks(gFrame, edges, false);
  reco->Draw("P"), unfolded->Draw("P");
  leg->AddEntry(reco, "Corrected reconstructed", "lp"), leg->AddEntry(unfolded, "Unfolded", "lp");
  leg->Draw();
  bottom->cd();
  std::vector<std::pair<double, double>> values;
  TGraphErrors *ratio = ratioPoints(u, c, x0, x1, gBlue, 24, values);
  edgeTicks(ratioFrame(x0, x1, values, "p_{T}^{#gamma} [GeV]", "Unfolded / reco"), edges, true);
  ratio->Draw("P");
  save(cv, name);
}

}  /* namespace */

#include "QAFigures.C"

/* results: the results file (or ""), fig: the output folder, systematics and qa: optional systematics and QA files.
   Draw saved quantities with their verified metadata; FIGURES.txt records written and unavailable figures. */
void IANFigures(const char *results, const char *fig, const char *systematics = "", const char *qa = "") {
  gWritten.clear(); gNotDrawn.clear(); gSystematics.reset(); gQA.reset(); gLegacy = false;
  gResultMetadata = gQAMetadata = PlotJson(); gRenderInputs = gStylePins = PlotJson::object(); gVisibleText.clear();
  gSemanticPackages.clear(); gProjectionSteps = PlotJson::array();
  const auto source = std::filesystem::path(__FILE__).parent_path();
  for (const auto &name : {"IANFigures.C", "QAFigures.C", "sPhenixStyle.C", "sPhenixStyle.h"})
    gStylePins[name] = photonjet_render::sha256((source/name).string(), true);
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  ianPalette();
  gFig = fig;
  auto open = [](const char *path, std::unique_ptr<TFile> &file) {
    file.reset(TFile::Open(path, "READ"));
    if (file && !file->IsZombie()) return true;
    throw std::runtime_error(std::string("IANFigures: cannot open ") + path);
  };
  if (systematics && *systematics && !open(systematics, gSystematics)) return;
  if (qa && *qa && !open(qa, gQA)) return;
  std::unique_ptr<TFile> f;
  if (results && *results && !open(results, f)) return;
  if (f) gResultMetadata = loadRenderMetadata(*f, results);
  if (gQA) gQAMetadata = loadRenderMetadata(*gQA, qa);
  if (gSystematics) gRenderInputs[systematics] = photonjet_render::sha256(systematics, true);
  if (f && gSystematics) {
    auto *nominal = dynamic_cast<TObjString *>(gSystematics->Get("nominal_input_sha256"));
    if (nominal && std::string(nominal->GetString().Data())!=gRenderInputs.at(results).get<std::string>())
      throw std::runtime_error("systematic bands belong to a different nominal result");
    if (!nominal) gLegacy = true;
  }
  if (!gResultMetadata.is_null() && !gQAMetadata.is_null())
    for (const auto &[package, metadata] : gResultMetadata.at("packages").items())
      if (gQAMetadata.at("packages").contains(package) && metadata != gQAMetadata.at("packages").at(package))
        throw std::runtime_error("QA and measurement refer to different sample populations or selections: "+package);
  for (const auto &metadata : {gResultMetadata,gQAMetadata}) if (!metadata.is_null())
    for (const auto &package : metadata.at("packages")) if (package.at("selection").is_null()) gLegacy = true;
  if ((f && gResultMetadata.is_null()) || (gQA && gQAMetadata.is_null())) gLegacy = true;
  if (auto *name = f ? dynamic_cast<TObjString *>(f->Get("analysis")) : nullptr) gLegacy = gLegacy || std::string(name->GetString().Data()).find("legacy") == 0;
  if (auto *coverage = f ? dynamic_cast<TObjString *>(f->Get("input_coverage")) : gQA ? dynamic_cast<TObjString *>(gQA->Get("input_coverage")) : nullptr) {
    gPopulation = nlohmann::json::parse(coverage->GetString().Data()).value("status", "UNKNOWN");
    if (gPopulation == "REFERENCE_INPUTS") gLegacy = true;
  } else gPopulation = "UNKNOWN";
  const std::vector<Result> all = f ? resultsOf(*f) : std::vector<Result>();
  if (f && all.empty()) throw std::runtime_error("IANFigures: no results in the supplied file");
  gSemanticPackages.clear();
  for (const auto &r : all) {
    const bool pp = r.tag == "pp";
    purityCorrection(*f, r, "purityAndCombSub/" + r.tag + "_xj_purity_correction.png");
    combinatoric(*f, r, "purityAndCombSub/" + r.tag + (pp ? "_xj_combinatoric_crosscheck.png" : "_xj_combinatoric_subtraction.png"));
    photonResponse(*f, r, "Unfolding/" + r.longTag + "_photon_response.png");
    jointResponse(*f, r, "Unfolding/" + r.longTag + "_joint_response.png");
    iterationRule(*f, r, "iteration_rule_photons", {systemLine(r), photonLine(*f, r)},
                  "Unfolding/" + r.tag + (pp ? "_photon_yield_iteration_stability_publication.png" : "_photon_yield_iteration_stability.png"));
    iterationRule(*f, r, "iteration_rule_pairs", {systemLine(r), photonLine(*f, r), "Joint (p_{T}^{#gamma}, x_{J#gamma}) unfolding",
                                                  Form("%g #leq x_{J#gamma} #leq %g", r.low, r.high)},
                  "Unfolding/" + r.tag + "_joint_iteration_stability.png");
    photonClosure(*f, r, "Unfolding/" + r.tag + "_photon_closure_overlay.png");
    photonYield(*f, r, "Unfolding/" + r.tag + "_photon_yield_unfolding_nominal_2bin_k" + std::to_string(r.photonIterations) + ".png");
    priorReweighting(*f, r, "Unfolding/" + r.tag + "_xj_prior_reweighting.png");
    beforeAfter(*f, r, "Unfolding/" + r.tag + "_xj_before_after_unfolding.png");
    if (pp || r.centralityLow == 0) atlasOverlay(*f, r, pp ? "pp" : "pbpb_0_10", "LHCComparison/" + r.tag + "_xj_atlas_overlay.png");
  }
  if (f) {
    ppVsCentral(*f, all, "Unfolding/pp_vs_auau_xj_unfolded.png");
    centralityOverlay(*f, all, "FinalResults/pp_vs_auau_centrality_xj_unfolded.png");
  }
  if (gQA) qaFigures();
  std::ofstream list((gFig + "/FIGURES.txt").c_str());
  list << "# Figures written by FinalAnalysis/plots/IANFigures.C from " << (f ? results : "no results file") << (gQA ? std::string(" and ") + qa : "") << "\n";
  for (const auto &w : gWritten) list << "fig/" << w << "\n";
  for (const auto &n : gNotDrawn) list << "# not drawn: fig/" << n << "\n";
  if (!list) throw std::runtime_error("IANFigures: cannot write FIGURES.txt");
  printf("IANFigures: %zu figures in %s, %zu not drawn\n", gWritten.size(), fig, gNotDrawn.size());
}
