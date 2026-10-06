/* IANFigures: every FinalAnalysis figure, in sPHENIX style, under the IAN's own file names. */

/* ./FinalAnalysis/final_analysis.sh draw --results results.root --qa qa.root --fig FinalAnalysis/fig */

/* The measurement figures draw results.root; the photon and event QA figures (QAFigures.C) draw the QA file. Each
   figure draws saved quantities and computes nothing beyond projections, ratios and normalisations for display.
   The results, their labels and their reported xJ ranges all come from results.root, so a new binning or centrality
   class needs no change here. Every xJ figure shows only the reported range, where the distribution is fully
   populated (results.root derives it from the jet threshold and the lowest photon pT); its axis starts at 0 and ends
   at the range's upper edge rounded up to 0.1. The output folder mirrors the IAN's fig/ folder. */
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
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cctype>

namespace {

/* One result of results.root, its IAN file-name tags and its reported xJ range. */
struct Result {
  std::string key, package, tag, longTag, label, energy;
  double low = 0, high = 0, jetThreshold = 0, centralityLow = -1;  /* centralityLow: Au+Au only */
  int pairIterations = 0, photonIterations = 0;
  std::string combinatoric;  /* none, simulation, or data_dphi_<shape>_<scope> */
};
const int kColors[] = {kBlack, kAzure + 2, kOrange + 7, kRed + 1, kGreen + 3, kViolet + 1, kGray + 2, kCyan + 2};
const int kMarkers[] = {20, 21, 33, 34, 24, 25, 27, 28};

std::string gFig;
std::vector<std::string> gWritten;
bool gLegacy = false;
std::string gPopulation;
std::unique_ptr<TFile> gSystematics;  /* systematics.root of final_analysis.sh systematics, when given */
using PlotJson = nlohmann::json;
PlotJson gResultMetadata, gQAMetadata, gRenderInputs, gStylePins;
std::vector<std::string> gVisibleText;
std::set<std::string> gSemanticPackages;
PlotJson gProjectionSteps = PlotJson::array();
const Result *gCurrentResult = nullptr;
bool gCurrentUsesData = true;
bool gFullPhotonResponse = false;

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
  static const PlotJson missing;
  const auto &metadata = qa ? gQAMetadata : gResultMetadata;
  gSemanticPackages.insert((qa ? "qa/" : "results/")+pkg);
  if (metadata.is_null() || !metadata.at("packages").contains(pkg)) return missing;
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
std::string jetSelectionLabel(const Result &r) {
  const auto &selection = renderSelection(r.package);
  if (selection.is_null()) return "Jet-selection metadata unavailable";
  const auto &j = selection.at("jet");
  const auto &fraction = j.at("delta_phi_pi_fraction");
  const std::string angle = fraction.is_null() ? Form("%.15g rad",j.at("delta_phi_min_rad").get<double>()) :
    Form("%d#pi/%d",fraction.at(0).get<int>(),fraction.at(1).get<int>());
  auto op = [](const PlotJson &value) {const auto text=value.get<std::string>();
    if(text==">=")return std::string("#geq"); if(text==">")return std::string(">"); throw std::runtime_error("unsupported jet-cut operator");};
  return Form("anti-k_{T} R = %g, #left|#eta^{jet}#right| < %g, p_{T}^{jet} %s %g GeV, #Delta#phi %s %s, #DeltaR_{#gamma j} %s %g",
              j.at("radius").get<double>(), j.at("abs_eta_max").get<double>(),op(j.at("pt_operator")).c_str(), r.jetThreshold,
              op(j.at("delta_phi_operator")).c_str(),angle.c_str(),op(j.at("overlap_operator")).c_str(),j.at("overlap_delta_r").get<double>());
}
std::vector<std::string> photonSelectionLines(const Result &r) {
  const auto &selection = renderSelection(r.package);
  if (selection.is_null()) return {};
  const auto &result = gResultMetadata.at("results").at(r.key);
  const auto &edges = result.at("photon_edges_gev");
  const double low = gFullPhotonResponse ? selection.at("photon").at("et_min_gev").get<double>() : edges.front().get<double>();
  const double high = gFullPhotonResponse ? selection.at("photon").at("et_max_gev").get<double>() : edges.back().get<double>();
  const auto &event = selection.at("event");
  std::vector<std::string> lines{Form("%g #leq E_{T}^{#gamma,reco} < %g GeV, #left|#eta^{#gamma}#right| < %g, #left|z_{vtx}#right| < %g cm",
    low, high, selection.at("photon").at("abs_eta_max").get<double>(), event.at("vertex_abs_z_max_cm").get<double>())};
  const auto &id = selection.at("photon_id");
  std::string region = Form("Reco: tight + isolated (R_{iso} = %g)", id.at("isolation").at("radius").get<double>());
  if (gCurrentUsesData) region += "; DATA: "+id.at("data_trigger").get<std::string>()+" ("+event.at("trigger_word").get<std::string>()+")";
  lines.push_back(region);
  if (gFullPhotonResponse) {
    const auto &truth = selection.at("truth_photon");
    lines.push_back(Form("Truth: #left|#eta^{#gamma}#right| < %g, E_{T}^{iso} (R=%g) < %g GeV",
      truth.at("abs_eta_max").get<double>(),truth.at("isolation_radius").get<double>(),truth.at("isolation_max_gev").get<double>()));
  }
  return lines;
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
  h = static_cast<TH1 *>(h->Clone((r.key + name).c_str()));
  h->SetDirectory(nullptr);
  return h;
}

/* The header: sPHENIX Internal, the lines below it and the legacy-input note. Returns the baseline of its lowest line;
   draw = false only measures. */
std::vector<std::string> compactLines(const std::vector<std::string> &input) {
  std::vector<std::string> out;
  auto width = [](const std::string &text) {
    double units=0;
    for (size_t i=0;i<text.size();++i) {
      if (text[i]=='#') {
        const auto begin=++i;
        while (i<text.size() && std::isalpha(static_cast<unsigned char>(text[i]))) ++i;
        const auto command=text.substr(begin,i-begin); --i;
        if (command!="left" && command!="right" && command!="it" && command!="bf") ++units;
      } else if (text[i]!='{' && text[i]!='}' && text[i]!='^' && text[i]!='_') units += text[i]==' '?0.5:1;
    }
    return units;
  };
  for (const auto &line : input) {
    size_t begin=0; int depth=0; std::string row;
    for (size_t at=0; at<=line.size(); ++at) {
      if (at<line.size() && line[at]=='{') ++depth;
      if (at<line.size() && line[at]=='}') --depth;
      if (at==line.size() || (!depth && (line[at]==',' || line[at]==';'))) {
        auto part=line.substr(begin,at-begin); const auto first=part.find_first_not_of(' '), last=part.find_last_not_of(' ');
        if (first!=std::string::npos) {
          part=part.substr(first,last-first+1);
          if (!row.empty() && width(row+", "+part)>34) {out.push_back(row); row.clear();}
          row += (row.empty()?"":", ")+part;
        }
        begin=at+1;
      }
    }
    if (!row.empty()) out.push_back(row);
  }
  return out;
}
double header(double x, double y, const std::vector<std::string> &lines, bool draw = true, bool compact = false) {
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.042);
  if (draw) { t.DrawLatex(x, y, sPhenixInternalLabel()); gVisibleText.push_back(sPhenixInternalLabel()); }
  t.SetTextSize(0.030);
  std::vector<std::string> text = lines;
  if (gCurrentResult) {
    const auto cuts = photonSelectionLines(*gCurrentResult);
    text.insert(text.end(), cuts.begin(), cuts.end());
  }
  if (compact) text=compactLines(text);
  double last = y;
  for (size_t i = 0; i < text.size(); ++i) {
    last -= i ? 0.041 : 0.050;
    if (draw) { t.DrawLatex(x, last, text[i].c_str()); gVisibleText.push_back(text[i]); }
  }
  if (gLegacy) {
    last -= 0.05;
    t.SetTextColor(kGray + 2);
    if (draw) { t.DrawLatex(x, last, "Reference inputs (development)"); gVisibleText.push_back("Reference inputs (development)"); }
  } else if (gPopulation == "PARTIAL" || gPopulation == "UNKNOWN") {
    last -= 0.05;
    t.SetTextColor(kGray + 2);
    if (draw) { const std::string note = gPopulation == "PARTIAL" ? "Partial input population" : "Input coverage unverified";
      t.DrawLatex(x, last, note.c_str()); gVisibleText.push_back(note); }
  }
  return last;
}

void label(const Result &r, double x = 0.2, double y = 0.90) { header(x, y, {r.label + (r.energy.empty() ? "" : ", " + r.energy)}); }

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

/* the legend of every figure: right column, below the sPHENIX header */
TLegend *legend(double y1 = 0.56, double y2 = 0.78, double x1 = 0.52) {
  auto *l = new TLegend(x1, y1, 0.93, y2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextSize(0.034);
  return l;
}

template <class T> void style(T *h, int color, int marker, double size = 1.1) {
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(size);
  h->SetLineWidth(2);
}

/* pairs summed over photon bins, per photon and per unit xJ */
TH1D *perPhoton(const TH2 *pairs, double photons, const char *name) {
  auto *h = pairs->ProjectionY(name, 1, pairs->GetNbinsX(), "e");
  h->SetDirectory(nullptr);
  if (photons > 0) h->Scale(1.0 / photons, "width");
  return h;
}

/* The bins of an xJ distribution inside [low, high], as points; widths: horizontal bars over each bin. */
TGraphErrors *points(const TH1 *h, double low, double high, int color, int marker, bool widths = true, double size = 1.1) {
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

double highest(const std::vector<TGraphErrors *> &graphs) {
  double y = 0;
  for (auto *g : graphs)
    for (int i = 0; i < g->GetN(); ++i) y = std::max(y, g->GetY()[i] + g->GetEY()[i]);
  return y;
}

/* The frame of every xJ figure: x from 0 to the reported upper edge rounded up to 0.1. */
void frameXj(double high, double ymax, const char *ytitle = "1/N_{#gamma} dN_{J#gamma}/dx_{J#gamma}") {
  TH1F *frame = gPad->DrawFrame(0.0, 0.0, std::ceil(high * 10 - 1e-9) / 10, ymax);
  frame->GetXaxis()->SetTitle("x_{J#gamma}");
  frame->GetYaxis()->SetTitle(ytitle);
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

/* Every xJ figure: the header at the top left, the legend below it (two columns past four entries), and the y range
   chosen so the highest point stays below the legend; systematic boxes, when given, sit under the points.
   entries: graph, legend text. */
void drawXj(const std::vector<std::pair<TGraphErrors *, std::string>> &entries, double high, const std::string &second,
            const char *third = nullptr, const std::vector<TGraphAsymmErrors *> &boxes = {},
            const std::vector<std::string> &extra = {}) {
  const double y = 0.90;
  std::vector<std::string> lines = {second};
  if (third) lines.push_back(third);
  lines.insert(lines.end(),extra.begin(),extra.end());
  if (gCurrentResult) lines.push_back(jetSelectionLabel(*gCurrentResult));
  const int columns = entries.size() > 4 ? 2 : 1, rows = (entries.size() + columns - 1) / columns;
  std::vector<TGraphErrors *> graphs;
  for (const auto &e : entries) graphs.push_back(e.first);
  const double bottom = gStyle->GetPadBottomMargin(), plot = 1 - bottom - gStyle->GetPadTopMargin();
  double top = highest(graphs);
  for (auto *b : boxes)
    for (int i = 0; b && i < b->GetN(); ++i) top = std::max(top, b->GetY()[i] + b->GetEYhigh()[i]);
  const double rightX=0.52, leftX=0.16;
  const double split=(rightX-gStyle->GetPadLeftMargin())/(1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())*std::ceil(high*10-1e-9)/10;
  const double leftEnd=(0.59-gStyle->GetPadLeftMargin())/(1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())*std::ceil(high*10-1e-9)/10;
  double left=0,right=0;
  for (auto *graph : graphs) for (int i=0;i<graph->GetN();++i) {
    const double value=graph->GetY()[i]+graph->GetEY()[i];
    if (graph->GetX()[i]-graph->GetEX()[i]<=leftEnd) left=std::max(left,value);
    if (graph->GetX()[i]+graph->GetEX()[i]>=split) right=std::max(right,value);
  }
  for (auto *graph : boxes) for (int i=0;graph && i<graph->GetN();++i) {
    const double value=graph->GetY()[i]+graph->GetEYhigh()[i];
    if (graph->GetX()[i]-graph->GetEXlow()[i]<=leftEnd) left=std::max(left,value);
    if (graph->GetX()[i]+graph->GetEXhigh()[i]>=split) right=std::max(right,value);
  }
  const bool onRight=right<=left; const double x=onRight?rightX:leftX;
  const double y2 = header(x,y,lines,false,true)-0.02, y1=y2-0.040*entries.size();
  if (y1<=bottom+0.05) throw std::runtime_error("xJ annotations do not fit the canvas");
  const double fraction=(y1-0.03-bottom)/plot, quiet=onRight?right:left;
  frameXj(high,std::max(top*1.12,quiet*1.06/std::max(fraction,0.1)));
  for (auto *b : boxes) if (b) b->Draw("2");
  auto *leg = new TLegend(x, y1, onRight?0.95:0.59, y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  leg->SetNColumns(1);
  for (const auto &e : entries) {
    e.first->Draw("P");
    leg->AddEntry(e.first, e.second.c_str(), "lp");
  }
  leg->Draw();
  header(x, y, lines, true, true);
}

/* Photon-purity correction of the reconstructed xJ distribution: region A, the transferred region-C background and
   the leakage-restored, purity-corrected distribution, each per photon. */
void purityCorrection(TFile &f, const Result &r, const std::string &name) {
  auto *a = dynamic_cast<TH2 *>(get(f, r, "pairs_A"));
  auto *c = dynamic_cast<TH2 *>(get(f, r, "pairs_C_scaled"));
  auto *p = dynamic_cast<TH2 *>(get(f, r, "pairs_purity_corrected"));
  auto *abcd = dynamic_cast<TH2 *>(get(f, r, "abcd"));
  auto *photons = get(f, r, "photons_corrected");
  if (!a || !c || !p || !abcd || !photons) return;
  const double nA = abcd->Integral(1, 1, 1, abcd->GetNbinsY());
  TCanvas cv("cv", "", 700, 600);
  drawXj({{points(perPhoton(a, nA, "ha"), r.low, r.high, kBlack, 20), "Region A"},
          {points(perPhoton(c, nA, "hc"), r.low, r.high, kOrange + 7, 21), "k #times region C"},
          {points(perPhoton(p, photons->Integral(), "hp"), r.low, r.high, kAzure + 2, 24), "Purity corrected"}},
         r.high, r.label + ", " + r.energy);
  save(cv, name);
}

/* Combinatoric-jet subtraction: the purity-corrected distribution, the scaled simulated template and the result. */
void combinatoric(TFile &f, const Result &r, const std::string &name) {
  auto *p = dynamic_cast<TH2 *>(get(f, r, "pairs_purity_corrected"));
  auto *k = dynamic_cast<TH2 *>(get(f, r, "pairs_combinatoric"));
  auto *s = dynamic_cast<TH2 *>(get(f, r, "pairs_corrected"));
  auto *photons = get(f, r, "photons_corrected");
  if (!p || !k || !s || !photons) return;
  const double n = photons->Integral();
  TCanvas cv("cv", "", 700, 600);
  drawXj({{points(perPhoton(p, n, "hp2"), r.low, r.high, kBlack, 24), "Purity corrected"},
          {points(perPhoton(k, n, "hk2"), r.low, r.high, kRed + 1, 21),
           r.combinatoric.rfind("data_dphi", 0) == 0 ? "Combinatoric jets (data, #Delta#phi sideband)"
           : k->Integral() > 0 ? "Combinatoric template" : "No combinatoric subtraction"},
          {points(perPhoton(s, n, "hs2"), r.low, r.high, kAzure + 2, 20), "Combinatoric subtracted"}},
         r.high, r.label + ", " + r.energy);
  save(cv, name);
}

/* Legible numbers on a coloured response matrix, using the same palette as ROOT's COLZ draw. */
void matrixText(TH2 &h) {
  TLatex text; text.SetTextAlign(22); text.SetTextFont(42); text.SetTextSize(0.022);
  const int colors = gStyle->GetNumberOfColors();
  const double low=h.GetMinimum(), high=h.GetMaximum();
  for (int x=1;x<=h.GetNbinsX();++x) for (int y=1;y<=h.GetNbinsY();++y) {
    const double value=h.GetBinContent(x,y);
    if (!std::isfinite(value) || value==0) continue;
    const double fraction=high>low?std::clamp((value-low)/(high-low),0.0,1.0):0.5;
    auto *color=gROOT->GetColor(TColor::GetColorPalette(std::clamp(int(fraction*(colors-1)),0,colors-1)));
    const double light=color?0.2126*color->GetRed()+0.7152*color->GetGreen()+0.0722*color->GetBlue():1;
    text.SetTextColor(light<0.5?kWhite:kBlack);
    text.DrawLatex(h.GetXaxis()->GetBinCenter(x),h.GetYaxis()->GetBinCenter(y),Form("%.2f",value));
  }
}

/* Photon response on the response grid, each truth bin normalised to unit probability. */
void photonResponse(TFile &f, const Result &r, const std::string &name) {
  auto *m = dynamic_cast<TH2 *>(get(f, r, "photon_response_grid"));
  if (!m) return;
  for (int i = 1; i <= m->GetNbinsX(); ++i) {
    double row = 0;
    for (int j = 1; j <= m->GetNbinsY(); ++j) row += m->GetBinContent(i, j);
    for (int j = 1; j <= m->GetNbinsY(); ++j) m->SetBinContent(i, j, row > 0 ? m->GetBinContent(i, j) / row : 0);
  }
  TCanvas cv("cv", "", 760, 700);
  cv.SetRightMargin(0.16);
  cv.SetTopMargin(1-header(0.16, 0.95, {r.label + (r.energy.empty() ? "" : ", " + r.energy)}, false)+0.035);
  gStyle->SetPaintTextFormat(".2f");
  m->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
  m->GetYaxis()->SetTitle("E_{T}^{#gamma, reco} [GeV]");
  m->GetZaxis()->SetTitle("P(reco | truth)");
  m->SetMinimum(0);
  m->SetMaximum(1);
  m->SetMarkerSize(1.2);
  m->Draw("COLZ"); matrixText(*m);
  header(0.16, 0.95, {r.label + (r.energy.empty() ? "" : ", " + r.energy)});
  save(cv, name);
}

/* Joint (photon pT, xJ) response of the measurement window: truth cell x reco cell. */
void jointResponse(TFile &f, const Result &r, const std::string &name) {
  auto *m = dynamic_cast<TH2 *>(get(f, r, "pair_response"));
  if (!m) return;
  TCanvas cv("cv", "", 760, 700);
  cv.SetRightMargin(0.16);
  cv.SetTopMargin(1-header(0.16, 0.95, {r.label + (r.energy.empty() ? "" : ", " + r.energy)}, false)+0.035);
  cv.SetLogz();
  m->GetXaxis()->SetTitle("truth cell (p_{T}^{#gamma} bin #times x_{J#gamma} bin)");
  m->GetYaxis()->SetTitle("reco cell (E_{T}^{#gamma} bin #times x_{J#gamma} bin)");
  m->Draw("COLZ");
  header(0.16, 0.95, {r.label + (r.energy.empty() ? "" : ", " + r.energy)});
  save(cv, name);
}

/* The per-photon xJ distribution unfolded with 1 to N iterations. */
void jointIterations(TFile &f, const Result &r, const std::string &name) {
  auto *scan = dynamic_cast<TH2 *>(get(f, r, "scan_xj_per_photon"));
  if (!scan) return;
  TCanvas cv("cv", "", 700, 600);
  std::vector<std::pair<TGraphErrors *, std::string>> rows;
  for (int k = 1; k <= scan->GetNbinsX(); ++k) {
    auto *h = scan->ProjectionY(Form("k%d", k), k, k);
    h->SetDirectory(nullptr);
    rows.push_back({points(h, r.low, r.high, kColors[(k - 1) % 8], kMarkers[(k - 1) % 8], false, 0.9),
                    Form("%d iteration%s", k, k > 1 ? "s" : "")});
  }
  drawXj(rows, r.high, r.label + ", " + r.energy, Form("Nominal: %d iterations", r.pairIterations));
  save(cv, name);
}

/* Unfolded photon yield of every photon bin against the iteration count, relative to the nominal count. */
void photonIterations(TFile &f, const Result &r, const std::string &name) {
  auto *scan = dynamic_cast<TH2 *>(get(f, r, "scan_photons"));
  if (!scan) return;
  const int nominal = std::max(1, r.photonIterations);
  TCanvas cv("cv", "", 700, 600);
  TLegend &leg = *legend();
  for (int b = 1; b <= scan->GetNbinsY(); ++b) {
    auto *h = scan->ProjectionX(Form("pt%d", b), b, b);
    h->SetDirectory(nullptr);
    const double ref = h->GetBinContent(nominal);
    if (ref > 0) h->Scale(1.0 / ref);
    style(h, kColors[b % 8], 20 + b % 4);
    h->GetXaxis()->SetTitle("iterations");
    h->GetYaxis()->SetTitle(Form("N_{#gamma}(k) / N_{#gamma}(k = %d)", nominal));
    h->SetMinimum(0.9);
    h->SetMaximum(1.1);
    h->Draw(b == 1 ? "P" : "P SAME");
    leg.AddEntry(h, Form("%g #leq p_{T}^{#gamma} < %g GeV", scan->GetYaxis()->GetBinLowEdge(b), scan->GetYaxis()->GetBinUpEdge(b)), "p");
  }
  leg.Draw();
  label(r);
  save(cv, name);
}

/* Simulation closure of the photon-yield unfolding: unfolded simulated reconstruction against its truth. */
void photonClosure(TFile &f, const Result &r, const std::string &name) {
  auto *u = get(f, r, "closure_photons");
  auto *t = get(f, r, "photon_truth");
  if (!u || !t) return;
  TCanvas cv("cv", "", 700, 700);
  TPad top("top", "", 0, 0.3, 1, 1), bottom("bottom", "", 0, 0, 1, 0.3);
  top.SetBottomMargin(0.02); bottom.SetTopMargin(0.02); bottom.SetBottomMargin(0.35);
  top.Draw(); bottom.Draw();
  top.cd();
  top.SetLogy();
  style(t, kBlack, 1); t->SetLineWidth(3);
  style(u, kAzure + 2, 24);
  t->GetYaxis()->SetTitle("photons (simulation, weighted)");
  t->GetXaxis()->SetLabelSize(0);
  t->Draw("HIST");
  u->Draw("P SAME");
  TLegend &leg = *legend();
  leg.AddEntry(t, "Truth", "l");
  leg.AddEntry(u, "Unfolded reconstruction", "p");
  leg.Draw();
  label(r);
  bottom.cd();
  auto *ratio = static_cast<TH1 *>(u->Clone("ratio"));
  ratio->Divide(t);
  ratio->SetMinimum(0.9); ratio->SetMaximum(1.1);
  ratio->GetYaxis()->SetTitle("unfolded / truth");
  ratio->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  ratio->GetYaxis()->SetTitleSize(0.1); ratio->GetYaxis()->SetLabelSize(0.09); ratio->GetYaxis()->SetTitleOffset(0.6);
  ratio->GetXaxis()->SetTitleSize(0.12); ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->Draw("P");
  TLine one(ratio->GetXaxis()->GetXmin(), 1, ratio->GetXaxis()->GetXmax(), 1);
  one.SetLineStyle(2);
  one.Draw();
  save(cv, name);
}

/* Purity-corrected reconstructed photon yields and the unfolded yields. */
void photonYield(TFile &f, const Result &r, const std::string &name) {
  auto *c = get(f, r, "photons_corrected");
  auto *u = get(f, r, "photons_unfolded");
  if (!c || !u) return;
  style(c, kBlack, 20); style(u, kAzure + 2, 24);
  TCanvas cv("cv", "", 700, 600);
  c->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  c->GetYaxis()->SetTitle("photons");
  c->SetMinimum(0);
  c->SetMaximum(1.5 * std::max(c->GetMaximum(), u->GetMaximum()));
  c->Draw("P");
  u->Draw("P SAME");
  TLegend &leg = *legend();
  leg.AddEntry(c, "Purity-corrected reco.", "p");
  leg.AddEntry(u, Form("Unfolded, %d iterations", r.photonIterations), "p");
  leg.Draw();
  label(r);
  save(cv, name);
}

/* The joint unfolding's prior: the simulated truth, the truth reweighted to the data, and the unfolded data from each. */
void priorReweighting(TFile &f, const Result &r, const std::string &name) {
  auto *simulation = get(f, r, "xj_prior_simulation");
  auto *reweighted = get(f, r, "xj_prior_reweighted");
  auto *nominal = get(f, r, "xj_per_photon");
  auto *fromSimulation = get(f, r, "xj_per_photon_simulation_prior");
  if (!simulation || !reweighted || !nominal || !fromSimulation) return;
  TCanvas cv("cv", "", 700, 600);
  drawXj({{points(simulation, r.low, r.high, kGray + 2, 25), "Simulation truth"},
          {points(reweighted, r.low, r.high, kOrange + 7, 26), "Truth reweighted to the data"},
          {points(fromSimulation, r.low, r.high, kRed + 1, 24), "Data, simulation prior"},
          {points(nominal, r.low, r.high, kAzure + 2, 20), "Data, reweighted prior"}},
         r.high, r.label + ", " + r.energy);
  save(cv, name);
}

/* Per-photon xJ before (purity- and combinatoric-corrected) and after unfolding. */
void beforeAfter(TFile &f, const Result &r, const std::string &name) {
  auto *b = get(f, r, "xj_per_photon_reco");
  auto *a = get(f, r, "xj_per_photon");
  if (!a || !b) return;
  TCanvas cv("cv", "", 700, 600);
  drawXj({{points(b, r.low, r.high, kBlack, 24), "Corrected, before unfolding"},
          {points(a, r.low, r.high, kAzure + 2, 20), "Unfolded (statistical)"}},
         r.high, r.label + ", " + r.energy);
  save(cv, name);
}

/* The unfolded per-photon xJ of p+p (atlas: "pp") or of the most central Au+Au class ("pbpb_0_10") against ATLAS
   (Phys. Lett. B 789 (2019) 167, HEPData ins1694678 Table 1: p+p and 0-10% Pb+Pb at 5.02 TeV, plots/inputs/
   atlas_reference.json), the ATLAS systematic uncertainty shaded. */
void atlasOverlay(TFile &f, const Result &r, const std::string &atlas, const std::string &name) {
  auto *h = get(f, r, "xj_per_photon");
  std::ifstream in(std::string(gSystem->GetDirName(__FILE__).Data()) + "/inputs/atlas_reference.json");
  if (!h || !in) return;
  const nlohmann::json ref = nlohmann::json::parse(in), &sel = ref.at("selection");
  auto *points_ = new TGraphErrors;
  auto *boxes = new TGraphAsymmErrors;
  for (const auto &b : ref.at("bins")) {
    const double lo = b.at("low"), hi = b.at("high"), y = b.at(atlas).at("value"), stat = b.at(atlas).at("stat"), sys = b.at(atlas).at("sys");
    const int n = points_->GetN();
    points_->SetPoint(n, 0.5 * (lo + hi), y), points_->SetPointError(n, 0.5 * (hi - lo), stat);
    boxes->SetPoint(n, 0.5 * (lo + hi), y), boxes->SetPointError(n, 0.5 * (hi - lo), 0.5 * (hi - lo), sys, sys);
  }
  style(points_, kOrange + 7, 25);
  boxes->SetFillColorAlpha(kOrange + 7, 0.2);
  boxes->SetLineWidth(0);
  const auto photons = sel.at("photon_pt_gev");
  const std::string system = atlas == "pp" ? "pp" : Form("Pb+Pb %g-%g%%", sel.at("pbpb_centrality_percent")[0].get<double>(),
                                                          sel.at("pbpb_centrality_percent")[1].get<double>());
  TCanvas cv("cv", "", 700, 600);
  std::vector<TGraphAsymmErrors *> shaded = {boxes};
  if (auto *ours = systematicBoxes(r.key, r.low, r.high, kAzure + 2)) shaded.push_back(ours);
  drawXj({{points(h, r.low, r.high, kAzure + 2, 20), r.label + ", " + r.energy},
          {points_, Form("ATLAS %s, %g TeV", system.c_str(), sel.at("sqrt_s_nn_gev").get<double>() / 1000)}},
         r.high, "Unfolded measurements",
         Form("ATLAS: %g < p_{T}^{#gamma} < %g GeV, p_{T}^{jet} #geq %g GeV, shading syst.", photons[0].get<double>(),
              photons[1].get<double>(), sel.at("jet_pt_min_gev").get<double>()),
         shaded);
  save(cv, name);
}

/* The unfolded per-photon xJ of every result on one plot, over the range all of them report, with the systematic
   boxes when a systematics file is given. */
void centralityOverlay(TFile &f, const std::vector<Result> &results, const std::string &name) {
  if (results.empty()) return;
  double low = 0, high = 1e9;
  for (const auto &r : results) { low = std::max(low, r.low); high = std::min(high, r.high); }
  TCanvas cv("cv", "", 760, 620);
  std::vector<std::pair<TGraphErrors *, std::string>> rows;
  std::vector<TGraphAsymmErrors *> boxes;
  int i = 0;
  for (const auto &r : results)
    if (auto *h = get(f, r, "xj_per_photon")) {
      rows.push_back({points(h, low, high, kColors[i % 8], kMarkers[i % 8], true, 1.2), r.label});
      boxes.push_back(systematicBoxes(r.key, low, high, kColors[i % 8]));
      ++i;
    }
  std::vector<std::string> cuts;
  std::set<std::string> shown;
  for (const auto &r : results) {
    const auto &selection = renderSelection(r.package);
    if (selection.is_null()) continue;
    const auto &edges = gResultMetadata.at("results").at(r.key).at("photon_edges_gev");
    const auto &id = selection.at("photon_id");
    const auto &event = selection.at("event");
    for (const auto &line : std::vector<std::string>{Form("%s, %s", r.package=="PP" ? "p+p" : "Au+Au",energyLabel(r.package).c_str()),
             Form("%g #leq E_{T}^{#gamma} < %g GeV, #left|#eta^{#gamma}#right| < %g",edges.front().get<double>(),edges.back().get<double>(),
                  selection.at("photon").at("abs_eta_max").get<double>()),
             Form("%s: #left|z_{vtx}#right| < %g cm, DATA %s, tight + isolated (R=%g)",r.package=="PP" ? "p+p" : "Au+Au",
                  event.at("vertex_abs_z_max_cm").get<double>(),id.at("data_trigger").get<std::string>().c_str(),
                  id.at("isolation").at("radius").get<double>())})
      if (shown.insert(line).second) cuts.push_back(line);
    const auto jet = jetSelectionLabel(r);
    if (shown.insert(jet).second) cuts.push_back(jet);
  }
  drawXj(rows, high, gSystematics ? "Unfolded; bars statistical, boxes systematic" : "Unfolded, statistical errors", nullptr, boxes, cuts);
  save(cv, name);
}

}  /* namespace */

#include "QAFigures.C"

/* results: the results file (or ""), fig: the output folder, systematics and qa: optional systematics and QA files.
   Draw saved quantities with their verified metadata; FIGURES.txt records written and unavailable figures. */
void IANFigures(const char *results, const char *fig, const char *systematics = "", const char *qa = "") {
  gWritten.clear(); gNotDrawn.clear(); gSystematics.reset(); gQA.reset(); gLegacy = false;
  gResultMetadata = gQAMetadata = PlotJson(); gRenderInputs = gStylePins = PlotJson::object(); gVisibleText.clear();
  gSemanticPackages.clear(); gProjectionSteps = PlotJson::array(); gCurrentResult = nullptr;
  const auto source = std::filesystem::path(__FILE__).parent_path();
  for (const auto &name : {"IANFigures.C", "QAFigures.C", "sPhenixStyle.C", "sPhenixStyle.h"})
    gStylePins[name] = photonjet_render::sha256((source/name).string(), true);
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
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
    gCurrentResult = &r; gCurrentUsesData = true;
    const bool pp = r.tag == "pp";
    purityCorrection(*f, r, "purityAndCombSub/" + r.tag + "_xj_purity_correction.png");
    combinatoric(*f, r, "purityAndCombSub/" + r.tag + (pp ? "_xj_combinatoric_crosscheck.png" : "_xj_combinatoric_subtraction.png"));
    gCurrentUsesData = false;
    gFullPhotonResponse = true;
    photonResponse(*f, r, "Unfolding/" + r.longTag + "_photon_response.png");
    gFullPhotonResponse = false;
    jointResponse(*f, r, "Unfolding/" + r.longTag + "_joint_response.png");
    gCurrentUsesData = true;
    photonIterations(*f, r, "Unfolding/" + r.tag + (pp ? "_photon_yield_iteration_stability_publication.png" : "_photon_yield_iteration_stability.png"));
    jointIterations(*f, r, "Unfolding/" + r.tag + "_joint_iteration_stability.png");
    gCurrentUsesData = false;
    photonClosure(*f, r, "Unfolding/" + r.tag + "_photon_closure_overlay.png");
    gCurrentUsesData = true;
    photonYield(*f, r, "Unfolding/" + r.tag + "_photon_yield_unfolding_nominal_2bin_k" + std::to_string(r.photonIterations) + ".png");
    priorReweighting(*f, r, "Unfolding/" + r.tag + "_xj_prior_reweighting.png");
    beforeAfter(*f, r, "FinalResults/" + r.tag + "_xj_before_after_unfolding.png");
    if (pp || r.centralityLow == 0) atlasOverlay(*f, r, pp ? "pp" : "pbpb_0_10", "LHCComparison/" + r.tag + "_xj_atlas_overlay.png");
  }
  gCurrentResult = nullptr;
  if (f) centralityOverlay(*f, all, "FinalResults/pp_vs_auau_centrality_xj_unfolded.png");
  if (gQA) qaFigures();
  std::ofstream list((gFig + "/FIGURES.txt").c_str());
  list << "# Figures written by FinalAnalysis/plots/IANFigures.C from " << (f ? results : "no results file") << (gQA ? std::string(" and ") + qa : "") << "\n";
  for (const auto &w : gWritten) list << "fig/" << w << "\n";
  for (const auto &n : gNotDrawn) list << "# not drawn: fig/" << n << "\n";
  if (!list) throw std::runtime_error("IANFigures: cannot write FIGURES.txt");
  printf("IANFigures: %zu figures in %s, %zu not drawn\n", gWritten.size(), fig, gNotDrawn.size());
}
