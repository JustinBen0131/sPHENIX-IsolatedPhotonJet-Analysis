/* QAFigures: the photon and event QA figures of the IAN, from the QA file of final_analysis.sh qa. */

/* Included by IANFigures.C, whose helpers it uses (header, save, style), and called from its entry point with the QA
   file. Each figure draws the saved sums and numerics under qa/<package>/ of that file: the DATA and simulated QA
   histograms (data/, phosim/, jetsim/; the simulation summed over samples, each weighted by its cross section per
   generated event) and the ABCD counts, purity, purity closure, stitching fits, working points, trigger study and
   centrality target of the qa step. Drawing projects, rebins and normalises for display only. The classes, samples and
   triggers are the ones in the file. A figure whose inputs are absent or empty is not drawn and is listed in
   FIGURES.txt with the reason; nothing is invented. */
#include <TColor.h>
#include <TF1.h>
#include <TH3.h>
#include <TParameter.h>

#include <algorithm>
#include <array>
#include <functional>
#include <limits>
#include <tuple>

namespace {

std::unique_ptr<TFile> gQA;
std::vector<std::string> gNotDrawn;
int gUnique = 0;
TH1F *gFrame = nullptr;                /* the last single-panel frame and its y range */
double gFrameLow = 0, gFrameHigh = 1;
const int cBlack = kBlack, cBlue = kAzure + 2, cOrange = kOrange + 7, cGreen = kGreen + 3, cRed = kRed + 1, cPurple = kViolet + 1;

std::string fresh() { return Form("qa%d", ++gUnique); }

void missing(const std::string &name, const std::string &why) {
  gNotDrawn.push_back(name + ": " + why);
  printf("IANFigures: not drawn %s (%s)\n", name.c_str(), why.c_str());
}

bool empty(const TH1 *h) {
  for (int b = 0; b < h->GetNcells(); ++b)
    if (!h->IsBinUnderflow(b) && !h->IsBinOverflow(b) && (h->GetBinContent(b) != 0 || h->GetBinError(b) != 0)) return false;
  return true;
}

/* qa/<path> of the QA file, detached; null when absent or empty. */
TH1 *qa(const std::string &path) {
  auto *h = dynamic_cast<TH1 *>(gQA->Get(("qa/" + path).c_str()));
  if (!h || empty(h)) return nullptr;
  renderSelection(path.substr(0,path.find('/')), true);
  h = static_cast<TH1 *>(h->Clone(fresh().c_str()));
  h->SetDirectory(nullptr);
  return h;
}

TH1 *joined(TH1 *a, TH1 *b) {
  if (!a) return b;
  if (b) a->Add(b);
  return a;
}

/* ---- classes and candidates ------------------------------------------------------------------------------------ */

PlotJson classRange(const std::string &c) {
  if (!gQAMetadata.is_null()) for (const auto &package : gQAMetadata.at("packages")) {
    const auto &selection = package.at("selection");
    const auto &classes = selection.is_null() ? package.at("classes") : selection.at("qa").at("classes");
    for (const auto &item : classes) if (item.at("id").get<std::string>()==c) return item.at("centrality_percent");
  }
  return PlotJson();
}
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

/* The classes of a package (its DATA photon directories), most central first: p+p "all", Au+Au centA_B. */
std::vector<std::string> classesOf(const std::string &pkg) {
  std::vector<std::string> out;
  if (auto *d = gQA->GetDirectory(("qa/" + pkg + "/data/photons").c_str())) {
    TIter next(d->GetListOfKeys());
    while (auto *k = static_cast<TKey *>(next()))
      if (std::find(out.begin(), out.end(), k->GetName()) == out.end()) out.push_back(k->GetName());
  }
  std::sort(out.begin(), out.end(), [](const std::string &a, const std::string &b) { return lowOf(a) < lowOf(b); });
  return out;
}

std::string percent(double a, double b) { return Form("%g-%g%%", a, b); }
std::string classLabel(const std::string &c) { return c.rfind("cent", 0) == 0 ? "Au+Au " + percent(lowOf(c), highOf(c)) : "p+p"; }
std::string classTag(const std::string &c) { return Form("%gto%g", lowOf(c), highOf(c)); }
std::string energyOf(const std::string &pkg) { return energyLabel(pkg, true); }
std::pair<double, double> photonSupport(const std::string &pkg) {
  const auto &m = renderSelection(pkg, true);
  if (m.is_null()) return {15, 35};  /* historical diagnostic projection, never qualified annotations */
  return {m.at("photon").at("et_min_gev"), m.at("photon").at("et_max_gev")};
}
std::string photonWindowLabel(double low, double high) { return Form("%g #leq E_{T}^{#gamma,reco} < %g GeV", low, high); }
std::string photonEtaLabel(const std::string &pkg) {
  const auto &m = renderSelection(pkg, true);
  return m.is_null() ? "Selection metadata unavailable" : Form("#left|#eta^{#gamma}#right| < %g", m.at("photon").at("abs_eta_max").get<double>());
}
double selectionValue(const std::string &pkg, const char *section, const char *key, double legacy) {
  const auto &m = renderSelection(pkg, true);
  return m.is_null() ? legacy : m.at(section).at(key).get<double>();
}
double isolationRadius(const std::string &pkg) {
  const auto &m = renderSelection(pkg, true);
  return m.is_null() ? (pkg=="PP" ? 0.4 : 0.3) : m.at("photon_id").at("isolation").at("radius").get<double>();
}
std::string truthSelectionLabel(const std::string &pkg) {
  const double radius = selectionValue(pkg,"truth_photon","isolation_radius",0.3);
  const double eta = selectionValue(pkg,"qa","stored_abs_eta_max",0.7);
  return Form("Truth isolation R = %g, #left|#eta^{#gamma}_{truth}#right| < %g",radius,eta);
}
/* All Au+Au classes together: "Au+Au 0-80%" and "0to80". */
std::string spanLabel() { const auto c = classesOf("AuAu"); return c.empty() ? "Au+Au" : "Au+Au " + percent(lowOf(c.front()), highOf(c.back())); }
std::string spanTag() { const auto c = classesOf("AuAu"); return c.empty() ? "" : Form("%gto%g", lowOf(c.front()), highOf(c.back())); }

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

std::vector<std::string> labelsOf(const std::string &family) {
  if (family == "data") return {"data"};
  if (family == "phosim") return {"sig", "phoother"};
  return {"bkg", "jetprompt"};
}

/* The candidates of a family in one class: photons/<class>/<label>/<stem>, summed over the family's labels. */
TH1 *candidates(const std::string &pkg, const std::string &family, const std::string &c, const std::string &stem) {
  TH1 *total = nullptr;
  for (const auto &l : labelsOf(family)) total = joined(total, qa(pkg + "/" + family + "/photons/" + c + "/" + l + "/" + stem));
  return total;
}

TH1 *overClasses(const std::string &pkg, const std::function<TH1 *(const std::string &)> &fetch) {
  TH1 *total = nullptr;
  for (const auto &c : classesOf(pkg)) total = joined(total, fetch(c));
  return total;
}

/* ---- projections ----------------------------------------------------------------------------------------------- */

/* The bins of one axis: values and variances. */
struct Bins {
  std::vector<double> edges, v, e2;
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
struct Cut { int axis; double lo, hi; };

/* The sum of a histogram over every axis but `axis` (no flow bins), the other axes restricted to the bins inside the
   windows given. */
Bins project(const TH1 *h, int axis, const std::vector<Cut> &windows = {}) {
  const TAxis *axes[3] = {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()};
  const int dim = h->GetDimension();
  PlotJson cuts = PlotJson::array();
  for (const auto &w : windows) cuts.push_back({{"axis",w.axis},{"low",w.lo},{"high",w.hi},{"interval","[low,high)"}});
  gProjectionSteps.push_back({{"histogram",h->GetTitle()},{"retained_axis",axis},{"windows",cuts},{"flow_bins","excluded"}});
  Bins out;
  for (int b = 1; b <= axes[axis]->GetNbins() + 1; ++b) out.edges.push_back(axes[axis]->GetBinLowEdge(b));
  out.v.assign(out.edges.size() - 1, 0.0);
  out.e2 = out.v;
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
        out.v[at[axis] - 1] += h->GetBinContent(bin);
        out.e2[at[axis] - 1] += e * e;
      }
  return out;
}

Bins scaled(Bins b, double f) {
  for (int i = 0; i < b.n(); ++i) b.v[i] *= f, b.e2[i] *= f * f;
  return b;
}
Bins unit(const Bins &b) { return b.sum() > 0 ? scaled(b, 1 / b.sum()) : b; }
Bins perWidth(Bins b) {
  for (int i = 0; i < b.n(); ++i) b.v[i] /= b.width(i), b.e2[i] /= b.width(i) * b.width(i);
  return b;
}
/* The bins inside [lo, hi]. */
Bins within(const Bins &b, double lo, double hi) {
  Bins out;
  for (int i = 0; i < b.n(); ++i)
    if (b.inside(i, lo, hi)) {
      if (out.edges.empty()) out.edges.push_back(b.edges[i]);
      out.edges.push_back(b.edges[i + 1]), out.v.push_back(b.v[i]), out.e2.push_back(b.e2[i]);
    }
  return out;
}
/* Groups of k bins (a remainder is dropped). */
Bins rebin(const Bins &b, int k) {
  Bins out;
  const int n = b.n() / k * k;
  for (int i = 0; i <= n; i += k) out.edges.push_back(b.edges[i]);
  for (int i = 0; i < n; i += k) {
    double v = 0, e2 = 0;
    for (int j = i; j < i + k; ++j) v += b.v[j], e2 += b.e2[j];
    out.v.push_back(v), out.e2.push_back(e2);
  }
  return out;
}
Bins sum(Bins a, const Bins &b) {
  for (int i = 0; i < a.n(); ++i) a.v[i] += b.v[i], a.e2[i] += b.e2[i];
  return a;
}
/* a / b per bin (NaN where b is not positive), without errors. */
Bins over(const Bins &a, const Bins &b) {
  Bins out = a;
  for (int i = 0; i < a.n(); ++i) out.v[i] = b.v[i] > 0 ? a.v[i] / b.v[i] : std::nan(""), out.e2[i] = 0;
  return out;
}
/* The bins of a saved one-dimensional result (purity, thresholds, spectra). */
Bins binsOf(const TH1 *h) {
  Bins out = project(h, 0);
  return out;
}

std::vector<double> edgesOf(const TAxis *a) {
  std::vector<double> e;
  for (int b = 1; b <= a->GetNbins() + 1; ++b) e.push_back(a->GetBinLowEdge(b));
  return e;
}

/* ---- drawing --------------------------------------------------------------------------------------------------- */

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
           reaches zero up to rounding (a single weighted entry) runs to the frame instead of stretching it over decades. */
        if (lower > 1e-3 * b.v[i]) add(lower);
        else if (std::isfinite(lower)) low = std::min(low, lower);
      }
  }
  double lowFor(bool logy) const { return logy ? positive : low; }
};

TH1D *line(const Bins &b, int color, int style = 1, int width = 2) {
  auto *h = new TH1D(fresh().c_str(), "", b.n(), b.edges.data());
  h->SetDirectory(nullptr);
  for (int i = 0; i < b.n(); ++i) h->SetBinContent(i + 1, b.v[i]);
  h->SetLineColor(color);
  h->SetLineStyle(style);
  h->SetLineWidth(width);
  return h;
}

/* Points at the bin centres (shifted by `shift`) with their errors; x errors: half the bin widths. Non-finite values
   are skipped. */
TGraphErrors *markers(const Bins &b, int color, int marker, bool errors = true, bool xerrors = false, double shift = 0, double size = 1.1) {
  auto *g = new TGraphErrors;
  for (int i = 0; i < b.n(); ++i) {
    if (!std::isfinite(b.v[i])) continue;
    const int n = g->GetN();
    g->SetPoint(n, b.center(i) + shift, b.v[i]);
    g->SetPointError(n, xerrors ? 0.5 * b.width(i) : 0.0, errors && std::isfinite(b.e2[i]) ? std::sqrt(b.e2[i]) : 0.0);
  }
  style(g, color, marker, size);
  return g;
}

void vline(double x, double y0, double y1, int style = 2) {
  TLine l;
  l.SetLineStyle(style);
  l.SetLineWidth(1);
  l.DrawLine(x, y0, x, y1);
}
void hline(double y, double x0, double x1, int style = 2) {
  TLine l;
  l.SetLineStyle(style);
  l.SetLineWidth(1);
  l.DrawLine(x0, y, x1, y);
}

/* A single-panel frame: the header at the top left and a legend of `entries` below it; the y range keeps the drawn
   values [low, high] (on a log scale low is the smallest positive value) below the legend. Returns the legend. */
TLegend *frame(double x0, double x1, double low, double high, bool logy, const char *xtitle, const char *ytitle,
               const std::vector<std::string> &lines, int entries, double x = 0.16, double quietHigh = std::nan("")) {
  const double y = 0.90;
  const int columns = entries > 4 ? 2 : 1, rows = (entries + columns - 1) / columns;
  const bool compact = std::isfinite(quietHigh);
  const double y2 = header(x, y, lines, false, compact) - 0.02, y1 = y2 - 0.042 * rows;
  const double bottom = gPad->GetBottomMargin(), plot = 1 - bottom - gPad->GetTopMargin();
  const double f = std::max(0.25, (y1 - 0.03 - bottom) / plot);
  if (!std::isfinite(high)) high = 1;
  if (logy) {
    gFrameLow = std::isfinite(low) && low > 0 ? low / 2 : std::max(high, 1e-9) * 1e-4;
    gFrameHigh = compact ? std::max(high*1.15,gFrameLow*std::pow(std::max(quietHigh,2*gFrameLow)/gFrameLow,1/f)) :
                          gFrameLow * std::pow(std::max(high, 2 * gFrameLow) / gFrameLow, 1 / f);
  } else {
    gFrameLow = std::isfinite(low) ? std::min(0.0, low) : 0.0;
    if (gFrameLow<0) gFrameLow-=0.04*std::max(high-gFrameLow,1e-9);
    gFrameHigh = compact ? std::max(high+(high-gFrameLow)*0.12,gFrameLow+(quietHigh-gFrameLow)/f) : gFrameLow + (high - gFrameLow) / f;
    if (!(gFrameHigh > gFrameLow)) gFrameHigh = gFrameLow + 1;
  }
  gPad->SetLogy(logy);
  gFrame = gPad->DrawFrame(x0, gFrameLow, x1, gFrameHigh);
  gFrame->GetXaxis()->SetTitle(xtitle);
  gFrame->GetYaxis()->SetTitle(ytitle);
  header(x, y, lines, true, compact);
  auto *leg = new TLegend(x, y1, compact && x<0.4 ? 0.59 : 0.95, y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  leg->SetNColumns(columns);
  return leg;
}

/* One curve of a single-panel figure: a line (simulation) or points (DATA). */
struct Curve {
  Bins bins;
  int color = cBlack, marker = 20;
  std::string text;
  bool points = false, errors = true, xerrors = false;
  double shift = 0;
};
Curve lineCurve(const Bins &b, int color, const std::string &text) { return Curve{b, color, 1, text}; }
Curve pointCurve(const Bins &b, int color, int marker, const std::string &text, bool errors = true) {
  Curve c{b, color, marker, text, true, errors};
  return c;
}

/* A single-panel figure of curves: lines under points, x from x0 to x1, the y range from the drawn values. */
void drawCurves(const std::vector<Curve> &curves, double x0, double x1, bool logy, const char *xtitle, const char *ytitle,
                const std::vector<std::string> &lines, const std::string &name, const std::function<void()> &extra = nullptr) {
  TCanvas cv("cv", "", 700, 600);
  Span s;
  for (const auto &c : curves) s.add(c.bins, c.points && c.errors, x0, x1);
  const double split=(0.52-gStyle->GetPadLeftMargin())/(1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())*(x1-x0)+x0;
  const double leftEnd=(0.59-gStyle->GetPadLeftMargin())/(1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())*(x1-x0)+x0;
  Span left,right;
  for (const auto &c : curves) {
    left.add(c.bins,c.points && c.errors,x0,leftEnd); right.add(c.bins,c.points && c.errors,split,x1);
  }
  const bool onRight=right.high<=left.high;
  const double quiet=std::isfinite(onRight?right.high:left.high)?(onRight?right.high:left.high):0;
  TLegend *leg = frame(x0, x1, s.lowFor(logy), s.high, logy, xtitle, ytitle, lines, curves.size(),onRight?0.52:0.16,quiet);
  for (const auto &c : curves)
    if (!c.points) {
      TH1D *h = line(c.bins, c.color);
      h->Draw("HIST SAME");
      leg->AddEntry(h, c.text.c_str(), "l");
    }
  for (const auto &c : curves)
    if (c.points) {
      TGraphErrors *g = markers(c.bins, c.color, c.marker, c.errors, c.xerrors, c.shift);
      g->Draw("P");
      leg->AddEntry(g, c.text.c_str(), "p");
    }
  if (extra) extra();
  leg->Draw();
  save(cv, name);
}

/* A map: the header above the frame, the palette on the right. */
void drawMap(TH2 *m, const char *xtitle, const char *ytitle, const char *ztitle, const std::vector<std::string> &lines, bool logz, const std::string &name) {
  TCanvas cv("cv", "", 760, 700);
  cv.SetRightMargin(0.2);
  cv.SetTopMargin(0.09 + 0.05 * lines.size());
  cv.SetLogz(logz);
  m->GetXaxis()->SetTitle(xtitle);
  m->GetYaxis()->SetTitle(ytitle);
  m->GetZaxis()->SetTitle(ztitle);
  m->GetZaxis()->SetTitleOffset(1.3);
  if (logz) {
    double low = std::numeric_limits<double>::infinity();
    for (int b = 0; b < m->GetNcells(); ++b) if (m->GetBinContent(b) > 0) low = std::min(low, m->GetBinContent(b));
    if (std::isfinite(low)) m->SetMinimum(low);
  }
  m->Draw("COLZ");
  header(0.16, 0.95, lines);
  save(cv, name);
}

/* The x-y plane of a three-dimensional histogram at one z bin (a selection stage), detached. */
TH2D *plane(const TH1 *h, int zbin) {
  const auto x = edgesOf(h->GetXaxis()), y = edgesOf(h->GetYaxis());
  auto *m = new TH2D(fresh().c_str(), "", x.size() - 1, x.data(), y.size() - 1, y.data());
  m->SetDirectory(nullptr);
  for (int i = 1; i < static_cast<int>(x.size()); ++i)
    for (int j = 1; j < static_cast<int>(y.size()); ++j) m->SetBinContent(i, j, h->GetBinContent(i, j, zbin));
  return m;
}

/* A canvas of nx x ny panels below a title line (sPHENIX Internal and the text); returns the pad that holds them. */
TPad *panels(TCanvas &cv, int nx, int ny, const std::string &title) {
  cv.cd();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(43);
  t.SetTextSize(21);
  const double h = cv.GetWh();
  const std::string line = std::string(sPhenixInternalLabel()) + "      " + title;
  t.DrawLatex(0.012, 1 - 23 / h, line.c_str()); gVisibleText.push_back(line);
  const std::string note = gLegacy ? "Reference inputs (development)" : gPopulation == "PARTIAL" ? "Partial input population" :
      gPopulation == "UNKNOWN" ? "Input coverage unverified" : "";
  if (!note.empty()) {
    t.SetTextSize(16); t.SetTextColor(kGray+2); t.DrawLatex(0.012, 1-44/h, note.c_str()); gVisibleText.push_back(note);
  }
  auto *body = new TPad(fresh().c_str(), "", 0, 0, 1, 1 - (note.empty() ? 30 : 51) / h);
  body->Draw();
  body->cd();
  body->Divide(nx, ny, 0.001, 0.001);
  return body;
}

/* The frame of one panel, with fonts and margins in pixels so every panel of a figure reads alike; the title (bold)
   sits above the frame at its left edge. */
TH1F *panelFrame(double x0, double y0, double x1, double y1, const char *xtitle, const char *ytitle, const std::string &title = "",
                 bool logy = false, double right = 16) {
  const double w = gPad->GetWw() * gPad->GetAbsWNDC(), h = gPad->GetWh() * gPad->GetAbsHNDC();
  gPad->SetLeftMargin(80 / w);
  gPad->SetRightMargin(right / w);
  gPad->SetBottomMargin(62 / h);
  gPad->SetTopMargin(34 / h);
  gPad->SetLogy(logy);
  TH1F *f = gPad->DrawFrame(x0, y0, x1, y1);
  f->GetXaxis()->SetTitle(xtitle);
  f->GetYaxis()->SetTitle(ytitle);
  for (TAxis *a : {f->GetXaxis(), f->GetYaxis()}) {
    a->SetLabelFont(43), a->SetLabelSize(16), a->SetTitleFont(43), a->SetTitleSize(18);
  }
  f->GetXaxis()->SetTitleOffset(1.4);
  f->GetYaxis()->SetTitleOffset(1.9);
  f->GetYaxis()->SetNdivisions(505);
  if (!title.empty()) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(63);
    t.SetTextSize(16);
    t.DrawLatex(gPad->GetLeftMargin(), 1 - 24 / h, title.c_str());
  }
  return f;
}

void panelText(double x, double y, const std::string &text, int size = 17, bool bold = true, int color = cBlack) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(bold ? 63 : 43);
  t.SetTextSize(size);
  t.SetTextColor(color);
  t.DrawLatex(x, y, text.c_str());
}

TLegend *panelLegend(double x1, double y1, double x2, double y2) {
  auto *l = new TLegend(x1, y1, x2, y2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(43);
  l->SetTextSize(15);
  return l;
}

/* A figure with a ratio panel below: the upper (70%) and lower (30%) pads. */
std::pair<TPad *, TPad *> ratioPads(TCanvas &cv) {
  cv.cd();
  auto *top = new TPad(fresh().c_str(), "", 0, 0.3, 1, 1), *bottom = new TPad(fresh().c_str(), "", 0, 0, 1, 0.3);
  top->SetBottomMargin(0.02);
  bottom->SetTopMargin(0.02);
  bottom->SetBottomMargin(0.35);
  top->Draw();
  bottom->Draw();
  return {top, bottom};
}

TH1F *ratioFrame(double x0, double x1, double y0, double y1, const char *xtitle, const char *ytitle) {
  TH1F *f = gPad->DrawFrame(x0, y0, x1, y1);
  f->GetXaxis()->SetTitle(xtitle);
  f->GetYaxis()->SetTitle(ytitle);
  f->GetYaxis()->SetTitleSize(0.1), f->GetYaxis()->SetLabelSize(0.09), f->GetYaxis()->SetTitleOffset(0.65);
  f->GetXaxis()->SetTitleSize(0.12), f->GetXaxis()->SetLabelSize(0.1), f->GetXaxis()->SetTitleOffset(1.2);
  f->GetYaxis()->SetNdivisions(504);
  return f;
}

bool anyFinite(const Bins &b) {
  for (double v : b.v) if (std::isfinite(v)) return true;
  return false;
}

/* ---- photon selection ------------------------------------------------------------------------------------------ */

/* NPB (non-collision background) classifier score of the candidates before the preselection, 16-35 GeV. */
void ncbScore(const std::string &pkg) {
  const std::string name = "photonOnly/" + (pkg == "PP" ? std::string("pp") : "auau_" + spanTag()) + "_ncb_score_data_mc.png";
  auto of = [&](const char *family) { return overClasses(pkg, [&](const std::string &c) { return candidates(pkg, family, c, "ncb"); }); };
  TH1 *data = of("data");
  if (!data) return missing(name, "no DATA QA");
  const double low = 16, high = photonSupport(pkg).second;
  auto spectrum = [&](const TH1 *h) { return unit(within(project(h, 0, {{1, low, high}}), 0, 1e300)); };
  std::vector<Curve> curves;
  if (TH1 *h = of("phosim")) curves.push_back(lineCurve(spectrum(h), cRed, "Photon+jet simulation"));
  if (TH1 *h = of("jetsim")) curves.push_back(lineCurve(spectrum(h), cBlue, "Inclusive-jet simulation"));
  curves.push_back(pointCurve(spectrum(data), cBlack, 20, "DATA"));
  drawCurves(curves, 0, 1, true, "NCB classifier score", "Normalized entries",
             {(pkg == "PP" ? std::string("p+p, ") : spanLabel() + ", ") + energyOf(pkg), photonWindowLabel(low,high)}, name,
             [] { vline(0.5, gFrameLow, gFrameHigh); });
}

/* Cluster time - MBD time against cluster eta of the p+p DATA candidates before the NPB requirement. */
void ncbEtaTime() {
  const std::string name = "photonOnly/pp_ncb_eta_time_before.png";
  TH1 *h = overClasses("PP", [](const std::string &c) { return candidates("PP", "data", c, "etatime"); });
  if (!h) return missing(name, "no p+p DATA QA");
  drawMap(static_cast<TH2 *>(h), "Cluster #eta", "Cluster time #minus MBD time [ns]", "Clusters",
          {"p+p, " + energyOf("PP"), Form("%g #leq E_{T}^{cluster} < %g GeV, before NCB preselection",
            selectionValue("PP","qa","et_min_gev",10),selectionValue("PP","qa","et_max_gev",40))}, true, name);
}

/* Shower shapes and the photon-ID score after the preselection, 15-20 GeV: DATA against the prompt photons of the
   photon simulation and the inclusive-jet candidates, one figure per Au+Au class (dis_<shape>_cent<k>_pt0_cut1). */
void showerShapes() {
  const std::vector<std::array<const char *, 3>> shapes = {
      {"e11_to_e33", "e11_over_e33", "E_{11}/E_{33}"}, {"e32_to_e35", "e32_over_e35", "E_{32}/E_{35}"},
      {"bdt", "bdt", "photon-ID score s_{BDT}"}, {"weta_cogx", "weta_cogx", "w_{#eta}^{COG}"},
      {"wphi_cogx", "wphi_cogx", "w_{#phi}^{COG}"}, {"et1", "et1", "E_{T,1} fraction"}, {"et2", "et2", "E_{T,2} fraction"},
      {"et3", "et3", "E_{T,3} fraction"}, {"et4", "et4", "E_{T,4} fraction"}};
  const auto classes = classesOf("AuAu");
  for (size_t k = 0; k < classes.size(); ++k)
    for (const auto &[slot, var, xtitle] : shapes) {
      const std::string name = Form("shower_shapes/dis_%s_cent%zu_pt0_cut1.pdf", slot, k), c = classes[k];
      TH1 *data = qa("AuAu/data/photons/" + c + "/data/shape_" + var);
      if (!data) { missing(name, "no Au+Au DATA QA"); continue; }
      const double low = 15, high = 20;
      auto spectrum = [&](const TH1 *h) { return unit(project(h, 0, {{1, low, high}, {2, 1, 2}})); };
      std::vector<Curve> curves;
      if (TH1 *h = qa("AuAu/phosim/photons/" + c + "/sig/shape_" + var)) curves.push_back(lineCurve(spectrum(h), cRed, "Photon+jet simulation (prompt)"));
      if (TH1 *h = qa("AuAu/jetsim/photons/" + c + "/bkg/shape_" + var)) curves.push_back(lineCurve(spectrum(h), cBlue, "Inclusive-jet simulation"));
      curves.push_back(pointCurve(spectrum(data), cBlack, 20, "DATA"));
      const Bins &b = curves.back().bins;
      drawCurves(curves, b.edges.front(), b.edges.back(), false, xtitle, "Normalized counts",
                 {classLabel(c) + ", " + energyOf("AuAu"), photonWindowLabel(low,high)+", "+photonEtaLabel("AuAu"), "after preselection"}, name);
    }
}

/* Cluster widths of the preselected prompt photons (photon simulation) and inclusive-jet candidates per Au+Au class. */
void widths() {
  for (const auto &[label, family, who] : {std::make_tuple("sig", "phosim", "signal"), std::make_tuple("bkg", "jetsim", "background")})
    for (const auto &[var, slot, xtitle] : {std::make_tuple("weta_cogx", "full_cluster", "w_{#eta}^{COG}"),
                                            std::make_tuple("weta33_cogx", "3x3", "w_{#eta}^{3#times3, COG}")}) {
      const std::string name = std::string("photonOnly/auau_") + who + "_" + slot + "_width.png";
      TH1 *h = qa(std::string("AuAu/") + family + "/slices/" + label + "/width_" + var);
      if (!h) { missing(name, std::string("no Au+Au ") + family + " simulation QA"); continue; }
      std::vector<Curve> curves;
      const int colors[] = {cBlue, cOrange, cGreen, cRed, cPurple};
      int i = 0;
      const auto domain = photonSupport("AuAu");
      for (const auto &c : classesOf("AuAu"))
        curves.push_back(lineCurve(unit(project(h, 0, {{1, domain.first, domain.second}, {2, lowOf(c) / 5, highOf(c) / 5}})), colors[i++ % 5], percent(lowOf(c), highOf(c))));
      drawCurves(curves, 0, 0.7, false, xtitle, "Unit-normalized candidates / 0.01",
                 {"Au+Au, " + energyOf("AuAu"), std::string(label == std::string("sig") ? "Truth-matched prompt photons" : "Inclusive-jet candidates") + ", preselected",
                  photonWindowLabel(domain.first,domain.second)}, name);
    }
}

/* The photon-ID score of the preselected prompt photons and inclusive-jet candidates, one panel per Au+Au class. */
void bdtPerClass() {
  const std::string name = "photonOnly/auau_bdt_score_perCent.png";
  TH1 *sig = qa("AuAu/phosim/slices/sig/bdt"), *bkg = qa("AuAu/jetsim/slices/bkg/bdt");
  const auto classes = classesOf("AuAu");
  if (!sig || !bkg || classes.empty()) return missing(name, "needs Au+Au photon and inclusive simulation QA");
  std::vector<std::pair<Bins, Bins>> spectra;
  const auto domain = photonSupport("AuAu");
  Span s;
  for (const auto &c : classes) {
    const std::vector<Cut> w = {{1, domain.first, domain.second}, {2, lowOf(c) / 5, highOf(c) / 5}};
    spectra.push_back({unit(project(sig, 0, w)), unit(project(bkg, 0, w))});
    s.add(spectra.back().first, false, 0, 1), s.add(spectra.back().second, false, 0, 1);
  }
  TCanvas cv("cv", "", 520 * classes.size(), 560);
  TPad *body = panels(cv, classes.size(), 1, "Au+Au SIM, " + energyOf("AuAu") + ", " + photonWindowLabel(domain.first, domain.second) +
    ", "+photonEtaLabel("AuAu")+Form(", #left|z_{vtx}#right| < %g cm, preselected, no isolation",selectionValue("AuAu","event","vertex_abs_z_max_cm",10)));
  /* Place the legend above the right-hand score range instead of inflating the entire y axis by a fixed factor. */
  double underneath = 0;
  for (const auto &[signal, background] : spectra) for (const auto &b : {signal, background})
    for (size_t i=0; i<b.v.size(); ++i) if (b.edges[i+1] > 0.5) underneath = std::max(underneath, b.v[i]);
  const double floor = s.positive/2, top = std::max(s.high*1.12, floor*std::pow(std::max(underneath/floor, 1.0), 1.0/0.70));
  for (size_t k = 0; k < classes.size(); ++k) {
    body->cd(k + 1);
    panelFrame(0, floor, 1, top, "BDT score, s_{BDT}", k == 0 ? "Fraction of candidates" : "", classLabel(classes[k]), true);
    TH1D *a = line(spectra[k].first, cBlue), *b = line(spectra[k].second, cOrange);
    a->Draw("HIST SAME");
    b->Draw("HIST SAME");
    if (k == 0) {
      TLegend *l = panelLegend(0.50, 0.74, 0.97, 0.9);
      l->AddEntry(a, "Prompt photons", "l");
      l->AddEntry(b, "Inclusive-jet candidates", "l");
      l->Draw();
    }
  }
  save(cv, name);
}

/* E11/E33 through the selection (before and after the preselection, after the tight ID), 15-35 GeV, one row per
   class; DATA with cluster time - MBD time in [-40, -5) ns (NPB timing tag, before the preselection) scaled to its panel. */
void e11Flow() {
  const std::string name = "shower_shapes/e11_e33_selection_flow.png";
  std::vector<std::tuple<std::string, std::string, std::string>> rows;
  for (const auto &c : classesOf("AuAu")) rows.push_back({"AuAu", c, classLabel(c)});
  for (const auto &c : classesOf("PP")) rows.push_back({"PP", c, "p+p reference"});
  bool drawn = false;
  for (const auto &[pkg, c, text] : rows) drawn = drawn || qa(pkg + "/data/photons/" + c + "/data/shape_e11_over_e33");
  if (!drawn) return missing(name, "no DATA QA");
  TCanvas cv("cv", "", 1260, 40 + 310 * rows.size());
  TPad *body = panels(cv, 3, rows.size(), "E_{11}/E_{33}: DATA and SIM selection stages");
  const char *titles[] = {"Before preselection", "After preselection", "After tight ID"};
  for (size_t r = 0; r < rows.size(); ++r) {
    const auto &[pkg, c, text] = rows[r];
    const auto domain = photonSupport(pkg);
    TH1 *data = qa(pkg + "/data/photons/" + c + "/data/shape_e11_over_e33");
    for (int stage = 0; stage < 3; ++stage) {
      body->cd(r * 3 + stage + 1);
      const std::vector<Cut> w = {{1, domain.first, domain.second}, {2, double(stage), stage + 1.0}};
      std::vector<std::pair<Bins, int>> lines;
      if (TH1 *h = qa(pkg + "/phosim/photons/" + c + "/sig/shape_e11_over_e33")) lines.push_back({unit(project(h, 0, w)), cRed});
      if (TH1 *h = qa(pkg + "/jetsim/photons/" + c + "/bkg/shape_e11_over_e33")) lines.push_back({unit(project(h, 0, w)), cBlue});
      Bins d, tagged;
      double total = 0;
      if (data) {
        const Bins raw = project(data, 0, w);
        total = raw.sum();
        d = unit(raw);
      }
      if (data && stage == 0)
        if (TH1 *t = qa(pkg + "/data/photons/" + c + "/data/e11e33_dt")) tagged = project(t, 0, {{1, domain.first, domain.second}, {2, -40, -5}});
      Span s;
      for (const auto &l : lines) s.add(l.first, false, 0, 1);
      s.add(d, true, 0, 1);
      panelFrame(0, 0, 1, 1.35 * (std::isfinite(s.high) ? s.high : 1), r + 1 == rows.size() ? "E_{11}/E_{33}" : "",
                 r == rows.size() / 2 && stage == 0 ? "Fraction of candidates" : "", r == 0 ? titles[stage] : "");
      for (const auto &l : lines) line(l.first, l.second)->Draw("HIST SAME");
      if (!tagged.v.empty() && total > 0) {
        line(scaled(tagged, 1 / total), cGreen)->Draw("HIST SAME");
        panelText(gPad->GetLeftMargin() + 0.03, 1 - gPad->GetTopMargin() - 0.17, Form("NCB timing tag: %.0f", tagged.sum()), 13, false, cGreen);
      }
      if (data) markers(d, cBlack, 20, true, false, 0, 0.6)->Draw("P");
      if (stage == 0) {
        panelText(gPad->GetLeftMargin() + 0.03, 1 - gPad->GetTopMargin() - 0.08, text, 15);
        panelText(gPad->GetLeftMargin() + 0.03, 1 - gPad->GetTopMargin() - 0.13, photonWindowLabel(domain.first,domain.second), 13, false);
      }
      if (r == 0 && stage == 2) {
        TLegend *l = panelLegend(gPad->GetLeftMargin() + 0.02, 0.5, 0.98, 1 - gPad->GetTopMargin() - 0.02);
        l->SetTextSize(14);
        l->AddEntry(markers(Bins{{0, 1}, {0}, {0}}, cBlack, 20, false, false, 0, 0.6), "DATA", "p");
        l->AddEntry(line(Bins{{0, 1}, {0}, {0}}, cRed), "Photon+jet simulation", "l");
        l->AddEntry(line(Bins{{0, 1}, {0}, {0}}, cBlue), "Inclusive-jet simulation", "l");
        l->AddEntry(line(Bins{{0, 1}, {0}, {0}}, cGreen), "NCB timing tag (scaled)", "l");
        l->Draw();
      }
    }
  }
  save(cv, name);
}

/* Truth isolation of the prompt photons of the photon simulation: the fraction below a threshold, direct and
   fragmentation photons, 15-20 and 25-30 GeV; one panel per class. */
void truthIsolation() {
  const std::string name = "photonOnly/pp_auau_truth_isolation_direct_fragmentation_2x2.png";
  const auto order = panelOrder();
  std::vector<TH1 *> hs;
  bool any = false;
  for (const auto &[pkg, c] : order) hs.push_back(qa(pkg + "/phosim/truth/" + c + "/truth_iso")), any = any || hs.back();
  if (!any) return missing(name, "no photon-family simulation QA");
  const int nx = 2, ny = (order.size() + 1) / 2;
  TCanvas cv("cv", "", 1100, 40 + 400 * ny);
  TPad *body = panels(cv, nx, ny, "Photon+jet SIM; "+truthSelectionLabel(order.front().first));
  for (size_t p = 0; p < order.size(); ++p) {
    body->cd(p + 1);
    const bool bottom = static_cast<int>(p) >= (ny - 1) * nx, left = p % nx == 0;
    panelFrame(0.5, 0, 10.5, 1.1, bottom ? "Truth isolation threshold E_{T,cut}^{iso} [GeV]" : "",
               left ? "Fraction with E_{T}^{iso,truth} < E_{T,cut}^{iso}" : "", classLabel(order[p].second) + ", " + energyOf(order[p].first));
    vline(selectionValue(order[p].first,"truth_photon","isolation_max_gev",4), 0, 1.1);
    if (!hs[p]) continue;
    TLegend *l = p == 0 ? panelLegend(0.5, 0.2, 0.97, 0.5) : nullptr;
    for (const auto &[kind, color, who] : {std::make_tuple(1, cBlue, "Direct"), std::make_tuple(2, cOrange, "Fragmentation")})
      for (const auto &[lo, hi, marker] : {std::make_tuple(15.0, 20.0, 20), std::make_tuple(25.0, 30.0, 21)}) {
        const Bins v = project(hs[p], 0, {{1, lo, hi}, {2, double(kind), kind + 1.0}});
        const double total = v.sum();
        if (!(total > 0)) continue;
        auto *g = new TGraphErrors;
        for (int cut = 1; cut <= 10; ++cut) {
          double below = 0;
          for (int i = 0; i < v.n(); ++i) if (v.edges[i + 1] <= cut + 1e-9) below += v.v[i];
          const double f = below / total;
          g->SetPoint(cut - 1, cut, f);
          g->SetPointError(cut - 1, 0, std::sqrt(std::max(f * (1 - f), 0.0) / std::max(total, 1.0)));
        }
        style(g, color, marker, 0.9);
        g->Draw("P");
        if (l) l->AddEntry(g, Form("%s, %g-%g GeV", who, lo, hi), "p");
      }
    if (l) l->Draw();
  }
  save(cv, name);
}

/* The prompt photons below 4 GeV of truth isolation against truth pT: all, direct, fragmentation; one per class. */
void truthPtAfterIsolation() {
  for (const auto &[pkg, c] : figureOrder()) {
    const std::string name = "photonOnly/" + (pkg == "PP" ? std::string("pp") : "auau_" + classTag(c)) + "_truth_photon_pt_after_isolation.png";
    TH1 *h = qa(pkg + "/phosim/truth/" + c + "/truth_iso");
    if (!h) { missing(name, "no photon-family simulation QA"); continue; }
    const double isoMax = selectionValue(pkg,"truth_photon","isolation_max_gev",4);
    const auto domain = photonSupport(pkg);
    const Bins direct = within(project(h, 1, {{0, 0, isoMax}, {2, 1, 2}}), domain.first, domain.second),
               fragmentation = within(project(h, 1, {{0, 0, isoMax}, {2, 2, 3}}), domain.first, domain.second);
    const Bins all = sum(direct, fragmentation);
    TCanvas cv("cv", "", 700, 700);
    auto [top, bottom] = ratioPads(cv);
    top->cd();
    Span s;
    for (const Bins *b : {&all, &direct, &fragmentation}) s.add(perWidth(*b), true, domain.first, domain.second);
    TLegend *leg = frame(domain.first, domain.second, s.positive, s.high, true, "", "Weighted photons / GeV",
                         {classLabel(c) + ", " + energyOf(pkg), Form("0 #leq E_{T,truth}^{iso} < %g GeV; ",isoMax)+truthSelectionLabel(pkg)}, 3);
    gFrame->GetXaxis()->SetLabelSize(0);
    for (const auto &[b, color, marker, text] : {std::make_tuple(&all, cBlack, 20, "All prompt photons"), std::make_tuple(&direct, cRed, 21, "Direct"),
                                                 std::make_tuple(&fragmentation, cBlue, 22, "Fragmentation")}) {
      TGraphErrors *g = markers(perWidth(*b), color, marker);
      g->Draw("P");
      leg->AddEntry(g, text, "p");
    }
    leg->Draw();
    bottom->cd();
    ratioFrame(domain.first, domain.second, 0, 1.1, "Truth photon p_{T} [GeV]", "Direct / total");
    markers(over(direct, all), cRed, 21, false)->Draw("P");
    save(cv, name);
  }
}

/* Reconstructed isolation of Au+Au DATA candidates at 14-16 GeV in five centrality ranges, with and without the
   underlying-event subtraction. */
void isolationOverlays() {
  for (const auto &[kind, stem, note] : {std::make_tuple("noUE", "isonoueslice", "without UE subtraction"), std::make_tuple("withUE", "isoslice", "with SUB1 UE subtraction")}) {
    const std::string name = std::string("isolation/isoCentOverlay_") + kind + "_pT_14_16.png";
    TH1 *h = qa(std::string("AuAu/data/slices/data/") + stem);
    if (!h) { missing(name, "no Au+Au DATA QA"); continue; }
    std::vector<Curve> curves;
    const double low = 14, high = 16;
    for (const auto &[lo, hi, color] : {std::make_tuple(0, 2, cRed), std::make_tuple(2, 4, cBlue), std::make_tuple(4, 8, cGreen),
                                         std::make_tuple(8, 12, cPurple), std::make_tuple(12, 16, cOrange)}) {
      Bins b = unit(rebin(project(h, 0, {{1, low, high}, {2, double(lo), double(hi)}}), 5));
      for (int i = 0; i < b.n(); ++i) if (!(b.v[i] > 0)) b.v[i] = std::nan("");
      curves.push_back(pointCurve(b, color, 20, percent(5 * lo, 5 * hi)));
    }
    drawCurves(curves, -10, 40, false, "E_{T}^{iso} [GeV]", "Normalized to unit area",
               {"Au+Au, " + energyOf("AuAu"), photonWindowLabel(low,high)+Form(", R_{iso} = %g",isolationRadius("AuAu")), note}, name);
  }
}

/* Reconstructed isolation of the tight DATA candidates, 15-35 GeV, against the non-tight DATA scaled to it above
   6 GeV and the prompt photons of the photon simulation scaled to the remainder; one per class. */
void isolationStacks() {
  const double sideband = 6.0;
  for (const auto &[pkg, c] : figureOrder()) {
    const std::string name = "isolation/" + (pkg == "PP" ? std::string("pp_nominal_topo") : "auau_" + classTag(c)) + "_isoStack.png";
    TH1 *data = qa(pkg + "/data/photons/" + c + "/data/iso"), *sig = qa(pkg + "/phosim/photons/" + c + "/sig/iso");
    if (!data || !sig) { missing(name, "needs DATA and photon-family simulation QA"); continue; }
    const auto domain = photonSupport(pkg);
    const Bins tight = project(data, 0, {{1, domain.first, domain.second}, {2, 0, 1}}), nonTight = project(data, 0, {{1, domain.first, domain.second}, {2, 1, 2}});
    Bins signal = project(sig, 0, {{1, domain.first, domain.second}, {2, 0, 1}});
    double tightSide = 0, nonTightSide = 0;
    for (int i = 0; i < tight.n(); ++i) if (tight.edges[i] >= sideband) tightSide += tight.v[i], nonTightSide += nonTight.v[i];
    const Bins background = scaled(nonTight, nonTightSide > 0 ? tightSide / nonTightSide : 0.0);
    signal = scaled(signal, signal.sum() > 0 ? (tight.sum() - background.sum()) / signal.sum() : 0.0);
    const int k = std::lround(1.0 / tight.width(0));
    const Bins t = rebin(tight, k), b = rebin(background, k), stacked = sum(rebin(signal, k), b);
    const double x0 = pkg == "PP" ? -5 : -20;
    Span s;
    s.add(t, true, x0, 20), s.add(stacked, false, x0, 20);
    TCanvas cv("cv", "", 700, 600);
    TLegend *leg = frame(x0, 20, s.low, s.high, false, "E_{T}^{iso,reco} [GeV]", "Counts / 1 GeV",
                         {classLabel(c) + ", " + energyOf(pkg), photonWindowLabel(domain.first,domain.second)+", "+photonEtaLabel(pkg),
                          Form("non-tight normalized at E_{T}^{iso} #geq %g GeV", sideband)}, 3);
    TH1D *all = line(stacked, cBlue, 1, 2), *under = line(b, TColor::GetColor("#F4B6B6"), 1, 1);
    all->SetFillColor(TColor::GetColor("#B8C0EE"));
    under->SetFillColor(TColor::GetColor("#F4B6B6"));
    all->Draw("HIST SAME");
    under->Draw("HIST SAME");
    line(stacked, cBlue, 1, 2)->Draw("HIST SAME");
    TGraphErrors *g = markers(t, cBlack, 20);
    g->Draw("P");
    gPad->RedrawAxis();
    leg->AddEntry(under, "Non-tight data control", "f");
    leg->AddEntry(all, "Signal MC", "f");
    leg->AddEntry(g, "Tight-photon data", "p");
    leg->Draw();
    save(cv, name);
  }
}

/* DATA event-leading candidates per region (A-D) and GeV on the ABCD photon bins; one panel per class. */
void abcdCounts() {
  const std::string name = "photonOnly/pp_auau_abcd_region_counts_1x4.png";
  const auto order = panelOrder();
  std::vector<TH1 *> hs;
  Span s;
  for (const auto &[pkg, c] : order) {
    hs.push_back(qa(pkg + "/abcd/" + c));
    if (hs.back())
      for (int r = 0; r < 4; ++r) s.add(perWidth(project(hs.back(), 1, {{0, double(r), r + 1.0}})), true, 15, 35);
  }
  if (!std::isfinite(s.high)) return missing(name, "no DATA analysis fills");
  TCanvas cv("cv", "", 420 * order.size(), 480);
  TPad *body = panels(cv, order.size(), 1, "DATA event-leading photon candidates per ABCD region");
  const std::array<std::pair<int, const char *>, 4> regions = {{{cBlack, "A: tight, isolated"}, {cRed, "B: tight, non-isolated"},
                                                                {cBlue, "C: non-tight, isolated"}, {cGreen, "D: non-tight, non-isolated"}}};
  for (size_t p = 0; p < order.size(); ++p) {
    body->cd(p + 1);
    panelFrame(15, s.positive / 3, 35, s.high * 400, "E_{T}^{#gamma,reco} [GeV]", p == 0 ? "Candidates / GeV" : "", classLabel(order[p].second), true);
    if (!hs[p]) continue;
    TLegend *l = p + 1 == order.size() ? panelLegend(0.36, 0.66, 0.97, 1 - gPad->GetTopMargin() - 0.02) : nullptr;
    for (int r = 0; r < 4; ++r) {
      TGraphErrors *g = markers(perWidth(project(hs[p], 1, {{0, double(r), r + 1.0}})), regions[r].first, 20, true, true, 0, 0.9);
      g->Draw("P");
      if (l) l->AddEntry(g, regions[r].second, "p");
    }
    if (l) l->Draw();
  }
  save(cv, name);
}

/* Two rows of panels per class: DATA or simulated purity above, a second quantity below (the purity figures). */
void purityPanels(const std::string &name, const std::string &title, const char *what, const char *upperTitle, double upperHigh,
                  const char *lowerTitle, double lowerLow, double lowerHigh) {
  const auto order = panelOrder();
  const bool closure = std::string(what) == "closure";
  bool any = false;
  for (const auto &[pkg, c] : order)
    if (TH1 *raw = qa(pkg + "/" + what + "/" + c + "/raw")) any = any || anyFinite(binsOf(raw));
  if (!any) return missing(name, closure ? "needs both simulation families and the responses" : "needs DATA analysis fills and the simulation responses");
  TCanvas cv("cv", "", 420 * order.size(), 820);
  TPad *body = panels(cv, order.size(), 2, title);
  for (size_t p = 0; p < order.size(); ++p) {
    const std::string dir = order[p].first + "/" + what + "/" + order[p].second + "/";
    TH1 *raw = qa(dir + "raw"), *corrected = qa(dir + "corrected");
    const bool drawn = raw && anyFinite(binsOf(raw));
    body->cd(p + 1);
    panelFrame(15, 0, 35, upperHigh, "", p == 0 ? upperTitle : "", classLabel(order[p].second));
    TLegend *upper = p == 0 ? panelLegend(0.35, 0.08 + gPad->GetBottomMargin(), 0.97, 0.3 + gPad->GetBottomMargin()) : nullptr;
    if (drawn) {
      TGraphErrors *a = markers(binsOf(raw), cBlack, 33, false, true, 0, 1.3), *b = corrected ? markers(binsOf(corrected), cBlue, 27, false, true, 0, 1.3) : nullptr;
      a->Draw("P");
      if (b) b->Draw("P");
      if (upper) upper->AddEntry(a, "Raw ABCD", "p");
      if (upper && b) upper->AddEntry(b, "Leakage-corrected ABCD", "p");
      if (closure)
        if (TH1 *truth = qa(dir + "truth")) {
          TGraphErrors *t = markers(binsOf(truth), cRed, 20, false, true);
          t->Draw("P");
          if (upper) upper->AddEntry(t, "Truth purity", "p");
        }
    }
    if (upper) upper->Draw();
    body->cd(order.size() + p + 1);
    panelFrame(15, lowerLow, 35, lowerHigh, "E_{T}^{#gamma,reco} [GeV]", p == 0 ? lowerTitle : "");
    if (!drawn) continue;
    if (closure) {
      hline(1, 15, 35);
      if (TH1 *truth = qa(dir + "truth"))
        if (corrected) markers(over(binsOf(corrected), binsOf(truth)), cBlue, 33, false, true, 0, 1.3)->Draw("P");
    } else {
      TLegend *lower = p == 0 ? panelLegend(0.3, 0.62, 0.8, 0.9) : nullptr;
      for (const auto &[r, color, text] : {std::make_tuple("B", cBlack, "S_{B}/S_{A}"), std::make_tuple("C", cRed, "S_{C}/S_{A}"), std::make_tuple("D", cBlue, "S_{D}/S_{A}")})
        if (TH1 *h = qa(dir + "leakage_" + r)) {
          TGraphErrors *g = markers(binsOf(h), color, 20, false, true, 0, 0.9);
          g->Draw("P");
          if (lower) lower->AddEntry(g, text, "p");
        }
      if (lower) lower->Draw();
    }
  }
  save(cv, name);
}

/* Photon efficiency of the analysis-signal truth photons of the photon simulation against truth pT. */
void efficiency() {
  const std::string name = "photonOnly/pp_si_di_auau_photon_efficiency_1x4.png";
  const auto order = panelOrder();
  std::vector<TH1 *> hs;
  bool any = false;
  for (const auto &[pkg, c] : order) hs.push_back(qa(pkg + "/phosim/truth/" + c + "/eff")), any = any || hs.back();
  if (!any) return missing(name, "no photon-family simulation QA");
  TCanvas cv("cv", "", 420 * order.size(), 480);
  TPad *body = panels(cv, order.size(), 1, "Photon+jet simulation, analysis-signal truth photons");
  const int magenta = TColor::GetColor("#C2185B");
  for (size_t p = 0; p < order.size(); ++p) {
    body->cd(p + 1);
    panelFrame(15, 0, 35, 1.1, "Truth photon p_{T} [GeV]", p == 0 ? "Photon efficiency" : "", classLabel(order[p].second));
    if (!hs[p]) continue;
    std::vector<Bins> stage;
    for (int k = 0; k < 6; ++k) stage.push_back(within(project(hs[p], 0, {{1, double(k), k + 1.0}}), 15, 35));
    Bins reco = over(stage[1], stage[0]), idIso = over(stage[5], stage[3]), selected = over(stage[5], stage[0]);
    Bins both = reco;
    for (int i = 0; i < both.n(); ++i) both.v[i] = reco.v[i] * idIso.v[i];
    TLegend *l = p == 0 ? panelLegend(gPad->GetLeftMargin() + 0.03, 0.08 + gPad->GetBottomMargin(), 0.9, 0.34 + gPad->GetBottomMargin()) : nullptr;
    for (const auto &[b, color, marker, text] : {std::make_tuple(&reco, cBlack, 21, "#varepsilon_{reco}"),
                                                 std::make_tuple(&both, magenta, 20, "#varepsilon_{reco} #varepsilon_{ID|iso}"),
                                                 std::make_tuple(&selected, cGreen, 34, "#varepsilon_{sel}")}) {
      TGraphErrors *g = markers(*b, color, marker, false, false, 0, 1.0);
      g->Draw("P");
      if (l) l->AddEntry(g, text, "p");
    }
    if (l) l->Draw();
  }
  save(cv, name);
}

/* Total calorimeter energy against centrality of the embedded photon simulation, before and after the minimum-bias
   classifier. */
void minimumBias() {
  const std::string name = "photonOnly/auau_minbias_pretraining_before_after.png";
  TH1 *h = qa("AuAu/phosim/events/calo_total");
  if (!h) return missing(name, "no Au+Au photon-family simulation QA");
  const auto energy = edgesOf(h->GetXaxis()), centrality = edgesOf(h->GetYaxis());
  TH2D *maps[2];
  for (int m = 0; m < 2; ++m) {
    maps[m] = new TH2D(fresh().c_str(), "", centrality.size() - 1, centrality.data(), energy.size() - 1, energy.data());
    maps[m]->SetDirectory(nullptr);
    for (int e = 1; e < static_cast<int>(energy.size()); ++e)
      for (int c = 1; c < static_cast<int>(centrality.size()); ++c)
        for (int state = 1; state <= 4; ++state)
          if (m == 0 || state == 2 || state == 4) maps[m]->AddBinContent(maps[m]->GetBin(c, e), h->GetBinContent(e, c, state));
  }
  const double top = maps[0]->GetMaximum();
  TCanvas cv("cv", "", 1300, 600);
  TPad *body = panels(cv, 2, 1, "Au+Au photon+jet SIM, " + energyOf("AuAu"));
  const char *titles[] = {"Before: no classifier requirement", "After: minimum-bias pass"};
  for (int m = 0; m < 2; ++m) {
    body->cd(m + 1);
    panelFrame(0, energy.front(), 80, energy.back(), "Centrality percentile [%]", m == 0 ? "E^{total}_{calo} [GeV]" : "", titles[m], false, m == 1 ? 110 : 16);
    gPad->SetLogz();
    maps[m]->SetMinimum(top * 1e-5);
    maps[m]->SetMaximum(top);
    maps[m]->GetZaxis()->SetTitle("Weighted events / bin");
    maps[m]->GetZaxis()->SetLabelFont(43), maps[m]->GetZaxis()->SetLabelSize(15);
    maps[m]->GetZaxis()->SetTitleFont(43), maps[m]->GetZaxis()->SetTitleSize(17), maps[m]->GetZaxis()->SetTitleOffset(2.2);
    maps[m]->Draw(m == 1 ? "COLZ SAME" : "COL SAME");
    gPad->RedrawAxis();
  }
  save(cv, name);
}

/* ---- events and samples ---------------------------------------------------------------------------------------- */

/* The DATA centrality of the package trigger (the centrality weights' target) against the embedded simulation before
   the centrality weights, in 2% bins of 0-80%. */
void centralityTarget() {
  const std::string name = "simulation_stitching/auau_centrality_data_embeddedMC.png";
  TH1 *data = qa("AuAu/centrality_target");
  TH1 *sims[2] = {qa("AuAu/phosim/events/centrality"), qa("AuAu/jetsim/events/centrality")};
  if (!data && !sims[0] && !sims[1]) return missing(name, "no DATA target events and no simulation QA");
  auto twoPercent = [](const Bins &b) { return unit(within(rebin(b, std::max(1L, std::lround(2.0 / b.width(0)))), 0, 80)); };
  std::vector<Curve> curves;
  if (data) curves.push_back(pointCurve(twoPercent(binsOf(data)), cBlack, 33, "DATA (centrality-weight target)", false));
  const std::array<std::tuple<int, int, const char *, double>, 2> looks = {{{cBlue, 21, "Photon+jet SIM", 0.0}, {cOrange, 22, "Inclusive-jet SIM", 0.3}}};
  for (int f = 0; f < 2; ++f)
    if (sims[f]) {
      Curve c = pointCurve(twoPercent(project(sims[f], 0, {{1, 0, 1}})), std::get<0>(looks[f]), std::get<1>(looks[f]), std::get<2>(looks[f]), false);
      c.shift = std::get<3>(looks[f]);
      curves.push_back(c);
    }
  drawCurves(curves, 0, 80, false, "Centrality percentile [%]", "Fraction per 2% bin",
             {"Au+Au, " + energyOf("AuAu"), "simulation before centrality reweighting"}, name);
}

/* Au+Au event properties (all classes together): DATA against the photon and inclusive-jet simulation. */
void eventProperties() {
  const std::string span = spanLabel() + ", " + energyOf("AuAu");
  struct Look { std::string family; int color; std::string text; };
  const std::vector<Look> sims = {{"phosim", cRed, "Photon+jet simulation"}, {"jetsim", cBlue, "Inclusive-jet simulation"}};
  /* centrality: DATA with w_event, the simulation with the full weight (after the centrality weights) */
  if (TH1 *data = qa("AuAu/data/events/centrality")) {
    std::vector<Curve> curves;
    for (const auto &s : sims)
      if (TH1 *h = qa("AuAu/" + s.family + "/events/centrality")) curves.push_back(lineCurve(unit(project(h, 0, {{1, 1, 2}})), s.color, s.text + " (reweighted)"));
    curves.push_back(pointCurve(unit(project(data, 0, {{1, 0, 1}})), cBlack, 20, "DATA", false));
    drawCurves(curves, 0, 100, false, "Centrality [%]", "Normalized events", {"Au+Au, " + energyOf("AuAu")}, "event_obs/centrality.pdf");
  } else missing("event_obs/centrality.pdf", "no Au+Au DATA QA");
  /* vertex */
  if (TH1 *data = overClasses("AuAu", [](const std::string &c) { return qa("AuAu/data/events/" + c + "/evt_vz"); })) {
    std::vector<Curve> curves;
    for (const auto &s : sims)
      if (TH1 *h = overClasses("AuAu", [&](const std::string &c) { return qa("AuAu/" + s.family + "/events/" + c + "/evt_vz"); }))
        curves.push_back(lineCurve(unit(binsOf(h)), s.color, s.text));
    curves.push_back(pointCurve(unit(binsOf(data)), cBlack, 20, "DATA", false));
    drawCurves(curves, -30, 30, false, "v_{z} [cm]", "Normalized events", {span}, "event_obs/vertexz.pdf");
  } else missing("event_obs/vertexz.pdf", "no Au+Au DATA QA");
  /* cluster ET of the preselected candidates; the simulation scaled to the DATA */
  if (TH1 *data = overClasses("AuAu", [](const std::string &c) { return candidates("AuAu", "data", c, "clu_et"); })) {
    const Bins d = project(data, 0, {{1, 1, 2}});
    std::vector<Curve> curves;
    for (const auto &s : sims)
      if (TH1 *h = overClasses("AuAu", [&](const std::string &c) { return candidates("AuAu", s.family, c, "clu_et"); })) {
        const Bins b = project(h, 0, {{1, 1, 2}});
        curves.push_back(lineCurve(scaled(b, b.sum() > 0 ? d.sum() / b.sum() : 0.0), s.color, s.text));
      }
    curves.push_back(pointCurve(d, cBlack, 20, "DATA"));
    drawCurves(curves, 5, 50, true, "Cluster E_{T} [GeV]", "Clusters", {span, "preselected candidates"}, "event_obs/cluster_pt.pdf");
  } else missing("event_obs/cluster_pt.pdf", "no Au+Au DATA QA");
  /* cluster and jet maps */
  for (const auto &m : std::vector<Look>{{"data", 0, "DATA"}, {"phosim", 0, "Photon+jet simulation"}}) {
    const std::string tag = m.family == "data" ? "data" : "photon";
    const std::string clusters = "event_obs/cluster_eta_phi_" + tag + ".pdf", jets = "event_obs/jet_eta_phi_" + tag + ".pdf";
    if (TH1 *h = overClasses("AuAu", [&](const std::string &c) { return candidates("AuAu", m.family, c, "clu_etaphi"); }))
      drawMap(plane(h, 2), "Cluster #eta", "Cluster #phi [rad]", "Clusters", {m.text + ", " + spanLabel(), "E_{T}^{cluster} > 10 GeV, preselected"}, false, clusters);
    else missing(clusters, "no Au+Au " + m.family + " QA");
    if (TH1 *h = overClasses("AuAu", [&](const std::string &c) { return qa("AuAu/" + m.family + "/events/" + c + "/jet_etaphi"); }))
      drawMap(plane(h, 2), "Jet #eta", "Jet #phi [rad]", "Jets", {m.text + ", " + spanLabel(),
        Form("p_{T}^{jet} #geq %g GeV, events with a Region-A photon",selectionValue("AuAu","qa","jet_map_pt_min_gev",5))}, false, jets);
    else missing(jets, "no Au+Au " + m.family + " QA");
  }
  /* jets of the events with a Region-A photon; the simulation scaled to the DATA */
  if (TH1 *data = overClasses("AuAu", [](const std::string &c) { return qa("AuAu/data/events/" + c + "/jet_pt"); })) {
    const Bins d = project(data, 0, {{1, 1, 2}});
    std::vector<Curve> curves;
    for (const auto &s : sims)
      if (TH1 *h = overClasses("AuAu", [&](const std::string &c) { return qa("AuAu/" + s.family + "/events/" + c + "/jet_pt"); })) {
        const Bins b = project(h, 0, {{1, 1, 2}});
        curves.push_back(lineCurve(scaled(b, b.sum() > 0 ? d.sum() / b.sum() : 0.0), s.color, s.text));
      }
    curves.push_back(pointCurve(d, cBlack, 20, "DATA"));
    drawCurves(curves, 0, 60, true, "Jet p_{T} [GeV]", "Jets", {span,
      Form("anti-k_{T} R = %g, events with a Region-A photon",selectionValue("AuAu","jet","radius",0.4))}, "event_obs/jet_pt.pdf");
  } else missing("event_obs/jet_pt.pdf", "no Au+Au DATA QA");
}

/* The stitched simulated spectrum of a family: every sample, their sum and its descriptive fit, and sum / fit. */
void stitching(const std::string &pkg, const std::string &family, bool generator) {
  const bool pp = pkg == "PP", photons = family == "phosim";
  const std::string name = std::string("simulation_stitching/") + (pp ? "pp" : "auau") + "_" + (photons ? "photonjet" : "inclusivejet") + "_" +
                           (generator ? "generator" : "reco_tag") + "_stitching.png";
  TDirectory *d = gQA->GetDirectory(("qa/" + pkg + "/stitching/" + family + (generator ? "_generator" : "_reco")).c_str());
  auto *total = d ? dynamic_cast<TH1 *>(d->Get("total")) : nullptr;
  if (!total) return missing(name, "no " + pkg + " " + family + " simulation QA");
  const double lo = generator ? (photons ? 15 : 10) : 15, hi = generator ? (photons ? 40 : 60) : 35;
  std::vector<std::pair<std::string, Bins>> samples;
  TIter next(d->GetListOfKeys());
  while (auto *k = static_cast<TKey *>(next()))
    if (std::string(k->GetName()).rfind("sample_", 0) == 0) samples.push_back({std::string(k->GetName()).substr(7), within(perWidth(binsOf(static_cast<TH1 *>(d->Get(k->GetName())))), lo, hi)});
  std::sort(samples.begin(), samples.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
  const Bins all = within(perWidth(binsOf(total)), lo, hi);
  /* An empty stitched sum is listed, not drawn as an empty canvas. */
  int first = -1, last = -1;
  for (int i = 0; i < all.n(); ++i)
    if (all.v[i] > 0) { if (first < 0) first = i; last = i; }
  if (first < 0) return missing(name, "no " + pkg + " " + family + " entries in the displayed range");
  auto *fit = dynamic_cast<TH1 *>(d->Get("fit"));
  auto *rms = dynamic_cast<TParameter<double> *>(d->Get("rms"));
  TF1 *model = nullptr;
  if (fit && std::isfinite(fit->GetBinContent(1))) {  /* the descriptive fit, drawn only over the bins it describes */
    model = new TF1(fresh().c_str(), "exp([0] - ([1] + [2]*log(x) + [3]*x)*log(x))", all.edges[first], all.edges[last + 1]);
    for (int p = 0; p < 4; ++p) model->SetParameter(p, fit->GetBinContent(p + 1));
    model->SetLineColor(cBlack), model->SetLineStyle(2), model->SetLineWidth(1), model->SetNpx(200);
  }
  TCanvas cv("cv", "", 700, 700);
  auto [top, bottom] = ratioPads(cv);
  top->cd();
  Span s;
  s.add(all, true, lo, hi);
  for (const auto &sample : samples) s.add(sample.second, true, lo, hi);
  const std::string central = classesOf(pkg).empty() ? "" : classesOf(pkg).front();
  TLegend *leg = frame(lo, hi, s.positive, s.high, true, "", "Weighted spectrum / GeV",
                       {(pp ? "p+p SIM, " : "Au+Au SIM, ") + energyOf(pkg), photons ? "Photon+jet samples" : "Inclusive-jet samples",
                        generator ? "Generator truth" : "Leading tight, isolated photon" + (pp ? std::string() : ", " + percent(lowOf(central), highOf(central)))},
                       samples.size() + 1);
  gFrame->GetXaxis()->SetLabelSize(0);
  const int colors[] = {cOrange, cPurple, cGreen, cRed, cBlue, kOrange + 3, kCyan + 2, kSpring + 5};
  const int open[] = {24, 25, 26, 27, 28, 30, 32, 46};
  for (size_t i = 0; i < samples.size(); ++i) {
    TGraphErrors *g = markers(samples[i].second, colors[i % 8], open[i % 8], true, false, 0, 1.0);
    g->Draw("P");
    leg->AddEntry(g, samples[i].first.c_str(), "p");
  }
  TGraphErrors *g = markers(all, cBlack, 33, true, false, 0, 1.3);
  g->Draw("P");
  leg->AddEntry(g, "Stitched sum", "p");
  if (model) model->Draw("SAME");
  leg->Draw();
  bottom->cd();
  ratioFrame(lo, hi, 0.5, 1.5, generator ? (photons ? "Maximum truth photon p_{T} [GeV]" : "Maximum truth jet p_{T} [GeV]") : "Reconstructed photon E_{T} [GeV]",
             "Spectrum / fit");
  hline(1, lo, hi);
  if (model) {
    Bins r = all;
    for (int i = 0; i < r.n(); ++i) {
      const double m = model->Eval(r.center(i));
      r.v[i] = all.v[i] > 0 ? all.v[i] / m : std::nan(""), r.e2[i] = all.e2[i] / (m * m);
    }
    markers(r, cBlack, 20, true, false, 0, 0.8)->Draw("P");
    if (rms && std::isfinite(rms->GetVal())) {
      TLatex t;
      t.SetNDC();
      t.SetTextSize(0.1);
      t.SetTextAlign(31);
      t.DrawLatex(0.93, 0.8, Form("RMS %.1f%%", 100 * rms->GetVal()));
    }
  }
  save(cv, name);
}

/* The photon-ID score (bdt) or isolation (isolation) thresholds at 70, 80 and 90% signal efficiency: per 5% slice
   against photon pT with the slice medians, and the medians against centrality with their lines (the qa step's). */
void workingPoints(const std::string &what) {
  const bool bdt = what == "bdt";
  const double radius = isolationRadius("AuAu");
  const std::string radiusTag = Form("%02d",int(std::round(radius*10)));
  const std::string grid = bdt ? "photonOnly/auau_h70_bdt_pt_centrality_grid.png" : "isolation/auau_r"+radiusTag+"_isolation_pt_centrality_grid.png";
  const std::string fits = bdt ? "photonOnly/auau_h70_bdt_centrality_fits.png" : "isolation/auau_r"+radiusTag+"_isolation_centrality_fits.png";
  TDirectory *d = gQA->GetDirectory(("qa/AuAu/working_points/" + what).c_str());
  if (!d) { missing(grid, "no Au+Au photon-family simulation QA"); return missing(fits, "no Au+Au photon-family simulation QA"); }
  const std::string ytitle = bdt ? "BDT score threshold" : Form("R = %g isolation threshold [GeV]",radius);
  const std::array<std::pair<int, int>, 3> levels = {{{70, cBlue}, {80, cGreen}, {90, cOrange}}};
  auto *first = dynamic_cast<TH2 *>(d->Get("thresholds_70"));
  const int slices = first ? first->GetNbinsX() : 0;
  if (!slices) { missing(grid,"missing threshold grid"); return missing(fits,"missing threshold grid"); }
  const double ptLow = first->GetYaxis()->GetXmin(), ptHigh = first->GetYaxis()->GetXmax();
  auto *centrality = dynamic_cast<TH1 *>(d->Get("median_70"));
  if (!centrality || centrality->GetNbinsX()!=slices) throw std::runtime_error("working-point centrality axes differ");
  /* per slice */
  {
    TCanvas cv("cv", "", 1400, 40 + 290 * ((slices + 3) / 4));
    TPad *body = panels(cv, 4, (slices + 3) / 4, std::string("Au+Au photon+jet SIM, ") + ytitle + " at 70, 80, 90% signal efficiency");
    const int rows = (slices + 3) / 4;
    for (int row = 0; row < rows; ++row) {
      Span s;
      for (const auto &[level, color] : levels)
        if (auto *t = dynamic_cast<TH2 *>(d->Get(Form("thresholds_%d", level))))
          for (int x = 4 * row + 1; x <= std::min(slices, 4 * row + 4); ++x)
            for (int y = 1; y <= t->GetNbinsY(); ++y) s.add(t->GetBinContent(x, y));
      const double pad = std::isfinite(s.high) ? 0.15 * std::max(s.high - s.low, 1e-3) : 1;
      for (int col = 0; col < 4 && 4 * row + col < slices; ++col) {
        const int x = 4 * row + col + 1;
        body->cd(x);
        panelFrame(ptLow, (std::isfinite(s.low) ? s.low : 0) - pad, ptHigh, (std::isfinite(s.high) ? s.high : 1) + 2.5 * pad,
                   row + 1 == rows ? "Photon E_{T}^{reco} [GeV]" : "", col == 0 ? (bdt ? "Score threshold" : "Threshold [GeV]") : "",
                   percent(centrality->GetXaxis()->GetBinLowEdge(x),centrality->GetXaxis()->GetBinUpEdge(x)));
        TLegend *l = x == 4 ? panelLegend(0.55, 0.62, 0.97, 0.88) : nullptr;
        for (const auto &[level, color] : levels) {
          auto *t = dynamic_cast<TH2 *>(d->Get(Form("thresholds_%d", level)));
          auto *median = dynamic_cast<TH1 *>(d->Get(Form("median_%d", level)));
          if (!t) continue;
          auto *g = new TGraphErrors;
          for (int y = 1; y <= t->GetNbinsY(); ++y)
            if (std::isfinite(t->GetBinContent(x, y))) g->SetPoint(g->GetN(), t->GetYaxis()->GetBinCenter(y), t->GetBinContent(x, y));
          style(g, color, 20, 0.7);
          if (g->GetN()) g->Draw("P");
          if (median && std::isfinite(median->GetBinContent(x))) {
            TLine m;
            m.SetLineColor(color), m.SetLineStyle(2), m.SetLineWidth(1);
            m.DrawLine(ptLow, median->GetBinContent(x), ptHigh, median->GetBinContent(x));
          }
          if (l) l->AddEntry(g, Form("%d%%", level), "p");
        }
        if (l) l->Draw();
      }
    }
    save(cv, grid);
  }
  /* the slice medians against centrality */
  std::vector<Curve> curves;
  std::vector<std::pair<TF1 *, int>> lines;
  for (const auto &[level, color] : levels) {
    auto *median = dynamic_cast<TH1 *>(d->Get(Form("median_%d", level)));
    auto *fit = dynamic_cast<TH1 *>(d->Get(Form("line_%d", level)));
    if (!median) continue;
    const double a = fit ? fit->GetBinContent(1) : std::nan(""), b = fit ? fit->GetBinContent(2) : std::nan("");
    curves.push_back(pointCurve(binsOf(median), color, 20, std::isfinite(a) ? Form("%d%%: %.4f %+.5f C", level, a, b) : Form("%d%%", level), false));
    if (std::isfinite(a)) {
      auto *f = new TF1(fresh().c_str(), "[0] + [1]*x", centrality->GetXaxis()->GetXmin(), centrality->GetXaxis()->GetXmax());
      f->SetParameters(a, b);
      f->SetLineColor(color), f->SetLineWidth(2);
      lines.push_back({f, color});
    }
  }
  drawCurves(curves, centrality->GetXaxis()->GetXmin(), centrality->GetXaxis()->GetXmax(), false, "Centrality percentile [%]", ytitle.c_str(),
             {"Au+Au, " + energyOf("AuAu")+", "+photonWindowLabel(ptLow,ptHigh), "QA threshold fits; descriptive lines"}, fits,
             [&] { for (auto &l : lines) l.first->Draw("SAME"); });
}

/* The maximum-cluster-energy spectra of the reference and photon triggers and the photon-trigger turn-ons against the
   reference, over the runs in which every trigger of the study is enabled (the qa step's selection). */
void triggerStudy(const std::string &pkg) {
  const bool pp = pkg == "PP";
  const std::string spectra = pp ? "event_obs/pp_max_cluster_energy_overlay.png" : "event_obs/hMaxClusterEnergy_groupOverlay.png";
  const std::string turnOn = pp ? "event_obs/pp_max_cluster_energy_turnon.png" : "event_obs/hMaxClusterEnergy_groupTurnOnOverlay.png";
  TDirectory *d = gQA->GetDirectory(("qa/" + pkg + "/triggers").c_str());
  auto *study = d ? dynamic_cast<TObjString *>(d->Get("study")) : nullptr;
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
  TIter next(d->GetListOfKeys());
  while (auto *k = static_cast<TKey *>(next()))
    if (std::string(k->GetName()).rfind("spectrum_", 0) == 0) photons.push_back(std::string(k->GetName()).substr(9));
  std::sort(photons.begin(), photons.end(), [](const std::string &a, const std::string &b) {
    return std::atof(a.c_str() + a.find_first_of("0123456789")) < std::atof(b.c_str() + b.find_first_of("0123456789"));
  });
  const std::vector<std::string> lines = {(pp ? "p+p, " : "Au+Au, ") + energyOf(pkg), Form("scaled-trigger study, %d runs", runs)};
  const int colors[] = {cRed, cBlue, cGreen, cPurple};
  std::vector<Curve> curves = {lineCurve(binsOf(static_cast<TH1 *>(d->Get("reference"))), cBlack, info["reference"])};
  for (size_t i = 0; i < photons.size(); ++i) curves.push_back(lineCurve(binsOf(static_cast<TH1 *>(d->Get(("spectrum_" + photons[i]).c_str()))), colors[i % 4], photons[i]));
  const Bins &reference = curves.front().bins;
  drawCurves(curves, reference.edges.front(), reference.edges.back(), true, "Maximum cluster energy [GeV]", "Scaled-trigger events", lines, spectra);
  const Bins denominator = binsOf(static_cast<TH1 *>(d->Get("reference_denominator")));
  std::vector<Curve> efficiencies;
  for (size_t i = 0; i < photons.size(); ++i) {
    Bins eff = binsOf(static_cast<TH1 *>(d->Get(("numerator_" + photons[i]).c_str())));
    double crossing = std::nan("");
    for (int b = 0; b < eff.n(); ++b) {
      const double n = denominator.v[b], e = n > 0 ? eff.v[b] / n : std::nan("");
      eff.v[b] = e, eff.e2[b] = n > 0 ? std::max(e * (1 - e), 0.0) / n : 0;
      if (std::isnan(crossing) && std::isfinite(e) && e >= 0.9) crossing = eff.center(b);
    }
    efficiencies.push_back(pointCurve(eff, colors[i % 4], 20, photons[i] + "/MBD" + (std::isfinite(crossing) ? Form(" (90%%: %.1f GeV)", crossing) : "")));
  }
  drawCurves(efficiencies, 0, pp ? 20 : 30, false, "Maximum cluster energy [GeV]", "Trigger / MBD efficiency", lines, turnOn,
             [pp] { hline(1, 0, pp ? 20 : 30); });
}

/* Draw the supported QA families from saved quantities. Individual helpers report missing inputs through missing(). */
void qaFigures() {
  for (const std::string pkg : {"AuAu", "PP"}) triggerStudy(pkg);
  for (const std::string pkg : {"PP", "AuAu"})
    for (const std::string family : {"phosim", "jetsim"})
      for (bool generator : {true, false}) stitching(pkg, family, generator);
  centralityTarget();
  eventProperties();
  minimumBias();
  ncbScore("PP");
  ncbScore("AuAu");
  ncbEtaTime();
  showerShapes();
  widths();
  workingPoints("bdt");
  bdtPerClass();
  e11Flow();
  truthIsolation();
  truthPtAfterIsolation();
  isolationOverlays();
  workingPoints("isolation");
  isolationStacks();
  abcdCounts();
  purityPanels("photonOnly/data_purity_and_prompt_leakage_2x4_style.png", "DATA photon purity and the prompt-signal leakage of the photon simulation",
               "purity", "DATA photon purity", 1.2, "Prompt-signal leakage fraction", 0, 0.6);
  purityPanels("photonOnly/purity_qa_pp_auau_simCLOSURE.png", "Purity closure: every simulated candidate as pseudo-data", "closure",
               "Photon purity (simulation)", 1.3, "Corrected ABCD / truth", 0.6, 1.6);
  efficiency();
}

}  /* namespace */
