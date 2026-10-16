/* Shared reconstructed-energy variations for measurement and response.
   Counter-based draws preserve original source/entry/object identities.
   Truth-based resolutions use the existing exact reco-truth matches. */
#pragma once
#include "PhotonJet.h"
namespace photonjet { namespace response {
/* A variation of the reco energy scale and resolution, for the systematic uncertainties: every reco jet of the analysis collection and
   photon of the event is taken at value x scale x (1 + smear g) before any selection, g a standard Gaussian drawn per
   object from a counter-based generator keyed by (seed, the file's production input, entry, object), so the draws of
   one file do not depend on the other inputs or their order. The reco-truth jet match stays the producer's, and the
   photon classes and isolation thresholds are the stored ones. The text is "" (nominal) or whitespace-separated
   key=value pairs of jet_scale, jet_smear, photon_scale, photon_smear, seed. */
struct Variation {
  bool centralityOff = false;
  double jetScale = 1.0, jetSmear = 0.0, photonScale = 1.0, photonSmear = 0.0;
  double emScale = 1.0;
  std::string photonResolution = "off";
  std::string photonSelection = "varied";
  std::string jetResolution = "off";
  std::uint64_t seed = 0;
};

inline bool responseOnly(const Variation& v) { return v.photonSelection == "response_only"; }
inline bool emsSelectionResponse(const Variation& v) { return v.photonSelection == "ems_selection_response"; }

// Native EMS requires its explicit source/selection/model/cone contract.
// Retained-record replay remains unqualified and supplies no such context.
inline void requireRunnableVariation(const Variation& v,bool hasExplicitContract=false) {
  if(emsSelectionResponse(v) && !hasExplicitContract)
    throw std::runtime_error("EMS composition requires an explicit source-bound selection contract");
}

inline void validatePhotonSelection(const Variation& v) {
  if (v.photonSelection != "varied" && !responseOnly(v) && !emsSelectionResponse(v)) throw std::runtime_error("unknown photon selection policy");
  if (responseOnly(v) && (v.photonResolution == "off" || v.photonScale != 1. || v.emScale != 1. ||
      v.photonSmear != 0. || v.jetSmear != 0.))
    throw std::runtime_error("response-only photon resolution cannot be combined with another energy effect");
}
inline Variation parseVariation(const char* text) {
  Variation v;
  std::vector<std::string> seen;
  std::istringstream stream(text ? text : "");
  for (std::string item; stream >> item;) {
    const size_t at = item.find('=');
    if (at == std::string::npos) throw std::runtime_error("a variation item is key=value: " + item);
    const std::string key = item.substr(0, at), value = item.substr(at + 1);
    if (std::find(seen.begin(), seen.end(), key) != seen.end()) throw std::runtime_error("a variation key appears twice: " + key);
    seen.push_back(key);
    if (key == "centrality_weight") {
      if (value != "on" && value != "off") throw std::runtime_error("centrality_weight requires on or off");
      v.centralityOff = value == "off"; continue;
    }
    if (key == "photon_selection") {
      if (value != "response_only" && value != "ems_selection_response") throw std::runtime_error("unknown photon selection policy");
      v.photonSelection = value; continue;
    }
    if (key == "photon_resolution") {
      if (value!="off" && value!="ppg18_nominal" && value!="ppg18_high" && value!="ppg18_low")
        throw std::runtime_error("unknown photon resolution prescription: " + value);
      v.photonResolution=value;
      continue;
    }
    if(key=="jet_resolution") {
      if(value!="off" && value!="sphenix_pp2026_nominal" && value!="sphenix_pp2026_high" && value!="sphenix_pp2026_low")
        throw std::runtime_error("unknown jet resolution prescription: "+value);
      v.jetResolution=value;continue;
    }
    size_t used = 0;
    if (key == "seed") {
      if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("a seed is a non-negative integer: " + item);
      v.seed = std::stoull(value);
      continue;
    }
    double x = 0.0;
    try {
      x = std::stod(value, &used);
    } catch (const std::exception&) {
      used = 0;
    }
    if (used == 0 || used != value.size()) throw std::runtime_error("a variation value is one number: " + item);
    if (!std::isfinite(x)) throw std::runtime_error("a variation value must be finite: " + item);
    if (key == "em_scale") {
      if (!(x > 0.0)) throw std::runtime_error("an EM scale must be positive: " + item);
      v.emScale=x;
    } else if (key == "jet_scale" || key == "photon_scale") {
      if (!(x > 0.0)) throw std::runtime_error("a scale must be positive: " + item);
      (key == "jet_scale" ? v.jetScale : v.photonScale) = x;
    } else if (key == "jet_smear" || key == "photon_smear") {
      if (x < 0.0) throw std::runtime_error("a smear must not be negative: " + item);
      (key == "jet_smear" ? v.jetSmear : v.photonSmear) = x;
    } else {
      throw std::runtime_error("unknown variation key " + key);
    }
  }
  if (v.emScale!=1.0 && v.photonScale!=1.0) throw std::runtime_error("em_scale already scales photons; do not also set photon_scale");
  if (v.photonResolution!="off" && v.photonSmear!=0.0) throw std::runtime_error("photon_resolution cannot be combined with constant photon_smear");
  if(v.jetResolution!="off" && v.jetSmear!=0.)throw std::runtime_error("jet_resolution cannot be combined with legacy reco-scaled jet_smear");
  validatePhotonSelection(v);
  if(emsSelectionResponse(v) && (v.photonResolution=="off" || v.jetResolution=="off" || v.photonScale!=1. || v.jetScale!=1. || v.photonSmear!=0. || v.jetSmear!=0.))
    throw std::runtime_error("EMS composition requires named photon/jet resolution and em_scale only");
  const bool centralityNominalResolution = responseOnly(v) &&
      v.photonResolution == "ppg18_nominal" && v.jetResolution == "sphenix_pp2026_nominal";
  if (v.centralityOff && (v.jetScale != 1. || v.photonScale != 1. || v.emScale != 1. ||
      v.jetSmear != 0. || v.photonSmear != 0. ||
      (!centralityNominalResolution && (v.photonResolution != "off" || v.jetResolution != "off" || responseOnly(v))) ||
      emsSelectionResponse(v)))
    throw std::runtime_error("centrality off requires a pure weight alternative or the nominal response-only resolutions");
  return v;
}

inline std::string variationJson(const Variation& v) {
  std::ostringstream s;
  s.precision(17);
  s << "{\"jet_scale\":" << v.jetScale << ",\"jet_smear\":" << v.jetSmear << ",\"photon_scale\":" << v.photonScale
    << ",\"photon_smear\":" << v.photonSmear << ",\"em_scale\":" << v.emScale
    << ",\"photon_resolution\":\"" << v.photonResolution << "\",\"jet_resolution\":\"" << v.jetResolution << "\",\"seed\":" << v.seed << '}';
  if (emsSelectionResponse(v)) { auto text=s.str();text.pop_back();return text+",\"photon_selection\":\"ems_selection_response\"}"; }
  if (responseOnly(v)) {
    auto text=s.str();text.pop_back();text+=",\"photon_selection\":\"response_only\"";
    if(v.centralityOff)text+=",\"centrality_weight\":\"off\"";
    return text+'}';
  }
  if (v.centralityOff) { auto text=s.str(); text.pop_back(); return text+",\"centrality_weight\":\"off\"}"; }
  return s.str();
}

/* Centrality off is an AuAu SIM response alternative. DATA measurement and
   the frozen centrality map used to check stored tree factors stay unchanged. */
inline void validateVariationSystem(const Variation& variation, System system) {
  if (variation.centralityOff && system != System::auau)
    throw std::runtime_error("centrality_weight=off applies only to AuAu simulation responses");
}

/* Evaluate the pre-centrality factors directly. Dividing an already weighted
   event loses zero-weight events and can apply the truth-vertex fit twice.
   Source/stitching normalization remains at the package consumer's stage. */
inline double responseEventWeight(const Variation& variation, double base, double vertex, double centrality) {
  const double product = base * vertex * (variation.centralityOff ? 1.0 : centrality);
  if (strictReduction && !std::isfinite(product)) throw std::runtime_error("nonfinite response weight");
  return std::isfinite(product) ? product : 1.0;
}


inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

inline std::uint64_t textKey(const std::string& text) {  /* FNV-1a */
  std::uint64_t h = 14695981039346656037ULL;
  for (unsigned char c : text) h = (h ^ c) * 1099511628211ULL;
  return h;
}

/* Standard Gaussian of one key: Box-Muller on two 53-bit uniforms. */
inline double gaussian(std::uint64_t key) {
  const std::uint64_t a = splitmix64(key), b = splitmix64(a);
  const double unit = 1.0 / 9007199254740992.0;  /* 2^-53 */
  const double u1 = ((a >> 11) + 1.0) * unit, u2 = (b >> 11) * unit;
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

/* PPG18 section 6.4: additional Gaussian width at truth pT, in GeV.
   Keep the same object draw for nominal/high; low has no additional smear.
   The mode is explicit: historical nominal never silently changes. */
inline double photonResolutionWidth(double truthPt, const std::string &mode) {
  if (mode=="off" || mode=="ppg18_low") return 0.;
  if (!(std::isfinite(truthPt) && truthPt>0.)) throw std::runtime_error("EM resolution requires positive finite truth pT");
  const auto variance=[truthPt](double p0,double p1,double p2) {
    return p0*p0/truthPt+p1*p1/(truthPt*truthPt)+p2*p2;
  };
  const double simulation=variance(.185,0.,.040);
  double data;
  if (mode=="ppg18_nominal") data=variance(.15,.05,.05);
  else if (mode=="ppg18_high") data=variance(.13,.08,.08);
  else throw std::runtime_error("unknown photon resolution mode");
  return truthPt*std::sqrt(std::max(0.,data-simulation));
}

/* sPH-JET-2026-01 / arXiv:2606.17184v1, supplemental Table 1:
   additional independent resolution, NOT the full detector resolution.
   Width is evaluated at the matched particle jet pT. R02 uses R03 as
   stated in sPH-CONF-JET-2026-01 (June 2026), section 3.
   This is the pp detector component; AuAu UE/fragmentation are separate
   calibration inputs. An explicit mode never silently changes nominal. */
inline double jetResolutionFraction(double radius,const std::string &mode) {
  if(mode=="off")return 0.;
  const int r=std::lround(radius*10.);
  if(std::abs(radius-r/10.)>1e-9)throw std::runtime_error("JER table requires an exact calibrated jet radius");
  const std::map<int,std::pair<double,double>> table={{2,{.107,.024}},{3,{.107,.024}},{4,{.100,.018}},
    {5,{.097,.021}},{6,{.087,.015}},{7,{.087,.012}},{8,{.091,.028}}};
  if(!table.count(r))throw std::runtime_error("jet radius is absent from the JER prescription");
  const auto [nominal,uncertainty]=table.at(r);
  if(mode=="sphenix_pp2026_nominal")return nominal;
  if(mode=="sphenix_pp2026_high")return nominal+uncertainty;
  if(mode=="sphenix_pp2026_low")return nominal-uncertainty;
  throw std::runtime_error("unknown jet resolution mode");
}

/* The reco values of one collection of one entry under a variation; scale 1 without smearing keeps them exactly. */
class Varied {
 public:
  template <class Values>
  Varied(const Values& values, double scale, double smear, std::uint64_t seed, std::uint64_t source, std::uint64_t entry,
         std::uint64_t collection)
      : v_(values.GetSize()) {
    for (size_t i = 0; i < v_.size(); ++i) {
      double factor = scale;
      if (smear != 0.0) factor *= 1.0 + smear * gaussian(splitmix64(splitmix64(splitmix64(seed ^ source) ^ entry) ^ collection) ^ i);
      v_[i] = factor == 1.0 ? values[i] : values[i] * factor;
    }
  }
  size_t GetSize() const { return v_.size(); }
  double operator[](size_t i) const { return v_[i]; }
  void scaleAll(double factor) {
    if(!(std::isfinite(factor)&&factor>0.))throw std::runtime_error("invalid response scale");
    if(factor!=1.)for(auto &value:v_)value*=factor;
  }

  template <class Values>
  void scaleEmFraction(const Values &fraction,double scale) {
    if (fraction.GetSize()!=v_.size()) throw std::runtime_error("jet EM fractions are not aligned");
    if (!std::isfinite(scale) || scale<=0.) throw std::runtime_error("invalid EM scale");
    for (size_t i=0;i<v_.size();++i) {
      const double factor=1.+(scale-1.)*fraction[i];
      // Subtracted component fractions need not lie within [0,1]. A finite
      // variation can cross zero for a soft jet; preserve its signed value
      // so the ordinary positive-pT selection handles the crossing.
      if (!std::isfinite(fraction[i]) || !std::isfinite(factor) || !std::isfinite(v_[i]*factor))
        throw std::runtime_error("invalid jet EM fraction scale");
      if (factor!=1.) v_[i]*=factor;
    }
  }
  template <class Truth>
  void smearPhotons(const Truth &truth,const std::vector<long> &matches,const Variation &variation,
                    std::uint64_t source,std::uint64_t entry) {
    if (matches.size()!=v_.size()) throw std::runtime_error("photon resolution matches are not aligned");
    for (size_t i=0;i<v_.size();++i) {
      if (matches[i]<0) continue; // No invented truth energy for unmatched photons.
      if (static_cast<size_t>(matches[i])>=truth.GetSize()) throw std::runtime_error("photon resolution truth index outside array");
      const double width=photonResolutionWidth(truth[matches[i]],variation.photonResolution);
      if (width==0.) continue;
      const auto key=splitmix64(splitmix64(splitmix64(variation.seed ^ source) ^ entry) ^ 0ULL) ^ i;
      v_[i]+=width*gaussian(key);
    }
  }
  template <class Truth>
  void smearJets(const Truth &truth,const std::vector<long> &matches,const Variation &variation,
                 double radius,std::uint64_t source,std::uint64_t entry) {
    if(matches.size()!=v_.size())throw std::runtime_error("JER matches are not aligned");
    const double fraction=jetResolutionFraction(radius,variation.jetResolution);
    if(fraction==0.)return;
    for(size_t i=0;i<v_.size();++i) {
      if(matches[i]<0)continue; // No invented particle momentum for unmatched jets.
      if(static_cast<size_t>(matches[i])>=truth.GetSize())throw std::runtime_error("JER truth index outside array");
      const double pt=truth[matches[i]];
      if(!(std::isfinite(pt)&&pt>0.))throw std::runtime_error("JER requires positive finite truth pT");
      const auto key=splitmix64(splitmix64(splitmix64(variation.seed ^ source) ^ entry) ^ 1ULL) ^ i;
      v_[i]+=fraction*pt*gaussian(key);
    }
  }

 private:
  std::vector<double> v_;
};


} }
