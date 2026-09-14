#ifndef CALOANA_RJFASTJETALGOSUBAREA_H
#define CALOANA_RJFASTJETALGOSUBAREA_H

#include <jetbase/FastJetOptions.h>
#include <jetbase/Jet.h>
#include <jetbase/JetAlgo.h>
#include <jetbase/JetContainer.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

// FastJetAlgoSub preserves negative-energy tower constituents by mapping them
// to positive 1 MeV pseudojets during clustering and then restoring the
// original constituent four-vectors in the output jet.  The release class does
// not implement FastJetOptions::calc_area.  This project-local extension keeps
// that established clustering behavior byte-for-byte at the algorithm level
// while adding FastJet active-area ghosts and the standard prop_area field.
class RJFastJetAlgoSubArea final : public JetAlgo
{
 public:
  explicit RJFastJetAlgoSubArea(const FastJetOptions& options)
    : m_opt(options)
  {
    fastjet::ClusterSequence banner;
    if (m_opt.verbosity > 0)
    {
      fastjet::ClusterSequence::print_banner();
    }
    else
    {
      std::ostringstream sink;
      fastjet::ClusterSequence::set_fastjet_banner_stream(&sink);
      fastjet::ClusterSequence::print_banner();
      fastjet::ClusterSequence::set_fastjet_banner_stream(&std::cout);
    }
  }

  ~RJFastJetAlgoSubArea() override = default;

  void identify(std::ostream& os = std::cout) override
  {
    os << "   RJFastJetAlgoSubArea: ";
    if (m_opt.algo == Jet::ANTIKT) os << "ANTIKT";
    else if (m_opt.algo == Jet::KT) os << "KT";
    else if (m_opt.algo == Jet::CAMBRIDGE) os << "CAMBRIDGE";
    else os << "UNKNOWN";
    os << " r=" << m_opt.jet_R
       << " active_area=" << (m_opt.calc_area ? 1 : 0) << std::endl;
  }

  Jet::ALGO get_algo() override { return m_opt.algo; }
  float get_par() override { return m_opt.jet_R; }

  void cluster_and_fill(std::vector<Jet*>& particles,
                        JetContainer* jetcont) override
  {
    if (!jetcont) return;
    initializeContainer(jetcont);

    std::vector<fastjet::PseudoJet> pseudojets;
    pseudojets.reserve(particles.size());
    for (std::size_t ipart = 0; ipart < particles.size(); ++ipart)
    {
      Jet* particle = particles[ipart];
      if (!particle) continue;

      float energy = particle->get_e();
      if (energy == 0.0F) continue;
      float px = particle->get_px();
      float py = particle->get_py();
      float pz = particle->get_pz();

      // Preserve the production FastJetAlgoSub negative-energy contract.
      if (energy < 0.0F)
      {
        const float ratio = 0.001F / energy;
        energy *= ratio;
        px *= ratio;
        py *= ratio;
        pz *= ratio;
      }

      fastjet::PseudoJet pseudojet(px, py, pz, energy);
      pseudojet.set_user_index(static_cast<int>(ipart));
      pseudojets.push_back(pseudojet);
    }

    fastjet::JetDefinition definition = jetDefinition();
    std::unique_ptr<fastjet::ClusterSequence> plainSequence;
    std::unique_ptr<fastjet::ClusterSequenceArea> areaSequence;
    std::vector<fastjet::PseudoJet> fastjets;
    if (m_opt.calc_area)
    {
      const fastjet::AreaDefinition areaDefinition(
          fastjet::active_area_explicit_ghosts,
          fastjet::GhostedAreaSpec(
              m_opt.ghost_max_rap, 1, m_opt.ghost_area));
      areaSequence = std::make_unique<fastjet::ClusterSequenceArea>(
          pseudojets, definition, areaDefinition);
      fastjets = areaSequence->inclusive_jets();
    }
    else
    {
      plainSequence = std::make_unique<fastjet::ClusterSequence>(
          pseudojets, definition);
      fastjets = plainSequence->inclusive_jets();
    }

    unsigned int outputId = 0;
    for (const fastjet::PseudoJet& fastjet : fastjets)
    {
      if (m_opt.calc_area && fastjet.is_pure_ghost()) continue;

      Jet* jet = jetcont->add_jet();
      if (!jet) continue;

      float totalPx = 0.0F;
      float totalPy = 0.0F;
      float totalPz = 0.0F;
      float totalEnergy = 0.0F;
      float weightedTime = 0.0F;
      float weightedEnergy = 0.0F;

      for (const fastjet::PseudoJet& constituent : fastjet.constituents())
      {
        if (m_opt.calc_area && constituent.is_pure_ghost()) continue;
        const int inputIndex = constituent.user_index();
        if (inputIndex < 0 ||
            inputIndex >= static_cast<int>(particles.size())) continue;
        Jet* particle = particles[static_cast<std::size_t>(inputIndex)];
        if (!particle) continue;

        totalPx += particle->get_px();
        totalPy += particle->get_py();
        totalPz += particle->get_pz();
        totalEnergy += particle->get_e();
        if (particle->size_properties() > Jet::PROPERTY::prop_t &&
            !std::isnan(particle->get_property(Jet::PROPERTY::prop_t)))
        {
          weightedTime +=
              particle->get_property(Jet::PROPERTY::prop_t) *
              particle->get_e();
          weightedEnergy += particle->get_e();
        }
        jet->insert_comp(particle->get_comp_vec(), true);
      }

      if (jet->size_properties() < Jet::PROPERTY::prop_t + 1)
        jet->resize_properties(Jet::PROPERTY::prop_t + 1);
      jet->set_property(
          Jet::PROPERTY::prop_t, weightedTime / weightedEnergy);
      if (m_opt.calc_area)
        jet->set_property(m_areaIndex, fastjet.area());

      jet->set_comp_sort_flag();
      jet->set_px(totalPx);
      jet->set_py(totalPy);
      jet->set_pz(totalPz);
      jet->set_e(totalEnergy);
      jet->set_id(outputId++);
    }
  }

 private:
  void initializeContainer(JetContainer* jetcont)
  {
    if (m_initialized) return;
    m_opt.initialize();
    jetcont->set_algo(m_opt.algo);
    jetcont->set_jetpar_R(m_opt.jet_R);
    if (m_opt.calc_area)
    {
      jetcont->add_property(Jet::PROPERTY::prop_area);
      m_areaIndex =
          jetcont->property_index(Jet::PROPERTY::prop_area);
    }
    m_initialized = true;
  }

  fastjet::JetDefinition jetDefinition() const
  {
    if (m_opt.algo == Jet::KT)
      return fastjet::JetDefinition(
          fastjet::kt_algorithm, m_opt.jet_R,
          fastjet::E_scheme, fastjet::Best);
    if (m_opt.algo == Jet::CAMBRIDGE)
      return fastjet::JetDefinition(
          fastjet::cambridge_algorithm, m_opt.jet_R,
          fastjet::E_scheme, fastjet::Best);
    return fastjet::JetDefinition(
        fastjet::antikt_algorithm, m_opt.jet_R,
        fastjet::E_scheme, fastjet::Best);
  }

  FastJetOptions m_opt{};
  bool m_initialized{false};
  Jet::PROPERTY m_areaIndex{Jet::PROPERTY::no_property};
};

#endif
