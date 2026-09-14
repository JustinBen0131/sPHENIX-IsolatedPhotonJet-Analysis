#ifndef RJ_CEMC_TOWER_STATUS_GUARD_H
#define RJ_CEMC_TOWER_STATUS_GUARD_H

#include "RJCemcStatusRules.h"
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <cdbobjects/CDBTTree.h>
#include <ffaobjects/CdbUrlSave.h>
#include <ffaobjects/RunHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/phool.h>
#include <TFile.h>
#include <TSystem.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <tuple>
#include <utility>

namespace rj_cemc_status
{
// One state shared by the raw-input and pre-clustering calibrated checks.
// The worker packet must pin the recovery payload bytes. This code never
// searches CVMFS or silently falls back to a current CDB globaltag lookup.
struct State
{
  int expected_run = 0;
  std::string original_calo_input;
  bool recovery_authorized = false;
  std::string recovery_payload;
  std::string selected_payload;
  bool restoring_missing_status = false;
  bool initialized = false;
  std::vector<std::uint8_t> mask;
  unsigned long long raw_events = 0;
  unsigned long long calibrated_events = 0;
};

class Guard final : public SubsysReco
{
 public:
  Guard(std::shared_ptr<State> state, bool calibrated)
    : SubsysReco(calibrated ? "CEMCStatusBeforeClustering" : "CEMCStatusBeforeCalibration")
    , m_state(std::move(state)), m_calibrated(calibrated) {}

  int InitRun(PHCompositeNode* top) override
  {
    try
    {
      require(bool(m_state) && m_state->expected_run > 0, "missing run binding");
      auto* runHeader = findNode::getClass<RunHeader>(top, "RunHeader");
      require(runHeader && runHeader->get_RunNumber() == m_state->expected_run,
              "input RunHeader does not match the bound run");
      if (m_calibrated)
      {
        require(m_state->initialized, "calibrated check registered before raw check");
        return Fun4AllReturnCodes::EVENT_OK;
      }
      auto* towers = findNode::getClass<TowerInfoContainer>(top, "TOWERS_CEMC");
      require(towers && towers->size() == channel_count, "TOWERS_CEMC input missing");
      std::set<std::string> saved;
      // Paired inputs and the geometry input share a live RUN node. Its last
      // loaded CdbUrl is not necessarily the original calorimeter provenance.
      // Read only the bound source's RUN metadata into a separate node.
      require(!m_state->original_calo_input.empty(), "missing exact original calorimeter source");
      auto originalRun = std::make_unique<PHCompositeNode>("ORIGINAL_CALO_RUN");
      {
        PHNodeIOManager original(m_state->original_calo_input, PHReadOnly, PHRunTree);
        require(original.isFunctional() && original.read(originalRun.get()),
                "unreadable original calorimeter RUN metadata");
      }
      auto* originalHeader = findNode::getClass<RunHeader>(originalRun.get(), "RunHeader");
      require(originalHeader && originalHeader->get_RunNumber() == m_state->expected_run,
              "original calorimeter metadata run mismatch");
      auto* urls = findNode::getClass<CdbUrlSave>(originalRun.get(), "CdbUrl");
      if (urls)
        for (auto it = urls->begin(); it != urls->end(); ++it)
          if (std::get<0>(*it) == "CEMC_BadTowerMap") saved.insert(std::get<1>(*it));
      require(saved.size() <= 1, "ambiguous saved CEMC map provenance");
      if (saved.empty())
      {
        require(m_state->recovery_authorized && !m_state->recovery_payload.empty(),
                "no saved CEMC map; explicit validated recovery payload required");
        // Limit this repair to the demonstrated missing-map/zero-HOT case.
        for (unsigned i = 0; i < towers->size(); ++i)
        {
          auto* t = towers->get_tower_at_channel(i);
          require(t && !t->get_isHot(), "missing provenance with existing HOT state requires review");
        }
        m_state->selected_payload = m_state->recovery_payload;
        m_state->restoring_missing_status = true;
      }
      else
      {
        // An explicit repair payload must never refresh already-valid inputs.
        m_state->selected_payload = *saved.begin();
        m_state->restoring_missing_status = false;
      }
      const auto& path = m_state->selected_payload;
      const auto suffix = "_" + std::to_string(m_state->expected_run) + "cdb.root";
      require(path.find("/cdb/CEMC_BadTowerMap/") != std::string::npos &&
              path.size() > suffix.size() &&
              path.compare(path.size()-suffix.size(), suffix.size(), suffix) == 0,
              "payload path is not bound to this run/CEMC domain");
      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      require(file && !file->IsZombie(), "unreadable CEMC map file");
      file->Close();
      CDBTTree map(path);
      std::vector<int> values;
      values.reserve(channel_count);
      for (unsigned i = 0; i < towers->size(); ++i)
        values.push_back(map.GetIntValue(towers->encode_key(i), "status", 0));
      m_state->mask = validated_mask(values);
      m_state->initialized = true;
      std::cout << "CEMC_STATUS_CONTRACT_INIT run=" << m_state->expected_run
                << " original_input=" << m_state->original_calo_input
                << " payload=" << path
                << " action=" << (m_state->restoring_missing_status ? "RESTORE_MISSING_HOT_ONLY" : "PRESERVE_VALIDATE")
                << " map_rejected=" << std::count(m_state->mask.begin(), m_state->mask.end(), 1)
                << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    catch (const std::exception& e) { return fatal(e.what()); }
  }

  int process_event(PHCompositeNode* top) override
  {
    try
    {
      require(m_state && m_state->initialized, "uninitialized status guard");
      auto* towers = findNode::getClass<TowerInfoContainer>(top,
          m_calibrated ? "TOWERINFO_CALIB_CEMC" : "TOWERS_CEMC");
      require(towers && towers->size() == channel_count, "required CEMC node missing");
      std::vector<std::uint8_t> flags;
      flags.reserve(channel_count);
      for (unsigned i = 0; i < towers->size(); ++i)
      {
        auto* t = towers->get_tower_at_channel(i);
        require(t != nullptr, "null tower");
        flags.push_back(t->get_status());
      }
      if (!m_calibrated && m_state->restoring_missing_status)
      {
        add_missing_hot_flags(flags, m_state->mask);
        for (unsigned i = 0; i < towers->size(); ++i)
          if (m_state->mask[i]) towers->get_tower_at_channel(i)->set_isHot(true);
      }
      validate_flags(flags, m_state->mask);
      for (unsigned i = 0; i < towers->size(); ++i)
        require(!m_state->mask[i] ||
                (towers->get_tower_at_channel(i)->get_isHot() &&
                 !towers->get_tower_at_channel(i)->get_isGood()),
                "map-rejected tower must be HOT and not GOOD in the actual node");
      if (m_calibrated)
      {
        ++m_state->calibrated_events;
        require(m_state->calibrated_events == m_state->raw_events,
                "raw/calibrated status event ordering mismatch");
      }
      else ++m_state->raw_events;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    catch (const std::exception& e) { return fatal(e.what()); }
  }

  int End(PHCompositeNode*) override
  {
    if (m_calibrated)
    {
      if (!m_state || !m_state->initialized ||
          m_state->raw_events != m_state->calibrated_events)
        return fatal("incomplete paired CEMC status audit");
      std::cout << "CEMC_STATUS_CONTRACT_END run=" << m_state->expected_run
                << " raw_events=" << m_state->raw_events
                << " calibrated_events=" << m_state->calibrated_events
                << " PASS" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  static int fatal(const std::string& reason)
  {
    std::cerr << "CEMC_STATUS_CONTRACT_FATAL " << reason << std::endl;
    // Existing steering can ignore se->run's return. Nonzero process exit is
    // mandatory: an ABORTRUN alone must not become a published partial tree.
    gSystem->Exit(1);
    return Fun4AllReturnCodes::ABORTRUN;
  }
  std::shared_ptr<State> m_state;
  bool m_calibrated;
};
}
#endif
