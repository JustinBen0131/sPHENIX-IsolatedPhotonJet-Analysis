// ProductionConditions.cc
//
// Which external conditions and immutable payloads does this job bind?
// Owns CDB setup, resolved payloads and digest checks.
// Does not own the reconstruction module graph.
// Called through Production.h before PhotonJetTree captures the event.
//
#include "Production.h"

#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>
#include <calobase/TowerInfoDefs.h>
#include <phool/recoConsts.h>
#include <TDirectory.h>
#include <TFile.h>

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace photonjet
{
namespace production
{
using detail::fail;
using detail::lower;
using detail::trim;
using detail::validSha256;
using detail::fileSha256;
using detail::requireStatus;
using detail::kArchivedDoubleInteractionRun;
using detail::kCemcChannelCount;

namespace
{
bool fileReadable(const std::string& path)
{
  std::ifstream input(path, std::ios::binary);
  return input.good();
}

bool safeHashPath(const std::string& path)
{
  return !path.empty() &&
         std::all_of(path.begin(), path.end(), [](const unsigned char c) {
           return std::isalnum(c) || c == '/' || c == '.' || c == '_' ||
                  c == '-' || c == '+';
         });
}

void verifySha256(const std::string& path, const std::string& expected)
{
  if (path.empty())
  {
    fail("cannot verify an empty payload path");
  }
  if (!fileReadable(path))
  {
    fail("payload is not readable: " + path);
  }
  if (!validSha256(expected))
  {
    fail("invalid expected SHA256 for: " + path);
  }

  const std::string observed = fileSha256(path);
  if (observed != lower(expected))
  {
    fail("SHA256 mismatch for " + path + " expected=" + expected +
         " observed=" + observed);
  }
}

// The CEMC zero-suppression cross calibration interprets a stored ratio of
// zero as unity and a missing field as NaN. Validate every channel of the
// pinned payload before the calibrator consumes it.
unsigned validateZeroSuppressionPayload(const std::string& path)
{
  requireStatus(!path.empty() && fileReadable(path), "missing zero-suppression payload");
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  requireStatus(file && !file->IsZombie(), "unreadable zero-suppression payload");
  file->Close();

  CDBTTree payload(path);
  unsigned zeroSentinels = 0;
  for (unsigned channel = 0; channel < kCemcChannelCount; ++channel)
  {
    const unsigned key = TowerInfoDefs::encode_emcal(channel);
    const float ratio = payload.GetFloatValue(key, "ratio", 0);
    requireStatus(std::isfinite(ratio) && ratio >= 0,
                  "missing or negative zero-suppression ratio at channel " +
                      std::to_string(channel));
    zeroSentinels += ratio == 0;
  }
  return zeroSentinels;
}
}  // namespace

namespace detail
{
// The file actually resolved at run time is hashed, never a label or a path
// that merely claims to be the payload. Plan loading uses the same hash route.
bool validSha256(const std::string& digest)
{
  return digest.size() == 64 &&
         std::all_of(digest.begin(), digest.end(),
                     [](const unsigned char c) { return std::isxdigit(c) != 0; });
}

std::string fileSha256(const std::string& path)
{
  if (!safeHashPath(path))
  {
    fail("unsafe path supplied for SHA256 verification: " + path);
  }

#ifdef __APPLE__
  const std::string command = "/usr/bin/shasum -a 256 -- " + path;
#else
  const std::string command = "sha256sum -- " + path;
#endif

  FILE* pipe = ::popen(command.c_str(), "r");
  if (!pipe)
  {
    fail("failed to launch the SHA256 tool for: " + path);
  }

  char buffer[256] = {};
  const bool readOk = std::fgets(buffer, sizeof(buffer), pipe) != nullptr;

  errno = 0;
  const int status = ::pclose(pipe);
  const int closeErrno = errno;

  if (!readOk)
  {
    fail("failed to read the SHA256 result for: " + path);
  }

  std::string digest(buffer);
  digest = lower(trim(digest.substr(0, digest.find_first_of(" \t\r\n"))));

  if (!validSha256(digest))
  {
    fail("invalid SHA256 result for: " + path);
  }

  // ROOT or Fun4All can reap the short-lived child before pclose waits for
  // it. That is not a hash failure; the digest comparison still decides.
  const bool externallyReaped = status == -1 && closeErrno == ECHILD;
  if (status != 0 && !externallyReaped)
  {
    fail("SHA256 tool failed for: " + path);
  }

  return digest;
}

void requireStatus(const bool condition, const std::string& reason)
{
  if (!condition)
  {
    throw std::runtime_error("CEMC status contract: " + reason);
  }
}
}  // namespace detail

// ConfigureConditions
//
void ConfigureConditions(Plan& plan)
{
  Conditions& cond = plan.conditions;
  const Reconstruction& reco = plan.reconstruction;

  recoConsts* rc = recoConsts::instance();
  const int run = reco.isArchivedDoubleInteraction ? kArchivedDoubleInteractionRun : plan.job.run;
  rc->set_IntFlag("RUNNUMBER", run);
  rc->set_uint64Flag("TIMESTAMP", cond.cdbTimestamp);
  if (!cond.cdbGlobalTag.empty())
  {
    rc->set_StringFlag("CDB_GLOBALTAG", cond.cdbGlobalTag);
  }

  if (reco.reconstructCentrality && cond.centralityLocalPayloads)
  {
    verifySha256(cond.centralityDivisions, cond.centralityDivisionsSha256);
    verifySha256(cond.centralityRunScale, cond.centralityRunScaleSha256);
    verifySha256(cond.centralityVertexScale, cond.centralityVertexScaleSha256);
  }

  if (cond.cemcStatusRecoveryAuthorized)
  {
    verifySha256(cond.cemcStatusRecoveryPayload, cond.cemcStatusRecoveryPayloadSha256);
  }

  if (cond.cemcZeroSuppressionCrossCalibration)
  {
    cond.cemcZeroSuppressionPayload = CDBInterface::instance()->getUrl("CEMC_ZSCrossCalib");
    const unsigned zeroSentinels = validateZeroSuppressionPayload(cond.cemcZeroSuppressionPayload);
    std::cout << "CEMC_ZS_PAYLOAD path=" << cond.cemcZeroSuppressionPayload
              << " zero_means_unity_channels=" << zeroSentinels << std::endl;
  }

  if (cond.applyJetEnergyScale)
  {
    // Resolve the same key the legacy method consumes, and verify those bytes.
    const std::string payload = CDBInterface::instance()->getUrl(cond.jetEnergyScaleCdbKey);
    if (payload.empty())
    {
      fail("the conditions database returned no jet energy-scale payload for key " + cond.jetEnergyScaleCdbKey);
    }
    verifySha256(payload, cond.jetEnergyScalePayloadSha256);
    {
      TDirectory::TContext context;
      std::unique_ptr<TFile> file(TFile::Open(payload.c_str(), "READ"));
      if (!file || file->IsZombie() || !file->IsOpen())
      {
        fail("resolved jet energy-scale payload cannot be opened: " + payload);
      }
    }
    cond.jetEnergyScalePayload = payload;
    plan.tree.jets.energyScalePayloadPath = payload;
  }

  if (cond.applyVertexReweight && !cond.vertexReweightFileSha256.empty())
  {
    verifySha256(cond.vertexReweightFile, cond.vertexReweightFileSha256);
  }
}

}  // namespace production
}  // namespace photonjet
