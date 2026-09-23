// ProductionInputs.cc
//
// Which physical source files enter Fun4All, and how are they synchronized?
// Owns input-list interpretation and input-manager registration.
// Does not own detector reconstruction.
// Called through Production.h before PhotonJetTree captures the event.
//
#include "Production.h"

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace photonjet
{
namespace production
{
using detail::fail;
using detail::trim;

namespace
{
//
// Input-list parsing
//
// One line per source. Columns are whitespace separated and "NONE" or "-"
// marks an absent lane. This is the only code that knows the column layout.
//
//   paired data                     CALO  AUXILIARY
//   reconstructed simulation        CALO  G4  TRUTH_JETS  GLOBAL  MBD
//   archived double interaction     G4  TRUTH_JETS
//
std::vector<std::string> splitColumns(const std::string& line)
{
  std::istringstream stream(line);
  std::vector<std::string> columns;
  std::string token;
  while (stream >> token)
  {
    if (token == "NONE" || token == "none" || token == "-")
    {
      columns.emplace_back();
    }
    else
    {
      columns.push_back(token);
    }
  }
  return columns;
}

void appendColumn(std::vector<std::string>& destination,
                  const std::vector<std::string>& columns,
                  const std::size_t index)
{
  if (index < columns.size() && !columns[index].empty())
  {
    destination.push_back(columns[index]);
  }
}

// Preserve the resolved stream policy: ordinary DST streams synchronize;
// embedded and archived-DI streams use the existing NoSync path.
Fun4AllInputManager* makeInputManager(const InputStream& stream)
{
  const std::string name = "DST_" + stream.name + "_IN";
  if (stream.synchronized)
  {
    return new Fun4AllDstInputManager(name);
  }
  return new Fun4AllNoSyncDstInputManager(name);
}
}  // namespace

namespace detail
{
std::vector<InputStream> loadInputs(const Job& job, const Reconstruction& reconstruction)
{
  std::ifstream input(job.inputList);
  if (!input)
  {
    fail("cannot read input list: " + job.inputList);
  }

  std::vector<std::string> primary, auxiliary, g4, truthJets, global, mbd;

  std::string line;
  std::size_t sourceRows = 0;
  while (std::getline(input, line))
  {
    line = trim(line);
    if (line.empty() || line.front() == '#')
    {
      continue;
    }
    const std::vector<std::string> columns = splitColumns(line);
    if (columns.empty())
    {
      continue;
    }

    if (++sourceRows > 1) fail("one physical source bundle per invocation is required");
    const std::size_t expected = reconstruction.inputLayout == InputLayout::ReconstructedSimulation ? 5 : 2;
    if (columns.size() != expected) fail("input-list row has an unexpected column count");
    switch (reconstruction.inputLayout)
    {
      case InputLayout::PairedData:
        appendColumn(primary, columns, 0);
        appendColumn(auxiliary, columns, 1);
        break;

      case InputLayout::ReconstructedSimulation:
        appendColumn(primary, columns, 0);
        appendColumn(g4, columns, 1);
        appendColumn(truthJets, columns, 2);
        appendColumn(global, columns, 3);
        appendColumn(mbd, columns, 4);
        break;

      case InputLayout::ArchivedDoubleInteraction:
        appendColumn(g4, columns, 0);
        appendColumn(truthJets, columns, 1);
        break;
    }
  }

  std::vector<InputStream> streams;

  auto add = [&](const char* name, std::vector<std::string> files,
                 const bool requiredStream, const bool synchronized)
  {
    if (files.empty() && !requiredStream)
    {
      return;
    }
    InputStream stream;
    stream.name = name;
    stream.files = std::move(files);
    stream.required = requiredStream;
    stream.synchronized = synchronized;
    streams.push_back(std::move(stream));
  };

  // Embedded and archived inputs are read without Fun4All synchronisation,
  // as the historical production did; paired data is synchronised.
  const bool sync = !(reconstruction.isEmbedded || reconstruction.isArchivedDoubleInteraction);

  switch (reconstruction.inputLayout)
  {
    case InputLayout::PairedData:
      add("PRIMARY", std::move(primary), true, true);
      add("AUXILIARY", std::move(auxiliary), false, true);
      break;

    case InputLayout::ArchivedDoubleInteraction:
      add("G4HITS", std::move(g4), true, false);
      add("TRUTH_JETS", std::move(truthJets), true, false);
      break;

    case InputLayout::ReconstructedSimulation:
      add("CALO", std::move(primary), true, sync);
      add("G4HITS", std::move(g4), true, sync);
      add("TRUTH_JETS", std::move(truthJets), false, sync);
      add("GLOBAL", std::move(global), false, sync);
      add("MBD", std::move(mbd), false, sync);
      break;
  }

  return streams;
}

bool hasInputStream(const Plan& plan, const std::string& name)
{
  return std::any_of(plan.inputs.begin(), plan.inputs.end(),
                     [&](const InputStream& stream) {
                       return stream.name == name && !stream.files.empty();
                     });
}
}  // namespace detail

// RegisterInputs
//
void RegisterInputs(Fun4AllServer* server, const Plan& plan)
{
  if (!server) fail("null Fun4AllServer");

  // Reconstructing the calorimeter needs the tower geometry run node.
  if (plan.reconstruction.reconstructCalorimeter)
  {
    auto* geometry = new Fun4AllRunNodeInputManager("DST_GEO");
    geometry->AddFile(CDBInterface::instance()->getUrl("calo_geo"));
    server->registerInputManager(geometry);
  }

  for (const InputStream& stream : plan.inputs)
  {
    if (stream.files.empty())
    {
      if (stream.required)
      {
        fail("required input stream has no files: " + stream.name);
      }
      continue;
    }
    Fun4AllInputManager* manager = makeInputManager(stream);
    manager->Verbosity(plan.job.verbosity);
    for (const std::string& file : stream.files)
    {
      manager->AddFile(file);
    }
    server->registerInputManager(manager);
  }
}

}  // namespace production
}  // namespace photonjet
