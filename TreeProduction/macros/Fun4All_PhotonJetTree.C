//
// Fun4All_PhotonJetTree.C
//
// Steering macro for the isolated-photon + jet tree production.
//
//   root -l -b -q 'Fun4All_PhotonJetTree.C("pp_data", "inputs.list",
//                                           "out.root", 47123, 0)'
//
// It loads and resolves the configuration, creates the server, configures
// conditions and reconstruction, registers the producer and the inputs, runs,
// and ends. Physics lives in the producer library; nothing is computed here.
//
#include <TreeProduction/PhotonJetTree.h>
#include <TreeProduction/internal/Production.h>

#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(libPhotonJetTree.so)

int Fun4All_PhotonJetTree(const std::string& profile,
                          const std::string& inputList,
                          const std::string& outputFile,
                          const int run,
                          const int segment,
                          const std::int64_t nEvents = -1,
                          const std::int64_t firstEntry = 0,
                          const std::string& configFile = "config/tree_production.yaml",
                          const std::int64_t sourceFileOrdinal = -1,
                          const std::string& inputFileSha256 = "",
                          const std::string& sourceManifestSha256 = "",
                          const std::string& sample = "",
                          const std::string& period = "",
                          const std::string& siDiRole = "",
                          const int verbosity = 0)
{
  using photonjet::production::Job;
  using photonjet::production::Plan;

  try
  {
    Job job;
    job.profile = profile;
    job.inputList = inputList;
    job.outputFile = outputFile;
    job.run = run;
    job.segment = segment;
    job.firstEntry = firstEntry;
    job.nEvents = nEvents;
    job.sourceFileOrdinal = sourceFileOrdinal;
    job.inputFileSha256 = inputFileSha256;
    job.sourceManifestSha256 = sourceManifestSha256;
    job.sample = sample;
    job.period = period;
    job.siDiRole = siDiRole;
    job.verbosity = verbosity;

    Plan plan = photonjet::production::LoadPlan(configFile, job);
    photonjet::production::ConfigureConditions(plan);
    photonjet::production::Print(plan);

    Fun4AllServer* server = Fun4AllServer::instance();
    server->Verbosity(verbosity);

    auto* producer = new PhotonJetTree(plan.tree, "PhotonJetTree");

    photonjet::production::RegisterReconstruction(server, plan, producer);
    photonjet::production::RegisterInputs(server, plan);

    if (firstEntry > 0)
    {
      if (server->skip(static_cast<int>(firstEntry)) != 0)
      {
        producer->markAborted();
        server->End();
        delete server;
        gSystem->Exit(1);
        return 1;
      }
    }

    const int status = nEvents < 0 ? server->run() : server->run(static_cast<int>(nEvents));
    if (status != 0)
    {
      std::cerr << "Fun4All_PhotonJetTree: the event loop stopped with status "
                << status << "; the output is marked aborted" << std::endl;
      producer->markAborted();
    }

    const int endStatus = server->End();
    delete server;

    if (status != 0 || endStatus != 0)
    {
      gSystem->Exit(1);
      return 1;
    }

    std::cout << "Fun4All_PhotonJetTree: complete" << std::endl;
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "Fun4All_PhotonJetTree failed: " << error.what() << std::endl;
    gSystem->Exit(1);
    return 1;
  }
}
