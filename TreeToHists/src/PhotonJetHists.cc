// Run one bounded source reduction or one checked sample consolidation.
#include "PhotonJetHists.h"
#include "internal/Types.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
namespace photonjet {
int runHistogramJob(const std::string& path) {
  try {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open histogram job: " + path);
  HistogramJob job; input >> job;
  if (job.at("operation") == "package") packageHistogramJob(job);
  else { configureHistogramJob(job); fillHistogramJob(job); }
  return 0;
  } catch(const std::exception& error) {
    std::cerr << "PhotonJetHists: " << error.what() << std::endl;
    return 1;
  }
}
}
