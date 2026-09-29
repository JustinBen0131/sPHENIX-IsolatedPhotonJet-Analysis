// Validate storage axes and activate strict checks before any ROOT booking.
#include "Types.h"
#include "../SelectionInput.h"
namespace photonjet {
void configureHistogramJob(const HistogramJob& job) {
  configuredXjEdges = job.at("histograms").at("native_xj_edges").get<std::vector<double>>();
  if (configuredXjEdges.size() < 2) throw std::runtime_error("empty xJ axis");
  for (size_t i=0; i<configuredXjEdges.size(); ++i)
    if (!std::isfinite(configuredXjEdges[i]) || (i && configuredXjEdges[i]<=configuredXjEdges[i-1]))
      throw std::runtime_error("xJ axis must be finite and increasing");
  if (configuredXjEdges.front()!=0 || configuredXjEdges.back()!=3)
    throw std::runtime_error("xJ storage must retain support [0,3]");
  strictReduction = true;
  progressPath = job.value("progress", "");
  progressTotal = job.value("entries", 0LL);
  configureVertexFit(job.value("vertex_reweighting",HistogramJob::object()));
}
}
