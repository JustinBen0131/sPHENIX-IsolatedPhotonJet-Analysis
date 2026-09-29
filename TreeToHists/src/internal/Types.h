// Shared job representation and the small interfaces between lifecycle stages.
#pragma once
#include <nlohmann/json.hpp>
namespace photonjet {
using HistogramJob = nlohmann::json;
void configureHistogramJob(const HistogramJob&);
void fillHistogramJob(const HistogramJob&);
void packageHistogramJob(const HistogramJob&);
void fillExtensions(const HistogramJob&);
}
