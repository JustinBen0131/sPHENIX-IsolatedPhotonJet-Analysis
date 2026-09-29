// Public reduction lifecycle. A resolved JSON job binds inputs, definitions and
// output paths; the Python command owns planning and immutable publication.
#pragma once
#include <string>
namespace photonjet {
int runHistogramJob(const std::string& jobPath);
}
