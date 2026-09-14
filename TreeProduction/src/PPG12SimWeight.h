#ifndef PPG12_SIM_WEIGHT_H
#define PPG12_SIM_WEIGHT_H

#include <array>
#include <cmath>

namespace PPG12SimWeight
{
  struct Factors
  {
    double slice = 1.0;
    double vertex = 1.0;
    double mix = 1.0;
    double period = 1.0;
    double final = 1.0;
  };

  inline bool Build(const double slice,
                    const double vertex,
                    const double mix,
                    const double period,
                    Factors& factors)
  {
    const std::array<double, 4> inputs{{slice, vertex, mix, period}};
    for (const double value : inputs)
    {
      if (!std::isfinite(value) || value < 0.0) return false;
    }

    factors.slice = slice;
    factors.vertex = vertex;
    factors.mix = mix;
    factors.period = period;
    factors.final = slice * vertex * mix * period;
    return std::isfinite(factors.final) && factors.final >= 0.0;
  }
}

#endif
