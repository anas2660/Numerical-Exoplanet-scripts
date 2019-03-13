#ifndef LANEEMDEN_H
#define LANEEMDEN_H
#include "definitions.h"
class LaneEmden {  // Lane-Emden system
 public:           // Initial conditions
  LaneEmden(ℝ initialξ) : t(initialξ), x({0, 1}), dxdt(2) {}
  ~LaneEmden() {}
  ℝ t;                  // ξ
  valarray<ℝ> x, dxdt;  // 0 : v, 1 : θ
  valarray<ℝ> &operator()(ℝ ξ, valarray<ℝ> x) {
    dxdt[0] = -((ξ * ξ) * pow(x[1], N));
    dxdt[1] = x[0] * (1.0 / (ξ * ξ));
    return dxdt;
  }
  // Analytical solutions to calculate error
#if N == 1
  ℝ GetError() { return (sin(t) / t) - x[1]; }
#elif N == 5
  ℝ GetError() { return (1.0 / sqrt(1.0 + (t * t / 3.0))) - x[1]; }
#else
  ℝ GetError() { return 0; }
#endif
};
#endif
