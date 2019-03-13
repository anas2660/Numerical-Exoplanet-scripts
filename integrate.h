#ifndef INTEGRATE_H
#define INTEGRATE_H
#include "definitions.h"
namespace ODE {
enum StepType { Euler, RK2, RK4 };
template <class System> System &EulerStep(System &s, ℝ dt) {
  s.x += dt * s(s.t, s.x);
  s.t += dt;
  return s;
}
template <class System> System &RK2Step(System &s, ℝ dt) {
  valarray<ℝ> k[2];
  k[0] = dt * s(s.t, s.x);
  k[1] = dt * s(s.t + ½ dt, s.x + ½ k[0]);
  s.x += k[1];
  s.t += dt;
  return s;
}
template <class System> System &RK4Step(System &s, ℝ dt) {
  valarray<ℝ> k[4];
  k[0] = dt * s(s.t, s.x);
  k[1] = dt * s(s.t + ½ dt, s.x + ½ k[0]);
  k[2] = dt * s(s.t + ½ dt, s.x + ½ k[1]);
  k[3] = dt * s(s.t + dt, s.x + k[2]);
  s.x += (k[1] + k[2] + ½(k[0] + k[3])) / 3.0;
  s.t += dt;
  return s;
}
template <class System> inline System &Step(StepType type, System &s, ℝ dt) {
  switch (type) {
    case Euler: return EulerStep(s, dt);
    case RK2: return RK2Step(s, dt);
    case RK4: return RK4Step(s, dt);
  }
}
}  // namespace ODE
#endif
