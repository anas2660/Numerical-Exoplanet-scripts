#ifndef EXOEQNEW_H
#define EXOEQNEW_H
#include "data.h"
#include "definitions.h"
#include "integrate.h"
/*
  d: distance, inm: initial mass
  inP: initial pressure, a: alpha, b: beta
 */
class ExoEq {
 public:  // Initial conditions
  ExoEq(ℝ ind, ℝ inm, ℝ inP, ℝ a, ℝ b, ℝ M, ℝ R)
      : x({inm, inP}), dxdt(2), t(ind), a(a), b(b), M(M), R(R) {}
  ~ExoEq() {}
  valarray<ℝ> x, dxdt;  // 0 : m, 1 : P
  ℝ t, a, b, M, R, r;   // t : d
  valarray<ℝ> &operator()(ℝ d, valarray<ℝ> x) {
    r = R - d;
    if (x[0] > (a + b) * M)
      return dxdt = {-4.0 * PI * r * r * C->GetDensity(x[1]),
                     G * x[0] * C->GetDensity(x[1]) / r / r};
    else if (x[0] > a * M)
      return dxdt = {-4.0 * PI * r * r * B->GetDensity(x[1]),
                     G * x[0] * B->GetDensity(x[1]) / r / r};
    else if (x[0] > 0)
      return dxdt = {-4.0 * PI * r * r * A->GetDensity(x[1]),
                     G * x[0] * A->GetDensity(x[1]) / r / r};
    return dxdt = {-4.0 * PI * r * r * B->GetDensity(x[1]), 0};
  }
  static void LoadData() {  // statically load data
    A = new DATA("iron.txt");
    B = new DATA("silicate.txt");
    C = new DATA("water.txt");
  }
  static DATA *A, *B, *C;
};
DATA *ExoEq::A, *ExoEq::B, *ExoEq::C;
ℝ ode4(ℝ di, ℝ df, ℝ inm, ℝ inP, ℝ a, ℝ b, ℝ M, ℝ R) {
  const ℝ dd = (df - di) / 2500.0;
  ExoEq eq(di, inm, inP, a, b, M, R);
  while (eq.t <= df) ODE::RK4Step(eq, dd);
  return eq.x[0];  // return m
}
#endif
