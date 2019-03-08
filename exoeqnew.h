#ifndef EXOEQNEW_H
#define EXOEQNEW_H
#include "data.h"
#include "definitions.h"

struct ydot {
  real m, P;
};

// d: distance
struct odeout {
  real d, m, P;
};

/*
  d: distance
  m: initial mass
  P: initial pressure
  a: control parameter
  b: control parameter
 */

ydot Exoeqnew(real d, real m, real P, real a, real b, real M, real R, DATA A,
              DATA B, DATA C) {
  if (m > (a + b) * M) {
    return {-4.0 * PI * (R - d) * (R - d) * C.GetDensity(P),
            G * m * C.GetDensity(P) / (R - d) / (R - d)};
  } else if (m > a * M) {
    return {-4.0 * PI * (R - d) * (R - d) * B.GetDensity(P),
            G * m * B.GetDensity(P) / (R - d) / (R - d)};
  } else if (m > 0) {
    return {-4.0 * PI * (R - d) * (R - d) * A.GetDensity(P),
            G * m * A.GetDensity(P) / (R - d) / (R - d)};
  } else {
    return {-4.0 * PI * (R - d) * (R - d) * B.GetDensity(P), 0};
  }
}

/*
  di: distance initial
  df: distance final
 */

vector<odeout> ode45(real di, real df, real inm, real inP, real a, real b,
                     real M, real R, DATA A, DATA B, DATA C) {
  // CSVWriter csv("rungekuttaint.csv");
  real m = inm, P = inP;
  vector<odeout> out = vector<odeout>();
  real d = di, dd = 256, km[4], kP[4];
  out.push_back({d, m, P});
  ydot result;
  for (d += dd; d < df; d += dd) {
    // Runge-Kutta step for v and theta
    result = Exoeqnew(d, m, P, a, b, M, R, A, B, C);
    km[0] = dd * result.m;
    kP[0] = dd * result.P;

    result = Exoeqnew(d + ½ dd, m + ½ km[0], P + ½ kP[0], a, b, M, R, A, B, C);
    km[1] = dd * result.m;
    kP[1] = dd * result.P;

    result = Exoeqnew(d + ½ dd, m + ½ km[1], P + ½ kP[1], a, b, M, R, A, B, C);
    km[2] = dd * result.m;
    kP[2] = dd * result.P;

    result = Exoeqnew(d + dd, m + km[2], P + kP[2], a, b, M, R, A, B, C);
    km[3] = dd * result.m;
    kP[3] = dd * result.P;

    m += (km[1] + km[2] + ½(km[0] + km[3])) / 3.0;
    P += (kP[1] + kP[2] + ½(kP[0] + kP[3])) / 3.0;
    out.push_back({d, m, P});
  }
  return out;
}

// #if USING_BOOST == 0
/*#else
#include <boost/numeric/odeint/

#endif
*/
#endif
