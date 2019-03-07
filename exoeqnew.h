#ifndef EXOEQNEW_H
#define EXOEQNEW_H
#include "data.h"
#include "definitions.h"

struct ydot {
  real m, P;
};

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

#endif
