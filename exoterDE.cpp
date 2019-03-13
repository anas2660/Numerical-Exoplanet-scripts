#include "csvwriter.h"
#include "data.h"
#include "definitions.h"
#include "exoeqnew.h"
#include "polyfit.h"
inline ℝ linear(ℝ x1, ℝ y1, ℝ x2, ℝ y2, ℝ x) {
  return ((y2 - y1) / (x2 - x1)) * (x - x1) + y1;
}
// Linear instead of pchip, extrapolates by default
ℝ interpolate(const vector<ℝ> &x, const vector<ℝ> &y, ℝ point) {
  assert(x.size() == y.size());  // Extrapolate linearly
  if (x.front() > point) return linear(x[0], y[0], x[1], y[1], point);
  if (x.back() < point)
    return linear(x.end()[-2], y.end()[-2], x.back(), y.back(), point);
  array<ℝ, 2> last = {x.front(), y.front()};  // Interpolate
  for (size_t i = 1; i < x.size(); i++) {
    if (x[i] > point) return linear(last[0], last[1], x[i], y[i], point);
    last = {x[i], y[i]};
  }
  throw;
}
void CalcCurve(const string &title, ℝ M, ℝ R, int v) {
  vector<ℝ> Fe = vector<ℝ>(), Si = vector<ℝ>();
  ℝ P0 = 0, a = 0, b, l, i = 0;
  bool u1, u2, j = false, k = false, sign = false;
  while ((a <= 1.0) && !(j && k)) {
    l = (b = (1.0 - a) / 2.0) / 2.0;
    u1 = u2 = false;
    while ((l > 5e-4)) {
      sign = ode4(0, R, M, P0, a, b, M, R) / M > 0;
      b += (sign ? l : -l);
      (sign ? u1 : u2) = true;
      l /= 2.0;
    }
    if (u1 && u2) {
      i++;
      Fe.push_back(a);
      Si.push_back(b);
      k = !(j = true);
    } else
      k = true;
    a += ASTEP;
  }
  printf("i : %f\n", i);
  CSVWriter csv("out/" + title + to_string(v) + ".csv", 2);
  if (i > 0.5)
    for (int n = 0; n < Fe.size(); n++)
      csv.WritePoint(
          {Fe[n], i > 3 ? polyval(polyfit(Fe, Si, 3), {Fe[n]})[0] : Si[n]});
  else
    printf("Curve for (m:%f, r:%f) does't exist.\n", M / Mearth, R / Rearth);
}
void ExoterDE(const string &title, ℝ M0, ℝ sigmaMun, ℝ R0, ℝ sigmaRun) {
  printf("Calculating curves for planet: %s\n", title.c_str());
  const ℝ ratio = interpolate({0.5, 1, 2, 4, 8, 16, 20},
                              {3.3, 3.4, 3.6, 3.9, 4.1, 4.4, 4.5}, M0),
          correction =
              sqrt(pow(sigmaMun / M0, 2) + pow(ratio * sigmaRun / R0, 2)) /
              (abs(sigmaMun / M0) + ratio * abs(sigmaRun / R0)),
          sigmaM = correction * sigmaMun, sigmaR = correction * sigmaRun;
  thread threads[7];
  for (int v = 0; v < 7; v++)
    threads[v] =
        thread(CalcCurve, title, M0 * (1.0 + (v - 3) * sigmaM / M0) * Mearth,
               R0 * (1.0 - (v - 3) * sigmaR / R0) * Rearth, v);
  for (int i = 0; i < 7; i++) threads[i].join();
}
