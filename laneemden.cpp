#include "laneemden.h"
#include "csvwriter.h"
#include "definitions.h"
#include "integrate.h"
void Integrate(const string &title, ODE::StepType type) {
  ℝ ξmax = 10.0, dξ = DKSI, error = 0;
  CSVWriter csv(title + "int.csv", 3);
  LaneEmden le(dξ);
  csv.WritePoint({0.0, 1.0, 0.0}, ξmax, 200.0);
  while (le.t < ξmax) {
    error += abs(ODE::Step(type, le, dξ).GetError());
    csv.WritePoint({le.t, le.x[1], le.GetError()}, ξmax, 200.0);
  }
  printf("%s error: %f\n", title.c_str(), dξ * error / ξmax);
}
int main() {
  vector<thread *> t = {new thread(Integrate, "euler", ODE::Euler),
                        new thread(Integrate, "rungekutta", ODE::RK4)};
  for (auto thr : t) thr->join();
  if (PLOT) system("python intplot.py");
}
