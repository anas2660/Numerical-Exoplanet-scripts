#include "csvwriter.h"
#include "definitions.h"
#include "integrate.h"
#include "laneemden.h"
ℝ Integrate(ODE::StepType type, unsigned int steps) {
  ℝ error = 0, ξmax = 10.0, dξ = ξmax / (ℝ)steps;
  LaneEmden le(dξ);
  while (le.t < ξmax) error += abs(ODE::Step(type, le, dξ).GetError());
  return error / (ℝ)steps;
}
int main(int argc, char *argv[]) {
  assert(argc > 1);
  const unsigned int points = stoi(argv[1]);
  future<ℝ> eule[points], rk2e[points], rk4e[points];
  CSVWriter *csv = new CSVWriter("error.csv", 4);
  for (unsigned int i = 0; i < points; i++) {
    eule[i] = async(launch::async, Integrate, ODE::Euler, i + 1);
    rk2e[i] = async(launch::async, Integrate, ODE::RK2, i + 1);
    rk4e[i] = async(launch::async, Integrate, ODE::RK4, i + 1);
  }
  for (size_t i = 0; i < points; i++)
    csv->WritePoint({(ℝ)i + 1, eule[i].get(), rk2e[i].get(), rk4e[i].get()});
  delete csv;  // make sure its closed
  if (PLOT) system("python errplot.py");
}
