#include "csvwriter.h"
#include "definitions.h"
#include <fstream>
#include <future>

inline ℝ dv(ℝ ξ, ℝ θ) { return -((ξ * ξ) * pow(θ, N)); }
inline ℝ dθ(ℝ ξ, ℝ v) { return v * (1.0 / (ξ * ξ)); }

ℝ EulerIntegration(unsigned int steps) {
  ℝ θ = 1.0, ξ = 0, ξmax = 10.0, v = 0, error = 0;
  ℝ dξ = ξmax / (ℝ)steps;
  error += abs(ERROR);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    v += dv(ξ, θ) * dξ;
    θ += dθ(ξ, v) * dξ;
    error += abs(ERROR);
  }
  return error;
}

ℝ RungeKuttaIntegration(unsigned int steps) {
  ℝ θ = 1.0, ξ = 0, ξmax = 10.0, v = 0, error = 0, kθ[4], kv[4];
  ℝ dξ = ξmax / (ℝ)steps;
  error += abs(ERROR);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    // Runge-Kutta step for v and theta
    kv[0] = dξ * dv(ξ, θ);
    kθ[0] = dξ * dθ(ξ, v);
    kv[1] = dξ * dv(ξ + ½ dξ, θ + ½ kθ[0]);
    kθ[1] = dξ * dθ(ξ + ½ dξ, v + ½ kv[0]);
    kv[2] = dξ * dv(ξ + ½ dξ, θ + ½ kθ[1]);
    kθ[2] = dξ * dθ(ξ + ½ dξ, v + ½ kv[1]);
    kv[3] = dξ * dv(ξ + dξ, θ + kθ[2]);
    kθ[3] = dξ * dθ(ξ + dξ, v + kv[2]);
    v += (kv[1] + kv[2] + ½(kv[0] + kv[3])) / 3.0;
    θ += (kθ[1] + kθ[2] + ½(kθ[0] + kθ[3])) / 3.0;
    error += abs(ERROR);
  }
  return error;
}

int main(int argc, char *argv[]) {
  if (argc > 1) {
    const unsigned int points = stoi(argv[1]);
    future<ℝ> eulererrors[points], rkerrors[points];
    for (unsigned int i = 0; i < points; i++) {
      eulererrors[i] = async(launch::async, EulerIntegration, i + 1);
      rkerrors[i] = async(launch::async, RungeKuttaIntegration, i + 1);
    }
    CSVWriter csv("error.csv", 3); // Write to csv file
    for (unsigned int i = 0; i < points; i++)
      csv.WritePoint({(double)i + 1, eulererrors[i].get(), rkerrors[i].get()});
  }
  if (PLOT)
    system("python errplot.py");
}
