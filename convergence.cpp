#include "definitions.h"
#include <fstream>
#include <future>
#include <iostream>
#include <numeric>

inline ℝ dv(ℝ ξ, ℝ θ) { return -((ξ * ξ) * pow(θ, N)); }

inline ℝ dθ(ℝ ξ, ℝ v) { return v * (1.0 / (ξ * ξ)); }

// Analytical solutions to calculate error
#if N == 1
inline ℝ θf1(ℝ ξ) { return sin(ξ) / ξ; }
inline ℝ errorf(ℝ ξ, ℝ θ) { return θf1(ξ) - θ; }
#elif N == 5
inline ℝ θf5(ℝ ξ) { return 1.0 / sqrt(1.0 + (ξ * ξ / 3.0)); }
inline ℝ errorf(ℝ ξ, ℝ θ) { return θf5(ξ) - θ; }
#else
inline ℝ errorf(ℝ ξ, ℝ θ) { return 0; }
#endif

real EulerIntegration(unsigned int steps) {
  ℝ θ = 1.0, ξ = 0, ξmax = 10.0, v = 0, error = 0;
  real dξ = ξmax / (real)steps;
  error += abs(ERROR);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    v += dv(ξ, θ) * dξ;
    θ += dθ(ξ, v) * dξ;
    error += abs(ERROR);
  }
  return error / (real)steps;
}

real RungeKuttaIntegration(unsigned int steps) {
  ℝ θ = 1.0, ξ = 0, ξmax = 10.0, v = 0, error = 0, kθ[4], kv[4];
  real dξ = ξmax / (real)steps;
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
  return error / (real)steps;
}

int main() {
  int a;
  cout << "Antal punkter: ";
  cin >> a;
  const unsigned int points = a;
  vector<unsigned int> steps(points);
  std::iota(steps.begin(), steps.end(), 1);
  future<real> eulererrors[points];
  future<real> rkerrors[points];
  for (unsigned int i = 0; i < points; i++) {
    eulererrors[i] = async(launch::async, EulerIntegration, steps[i]);
    rkerrors[i] = async(launch::async, RungeKuttaIntegration, steps[i]);
  }
  ofstream file("error.csv");
  for (unsigned int i = 0; i < points; i++) {
    file << steps[i] << "," << eulererrors[i].get() << "," << rkerrors[i].get()
         << endl;
  }
  file.close();

  if (PLOT)
    system("python errplot.py");
}
