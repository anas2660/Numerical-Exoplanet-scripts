#include "csvwriter.h"
#include "definitions.h"

inline ℝ dv(ℝ ξ, ℝ θ, ℝ v) { return -(pow(θ, N)) - 2 * v / ξ; }

inline ℝ dθ(ℝ v) { return v; }

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

void EulerIntegration() {
  CSVWriter csv("eulerint.csv");
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = -1.0, error = 0;
  csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
  error += abs(ERROR);
  for (ξ = 5 * dξ; ξ < ξmax; ξ += dξ) {
    v += dv(ξ, θ, v) * dξ;
    θ += dθ(v) * dξ;
    csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
    error += abs(ERROR);
  }
  printf("Euler error: %f\n", error);
}

void RungeKuttaIntegration() {
  CSVWriter csv("rungekuttaint.csv");
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = 0, error = 0, kθ[4], kv[4];
  csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
  error += abs(ERROR);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    // Runge-Kutta step for v and theta
    kv[0] = dξ * dv(ξ, θ, v);
    kθ[0] = dξ * dθ(v);
    kv[1] = dξ * dv(ξ + ½ dξ, θ + ½ kθ[0], v + ½ kv[0]);
    kθ[1] = dξ * dθ(v + ½ kv[0]);
    kv[2] = dξ * dv(ξ + ½ dξ, θ + ½ kθ[1], v + ½ kv[1]);
    kθ[2] = dξ * dθ(v + ½ kv[1]);
    kv[3] = dξ * dv(ξ + dξ, θ + kθ[2], v + kv[2]);
    kθ[3] = dξ * dθ(v + kv[2]);
    v += (kv[1] + kv[2] + ½(kv[0] + kv[3])) / 3.0;
    θ += (kθ[1] + kθ[2] + ½(kθ[0] + kθ[3])) / 3.0;
    csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
    error += abs(ERROR);
  }
  printf("Runge-Kutta error: %f\n", error);
}

int main() {
  thread euler(EulerIntegration);
  thread rungekutta(RungeKuttaIntegration);
  euler.join();
  rungekutta.join();
  if (PLOT)
    system("python intplot.py");
}
