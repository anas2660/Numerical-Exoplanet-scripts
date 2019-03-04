#include "csvwriter.h"
#include "definitions.h"

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

void EulerIntegration() {
  CSVWriter csv("eulerint.csv");
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = 0, error = 0, m = 0, P = 0.0;
  csv.WritePointSimplify(log10(P), θ, 0, ξmax, 500);
  error += abs(ERROR);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    v += dv(ξ, θ) * dξ;
    θ += dθ(ξ, v) * dξ;        // density
    m += 4 * PI * (ξ * ξ) * θ; // mass
    P += -G * m * θ / (ξ * ξ); // pressure
    csv.WritePointSimplify(log10(P), θ, 0, ξmax, 500);
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
    csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
    error += abs(ERROR);
  }
  printf("Runge-Kutta error: %f\n", error);
}

int main() {
  EulerIntegration();
  if (PLOT)
    system("python intplot.py");
}
