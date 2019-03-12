#include "csvwriter.h"
#include "definitions.h"

inline ℝ dv(ℝ ξ, ℝ θ, ℝ v) { return -(pow(θ, N)) - 2 * v / ξ; }
inline ℝ dθ(ℝ v) { return v; }

void EulerIntegration() {
  CSVWriter csv("eulerint.csv", 3);
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = 0, error = 0;
  csv.WritePoint({ξ, θ, ERROR}, ξmax, 200);
  error += abs(ERROR);
  for (ξ = 5 * dξ; ξ < ξmax; ξ += dξ) {
    v += dv(ξ, θ, v) * dξ;
    θ += dθ(v) * dξ;
    csv.WritePoint({ξ, θ, ERROR}, ξmax, 200);
    error += abs(ERROR);
  }
  printf("Euler error: %f\n", dξ * error / ξmax);
}
void RungeKuttaIntegration() {
  CSVWriter csv("rungekuttaint.csv", 3);
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = 0, error = 0, kθ[4], kv[4];
  csv.WritePoint({ξ, θ, ERROR}, ξmax, 200);
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
    csv.WritePoint({ξ, θ, ERROR}, ξmax, 200);
    error += abs(ERROR);
  }
  printf("Runge-Kutta error: %f\n", dξ * error / ξmax);
}

int main() {
  thread euler(EulerIntegration);
  thread rungekutta(RungeKuttaIntegration);
  euler.join();
  rungekutta.join();
  if (PLOT)
    system("python intplot.py");
}
