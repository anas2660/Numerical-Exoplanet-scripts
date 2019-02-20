#include "definitions.h"
#include "csvwriter.h"
#include <thread>


// Analytical solution for n = 5
ℝ θf5(ℝ ξ){
  return 1.0/sqrt(1.0+(ξ*ξ/3.0));
}

// Analytical solution for n = 1
ℝ θf1(ℝ ξ){
  return sin(ξ)/ξ;
}


ℝ dv(ℝ ξ, ℝ θ){
  return -((ξ*ξ)*pow(θ, N));
}

ℝ dθ(ℝ ξ, ℝ v) {
  return v*(1.0/(ξ*ξ));
}

#if N == 1
ℝ errorf(ℝ ξ, ℝ θ) { return θf1(ξ) - θ; }
#elif N == 5
ℝ errorf(ℝ ξ, ℝ θ) { return θf5(ξ) - θ; }
#else
ℝ errorf(ℝ ξ, ℝ θ) { return 0; }
#endif

void EulerIntegration(){
  CSVWriter  csv("eulerint.csv");
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = 0, error = 0;
  csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
  error += abs(ERROR);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    v += dv(ξ, θ) * dξ;
    θ += dθ(ξ, v) * dξ;
    csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
    error += abs(ERROR);
  }
  printf("Euler error: %f\n", error);
}


void RungeKuttaFourthOrderIntegration() {
  CSVWriter  csv("rungekuttaint.csv");
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = 0, error = 0;
  ℝ k1, k2, k3, k4, kv1, kv2, kv3, kv4;
  csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
  error += abs(ERROR);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    // Runge-Kutta step for v and theta
    kv1 = dξ * dv(ξ, θ);
    k1 = dξ * dθ(ξ, v);
    kv2 = dξ * dv(ξ + ½ dξ, θ + ½ k1);
    k2 = dξ * dθ(ξ + ½ dξ, v + ½ kv1);
    kv3 = dξ * dv(ξ + ½ dξ, θ + ½ k2);
    k3 = dξ * dθ(ξ + ½ dξ, v+ ½ kv2);
    kv4 = dξ * dv(ξ + dξ, θ + k3);
    k4 = dξ * dθ(ξ + dξ, v + kv3);
    v += (kv2 + kv3 + ½ (kv1 + kv4)) / 3.0;
    θ += (k2 + k3 + ½ (k1 + k4)) / 3.0;

    csv.WritePointSimplify(ξ, θ, ERROR, ξmax, 200);
    error += abs(ERROR);
  }
  printf("Runge-Kutta error: %f\n", error);
}



int main(){
  thread euler(EulerIntegration);
  thread rungekutta(RungeKuttaFourthOrderIntegration);
  euler.join();
  rungekutta.join();
  if (PLOT) system("python intplot.py");
}
