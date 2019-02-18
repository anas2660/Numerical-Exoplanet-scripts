#define PLOT true
#define ½ (1.0 / 2.0) *
#include <assert.h>
#include <cmath>
#include <fstream>
using namespace std;
typedef double real;
typedef double ℝ;

void WritePoint(ofstream &outfile, ℝ x, ℝ y, ℝ dy);

ℝ n = 1;

ofstream eulercsv;
ofstream rungekuttacsv;

// Analytical solution for n = 5
ℝ θf5(ℝ ξ) { return 1.0 / sqrt(1.0 + (ξ * ξ / 3.0)); }

// Analytical solution for n = 1
ℝ θf1(ℝ ξ) { return sin(ξ) / ξ; }

ℝ dv(ℝ ξ, ℝ θ, ℝ v) { return -((ξ * ξ) * pow(θ, n)) - 2 * ξ * v; }

ℝ dθ(ℝ v) { return v; }

void EulerIntegration() {
  ℝ θ = 1.0, ξ = 0, dξ = 0.05, ξmax = 10.0, v = -1.0, error = 0;
  WritePoint(eulercsv, ξ, θ, θf1(ξ + dξ) - θ);
  error += abs(θf1(ξ + dξ) - θ);
  for (ξ = 5 * dξ; ξ < ξmax; ξ += dξ) {
    v += dv(ξ, θ, v) * dξ;
    θ += dθ(v) * dξ;
    WritePoint(eulercsv, ξ, θ, θf1(ξ) - θ);
    error += abs(θf1(ξ + dξ) - θ);
  }
  printf("Euler error: %f\n", error);
}

void RungeKuttaFourthOrderIntegration() {
  ℝ θ = 1.0, ξ = 0, dξ = 0.05, ξmax = 10.0, v = 0, k1, k2, k3, k4, error = 0, k11, k12,k13,k14;
  

  WritePoint(rungekuttacsv, ξ, θ, θf1(ξ + dξ) - θ);
  error += abs(θf1(ξ + dξ) - θ);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    // Runge-Kutta step for v
    /*k1 = dξ * dv(ξ, θ, v);
    k11 = dξ * dθ(v);
    k2 = dξ * dv(ξ + ½ dξ, θ + ½ k11, v + ½ k1);
    k12 =  dξ * dθ(v + ½ k1);
    k3 = dξ * dv(ξ + ½ dξ, θ + ½ k12, v + ½ k2);
    k13 = dξ * dθ(v + ½ k2);
    k4 = dξ * dv(ξ + dξ, θ + k13, v + k3);
    k14 = dξ * dθ(v + k3);
    v += (k2 + k3 + ½(k1 + k4)) / 3.0;
    θ += (k12 + k13 + ½(k11 + k14)) / 3.0;*/
    k1 = dξ * dv(ξ, θ, v);
    k2 = dξ * dv(ξ + ½ dξ, θ + ½ k1, v);
    k3 = dξ * dv(ξ + ½ dξ, θ + ½ k2, v);
    k4 = dξ * dv(ξ + dξ, θ + k3, v);
    v += (k2 + k3 + ½(k1 + k4)) / 3.0;

    /*
     */
    // Runge-Kutta step for θ
    /*
      k1 = dξ * dθ(ξ, v);
      k2 = dξ * dθ(ξ + ½ dξ, v + ½ k1);
      k3 = dξ * dθ(ξ + ½ dξ, v + ½ k2);
      k4 = dξ * dθ(ξ + dξ, v + k3);
      θ += (k2 + k3 + ½(k1 + k4)) / 3.0;*/
    θ += dθ(v) * dξ;
    WritePoint(rungekuttacsv, ξ, θ, θf1(ξ) - θ);
    error += abs(θf1(ξ + dξ) - θ);
  }
  printf("Runge-Kutta error: %f\n", error);
}

int main() {
  eulercsv = ofstream("eulerint.csv");
  rungekuttacsv = ofstream("rungekuttaint.csv");
  assert(eulercsv.is_open());
  assert(rungekuttacsv.is_open());
  EulerIntegration();
  RungeKuttaFourthOrderIntegration();
  eulercsv.close();
  rungekuttacsv.close();
  if (PLOT) {
    system("python intplot.py");
  }
}

void WritePoint(ofstream &outfile, ℝ x, ℝ y, ℝ dy) {
  outfile << x << "," << y << "," << dy << endl;
}
