#define PLOT true
#define ERROR errorf(ξ + dξ, θ)
#define ½ (1.0/2.0)*
#define N 5
#define DKSI 0.4
#include <cmath>
#include <fstream>
#include <assert.h>
using namespace std;
typedef double real;
typedef double ℝ;

void WritePoint(ofstream& outfile, ℝ x, ℝ y, ℝ dy);

ofstream  eulercsv;
ofstream  rungekuttacsv;

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
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = 0, error = 0;
  WritePoint(eulercsv, ξ, θ, ERROR);
  error += abs(ERROR);
  for (ξ = dξ; ξ < ξmax; ξ += dξ) {
    v += dv(ξ, θ) * dξ;
    θ += dθ(ξ, v) * dξ;
    WritePoint(eulercsv, ξ, θ, ERROR);
    error += abs(ERROR);
  }
  printf("Euler error: %f\n", error);
}


void RungeKuttaFourthOrderIntegration() {
  ℝ θ = 1.0, ξ = 0, dξ = DKSI, ξmax = 10.0, v = 0, error = 0;
  ℝ k1, k2, k3, k4, kv1, kv2, kv3, kv4;
  WritePoint(rungekuttacsv, ξ, θ, ERROR);
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

    WritePoint(rungekuttacsv, ξ, θ, ERROR);
    error += abs(ERROR);
  }
  printf("Runge-Kutta error: %f\n", error);
}





int main(){
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


void WritePoint(ofstream& outfile, ℝ x, ℝ y, ℝ dy){
  outfile << x << "," << y << "," << dy <<  endl;
}
