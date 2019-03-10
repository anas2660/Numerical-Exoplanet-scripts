#include "csvwriter.h"
#include "data.h"
#include "definitions.h"
#include "exoeqnew.h"
#include "polyfit.h"

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
  CSVWriter csv("euler.csv");
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

real linear(real x1, real y1, real x2, real y2, real x) {
  const real a = (y2 - y1) / (x2 - x1);
  const real b = y1 - a * x1;
  return (a * x + b);
}

// Linear instead of pchip, extrapolates by default
real interpolate(const vector<real> &x, const vector<real> &y, real point) {
  if (x.size() != y.size()) {
    printf("size of x does not match size of y\n");
    exit(-1);
  }

  // Interpolate linearly before and after x vector
  if (x.front() > point) {
    return linear(x[0], y[0], x[1], y[1], point);
  } else if (x.back() < point) {
    return linear(x.end()[-2], y.end()[-2], x.back(), y.back(), point);
  }
  real lastx = x.front(), lasty = y.front();
  for (size_t i = 1; i < x.size(); i++) {
    if (x[i] > point) {
      return linear(lastx, lasty, x[i], y[i], point);
    }
    lastx = x[i];
    lasty = y[i];
  }
  return 1;
}

struct PlotPoint {
  real x, y;
};

vector<PlotPoint> ternplot(const vector<real> &c1, const vector<real> &c2) {
  vector<PlotPoint> out(c1.size());
  for (int i = 0; i < c1.size(); i++) {
    out[i] = {0.5 * (1.0 + c1[i] - c2[i]),
              (sqrt(3.0) / 2.0) * (1.0 - c1[i] - c2[i])};
  }
  return out;
}

// outer index : v, inner index : i
vector<vector<real>> iron(7);
vector<vector<real>> silicate(7);
vector<vector<real>> silicatenew(7);
vector<vector<PlotPoint>> hd(7);

void ExoterDE(real M0, real sigmaMun, real R0, real sigmaRun) {
  // Initialize vectors
  for (int i = 0; i < iron.size(); i++) {
    iron[i] = vector<real>(5);
    silicate[i] = vector<real>(5);
    // silicatenew[i] = vector<real>(5);
  }
  DensityPressureData A("iron.txt"), B("silicate.txt"), C("water.txt");
  const real step = 0.005;

  vector<real> x = {0.5, 1, 2, 4, 8, 16, 20},
               y = {3.3, 3.4, 3.6, 3.9, 4.1, 4.4, 4.5};
  const real ratio = interpolate(x, y, M0);
  const real correction =
      sqrt(pow((sigmaMun / M0), 2) + pow(ratio * (sigmaRun / R0), 2)) /
      (abs(sigmaMun / M0) + ratio * abs(sigmaRun / R0));
  const real sigmaM = correction * sigmaMun;
  const real sigmaR = correction * sigmaRun;

  for (int v = -3; v < 4; v++) {

    real M = M0 * (1.0 + v * sigmaM / M0) * Mearth;
    real R = R0 * (1.0 - v * sigmaR / R0) * Rearth;

    int i = 0;
    bool u1 = false, u2 = false, j = false, k = false;
    real m0 = M0, P0 = 0, a = 0.0, b, mid, l, sign;
    // real w[7]{0, 0, 0, 0, 0, 0, 0};
    while (((a >= 0.0) && (a <= 1.0)) && !(j && k)) {

      mid = (1.0 - a) / 2.0;
      b = mid;
      l = mid / 2;
      //      printf("l : %f\n", l);
      u1 = u2 = false;
      while (l > 5e-4) {
        // Solve the ODEs :
        vector<odeout> out = ode45(0, 0.9999 * R, m0, P0, a, b, m0, R, A, B, C);
        sign = out.back().m / M;

        // Improve to Newton's method
        if (sign > 0) {
          b += l;
          u1 = true;
        } else {
          b -= l;
          u2 = true;
        }
        l /= 2.0;
        out.clear(); // no need since it leaves scope anyways
      }
      if (u1 && u2) {
        i++;
        iron[v + 3][i - 1] = a;
        silicate[v + 3][i - 1] = b;
        j = true;
        k = false;
      } else {
        k = true;
      }
      a += step;
    }
    printf("i : %d\n", i);

    if (i > 3) {
      vector<real> ta(iron[v + 3].begin(), iron[v + 3].begin() + i),
          tb(silicate[v + 3].begin(), silicate[v + 3].begin() + i);
      vector<real> p1 = polyfit(ta, tb, 3);
      silicatenew[v + 3] = polyval(p1, ta);
      hd[v + 3] = ternplot(ta, silicatenew[v + 3]);
    } else if (i > 0) {
      vector<real> ta(iron[v + 3].begin(), iron[v + 3].begin() + i),
          tb(silicate[v + 3].begin(), silicate[v + 3].begin() + i);
      hd[v + 3] = ternplot(ta, tb);
    } else {
      hd[v + 3] = vector<PlotPoint>();
      printf("Curve for mass=%f and radius=%f does not exist.\n",
             M0 * (1 + v * sigmaM), R0 * (1 - v * sigmaR));
    }
  }
}

int main() {
  /*  EulerIntegration();
  if (PLOT)
    system("python intplot.py");*/

  // Planet 1
  // ExoterDE(4.8, 0.8, 1.68, 0.09);
  ExoterDE(6.55, 0.98, 2.678, 0.13);

  // Output plot
  for (int i = 0; i < hd.size(); i++) {
    string filepath = "out/plot" + to_string(i) + ".csv";
    ofstream file(filepath);
    for (auto it = hd[i].begin(); it != hd[i].end(); ++it) {
      file << it->x << "," << it->y << endl;
    }
    file.close();
  }
}
