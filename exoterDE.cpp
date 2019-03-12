#include "csvwriter.h"
#include "data.h"
#include "definitions.h"
#include "exoeqnew.h"
#include "polyfit.h"

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
    out[i] = {c1[i], c2[i]};
  }
  return out;
}

void ExoterDE(const string &title, real M0, real sigmaMun, real R0,
              real sigmaRun) {
  printf("Calculating curves for planet: %s\n", title.c_str());

  // Initialize vectors
  vector<vector<PlotPoint>> hd(7);

  vector<real> iron[7];
  vector<real> silicate[7];
  vector<real> silicatenew[7];
  for (int i = 0; i < 7; i++) {
    iron[i] = vector<real>();
    silicate[i] = vector<real>();
    silicatenew[i] = vector<real>();
  }

  DensityPressureData A("iron.txt"), B("silicate.txt"), C("water.txt");
  const real step = 0.005;
  const real ratio = interpolate({0.5, 1, 2, 4, 8, 16, 20},
                                 {3.3, 3.4, 3.6, 3.9, 4.1, 4.4, 4.5}, M0);
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
    real P0 = 0, a = 0.0, b, mid, l, sign, lastm, lastsign;
    // real w[7]{0, 0, 0, 0, 0, 0, 0};

    while (((a >= 0.0) && (a <= 1.0)) && !(j && k)) {
      mid = (1.0 - a) / 2.0;
      b = mid;
      l = mid / 2;
      //      printf("l : %f\n", l);
      u1 = u2 = false;
      while ((l > 5e-4)) {
        // Solve the ODEs : // 0.9999 *
        odeout out = ode45(0, R, M, P0, a, b, M, R, A, B, C);
        sign = out.m / M;
        lastm = out.m;
        lastsign = sign;

        // Improve to Newton's method
        if (sign > 0) {
          b += l;
          u1 = true;
        } else {
          b -= l;
          u2 = true;
        }
        l /= 2.0;
        // out.clear(); // no need since it leaves scope anyways
      }

      if (u1 && u2) {
        i++;
        iron[v + 3].push_back(a);
        silicate[v + 3].push_back(b);
        j = true;
        k = false;
      } else {
        k = true;
      }
      a += step;
    }
    printf("i : %d\tm : %f\t ratio : %f\n", i, lastm, lastsign);

    if (i > 3) {
      vector<real> p1 = polyfit(iron[v + 3], silicate[v + 3], 3);
      silicatenew[v + 3] = polyval(p1, iron[v + 3]);
      hd[v + 3] = ternplot(iron[v + 3], silicatenew[v + 3]);
    } else if (i > 0) {
      hd[v + 3] = ternplot(iron[v + 3], silicate[v + 3]);
    } else {
      printf("Curve for mass=%f and radius=%f does not exist.\n",
             M0 * (1 + v * sigmaM), R0 * (1 - v * sigmaR));
    }
  }
  // Output plot
  for (int i = 0; i < hd.size(); i++) {
    string filepath = "out/" + title + to_string(i) + ".csv";
    ofstream file(filepath);
    for (auto it = hd[i].begin(); it != hd[i].end(); ++it) {
      file << it->x << "," << it->y << endl;
    }
    file.close();
  }
}
