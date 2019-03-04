#include "data.h"
#include <fstream>

DensityPressureData::DensityPressureData(const char *filename) {
  ifstream file(filename);
  if (!file.is_open())
    throw;
  real d, p;
  while (file >> d >> p) {
    Density.emplace_back(d);
    Pressure.emplace_back(p);
  }
  file.close();
}

DensityPressureData::~DensityPressureData() {}

real DensityPressureData::GetDensity(real P) {
  int i;
  const double pace = 0.01;
  double x1, x2;
  if ((P < 1) || (log10(P) < 1)) {
    return Density[0];
  } else if (P < Pressure.back()) {
    i = trunc((log10(P) - 1) * 100);
    x1 = log10(P) - log10(Pressure[i]);
    x2 = pace - x1;
    return Density[i] * x2 / pace + Density[i + 1] * x1 / pace;
  } else {
    return Density.back();
  }
}
