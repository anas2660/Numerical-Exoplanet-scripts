#include "data.h"
#include <fstream>
#include <iostream>
#include <string>

DensityPressureData::DensityPressureData(const string &filename) {
  ifstream file("data/" + filename);
  assert(file.is_open());
  ℝ d, p;
  while (file >> d >> p) {
    Density.emplace_back(d);
    Pressure.emplace_back(p);
  }
  file.close();
}
DensityPressureData::~DensityPressureData() {}
ℝ DensityPressureData::GetDensity(ℝ P) {
  int i;
  ℝ step = 0.01, x1;
  if ((P < 1) || (log10(P) < 1))
    return Density[0];
  else if (P < Pressure.back()) {
    i = trunc((log10(P) - 1.0) * 100.0);
    x1 = log10(P) - log10(Pressure[i]);
    return Density[i] + (Density[i + 1] - Density[i]) * x1 / step;
  } else
    return Density.back();
}
