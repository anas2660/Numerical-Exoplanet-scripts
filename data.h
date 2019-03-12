#ifndef DATA_H
#define DATA_H
#include "definitions.h"
#include <vector>

class DensityPressureData {
public:
  DensityPressureData(const string &filename);
  ~DensityPressureData();
  ℝ GetDensity(ℝ P);

private:
  vector<ℝ> Density, Pressure;
};
#define DATA DensityPressureData &
#endif
