#ifndef DATA_H
#define DATA_H
#include "definitions.h"
#include <vector>

class DensityPressureData {
public:
  DensityPressureData(const char *filename);
  ~DensityPressureData();

  real GetDensity(real P);
  // real GetPressure(real d);

private:
  vector<real> Density, Pressure;
};

#define DATA DensityPressureData &
#endif
