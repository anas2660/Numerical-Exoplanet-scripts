#ifndef CSVWRITER_H
#define CSVWRITER_H

#include "definitions.h"
#include <fstream> // File streams
#include <vector> // Dynamic lists

class CSVWriter
{
public:
  CSVWriter(const char* filename);
  ~CSVWriter();

  void WritePointSimplify(ℝ x, ℝ y, ℝ dy, ℝ xmax, ℝ numpoints);

private:
  struct Point{ real x, y, z; }; // 3D point structure

  ofstream file; // Out-filestream
  vector<Point> buffer; // Point buffer
  real lastx;

  Point getAverage();

};

#endif
