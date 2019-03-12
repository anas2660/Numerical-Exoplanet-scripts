#ifndef CSVWRITER_H
#define CSVWRITER_H
#include "definitions.h"
#include <fstream>
#include <valarray>
#include <vector>

class CSVWriter {
public:
  typedef valarray<ℝ> row;
  CSVWriter(const string &filename, size_t columns);
  ~CSVWriter();
  void WritePoint(row Row, ℝ xmax, ℝ numpoints);
  void WritePoint(row Row);

private:
  ofstream file;      // Out filestream
  vector<row> buffer; // Simplification buffer
  ℝ lastx;
  size_t rowsize;
  row getAverage();
};
#endif
