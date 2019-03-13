#ifndef CSVWRITER_H
#define CSVWRITER_H
#include <fstream>
#include "definitions.h"
class CSVWriter {
 public:
  typedef valarray<ℝ> Row;
  CSVWriter(const string &filename, size_t columns);
  ~CSVWriter();
  void WritePoint(Row row, ℝ xmax, ℝ numpoints);
  void WritePoint(Row row);

 private:
  ofstream file;       // Out filestream
  vector<Row> buffer;  // Simplification buffer
  ℝ lastx;
  size_t rowsize;
  Row getAverage();
};
#endif
