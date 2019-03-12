#include "csvwriter.h"
#include <algorithm>
#include <assert.h>
#include <functional>
CSVWriter::CSVWriter(const string &filename, size_t columns)
    : file(filename), lastx(0), rowsize(columns) {
  assert(file.is_open());
}
CSVWriter::~CSVWriter() {
  if (buffer.size() > 0) // Write whats left
    WritePoint(getAverage());
  file.close();
}
CSVWriter::row CSVWriter::getAverage() {
  row avg(rowsize);
  for (auto row : buffer)
    avg += row;
  avg /= (ℝ)buffer.size();
  buffer.clear();
  return avg;
}
void CSVWriter::WritePoint(row Row, ℝ xmax, ℝ numpoints) {
  buffer.push_back(Row);
  const ℝ dx = xmax / numpoints;
  if (Row[0] >= lastx + dx) {
    WritePoint(getAverage());
    lastx += dx;
  }
}
void CSVWriter::WritePoint(row Row) {
  file << Row[0];
  if (rowsize > 1)
    for_each(begin(Row) + 1, end(Row), [=](ℝ &n) { file << "," << n; });
  file << endl;
}
