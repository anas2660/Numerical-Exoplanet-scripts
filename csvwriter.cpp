#include "csvwriter.h"
CSVWriter::CSVWriter(const string &filename, size_t columns)
    : file(filename), lastx(0), rowsize(columns) {
  assert(file.is_open());
}
CSVWriter::~CSVWriter() {  // Write whats left
  if (buffer.size() > 0) WritePoint(getAverage());
  file.close();
}
CSVWriter::Row CSVWriter::getAverage() {
  Row avg(rowsize);
  for (auto row : buffer) avg += row;
  avg /= (ℝ)buffer.size();
  buffer.clear();
  return avg;
}
void CSVWriter::WritePoint(Row row, ℝ xmax, ℝ numpoints) {
  buffer.push_back(row);
  const ℝ dx = xmax / numpoints;
  if (row[0] >= lastx + dx) {
    WritePoint(getAverage());
    lastx += dx;
  }
}
void CSVWriter::WritePoint(Row row) {
  file << row[0];
  if (rowsize > 1)
    for (int n = 1; n < row.size(); n++) file << "," << row[n];
  file << endl;
}
