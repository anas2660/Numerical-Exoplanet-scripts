#include "csvwriter.h"
#include <assert.h>

CSVWriter::CSVWriter(const char* filename) : file(filename), lastx(0){
  assert(file.is_open());
}

CSVWriter::~CSVWriter(){
  // Write whats left
  // Calculate average
  Point avg = getAverage();
  // Write average
  file << avg.x << "," << avg.y << "," << avg.z << std::endl;
  buffer.clear();
  file.close();
}


CSVWriter::Point CSVWriter::getAverage(){
  Point avg = {0, 0, 0};
  for (auto it = buffer.begin(); it != buffer.end(); ++it) {
    avg.x += it->x;
    avg.y += it->y;
    avg.z += it->z;
  }
  avg.x /= (real)buffer.size();
  avg.y /= (real)buffer.size();
  avg.z /= (real)buffer.size();
  return avg;
}


void CSVWriter::WritePointSimplify(ℝ x, ℝ y, ℝ dy, ℝ xmax, ℝ numpoints){
  buffer.push_back({x, y, dy});
  const real dx = xmax / numpoints;
  if(x >= lastx+dx){
    // Calculate average
    Point avg = getAverage();

    // Write average
    file << avg.x << "," << avg.y << "," << avg.z << std::endl;
    buffer.clear();
    lastx += dx;
  }
}
