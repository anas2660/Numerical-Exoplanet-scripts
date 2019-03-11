#include "definitions.h"
#include "exoterDE.cpp"

int main() {
  // Planet 1
  ExoterDE("CoRoT7b", 4.8, 0.8, 1.68, 0.09);
  // Planet 2
  ExoterDE("GJ1214b", 6.55, 0.98, 2.678, 0.13);
  if (PLOT)
    system("python ternaryplot.py CoRoT7b GJ1214b");
}
