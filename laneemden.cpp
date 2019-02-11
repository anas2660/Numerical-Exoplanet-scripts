#include <cmath>
#include <fstream>
using namespace std;
typedef double real;
void WritePoint(real x, real y, real dy);

real n = 5;

ofstream  outfile;

real thetaf(real xi){
  return 1.0/sqrt(1.0+(xi*xi/3.0));
}

void EulerIntegration(double n){
  real theta = 1.0, xi = 0, dxi = 0.01, ximax = 10.0, v = 0, dv = 0;
  WritePoint(xi, theta, thetaf(xi+dxi) - theta);
  for (xi = dxi; xi < ximax; xi += dxi) {
    v += -((xi*xi)*pow(theta, n)*dxi);
    theta += v*(1.0/(xi*xi))*dxi;
    WritePoint(xi, theta, thetaf(xi)-theta);
  }
}










int main(){
  outfile = ofstream("eulerint.csv");
  if(!outfile.is_open())
    throw "Could not open file.";
  try {
    EulerIntegration(n);
  } catch (const char* error) {
    printf("%s\n", error);
  }
  outfile.close();
}


void WritePoint(double x, double y, double dy){
  outfile << x << "," << y << "," << dy <<  endl;
}
