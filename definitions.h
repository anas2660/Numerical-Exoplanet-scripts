#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define PLOT true
#define ERROR errorf(ξ + dξ, θ)
#define ½ (1.0 / 2.0) *
#define N 5
#define DKSI 0.001
#define PI 3.1415926535897932384626434
#define G 0.0000000000667408
#define ASTEP 0.005
#define Mearth 5.9736e24
#define Rearth 6372797.0
#define USING_BOOST 1
#include <cassert>
#include <cmath>
#include <thread>
#include <vector>
using namespace std;
typedef double real;
typedef double ℝ;

// Analytical solutions to calculate error
#if N == 1
inline ℝ θf1(ℝ ξ) { return sin(ξ) / ξ; }
inline ℝ errorf(ℝ ξ, ℝ θ) { return θf1(ξ) - θ; }
#elif N == 5
inline ℝ θf5(ℝ ξ) { return 1.0 / sqrt(1.0 + (ξ * ξ / 3.0)); }
inline ℝ errorf(ℝ ξ, ℝ θ) { return θf5(ξ) - θ; }
#else
inline ℝ errorf(ℝ ξ, ℝ θ) { return 0; }
#endif

#endif
