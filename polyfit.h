#ifndef POLYFIT_H
#define POLYFIT_H
#include "definitions.h"
template <typename T>
vector<T> polyfit(const vector<T> &oX, const vector<T> &oY, int nDegree);
template <typename T>
vector<T> polyval(const vector<T> &oCoeff, const vector<T> &oX);
#endif
