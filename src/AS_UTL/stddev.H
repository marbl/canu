#ifndef STDDEV_H
#define STDDEV_H

#include "AS_global.H"

#include <vector>

using namespace std;

void
computeStdDev(vector<double> dist, double &mean, double &stddev);

double
computeExponentialMovingAverage(double alpha, double ema, double value);

#endif

