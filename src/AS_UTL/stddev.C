
static const char *rcsid = "$Id$";

#include "stddev.H"

#include <algorithm>

void
computeStdDev(vector<double> dist, double &mean, double &stddev) {
  mean   = 0.0;
  stddev = 0.0;

  if (dist.size() == 0) {
    //fprintf(stderr, "NO SAMPLES\n");
    return;
  }

  sort(dist.begin(), dist.end());

  double median     = dist[dist.size() * 1 / 2];
  double oneThird   = dist[dist.size() * 1 / 3];
  double twoThird   = dist[dist.size() * 2 / 3];

  double approxStd  = max(median - oneThird, twoThird - median);

  double biggest    = median + approxStd * 5;
  double smallest   = median - approxStd * 5;

  size_t numSamples = 0;

  for (size_t x=0; x<dist.size(); x++)
    if ((smallest  <= dist[x]) &&
        (dist[x]   <= biggest)) {
      numSamples += 1;
      mean       += dist[x];
    }

  if (numSamples == 0) {
    fprintf(stderr, "NO SAMPLES in range %f - %f\n", smallest, biggest);
    return;
  }

  mean   = mean / numSamples;
  stddev = 0.0;

  for (uint64 x=0; x<dist.size(); x++)
    if ((smallest  <= dist[x]) &&
        (dist[x]   <= biggest))
      stddev += (dist[x] - mean) * (dist[x] - mean);
          
  if (numSamples > 1)
    stddev = sqrt(stddev / (numSamples - 1));
}


double
computeExponentialMovingAverage(double alpha, double ema, double value) {

  assert(0.0   <= alpha);
  assert(alpha <= 1.0);

  return(alpha * value + (1 - alpha) * ema);
}
