
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-FEB-09 to 2015-JUN-02
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "Binomial_Bound.H"
#include "gkStore.H"

#undef COMPUTE_IN_LOG_SPACE

//  Determined by  EDIT_DIST_PROB_BOUND
#define  NORMAL_DISTRIB_THOLD    3.62

//  Probability limit to "band" edit-distance calculation
//  Determines  NORMAL_DISTRIB_THOLD
#define  EDIT_DIST_PROB_BOUND    1e-4



//  Return the smallest  n >= Start  s.t.
//    prob [>= e  errors in  n  binomial trials (p = error prob)] > EDIT_DIST_PROB_BOUND
//
int
Binomial_Bound(int e, double p, int Start) {
  double  Normal_Z, Mu_Power, Factorial, Poisson_Coeff;
  double  q, Sum, P_Power, Q_Power, X;
  int  k, n, Bin_Coeff, Ct;

  q = 1.0 - p;
  if (Start < e)
    Start = e;

  for (n = Start;  n < AS_MAX_READLEN;  n ++) {
    if (n <= 35) {
      Sum = 0.0;
      Bin_Coeff = 1;
      Ct = 0;
      P_Power = 1.0;
      Q_Power = pow (q, n);

      for (k = 0;  k < e && 1.0 - Sum > EDIT_DIST_PROB_BOUND;  k ++) {
        X = Bin_Coeff * P_Power * Q_Power;
        Sum += X;
        Bin_Coeff *= n - Ct;
        Bin_Coeff /= ++ Ct;
        P_Power *= p;
        Q_Power /= q;
      }

      if (1.0 - Sum > EDIT_DIST_PROB_BOUND)
        return(n);

    } else {
      Normal_Z = (e - 0.5 - n * p) / sqrt (n * p * q);
      if (Normal_Z <= NORMAL_DISTRIB_THOLD)
        return  n;

#ifndef COMPUTE_IN_LOG_SPACE
      Sum = 0.0;
      Mu_Power = 1.0;
      Factorial = 1.0;
      Poisson_Coeff = exp (- n * p);
      for (k = 0;  k < e;  k ++) {
        Sum += Mu_Power * Poisson_Coeff / Factorial;
        Mu_Power *= n * p;
        Factorial *= k + 1;
      }
#else
      Sum = 0.0;
      Mu_Power = 0.0;
      Factorial = 0.0;
      Poisson_Coeff = - n * p;
      for (k = 0;  k < e;  k ++) {
        Sum += exp(Mu_Power + Poisson_Coeff - Factorial);
        Mu_Power += log(n * p);
        Factorial = lgamma(k + 1);
      }
#endif

      if (1.0 - Sum > EDIT_DIST_PROB_BOUND)
        return(n);
    }
  }

  return(AS_MAX_READLEN);
}




void
Initialize_Match_Limit(int32 *ml, double maxErate, int32 maxErrors) {
  int32 e = 0;
  int32 s = 1;
  int32 l = MIN(maxErrors, 2000);  //  Compute the first 2000 values; set to maxErrors to do no estimation

  //  The number of errors that are ignored in setting probability bound for terminating alignment
  //  extensions in edit distance calculations
  int32 ERRORS_FOR_FREE = 1;

  //  Free errors.

  while (e <= ERRORS_FOR_FREE)
    ml[e++] = 0;

  //  Compute the actual limits.  This is _VERY_ expensive for longer reads.  BITS=17 is about all
  //  it can support.

#ifdef DUMP_MATCH_LIMIT
  l = maxErrors;  //  For showing deviations
#endif

  while (e < l) {
    s = Binomial_Bound(e - ERRORS_FOR_FREE, maxErate, s);
    ml[e] = s - 1;

    assert(ml[e] >= ml[e-1]);

    //if ((e % 100) == 0)
    //  fprintf(stderr, " %8.4f%% - %8d / %8d\r", 100.0 * e / maxErrors, e, maxErrors);

    e++;
  }

  //  Estimate the remaining limits.  using a linear function based on a precomputed slope.
  //  prefixEditDistance-matchLimitGenerate computes the data values for a bunch of error rates.
  //  These are used to compute the slope of a line from the [2000] point through the [max] point.
  //  These slopes fit, almost exactly, an a/x+b curve, and that curve is used to compute the slope
  //  for any error rate.
  //
#if AS_MAX_READLEN_BITS == 17
  double sl = 0.962830901135531 / maxErate + 0.096810267016486;
#endif

#if AS_MAX_READLEN_BITS == 18
  double sl = 0.964368146781421 / maxErate + 0.118101522100597;
#endif

#if AS_MAX_READLEN_BITS == 19
  double sl = 0.963823337297648 / maxErate + 0.156091528250625;
#endif

#if AS_MAX_READLEN_BITS == 20
  double sl = 0.971023157863311 / maxErate + 0.154425731746994;
#endif

#if AS_MAX_READLEN_BITS == 21
  double sl = 0.982064188397525 / maxErate + 0.067835741959926;
#endif

#if AS_MAX_READLEN_BITS == 22
  double sl = 0.986446300363063 / maxErate + 0.052358358862826;
#endif

  //  And the first value.
  double  vl = ml[e-1] + sl;

  while (e < maxErrors) {
    ml[e] = (int32)ceil(vl);
    vl += sl;
    e++;
  }

#ifdef DUMP_MATCH_LIMIT
  FILE *F = fopen("values-new.dat", "w");
  for (int32 e=0; e<maxErrors; e++)
    fprintf(F, "%d %d\n", e, Edit_Match_Limit[e]);
  fclose(F);

  fprintf(stderr, "values-orig.dat and values-new.dat dumped.  exiting.\n");
  exit(1);
#endif

}

