
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

