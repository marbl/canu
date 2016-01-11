
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
 *    Brian P. Walenz from 2008-JUN-27 to 2013-AUG-01
 *      are Copyright 2008,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz beginning on 2015-OCT-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "AS_UTL_rand.H"

#define ITERATIONS 1000000

int main(void){
  double
    sumx = 0.0,
    sumx2 = 0.0,
    sumx3 = 0.0;
  int i;
  int Seed = getpid();
  double num;
  srand48(Seed); /* Initialize Random Number Generator */

  for(i = 0; i < ITERATIONS; i++){
    double rand = GaussRandomNormalized_AS();
    double mult = rand;
    sumx += mult;
    mult *= rand;
    sumx2 += mult;
    mult *= rand;
    sumx3 += mult;

  }

  num = (double)ITERATIONS;

  fprintf(stderr,"***** Normalized Gaussian *****");
  fprintf(stderr,"* avg = %f  avg2 = %f  avg3 = %f\n",
          sumx/num, sumx2/num, sumx3/num);



  sumx = 0.0;
  sumx2 = 0.0;
  sumx3 = 0.0;


  for(i = 0; i < ITERATIONS; i++){
    double rand = GaussRandom_AS(0.0, 5.0);
    double mult = rand;
    sumx += mult;
    mult *= rand;
    sumx2 += mult;
    mult *= rand;
    sumx3 += mult;

  }

  num = (double)ITERATIONS;
  sumx = sumx/num;
  sumx2 = sqrt(sumx2/num);
  sumx3 = pow(sumx3/num, 0.33);
  fprintf(stderr,"*****  Gaussian with STDEV 5.0 *****");
  fprintf(stderr,"* avg = %f  sqrt avg2 = %f  curt avg3 = %f\n",
          sumx, sumx2, sumx3);


  return 0;
}
