
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "AS_UTL_rand.h"

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
