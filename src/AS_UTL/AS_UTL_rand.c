
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
static char CM_ID[] = "$Id: AS_UTL_rand.c,v 1.4 2005-03-22 19:49:29 jason_miller Exp $";
#include <stdio.h>
#include <math.h>
#include "AS_UTL_rand.h"

static int32 iset;
static double gset;


double GaussRandomNormalized_AS(void){
  double fac, r, v1, v2;
  double ret;

  if(iset == 0){
    do{
      v1 = 2.0 * drand48()  -1.0; /* Uniform on -1,1 */
      v2 = 2.0 * drand48()  -1.0; /* Uniform on -1,1 */
      /* v1 and v2 are random in the square (-1,-1) (1,1) */
      /* If the sum of their squares is < 1 then they are in the unit circle */
      r = v1 * v1  + v2 * v2;
    } while(r >= 1.0);
    //    fprintf(stderr,"* v1 = %g v2 = %g r = %g\n",
    //    v1,v2,r);
    if(r > 1.0e-09)
      fac = sqrt(-2 * log(r)/r);
    else
      fac = 0.0;
    /* Now make the Box-Muller transformation to get two normal deviates.  Save
       one and return the other */
    gset = v1 * fac;
    iset = 1;
    ret =  v2 * fac;
  }else{
    iset = 0;
    ret = gset;
  }
  //  fprintf(stderr,"* Returning %g\n", ret);
  return ret;
}


/* **************************************************** */

int GetRand_AS(int min, int max, int uniform){
  // int retValue;

  if(min == max)
    return min;

  if(uniform)
    return  (int)(min + (max - min + 0.9999) * drand48());

  {
    double mean = (double)(min + max)/2.0;
    double std = (double)(max - min)/6.0; /* Use 99% interval, roughly 3 std devs, symmetrically */
    double randnum;

    
    /* Filter out values that are > 5.0 std from mean */
    randnum = GaussRandom_AS(mean,std);
    while(fabs((randnum-mean)/std) > 5.0){
      fprintf(stderr,"* rand %g was > 5stds from mean ...trying again\n",
	      randnum);
	randnum = GaussRandom_AS(mean,std);
    }
    return (int)randnum;
  }
}



double GetDrand_AS(double min, double max)
{
  double l = max-min;

  if( min == max )
    return min;
  
  return(drand48()*l+min);
}
  
