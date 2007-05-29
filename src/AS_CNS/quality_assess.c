
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
/* RCS info
 * $Id: quality_assess.c,v 1.7 2007-05-29 10:54:28 brianwalenz Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "assert.h"
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "MultiAlignStore_CNS.h"
#include "AS_PER_gkpStore.h"

int main(int argc, char *argv[])
{
FILE * hist=fopen("coverage_histogram.cgm","w");
int i;
int contig;
uint32 column;
int cov;
int bcov;
int ecov;
char call;
int qv;
int mm;
uint32 ungapped;
double expected=0.0; 
double ext_expected=0.0; 
uint32       zbases=0;
uint32       vbases=0;
uint32       ebases=0;
double EPROB[61]; // prob of error for each quality value
double PROB[61];  // prob of correct call for each quality value (1-eprob)

{ int qv=0;
  for (i=0;i<61;i++) {
     EPROB[i]= pow(10,-qv/10.);
     PROB[i] = (1.0 - EPROB[i]);
     //fprintf(stderr,"%d %g %g\n",i,EPROB[i],PROB[i]);
     qv++;
   }
}

fprintf(hist,"Coverage depth\n");
while( stdin ) {
        if ( fscanf(stdin,"%d %u %d %d %d %c %d %d %u", &contig,
              &column,&cov,&bcov,&ecov,&call,&qv,&mm,&ungapped) != 9 ) break;
        if ( qv == 0 ) {
          zbases++; // zero qv 
        } else if ( cov - ( bcov+ecov) == 0 )  {
          ext_expected+=EPROB[qv];    
          ebases++; // external data only
        } else {
          expected+=EPROB[qv];    
          fprintf(hist,"%d\n",cov);
          vbases++;
        } 
      }

 fprintf(stdout," expected errors = %g\n", expected);
 fprintf(stdout," expected errors in non-Celera columns = %g\n", ext_expected);
 fprintf(stdout," zero quality columns = %d\n", zbases);
 fprintf(stdout," external only columns = %d\n", ebases);
 fprintf(stdout," columns counting to erate = %d\n", vbases);
 fprintf(stdout," celera_erate = (expected errors) / (bases) = (%g) / (%d) = %.6f\n\n",expected,vbases,expected/vbases);
 fprintf(stdout," all_erate = (expected errors + ext_expected errors) / (bases+ebases) = (%g + %g) / (%d + %d) = %.6f\n\n",expected,ext_expected,vbases,ebases,(expected+ext_expected)/(vbases+ebases));

 exit (0);
}
