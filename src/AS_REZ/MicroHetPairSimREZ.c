
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
/**********************************************************************
 Module:      The functions in this file are only for testing purposes
 Description:
 Assumptions:
**********************************************************************/

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "UtilsREZ.h"
#include "AS_UTL_rand.h"
#include "AS_UTL_skiplist.h"
#include "MicroHetREZ.h"
#include "MicroHetScoreREZ.h"
#include "MicroHetPartitionsREZ.h"


/* poor man's version of Bernard Sulzer's mbs macros */

/*********************************************************************
**
**   MACROS:        
**
**   ROUTINES:
**                  mbs_2D_alloc_array
**                  mbs_2D_free_array
**                 
**********************************************************************/
double RANDOMMAX;

double getrandomrange(void);

double getrandomrange(void){
  return( ((double)random())/RANDOMMAX);
}

static void **
mbs_2D_alloc_array(int dim1, int dim2, int ob_type)
{
        int            i;
	void           **ptr_2, *ptr_1;
	
	if (!(ptr_1 = (void *)calloc(dim1*dim2, ob_type))){
	    fprintf(stderr,"%s\n\t%s\n","MBS ALLOC ARRAY (2D)", "level 0 data");
	    exit(-1);
	}
	if (!(ptr_2 = (void **)calloc(dim1, sizeof(void *)))){
	    fprintf(stderr,"%s\n\t%s\n","MBS ALLOC ARRAY (2D)", "level 1 pointers");
	    exit(-1);
	}

	for (i=0 ; i<dim1 ; i++)
	{
	    ptr_2[i] = ptr_1;
	    ptr_1 = (void *) (dim2 * ob_type + (long) ptr_1);
	}
	return (ptr_2);
}

static void    
mbs_2D_free_array(void **ptr)
{
	free(*ptr);
	free(ptr);
}


int main(int argc,char *argv[])
{
  int rows=10,cols=60;
  double erate=.03;
  int i,j;
  double l;
  char **ali;
  Alignment_t aln;
  Marker_t m;

  ali=(char**)mbs_2D_alloc_array(rows,cols+1,sizeof(char));

  RANDOMMAX= pow(2,31);
  srandom(atoi(argv[1]));

  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      ali[i][j]='A';
    }
    ali[i][j]='\0';
  }
  for(i=0;i<rows;i++){
      for(j=0;j<cols;j++){
	l = getrandomrange();
	if( l < erate ){
	  l = getrandomrange();
	  if ( l < .25 )
	    ali[i][j]='C';
	  else if ( l < .5)
	    ali[i][j]='G';
	  else if ( l < .75)
	    ali[i][j]='T';
	  else
	    ali[i][j]='-';
	}
      }
      //      printf("%s\n",ali[i]);
  }
  aln.ali=ali;
  aln.rows=rows;
  aln.cols=cols;
  m.set=(int*)malloc(sizeof(int)*rows);
  if(m.set==NULL){
    fprintf(stderr,"Out of memory, pairwise test marker");
    exit(-1);
  }
  for(i=0;i<rows;i++){
    (m.set)[i]=TRUE;
  }
  m.len=rows;

  printf("%e\n",Pairwise_Only_Prob_Simple(&aln,0,cols-1,&m));

  free(m.set);
  free(ali);
  return(0);
}
