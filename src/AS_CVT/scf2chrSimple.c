
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
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"

int main(int argc, char ** argv)
{
  FILE * fp1;
  FILE * fp2;
  
  if(argc != 3)
  {
    fprintf(stderr, "Usage: %s fileWfield1ScfUIDs  fileWfield2ChrNums\n",
            argv[0]);
    exit(-1);
  }
  fp1 = fopen(argv[1], "r");
  fp2 = fopen(argv[2], "r");
  assert(fp1 != NULL && fp2 != NULL);


  {
    char line1[1024];
    CDS_UID_t uid1;
    CDS_UID_t scfUID1;
    CDS_COORD_t fiveP1;
    CDS_COORD_t threeP1;
    int   numInstances;
    
    char line2[1024];
    CDS_UID_t scfUID2 = 0;
    char chrom[1024];
    CDS_COORD_t start2;
    CDS_COORD_t end2;
    int orient2;

    scfUID2 = 0;
    while(fgets(line1, 1023, fp1))
    {
      sscanf(line1, F_UID " " F_UID " " F_COORD " " F_COORD " %d",
             &scfUID1, &uid1, &fiveP1, &threeP1, &numInstances);

      while(scfUID2 < scfUID1)
      {
        if(fgets(line2, 1023, fp2))
          sscanf(line2, F_UID " %s " F_COORD " " F_COORD " %d",
                 &scfUID2, chrom, &start2, &end2, &orient2);
        else
          return 0;
      }
      if(scfUID2 > scfUID1) continue;
      assert(scfUID2 == scfUID1);

      switch(orient2)
      {
        case 1: // scaffold is oriented A_B
          if(fiveP1 < threeP1) // interval is A_B
            fprintf(stdout, "%s " F_UID " " F_COORD " " F_COORD " 1\n",
                    chrom, uid1,
                    start2 + fiveP1, threeP1 - fiveP1);
          else
            fprintf(stdout, "%s " F_UID " " F_COORD " " F_COORD " -1\n",
                    chrom, uid1,
                    start2 + threeP1, fiveP1 - threeP1);
          break;
        case -1: // scaffold is oriented B_A
          if(fiveP1 < threeP1)
            fprintf(stdout, "%s " F_UID " " F_COORD " " F_COORD " -1\n",
                    chrom, uid1,
                    end2 - threeP1, threeP1 - fiveP1);
          else
            fprintf(stdout, "%s " F_UID " " F_COORD " " F_COORD " 1\n",
                    chrom, uid1,
                    end2 - fiveP1, fiveP1 - threeP1);
          break;
        default: // scaffold is unoriented
          break;
      }
    }
  }
  
  return 0;
}
