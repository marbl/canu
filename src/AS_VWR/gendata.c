
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
/* $Id: gendata.c,v 1.4 2005-03-22 19:49:31 jason_miller Exp $ */

#include <stdio.h>
#include <stdlib.h>

#define NUM_FRAGS     50
#define MAX_LEN       650
#define MIN_LEN       350
#define RANGE       10000

int main(int argc, char *argv[])
{ int i, beg, len;

  srand48(getpid());
  printf("0: CFF0000 T2  # Red\n");
  printf("1: C00AF00     # Green\n");
  printf("2: CA0A0FF T3  # Blue Mark\n");
  printf("3: CFFAA00 T1  # Linker\n");
  for (i = 0; i < NUM_FRAGS; i++)
    { len = drand48()*(MAX_LEN-MIN_LEN) + MIN_LEN;
      beg = drand48()*(RANGE - len); 
      if (drand48() > .75)
        printf("%d: %d A0 %d  # Comment %d\n",i,beg,beg+len,i);
      else if (drand48() > .66)
        printf("%d: %d A1 %d  # Comment %d\n",i,beg,beg+len,i);
      else if (drand48() > .50)
        printf("%d: %d A2 M%d A0 %d  # Comment %d\n",
               i,beg,beg+len/3,beg+len,i);
      else
        printf("%d: %d A2 M%d A2 M%d A1 %d  # Comment %d\n",i,
                beg,beg+len/3,beg+(2*len)/3,beg+len,i);
    }
  printf("ATT: 3\n");
  for (i = NUM_FRAGS*.5/3.5; i >= 0; i--)
    { int size, chosen[6];
      int k, j, min;

      size = drand48()*4 + 2;
      printf("LNK: ");
      for (j = 0; j < size; j++)
        chosen[j] = drand48()*(NUM_FRAGS-size);
      for (k = 0; k < size-1; k++)
        { min = k;
          for (j = k+1; j < size; j++)
            if (chosen[j] < chosen[min])
              min = j;
          if (min != k)
            { j = chosen[k];
              chosen[k] = chosen[min];
              chosen[min] = j;
            }
        }
      for (k = 0; k < size; k++)
        printf(" %d",chosen[k]+k);
      printf("   # Link set %d\n",k);
    }
}
