
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

/*
  Program to concatenate 9 output files from Brian & renumber in the process
  
run1.regions
run2.regions
run3.regions
run4.regions
run5.regions
run6.regions
run7.regions
run8a.regions
run8b.regions

1:         0  2690659 
2:   2690660  5381318 
3:   5381319  8071977 
4:   8071978 10762636 
5:  10762637 13453295 
6:  13453296 16143954 
7:  16143955 18834613 
8:  18834614 21525271 

8a:        0 1345328 
8b:  1345329 2690657 

8a: 18834614 20179942
8b: 20179943 21525271
*/

#define NUM_FILES   8
char * files[NUM_FILES] =
{
 "run1.regions",
 "run2.regions",
 "run3.regions",
 "run4.regions",
 "run5.regions",
 "run6.regions",
 "run7.regions",
 "run8.regions"
};

int offsets[NUM_FILES] =
{
       0,
 2690660,
 5381319,
 8071978,
10762637,
13453296,
16143955,
18834614
};

int main(int argc, char ** argv)
{
  int i;
  int id, bp, chr, matchbp, matchpct;
  long int left, right;
  char string[1000];
  int j1, j2, j3;
  char line[2048];
  FILE * fp;
  int numLines;
  
  for(i = 0; i < NUM_FILES; i++)
  {
    fprintf(stderr, "Working on file %s\n", files[i]);

    numLines = 0;
    fp = fopen(files[i], "r");
    while(fgets(line, 2047, fp))
    {
      sscanf(line, "%d[%d-%d-%d] %d[%ld-%ld] <%d-%d-%d-%s",
             &id, &bp, &j1, &j2,
             &chr, &left, &right,
             &matchbp, &j3, &matchpct, string);

      fprintf(stdout, "%d[%d-%d-%d] %d[%ld-%ld] <%d-%d-%d-%s\n",
              id + offsets[i], bp, j1, j2,
              chr, left, right,
              matchbp, j3, matchpct, string);
      if(numLines == 0)
      {
        fprintf(stderr, "First line input:\n%s\n", line);
        fprintf(stderr, "First line output:\n");
        fprintf(stderr, "%d[%d-%d-%d] %d[%ld-%ld] <%d-%d-%d-%s\n",
                id + offsets[i], bp, j1, j2,
                chr, left, right,
                matchbp, j3, matchpct, string);
      }
      numLines++;
    }
    fprintf(stderr, "Last line input:\n%s\n", line);
    fprintf(stderr, "Last line output:\n");
    fprintf(stderr, "%d[%d-%d-%d] %d[%ld-%ld] <%d-%d-%d-%s\n",
            id + offsets[i], bp, j1, j2,
            chr, left, right,
            matchbp, j3, matchpct, string);
    fclose(fp);
  }
  return 0;
}
