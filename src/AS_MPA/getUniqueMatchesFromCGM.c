
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
/* $Id: getUniqueMatchesFromCGM.c,v 1.6 2008-03-18 07:02:45 brianwalenz Exp $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HUMAN
#ifdef HUMAN
#define NUM_MATED_FRAGS 21525272
#else
#define NUM_MATED_FRAGS  2309852
#endif

int main(int argc, char ** argv)
{
  char line[1024], line2[1024];
  int id, bp, chr, matchbp, matchpct;
  long int left, right;
  char string[1000];
  FILE * fp, * fp2;
  unsigned long doThisID[NUM_MATED_FRAGS];
  int rejectedCount = 0;
  int acceptedCount = 0;
  int multiCount = 0;
  int singleCount = 0;
  int zeroCount = 0;

  if(argc != 4) exit(1);

  // celegram file
  fp = fopen(argv[1], "r");
  if(fp == NULL) exit(1);

  // uid file
  fp2 = fopen(argv[2], "r");
  if(fp2 == NULL) exit(1);

  fprintf(stderr, "Reading count file & uid file\n");
  for(id = 0, fgets(line, 1023, fp); fgets(line, 1023, fp); id++)
  {
    int mappingCount = atoi(line);
    fgets(line2, 1023, fp2);
    switch(mappingCount)
    {
      case 0:
        doThisID[id] = 0;
        zeroCount++;
        break;
      case 1:
        sscanf(line2, "%*d %lu", &(doThisID[id]));
        singleCount++;
        break;
      default:
        doThisID[id] = 0;
        multiCount++;
        break;
    }
  }
  fclose(fp);
  fclose(fp2);

  fp = fopen(argv[3], "r");
  if(fp == NULL) exit(1);
  fprintf(stderr, "Writing uniquely mapped fragment list\n");
  while(fgets(line, 1023, fp))
  {
    if(doThisID[atoi(line)] != 0)
    {
      sscanf(line, "%d[%d-%*d-%*d] %d[%ld-%ld] <%d-%*d-%d-%s",
             &id, &bp, &chr, &left, &right, &matchbp, &matchpct, string);
      if((right-left-matchbp)/((double) matchbp) < .1)
      {
        if(strncmp(string, "complement", 10) == 0)
        {
          fprintf(stdout, "%lu %d %d %ld %ld %d %d %d (%d)\n",
                  doThisID[id], id, chr, right, left, bp, matchbp, matchpct,
                  right-left-matchbp);
        }
        else if(strncmp(string, "forward", 7) == 0)
        {
          fprintf(stdout, "%lu %d %d %ld %ld %d %d %d (%d)\n",
                  doThisID[id], id, chr, left, right, bp, matchbp, matchpct,
                  right-left-matchbp);
        }
        else
        {
          fprintf(stderr, "ERROR!\n%s\n", string);
          exit(1);
        }
        acceptedCount++;
      }
      else
      {
        rejectedCount++;
      }
    }
  }
  fclose(fp);

  fprintf(stderr, "%d unique fragment mappings total.\n", singleCount);
  fprintf(stderr, "  %d accepted.\n", acceptedCount);
  fprintf(stderr, "  %d rejected - lengths don't match\n", rejectedCount);
  fprintf(stderr, "%d fragments had no mapping\n", zeroCount);
  fprintf(stderr, "%d fragments had multiple mappings\n", multiCount);
          
  return 0;
}
