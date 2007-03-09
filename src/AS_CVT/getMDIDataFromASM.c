
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
#include <string.h>
#include <unistd.h>
#include "AS_global.h"
#include "AS_MSG_pmesg.h"


void Usage(char * progName)
{
  fprintf(stderr, "Usage: %s  [-a filename]  -s   -h\n", progName);
  fprintf(stderr,
          "\t-a filename     .asm filename\n"
          "\t-s              use stddev as bucket size\n"
          "\t                (default is cgw's bucket size)\n"
          "\t-h              print usage\n"
          );
  exit(-1);
}


int main( int argc, char ** argv )
{
  GenericMesg * gen;
  int inMDIs = 0;
  int byBuckets = 1;
  char * asmFilename = NULL;
  FILE * fi;
  
  // parse command line
  {
    int ch, errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "a:sh")) != EOF))
    {
      switch(ch) 
      {
	case 'a':
          asmFilename = optarg;
          break;
        case 's':
          byBuckets = 0;
          break;
        case 'h':
          Usage(argv[0]);
          break;
        default:
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  errflg++;
          break;
      }
    }
  }
  
  if(asmFilename == NULL)
    Usage(argv[0]);

  fi = fopen(asmFilename, "r");
  if(fi == NULL)
  {
    fprintf(stderr, "Failed to open %s for reading.\n", asmFilename);
    exit(-1);
  }
  
  while(ReadProtoMesg_AS(fi, &gen) != EOF)
  {
    switch(gen->t)
    {
      case MESG_MDI:
      {
        SnapMateDistMesg * smdm = gen->m;

        if(smdm->num_buckets > 0)
        {
          FILE * fout;
          char filename[1024];
          int i;
          float bucketSize = smdm->stddev / 3.;

          // set up a file for this mdi
          sprintf(filename, F_UID ".mdi", smdm->erefines);
          fout = fopen(filename, "w");

          fprintf(fout, "(" F_UID "," F_IID ")\t%f\t%f\t(%d,%d)",
                  smdm->erefines, smdm->irefines,
                  smdm->mean, smdm->stddev,
                  smdm->min, smdm->max);

          if(byBuckets == 1)
          {
            
            fprintf(fout, "\t%d\n", smdm->num_buckets);
            
            for(i = 0; i < smdm->num_buckets; i++)
            {
              fprintf(fout, "%d\t%d\n",
                      (int) (i * bucketSize + smdm->min + .5),
                      smdm->histogram[i]);
            }
          }
          else
          {
            int * sc = NULL;
            float minStddev = (smdm->min - smdm->mean) / smdm->stddev;
            float maxStddev = (smdm->max - smdm->mean) / smdm->stddev;
            
            fprintf(fout, "\t%d\n", (int) (.5+maxStddev + 1 - minStddev));
            
            sc = (int *) safe_calloc((int) (maxStddev + 1.5 - minStddev), sizeof(int));
            for(i = 0; i < smdm->num_buckets; i++)
            {
              sc[(int) ((i * bucketSize + smdm->min - smdm->mean) /
                        smdm->stddev - minStddev)] += smdm->histogram[i];
            }
            
            for(i = 0; i < maxStddev + 1 - minStddev; i++)
              fprintf(fout, "%d\t%d\n",
                      (int) (.5 + smdm->mean + smdm->stddev * (minStddev + i)),
                      sc[i]);

            safe_free(sc);
          }

          fclose(fout);
        }
        inMDIs = 1;
      }
      break;
      default:
        if(inMDIs)
          exit(0);
        break;
    }
  }
  fclose(fi);

  return 0;
}
