/**************************************************************************
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
/* $Id: identifyMistrackedClones.c,v 1.4 2005-09-23 01:17:07 brianwalenz Exp $ */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define DEFAULT_NUM_SIGMAS       4


void usage(char * progname)
{
  fprintf(stderr, "Usage: %s [-h] -f lengthsFile  -l libFile  -n libNum\n"
          "\t[-s sigmas]  [-p]  [-z count]  [-o count]  [-t count]\n"
          "  -h               print help\n"
          "  -f lengthsFile   name of file containing clone length data\n"
          "                     default is stdin\n"
          "  -l libFile       file with mean & stdev of libs\n"
          "  -n libNum        UID of 'other' library in libFile to check against\n"
          "  -s sigmas        require lengths to be within sigmas from mean\n"
          "                     default is %u\n"
          "  -p               identify clones if 100%% agree with other lib\n"
          "  -z count         count instances of clone with zero disagreements identify it as mistracked\n"
          "  -o count         count instances of clone with one disagreement identify it as mistracked\n\n"
          "  -t count         count instances of clone with two disagreements identify it as mistracked\n\n"
          "Output is written to stdout.\n"
          "If clone has a disagreement, the disagreeing distance is written\n",

          progname, DEFAULT_NUM_SIGMAS);
  exit(1);
}


int main(int argc, char ** argv)
{
  char * dfn = NULL;
  char * lfn = NULL;
  FILE * dfp;
  FILE * lfp = stdin;
  char line[2048];
  unsigned long n;
  unsigned long long libNum;
  unsigned int minZCount = 0;
  unsigned int minOCount = 0;
  unsigned int minTCount = 0;
  int do100Pct = 0;
  double numSigmas = DEFAULT_NUM_SIGMAS;
  double val;
  double mean;
  double sigma;

  // scope for scope's sake
  {
    int ch;
    while((ch = getopt(argc, argv, "hf:l:n:s:z:o:t:p")) != EOF)
    {
      switch(ch)
      {
        case 'h':
          usage(argv[0]);
          break;
        case 'f':
          dfn = optarg;
          break;
        case 'l':
          lfn = optarg;
          break;
        case 'n':
          sscanf(optarg, "%llu", &libNum);
          break;
        case 's':
          numSigmas = atoi(optarg);
          break;
        case 'z':
          minZCount = atoi(optarg);
          break;
        case 'o':
          minOCount = atoi(optarg);
          break;
        case 't':
          minTCount = atoi(optarg);
          break;
        case 'p':
          do100Pct = 1;
          break;
        default:
          usage(argv[0]);
          break;
      }
    }
  }

  if(minZCount <= 1 || numSigmas <= 0 || minOCount <= 0)
    usage(argv[0]);
  
  // read lib info
  if(lfn == NULL || (lfp = fopen(lfn, "r")) == NULL)
    usage(argv[0]);

  while(fgets(line, 2047, lfp) != NULL)
  {
    unsigned long long ln;
    sscanf(line, "%llu %lf %lf", &ln, &mean, &sigma);
    if(ln == libNum)
      break;
  }
  fclose(lfp);

  // fprintf(stderr, "%llu %lf %lf\n", libNum, mean, sigma);

  // read lengths
  {
    double minLength = mean - numSigmas * sigma;
    double maxLength = mean + numSigmas * sigma;
    unsigned long long minUID;
    unsigned long long maxUID;
    unsigned long long lastMinUID = 0;
    unsigned long long lastMaxUID = 0;
    int numThisUID = 1;
    int numThisUIDOkay = 0;
    double printLength;
    double lastVal;

    fprintf(stderr, "Identifying clones between %lf and %lf bp with\n"
            "%d or more non-conflicting instances or with\n"
            "%d or more instances with at most one conflict or with\n"
            "%d or more instances with at most two conflicts\n",
            minLength, maxLength, minZCount, minOCount, minTCount);

    if(do100Pct)
      fprintf(stderr, "Or if 100%% of clones agree\n");
      
    if(dfn == NULL || (dfp = fopen(dfn, "r")) == NULL)
      usage(argv[0]);
    while(fgets(line, 2047, lfp) != NULL)
    {
      sscanf(line, "%lf %llu %llu", &val, &minUID, &maxUID);

      if(minUID == lastMinUID)
      {
        numThisUID++;
      }
      else
      {
        printLength = ((numThisUIDOkay != numThisUID) ? printLength : lastVal);
        if(((numThisUIDOkay >= minZCount || do100Pct) &&
            numThisUIDOkay == numThisUID) ||
           (numThisUIDOkay >= minOCount && numThisUID - numThisUIDOkay == 1) ||
           (numThisUIDOkay >= minTCount && numThisUID - numThisUIDOkay == 2))
        {
          fprintf(stdout, "%lf %llu %llu %d %d\n",
                  printLength, lastMinUID, lastMaxUID,
                  numThisUIDOkay, numThisUID - numThisUIDOkay);
        }
        numThisUID = 1;
        numThisUIDOkay = 0;
      }
      
      if(val < minLength || val > maxLength)
      {
        printLength = val;
      }
      else
      {
        numThisUIDOkay++;
      }

      lastVal = val;
      lastMinUID = minUID;
      lastMaxUID = maxUID;
    } // lines in file
    fclose(lfp);

    printLength = ((numThisUIDOkay != numThisUID) ? printLength : lastVal);
    if(((numThisUIDOkay >= minZCount || do100Pct) &&
        numThisUIDOkay == numThisUID) ||
       (numThisUIDOkay >= minOCount && numThisUID - numThisUIDOkay == 1) ||
       (numThisUIDOkay >= minTCount && numThisUID - numThisUIDOkay == 2))
    {
      fprintf(stdout, "%lf %llu %llu %d %d\n",
              printLength, lastMinUID, lastMaxUID,
              numThisUIDOkay, numThisUID - numThisUIDOkay);
    }
  } // scope

  return 0;
}

