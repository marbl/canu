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
/* $Id: reestimateLibs.c,v 1.3 2005-09-21 20:13:07 catmandew Exp $ */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <values.h>
#include <math.h>

#define DEFAULT_NUM_ITERATIONS   4
#define DEFAULT_NUM_SIGMAS       4


void usage(char * progname)
{
  fprintf(stderr, "Usage: %s [-h] -f lengthsFile  -l libFile  -n libNum  -s sigmas  -i iterations\n"
          "  -h               print help\n"
          "  -v               be verbose to stderr\n"
          "  -f lengthsFile   name of file containing clone length data\n"
          "                     default is stdin\n"
          "  -l libFile       file with initial mean & stdev of libs\n"
          "  -n libNum        UID of library (in libFile)\n"
          "  -s sigmas        exlude lengths more than sigmas from mean\n"
          "                     default is %u\n"
          "  -i iterations    number of iterations to perform\n"
          "                     default is %u\n"
          "\n\n"
          "New libFile-formatted output is written to stdout\n\n",
          progname, DEFAULT_NUM_SIGMAS, DEFAULT_NUM_ITERATIONS);
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
  unsigned int numIterations = DEFAULT_NUM_ITERATIONS;
  double numSigmas = DEFAULT_NUM_SIGMAS;
  double sumX;
  double sumX2;
  double val;
  double mean;
  double sigma;
  double meanOrig;
  double sigmaOrig;
  int verbose = 0;

  // scope for scope's sake
  {
    int ch;
    while((ch = getopt(argc, argv, "hvf:l:n:s:i:")) != EOF)
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
        case 'i':
          numIterations = atoi(optarg);
          break;
        case 'v':
          verbose = 1;
          break;
        default:
          usage(argv[0]);
          break;
      }
    }
  }

  // read lib info
  if(lfn == NULL || (lfp = fopen(lfn, "r")) == NULL)
    usage(argv[0]);

  while(fgets(line, 2047, lfp) != NULL)
  {
    unsigned long long ln;
    sscanf(line, "%llu %lf %lf", &ln, &meanOrig, &sigmaOrig);
    if(ln == libNum)
      break;
  }
  fclose(lfp);

  if(verbose)
    fprintf(stderr, "0 %llu %lf %lf\n", libNum, meanOrig, sigmaOrig);

  // read lengths
  {
    int i;
    mean = meanOrig;
    sigma = sigmaOrig;
    for(i = 0; i < numIterations; i++)
    {
      double minLength;
      double maxLength;
      
      minLength = mean - numSigmas * sigma;
      maxLength = mean + numSigmas * sigma;
      sumX = sumX2 = n = 0;
      
      if(dfn == NULL || (dfp = fopen(dfn, "r")) == NULL)
        usage(argv[0]);
      while(fgets(line, 2047, lfp) != NULL)
      {
        sscanf(line, "%lf", &val);

        if(val < minLength || val > maxLength) continue;
      
        n++;
        sumX += val;
        sumX2 += val * val;
      } // lines in file
      fclose(lfp);

      if(n == 0)
      {
        fprintf(stderr, "Not reestimating library %llu. "
                "It has no clones!\n", libNum);
        break;
      }
      mean = sumX / n;
      sigma = sqrt((sumX2 - n * mean * mean) / n);

      if(verbose)
        fprintf(stderr, "%d %llu %lf %lf\n", i + 1, libNum, mean, sigma);
      
    } // iterations
  } // scope

  fprintf(stdout, "%llu %lf %lf\n", libNum, mean, sigma);
      
  return 0;
}

