#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
//#include <values.h>
#include <math.h>

#define DEFAULT_NUM_SIGMAS       4


void usage(char * progname)
{
  fprintf(stderr, "Usage: %s [-h] -f lengthsFile  -l libFile  -n libNum  -s sigmas -c count\n"
          "  -h               print help\n"
          "  -f lengthsFile   name of file containing clone length data\n"
          "                     default is stdin\n"
          "  -l libFile       file with mean & stdev of libs\n"
          "  -n libNum        UID of 'other' library in libFile to check against\n"
          "  -s sigmas        require lengths to be within sigmas from mean\n"
          "                     default is %u\n"
          "  -c count         require count instances of clone to identify as mistracked\n",

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
  unsigned int minCount = 0;
  double numSigmas = DEFAULT_NUM_SIGMAS;
  double val;
  double mean;
  double sigma;

  // scope for scope's sake
  {
    int ch;
    while((ch = getopt(argc, argv, "hf:l:n:s:c:")) != EOF)
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
        case 'c':
          minCount = atoi(optarg);
          break;
        default:
          usage(argv[0]);
          break;
      }
    }
  }

  if(minCount <= 1 || numSigmas <= 0)
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
    int thisUIDCount = 1;

    fprintf(stderr, "Identifying clones between %lf and %lf bp with %d or more instances\n", minLength, maxLength, minCount);
      
    if(dfn == NULL || (dfp = fopen(dfn, "r")) == NULL)
      usage(argv[0]);
    while(fgets(line, 2047, lfp) != NULL)
    {
      sscanf(line, "%lf %llu %llu", &val, &minUID, &maxUID);

      if(val < minLength || val > maxLength) continue;

      if(minUID == lastMinUID)
      {
        thisUIDCount++;
      }
      else
      {
        if(thisUIDCount >= minCount)
        {
          fprintf(stdout, "%lf %llu %llu %d\n",
                  val, minUID, maxUID, thisUIDCount);
        }
        thisUIDCount = 1;
      }
      lastMinUID = minUID;
    } // lines in file
    fclose(lfp);
  } // scope

  return 0;
}

