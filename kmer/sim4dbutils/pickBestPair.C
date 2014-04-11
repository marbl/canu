#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "bio.h"
#include "sim4.H"

#include <vector>
using namespace std;

#define NAME_MAX  32

struct mapResult {
  uint32 seqIdx;
  char   seqName[NAME_MAX];

  uint32 refIdx;
  char   refName[NAME_MAX];
  uint32 refBgn;
  uint32 refEnd;

  bool   forward;
};




bool
readMR(FILE *in, mapResult &mr) {
  static  char           line[1024];
  static  splitToWords   W;

  fgets(line, 1024, in);

  if (feof(in))
    return(false);

  W.split(line);

  if (strcmp(W[0], "cDNAid") == 0)
    return(readMR(in, mr));

  if (strlen(W[0]) >= NAME_MAX)
    W[0][NAME_MAX-1] = 0;

  if (strlen(W[6]) >= NAME_MAX)
    W[6][NAME_MAX-1] = 0;

  assert(strlen(W[0]) < NAME_MAX);
  assert(strlen(W[6]) < NAME_MAX);

  mr.seqIdx = W(1);
  mr.refIdx = W(7);

  mr.refBgn = W(8);
  mr.refEnd = W(9);
      
  mr.forward = (W(4) < W(5)) ? true : false;

  strcpy(mr.seqName, W[0]);
  strcpy(mr.refName, W[6]);

  return(true);
}


int
main(int argc, char **argv) {
  char  *in1 = NULL;
  char  *in2 = NULL;
  char  *out = NULL;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-1") == 0)
      in1  = argv[++arg];

    else if (strcmp(argv[arg], "-2") == 0)
      in2  = argv[++arg];

    else if (strcmp(argv[arg], "-o") == 0)
      out  = argv[++arg];

    else
      err++;

    arg++;
  }
  if ((in1 == NULL) || (in2 == NULL) || (out == NULL)) {
    fprintf(stderr, "usage: %s -1 in1.extent -2 in2.extent -o prefix\n", argv[0]);
    exit(1);
  }

  vector<mapResult>   mr1;
  vector<mapResult>   mr2;

  mapResult           mr;

  fprintf(stderr, "Loading alignments from '%s'\n", in1);
  FILE  *IN1 = fopen(in1, "r");
  while (readMR(IN1, mr) == true)
    mr1.push_back(mr);
  fclose(IN1);

  fprintf(stderr, "Loading alignments from '%s'\n", in2);
  FILE  *IN2 = fopen(in2, "r");
  while (readMR(IN2, mr) == true)
    mr2.push_back(mr);
  fclose(IN2);

  FILE  *LOG = NULL;
  FILE  *STA = NULL;

  char   name[10240];

  sprintf(name, "%s.pairLog", out);
  LOG = fopen(name, "w");

  sprintf(name, "%s.stats", out);
  STA = fopen(name, "w");

  uint32   mr1bgn = 0;
  uint32   mr1end = 0;
  uint32   mr1END = mr1.size();

  uint32   mr2bgn = 0;
  uint32   mr2end = 0;
  uint32   mr2END = mr2.size();

  while ((mr1bgn < mr1END) && (mr2bgn < mr2END)) {

    if ((mr1[mr1bgn].seqIdx < mr2[mr2bgn].seqIdx) && (mr1bgn < mr1END))
      mr1bgn++;

    if ((mr2[mr2bgn].seqIdx < mr1[mr1bgn].seqIdx) && (mr2bgn < mr2END))
      mr2bgn++;

    if (mr1[mr1bgn].seqIdx != mr2[mr2bgn].seqIdx)
      //  SequenceA 1 3 5 7 8
      //  SequenceB  2 4 6  8
      //  1st pass, A increases to 3, B increases to 4
      //  2nd pass, A increases to 5, B increases to 6
      //  3rd pass, A increases to 7, B increases to 8
      //  4th pass, A increases to 8, B doesn't change.
      continue;

    assert(mr1[mr1bgn].seqIdx == mr2[mr2bgn].seqIdx);

    mr1end = mr1bgn + 1;
    mr2end = mr2bgn + 1;

    while (mr1[mr1bgn].seqIdx == mr1[mr1end].seqIdx)
      mr1end++;

    while (mr2[mr2bgn].seqIdx == mr2[mr2end].seqIdx)
      mr2end++;

    //  Group of reads from mr1bgn-mr1end and mr2bgn-mr2end need to be compared.

    for (uint32 i1=mr1bgn; i1<mr1end; i1++) {
      for (uint32 i2=mr2bgn; i2<mr2end; i2++) {
        if (mr1[i1].refIdx != mr2[i2].refIdx)
          continue;

        //validParis++;

        uint32  df  = 0;
        uint32  dr  = 0;
        char    ori = 'X';

        if (mr1[i1].refBgn < mr2[i2].refEnd)
          df = mr2[i2].refEnd - mr1[i1].refBgn;

        if (mr2[i2].refBgn < mr1[i1].refEnd)
          dr = mr1[i1].refEnd - mr2[i2].refBgn;

        assert(df + dr > 0);

        if (df > dr) {
          if ((mr1[i1].forward == true)  && (mr2[i2].forward == true))
            ori = 'N';
          if ((mr1[i1].forward == true)  && (mr2[i2].forward == false))
            ori = 'I';
          if ((mr1[i1].forward == false) && (mr2[i2].forward == true))
            ori = 'O';
          if ((mr1[i1].forward == false) && (mr2[i2].forward == false))
            ori = 'A';

          fprintf(LOG, "%c "uint32FMT" "uint32FMT" %s ("uint32FMT","uint32FMT") "uint32FMT" %s ("uint32FMT","uint32FMT") "uint32FMT" %s\n",
                  ori,
                  df,
                  mr1[i1].seqIdx, mr1[i1].seqName, mr1[i1].refBgn, mr1[i1].refEnd,
                  mr2[i2].seqIdx, mr2[i2].seqName, mr2[i2].refBgn, mr2[i2].refEnd,
                  mr1[i1].refIdx, mr1[i1].refName);

        } else {
          if ((mr2[i2].forward == true)  && (mr1[i1].forward == true))
            ori = 'N';
          if ((mr2[i2].forward == true)  && (mr1[i1].forward == false))
            ori = 'I';
          if ((mr2[i2].forward == false) && (mr1[i1].forward == true))
            ori = 'O';
          if ((mr2[i2].forward == false) && (mr1[i1].forward == false))
            ori = 'A';

          fprintf(LOG, "%c "uint32FMT" "uint32FMT" %s ("uint32FMT","uint32FMT") "uint32FMT" %s ("uint32FMT","uint32FMT") "uint32FMT" %s\n",
                  ori,
                  dr,
                  mr2[i2].seqIdx, mr2[i2].seqName, mr2[i2].refBgn, mr2[i2].refEnd,
                  mr1[i1].seqIdx, mr1[i1].seqName, mr1[i1].refBgn, mr1[i1].refEnd,
                  mr2[i2].refIdx, mr2[i2].refName);
        }
      }
    }

    mr1bgn = mr1end;
    mr2bgn = mr2end;
  }

  fprintf(STA, "alignments:  "uint32FMT" "uint32FMT"\n", mr1END, mr2END);

  fclose(LOG);
  fclose(STA);

  exit(0);
}
