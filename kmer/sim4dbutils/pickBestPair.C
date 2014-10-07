#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "bio.h"
#include "sim4.H"

#include <vector>
#include <map>
#include <string>
using namespace std;

#define SEQNAME_MAX  64

class mapResult {
public:
  uint32 seqIdx;
  char   seqName[SEQNAME_MAX];

  uint32 refIdx;
  char   refName[SEQNAME_MAX];
  uint32 refBgn;
  uint32 refEnd;

  bool   forward;
};

class readData {
public:
  readData() {
    cloneIndex  = 999999999;
    isFirstMate = 0;
  };
  readData(uint32 index, uint32 first) {
    cloneIndex  = index;
    isFirstMate = first;
  };

  uint32 cloneIndex  : 31;
  uint32 isFirstMate : 1;
};

map<string,readData>  nameToIndex;
uint32                nameToIndexIndex = 0;


bool
readMR(FILE *in, mapResult &mr) {
  static  char           line[1024];
  static  splitToWords   W;

  //  Skip header.
  if (ftell(in) == 0) {
    fgets(line, 1024, in);
  }

  fgets(line, 1024, in);

  if (feof(in))
    return(false);

  chomp(line);
  W.split(line);

  if (strlen(W[0]) >= SEQNAME_MAX)
    W[0][SEQNAME_MAX-1] = 0;

  if (strlen(W[6]) >= SEQNAME_MAX)
    W[6][SEQNAME_MAX-1] = 0;

  assert(strlen(W[0]) < SEQNAME_MAX);
  assert(strlen(W[6]) < SEQNAME_MAX);

  mr.seqIdx = W(1);
  mr.refIdx = W(7);

  mr.refBgn = W(8);
  mr.refEnd = W(9);

  mr.forward = (W(4) < W(5)) ? true : false;

  strcpy(mr.seqName, W[0]);
  strcpy(mr.refName, W[6]);

  return(true);
}


mapResult &
readMRsim4db(sim4polish *p, mapResult &mr) {

  if (strlen(p->_estDefLine) >= SEQNAME_MAX)
    p->_estDefLine[SEQNAME_MAX-1] = 0;

  if (strlen(p->_genDefLine) >= SEQNAME_MAX)
    p->_genDefLine[SEQNAME_MAX-1] = 0;

  assert(strlen(p->_estDefLine) < SEQNAME_MAX);
  assert(strlen(p->_genDefLine) < SEQNAME_MAX);

  mr.seqIdx = p->_estID;
  mr.refIdx = p->_genID;

  mr.refBgn = p->_exons[0]._genFrom - 1;
  mr.refEnd = p->_exons[0]._genTo;

  mr.forward = (p->_matchOrientation == SIM4_MATCH_FORWARD) ? true : false;

  strcpy(mr.seqName, p->_estDefLine);
  strcpy(mr.refName, p->_genDefLine);

  return(mr);
}


bool
readMRcoords(FILE *in, mapResult &mr) {
  static  char           line[1024];
  static  splitToWords   W;

  //  Skip header.
  if (ftell(in) == 0) {
    fgets(line, 1024, in);
    fgets(line, 1024, in);
    fgets(line, 1024, in);
    fgets(line, 1024, in);
  }

  fgets(line, 1024, in);

  if (feof(in))
    return(false);

  chomp(line);
  W.split(line);

  //  Since we don't have indexes in coords files, we must assign them based on
  //  object names.

  //  But we use "same index" to infer pairing.  This won't work.

  string refNam(W[9]);
  string seqNam(W[10]);

  if (nameToIndex.find(refNam) == nameToIndex.end()) {
    nameToIndex[refNam] = readData(nameToIndexIndex++, false);
  }

  if (nameToIndex.find(seqNam) == nameToIndex.end()) {
    fprintf(stderr, "1 failed to find mate index for read '%s'\n", W[9]);
  }

  uint32  seqIdx = nameToIndex[seqNam].cloneIndex;
  uint32  refIdx = nameToIndex[refNam].cloneIndex;


  if (strlen(W[9]) >= SEQNAME_MAX)
    W[9][SEQNAME_MAX-1] = 0;

  if (strlen(W[10]) >= SEQNAME_MAX)
    W[10][SEQNAME_MAX-1] = 0;

  assert(strlen(W[9]) < SEQNAME_MAX);
  assert(strlen(W[10]) < SEQNAME_MAX);

  mr.seqIdx = seqIdx;
  mr.refIdx = refIdx;

  mr.refBgn = W(0);
  mr.refEnd = W(1);

  mr.forward = (W(2) < W(3)) ? true : false;

  strcpy(mr.seqName, W[10]);
  strcpy(mr.refName, W[9]);

  return(true);
}



bool
readMRcoords(FILE *in, mapResult &mr, bool &is1) {
  static  char           line[1024];
  static  splitToWords   W;

  //  Skip header.
  if (ftell(in) == 0) {
    fgets(line, 1024, in);
    fgets(line, 1024, in);
    fgets(line, 1024, in);
    fgets(line, 1024, in);
  }

  fgets(line, 1024, in);

  if (feof(in))
    return(false);

  chomp(line);
  W.split(line);

  //  Since we don't have indexes in coords files, we must assign them based on
  //  object names.

  //  But we use "same index" to infer pairing.  This won't work.

  string refNam(W[9]);
  string seqNam(W[10]);

  if (nameToIndex.find(refNam) == nameToIndex.end()) {
    nameToIndex[refNam] = readData(nameToIndexIndex++, false);
  }

  if (nameToIndex.find(seqNam) == nameToIndex.end()) {
    fprintf(stderr, "2 failed to find mate index for read '%s'\n", W[10]);
    for (uint32 i=0; i<12; i++)
      fprintf(stderr, "%2d -- '%s'\n", i, W[i]);
    exit(1);
  }

  uint32  seqIdx = nameToIndex[seqNam].cloneIndex;
  uint32  refIdx = nameToIndex[refNam].cloneIndex;

  is1 = nameToIndex[seqNam].isFirstMate;


  if (strlen(W[9]) >= SEQNAME_MAX)
    W[9][SEQNAME_MAX-1] = 0;

  if (strlen(W[10]) >= SEQNAME_MAX)
    W[10][SEQNAME_MAX-1] = 0;

  assert(strlen(W[9]) < SEQNAME_MAX);
  assert(strlen(W[10]) < SEQNAME_MAX);

  mr.seqIdx = seqIdx;
  mr.refIdx = refIdx;

  mr.refBgn = W(0);
  mr.refEnd = W(1);

  mr.forward = (W(2) < W(3)) ? true : false;

  strcpy(mr.seqName, W[10]);
  strcpy(mr.refName, W[9]);

  return(true);
}






int
main(int argc, char **argv) {
  vector<char *>  in1extent, in1sim4db, in1coords, incoords;
  vector<char *>  in2extent, in2sim4db, in2coords;
  vector<char *>  mateMaps;
  char           *out         = NULL;
  char            orient      = 0;
  uint32          distMin     = 0;
  uint32          distMax     = uint32MAX;

  double          minIdent    = 0;
  double          minLength   = 0;
  double          minCoverage = 0;

  bool            allowDups   = false;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-1extent") == 0)
      while ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        in1extent.push_back(argv[++arg]);
    else if (strcmp(argv[arg], "-2extent") == 0)
      while ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        in2extent.push_back(argv[++arg]);

    else if (strcmp(argv[arg], "-1sim4db") == 0)
      while ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        in1sim4db.push_back(argv[++arg]);
    else if (strcmp(argv[arg], "-2sim4db") == 0)
      while ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        in2sim4db.push_back(argv[++arg]);

    else if (strcmp(argv[arg], "-1coords") == 0)
      while ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        in1coords.push_back(argv[++arg]);
    else if (strcmp(argv[arg], "-2coords") == 0)
      while ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        in2coords.push_back(argv[++arg]);

    else if (strcmp(argv[arg], "-coords") == 0)
      while ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        incoords.push_back(argv[++arg]);

    else if (strcmp(argv[arg], "-matemap") == 0)
      while ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        mateMaps.push_back(argv[++arg]);

    else if (strcmp(argv[arg], "-insert") == 0) {
      orient  = argv[++arg][0];
      distMin = atoi(argv[++arg]);
      distMax = atoi(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-minident") == 0)
      minIdent = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-minlength") == 0)
      minLength = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-mincoverage") == 0)
      minCoverage = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-allowduplicates") == 0)
      allowDups = true;

    else if (strcmp(argv[arg], "-o") == 0)
      out  = argv[++arg];

    else
      err++;

    arg++;
  }
  if (out == NULL) {
    fprintf(stderr, "usage: %s -1 in1.extent -2 in2.extent -o prefix\n", argv[0]);
    exit(1);
  }

  vector<mapResult>   mr1;
  vector<mapResult>   mr2;

  mapResult           mr;

  //  Load mate map if needed

  if (mateMaps.size() > 0) {
    for (uint32 mm=0; mm<mateMaps.size(); mm++) {
      uint32  numLoaded = 0;

      fprintf(stderr, "Loading mate pairings from '%s'.\n", mateMaps[mm]);

      errno = 0;
      FILE   *IN = fopen(mateMaps[mm], "r");
      if (errno)
        fprintf(stderr, "Failed to open mate map '%s': %s\n", mateMaps[mm], strerror(errno)), exit(1);

      char    LL[10240];
      fgets(LL, 10240, IN);

      while (!feof(IN)) {
        chomp(LL);
        splitToWords  W(LL);

        nameToIndex[string(W[0])] = readData(nameToIndexIndex,   true);
        nameToIndex[string(W[1])] = readData(nameToIndexIndex++, false);

        numLoaded++;

        fgets(LL, 10240, IN);
      }

      fprintf(stderr, "Loaded %u mate pairings from '%s', total %u.\n", numLoaded, mateMaps[mm], nameToIndexIndex);
    }
  }

  //  Load alignments

  for (uint32 ii=0; ii<in1extent.size(); ii++) {
    fprintf(stderr, "Loading alignments from '%s'\n", in1extent[ii]);
    FILE  *IN = fopen(in1extent[ii], "r");
    while (readMR(IN, mr) == true)
      mr1.push_back(mr);
    fclose(IN);
  }

  for (uint32 ii=0; ii<in1sim4db.size(); ii++) {
    fprintf(stderr, "Loading alignments from '%s'\n", in1sim4db[ii]);
    sim4polishReader  *IN = new sim4polishReader(in1sim4db[ii]);
    sim4polish        *p   = NULL;
    while (IN->nextAlignment(p)) {
      mr1.push_back(readMRsim4db(p, mr));
    }
    delete IN;
  }

  for (uint32 ii=0; ii<in1coords.size(); ii++) {
    fprintf(stderr, "Loading alignments from '%s'\n", in1coords[ii]);
    FILE  *IN = fopen(in1coords[ii], "r");
    while (readMRcoords(IN, mr) == true)
      mr1.push_back(mr);
    fclose(IN);
  }



  for (uint32 ii=0; ii<in2extent.size(); ii++) {
    fprintf(stderr, "Loading alignments from '%s'\n", in2extent[ii]);
    FILE  *IN = fopen(in2extent[ii], "r");
    while (readMR(IN, mr) == true)
      mr2.push_back(mr);
    fclose(IN);
  }

  for (uint32 ii=0; ii<in2sim4db.size(); ii++) {
    fprintf(stderr, "Loading alignments from '%s'\n", in2sim4db[ii]);
    sim4polishReader  *IN = new sim4polishReader(in2sim4db[ii]);
    sim4polish        *p   = NULL;
    while (IN->nextAlignment(p)) {
      mr2.push_back(readMRsim4db(p, mr));
    }
    delete IN;
  }

  for (uint32 ii=0; ii<in2coords.size(); ii++) {
    fprintf(stderr, "Loading alignments from '%s'\n", in2coords[ii]);
    FILE  *IN = fopen(in2coords[ii], "r");
    while (readMRcoords(IN, mr) == true)
      mr1.push_back(mr);
    fclose(IN);
  }



  for (uint32 ii=0; ii<incoords.size(); ii++) {
    fprintf(stderr, "Loading alignments from '%s'\n", incoords[ii]);
    FILE  *IN = fopen(incoords[ii], "r");
    bool   is1;

    while (readMRcoords(IN, mr, is1) == true)
      if (is1)
        mr1.push_back(mr);
      else
        mr2.push_back(mr);

    fclose(IN);
  }

  fprintf(stderr, "Loaded %lu '1' alignments.\n", mr1.size());
  fprintf(stderr, "Loaded %lu '2' alignments.\n", mr2.size());

  char   name[10240];

  sprintf(name, "%s.pairLog", out);
  FILE *LOG = fopen(name, "w");

  sprintf(name, "%s.duplicates", out);
  FILE *DUP = fopen(name, "w");

  sprintf(name, "%s.stats", out);
  FILE *STA = fopen(name, "w");

  uint32   mr1bgn = 0;
  uint32   mr1end = 0;
  uint32   mr1END = mr1.size();

  uint32   mr2bgn = 0;
  uint32   mr2end = 0;
  uint32   mr2END = mr2.size();

  map<char,uint32>   totalPairs;
  map<char,uint32>   sizedPairs;

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

    if ((mr1end - mr1bgn > 1) &&
        (mr2end - mr2bgn > 1)) {
      fprintf(DUP, "%s\t%u\t%s\t%u\n",
              mr1[mr1bgn].seqName, mr1end - mr1bgn,
              mr2[mr2bgn].seqName, mr2end - mr2bgn);
      if (allowDups == false) {
        mr1bgn = mr1end;
        mr2bgn = mr2end;
      }
    }

    //  Now find all possible pairs.

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

          totalPairs[ori]++;

          if ((orient == 0) ||
              ((ori == orient) && (distMin <= df) && (df <= distMax))) {
            sizedPairs[ori]++;
            fprintf(LOG, "%c "uint32FMT" "uint32FMT" %s ("uint32FMT","uint32FMT") "uint32FMT" %s ("uint32FMT","uint32FMT") "uint32FMT" %s\n",
                    ori,
                    df,
                    mr1[i1].seqIdx, mr1[i1].seqName, mr1[i1].refBgn, mr1[i1].refEnd,
                    mr2[i2].seqIdx, mr2[i2].seqName, mr2[i2].refBgn, mr2[i2].refEnd,
                    mr1[i1].refIdx, mr1[i1].refName);
          }

        } else {
          if ((mr2[i2].forward == true)  && (mr1[i1].forward == true))
            ori = 'N';
          if ((mr2[i2].forward == true)  && (mr1[i1].forward == false))
            ori = 'I';
          if ((mr2[i2].forward == false) && (mr1[i1].forward == true))
            ori = 'O';
          if ((mr2[i2].forward == false) && (mr1[i1].forward == false))
            ori = 'A';

          totalPairs[ori]++;

          if ((orient == 0) ||
              ((ori == orient) && (distMin <= dr) && (dr <= distMax))) {
            sizedPairs[ori]++;
            fprintf(LOG, "%c "uint32FMT" "uint32FMT" %s ("uint32FMT","uint32FMT") "uint32FMT" %s ("uint32FMT","uint32FMT") "uint32FMT" %s\n",
                    ori,
                    dr,
                    mr2[i2].seqIdx, mr2[i2].seqName, mr2[i2].refBgn, mr2[i2].refEnd,
                    mr1[i1].seqIdx, mr1[i1].seqName, mr1[i1].refBgn, mr1[i1].refEnd,
                    mr2[i2].refIdx, mr2[i2].refName);
          }
        }
      }
    }

    mr1bgn = mr1end;
    mr2bgn = mr2end;
  }

  fprintf(STA, "alignments:      "uint32FMT" "uint32FMT"\n", mr1END, mr2END);
  fprintf(STA, "totalPairs[%c]:  %u\n", 'N', totalPairs['N']);
  fprintf(STA, "totalPairs[%c]:  %u\n", 'I', totalPairs['I']);
  fprintf(STA, "totalPairs[%c]:  %u\n", 'O', totalPairs['O']);
  fprintf(STA, "totalPairs[%c]:  %u\n", 'A', totalPairs['A']);
  fprintf(STA, "sizedPairs[%c]:  %u\n", 'N', sizedPairs['N']);
  fprintf(STA, "sizedPairs[%c]:  %u\n", 'I', sizedPairs['I']);
  fprintf(STA, "sizedPairs[%c]:  %u\n", 'O', sizedPairs['O']);
  fprintf(STA, "sizedPairs[%c]:  %u\n", 'A', sizedPairs['A']);

  fclose(LOG);
  fclose(DUP);
  fclose(STA);

  exit(0);
}
