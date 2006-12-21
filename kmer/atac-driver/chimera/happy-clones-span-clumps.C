#include <stdio.h>
#include <stdlib.h>

using namespace std;
#include <map>

#include "util++.H"

//  Reads a clump-annotated atac file, builds a search tree of all the
//  matches in those clumps.  Then reads a list of happy clones mapped
//  to the sequence, figures out what clump each read in the clone is
//  in, and reports whenever the clone spans a clump.


//  Contains a list of intervalLists, one for each clump.  The
//  intervalList stores the positions of the matches in this clump.
//
class atacClumpCoordTreeScaffold {
public:
  atacClumpCoordTreeScaffold() {
    clumpsLen = 0;
    clumpsMax = 8;
    clumpID   = new u32bit [clumpsMax];
    clumps    = new intervalList * [clumpsMax];

    intervalsLen = 0;
    intervalsMax = 0;
    intervals    = 0L;
  };

  ~atacClumpCoordTreeScaffold() {
    for (u32bit i=0; i<clumpsLen; i++)
      delete clumps[i];
    delete [] clumpID;
    delete [] clumps;
    delete [] intervals;
  };

  //  Add a match to some clump.
  //
  void    addMatch(s32bit clumpid, u32bit begin, u32bit length) {

    //  Not in a clump, get the heck outta here!
    //
    if (clumpid < 0)
      return;

    //  Linear search through the clumps to find the correct id, we
    //  don't expect to have many clumps per scaffold.
    //
    for (u32bit i=0; i<clumpsLen; i++) {
      if (clumpID[i] == (u32bit)clumpid) {
        clumps[i]->add(begin, length);
        return;
      }
    }

    if (clumpsLen == clumpsMax) {
      fprintf(stderr, "ERROR: increase clumpsMax!\n");
      exit(1);
    }

    //  Didn't add to an existing clump, so must be a new clump.
    //
    clumpID[clumpsLen] = clumpid;
    clumps[clumpsLen]  = new intervalList;
    clumps[clumpsLen]->add(begin, length);
    clumpsLen++;
  };


  u32bit  getClumpID(u32bit begin, u32bit end) {
    u32bit  clumpid = 0;
    u32bit  numhits = 0;

    //  We can make this much quicker if we remember the extent of each interval list.

    for (u32bit i=0; i<clumpsLen; i++) {
      if (clumps[i]->overlapping(begin, end, intervals, intervalsLen, intervalsMax) > 0) {
        clumpid = clumpID[i];
        numhits++;
      }
    }

    if (numhits == 0)
      return(0);
    if (numhits == 1)
      return(clumpid);

    //fprintf(stderr, "FOUND MORE THAN ONE CLUMP MATCHING!\n");
    return(~u32bitZERO);
  };


  u32bit           clumpsLen;
  u32bit           clumpsMax;
  u32bit          *clumpID;
  intervalList   **clumps;

  u32bit           intervalsLen;
  u32bit           intervalsMax;
  u32bit          *intervals;
};


class atacClumpCoordTree {
public:
  atacClumpCoordTree() {
    scaffoldsMax = 262144;
    scaffolds    = new atacClumpCoordTreeScaffold * [scaffoldsMax];
    for (u32bit i=0; i<scaffoldsMax; i++)
      scaffolds[i] = 0L;
  };
  ~atacClumpCoordTree() {
    for (u32bit i=0; i<scaffoldsMax; i++)
      if (scaffolds[i])
        delete scaffolds[i];
    delete [] scaffolds;
  };

  void    addMatch(u32bit scaffoldid, s32bit clumpid, u32bit begin, u32bit length) {
    if (scaffoldid >= scaffoldsMax) {
      fprintf(stderr, "ERROR: increase scaffoldsMax "u32bitFMT"\n", scaffoldid);
      exit(1);
    }

    if (scaffolds[scaffoldid] == 0L)
      scaffolds[scaffoldid] = new atacClumpCoordTreeScaffold;

    scaffolds[scaffoldid]->addMatch(clumpid, begin, length);
  };
  
  u32bit  getClumpID(u32bit scaffoldid, u32bit begin, u32bit end) {
    if (scaffolds[scaffoldid])
      return(scaffolds[scaffoldid]->getClumpID(begin, end));
    return(0);
  };
  
  u32bit                       scaffoldsMax;
  atacClumpCoordTreeScaffold **scaffolds;
};





atacClumpCoordTree*
buildCoordTree(char *clumpFile) {
  atacClumpCoordTree  *ct = new atacClumpCoordTree;

  FILE                *inf;
  char                 inl[1024];

  //  We can't use the built-in atac reader, because it strips out
  //  clump information.  Bummer.

  errno = 0;
  inf = fopen(clumpFile, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", clumpFile, strerror(errno)), exit(1);

  fgets(inl, 1024, inf);

  if (feof(inf))
    return(0L);

  while (!feof(inf)) {
    if ((inl[0] == 'M') &&
        (inl[2] == 'u')) {
      splitToWords  S(inl);

      //fprintf(stderr, "%s", inl);

      if (S[12][0] != '#')
        fprintf(stderr, "no clump for '%s'\n", inl);

      if (S[13][0] != '-') {
        char *scfid = S[8];
        while (*scfid != ':')
          scfid++;

        ct->addMatch(atoi(scfid + 1),
                     atoi(S[13]),
                     atoi(S[9]),
                     atoi(S[10]));
      }
    }

    fgets(inl, 1024, inf);
  }

  fclose(inf);

  return(ct);
}








int
main(int argc, char **argv) {
  char  *clumpFile = 0L;
  char  *happyFile = 0L;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-clumps") == 0) {
      clumpFile = argv[++arg];
    } else if (strcmp(argv[arg], "-happy") == 0) {
      happyFile = argv[++arg];
    } else {
      err++;
    }

    arg++;
  }

  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    exit(1);
  }

  atacClumpCoordTree *ct = buildCoordTree(clumpFile);

  ////////////////////////////////////////
  //
  //  ugly hack -- read in the map from HUREF6A UID to scaffold
  //
  map<u64bit,u32bit>  UIDtoIID;

  {
    FILE *F = fopen("/home/work/chimera/HUREF6A.info", "r");
    char  L[1024];

    fgets(L, 1024, F);
    while (!feof(F)) {
      if (L[0] == 'G') {
        splitToWords  S(L);

        UIDtoIID[strtou64bit(S[13]+1, 0L)] = strtou32bit(S[10], 0L);
      }
      fgets(L, 1024, F);
    }
  }
  //
  ////////////////////////////////////////

  FILE                *inf;
  char                 ina[1024];
  char                 inb[1024];

  errno = 0;
  inf = fopen(happyFile, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", happyFile, strerror(errno)), exit(1);

  fgets(ina, 1024, inf);  chomp(ina);
  fgets(inb, 1024, inf);  chomp(inb);
  while (!feof(inf)) {
    splitToWords  A(ina);
    splitToWords  B(inb);

    //  Some sanity checking.
    if (strcmp(A[2], B[2]) != 0) {
      fprintf(stderr, "ERROR:  Different clone!\n%s\n%s\n", ina, inb);
    }


    u32bit  scfa = UIDtoIID[strtou64bit(A[7], 0L)];
    u32bit  scfb = UIDtoIID[strtou64bit(B[7], 0L)];

    u32bit  cla = ct->getClumpID(scfa,
                                 atoi(A[8]),
                                 atoi(A[9]));
    u32bit  clb = ct->getClumpID(scfb,
                                 atoi(B[8]),
                                 atoi(B[9]));

    if (cla == ~u32bitZERO) {
      fprintf(stdout, "%s spans clump in scaffold "u32bitFMT"\n", ina, scfa);
      cla = 0;
    }
    if (clb == ~u32bitZERO) {
      fprintf(stdout, "%s spans clump in scaffold "u32bitFMT" \n", inb, scfb);
      clb = 0;
    }

    if ((cla != 0) &&
        (clb != 0) &&
        (cla != clb)) {
      fprintf(stdout, "scaffold "u32bitFMT" clump "u32bitFMT" "u32bitFMT" confirmed by %s\n",
              scfa,
              (cla < clb) ? cla : clb,
              (cla < clb) ? clb : cla,
              A[2]);
    } else {
    }

    fgets(ina, 1024, inf);  chomp(ina);
    fgets(inb, 1024, inf);  chomp(inb);
  }

  fclose(inf);

  delete ct;
}
