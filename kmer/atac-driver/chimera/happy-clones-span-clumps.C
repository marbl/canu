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
    clumpsMax = 64;
    clumpID   = new u32bit [clumpsMax];
    clumps    = new intervalList * [clumpsMax];
    clumpmin  = new u32bit [clumpsMax];
    clumpmax  = new u32bit [clumpsMax];

    clumpconfirm = new u32bit [clumpsMax * clumpsMax];

    for (u32bit i=0; i<clumpsMax * clumpsMax; i++)
      clumpconfirm[i] = 0;

    intervalsLen = 0;
    intervalsMax = 0;
    intervals    = 0L;
  };

  ~atacClumpCoordTreeScaffold() {
    for (u32bit i=0; i<clumpsLen; i++)
      delete clumps[i];
    delete [] clumpID;
    delete [] clumps;
    delete [] clumpmin;
    delete [] clumpmax;
    delete [] clumpconfirm;
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
        if (clumpmin[i] > begin)
          clumpmin[i] = begin;
        if (clumpmax[i] < begin + length)
          clumpmax[i] = begin + length;
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
    clumpmin[clumpsLen] = begin;
    clumpmax[clumpsLen] = begin + length;
    clumpsLen++;
  };


  u32bit  getClumpID(u32bit begin, u32bit end) {
    u32bit  clumpid = 0;
    u32bit  numhits = 0;

    //  We can make this much quicker if we remember the extent of
    //  each interval list.
    //
    //  We want to allow partial matches, so check that the end is
    //  above the min, and the begin is before the max.
    //
    //         b-------e               b-----e
    //                 -------clump------
    //
    for (u32bit i=0; i<clumpsLen; i++) {
      if ((clumpmin[i] <= end) && (begin <= clumpmax[i])) {
        if (clumps[i]->overlapping(begin, end, intervals, intervalsLen, intervalsMax) > 0) {
          clumpid = clumpID[i];
          numhits++;
        }
      } else {
        //  If you really want to check....
        //if (clumps[i]->overlapping(begin, end, intervals, intervalsLen, intervalsMax) > 0)
        //  fprintf(stderr, "WARNING: Found overlapping clump outside extent!\n");
      }
    }

    if (numhits == 0)
      return(0);
    if (numhits == 1)
      return(clumpid);

    //fprintf(stderr, "FOUND MORE THAN ONE CLUMP MATCHING!\n");
    return(~u32bitZERO);
  };


  void    sortClumps(void) {
    u32bit         ciid;
    intervalList  *cptr;
    u32bit         cmin;
    u32bit         cmax;

    u32bit         i = 0;
    u32bit         j = 0;

    //  an insertion sort

    for (i=clumpsLen; i--; ) {
      ciid = clumpID[i];
      cptr = clumps[i];
      cmin = clumpmin[i];
      cmax = clumpmax[i];

      for (j=i+1; (j < clumpsLen) && (cmin > clumpmin[j]); j++) {
        clumpID[j-1]  = clumpID[j];
        clumps[j-1]   = clumps[j];
        clumpmin[j-1] = clumpmin[j];
        clumpmax[j-1] = clumpmax[j];
      }

      clumpID[j-1]  = ciid;
      clumps[j-1]   = cptr;
      clumpmin[j-1] = cmin;
      clumpmax[j-1] = cmax;
    }
  };


  void    confirm(u32bit ca, u32bit cb) {
    u32bit  caidx = 0;
    u32bit  cbidx = 0;
    for (u32bit i=0; i<clumpsLen; i++) {
      if (ca == clumpID[i])
        caidx = i;
      if (cb == clumpID[i])
        cbidx = i;
    }
    clumpconfirm[caidx * clumpsMax + cbidx]++;
  };


  u32bit           clumpsLen;
  u32bit           clumpsMax;
  u32bit          *clumpID;
  intervalList   **clumps;
  u32bit          *clumpmin;
  u32bit          *clumpmax;

  u32bit          *clumpconfirm;

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


  void    removeSingleClumpScaffolds(void) {
    u32bit deleted = 0;
    u32bit remain = 0;


    for (u32bit i=0; i<scaffoldsMax; i++) {
      if ((scaffolds[i]) && (scaffolds[i]->clumpsLen < 2)) {
        delete scaffolds[i];
        scaffolds[i] = 0L;
        deleted++;
      }
      if (scaffolds[i]) {
        scaffolds[i]->sortClumps();
        remain++;
      }
    }
    fprintf(stderr, "Deleted "u32bitFMT" scaffolds with less than 2 clumps.\n", deleted);
    fprintf(stderr, "Remain  "u32bitFMT" scaffolds with more than 2 clumps.\n", remain);
  };


  void    showMultipleClumpScaffolds(void) {

    for (u32bit i=0; i<scaffoldsMax; i++) {
      if ((scaffolds[i]) && (scaffolds[i]->clumpsLen >= 2)) {

        fprintf(stdout, "\n");

        for (u32bit j=0; j<scaffolds[i]->clumpsLen; j++) {
          bool  overlap = false;

          if ((j+1 < scaffolds[i]->clumpsLen) &&
              (scaffolds[i]->clumpmax[j] > scaffolds[i]->clumpmin[j+1]))
            overlap = true;

          fprintf(stdout, "scaffold "u32bitFMT" clump "u32bitFMT" begin "u32bitFMT" end "u32bitFMT"\n",
                  i,
                  scaffolds[i]->clumpID[j],
                  scaffolds[i]->clumpmin[j],
                  scaffolds[i]->clumpmax[j]);

          if (overlap)
            fprintf(stdout, "scaffold "u32bitFMT" clump "u32bitFMT" and clump "u32bitFMT" OVERLAP\n",
                    i,
                    scaffolds[i]->clumpID[j],
                    scaffolds[i]->clumpID[j+1]);

          for (u32bit b=0; b<scaffolds[i]->clumpsLen; b++) {
            u32bit cc = j * scaffolds[i]->clumpsMax + b;
            if (scaffolds[i]->clumpconfirm[cc]) {
              fprintf(stdout, "scaffold "u32bitFMT" clump "u32bitFMT" and "u32bitFMT" confirmed by "u32bitFMT" clones.\n",
                      i,
                      scaffolds[i]->clumpID[j],
                      scaffolds[i]->clumpID[b],
                      scaffolds[i]->clumpconfirm[cc]);
              
            }
          }
        }
      }
    }
  };


  u32bit  getClumpID(u32bit scaffoldid, u32bit begin, u32bit end) {
    if (scaffolds[scaffoldid])
      return(scaffolds[scaffoldid]->getClumpID(begin, end));
    return(0);
  };

  void    confirmClump(u32bit scaffoldid, u32bit ca, u32bit cb) {
    if (scaffolds[scaffoldid]) {
      scaffolds[scaffoldid]->confirm(ca, cb);
    }
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

  if (clumpFile == 0L)
    fprintf(stderr, "No -clumps supplied!\n"), err++;
  if (happyFile == 0L)
    fprintf(stderr, "No -happy clones supplied!\n"), err++;
  if (err)
    fprintf(stderr, "usage: %s ...\n", argv[0]), exit(1);

  atacClumpCoordTree *ct = buildCoordTree(clumpFile);

  ct->removeSingleClumpScaffolds();
  ct->showMultipleClumpScaffolds();

  ////////////////////////////////////////
  //
  //  ugly hack -- read in the map from HUREF6A UID to scaffold
  //
  map<u64bit,u32bit>  UIDtoIID;

  {
    char *uidmapName = "/project/huref6/assembly/fasta/HUREF6A.info";

    errno = 0;
    FILE *F = fopen(uidmapName, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", uidmapName, strerror(errno)), exit(1);

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

  speedCounter         S("%9.0f clones (%6.1f clones/sec)\r", 1, 4096, true);

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
      fprintf(stdout, "%s spans clump in scaffold %s,"u32bitFMT"\n", ina, A[7], scfa);
      cla = 0;
    }
    if (clb == ~u32bitZERO) {
      fprintf(stdout, "%s spans clump in scaffold %s,"u32bitFMT" \n", inb, B[7], scfb);
      clb = 0;
    }

    if ((cla != 0) &&
        (clb != 0) &&
        (cla != clb)) {
      ct->confirmClump(scfa,
                       (cla < clb) ? cla : clb,
                       (cla < clb) ? clb : cla);

      fprintf(stdout, "scaffold %s,"u32bitFMT" clump "u32bitFMT" "u32bitFMT" confirmed by %s\n",
              A[7], scfa,
              (cla < clb) ? cla : clb,
              (cla < clb) ? clb : cla,
              A[2]);
    }

    S.tick();

    fgets(ina, 1024, inf);  chomp(ina);
    fgets(inb, 1024, inf);  chomp(inb);
  }

  fclose(inf);

  S.finish();

  ct->showMultipleClumpScaffolds();

  delete ct;
}
