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
    clumpID   = new uint32 [clumpsMax];
    clumps    = new intervalList * [clumpsMax];
    clumpmin  = new uint32 [clumpsMax];
    clumpmax  = new uint32 [clumpsMax];

    clumpconfirm = new uint32 [clumpsMax * clumpsMax];

    for (uint32 i=0; i<clumpsMax * clumpsMax; i++)
      clumpconfirm[i] = 0;

    intervalsLen = 0;
    intervalsMax = 0;
    intervals    = 0L;
  };

  ~atacClumpCoordTreeScaffold() {
    for (uint32 i=0; i<clumpsLen; i++)
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
  void    addMatch(int32 clumpid, uint32 begin, uint32 length) {

    //  Not in a clump, get the heck outta here!
    //
    if (clumpid < 0)
      return;

    //  Linear search through the clumps to find the correct id, we
    //  don't expect to have many clumps per scaffold.
    //
    for (uint32 i=0; i<clumpsLen; i++) {
      if (clumpID[i] == (uint32)clumpid) {
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


  uint32  getClumpID(uint32 begin, uint32 end) {
    uint32  clumpid = 0;
    uint32  numhits = 0;

    //  We can make this much quicker if we remember the extent of
    //  each interval list.
    //
    //  We want to allow partial matches, so check that the end is
    //  above the min, and the begin is before the max.
    //
    //         b-------e               b-----e
    //                 -------clump------
    //
    for (uint32 i=0; i<clumpsLen; i++) {
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
    return(~uint32ZERO);
  };


  void    sortClumps(void) {
    uint32         ciid;
    intervalList  *cptr;
    uint32         cmin;
    uint32         cmax;

    uint32         i = 0;
    uint32         j = 0;

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


  void    confirm(uint32 ca, uint32 cb) {
    uint32  caidx = 0;
    uint32  cbidx = 0;
    for (uint32 i=0; i<clumpsLen; i++) {
      if (ca == clumpID[i])
        caidx = i;
      if (cb == clumpID[i])
        cbidx = i;
    }
    clumpconfirm[caidx * clumpsMax + cbidx]++;
  };


  uint32           clumpsLen;
  uint32           clumpsMax;
  uint32          *clumpID;
  intervalList   **clumps;
  uint32          *clumpmin;
  uint32          *clumpmax;

  uint32          *clumpconfirm;

  uint32           intervalsLen;
  uint32           intervalsMax;
  uint32          *intervals;
};


class atacClumpCoordTree {
public:
  atacClumpCoordTree() {
    scaffoldsMax = 262144;
    scaffolds    = new atacClumpCoordTreeScaffold * [scaffoldsMax];
    for (uint32 i=0; i<scaffoldsMax; i++)
      scaffolds[i] = 0L;
  };
  ~atacClumpCoordTree() {
    for (uint32 i=0; i<scaffoldsMax; i++)
      if (scaffolds[i])
        delete scaffolds[i];
    delete [] scaffolds;
  };


  void    addMatch(uint32 scaffoldid, int32 clumpid, uint32 begin, uint32 length) {
    if (scaffoldid >= scaffoldsMax) {
      fprintf(stderr, "ERROR: increase scaffoldsMax "uint32FMT"\n", scaffoldid);
      exit(1);
    }

    if (scaffolds[scaffoldid] == 0L)
      scaffolds[scaffoldid] = new atacClumpCoordTreeScaffold;

    scaffolds[scaffoldid]->addMatch(clumpid, begin, length);
  };


  void    removeSingleClumpScaffolds(void) {
    uint32 deleted = 0;
    uint32 remain = 0;


    for (uint32 i=0; i<scaffoldsMax; i++) {
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
    fprintf(stderr, "Deleted "uint32FMT" scaffolds with less than 2 clumps.\n", deleted);
    fprintf(stderr, "Remain  "uint32FMT" scaffolds with more than 2 clumps.\n", remain);
  };


  void    showMultipleClumpScaffolds(void) {

    for (uint32 i=0; i<scaffoldsMax; i++) {
      if ((scaffolds[i]) && (scaffolds[i]->clumpsLen >= 2)) {

        fprintf(stdout, "\n");

        for (uint32 j=0; j<scaffolds[i]->clumpsLen; j++) {
          bool  overlap = false;

          if ((j+1 < scaffolds[i]->clumpsLen) &&
              (scaffolds[i]->clumpmax[j] > scaffolds[i]->clumpmin[j+1]))
            overlap = true;

          fprintf(stdout, "scaffold "uint32FMT" clump "uint32FMT" begin "uint32FMT" end "uint32FMT"\n",
                  i,
                  scaffolds[i]->clumpID[j],
                  scaffolds[i]->clumpmin[j],
                  scaffolds[i]->clumpmax[j]);

          if (overlap)
            fprintf(stdout, "scaffold "uint32FMT" clump "uint32FMT" and clump "uint32FMT" OVERLAP\n",
                    i,
                    scaffolds[i]->clumpID[j],
                    scaffolds[i]->clumpID[j+1]);

          for (uint32 b=0; b<scaffolds[i]->clumpsLen; b++) {
            uint32 cc = j * scaffolds[i]->clumpsMax + b;
            if (scaffolds[i]->clumpconfirm[cc]) {
              fprintf(stdout, "scaffold "uint32FMT" clump "uint32FMT" and "uint32FMT" confirmed by "uint32FMT" clones.\n",
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


  uint32  getClumpID(uint32 scaffoldid, uint32 begin, uint32 end) {
    if (scaffolds[scaffoldid])
      return(scaffolds[scaffoldid]->getClumpID(begin, end));
    return(0);
  };

  void    confirmClump(uint32 scaffoldid, uint32 ca, uint32 cb) {
    if (scaffolds[scaffoldid]) {
      scaffolds[scaffoldid]->confirm(ca, cb);
    }
  };

  
  uint32                       scaffoldsMax;
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
  map<uint64,uint32>  UIDtoIID;

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

        UIDtoIID[strtouint64(S[13]+1, 0L)] = strtouint32(S[10], 0L);
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


    uint32  scfa = UIDtoIID[strtouint64(A[7], 0L)];
    uint32  scfb = UIDtoIID[strtouint64(B[7], 0L)];

    uint32  cla = ct->getClumpID(scfa,
                                 atoi(A[8]),
                                 atoi(A[9]));
    uint32  clb = ct->getClumpID(scfb,
                                 atoi(B[8]),
                                 atoi(B[9]));

    if (cla == ~uint32ZERO) {
      fprintf(stdout, "%s spans clump in scaffold %s,"uint32FMT"\n", ina, A[7], scfa);
      cla = 0;
    }
    if (clb == ~uint32ZERO) {
      fprintf(stdout, "%s spans clump in scaffold %s,"uint32FMT" \n", inb, B[7], scfb);
      clb = 0;
    }

    if ((cla != 0) &&
        (clb != 0) &&
        (cla != clb)) {
      ct->confirmClump(scfa,
                       (cla < clb) ? cla : clb,
                       (cla < clb) ? clb : cla);

      fprintf(stdout, "scaffold %s,"uint32FMT" clump "uint32FMT" "uint32FMT" confirmed by %s\n",
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
