// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Copyright (c) 2005 The J. Craig Venter Institute
// Author: Clark Mobarry
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>

#include "util++.H"
#include "heavychains.H"


extern "C" {
  void    *construct(char *options);
  void     destruct(void *handle);
  void     addHit(void *handle,
                  char    orientation,
                  u32bit  id1,
                  u32bit  pos1,
                  u32bit  len1,
                  u32bit  id2,
                  u32bit  pos2,
                  u32bit  len2,
                  u32bit  filled);
  void     filter(void *handle);
  u64bit   output(void *handle, FILE *file, u64bit matchid);

  void    *constructStats(char *options);
  void     destructStats(void *handle);
  void     addStats(void  *handle, void *sp);
  void     showStats(void *handle, FILE *file);
}



//  HeavyChains is implemented in the StrandPair class.  It takes all
//  hits for a single pair of sequences and....does something.  Seatac
//  gives the filterObj interface (aka, the interface in this file)
//  all hits for a single sequence to the whole genome (or part of).
//  So, the StrandPairManager acts as the, uhhh, manager for a bunch
//  of StrandPairs, ensuring that each StrandPair is in fact a pair.
//  
//  It is interface compatible with a StrandPair.
//
class StrandPairManager {
private:
  int          beVerbose;
  char         assemblyId1[32];
  char         assemblyId2[32];
  int          maxJump;          // Default maximum intra-run jump allowed in a good run.
  double       minScore;         // Default minimum of bp filled in a good run.

  bool         isForward;

  StrandPair  *P;
  StrandPair  *Proot;
public:
  StrandPairManager(bool   verbose,
                    char  *assemblyid1,
                    char  *assemblyid2,
                    int    maxjump,
                    double minscore) {
    beVerbose    = verbose;
    strncpy(assemblyId1, assemblyid1, 31);
    strncpy(assemblyId2, assemblyid2, 31);
    maxJump      = maxjump;
    minScore     = minscore;

    isForward    = true;

    Proot        = 0L;
    P            = 0L;
  };

  ~StrandPairManager(void) {
    P = Proot;
    while (Proot) {
      Proot = Proot->next();
      delete P;
      P = Proot;
    }
  };

  void addHit(char   direction,
              u32bit id1,
              u32bit xlo,
              u32bit xln,
              u32bit id2,
              u32bit ylo,
              u32bit yln,
              u32bit filled) {

    //  We're given hits for exactly one id2 and all id1, forward hits
    //  followed by reverse hits.  Which means that id1 makes two
    //  passes through, both passes are increasing (enforced by the
    //  seqStream used in seatac).
    //
    //  A linked list of strand pairs is kept (the links are built
    //  into StrandPair for convenience), each strand pair knows it's
    //  pair of ids.
    //

    //  No root?  Make one and add the hit.
    //
    if (Proot == 0L) {
      P = Proot = new StrandPair(beVerbose, assemblyId1, assemblyId2, maxJump, minScore);
      P->addHit(direction, id1, xlo, xln, id2, ylo, yln, filled);
      return;
    }

    //  Reset to the start if we just switched from forward to
    //  reverse.  This is also the only time that the sequence id can
    //  decrease, and we might have to make a new root.
    //
    if (isForward && (direction == 'r')) {
      isForward = false;

      if (id1 < Proot->sequenceIID1()) {
        StrandPair *N = new StrandPair(beVerbose, assemblyId1, assemblyId2, maxJump, minScore);
        N->addHit(direction, id1, xlo, xln, id2, ylo, yln, filled);
        N->addNext(Proot);
        P = Proot = N;
        return;
      }

      P = Proot;
    }

    //  Verify that id1 didn't decrease.
    //
    if (id1 < P->sequenceIID1()) {
      fprintf(stderr, "Why did the sequence id just decrease?  This should not have happened.\n");
      fprintf(stderr, "Crash.  %s at line %d\n", __FILE__, __LINE__ - 2);
      exit(1);
    }

    //  Move to the node just before, or exactly at, the one we want
    //  to add to.  Remember, id1 never decreases.
    //
    while ((P->next()) && (P->next()->sequenceIID1() <= id1))
      P = P->next();

    //  If we're not at the correct node, insert one after the
    //  current, and make it the correct one.
    //
    if (P->sequenceIID1() != id1) {
      StrandPair *NP = new StrandPair(beVerbose, assemblyId1, assemblyId2, maxJump, minScore);
      NP->addNext(P->next());
      P->addNext(NP);
      P = NP;  //  Hooray!
    }

    //  And now we can just add the hit.
    //
    P->addHit(direction, id1, xlo, xln, id2, ylo, yln, filled);
  };

  void process(void) {
    for (StrandPair *SP=Proot; SP; SP=SP->next())
      SP->process();
  };

  u64bit print(FILE *outF, u64bit matchid) {
    for (StrandPair *SP=Proot; SP; SP=SP->next())
      matchid = SP->print(outF, matchid);
    return(matchid);
  };

  void addStats(TheStats *ST) {
    for (StrandPair *SP=Proot; SP; SP=SP->next())
      ST->add(SP);
  };
};








void*
construct(char *options) {
  int    beVerbose       = 0;
  char   assemblyIdD[4]  = { 'U', 'N', 'K', 0 };
  char  *assemblyId1     = assemblyIdD;
  char  *assemblyId2     = assemblyIdD;
  double minScore        = 100.0;   // Default minimum of bp filled in a good run.
  int    maxJump         = 100000;  // Default maximum intra-run jump allowed in a good run.

  //  Parse the options to find the parameters
  //
  splitToWords  W(options);

  u32bit arg = 0;
  while (arg < W.numWords()) {
    if        (strcmp(W.getWord(arg), "-v") == 0) {
      beVerbose++;
    } else if (strcmp(W.getWord(arg), "-s") == 0) {
      minScore = atof(W.getWord(++arg));
    } else if (strcmp(W.getWord(arg), "-j") == 0) {
      maxJump = atoi(W.getWord(++arg));
    } else if (strcmp(W.getWord(arg), "-1") == 0) {
      assemblyId1 = W.getWord(++arg);
    } else if (strcmp(W.getWord(arg), "-2") == 0) {
      assemblyId2 = W.getWord(++arg);
    }

    arg++;
  }

  return((void *)(new StrandPairManager(beVerbose, assemblyId1, assemblyId2, maxJump, minScore)));
}

void
destruct(void *handle) {
  delete (StrandPairManager *)handle;
}

void
addHit(void   *handle,
       char    orientation,
       u32bit  id1,
       u32bit  pos1,
       u32bit  len1,
       u32bit  id2,
       u32bit  pos2,
       u32bit  len2,
       u32bit  filled) {
  ((StrandPairManager *)handle)->addHit(orientation, id1, pos1, len1, id2, pos2, len2, filled);
}

void
filter(void *handle) {
  ((StrandPairManager *)handle)->process();
}


u64bit
output(void *handle, FILE *file, u64bit matchid) {
  return(((StrandPairManager *)handle)->print(file, matchid));
}






void*
constructStats(char *options) {
  int    beVerbose       = 0;
  char   assemblyIdD[4]  = { 'U', 'N', 'K', 0 };
  char  *assemblyId1     = assemblyIdD;
  char  *assemblyId2     = assemblyIdD;
  double minScore        = 100.0;   // Default minimum of bp filled in a good run.
  int    maxJump         = 100000;  // Default maximum intra-run jump allowed in a good run.

  //  Parse the options to find the parameters
  //
  splitToWords  W(options);

  u32bit arg = 0;
  while (arg < W.numWords()) {
    if        (strcmp(W.getWord(arg), "-v") == 0) {
      beVerbose++;
    } else if (strcmp(W.getWord(arg), "-s") == 0) {
      minScore = atof(W.getWord(++arg));
    } else if (strcmp(W.getWord(arg), "-j") == 0) {
      maxJump = atoi(W.getWord(++arg));
    } else if (strcmp(W.getWord(arg), "-1") == 0) {
      assemblyId1 = W.getWord(++arg);
    } else if (strcmp(W.getWord(arg), "-2") == 0) {
      assemblyId2 = W.getWord(++arg);
    }

    arg++;
  }

  return((void *)(new TheStats(beVerbose, assemblyId1, assemblyId2, maxJump, minScore)));
}

void
destructStats(void *handle) {
  delete (TheStats *)handle;
}

void
addStats(void   *handle, void *sp) {

  //  We aren't getting a single StrandPair anymore, we're getting a StrandPairManager now.
  //
  //((TheStats *)handle)->add((StrandPair *)sp);
  //
  ((StrandPairManager *)sp)->addStats((TheStats *)handle);
}

void
showStats(void *handle, FILE *file) {
  ((TheStats *)handle)->show(file);
}
