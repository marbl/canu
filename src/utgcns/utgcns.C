
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

const char *mainid = "$Id$";

#include "AS_global.H"
#include "gkStore.H"
#include "tgStore.H"
#include "abAbacus.H"

#include "AS_UTL_decodeRange.H"

#include <map>
#include <algorithm>



//  Create a new f_list for the ma that has no contained reads.
//  The original f_list is returned.
//

class readLength {
public:
  uint32    idx;
  int32     len;

  bool operator<(const readLength &that) const {
    return(len < that.len);
  };
};


class savedChildren {
public:
  savedChildren(tgTig *tig) {
    childrenLen = tig->_childrenLen;
    childrenMax = tig->_childrenMax;
    children    = tig->_children;
  };

  uint32      childrenLen;
  uint32      childrenMax;
  tgPosition *children;
};



//  Replace the children list in tig with one that has fewer contains.  The original
//  list is returned.
savedChildren *
stashContains(tgTig       *tig,
              double       maxCov) {

  if (tig->numberOfChildren() == 1)
    return(NULL);

  //  Stats we report
  int32  nOrig     = tig->numberOfChildren();
  int32  nBack     = 0;
  int32  nCont     = 0;
  int32  nSave     = 0;
  int64  nBase     = 0;
  int64  nBaseDove = 0;
  int64  nBaseCont = 0;
  int64  nBaseSave = 0;

  //  Save the original children
  savedChildren   *saved = new savedChildren(tig);

  bool         *isBack   = new bool       [nOrig];   //  True, we save the child for processing
  readLength   *posLen   = new readLength [nOrig];   //  Sorting by length of child

  //  Sort the original children by position.

  std::sort(saved->children, saved->children + saved->childrenLen);

  //  The first read is always saved

  int32         loEnd = saved->children[0].min();
  int32         hiEnd = saved->children[0].max();

  isBack[0]      = 1;
  nBack          = 1;
  posLen[0].idx  = 0;
  posLen[0].len  = hiEnd - loEnd;
  nBaseDove     += posLen[0].len;
  nBase         += posLen[0].len;

  //  For the other reads, save it if it extends the backbone sequence.

  for (uint32 fi=1; fi<nOrig; fi++) {
    int32  lo = saved->children[fi].min();
    int32  hi = saved->children[fi].max();

    posLen[fi].idx  = fi;
    posLen[fi].len  = hi - lo;
    nBase          += posLen[fi].len;

    if (hi <= hiEnd) {
      isBack[fi] = false;
      nCont++;
      nBaseCont += posLen[fi].len;

    } else {
      isBack[fi] = true;
      nBack++;
      nBaseDove += posLen[fi].len;
    }

    hiEnd = MAX(hi, hiEnd);
  }

  //  Entertain the user with some statistics

  double percCont = 100.0 * nBaseCont / nBase;
  double percDove = 100.0 * nBaseDove / nBase;
  double totlCov  = (double)nBase / hiEnd;

  fprintf(stderr, "  unitig %d detected "F_S32" contains (%.2fx, %.2f%%) "F_S32" dovetail (%.2fx, %.2f%%)\n",
          tig->tigID(),
          nCont, (double)nBaseCont / hiEnd, percCont,
          nBack, (double)nBaseDove / hiEnd, percDove);

  //  If the tig has more coverage than allowed, throw out some of the contained reads.

  if ((totlCov  >= maxCov) &&
      (maxCov   > 0)) {
#warning is this really larger first?
    std::sort(posLen, posLen + nOrig);  //  Sort by length, larger first

    nBaseSave = 0.0;

    for (uint32 ii=0; ((ii < nOrig) && ((double)(nBaseSave + nBaseDove) / hiEnd < maxCov)); ii++) {
      if (isBack[posLen[ii].idx])
        //  Already a backbone read.
        continue;

      isBack[posLen[ii].idx] = true;  //  Save it.

      nSave++;
      nBaseSave += posLen[ii].len;
    }

    fprintf(stderr, "    unitig %d removing "F_S32" (%.2fx) contained reads; processing only "F_S32" contained (%.2fx) and "F_S32" dovetail (%.2fx) reads\n",
            tig->tigID(),
            nOrig - nBack - nSave,
            (double)(nBaseCont - nBaseSave) / hiEnd,
            nSave, (double)nBaseSave / hiEnd,
            nBack, (double)nBaseDove / hiEnd);

    //  For all the reads we saved, copy them to a new children list in the tig

    tig->_childrenLen = 0;
    tig->_childrenMax = nBack + nSave;
    tig->_children    = new tgPosition [tig->_childrenMax];

    for (uint32 fi=0; fi<nOrig; fi++) {
      if (isBack[fi] == false)
        continue;

      //fprintf(stderr, "    ident %9d position %6d %6d\n",
      //        saved->children[fi].ident(), saved->children[fi].bgn(), children[fi].end());

      tig->_children[tig->_childrenLen] = saved->children[fi];
    }
  }

  //  Else, the tig coverage is acceptable and we do no filtering.
  else {
    delete saved;
    saved = NULL;
  }

  delete [] isBack;
  delete [] posLen;

  return(saved);
}


//  Restores the f_list, and updates the position of non-contained reads.
//
void
unstashContains(tgTig                *tig,
                savedChildren        *saved) {

  if (saved == NULL)
    return;

  uint32   oldMax = 0;
  uint32   newMax = 0;

  //  For fragments not involved in the consensus computation, we'll scale their position linearly
  //  from the old max to the new max.
  //
  //  We probably should do an alignment to the consensus sequence to find the true location, but
  //  that's (a) expensive and (b) likely overkill for these unitigs.

  //  Find the oldMax
  for (uint32 fi=0, ci=0; fi<saved->childrenLen; fi++)
    if (oldMax < saved->children[fi].max())
      oldMax = saved->children[fi].max();

  //  Find the newMax
  //  We could have just done: newMax = tig->gappedLength();
  for (uint32 fi=0, ci=0; fi<tig->numberOfChildren(); fi++)
    if (newMax < tig->getChild(fi)->max())
      newMax = tig->getChild(fi)->max();

  double sf = (double)newMax / oldMax;

  //  First, we need a map from the child id to the location in the current tig

  map<int32, tgPosition *>   idmap;

  for (uint32 ci=0; ci < tig->numberOfChildren(); ci++)
    idmap[tig->getChild(ci)->ident()] = tig->getChild(ci);

  //  Now, over all the reads in the original saved fragment list, update the position.  Either from
  //  the computed result, or by extrapolating.

  for (uint32 fi=0; fi<saved->childrenLen; fi++) {
    uint32  iid = saved->children[fi].ident();

    //  Does the ID exist in the new positions?  Copy the new position to the original list.
    if (idmap.find(iid) != idmap.end()) {
      saved->children[fi] = *idmap[iid];
      idmap.erase(iid);
    }

    //  Otherwise, fudge the positions.
    else {
      saved->children[fi].bgn() = sf * saved->children[fi].bgn();
      saved->children[fi].end() = sf * saved->children[fi].end();

      if (saved->children[fi].bgn() > newMax)  saved->children[fi].bgn() = newMax;
      if (saved->children[fi].end() > newMax)  saved->children[fi].end() = newMax;
    }
  }

  if (idmap.empty() == false)
    fprintf(stderr, "Failed to unstash the contained reads.  Still have "F_SIZE_T" reads unplaced.\n",
            idmap.size());
  assert(idmap.empty() == true);

  //  Throw out the reduced list, and restore the original.

  delete [] tig->_children;

  tig->_childrenLen = saved->childrenLen;
  tig->_childrenMax = saved->childrenMax;
  tig->_children    = saved->children;
}



int
main (int argc, char **argv) {
  char  *gkpName = NULL;

  char  *tigName = NULL;
  int32  tigVers = -1;
  int32  tigPart = -1;

  int64  utgBgn = -1;
  int64  utgEnd = -1;
  char  *utgFile = NULL;

  bool   forceCompute = false;

  int32  numFailures = 0;
  int32  numSkipped  = 0;

  bool   showResult = false;

  double maxCov = 0.0;
  uint32 maxLen = UINT32_MAX;

  bool   inplace  = false;
  bool   loadall  = false;
  bool   doUpdate = true;

  uint32 verbosity = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);
      tigPart = atoi(argv[++arg]);

      if (tigVers <= 0)
        fprintf(stderr, "invalid tigStore version (-t store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);
      if ((tigPart <= 0) && (argv[arg][0] != '.'))
        fprintf(stderr, "invalid tigStore partition (-t store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);

    } else if (strcmp(argv[arg], "-u") == 0) {
      AS_UTL_decodeRange(argv[++arg], utgBgn, utgEnd);

    } else if (strcmp(argv[arg], "-T") == 0) {
      utgFile = argv[++arg];

    } else if (strcmp(argv[arg], "-f") == 0) {
      forceCompute = true;

    } else if (strcmp(argv[arg], "-v") == 0) {
      showResult = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      verbosity++;

    } else if (strcmp(argv[arg], "-maxcoverage") == 0) {
      maxCov   = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-maxlength") == 0) {
      maxLen   = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-inplace") == 0) {
      inplace = true;

    } else if (strcmp(argv[arg], "-loadall") == 0) {
      loadall = true;

    } else if (strcmp(argv[arg], "-n") == 0) {
      doUpdate = false;

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if (gkpName == NULL)
    err++;
  if ((utgFile == NULL) && (tigName == NULL))
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version partition [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    -u b            Compute only unitig ID 'b' (must be in the correct partition!)\n");
    fprintf(stderr, "    -u b-e          Compute only unitigs from ID 'b' to ID 'e'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -T file         Test the computation of the unitig layout in 'file'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -f              Recompute unitigs that already have a multialignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -v              Show multialigns.\n");
    fprintf(stderr, "    -V              Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  ADVANCED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -n              Do not update the store after computing consensus.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -maxcoverage c  Use non-contained reads and the longest contained reads, up to\n");
    fprintf(stderr, "                    C coverage, for consensus generation.  The default is 0, and will\n");
    fprintf(stderr, "                    use all reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -maxlength l    Do not compute consensus for unitigs longer than l bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -inplace        Write the updated unitig to the same version it was read from.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -t S V P        If 'partition' is '.', use an unpartitioned tigStore/gkpStore.\n");
    fprintf(stderr, "    -loadall        Load ALL reads into memory.  Ignores partition if it exists.\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR:  No gkpStore (-g) supplied.\n");

    if ((utgFile == NULL) && (tigName == NULL))
      fprintf(stderr, "ERROR:  No tigStore (-t) OR no test unitig (-T) supplied.\n");

    exit(1);
  }

  //  Open gatekeeper for read only, and load the partitioned data if tigPart > 0.

  fprintf(stderr, "Opening gkpStore.\n");
  gkStore   *gkpStore = new gkStore(gkpName, gkStore_readOnly, tigPart);

  //  Create a consensus object.

  fprintf(stderr, "Creating abacus.\n");
  abAbacus  *abacus   = new abAbacus(gkpStore);

  //  If we are testing a unitig, do that.

  if (utgFile != NULL) {
    fprintf(stderr, "utgFile not supported.\n");
    exit(1);
#if 0
    errno = 0;
    FILE  *F = fopen(utgFile, "r");
    if (errno)
      fprintf(stderr, "Failed to open input unitig file '%s': %s\n", utgFile, strerror(errno)), exit(1);

    MultiAlignT  *ma       = CreateEmptyMultiAlignT();
    bool          isUnitig = false;

    while (LoadMultiAlignFromHuman(tig, isUnitig, F) == true) {
      if (generateMultiAlignment(tig, gkpStore, NULL)) {
        if (showResult)
          abacus->getMultiAlignment()->printAlignment(abacus, stdout);

      } else {
        fprintf(stderr, "tig %d failed.\n", tig->tigID());
        numFailures++;
      }
    }

    DeleteMultiAlignT(ma);
#endif

    exit(0);
  }

  //  Otherwise, we're computing unitigs from the store.  Open it for read only.
  //  Outputs get written to a single output file.

  fprintf(stderr, "Opening tigStore.\n");
  tgStore *tigStore = new tgStore(tigName, tigVers);

  //  Decide on what to compute.  Either all unitigs, or a single unitig, or a special case test.

  uint32  b = 0;
  uint32  e = tigStore->numTigs();

  if (utgBgn != -1) {
    b = utgBgn;
    e = utgEnd + 1;
  }

  //  Reopen for writing, if we have work to do.

  //if (b < e) {
  //  delete tigStore;
  //  tigStore = new MultiAlignStore(tigName, tigVers, tigPart, 0, doUpdate, inplace, !inplace);
  //}

  fprintf(stderr, "Computing unitig consensus for b="F_U32" to e="F_U32"\n", b, e);

  //  Now the usual case.  Iterate over all unitigs, compute and update.

  for (uint32 ti=b; ti<e; ti++) {
    tgTig  *tig = tigStore->loadTig(ti);

    if (tig == NULL) {
      //  Not in our partition, or deleted.
      continue;
    }

    bool exists = (tig->gappedLength() > 0);

    if ((forceCompute == false) && (exists == true)) {
      //  Already finished unitig consensus.
      if (tig->numberOfChildren() > 1)
        fprintf(stderr, "Working on unitig %d of length %d (%d children) - already computed, skipped\n",
                tig->tigID(), tig->layoutLength(), tig->numberOfChildren());
      numSkipped++;
      continue;
    }

    if (tig->layoutLength() > maxLen) {
      fprintf(stderr, "SKIP unitig %d of length %d (%d children) - too long, skipped\n",
              tig->tigID(), tig->layoutLength(), tig->numberOfChildren());
      continue;
    }

    if (tig->numberOfChildren() > 1)
      fprintf(stderr, "Working on unitig %d of length %d (%d children)%s\n",
              tig->tigID(), tig->layoutLength(), tig->numberOfChildren(),
              (exists) ? " - already computed, recomputing" : "");

    //  Build a new ma if we're ignoring contains.  We'll need to put back the reads we remove
    //  before we add it to the store.

    savedChildren *origChildren = stashContains(tig, maxCov);

    tig->_utgcns_verboseLevel = verbosity;
    tig->_utgcns_smoothWindow = 11;
    tig->_utgcns_splitAlleles = false;
    tig->_utgcns_doPhasing    = false;

    if (generateMultiAlignment(tig, gkpStore, NULL)) {
      if (showResult)
        //abacus->getMultiAlign()->printAlignment(abacus, stdout);
        //  

      unstashContains(tig, origChildren);

      //if (doUpdate) {
      //  tigStore->insertMultiAlign(ma, true, true);
      //  tigStore->unloadMultiAlign(ma->maID, true, false);
      //} else {
      //  tigStore->unloadMultiAlign(ma->maID, true, true);
      //}

    } else {
      fprintf(stderr, "unitigConsensus()-- unitig %d failed.\n", tig->tigID());
      numFailures++;
    }
  }

 finish:
  delete tigStore;

#if 0
  fprintf(stderr, "\n");
  fprintf(stderr, "NumColumnsInUnitigs             = %d\n", NumColumnsInUnitigs);
  fprintf(stderr, "NumGapsInUnitigs                = %d\n", NumGapsInUnitigs);
  fprintf(stderr, "NumRunsOfGapsInUnitigReads      = %d\n", NumRunsOfGapsInUnitigReads);
  fprintf(stderr, "NumColumnsInContigs             = %d\n", NumColumnsInContigs);
  fprintf(stderr, "NumGapsInContigs                = %d\n", NumGapsInContigs);
  fprintf(stderr, "NumRunsOfGapsInContigReads      = %d\n", NumRunsOfGapsInContigReads);
  fprintf(stderr, "NumAAMismatches                 = %d\n", NumAAMismatches);
  fprintf(stderr, "NumVARRecords                   = %d\n", NumVARRecords);
  fprintf(stderr, "NumVARStringsWithFlankingGaps   = %d\n", NumVARStringsWithFlankingGaps);
  fprintf(stderr, "NumUnitigRetrySuccess           = %d\n", NumUnitigRetrySuccess);
  fprintf(stderr, "\n");
#endif

  if (numFailures) {
    fprintf(stderr, "WARNING:  Total number of unitig failures = %d\n", numFailures);
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus did NOT finish successfully.\n");
  } else {
    fprintf(stderr, "Consensus finished successfully.  Bye.\n");
  }

  return(0);
}
