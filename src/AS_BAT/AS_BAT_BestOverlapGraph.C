
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: AS_BAT_BestOverlapGraph.C,v 1.4 2010-12-07 00:25:56 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"

const uint64 ogMagicNumber   = 0x72476c764f747362llu;  //  'bstOvlGr'
const uint64 ogVersionNumber = 2;

//  Checkpointing will save the best overlap graph after it is loaded from the overlap store.  It
//  seems to be working, but is missing a few features:
//
//  1) It has no version number of any other verification of the data (like, is this a checkpoint
//     file).
//
//  2) It doesn't remember the erate/elimit used when reading overlaps.  Changing these command line
//     parameters with a checkpoint file will result in the original overlaps being used.
//
//  3) It doesn't checkpoint _enough_ (maybe).  It is storing only best overlaps, not the state of
//     bog at the time.  It might be useful to save the gatekeeper information too.


BestOverlapGraph::BestOverlapGraph(OverlapStore        *ovlStoreUniq,
                                   OverlapStore        *ovlStoreRept,
                                   double               AS_UTG_ERROR_RATE,
                                   double               AS_UTG_ERROR_LIMIT,
                                   const char          *prefix) {
  OVSoverlap olap;

  _best5 = new BestEdgeOverlap [FI->numFragments() + 1];
  _best3 = new BestEdgeOverlap [FI->numFragments() + 1];
  _bestC = new BestContainment [FI->numFragments() + 1];

  memset(_best5, 0, sizeof(BestEdgeOverlap) * (FI->numFragments() + 1));
  memset(_best3, 0, sizeof(BestEdgeOverlap) * (FI->numFragments() + 1));
  memset(_bestC, 0, sizeof(BestContainment) * (FI->numFragments() + 1));

  assert(AS_UTG_ERROR_RATE >= 0.0);
  assert(AS_UTG_ERROR_RATE <= AS_MAX_ERROR_RATE);

  assert(AS_CNS_ERROR_RATE >= 0.0);
  assert(AS_CNS_ERROR_RATE <= AS_MAX_ERROR_RATE);

  mismatchCutoff  = AS_OVS_encodeQuality(AS_UTG_ERROR_RATE);
  consensusCutoff = AS_OVS_encodeQuality(AS_CNS_ERROR_RATE);

  mismatchLimit   = AS_UTG_ERROR_LIMIT;

  if (load(prefix, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT)) {
    logFileOrder++;  //  To keep indices the same on log names
    setLogFile(prefix, NULL);
    return;
  }

  setLogFile(prefix, "bestoverlapgraph");

  _bestCscore = new uint64 [FI->numFragments() + 1];
  _best5score = new uint64 [FI->numFragments() + 1];
  _best3score = new uint64 [FI->numFragments() + 1];

  memset(_bestCscore, 0, sizeof(uint64) * (FI->numFragments() + 1));
  memset(_best5score, 0, sizeof(uint64) * (FI->numFragments() + 1));
  memset(_best3score, 0, sizeof(uint64) * (FI->numFragments() + 1));

  AS_OVS_resetRangeOverlapStore(ovlStoreUniq);
  while  (AS_OVS_readOverlapFromStore(ovlStoreUniq, &olap, AS_OVS_TYPE_OVL)) {
    scoreContainment(olap);
  }

  if (ovlStoreRept) {
    AS_OVS_resetRangeOverlapStore(ovlStoreRept);
    while  (AS_OVS_readOverlapFromStore(ovlStoreRept, &olap, AS_OVS_TYPE_OVL)) {
      scoreContainment(olap);
    }
  }

  //  Until we get a list of the fragments that are contained, we must make two passes
  //  through the ovlStore.  The first pass really just marks fragments as contained,
  //  the second pass can then load the overlaps.  The issue is that we cannot load
  //  as a best edge an overlap between a non-contained and a contained fragment.  The
  //  only way to see if that fragment (the B fragment, for argument) is contained is
  //  to load all its overlaps -- and so if A < B, we won't know if B is contained until
  //  too late.

  AS_OVS_resetRangeOverlapStore(ovlStoreUniq);
  while  (AS_OVS_readOverlapFromStore(ovlStoreUniq, &olap, AS_OVS_TYPE_OVL)) {
    //scoreContainment(olap);
    scoreEdge(olap);
  }

  if (ovlStoreRept) {
    AS_OVS_resetRangeOverlapStore(ovlStoreRept);
    while  (AS_OVS_readOverlapFromStore(ovlStoreRept, &olap, AS_OVS_TYPE_OVL)) {
      //scoreContainment(olap);
      scoreEdge(olap);
    }
  }

  delete [] _bestCscore;
  delete [] _best5score;
  delete [] _best3score;

  _bestCscore = NULL;
  _best5score = NULL;
  _best3score = NULL;

  //  Remove dovetail overlaps for contained fragments.

  for (uint32 id=1; id<FI->numFragments() + 1; id++) {
    if (isContained(id) == true) {
      getBestEdgeOverlap(id, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(id, true) ->set(0, 0, 0, 0);
    }
  }

  setLogFile(prefix, NULL);

  //  Diagnostic.  Dump the best edges, count the number of contained reads, etc.
  {
    FILE *BC = fopen("best.contains", "w");
    FILE *BE = fopen("best.edges", "w");
    FILE *BS = fopen("best.singletons", "w");

    if ((BC) && (BE)) {
      fprintf(BC, "#fragId\tlibId\tmated\tbestCont\n");
      fprintf(BE, "#fragId\tlibId\tbest5\tbest3\n");
      fprintf(BS, "#fragId\tlibId\tmated\n");

      for (uint32 id=1; id<FI->numFragments() + 1; id++) {
        BestContainment *bestcont  = getBestContainer(id);
        BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
        BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

        if (bestcont)
          fprintf(BC, "%u\t%u\t%c\t%u\n", id, FI->libraryIID(id), (FI->mateIID(id) > 0) ? 'm' : 'f', bestcont->container);
        else if ((bestedge5->fragId() > 0) || (bestedge3->fragId() > 0))
          fprintf(BE, "%u\t%u\t%u\t%c'\t%u\t%c'\n", id, FI->libraryIID(id),
                  bestedge5->fragId(), bestedge5->frag3p() ? '3' : '5',
                  bestedge3->fragId(), bestedge3->frag3p() ? '3' : '5');
        else
          fprintf(BS, "%u\t%u\t%c\n", id, FI->libraryIID(id), (FI->mateIID(id) > 0) ? 'm' : 'f');

      }

      fclose(BC);
      fclose(BE);
    }
  }

  save(prefix, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT);
}

BestOverlapGraph::~BestOverlapGraph(){
  delete[] _best5;
  delete[] _best3;
  delete[] _bestC;
}



void
BestOverlapGraph::scoreContainment(const OVSoverlap& olap) {

  if (isOverlapBadQuality(olap))
    //  Yuck.  Don't want to use this crud.
    return;

  if ((olap.dat.ovl.a_hang == 0) &&
      (olap.dat.ovl.b_hang == 0) &&
      (olap.a_iid > olap.b_iid))
    //  Exact!  Each contains the other.  Make the lower IID the container.
    return;

  if ((olap.dat.ovl.a_hang < 0) ||
      (olap.dat.ovl.b_hang > 0))
    //  We only save if A contains B.
    return;

  uint64           newScr = scoreOverlap(olap);
  BestContainment      *c = &_bestC[olap.b_iid];

  assert(newScr > 0);

  if (newScr > _bestCscore[olap.b_iid]) {
    c->container         = olap.a_iid;
    c->isContained       = true;
    c->sameOrientation   = olap.dat.ovl.flipped ? false : true;
    c->a_hang            = olap.dat.ovl.a_hang;
    c->b_hang            = olap.dat.ovl.b_hang;

    _bestCscore[olap.b_iid] = newScr;
  }
}



void
BestOverlapGraph::scoreEdge(const OVSoverlap& olap) {

  if (isOverlapBadQuality(olap))
    //  Yuck.  Don't want to use this crud.
    return;

  if (((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) ||
      ((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0)))
    //  Skip containment overlaps.
    return;

  if ((isContained(olap.a_iid) == true) ||
      (isContained(olap.b_iid) == true))
    //  Skip contained fragments.
    return;

  uint64           newScr = scoreOverlap(olap);
  bool             a3p    = AS_OVS_overlapAEndIs3prime(olap);
  BestEdgeOverlap *best   = getBestEdgeOverlap(olap.a_iid, a3p);
  uint64          *score  = (a3p) ? (_best3score + olap.a_iid) : (_best5score + olap.a_iid);

  assert(newScr > 0);

  if (newScr <= *score)
    return;

  best->set(olap);

  *score = newScr;
}



void
BestOverlapGraph::save(const char *prefix, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT) {
  char name[FILENAME_MAX];

  sprintf(name, "%s.bog", prefix);

  assert(_best5score == NULL);
  assert(_best3score == NULL);
  assert(_bestCscore == NULL);

  errno = 0;
  FILE *file = fopen(name, "w");
  if (errno) {
    fprintf(logFile, "BestOverlapGraph-- Failed to open '%s' for writing: %s\n", name, strerror(errno));
    fprintf(logFile, "BestOverlapGraph-- Will not save best overlap graph to cache.\n");
    return;
  }

  fprintf(logFile, "BestOverlapGraph()-- Saving overlap graph to '%s'.\n",
          name);

  AS_UTL_safeWrite(file, &ogMagicNumber,      "magicnumber",   sizeof(uint64),              1);
  AS_UTL_safeWrite(file, &ogVersionNumber,    "versionnumber", sizeof(uint64),              1);

  AS_UTL_safeWrite(file, &AS_UTG_ERROR_RATE,  "errorRate",     sizeof(double),              1);
  AS_UTL_safeWrite(file, &AS_UTG_ERROR_LIMIT, "errorLimit",    sizeof(double),              1);

  AS_UTL_safeWrite(file, _best5, "best overlaps 5", sizeof(BestEdgeOverlap), FI->numFragments() + 1);
  AS_UTL_safeWrite(file, _best3, "best overlaps 3", sizeof(BestEdgeOverlap), FI->numFragments() + 1);
  AS_UTL_safeWrite(file, _bestC, "best contains C", sizeof(BestContainment), FI->numFragments() + 1);

  fclose(file);
}

bool
BestOverlapGraph::load(const char *prefix, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT) {
  char name[FILENAME_MAX];

  sprintf(name, "%s.bog", prefix);

  errno = 0;
  FILE *file = fopen(name, "r");
  if (errno)
    return(false);

  assert(_best5 != NULL);
  assert(_best3 != NULL);
  assert(_bestC != NULL);

  uint64 magicNumber;
  uint64 versionNumber;

  AS_UTL_safeRead(file, &magicNumber,   "magicnumber",   sizeof(uint64), 1);
  AS_UTL_safeRead(file, &versionNumber, "versionnumber", sizeof(uint64), 1);

  if (magicNumber != ogMagicNumber) {
    fprintf(logFile, "BestOverlapGraph()-- File '%s' is not a best overlap graph; cannot load graph.\n", name);
    fclose(file);
    return(false);
  }
  if (versionNumber != ogVersionNumber) {
    fprintf(logFile, "BestOverlapGraph()-- File '%s' is version "F_U64", I can only read version "F_U64"; cannot load graph.\n",
            name, versionNumber, ogVersionNumber);
    fclose(file);
    return(false);
  }

  fprintf(logFile, "BestOverlapGraph()-- Loading overlap graph from '%s'.\n", name);

  double  eRate  = 0.0;
  double  eLimit = 0.0;

  AS_UTL_safeRead(file, &eRate,  "errorRate",     sizeof(double), 1);
  AS_UTL_safeRead(file, &eLimit, "errorLimit",    sizeof(double), 1);

  if (eRate  != AS_UTG_ERROR_RATE)
    fprintf(logFile, "BestOverlapGraph()-- Saved graph in '%s' has error rate %f, this run is expecting error rate %f; cannot load graph.\n",
            name, eRate, AS_UTG_ERROR_RATE);
  if (eLimit != AS_UTG_ERROR_LIMIT)
    fprintf(logFile, "BestOverlapGraph()-- Saved graph in '%s' has error limit %f, this run is expecting error limit %f; cannot load graph.\n",
            name, eLimit, AS_UTG_ERROR_LIMIT);
  if ((eRate  != AS_UTG_ERROR_RATE) ||
      (eLimit != AS_UTG_ERROR_LIMIT)) {
    fclose(file);
    return(false);
  }

  AS_UTL_safeRead(file, _best5, "best overlaps", sizeof(BestEdgeOverlap), FI->numFragments() + 1);
  AS_UTL_safeRead(file, _best3, "best overlaps", sizeof(BestEdgeOverlap), FI->numFragments() + 1);
  AS_UTL_safeRead(file, _bestC, "best contains", sizeof(BestContainment), FI->numFragments() + 1);

  fclose(file);

  return(true);
}
