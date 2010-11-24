
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

static const char *rcsid = "$Id: AS_BAT_BestOverlapGraph.C,v 1.1 2010-11-24 01:03:31 brianwalenz Exp $";

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
    logFileOrder += 2;  //  To keep indices the same on log names
    setLogFile(prefix, NULL);
    return;
  }

  //  Pass 1 through overlaps -- find the contained fragments.

  setLogFile(prefix, "bestoverlapgraph-containments");

  _bestCscore    = new uint64 [FI->numFragments() + 1];
  memset(_bestCscore,    0, sizeof(uint64) * (FI->numFragments() + 1));

  AS_OVS_resetRangeOverlapStore(ovlStoreUniq);
  while  (AS_OVS_readOverlapFromStore(ovlStoreUniq, &olap, AS_OVS_TYPE_OVL))
    scoreContainment(olap);

  if (ovlStoreRept) {
    AS_OVS_resetRangeOverlapStore(ovlStoreRept);
    while  (AS_OVS_readOverlapFromStore(ovlStoreRept, &olap, AS_OVS_TYPE_OVL))
      scoreContainment(olap);
  }

  delete [] _bestCscore;
  _bestCscore    = NULL;


  //  Report some statistics on overlaps
  {
    uint64  numContainsToSave = 0;

    for (uint32 i=0; i<FI->numFragments(); i++)
      numContainsToSave += _bestC[i].olapsLen;

    fprintf(logFile, "Need to save "F_U64" near-containment overlaps\n", numContainsToSave);
  }


  //  Pass 2 through overlaps -- find dovetails, build the overlap graph.  For each
  //  contained fragment, remember some of the almost containment overlaps.

  setLogFile(prefix, "bestoverlapgraph-dovetails");

  _best5score = new uint64 [FI->numFragments() + 1];
  _best3score = new uint64 [FI->numFragments() + 1];

  memset(_best5score, 0, sizeof(uint64) * (FI->numFragments() + 1));
  memset(_best3score, 0, sizeof(uint64) * (FI->numFragments() + 1));

  AS_OVS_resetRangeOverlapStore(ovlStoreUniq);
  while  (AS_OVS_readOverlapFromStore(ovlStoreUniq, &olap, AS_OVS_TYPE_OVL))
    scoreEdge(olap);

  if (ovlStoreRept) {
    AS_OVS_resetRangeOverlapStore(ovlStoreRept);
    while  (AS_OVS_readOverlapFromStore(ovlStoreRept, &olap, AS_OVS_TYPE_OVL))
      scoreEdge(olap);
  }

  delete [] _best5score;
  delete [] _best3score;

  _best5score = NULL;
  _best3score = NULL;

  setLogFile(prefix, NULL);


  //  Clean up our allocation.  We seem to over count the number of overlaps in the first pass,
  //  then don't add any overlaps in the second pass, leaving the count positive and the pointer
  //  NULL.  In the case where we just overcount, the pointer is valid, and the count is correct.
  //
  //  A better explanation is that we count near containment overlaps for ALL fragments in the
  //  first pass, but then only save near containment overlaps for contained fragments, leaving
  //  the dovetail fragments with a positive count and a NULL pointer.
  //
  for (uint32 i=0; i<FI->numFragments() + 1; i++)
    if (_bestC[i].olaps == NULL)
      _bestC[i].olapsLen = 0;


  //  Diagnostic.  Dump the best edges, count the number of contained
  //  reads, etc.
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



void BestOverlapGraph::scoreContainment(const OVSoverlap& olap) {

  if (isOverlapBadQuality(olap))
    return;

  //  Count the number of good containment and near-containment
  //  overlaps the B fragment has -- used by scoreEdge (see the
  //  comment there) to keep a list of dovetail overlaps to contained
  //  fragments.
  //
  if (((olap.dat.ovl.a_hang >= -10) && (olap.dat.ovl.b_hang <=  0)) ||
      ((olap.dat.ovl.a_hang >=   0) && (olap.dat.ovl.b_hang <= 10)))
    _bestC[olap.b_iid].olapsLen++;

  //  In the case of no hang, make the lower frag the container
  //
  if ((olap.dat.ovl.a_hang == 0) &&
      (olap.dat.ovl.b_hang == 0) &&
      (olap.a_iid > olap.b_iid))
    return;

  //  We only care if A contains B.

  if ((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) {
    uint64           newScr = scoreOverlap(olap);
    BestContainment      *c = &_bestC[olap.b_iid];

    if (newScr > _bestCscore[olap.b_iid]) {
      //  NOTE!  This is already initialized.  We do not need to, and
      //  it is an error to, initialize olaps to zero!  (We're
      //  counting olapsLen above, see?  This stupid bug took me about
      //  an hour to find, grrr.)
      c->container         = olap.a_iid;
      c->a_hang            = olap.dat.ovl.a_hang;
      c->b_hang            = olap.dat.ovl.b_hang;
      c->sameOrientation   = olap.dat.ovl.flipped ? false : true;
      c->isContained       = true;
      c->olapsSorted       = false;

      _bestCscore[olap.b_iid] = newScr;
    }
  }
}



void BestOverlapGraph::scoreEdge(const OVSoverlap& olap) {

  if (isOverlapBadQuality(olap))
    return;

  //  Store edges from contained frags to help with unhappy mate
  //  splitting.
  //
  //  From Eli: These are contained, but close either way.  We're
  //  storing the non-containment edges for this fragment, plus a few
  //  containment edges that are "close" to being dovetails.  "I think
  //  there are cases when a change in the alignemtn (consensus) will
  //  change which one is contained and screw up the order, so having
  //  this 10 base fudge factor helps things work out."
  //
  if (isContained(olap.b_iid)) {
    if (((olap.dat.ovl.a_hang >= -10) && (olap.dat.ovl.b_hang <=  0)) ||
        ((olap.dat.ovl.a_hang >=   0) && (olap.dat.ovl.b_hang <= 10))) {
      BestContainment *c = &_bestC[olap.b_iid];
      if (c->olaps == NULL) {
        c->olaps    = new uint32 [c->olapsLen];
        c->olapsLen = 0;
      }
      c->olaps[c->olapsLen++] = olap.a_iid;
    }
    return;
  }

  //  Skip contained fragments.
  if (isContained(olap.a_iid) || isContained(olap.b_iid))
    return;

  //  Skip containment overlaps.  Can this happen?  Yup.  How?
  //  The overlap could be above our allowed error.
  //
  if (((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) ||
      ((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0)))
    return;

  uint64 newScr = scoreOverlap(olap);

  //  If the score is 0, the overlap doesn't pass the scoring
  //  criteria at all so don't store the overlap whether or not
  //  it's dovetailing or containment.

  if (newScr == 0)
    return;

  //  Dove tailing overlap
  bool             a3p     = AS_OVS_overlapAEndIs3prime(olap);
  BestEdgeOverlap *best    = getBestEdgeOverlap(olap.a_iid, a3p);
  uint64           score   = 0;

  // Store the overlap if:
  //   1.)  The score is better than what is already in the graph
  //   2.)  If the scores are identical, the one with the longer length
  //
  // Since the order of how the overlaps are read in from the overlap
  // store are by A's increasing uint32, by default, if the score and
  // length are the same, the uint32 of the lower value will be kept.

  if (a3p)
    score = _best3score[olap.a_iid];
  else
    score = _best5score[olap.a_iid];

  if (newScr > score) {
    best->set(olap);

    if (a3p)
      _best3score[olap.a_iid] = newScr;
    else
      _best5score[olap.a_iid] = newScr;
  }
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
    fprintf(logFile, "BestOverlapGraph-- Will not save best overlap graph to cache.\n", strerror(errno));
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

  for (uint32 i=0; i<FI->numFragments() + 1; i++)
    if (_bestC[i].olaps != NULL)
      AS_UTL_safeWrite(file, _bestC[i].olaps, "best contains olaps", sizeof(uint32), _bestC[i].olapsLen);

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
            versionNumber, ogVersionNumber, name);
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

  for (uint32 i=0; i<FI->numFragments() + 1; i++) {
    if (_bestC[i].olapsLen > 0) {
      _bestC[i].olaps = new uint32 [_bestC[i].olapsLen];
      AS_UTL_safeRead(file, _bestC[i].olaps, "best contains olaps", sizeof(uint32), _bestC[i].olapsLen);
    } else {
      assert(_bestC[i].olaps == NULL);
    }
  }

  fclose(file);

  return(true);
}
