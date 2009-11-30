
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

static const char *rcsid = "$Id: AS_BOG_UnitigGraph.cc,v 1.119 2009-11-30 20:37:39 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"

#undef max

#define VERBOSEBUILD
#define VERBOSEBREAK


void
UnitigGraph::checkUnitigMembership(void) {
  int nutg = 0;
  int nfrg = 0;

  fprintf(stderr, "checkUnitigMembership()--  numfrags=%d\n", _fi->numFragments());

  uint32 *inUnitig = new uint32 [_fi->numFragments()+1];

  for (uint32 i=0; i<_fi->numFragments()+1; i++)
    inUnitig[i] = 987654321;

  for (uint32 ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg) {
      nutg++;
      for (DoveTailIter it=utg->dovetail_path_ptr->begin(); it != utg->dovetail_path_ptr->end(); it++) {
        if (it->ident > _fi->numFragments())
          fprintf(stderr, "HUH?  ident=%d numfrags=%d\n", it->ident, _fi->numFragments());
        inUnitig[it->ident] = ti;
        nfrg++;
      }
    }
  }

  int lost = 0;
  int found = 0;

  for (uint32 i=0; i<_fi->numFragments()+1; i++) {
    if (_fi->fragmentLength(i) > 0) {
      if (inUnitig[i] == 0) {
        fprintf(stderr, "ERROR frag %d is in unitig 0!\n", i);
      } else if (inUnitig[i] != 987654321) {
        found++;
      } else {
        fprintf(stderr, "ERROR frag %d disappeared!\n", i);
        lost++;
      }
    }
  }

  fprintf(stderr, "checkUnitigMembership()-- nutg=%d nfrg=%d lost=%d found=%d\n", nutg, nfrg, lost, found);

  assert(lost == 0);
};


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
UnitigGraph::reportOverlapsUsed(const char *filename) {

  return;

  FILE *F = fopen(filename, "w");

  if (F == NULL)
    return;

  for (uint32  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode  *frg = &(*utg->dovetail_path_ptr)[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);

      uint32           bestident5 = 0;
      uint32           bestident3 = 0;

      if (bestedge5)
        bestident5 = bestedge5->frag_b_id;

      if (bestedge3)
        bestident3 = bestedge3->frag_b_id;

      //  Now search ahead, reporting any overlap to any fragment.
      //
      for (uint32 oi=fi+1; oi<utg->dovetail_path_ptr->size(); oi++) {
        DoveTailNode  *ooo = &(*utg->dovetail_path_ptr)[oi];

        int frgbgn = MIN(frg->position.bgn, frg->position.end);
        int frgend = MAX(frg->position.bgn, frg->position.end);

        int ooobgn = MIN(ooo->position.bgn, ooo->position.end);
        int oooend = MAX(ooo->position.bgn, ooo->position.end);

        if ((frgbgn <= ooobgn) && (ooobgn + 40 < frgend)) {
          BestContainment *bestcont  = bog_ptr->getBestContainer(ooo->ident);

          uint32           bestident  = 0;
          if (bestcont)
            bestident = bestcont->container;

          bool isBest = ((frg->ident == bestident) ||
                         (ooo->ident == bestident5) ||
                         (ooo->ident == bestident3));

          fprintf(F, "%d\t%d%s\n", frg->ident, ooo->ident, (isBest) ? ((bestident) ? "\tbc" : "\tbe") : "");
        }

        if (frgend < ooobgn)
          break;
      }
    }
  }

  fclose(F);
}


UnitigGraph::UnitigGraph(FragmentInfo *fi, BestOverlapGraph *bp) {
  unitigs = new UnitigVector;
  _fi     = fi;
  bog_ptr = bp;
}

UnitigGraph::~UnitigGraph() {
  UnitigVector::iterator utg_itr;
  for(utg_itr=unitigs->begin(); utg_itr!=unitigs->end(); utg_itr++)
    delete *utg_itr;
  delete unitigs;
}


void UnitigGraph::build(ChunkGraph *cg_ptr, bool enableIntersectionBreaking, bool enableBubblePopping, char *output_prefix) {
  bool verbose = false;

  // Initialize where we've been to nowhere
  Unitig::resetFragUnitigMap( _fi->numFragments() );

  // Invert the containment map to key by container, instead of containee
  ContainerMap      cMap;

  for (uint32 c=0; c<=_fi->numFragments(); c++)
    if (bog_ptr->isContained(c))
      cMap[bog_ptr->getBestContainer(c)->container].push_back(c);

  // Step through all the fragments

  fprintf(stderr, "==> BUILDING UNITIGS from %d fragments.\n", _fi->numFragments());

  //  There is no 0th unitig.
  unitigs->push_back(NULL);

  uint32 frag_idx;
  while ((frag_idx = cg_ptr->nextFragByChunkLength()) > 0) {
    if ((_fi->fragmentLength(frag_idx) == 0) ||
        (Unitig::fragIn(frag_idx) != 0) ||
        (bog_ptr->isContained(frag_idx) == true))
      //  Skip deleted fragments, already placed fragments, and contained fragments.
      continue;

    populateUnitig(frag_idx, false);
  }

  reportOverlapsUsed("overlaps.afterbuild");

  fprintf(stderr, "==> BUILDING UNITIGS catching missed fragments.\n");

  //  Pick up frags missed above, leftovers from possibly circular unitigs

  for (frag_idx=1; frag_idx <= _fi->numFragments(); frag_idx++) {
    if ((_fi->fragmentLength(frag_idx) == 0) ||
        (Unitig::fragIn(frag_idx) != 0) ||
        (bog_ptr->isContained(frag_idx) == true))
      continue;

    populateUnitig(frag_idx, false);
  }

  reportOverlapsUsed("overlaps.afterbuild2");

  if (enableIntersectionBreaking)
    breakUnitigs(cMap, output_prefix);

  reportOverlapsUsed("overlaps.afterbreak");

  placeContains();

  reportOverlapsUsed("overlaps.aftercontains");

  placeZombies();

  checkUnitigMembership();

  if (enableBubblePopping)
    popBubbles();

  checkUnitigMembership();
}


void
UnitigGraph::setParentAndHang(ChunkGraph *cg) {

  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig        *utg = (*unitigs)[ti];
    DoveTailNode  *frg;

    if (utg == NULL)
      continue;

    if (utg->dovetail_path_ptr->size() == 0)
      continue;

    //  First fragment has no parent or hangs set.
    frg = &(*utg->dovetail_path_ptr)[0];

    frg->parent       = 0;
    frg->ahang        = 0;
    frg->bhang        = 0;

    //  All other fragments can have a parent, if the best overlap is
    //  in our unitig.
    //
    for (int fi=1; fi<utg->dovetail_path_ptr->size(); fi++) {
      frg = &(*utg->dovetail_path_ptr)[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);
      BestContainment *bestcont  = bog_ptr->getBestContainer(frg->ident);

      if (frg->contained) {
        //  Fragment is contained.  Don't use the edges off the ends.
        bestedge5 = NULL;
        bestedge3 = NULL;
      } else {
        //  Fragment is dovetail.  Don't use the best containment.
        bestcont == NULL;
      }

      if (frg->position.bgn < frg->position.end) {
        //  Fragment is forward.  Looking for the best edge off of our 5' end.
        bestedge3 = NULL;
      } else {
        //  Fragment is reverse.  Looking for the best edge off of our 3' end.
        bestedge5 = NULL;
      }

      //  Initialize to 'unset'.
      frg->parent       = 0;
      frg->ahang        = 0;
      frg->bhang        = 0;

      //  Now search for the correct overlap
      //
      for (uint32 oi=fi; oi--; ) {
        DoveTailNode  *ooo = &(*utg->dovetail_path_ptr)[oi];

        if (bestcont) {
          //  Frag is contained.

          if (ooo->ident == bestcont->container) {
            //  And the container is here too!
            frg->parent = bestcont->container;

            //  The hangs assume the container is forward; adjust if not so.
            if (ooo->position.bgn < ooo->position.end) {
              frg->ahang  = bestcont->a_hang;
              frg->bhang  = bestcont->b_hang;
            } else {
              frg->ahang  = -bestcont->b_hang;
              frg->bhang  = -bestcont->a_hang;
            }
          }

        } else {
          //  Frag is not contained.
          if ((bestedge5) && (ooo->ident == bestedge5->frag_b_id)) {
            //  Our fragment is forward.  Reverse the overlap.
            //
            //  reality  A -----------      bestedge  B    -------->
            //                B --------->            A -------
            //  
            frg->parent = bestedge5->frag_b_id;
            frg->ahang  = -bestedge5->ahang;
            frg->bhang  = -bestedge5->bhang;

            assert(frg->ahang >= 0);
          }
          if ((bestedge3) && (ooo->ident == bestedge3->frag_b_id)) {
            //  Our fragment is backward.  Reverse the overlap.
            //
            //  reality  A ---------       bestedge  B ------->
            //                 B <-------            A    --------->
            //
            frg->parent = bestedge3->frag_b_id;
            frg->ahang  = bestedge3->bhang;
            frg->bhang  = bestedge3->ahang;

            assert(frg->ahang >= 0);
          }
        }
      }
    }
  }
}




void
UnitigGraph::populateUnitig(Unitig           *unitig,
                            BestEdgeOverlap  *bestnext,
                            bool              verbose){

  assert(unitig->getLength() > 0);

  if ((bestnext == NULL) || (bestnext->frag_b_id == 0))
    //  Nothing to add!
    return;

  DoveTailNode frag = unitig->dovetail_path_ptr->back();

  //  The ID of the last fragment in the unitig, and the end we should walk off of it.
  int32 lastID  = frag.ident;
  int32 lastEnd = (frag.position.bgn < frag.position.end) ? THREE_PRIME : FIVE_PRIME;

                   
  if (Unitig::fragIn(bestnext->frag_b_id) == unitig->id())
    //  Cicrular unitig.  Deal with later.
    return;

  if (Unitig::fragIn(bestnext->frag_b_id) != 0) {
    //  Intersection.  Remember.
    fprintf(stderr,"unitigIntersect: unitig %5d frag %7d -> unitig %5d frag %7d (before construction)\n",
            unitig->id(), lastID, Unitig::fragIn(bestnext->frag_b_id), bestnext->frag_b_id);
    unitigIntersect[bestnext->frag_b_id].push_back(lastID);
    return;
  }

  //  While there are fragments to add, construct a reverse-edge, and add the fragment.

  while (bestnext->frag_b_id != 0) {
    BestEdgeOverlap  bestprev;

    //  Reverse nextedge (points from the unitig to the next fragment to add) so that it points from
    //  the next fragment to add back to something in the unitig.  If the fragments are
    //  innie/outtie, we need to reverse the overlap to maintain that the A fragment is forward.

    bestprev.frag_b_id = lastID;
    bestprev.bend      = lastEnd;
    bestprev.ahang     = -bestnext->ahang;
    bestprev.bhang     = -bestnext->bhang;

    if (bestprev.bend == bestnext->bend) {
      bestprev.ahang = bestnext->bhang;
      bestprev.bhang = bestnext->ahang;
    }

    //  Call the usual placement routine to place the next fragment relative to the last one.  This
    //  call depends on which end of the frag-to-be-added we are working with.

    frag.ident = bestnext->frag_b_id;

    int32  bidx5 = -1, bidx3 = -1;

    if (unitig->placeFrag(frag, bidx5, (bestnext->bend == THREE_PRIME) ? NULL : &bestprev,
                          frag, bidx3, (bestnext->bend == THREE_PRIME) ? &bestprev : NULL)) {
      unitig->addFrag(frag, 0, true);

    } else {

      fprintf(stderr, "ERROR:  Failed to place frag %d into BOG path.\n", frag.ident);
      assert(0);
    }

    //  Set up for the next fragmnet

    lastID  = frag.ident;
    lastEnd = (frag.position.bgn < frag.position.end) ? THREE_PRIME : FIVE_PRIME;

    bestnext = bog_ptr->getBestEdgeOverlap(lastID, lastEnd);

    //  Abort if we intersect, or are circular

    if (Unitig::fragIn(bestnext->frag_b_id) == unitig->id()) {
      //fprintf(stderr,"unitigIntersect: unitig %5d frag %7d -> unitig %5d frag %7d\n",
      //        unitig->id(), lastID, Unitig::fragIn(bestnext->frag_b_id), bestnext->frag_b_id);
      //unitigIntersect[bestnext->frag_b_id].push_back(lastID);
      //bestnext->frag_b_id = 0;BAD
      break;
    }

    if (Unitig::fragIn(bestnext->frag_b_id) != 0) {
      fprintf(stderr,"unitigIntersect: unitig %5d frag %7d -> unitig %5d frag %7d (during construction)\n",
              unitig->id(), lastID, Unitig::fragIn(bestnext->frag_b_id), bestnext->frag_b_id);
      unitigIntersect[bestnext->frag_b_id].push_back(lastID);
      //bestnext->frag_b_id = 0;
      break;
    }
  }
}




void
UnitigGraph::populateUnitig(int32 frag_idx,
                            bool  verbose) {

  Unitig *utg = new Unitig(verbose);

  unitigs->push_back(utg);

  //  Add a first fragment -- to be 'compatable' with the old code, the first fragment is added
  //  reversed, we walk off of its 5' end, flip it, and add the 3' walk.

  DoveTailNode  frag;

  frag.type         = AS_READ;
  frag.ident        = frag_idx;
  frag.contained    = 0;
  frag.parent       = 0;
  frag.ahang        = 0;
  frag.bhang        = 0;
  frag.position.bgn = _fi->fragmentLength(frag_idx);
  frag.position.end = 0;
  frag.delta_length = 0;
  frag.delta        = NULL;

  utg->addFrag(frag, 0, true);

  //  Add fragments as long as there is a path to follow...from the 3' end of the first fragment.

  BestEdgeOverlap  *bestedge5 = bog_ptr->getBestEdgeOverlap(frag_idx, FIVE_PRIME);
  BestEdgeOverlap  *bestedge3 = bog_ptr->getBestEdgeOverlap(frag_idx, THREE_PRIME);

  fprintf(stderr, "Adding 5' edges off of fragment %d in unitig %d\n",
          utg->dovetail_path_ptr->back().ident, utg->id());

  if (bestedge5->frag_b_id)
    populateUnitig(utg, bestedge5, verbose);

  utg->reverseComplement(false);

  fprintf(stderr, "Adding 3' edges off of fragment %d in unitig %d\n",
          utg->dovetail_path_ptr->back().ident, utg->id());

  if (bestedge3->frag_b_id)
    populateUnitig(utg, bestedge3, verbose);

  utg->reverseComplement(false);
}




void UnitigGraph::breakUnitigs(ContainerMap &cMap, char *output_prefix) {
  FILE *breakFile;

  {
    char name[FILENAME_MAX];
    sprintf(name, "%s.breaks.ovl", output_prefix);

    errno = 0;
    breakFile = fopen(name, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' to write unitig breaks.\n", name), exit(1);
  }

  fprintf(stderr, "==> BREAKING UNITIGS.\n");

  //  Debug output
  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig        *tig   = (*unitigs)[ti];

    if (tig == NULL)
      continue;

    DoveTailNode   first = tig->dovetail_path_ptr->front();
    uint32         prev  = 0;
    DoveTailNode   last  = tig->getLastBackboneNode(prev);

    if ( prev == 0 )
      continue; // skip singletons

    uint32 id1 = first.ident;
    uint32 id2 = last.ident;

    assert(prev != id2);
    assert(last.contained == 0);

    BestEdgeOverlap *firstBest = bog_ptr->getBestEdgeOverlap(id1, (first.position.bgn < first.position.end) ? FIVE_PRIME  : THREE_PRIME);
    BestEdgeOverlap *lastBest  = bog_ptr->getBestEdgeOverlap(id2, (last.position.bgn  < last.position.end)  ? THREE_PRIME : FIVE_PRIME);

    // Skip non bubbles, anchored ends
    if (lastBest->frag_b_id == 0 || firstBest->frag_b_id == 0)
      continue;

    if (firstBest->frag_b_id == lastBest->frag_b_id)
      fprintf(stderr, "Bubble %d to %d len %d (self bubble %d to %d)\n", firstBest->frag_b_id, lastBest->frag_b_id, tig->getLength(), id1, id2);
    else
      fprintf(stderr, "Bubble %d to %d len %d\n", firstBest->frag_b_id, lastBest->frag_b_id, tig->getLength());
  }


  //  Stop when we've seen all current unitigs.  Replace tiMax
  //  in the for loop below with unitigs->size() to recursively
  //  split unitigs.
  int  tiMax = unitigs->size();

  for (int  ti=0; ti<tiMax; ti++) {
    Unitig             *tig = (*unitigs)[ti];

    if (tig == NULL)
      continue;

    UnitigBreakPoints   breaks;
    DoveTailNode        lastBackbone;

    int                 numFragsInUnitig = 0;
    int                 fragCount        = 0;
    int                 fragIdx;

    tig->sort();

    //  Count the number of fragments in this unitig, including
    //  yet-to-be-placed contained fragments.
    for (fragIdx=0; fragIdx<tig->dovetail_path_ptr->size(); fragIdx++) {
      DoveTailNode  *f = &(*tig->dovetail_path_ptr)[fragIdx];

      if (cMap.find(f->ident) != cMap.end())
        numFragsInUnitig += cMap[f->ident].size();

      numFragsInUnitig++;

      //  Contained fragments should not be placed yet.
      assert(f->contained == 0);
    }


    for (fragIdx=0; fragIdx<tig->dovetail_path_ptr->size(); fragIdx++) {
      DoveTailNode  *f = &(*tig->dovetail_path_ptr)[fragIdx];

      fragCount++;
      if (cMap.find(f->ident) != cMap.end())
        fragCount += cMap[f->ident].size();

#ifdef VERBOSEBREAK
      if (fragIdx == 0) {
        uint32             dtEnd = (isReverse(f->position)) ? THREE_PRIME : FIVE_PRIME;
        BestEdgeOverlap   *bEdge = bog_ptr->getBestEdgeOverlap(f->ident, dtEnd);

        if ((bEdge) && (bEdge->frag_b_id > 0))
          fprintf(stderr,"unitig %d %c' frag %d points to unitig %d frag %d\n",
                  tig->id(), (dtEnd == THREE_PRIME) ? '3' : '5', f->ident,
                  Unitig::fragIn(bEdge->frag_b_id), bEdge->frag_b_id);
      }

      if (fragIdx + 1 == tig->dovetail_path_ptr->size()) {
        uint32             dtEnd = (isReverse(f->position)) ? FIVE_PRIME : THREE_PRIME;
        BestEdgeOverlap   *bEdge = bog_ptr->getBestEdgeOverlap(f->ident, dtEnd);

        if ((bEdge) && (bEdge->frag_b_id > 0))
          fprintf(stderr,"unitig %d %c' frag %d points to unitig %d frag %d\n",
                  tig->id(), (dtEnd == THREE_PRIME) ? '3' : '5', f->ident,
                  Unitig::fragIn(bEdge->frag_b_id), bEdge->frag_b_id);
      }
#endif

      FragmentEdgeList::const_iterator edge_itr = unitigIntersect.find(f->ident);

      if (edge_itr == unitigIntersect.end())
        continue;

      // We have a set of best edges incoming from other unitigs
      for (FragmentList::const_iterator fragItr = edge_itr->second.begin();
           fragItr != edge_itr->second.end();
           fragItr++) {
        uint32 inFrag = *fragItr;

        // check if it's incoming frag's 5' best edge.  If not, it must be the 3' edge.
        uint32            bestEnd  = FIVE_PRIME;
        BestEdgeOverlap  *bestEdge = bog_ptr->getBestEdgeOverlap(inFrag, bestEnd);

        if (bestEdge->frag_b_id != f->ident) {
          bestEnd  = THREE_PRIME;
          bestEdge = bog_ptr->getBestEdgeOverlap(inFrag, bestEnd);
          assert(bestEdge->frag_b_id == f->ident);
        }

        int pos = (bestEdge->bend == FIVE_PRIME) ? f->position.bgn : f->position.end;

        Unitig *inTig = (*unitigs)[Unitig::fragIn(inFrag)];

        assert(inTig->id() == Unitig::fragIn(inFrag));

        if (inTig) {
          UnitigBreakPoint breakPoint(f->ident, bestEdge->bend);

          breakPoint.fragPos     = f->position;
          breakPoint.fragsBefore = fragCount;
          breakPoint.fragsAfter  = numFragsInUnitig - fragCount;
          breakPoint.inSize      = inTig->getLength();
          breakPoint.inFrags     = inTig->getNumFrags();

          breaks.push_back(breakPoint);
        } else {
          //fprintf(stderr, "  NullBreak tig %5d at frag %7d\n", tig->id(), f->ident);
        }

#ifdef VERBOSEBREAK
        fprintf(stderr, "unitig %d (%d frags, len %d) frag %d end %c' into unitig %d frag %d end %c' pos %d\n",
                Unitig::fragIn(inFrag),
                inTig->getNumFrags(),
                inTig->getLength(),
                inFrag,
                (bestEnd == FIVE_PRIME) ? '5' : '3',
                tig->id(),
                f->ident,
                (bestEdge->bend == FIVE_PRIME) ? '5' : '3',
                pos);
#endif
        {
          GenericMesg  pmesg;
          OverlapMesg  omesg;

          omesg.aifrag          = inFrag;
          omesg.bifrag          = f->ident;
          omesg.ahg             = bestEdge->ahang;
          omesg.bhg             = bestEdge->bhang;
          omesg.orientation     = AS_UNKNOWN;
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          if ((bestEnd == FIVE_PRIME) && (bestEdge->bend == FIVE_PRIME))
            omesg.orientation = AS_OUTTIE;
          if ((bestEnd == FIVE_PRIME) && (bestEdge->bend == THREE_PRIME))
            omesg.orientation = AS_ANTI;
          if ((bestEnd == THREE_PRIME) && (bestEdge->bend == FIVE_PRIME))
            omesg.orientation = AS_NORMAL;
          if ((bestEnd == THREE_PRIME) && (bestEdge->bend == THREE_PRIME))
            omesg.orientation = AS_INNIE;

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(breakFile, &pmesg);
        }
      }  //  Over all fragments that intersect with us
    }  //  Over all fragments in the unitig

    if (breaks.empty() == false) {
      DoveTailNode  *f = &tig->dovetail_path_ptr->back();

      //  create a final fake bp for the last frag so we
      //  have a reference point to the end of the tig for
      //  filtering the breakpoints.  fakeEnd, for searching.
      //
      UnitigBreakPoint breakPoint(f->ident, (isReverse(f->position)) ? FIVE_PRIME : THREE_PRIME);

      //  +1 below seems like a bug, but it reproduces what was here
      //  before.
#warning possible bug

      breakPoint.fragPos     = f->position;
      breakPoint.fragsBefore = fragCount + 1;
      breakPoint.fragsAfter  = 0;
      breakPoint.inSize      = std::numeric_limits<int>::max();
      breakPoint.inFrags     = 0;

      breaks.push_back(breakPoint);

      UnitigVector* newUs = breakUnitigAt(cMap, tig, breaks);

      if (newUs != NULL) {
        delete tig;
        (*unitigs)[ti] = NULL;
        unitigs->insert(unitigs->end(), newUs->begin(), newUs->end());
      }

      delete newUs;
    }
  }  //  Over all tigs

  fclose(breakFile);
}



void
UnitigGraph::placeContains(void) {
  uint32   fragsPlaced  = 1;
  uint32   fragsPending = 0;

  while (fragsPlaced > 0) {
    fragsPlaced  = 0;
    fragsPending = 0;

    fprintf(stderr, "==> PLACING CONTAINED FRAGMENTS\n");

    for (uint32 fid=0; fid<_fi->numFragments()+1; fid++) {
      BestContainment *bestcont = bog_ptr->getBestContainer(fid);
      Unitig          *utg;

      if (bestcont == NULL)
        //  Not a contained fragment.
        continue;

      if (bestcont->isPlaced == true)
        //  Containee already placed.
        continue;

      if (Unitig::fragIn(bestcont->container) == 0) {
        //  Container not placed (yet).
        fragsPending++;
        continue;
      }

      utg = (*unitigs)[Unitig::fragIn(bestcont->container)];
      utg->addContainedFrag(fid, bestcont, true);

      bestcont->isPlaced = true;

      fragsPlaced++;
    }

    fprintf(stderr, "==> PLACING CONTAINED FRAGMENTS - placed %d fragments; still need to place %d\n",
            fragsPlaced, fragsPending);

    if ((fragsPlaced == 0) && (fragsPending > 0)) {
      fprintf(stderr, "Stopping contained fragment placement due to zombies.\n");
      fragsPlaced  = 0;
      fragsPending = 0;
    }
  }

  for (int ti=0; ti<unitigs->size(); ti++) {
    Unitig *utg = (*unitigs)[ti];

    if (utg)
      utg->sort();
  }
}



//  This is a huge hack to get around a bug somewhere before us.  We
//  seem to not be placing all fragments.  So, run through all
//  fragments, and for anything not placed, toss it into a new
//  unitig.  Let scaffolder figure out where to put it.
//
//  Notice the similarity between this code and checkUnitigMembership().
//
void
UnitigGraph::placeZombies(void) {

  fprintf(stderr, "==> SEARCHING FOR ZOMBIES\n");

  uint32 *inUnitig   = new uint32 [_fi->numFragments()+1];
  int     numZombies = 0;

  //  Sorry, poor fella that has more than 987,654,321 unitigs.  I
  //  can't imagine the rest of the pipeline would run though.
  //
  //  Mark fragments as dead.
  //
  for (uint32 i=0; i<_fi->numFragments()+1; i++)
    inUnitig[i] = 987654321;

  //  ZZZzzzzaapaapppp!  IT'S ALIVE!
  //
  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg)
      for (DoveTailIter it=utg->dovetail_path_ptr->begin(); it != utg->dovetail_path_ptr->end(); it++)
        inUnitig[it->ident] = utg->id();
  }

  //  Anything still dead?
  //
  for (uint32 i=0; i<_fi->numFragments()+1; i++) {
    if (_fi->fragmentLength(i) > 0) {
      if (inUnitig[i] == 0) {
        //  We'll catch this error inna second in checkUnitigMembership().
      } else if (inUnitig[i] != 987654321) {
        //  We'll count this inna second there too.
      } else {
        //  Ha!  Gotcha!  You're now a resurrected brain eating
        //  zomibie?  Some day we'll figure out how to put you in
        //  properly.  For now, enjoy the ride.

        Unitig *utg = new Unitig(false);

        DoveTailNode frag;

        frag.type         = AS_READ;
        frag.ident        = i;
        frag.contained    = 0;
        frag.delta_length = 0;
        frag.delta        = NULL;

        frag.position.bgn = 0;
        frag.position.end = _fi->fragmentLength(i);

        utg->addFrag(frag, 0, false);
        unitigs->push_back(utg);

        numZombies++;
      }
    }
  }

  if (numZombies > 0)
    fprintf(stderr, "RESURRECTED %d ZOMBIE FRAGMENTS.\n", numZombies);

  delete inUnitig;
}


//  Sometimes, we don't include fragments into a unitig, even though the best overlap off of both
//  ends is to the same unitig.  This will include those.
//
void
UnitigGraph::popBubbles(void) {

  fprintf(stderr, "==> SEARCHING FOR BUBBLES\n");

  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig        *utg = (*unitigs)[ti];

    if ((utg == NULL) ||
        (utg->dovetail_path_ptr == NULL) ||
        (utg->dovetail_path_ptr->size() >= 30))
      continue;

    uint32         otherUtg  = 987654321;
    uint32         conflicts = 0;
    uint32         self      = 0;
    uint32         nonmated  = 0;
    uint32         matedcont = 0;
    uint32         spurs     = 0;

    for (int fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode *frg = &(*utg->dovetail_path_ptr)[fi];

      int32  frgID = frg->ident;
      int32  utgID = utg->id();

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);
      BestContainment *bestcont  = bog_ptr->getBestContainer(frg->ident);

      if (_fi->mateIID(frgID) > 0) {
        if (bestcont)
          matedcont++;
      } else {
        nonmated++;
      }

      if (bestcont)
        continue;

      if (bestedge5->frag_b_id == 0) {
        spurs++;
      } else {
        int32 ou5 = utg->fragIn(bestedge5->frag_b_id);

        assert(ou5 > 0);

        if (ou5 == utg->id()) {
          self++;
        } else {
          if ((otherUtg == 987654321) && (ou5 != 0))
            otherUtg = ou5;
          if (otherUtg != ou5)
            conflicts++;
        }
      }

      if (bestedge3->frag_b_id == 0) {
        spurs++;
      } else {
        int32 ou3 = utg->fragIn(bestedge3->frag_b_id);

        assert(ou3 > 0);

        if (ou3 == utg->id()) {
          self++;
        } else {
          if ((otherUtg == 987654321) && (ou3 != 0))
            otherUtg = ou3;
          if (otherUtg != ou3)
            conflicts++;
        }
      }
    }  //  Over all frags

    fprintf(stderr, "CONFLICTS %d SPURS %d SELF %d len %d frags %d matedcont %d nonmated %d\n",
            conflicts, spurs, self, utg->getLength(), utg->dovetail_path_ptr->size(), matedcont, nonmated);

    if ((conflicts > 0) ||
        (spurs     > 0) ||
        (self      > 6) ||
        (matedcont > 6) ||
        (otherUtg == 987654321))
      continue;

    //  Does the new length make sense?

    Unitig  *mergeTig    = (*unitigs)[otherUtg];

#if 0
    int32  minU = 0;
    int32  maxU = 0;

    for (int fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode *frg = &(*utg->dovetail_path_ptr)[fi];

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);

      DoveTailNode place5 = {0};
      DoveTailNode place3 = {0};

      mergeTig->placeFrag(place5, place3, bestedge5, bestedge3);

      int32  minF5 = 0, minF3 = 0;
      int32  maxF5 = 0, maxF3 = 0;

      if (place5.position.bgn != place5.position.end) {
        minF5 = MIN(place5.position.bgn, place5.position.end);
        maxF5 = MAX(place5.position.bgn, place5.position.end);
      }

      if (place3.position.bgn != place3.position.end) {
        minF3 = MIN(place3.position.bgn, place3.position.end);
        maxF3 = MAX(place3.position.bgn, place3.position.end);
      }
    }
#endif

    //  Merge this unitig into otherUtg.

    fprintf(stderr, "MERGE unitig %d (len %d) into unitig %d (len %d)\n",
            utg->id(), utg->getLength(),
            mergeTig->id(), mergeTig->getLength());

    for (int fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode  *frag = &(*utg->dovetail_path_ptr)[fi];

      if (bog_ptr->getBestContainer(frag->ident))
        mergeTig->addContainedFrag(frag->ident,
                                   bog_ptr->getBestContainer(frag->ident),
                                   true);
      else
        mergeTig->addAndPlaceFrag(frag->ident,
                                  bog_ptr->getBestEdgeOverlap(frag->ident, FIVE_PRIME),
                                  bog_ptr->getBestEdgeOverlap(frag->ident, THREE_PRIME),
                                  true);
    }

    mergeTig->sort();

    //  Mark this unitig as dead.

    (*unitigs)[ti] = NULL;
    delete utg;
  }  //  over all unitigs
}  //  orphan searching scope





float UnitigGraph::getGlobalArrivalRate(long total_random_frags_in_genome, long genome_size){

  float _globalArrivalRate;

  //  If the genome size has not been specified, estimate the GAR.
  if(genome_size == 0){

    float total_rho=0, avg_rho;
    float total_arrival_frags=0;
    size_t rho_gt_10000 = 0;

    // Go through all the unitigs to sum rho and unitig arrival frags
    UnitigVector::const_iterator iter;
    for(
        iter=unitigs->begin();
        iter!=unitigs->end();
        iter++){

      if (*iter == NULL)
        continue;

      avg_rho = (*iter)->getAvgRho(_fi);
      total_rho += avg_rho;
      if (avg_rho > 10000.0)
        rho_gt_10000 += (size_t)avg_rho / 10000;

      float unitig_random_frags = (*iter)->getNumRandomFrags();
      if (--unitig_random_frags < 0)
        unitig_random_frags = 0;

      total_arrival_frags += unitig_random_frags;
      (*iter)->setLocalArrivalRate( unitig_random_frags / avg_rho );
    }
    // Estimate GAR
    _globalArrivalRate = (total_rho > 0) ? (total_arrival_frags / total_rho): 0;

    std::cerr << "Calculated Global Arrival rate " << _globalArrivalRate <<
      std::endl;
    // Now recalculate based on big unitigs, copied from AS_CGB/AS_CGB_cgb.c
    if (rho_gt_10000 * 20000 > total_rho) {
      float min_10_local_arrival_rate          = _globalArrivalRate;
      float median_local_arrival_rate          = _globalArrivalRate;
      float max_local_arrival_rate             = _globalArrivalRate;
      float recalibrated_fragment_arrival_rate = _globalArrivalRate;
      size_t num_arrival_rates=0;
      int median_index;
      std::vector<float> arrival_rate_array(rho_gt_10000);

      for( iter=unitigs->begin(); iter!=unitigs->end(); iter++) {
        if (*iter == NULL)
          continue;
        avg_rho = (*iter)->getAvgRho(_fi);
        if (avg_rho > 10000.0) {
          const int num_10000 = (size_t)avg_rho / 10000;
          const float local_arrival_rate =
            (*iter)->getNumRandomFrags() / avg_rho;
          assert(num_10000 > 0);
          int i;
          for(i=0;i<num_10000;i++){
            assert(i < rho_gt_10000);
            arrival_rate_array[i] = local_arrival_rate;
            num_arrival_rates++;
          }
          if(num_arrival_rates > 0){
            float tmp_fragment_arrival_rate, max_diff_arrival_rate;
            float prev_arrival_rate, cur_arrival_rate, diff_arrival_rate;
            int max_diff_index;
            std::sort(arrival_rate_array.begin(),arrival_rate_array.end());
            min_10_local_arrival_rate = arrival_rate_array[num_arrival_rates / 10];
            median_index = (num_arrival_rates * 5) / 10;
            median_local_arrival_rate = arrival_rate_array[median_index];
            max_local_arrival_rate = arrival_rate_array[num_arrival_rates-1];
            recalibrated_fragment_arrival_rate =
              arrival_rate_array[(num_arrival_rates * 19) / 20];
            prev_arrival_rate = min_10_local_arrival_rate;
            max_diff_arrival_rate = 0.0;
            for(i=num_arrival_rates / 10;i<median_index;i++){
              cur_arrival_rate = arrival_rate_array[i];
              diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
              prev_arrival_rate = cur_arrival_rate;
              if(diff_arrival_rate > max_diff_arrival_rate){
                max_diff_arrival_rate = diff_arrival_rate;
              }
            }
            max_diff_arrival_rate *= 2.0;
            max_diff_index = num_arrival_rates - 1;
            for(i=median_index;i<num_arrival_rates;i++){
              cur_arrival_rate = arrival_rate_array[i];
              diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
              prev_arrival_rate = cur_arrival_rate;
              if(diff_arrival_rate > max_diff_arrival_rate){
                max_diff_arrival_rate = diff_arrival_rate;
                max_diff_index = i - 1;
                break;
              }
            }
            max_diff_arrival_rate = arrival_rate_array[max_diff_index];
            tmp_fragment_arrival_rate =MIN(min_10_local_arrival_rate * 2.0,
                                           median_local_arrival_rate * 1.25);
            if(tmp_fragment_arrival_rate < recalibrated_fragment_arrival_rate){
              recalibrated_fragment_arrival_rate = tmp_fragment_arrival_rate;
            }
            if(max_diff_arrival_rate < recalibrated_fragment_arrival_rate){
              recalibrated_fragment_arrival_rate = max_diff_arrival_rate;
            }
          }
          if(recalibrated_fragment_arrival_rate > _globalArrivalRate){
            _globalArrivalRate = recalibrated_fragment_arrival_rate;
            std::cerr <<
              "Used recalibrated global_fragment_arrival_rate="
                      << _globalArrivalRate << std::endl
                      << "Used recalibrated global_fragment_arrival_distance="
                      <<((_globalArrivalRate > 0.) ? 1./(_globalArrivalRate) : 0.)
                      << std::endl
                      << "Chunk arrival rates sorted at 1/100s"
                      << std::endl ;
            for(i=0;i<100;i++) {
              std::cerr<<arrival_rate_array[((num_arrival_rates * i) / 100)]
                       << std::endl;
            }
            std::cerr << max_local_arrival_rate << std::endl;
          }
        }
      }
    }

   std::cerr << 
      "Computed genome_size="
      << (_globalArrivalRate > 0.f ? ((float)total_random_frags_in_genome / _globalArrivalRate) : 0.f) 
      << std::endl;
  }else{
    // Compute actual GAR
    _globalArrivalRate = (float)total_random_frags_in_genome / (float)genome_size;
  }

  return(_globalArrivalRate);

}






UnitigBreakPoint UnitigGraph::selectSmall(ContainerMap &cMap,
                                          const Unitig *tig,
                                          const UnitigBreakPoints &smalls,
                                          const UnitigBreakPoint& big,
                                          int   &lastBPCoord,
                                          int   &lastBPFragNum) {
  UnitigBreakPoint selection;

  double difference = 0.0;
  int  rFrgs   = big.fragsBefore;
  bool bRev    = isReverse(big.fragPos);
  int bContain = 0;

  int right    = (big.fragEnd.fragEnd() == FIVE_PRIME) ? big.fragPos.bgn : big.fragPos.end;

  if (cMap.find(big.fragEnd.fragId()) != cMap.end())
    bContain = cMap[big.fragEnd.fragId()].size();

  if (((bRev == true)  && (big.fragEnd == THREE_PRIME)) ||
      ((bRev == false) && (big.fragEnd == FIVE_PRIME)))
    rFrgs -= 1 + bContain;

  for(UnitigBreakPoints::const_iterator sIter = smalls.begin(); sIter != smalls.end(); sIter++) {
    UnitigBreakPoint small = *sIter;

    if (small.fragEnd.fragId() == big.fragEnd.fragId())
      continue;

    bool rev     = isReverse(small.fragPos);
    int sContain = 0;
    int lFrgs    = small.fragsBefore;
    uint32 sid     = small.fragEnd.fragId();

    if (cMap.find(sid) != cMap.end())
      sContain = cMap[sid].size();

    // left side of the frag in the unitig, don't count it
    if (((rev == true)  && (small.fragEnd == THREE_PRIME)) ||
        ((rev == false) && (small.fragEnd == FIVE_PRIME)))
      lFrgs -= 1 + sContain;

    if (rFrgs - lFrgs == 1)
      continue;

    // use middle of frag instead of end, to get some overlap
    double bp    = (small.fragPos.end + small.fragPos.bgn) / 2.0;

    double lRate = (lFrgs - lastBPFragNum) / (bp - lastBPCoord);
    double rRate = (rFrgs - lFrgs) / (right - bp);
    double ratio = (lRate > rRate) ? lRate / rRate : rRate / lRate;
    double diff  = fabs( lRate - rRate);

    if ((ratio > 1.8) && (diff > difference)) {
      //fprintf(stderr, "Break frg %7d b %4d l %4d pos b %5d e %5.0f lRate %.4f\n", sid, lastBPFragNum, lFrgs, lastBPCoord, bp, lRate );
      //fprintf(stderr, "     diff %4d r %4d pos %5d rRate %.4f ratio %.2f to frag %7d\n", rFrgs - lFrgs, rFrgs, right, rRate, ratio, big.fragEnd.fragId());
      //fprintf(stderr,"     select frg %d for break on arrival rate diff %.4f\n", sid, diff);
      difference = diff;
      selection = small;
    }
  }

  lastBPCoord   = right;
  lastBPFragNum = rFrgs;

  return selection;
}



void UnitigGraph::filterBreakPoints(ContainerMap &cMap,
                                    Unitig *tig,
                                    UnitigBreakPoints &breaks) {

  int lastBPCoord = 0;
  int lastBPFragNum = 0;

  // Check for sentinal at end, not added by MateChecker
  UnitigBreakPoint fakeEnd = breaks.back();
  if (fakeEnd.inSize == std::numeric_limits<int>::max())
    breaks.pop_back();

  UnitigBreakPoints smallBPs;
  UnitigBreakPoints newBPs;

  for(UnitigBreakPoints::iterator iter = breaks.begin(); iter != breaks.end(); iter++) {
    UnitigBreakPoint nextBP = *iter;

    if ((nextBP.inFrags > UnitigGraph::MIN_BREAK_FRAGS) && (nextBP.inSize > UnitigGraph::MIN_BREAK_LENGTH)) {

      // big one, compare against smalls -- if we haven't
      // seen this break point already.
      //
      if (newBPs.empty() || (nextBP.fragEnd != newBPs.back().fragEnd)) {

        if (smallBPs.empty() == false) {
          //  smalls exist, select one.
          UnitigBreakPoint small = selectSmall(cMap, tig, smallBPs, nextBP, lastBPCoord, lastBPFragNum);
          if (small.fragEnd.fragId() > 0)
            newBPs.push_back(small);
          smallBPs.clear();

        } else {
          //  No smalls.  Update state to move past the current big breakpoint.
          lastBPCoord   = (nextBP.fragEnd.fragEnd() == FIVE_PRIME) ? nextBP.fragPos.bgn : nextBP.fragPos.end;
          lastBPFragNum = nextBP.fragsBefore;

          int bContain = 0;

          if (cMap.find(nextBP.fragEnd.fragId()) != cMap.end())
            bContain = cMap[nextBP.fragEnd.fragId()].size();

          bool bRev = isReverse(nextBP.fragPos);

          if (((bRev == true)  && (nextBP.fragEnd == THREE_PRIME)) ||
              ((bRev == false) && (nextBP.fragEnd == FIVE_PRIME)))
            lastBPFragNum -= 1 + bContain;
        }

        //  Haven't seen this break, so add it.
        newBPs.push_back( nextBP );
      }
    } else {
      //  Not a big breaker.  Save the small breaker if we've not seen it yet.
      if (newBPs.empty() || ((nextBP.fragEnd != newBPs.back().fragEnd) &&
                             (nextBP.fragEnd != smallBPs.back().fragEnd)))
        smallBPs.push_back(nextBP);
    }
  }

  //  If we've got small ones saved, select one.

  if (smallBPs.empty() == false) {
    UnitigBreakPoint small = selectSmall(cMap, tig, smallBPs, fakeEnd, lastBPCoord, lastBPFragNum);
    if (small.fragEnd.fragId() > 0)
      newBPs.push_back(small);
    smallBPs.clear();
  }

  breaks = newBPs;
}

// Doesn't handle contained frags yet
UnitigVector* UnitigGraph::breakUnitigAt(ContainerMap &cMap,
                                         Unitig *tig,
                                         UnitigBreakPoints &breaks) {

  if (breaks.empty())
    return NULL;

  // remove small break points
  filterBreakPoints(cMap, tig, breaks);

  // if we filtered all the breaks out return an empty list
  if (breaks.empty())
    return NULL;

  //  This is to catch a bug in....something before us.
  //  Sometimes the breakpoint list contains nonsense that tells
  //  us to split a unitig before the first fragment, or after
  //  the last.
  //
  //  That code is pretty dense (and big) so instead we'll make
  //  one pass through the list of breaks and see if we'd do any
  //  splitting.
  //
  //  It's a huge code duplication of the main loop.
  //
  {
    UnitigBreakPoints breakstmp = breaks;

    UnitigBreakPoint  breakPoint = breakstmp.front();
    breakstmp.pop_front();

    UnitigBreakPoint nextBP;

    int   newUnitigsConstructed = 0;
    int   newUnitigExists = 0;

    if (!breakstmp.empty())
      nextBP = breakstmp.front();

    for(DoveTailIter dtIter = tig->dovetail_path_ptr->begin(); dtIter != tig->dovetail_path_ptr->end(); dtIter++) {
      DoveTailNode frg = *dtIter;

      bool bothEnds = false;
      while ( !breakstmp.empty() && nextBP.fragEnd.fragId() == breakPoint.fragEnd.fragId() ) {
        if (nextBP.fragEnd.fragEnd() != breakPoint.fragEnd.fragEnd())
          bothEnds = true;
        breakstmp.pop_front();
        if (!breakstmp.empty())
          nextBP = breakstmp.front();
      }

      bool reverse = isReverse(frg.position);

      if (breakPoint.fragEnd.fragId() == frg.ident) {

        if (bothEnds) {
          newUnitigsConstructed++;
          newUnitigExists = 0;
        }

        else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && !reverse ||
                 breakPoint.fragEnd.fragEnd() == THREE_PRIME && reverse) {
          newUnitigsConstructed++;
          newUnitigExists = 1;
        }

        else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && reverse ||
                 breakPoint.fragEnd.fragEnd() == THREE_PRIME && !reverse ) {
          if (newUnitigExists == 0)
            newUnitigsConstructed++;
          newUnitigExists = 0;
        } else {
          assert(0);
        }

        if (breakstmp.empty()) {
          breakPoint.fragEnd = FragmentEnd();
        } else {
          breakPoint = breakstmp.front();
          breakstmp.pop_front();
          if (!breakstmp.empty())
            nextBP = breakstmp.front();
        }
      } else {
        if (newUnitigExists == 0) {
          newUnitigsConstructed++;
          newUnitigExists = 1;
        }
      }
    }

    if (newUnitigsConstructed < 2) {
      fprintf(stderr, "SPLITTING BUG DETECTED!  Adjusting for it.\n");
      return NULL;
    }

  }  //  end of explicit split test




  UnitigVector *splits = new UnitigVector();
  Unitig       *newTig = NULL;
  int           offset = 0;

  //  Emit a log of unitig creation and fragment moves
  //
  bool          verbose = false;

  UnitigBreakPoint breakPoint = breaks.front();
  breaks.pop_front();

  UnitigBreakPoint nextBP;

  if (!breaks.empty())
    nextBP = breaks.front();

  for(DoveTailIter dtIter = tig->dovetail_path_ptr->begin(); dtIter != tig->dovetail_path_ptr->end(); dtIter++) {
    DoveTailNode frg = *dtIter;

    // reduce multiple breaks at the same fragment end down to one
    bool bothEnds = false;
    while ( !breaks.empty() && nextBP.fragEnd.fragId() == breakPoint.fragEnd.fragId() ) {
      if (nextBP.fragEnd.fragEnd() != breakPoint.fragEnd.fragEnd())
        bothEnds = true;
      breaks.pop_front();
      if (!breaks.empty())
        nextBP = breaks.front();
    }

    bool reverse = isReverse(frg.position);

    if (breakPoint.fragEnd.fragId() == frg.ident) {

      //  There is method to this madness.  The three if
      //  blocks below have the same two basic code blocks
      //  in various orders.
      //
      //  o One code block makes a new unitig.
      //
      //  o The other adds a fragment to it, possibly
      //    remembering the beginning offset.
      //
      //  Finally, we delay making a new unitig until we are
      //  absolutely sure we're going to put fragments into
      //  it.

      //  Notes on fixing the break/join bug -- the bug is when we have a big breaker
      //  come in and chop off a small guy at the end.  We should be joining the two
      //  large unitigs.
      //
      //  To fix this, BPW was going to store enough stuff in the breakPoint so that
      //  the case can be detected here.  We need to know how much is before/after
      //  (depends on which end of the unitig we're at) our incoming point, (e.g.,
      //  fragsBefore and fragsAfter), and the fragments involved (fragEnd and inEnd).
      //
      //  When detected, push fragEnd and inEnd onto a list of joins.  After all breaks
      //  are finished, use the list of joins to get the unitigs to join and, well,
      //  join them.

      if (bothEnds) {
        //
        //  Break at both ends, create a singleton for this fragment.
        //
        if (verbose)
          fprintf(stderr,"  Break tig %d at both ends of %d num %d\n", tig->id(), breakPoint.fragEnd.fragId(), breakPoint.fragsBefore);

        newTig = new Unitig(verbose);  //  always make a new unitig, we put a frag in it now
        splits->push_back( newTig );

        if (newTig->dovetail_path_ptr->empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;
        newTig->addFrag( frg, offset, false );

        newTig = NULL;  //  delay until we need to make it
      }

      else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && !reverse ||
               breakPoint.fragEnd.fragEnd() == THREE_PRIME && reverse) {
        //
        //  Break at left end of frg, frg starts new tig
        //
        if (verbose)
          fprintf(stderr,"  Break tig %d before %d num %d\n", tig->id(), breakPoint.fragEnd.fragId(), breakPoint.fragsBefore);

        newTig = new Unitig(verbose);  //  always make a new unitig, we put a frag in it now
        splits->push_back( newTig );

        if (newTig->dovetail_path_ptr->empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;
        newTig->addFrag( frg, offset, false );
      }

      else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && reverse ||
               breakPoint.fragEnd.fragEnd() == THREE_PRIME && !reverse ) {
        //
        //  Break at right end of frg, frg goes in existing tig, then make new unitig
        //

        if (newTig == NULL) {  //  delayed creation?
          newTig = new Unitig(verbose);
          splits->push_back( newTig );
        }

        if (newTig->dovetail_path_ptr->empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;
        newTig->addFrag( frg, offset, false );

        if (verbose)
          fprintf(stderr,"  Break tig %d after %d num %d\n", tig->id(), breakPoint.fragEnd.fragId(), breakPoint.fragsBefore);

        //  Delay making a new unitig until we need it.
        newTig = NULL;
      } else {
        // logically impossible!
        assert(0);
      }


      // Done breaking, continue adding remaining frags to newTig
      if (breaks.empty()) {
        breakPoint.fragEnd = FragmentEnd();
      } else {
        breakPoint = breaks.front();
        breaks.pop_front();
        if (!breaks.empty())
          nextBP = breaks.front();
      }
    } else {
      //
      //  frag is not a unitig break point fragment, add to current new tig
      //

      if (newTig == NULL) {  //  delayed creation?
        newTig = new Unitig(verbose);
        splits->push_back( newTig );
      }
      if (newTig->dovetail_path_ptr->empty())
        offset = reverse ? -frg.position.end : -frg.position.bgn;
      newTig->addFrag( frg, offset, false );
    }
  }

  return splits;
}



void UnitigGraph::writeIUMtoFile(char *fileprefix, char *tigStorePath, int frg_count_target){
  int32       utg_count              = 0;
  int32       frg_count              = 0;
  int32       prt_count              = 1;
  char        filename[FILENAME_MAX] = {0};
  uint32     *partmap                = new uint32 [unitigs->size()];

  //  This code closely follows that in AS_CGB_unitigger.c::output_the_chunks()

  checkUnitigMembership();

  // Open up the initial output file

  sprintf(filename, "%s.iidmap", fileprefix);
  FILE *iidm = fopen(filename, "w");
  assert(NULL != iidm);

  sprintf(filename, "%s.partitioning", fileprefix);
  FILE *part = fopen(filename, "w");
  assert(NULL != part);

  sprintf(filename, "%s.partitioningInfo", fileprefix);
  FILE *pari = fopen(filename, "w");
  assert(NULL != pari);

  //  Step through all the unitigs once to build the partition mapping and IID mapping.

  for (uint32 iumiid=0, ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];
    uint32   nf  = (utg) ? utg->getNumFrags() : 0;

    if ((utg == NULL) || (nf == 0))
      continue;

    assert(utg->getLength() > 0);
    assert(nf == utg->dovetail_path_ptr->size());

    if ((0              <= frg_count_target) &&
        (frg_count + nf >= frg_count_target) &&
        (frg_count                      >  0)) {
      fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
              prt_count, utg_count, frg_count);

      prt_count++;
      utg_count = 0;
      frg_count = 0;
    }

    partmap[iumiid] = prt_count;

    fprintf(iidm, "Unitig "F_U32" == IUM "F_U32" (in partition "F_U32" with "F_S64" frags)\n",
            utg->id(), iumiid, partmap[iumiid], nf);

    for (int32 fragIdx=0; fragIdx<nf; fragIdx++) {
      DoveTailNode  *f = &(*utg->dovetail_path_ptr)[fragIdx];

      fprintf(part, "%d\t%d\n", prt_count, f->ident);

      //  We abused the delta_length field earlier.  Make sure it's sane.  If not, we assert when
      //  writing the tig.
      f->delta_length = 0;
      f->delta        = NULL;
    }

    utg_count += 1;
    frg_count += nf;

    iumiid++;
  }

  fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
          prt_count, utg_count, frg_count);

  fclose(pari);
  fclose(part);
  fclose(iidm);

  //  Step through all the unitigs once to build the partition mapping and IID mapping.

  MultiAlignStore  *MAS = new MultiAlignStore(tigStorePath);
  MultiAlignT      *ma  = CreateEmptyMultiAlignT();

  MAS->writeToPartitioned(partmap, NULL);

  for (uint32 iumiid=0, ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];
    uint32   nf  = (utg) ? utg->getNumFrags() : 0;

    if ((utg == NULL) || (nf == 0))
      continue;

    //  Massage the Unitig into a MultiAlignT (also used in SplitChunks_CGW.c)

    ma->maID                      = iumiid;
    ma->data.unitig_coverage_stat = utg->getCovStat(_fi);
    ma->data.unitig_microhet_prob = 1.0;  //  Default to 100% probability of unique

    ma->data.unitig_status        = AS_UNASSIGNED;
    ma->data.unitig_unique_rept   = AS_FORCED_NONE;

    ma->data.contig_status        = AS_UNPLACED;

    //  Add the fragments

    ResetVA_IntMultiPos(ma->f_list);
    SetRangeVA_IntMultiPos(ma->f_list, 0, nf, &(*utg->dovetail_path_ptr)[0]);

    //  NOTE!  This is not currently a valid multialign as it has NO IntUnitigPos.  That is
    //  added during consensus.  CGW will correctly assert that it reads in unitigs with
    //  exactly one IUP.

    //  Stash the unitig in the store

    MAS->insertMultiAlign(ma, TRUE, FALSE);

    iumiid++;
  }

  delete    MAS;
  delete [] partmap;
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
UnitigGraph::writeOVLtoFile(char *fileprefix) {
  char         filename[FILENAME_MAX] = {0};
  GenericMesg  pmesg;
  OverlapMesg  omesg;

  sprintf(filename, "%s.unused.ovl", fileprefix);
  FILE *file = fopen(filename, "w");
  assert(file != NULL);

  for (uint32  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode  *frg = &(*utg->dovetail_path_ptr)[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);

      int              bestident5 = 0;
      int              bestident3 = 0;

      if (bestedge5) {
        bestident5 = bestedge5->frag_b_id;

        if ((bestident5 > 0) && (utg->fragIn(bestident5) != utg->id())) {
          omesg.aifrag          = frg->ident;
          omesg.bifrag          = bestident5;
          omesg.ahg             = bestedge5->ahang;
          omesg.bhg             = bestedge5->bhang;
          omesg.orientation     = AS_UNKNOWN;
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          //  This overlap is off of the 5' end of this fragment.
          if (bestedge5->bend == FIVE_PRIME)
            omesg.orientation = AS_OUTTIE;
          if (bestedge5->bend == THREE_PRIME)
            omesg.orientation = AS_ANTI;

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(file, &pmesg);
        }
      }

      if (bestedge3) {
        bestident3 = bestedge3->frag_b_id;

        if ((bestident3 > 0) && (utg->fragIn(bestident3) != utg->id())) {
          omesg.aifrag          = frg->ident;
          omesg.bifrag          = bestident3;
          omesg.ahg             = bestedge3->ahang;
          omesg.bhg             = bestedge3->bhang;
          omesg.orientation     = AS_UNKNOWN;
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          //  This overlap is off of the 3' end of this fragment.
          if (bestedge3->bend == FIVE_PRIME)
            omesg.orientation = AS_NORMAL;
          if (bestedge3->bend == THREE_PRIME)
            omesg.orientation = AS_INNIE;

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(file, &pmesg);
        }
      }
    }
  }

  fclose(file);
}
