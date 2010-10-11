
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

static const char *rcsid = "$Id: AS_BOG_Unitig.cc,v 1.33 2010-10-11 03:43:44 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_Unitig.hh"
#include "AS_BOG_BestOverlapGraph.hh"

static std::map<uint32,int>* containPartialOrder;

uint32  Unitig::_nextId       = 1;
uint32* Unitig::_inUnitig     = NULL;
uint32* Unitig::_pathPosition = NULL;

float   Unitig::_globalArrivalRate = -1;


Unitig::Unitig(bool report){
  _localArrivalRate = -1;
  _length           =  0;
  _avgRho           = -1;
  _id               = _nextId++;

  if (report)
    fprintf(logFile, "Creating Unitig %d\n", _id);
}


Unitig::~Unitig(void){
}


#warning WHAT REALLLY HAPPENS IF NO BACKBONE NODE, OR NO PREVIOUS BACKBONE NODE

ufNode Unitig::getLastBackboneNode(void) {
  for (int32 fi=ufpath.size()-1; fi >= 0; fi--) {
    ufNode  &node = ufpath[fi];

    if (node.contained)
      continue;

    return(node);
  }

  fprintf(logFile, "Unitig::getLastBackboneNode()--  WARNING:  unitig %d has no backbone nodes, all contained!\n", id());
  ufNode last;
  memset(&last, 0, sizeof(ufNode));
  return(last);
}


ufNode Unitig::getLastBackboneNode(uint32 &prevID) {
  ufNode last;

  memset(&last, 0, sizeof(ufNode));

  prevID = 0;

  for (int32 fi=ufpath.size()-1; (fi >= 0) && (prevID == 0); fi--) {
    ufNode  *node = &ufpath[fi];

    if (node->contained)
      continue;

    if (last.ident == 0)
      //  Save the last dovetail node, but keep looking....
      last = *node;
    else
      //  ...for the next to last ID.
      prevID = node->ident;
  }

  return(last);
}






float Unitig::getAvgRho(void){

  if(ufpath.size() == 1)
    _avgRho = 1;
  if(_avgRho!=-1)
    return(_avgRho);


  // We will compute the average rho.
  //
  // Since rho is the length(unitig) - length(last fragment),
  //   and the length(last fragment) is ambiguous depending on which
  //   direction we are walking the unitig from.  We will take the average
  //   of the rhos through both directions.

  if (ufpath.empty())
    return 1;

  // Get first fragment's length
  int32 ident1         = ufpath[0].ident;
  int32 first_frag_len = FI->fragmentLength(ident1);
  assert(first_frag_len > 0);

  // Get last fragment's length
  int32 ident2        = ufpath[ufpath.size()-1].ident;
  int32 last_frag_len = FI->fragmentLength(ident2);
  assert(last_frag_len > 0);

  // Get average of first and last fragment lengths
  double avg_frag_len = (last_frag_len + first_frag_len)/2.0;

  // Compute average rho
  long unitig_length=getLength();
  _avgRho = unitig_length - avg_frag_len;

  if (_avgRho <= 0 ) {
    //fprintf(logFile, "Negative Rho ident1 "F_IID" ident2 "F_IID" unitig_length %d first_frag_len %d last_frag_len %d avg_frag_len %f\n",
    //        ident1, ident2, unitig_length, first_frag_len, last_frag_len, avg_frag_len);
    _avgRho = 1;
  }
  return(_avgRho);
}



void Unitig::setGlobalArrivalRate(float global_arrival_rate){
  _globalArrivalRate = global_arrival_rate;
}

void Unitig::setLocalArrivalRate(float local_arrival_rate){

  if ( local_arrival_rate < std::numeric_limits<float>::epsilon())
    _localArrivalRate = 0;
  else
    _localArrivalRate = local_arrival_rate;
}

float Unitig::getLocalArrivalRate(void){
  if (_localArrivalRate != -1 )
    return _localArrivalRate;
  setLocalArrivalRate((getNumFrags() - 1) / getAvgRho());
  return _localArrivalRate;
}


float Unitig::getCovStat(void){
  const float ln2=0.69314718055994530941723212145818;

  // Note that we are using numFrags in this calculation.
  //   If the fragments in the unitig are not randomly sampled
  //   from the genome, this calculation will be wrong.
  //   Ie. if fragments being assembled are from a locally
  //   sequenced batch, the region may look artificially repetitive.
  //   The value should really be "number of randomly sampled
  //   fragments in the unitig".

  //if(_globalArrivalRate == -1)
  //  fprintf(logFile, "You have not set the _globalArrivalRate variable.\n");

  float covStat = 0.0;

  if (_globalArrivalRate > 0.0)
    covStat = (getAvgRho() * _globalArrivalRate) - (ln2 * (getNumFrags() - 1));

  return(covStat);
}


void Unitig::reverseComplement(bool doSort) {

  //  If there are contained fragments, we need to sort by position to place them correctly after
  //  their containers.  If there are no contained fragments, sorting can break the initial unitig
  //  building.  When two frags start at position zero, we'll exchange the order.  Initial unitig
  //  building depends on having the first fragment added become the last fragment in the unitig
  //  after reversing.

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode  *frg = &ufpath[fi];

    frg->position.bgn = getLength() - frg->position.bgn;
    frg->position.end = getLength() - frg->position.end;

    //if (frg->contained != 0)
    //  doSort = true;

    assert(frg->position.bgn >= 0);
    assert(frg->position.end >= 0);
  }

  //  We've updated the positions of everything.  Now, sort or reverse the list, and rebuild the
  //  pathPosition map.

  if (doSort) {
    assert(0);
    sort();
  } else {
    std::reverse(ufpath.begin(), ufpath.end());

    for (int fi=0; fi<ufpath.size(); fi++)
      _pathPosition[ufpath[fi].ident] = fi;
  }
}



int
ufNodeCmp(const void *a, const void *b){
  ufNode *impa = (ufNode *)a;
  ufNode *impb = (ufNode *)b;

  int32 abgn = (impa->position.bgn < impa->position.end) ? impa->position.bgn : impa->position.end;
  int32 aend = (impa->position.bgn < impa->position.end) ? impa->position.end : impa->position.bgn;

  int32 bbgn = (impb->position.bgn < impb->position.end) ? impb->position.bgn : impb->position.end;
  int32 frag3p = (impb->position.bgn < impb->position.end) ? impb->position.end : impb->position.bgn;

  //  NEWSORT does not work.  When bubbles are popped, we add non-contained fragments to
  //  a unitig, but just stick them at the end of the list.  NEWSORT would then maintain
  //  this ordering, which is an error.
  //
#undef NEWSORT

#ifdef NEWSORT
  bool  aIsCont = OG->isContained(impa->ident);
  bool  bIsCont = OG->isContained(impb->ident);

  if ((aIsCont == false) && (bIsCont == false))
    //  Both dovetail nodes, keep same order
    return((int)impa->containment_depth - (int)impb->containment_depth);
#endif

  if (abgn != bbgn)
    //  Return negative for the one that starts first.
    return(abgn - bbgn);

  if (aend != frag3p)
    //  Return negative for the one that ends last.
    return(frag3p - aend);

#ifdef NEWSORT
  if (bIsCont == true)
    //  b is contained in a, so it comes after a.
    return(-1);

 if (aIsCont == true)
    //  a is contained in b, so it comes after b.
    return(1);
#endif

  //  Both contained, fallback on depth added, negative for earliest added
  return((int)impa->containment_depth - (int)impb->containment_depth);
}


void
Unitig::sort(void) {

#ifdef NEWSORT
  for (int fi=0; fi<ufpath.size(); fi++) {
    ufNode *f = &(ufpath[fi]);

    if (OG->isContained(f->ident) == false)
      f->containment_depth = fi;
  }
#endif

  qsort( &(ufpath.front()), getNumFrags(), sizeof(ufNode), &ufNodeCmp );

  for (int fi=0; fi<ufpath.size(); fi++)
    _pathPosition[ufpath[fi].ident] = fi;
}
