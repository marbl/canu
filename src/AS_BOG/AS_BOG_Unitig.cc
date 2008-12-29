
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

static const char *rcsid = "$Id: AS_BOG_Unitig.cc,v 1.9 2008-12-29 16:07:17 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_Unitig.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#undef max


extern FragmentInfo     *debugfi;


// various class static methods and variables
static std::map<iuid,int>* containPartialOrder;

iuid Unitig::nextId        = 1;
iuid* Unitig::_inUnitig = NULL;


Unitig::Unitig(bool report){
  _localArrivalRate = -1;
  _covStat          = FLT_MAX;
  _length           = -1;
  _avgRho           = -1;
  dovetail_path_ptr = new DoveTailPath;
  _id               = nextId++;

  if (report)
    fprintf(stderr, "Creating Unitig %d\n", _id);
}


Unitig::~Unitig(void){
  delete dovetail_path_ptr;
}



DoveTailNode Unitig::getLastBackboneNode(iuid &prevId) {
  DoveTailNode lastNonContain;

  memset(&lastNonContain, 0, sizeof(DoveTailNode));

  for(DoveTailPath::reverse_iterator rIter = dovetail_path_ptr->rbegin(); rIter != dovetail_path_ptr->rend(); rIter++) {
    if (rIter->contained == 0) {
      if (lastNonContain.ident == 0) {
        lastNonContain = *rIter;
      } else {
        prevId = rIter->ident;
        return lastNonContain;
      }
    }
  }
  return lastNonContain;
}


void Unitig::addFrag( DoveTailNode node, int offset, bool report) {

  node.position.bgn += offset;
  node.position.end += offset;

  // keep track of the unitig a frag is in
  _inUnitig[ node.ident ] = id();

  // keep track of max position in unitig
  int frgEnd = MAX( node.position.bgn, node.position.end);
  if ( frgEnd > _length)
    _length = frgEnd;

  dovetail_path_ptr->push_back(node);

  report = true;

  int32 len = debugfi->fragmentLength(node.ident);
  int32 pos = (node.position.end > node.position.bgn) ? (node.position.end - node.position.bgn) : (node.position.bgn - node.position.end);

  if ((report) || (node.position.bgn < 0) || (node.position.end < 0))
    if (node.contained)
      fprintf(stderr, "Added frag %d (len %d) to unitig %d at %d,%d (lendiff %d) (contained in %d)\n",
              node.ident, len, id(), node.position.bgn, node.position.end,
              pos - len,
              node.contained);
    else
      fprintf(stderr, "Added frag %d (len %d) to unitig %d at %d,%d (lendiff %d)\n",
              node.ident, len, id(), node.position.bgn, node.position.end,
              pos - len);

  assert(node.position.bgn >= 0);
  assert(node.position.end >= 0);
}

void Unitig::addContainedFrag(DoveTailNode node, BestContainment *bestcont, bool report) {

  //  This will add a contained fragment to a unitig, adjusting the
  //  position as needed.  It is only needed when moving a contained
  //  read from unitig A to unitig B.  It is NOT needed when
  //  rebuilding a unitig.

  assert(node.contained == bestcont->container);

  //  Apparently, no way to retrieve a single fragment from a
  //  unitig without searching for it.
  //
  //  Orientation should be OK.  All we've done since the
  //  unitig was built was to split at various spots.  But we
  //  need to adjust the location of the read.
  //
  for (DoveTailIter it=dovetail_path_ptr->begin(); it != dovetail_path_ptr->end(); it++)
    if (it->ident == node.contained) {
      int offset = MIN(it->position.bgn, it->position.end) + bestcont->a_hang;
      int adj    = MIN(node.position.bgn, node.position.end);

      node.position.bgn += offset - adj;
      node.position.end += offset - adj;
    }

  addFrag(node, 0, report);

  //  Bump that new fragment up to be in the correct spot -- we can't
  //  use the sort() method on Unitig, since we lost the
  //  containPartialOrder.
  //
  int             i = dovetail_path_ptr->size() - 1;
  DoveTailNode   *f = &dovetail_path_ptr->front();

  //  Only needed if the frag we just added (i) begins before the second to last frag.

  if (MIN(f[i].position.bgn, f[i].position.end) < MIN(f[i-1].position.bgn, f[i-1].position.end)) {
    DoveTailNode    containee    = f[i];
    int             containeeMin = MIN(containee.position.bgn, containee.position.end);

    while ((i > 0) &&
           (containee.contained != f[i-1].ident) &&
           (containeeMin < MIN(f[i-1].position.bgn, f[i-1].position.end))) {
      f[i] = f[i-1];
      i--;
    }

    f[i] = containee;
  }
}




float Unitig::getAvgRho(FragmentInfo *fi){

  if(dovetail_path_ptr->size() == 1)
    _avgRho = 1;
  if(_avgRho!=-1)
    return(_avgRho);


  // We will compute the average rho.
  //
  // Since rho is the length(unitig) - length(last fragment),
  //   and the length(last fragment) is ambiguous depending on which
  //   direction we are walking the unitig from.  We will take the average
  //   of the rhos through both directions.

  DoveTailPath::const_iterator dtp_iter;
  if (dovetail_path_ptr == NULL || dovetail_path_ptr->empty()) {
    //fprintf(stderr,"NULL dovetailpath\n");
    return 1;
  }

  // Get first fragment's length
  dtp_iter=dovetail_path_ptr->begin();
  int ident1 = dtp_iter->ident;
  long first_frag_len = fi->fragmentLength(dtp_iter->ident);
  assert(first_frag_len > 0);

  // Get last fragment's length
  dtp_iter=dovetail_path_ptr->end();
  dtp_iter--;
  int ident2 = dtp_iter->ident;
  long last_frag_len = fi->fragmentLength(dtp_iter->ident);
  assert(last_frag_len > 0);

  // Get average of first and last fragment lengths
  double avg_frag_len = (last_frag_len + first_frag_len)/2.0;

  // Compute average rho
  long unitig_length=getLength();
  _avgRho = unitig_length - avg_frag_len;

  if (_avgRho <= 0 ) {
    fprintf(stderr, "Negative Rho ident1 "F_IID" ident2 "F_IID" unitig_length %d first_frag_len %d last_frag_len %d avg_frag_len %f\n",
            ident1, ident2, unitig_length, first_frag_len, last_frag_len, avg_frag_len);
    _avgRho = 1;
  }
  return(_avgRho);
}


float Unitig::_globalArrivalRate = -1;

void Unitig::setGlobalArrivalRate(float global_arrival_rate){
  _globalArrivalRate=global_arrival_rate;
}
void Unitig::setLocalArrivalRate(float local_arrival_rate){

  if ( local_arrival_rate < std::numeric_limits<float>::epsilon())
    _localArrivalRate = 0;
  else
    _localArrivalRate = local_arrival_rate;
}
float Unitig::getLocalArrivalRate(FragmentInfo *fi){
  if (_localArrivalRate != -1 )
    return _localArrivalRate;
  setLocalArrivalRate((getNumFrags() - 1) / getAvgRho(fi));
  return _localArrivalRate;
}


float Unitig::getCovStat(FragmentInfo *fi){

  const float ln2=0.69314718055994530941723212145818;

  // Note that we are using numFrags in this calculation.
  //   If the fragments in the unitig are not randomly sampled
  //   from the genome, this calculation will be wrong.
  //   Ie. if fragments being assembled are from a locally
  //   sequenced batch, the region may look artificially repetitive.
  //   The value should really be "number of randomly sampled
  //   fragments in the unitig".

  if(_globalArrivalRate==-1)
    fprintf(stderr, "You have not set the _globalArrivalRate variable.\n");

  if(_covStat == FLT_MAX){
    if(_globalArrivalRate > 0.0){
      _covStat = (getAvgRho(fi) * _globalArrivalRate) - (ln2 * (getNumFrags() -1));
    }else{
      _covStat = 0.0;
    }
  }

  return(_covStat);
}


void Unitig::reverseComplement() {
  DoveTailIter iter  = dovetail_path_ptr->begin();

  for(; iter != dovetail_path_ptr->end(); iter++) {
    iter->position.bgn = getLength() - iter->position.bgn;
    iter->position.end = getLength() - iter->position.end;

    assert(iter->position.bgn >= 0);
    assert(iter->position.end >= 0);
  }

  std::reverse(dovetail_path_ptr->begin(), dovetail_path_ptr->end());
}

// Recursively place all contains under this one into the FragmentPositionMap
//
// Compute assuming that containee is the same orientation as container
//  if(cntnr_intvl.begin < cntnr_intvl.end)
//
// Container is in forward direction
//
// |=========F0============|
// 0          |=============CER=============>|
//                    |=====CEE=====>|
//
// |---Cro----|
//            |--Ceo--|
// |-------Cep--------|
//
// Cro = container offset from beginning of unitig = cntnr_intvl.begin
// Ceo = containee offset from 5' end of container = cntee->olap_offset
// Cep = containee offset from beginning of unitig = cntee_intvl.begin
// CEE fragment can be from either orientation since
//   definition of olap_offset is based on 3' origin.
//
// else if(cntnr_intvl.begin > cntnr_intvl.end)
//
// Container is in reverse direction
//
// |=========F0============|
// 0          |<============CER==============|
//                    |<====CEE======|
//
// |---Cro----|
//                                   |--Ceo--|
// |-------Cep-----------------------|
//
// Cro = container offset from beginning of unitig = cntnr_intvl.end
// Ceo = containee offset from 5' end of container = cntee->olap_offset
// Cep = containee offset from beginning of unitig = cntee_intvl.end
// CEE fragment can be from either orientation since
//   definition of olap_offset is based on 3' origin.

void Unitig::placeContains(const ContainerMap &cMap,
                           BestOverlapGraph *bog,
                           const iuid containerId,
                           const SeqInterval containerPos,
                           const int level) {
  if (cMap.size() == 0)
    return;

  ContainerMap::const_iterator ctmp_itr = cMap.find( containerId );

  if (ctmp_itr == cMap.end() )
    return;

  for(ContaineeList::const_iterator ci  = ctmp_itr->second.begin();
      ci != ctmp_itr->second.end();
      ci++) {

    iuid             fragId = *ci;
    BestContainment *best   =  bog->getBestContainer(fragId);

    if (best->isPlaced)
      continue;

    assert(best->container == containerId);
    assert(containerPos.bgn != containerPos.end);

    (*containPartialOrder)[ fragId ] = level;

    //if (level > maxContainDepth)
    //    maxContainDepth = level;

    DoveTailNode frag;

    //  The second condition of the position field looks strange.

    frag.type         = AS_READ;
    frag.ident        = fragId;
    frag.contained    = containerId;
    frag.parent       = containerId;
    frag.sourceInt    = -1;

    //  The 'best' overlap picture is always the same:
    //
    //  -------------------->    containerPos
    //  a>0   -------->   b<0    frag
    //
    if (containerPos.bgn < containerPos.end) {
      //  Container is oriented same as overlap.
      frag.ahang        = best->a_hang;
      frag.bhang        = best->b_hang;

      if (best->sameOrientation) {
        //  Containee too.  Straightforward.
        frag.position.bgn = containerPos.bgn + frag.ahang;
        frag.position.end = containerPos.end + frag.bhang;
      } else {
        //  But containee is flipped and is now 'reversed'.
        frag.position.end = containerPos.bgn + frag.ahang;
        frag.position.bgn = containerPos.end + frag.bhang;
      }
    } else {
      //  Container is backwards from overlap.
      frag.ahang        = -best->b_hang;
      frag.bhang        = -best->a_hang;

      if (best->sameOrientation) {
        //  Containee too.  The tricky case.
        //
        //  We've already adjusted frag.ahang and frag.bhang to
        //  account for the flipped overlap.  The subtraction below
        //  accounts for the flipped coordinates.  The thing called
        //  'bgn' is really the end of the overlap, and to make the
        //  frag position smaller, we need to subtract the (adjusted
        //  and now positive) ahang.
        //
        //  end               bgn
        //  <--------------------    containerPos
        //  bhang <-------- ahang    frag
        //
        //  Finally, because the frag is also reversed the 'end'
        //  variable is storing the start position, and we add the new
        //  ahang to the start.
        //  
        frag.position.end = containerPos.end + frag.ahang;
        frag.position.bgn = containerPos.bgn + frag.bhang;
      } else {
        //  But containee is flipped, and is now 'forward'.
        frag.position.bgn = containerPos.end + frag.ahang;
        frag.position.end = containerPos.bgn + frag.bhang;
      }
    }

    frag.delta_length = 0;
    frag.delta        = NULL;

    //  See AS_BOG_UnitigGraph.c::populateUnitig() around line 557 for
    //  why we reset the end coordinate.
    //
    //  Containments are particularily painful.  A beautiful example:
    //  a fragment of length 253bp is contained in a fragment of
    //  length 251bp (both hangs are zero).  In this case, the
    //  "ahang+length" method fails, placing the contained fragment
    //  outside the container (and if backwards oriented, _BEFORE_ the
    //  contained fragment).  The "ahang,bhang" method works here, but
    //  fails on other instances, shrinking deep containments to
    //  nothing.
    //
    //  We can use either method first, then adjust using the other method.
    //
    //  We'll use 'ahang,bhang' first (mostly because it was already done, and we need
    //  to compute those values anyway) then reset the end based on the length, limited
    //  to maintain a containment relationship.
    //
#warning not knowing the overlap length really hurts.
    if (frag.position.bgn < frag.position.end) {
      frag.position.end = frag.position.bgn + bog->fragmentLength(fragId);
      if (frag.position.end > MAX(containerPos.bgn, containerPos.end))
        frag.position.end = MAX(containerPos.bgn, containerPos.end);
    } else {
      frag.position.bgn = frag.position.end + bog->fragmentLength(fragId);
      if (frag.position.bgn > MAX(containerPos.bgn, containerPos.end))
        frag.position.bgn = MAX(containerPos.bgn, containerPos.end);
    }

    {
      uint32  cbgn = MIN(containerPos.bgn, containerPos.end);
      uint32  cend = MAX(containerPos.bgn, containerPos.end);

      uint32  fbgn = MIN(frag.position.bgn, frag.position.end);
      uint32  fend = MAX(frag.position.bgn, frag.position.end);

      if ((fbgn < cbgn) || (cend < fend)) {
        fprintf(stderr, "ERROR:  Fragment became uncontained in the layout.\n");
        fprintf(stderr, "containerPos "F_IID"\t%d %d\n", containerId, containerPos.bgn, containerPos.end);
        fprintf(stderr, "fragment     "F_IID"\t%d %d\n", fragId,      frag.position.bgn, frag.position.end);
      }

      assert(cbgn <= fbgn);
      assert(fend <= cend);
    }

    addFrag(frag, 0, false);
    best->isPlaced = true;
    placeContains(cMap, bog, frag.ident, frag.position, level+1);
  }
}

void Unitig::recomputeFragmentPositions(ContainerMap &cMap,
                                        BestOverlapGraph *bog_ptr) {

  if (dovetail_path_ptr == NULL)
    return;

  containPartialOrder = new std::map<iuid,int>;

  for (int i=0; i < dovetail_path_ptr->size(); i++) {
    DoveTailNode *dt = &(*dovetail_path_ptr)[i];

    fprintf(stderr, "PlaceContains in frag %d position %d,%d\n", dt->ident, dt->position.bgn, dt->position.end);

    placeContains(cMap, bog_ptr, dt->ident, dt->position, 1);
  }

  this->sort();

  delete containPartialOrder;
  containPartialOrder = NULL;
}


int IntMultiPosCmp(const void *a, const void *b){
  IntMultiPos *impa=(IntMultiPos*)a;
  IntMultiPos *impb=(IntMultiPos*)b;
  long aleft,aright,bleft,bright;
  if (impa->position.bgn < impa->position.end) {
    aleft  = impa->position.bgn;
    aright = impa->position.end;
  } else {
    aright = impa->position.bgn;
    aleft  = impa->position.end;
  }
  if (impb->position.bgn < impb->position.end) {
    bleft  = impb->position.bgn;
    bright = impb->position.end;
  } else {
    bright = impb->position.bgn;
    bleft  = impb->position.end;
  }
  if(aleft!=bleft) {
    return(aleft - bleft);
  }
  else if (aright != bright) {
    if (impa->contained==0 && impb->contained!=0)
      return(-1);
    else if (impa->contained!=0 && impb->contained==0)
      return(1);
    if (impa->contained!=0 || impb->contained!=0)
      return(bright - aright);
    else
      return(aright - bright);
  }
  else {
    if(impa->contained == impb->ident)
      return(1);
    if(impb->contained == impa->ident)
      return(-1);
    if(impa->contained!=0 && impb->contained!=0)
      if (containPartialOrder == NULL)
        return(0);
      else
        return((*containPartialOrder)[impa->ident] - (*containPartialOrder)[impb->ident]);
    if(impa->contained!=0)
      return(1);
    if(impb->contained!=0)
      return(-1);
    return(0);
  }
}

void Unitig::sort() {
  qsort( &(dovetail_path_ptr->front()), getNumFrags(), sizeof(IntMultiPos), &IntMultiPosCmp );
}
