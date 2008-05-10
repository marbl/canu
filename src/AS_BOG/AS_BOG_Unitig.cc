
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

#include "AS_BOG_Unitig.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include <algorithm>

#undef max


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


void Unitig::computeFragmentPositions(FragmentInfo *fi, BestOverlapGraph* bog_ptr) {
  if (dovetail_path_ptr == NULL || dovetail_path_ptr->empty())
    return;
  // we need to determine which orientation the first frag is in
  // so we peek ahead to the 2nd frag to find out
  DoveTailIter iter = dovetail_path_ptr->begin();
  int frag_begin,fragPrevEnd,fragNextEnd;
  frag_begin = fragPrevEnd = fragNextEnd = 0;
  iuid first = iter->ident;
  int frag_end = fi->fragmentLength(first);
  assert(frag_end > 0);
  fragPrevEnd = frag_end;
  if(  dovetail_path_ptr->size() == 1) {
    // still set singleton's coords
    dovetail_path_ptr->front().position.bgn = frag_begin;
    dovetail_path_ptr->front().position.end = frag_end;
    return;
  }
  BestEdgeOverlap* bestEdge = bog_ptr->getBestEdgeOverlap(first, FIVE_PRIME );
  BestEdgeOverlap* threeP = bog_ptr->getBestEdgeOverlap(first, THREE_PRIME );
  fprintf(stderr,"1stFrag %d beg %d end %d\n",
          first, frag_begin, frag_end );
  iter++;
  if (iter->ident == bestEdge->frag_b_id) {
    dovetail_path_ptr->front().position.bgn = frag_end;
    dovetail_path_ptr->front().position.end = frag_begin;
    fprintf(stderr,"Go off 5' edge and get frag %d\n",iter->ident);
  } else {
    // doesn't handle contianed, only backbone
    fprintf(stderr,"Go off 3' edge and get frag %d %d\n",iter->ident,threeP->frag_b_id);
    assert( iter->ident == threeP->frag_b_id );
    bestEdge = threeP;
    dovetail_path_ptr->front().position.bgn = frag_begin;
    dovetail_path_ptr->front().position.end = frag_end;
  }
  fprintf(stderr,"ahang %d bhang %d\n", bestEdge->ahang, bestEdge->bhang);
  if ( bestEdge->ahang > 0 )
    frag_begin = bestEdge->ahang;
  else 
    frag_begin = -bestEdge->bhang;

  fragment_end_type whichEnd = FIVE_PRIME; // the end to walk off for next frag
  if (bestEdge->bend == FIVE_PRIME)
    whichEnd = THREE_PRIME;

  for (; iter != dovetail_path_ptr->end(); iter++ ) {
    DoveTailNode* node = &*iter;
    //DoveTailNode* nextNode = &*(iter+1);
    bestEdge = bog_ptr->getBestEdgeOverlap( node->ident, whichEnd);
    int currLen = fi->fragmentLength(node->ident);
    assert(currLen > 0);
    // The end of the fragment can be calulated as start + length
    // or as end of previous frag + b_hang. They should be the same
    // if there are no gaps, but with gaps the 2nd method should be better
    frag_end = frag_begin + currLen;

    fprintf(stderr,"Frag %d len %d beg %d end %d ahang %d bhang %d end %s bfrg %d bend %s\n",
            node->ident, currLen, frag_begin, frag_end,
            bestEdge->ahang, bestEdge->bhang, whichEnd == FIVE_PRIME ? "5'":"3'",
            bestEdge->frag_b_id, bestEdge->bend == FIVE_PRIME ? "5'":"3'");

    int end;
    // pick the smallest end that's greater then the previous end
    // this is critical for preserving the correct frag order after
    // a reverse complement in the unitig merging code
    if (frag_end < fragPrevEnd) {
      end = fragNextEnd;
    } else if (fragNextEnd < fragPrevEnd) {
      end = frag_end;
    } else {
      end = fragNextEnd > frag_end ? frag_end : fragNextEnd;
    }

    if(whichEnd == FIVE_PRIME){
      node->position.bgn = end;
      node->position.end = frag_begin;
    }else {
      node->position.bgn = frag_begin;
      node->position.end = end;
    }
    fragNextEnd = frag_end;
    fragPrevEnd = end;

    // Prep the start position of the next fragment
    if (bestEdge->ahang < 0 && bestEdge->bhang < 0 ) {
      fragNextEnd -= bestEdge->ahang ;
      frag_begin  -= bestEdge->bhang ;
    } else {
      fragNextEnd += bestEdge->bhang ;
      frag_begin  += bestEdge->ahang ;
    }
    if ( bestEdge->bend == FIVE_PRIME ) {
      whichEnd = THREE_PRIME;
    } else {
      whichEnd = FIVE_PRIME;
    }
  }
}


void Unitig::addFrag( DoveTailNode node, int offset, bool report) {

  node.position.bgn += offset;
  node.position.end += offset;

  assert(node.position.bgn >= 0);
  assert(node.position.end >= 0);

  // keep track of the unitig a frag is in
  _inUnitig[ node.ident ] = id();

  // keep track of max position in unitig
  int frgEnd = MAX( node.position.bgn, node.position.end);
  if ( frgEnd > _length)
    _length = frgEnd;

  dovetail_path_ptr->push_back(node);

  if (report)
    if (node.contained)
      fprintf(stderr, "Added frag %d to unitig %d at %d,%d (contained in %d)\n", node.ident, id(), node.position.bgn, node.position.end, node.contained);
    else
      fprintf(stderr, "Added frag %d to unitig %d at %d,%d\n", node.ident, id(), node.position.bgn, node.position.end);
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


long Unitig::getLength(void){
  return _length;
}

long Unitig::getNumFrags(void){
  return(dovetail_path_ptr->size());
}

// This is a placeholder, random frags should not contain guides, or
// other fragments that are not randomly sampled across the whole
// genome.
//
long Unitig::getNumRandomFrags(void){
  return(getNumFrags());
}

void Unitig::shiftCoordinates(int offset) {
  //simple version
  if (dovetail_path_ptr->empty()) {
  } else {
    //        fprintf(stderr,"shift unitig by %d 1st frg %d\n",
    //                offset,dovetail_path_ptr->front().ident);
  }
  DoveTailIter iter = dovetail_path_ptr->begin();
  for(; iter != dovetail_path_ptr->end(); iter++) {
    iter->position.bgn += offset;
    assert( iter->position.bgn >= 0);
    iter->position.end += offset;
    assert( iter->position.end >= 0);
  }
}

void Unitig::reverseComplement() {
  int length = this->getLength();
  DoveTailNode first = dovetail_path_ptr->front();
  DoveTailNode last = dovetail_path_ptr->back();
  DoveTailIter iter = dovetail_path_ptr->begin();
  for(; iter != dovetail_path_ptr->end(); iter++) {
    iter->position.bgn = length - iter->position.bgn;
    assert( iter->position.bgn >= 0);
    iter->position.end = length - iter->position.end;
    assert( iter->position.end >= 0);
  }
  reverse(dovetail_path_ptr->begin(),dovetail_path_ptr->end());
  first = dovetail_path_ptr->front();
  last = dovetail_path_ptr->back();
}

void Unitig::reverseComplement(int offset, BestOverlapGraph *bog_ptr) {
  iuid notUsed;
  DoveTailNode last  = getLastBackboneNode(notUsed);
  int lastEnd = last.position.end > last.position.bgn ? 
    last.position.end : last.position.bgn;
  int prevBeg = 0;
  DoveTailPath contains;
  DoveTailPath *revP = new DoveTailPath;
  DoveTailPath::reverse_iterator rend = dovetail_path_ptr->rend();
  DoveTailPath::reverse_iterator addIter = dovetail_path_ptr->rbegin();
  for(; addIter != rend; addIter++) {
    addIter->position.bgn = lastEnd - addIter->position.bgn + offset;
    assert( addIter->position.bgn >= 0);
    addIter->position.end = lastEnd - addIter->position.end + offset;
    assert( addIter->position.end >= 0);
    if (addIter->contained != 0) {
#ifdef NEW_UNITIGGER_INTERFACE
      int tmp = addIter->ahang;
      addIter->ahang = -addIter->bhang;
      addIter->bhang = -tmp;
#endif
      contains.push_back( *addIter );
    } else {
#ifdef NEW_UNITIGGER_INTERFACE
      int left = addIter->position.end < addIter->position.bgn ? 
        addIter->position.end : addIter->position.bgn;
      assert( prevBeg <= left ); 
      DoveTailPath::reverse_iterator prev = addIter+1;
      while (prev != rend &&
             prev->contained != 0) {
        prev = prev+1;
      }
      if (prev == rend) {

        // After reverse, need to switch to other edge
        BestEdgeOverlap* other;
        iuid id = addIter->ident;
        if (addIter->position.bgn > addIter->position.end) {
          other = bog_ptr->getBestEdgeOverlap( id, FIVE_PRIME );
        } else {
          other = bog_ptr->getBestEdgeOverlap( id, THREE_PRIME );
        }
        addIter->ident2 = other->frag_b_id;
        if (other->ahang < 0 && other->bhang < 0 ) {
          addIter->ahang = -other->bhang;
          addIter->bhang = -other->ahang;
        } else {
          addIter->ahang = other->ahang;
          addIter->bhang = other->bhang;
        }
      } else {
        assert( prev->contained == 0 );
        addIter->ident2 = prev->ident;
        addIter->bhang = prev->ahang;
        addIter->ahang = prev->bhang;
      }
#endif
      revP->push_back( *addIter );
      if (contains.size() > 0) {
        revP->insert(revP->end(), contains.begin(), contains.end());
        contains.clear();
      }
    }
  }
  delete dovetail_path_ptr;
  dovetail_path_ptr = revP;
}

//
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
                           BestContainmentMap *bestCtn,
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
    BestContainment &best   = (*bestCtn)[ fragId ];

    if (best.isPlaced)
      continue;

    assert( best.container == containerId );

    (*containPartialOrder)[ fragId ] = level;

    //if (level > maxContainDepth)
    //    maxContainDepth = level;

    DoveTailNode frag;

    frag.type         = AS_READ;
    frag.ident        = fragId;
    frag.contained    = containerId;
    frag.delta_length = 0;
    frag.delta        = NULL;

    if(containerPos.bgn < containerPos.end) {
      //  Container is forward
      frag.position.bgn = containerPos.bgn + best.a_hang;  //  BPW says "looks ok"
      frag.position.end = containerPos.end + best.b_hang;
#ifdef NEW_UNITIGGER_INTERFACE
      frag.ident2       = containerId;
      frag.ahang        = best.a_hang;
      frag.bhang        = best.b_hang;
#endif

    } else if (containerPos.bgn > containerPos.end) {
      //  Container is reverse
      frag.position.bgn = containerPos.bgn - best.a_hang;  //  BPW says "suspicious"
      frag.position.end = containerPos.end - best.b_hang;
#ifdef NEW_UNITIGGER_INTERFACE
      frag.ident2       = containerId;
      frag.ahang        = - best.b_hang;   //  consensus seems to want these reversed
      frag.bhang        = - best.a_hang;
#endif
    }else{
      fprintf(stderr, "Container size is zero?\n");
      assert(containerPos.bgn != containerPos.end);
    }

    // Swap ends if containee is not same strand as container

    if(!best.sameOrientation){
      int tmp          = frag.position.bgn;
      frag.position.bgn = frag.position.end;
      frag.position.end = tmp;
    }

    addFrag( frag, 0, false );
    best.isPlaced = true;
    placeContains(cMap, bestCtn, frag.ident, frag.position, level+1);
  }
}

void Unitig::recomputeFragmentPositions(ContainerMap &cMap,
                                        BestContainmentMap *bestContain,
                                        BestOverlapGraph *bog_ptr) {
  iuid lastFrag = 0;
#ifdef NEW_UNITIGGER_INTERFACE
  iuid nextFrag = 0;
#endif

  // place dovetails in a row

  if (dovetail_path_ptr == NULL)
    return;

  containPartialOrder = new std::map<iuid,int>;

  for (int i=0; i < dovetail_path_ptr->size(); i++) {
    DoveTailNode *dt = &(*dovetail_path_ptr)[i];

#ifdef NEW_UNITIGGER_INTERFACE
    if ( nextFrag != 0 )
      assert( nextFrag == dt->ident);
    nextFrag = dt->ident2;
#endif
    lastFrag = dt->ident;

    placeContains(cMap, bestContain, dt->ident, dt->position, 1);
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
