
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-APR-09 to 2015-MAY-09
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-01
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "stashContains.H"

//  Replace the children list in tig with one that has fewer contains.  The original
//  list is returned.
savedChildren *
stashContains(tgTig       *tig,
              double       maxCov,
              bool         beVerbose) {

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

  double totlCov  = (double)nBase / hiEnd;

  saved->numContains = nCont;
  saved->covContain  = (double)nBaseCont / hiEnd;
  saved->percContain = 100.0 * nBaseCont / nBase;;

  saved->numDovetails = nBack;
  saved->covDovetail  = (double)nBaseDove / hiEnd;
  saved->percDovetail = 100.0 * nBaseDove / nBase;;

  if (beVerbose)
    saved->reportDetected(stderr, tig->tigID());

  //  If the tig has more coverage than allowed, throw out some of the contained reads.

  if ((totlCov  >= maxCov) &&
      (maxCov   > 0)) {
    std::sort(posLen, posLen + nOrig, greater<readLength>());  //  Sort by length, larger first

    nBaseSave = 0.0;

    for (uint32 ii=0; ii < nOrig; ii++) {

      if (ii > 0)
        assert(posLen[ii-1].len >= posLen[ii].len);

      if (isBack[posLen[ii].idx])
        //  Already a backbone read.
        continue;

      if ((double)(nBaseSave + nBaseDove) / hiEnd < maxCov) {
        isBack[posLen[ii].idx] = true;  //  Save it.

        nSave++;
        nBaseSave += posLen[ii].len;
      }
    }


    saved->numContainsRemoved = nOrig - nBack - nSave;
    saved->covContainsRemoved = (double)(nBaseCont - nBaseSave) / hiEnd;

    saved->numContainsSaved   = nSave;
    saved->covContainsSaved   = (double)nBaseSave / hiEnd;

    if (beVerbose)
      saved->reportRemoved(stderr, tig->tigID());

    //  For all the reads we saved, copy them to a new children list in the tig

    tig->_childrenLen = 0;
    tig->_childrenMax = nBack + nSave;
    tig->_children    = new tgPosition [tig->_childrenMax];  //  The original is in savedChildren now

    for (uint32 fi=0; fi<nOrig; fi++) {
      if (isBack[fi] == false)
        continue;

      //fprintf(stderr, "    ident %9d position %6d %6d\n",
      //        saved->children[fi].ident(), saved->children[fi].bgn(), children[fi].end());

      tig->_children[tig->_childrenLen++] = saved->children[fi];
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
      int32  nmin = sf * saved->children[fi].min();
      int32  nmax = sf * saved->children[fi].max();

      if (nmin > newMax)  nmin = newMax;
      if (nmax > newMax)  nmax = newMax;

      saved->children[fi].setMinMax(nmin, nmax);
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

