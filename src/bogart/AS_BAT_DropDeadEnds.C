
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
 *    Brian P. Walenz beginning on 2017-MAY-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_AssemblyGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_TigVector.H"

#include "AS_BAT_CreateUnitigs.H"




//  Find the next/previous read in the tig.  Skips contained reads if they are isolated.

ufNode *
findNextRead(Unitig *tig,
             ufNode *fn) {

  for (uint32 ni = tig->ufpathIdx(fn->ident)+1; ni < tig->ufpath.size(); ni++) {
    ufNode *nn = &tig->ufpath[ni];

    //  If nn is dovetail, return it.
    //      fn       -----------
    //      nn             ---------
    //
    if (fn->position.max() < nn->position.max())
      return(nn);

    //  Otherwise, if it intersects the next-next read, return it.
    //      fn       ----------------------
    //      nn             ---------
    //      next-next             -------
    //
    if ((ni + 1 < tig->ufpath.size()) &&
        (tig->ufpath[ni+1].position.min() < nn->position.max()))
      return(nn);
  }

  //  Otherwise, ran out of reads.

  return(NULL);
}



#if 0

//  Not used anymore.  Might be incorrect.

ufNode *
findPrevRead(Unitig *tig,
             ufNode *li) {

  //  A significant complication of working with reads on the 3' end is that they aren't sorted by
  //  their end position.  We get around this by saving a copy of the existing reads, reverse
  //  complementing that, and using the same method as in findNextRead().
  //
  //  Don't be clever and thing you can just reverse complement the tig; that can change order of
  //  reads, and we don't want to do that here.

  vector<ufNode>  ufcopy;

  ufcopy.resize(tig->ufpath.size());

  for (uint32 ii=0; ii<tig->ufpath.size(); ii++) {
    ufcopy[ii] = tig->ufpath[ii];

    ufcopy[ii].position.bgn = tig->getLength() - tig->ufpath[ii].position.bgn;
    ufcopy[ii].position.end = tig->getLength() - tig->ufpath[ii].position.end;
  }

  std::sort(ufcopy.begin(), ufcopy.end());

  //  ufpathIdx() won't work anymore, but li should be the first read.

  uint32 niPos=0;

  while (ufcopy[niPos].ident != li->ident)
    niPos++;

  //  Set 'fn' to that first node, and search for the next node.  This is nearly cut-n-paste from
  //  above (just replaced the return value with one that uses the ufpath).

  ufNode *fn = &ufcopy[niPos];

  for (uint32 ni = niPos+1; ni < tig->ufpath.size(); ni++) {
    ufNode *nn = &ufcopy[ni];

    if (fn->position.max() < nn->position.max())
      return(&tig->ufpath[ tig->ufpathIdx(nn->ident) ]);

    if ((ni + 1 < tig->ufpath.size()) &&
        (ufcopy[ni+1].position.min() < nn->position.max()))
      return(&tig->ufpath[ tig->ufpathIdx(nn->ident) ]);
  }

  //  Otherwise, ran out of reads.

  return(NULL);
}

#endif




uint32
dropDeadFirstRead(AssemblyGraph *AG,
                  Unitig        *tig) {

  ufNode *fn = tig->firstRead();
  ufNode *sn = findNextRead(tig, fn);

   //  No next read, keep fn in the tig.

  if (sn == NULL)
    return(0);

  //  Over all edges from the first read, look for any edge to something else.
  //
  //  If a contained edge to anything, read fn is good and should be kept.
  //
  //  Otherwise, decide which overlap we want to be using, based on the orientation of the read in
  //  the tig.  We assume that this is always the first read, which is OK, because the function name
  //  says so.  Any edge to anywhere means the read is good and should be kept.

  for (uint32 pp=0; pp<AG->getForward(fn->ident).size(); pp++) {
    BestPlacement  &pf = AG->getForward(fn->ident)[pp];

    if (pf.bestC.b_iid > 0)
      return(0);

    if (((fn->position.isForward() == true)  && (pf.best5.b_iid != 0)) ||
        ((fn->position.isForward() == false) && (pf.best3.b_iid != 0)))
      return(0);
  }

  //  But no edge means we need to check the second read.  If it has an edge, then we infer the
  //  first read is bogus and should be removed.  If it also has no edge (except to the first read,
  //  duh) then we know nothing: this could be novel sequence or it could be the same garbage that
  //  is infecting the first read.
  //
  //  This is basically the same as the previous loop, except we also need to exclude edges to the
  //  first read.  Well, and that if the second read has an edge we declare the first read to be
  //  junk.  That's also a bit of a difference from the previous loop.

  for (uint32 pp=0; pp<AG->getForward(sn->ident).size(); pp++) {
    BestPlacement  &pf = AG->getForward(sn->ident)[pp];

    if ((pf.bestC.b_iid > 0) && (pf.bestC.b_iid != fn->ident))
      return(fn->ident);

    if (((sn->position.isForward() == true)  && (pf.best5.b_iid != 0) && (pf.best5.b_iid != fn->ident)) ||
        ((sn->position.isForward() == false) && (pf.best3.b_iid != 0) && (pf.best3.b_iid != fn->ident)))
      return(fn->ident);
  }

  //  Otherwise, the second read had only edges to the first read, and we should keep the first
  //  read.

  return(0);
}



void
dropDeadEnds(AssemblyGraph  *AG,
             TigVector      &tigs) {

  uint32  numF = 0;  //  Number of first-read drops
  uint32  numL = 0;  //  Number of last-read drops
  uint32  numB = 0;  //  Number of both-first-and-last-read drops
  uint32  numT = 0;  //  Number of tigs mucked with

  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig  *tig = tigs[ti];

    if ((tig == NULL) || (tig->ufpath.size() < 2))  //  No tig, or don't care.
      continue;

    uint32  fn = dropDeadFirstRead(AG, tig);        //  Decide if the first read is junk.

    tig->reverseComplement();                       //  Flip.
    uint32  ln = dropDeadFirstRead(AG, tig);        //  Decide if the last (now first) read is junk.
    tig->reverseComplement();                       //  Flip back.

    if ((fn == 0) && (ln == 0))                     //  Nothing to remove, just get out of here.
      continue;

    //  At least one read needs to be kicked out.  Make new tigs for everything.
 
    char   fnMsg[80] = {0};   Unitig  *fnTig = NULL;
    char   nnMsg[80] = {0};   Unitig  *nnTig = NULL;  int32  nnOff = INT32_MAX;
    char   lnMsg[80] = {0};   Unitig  *lnTig = NULL;

    if (fn > 0)
      fnTig = tigs.newUnitig(false);

    if (tig->ufpath.size() > (fn > 0) + (ln > 0))
      nnTig = tigs.newUnitig(false);

    if (ln > 0)
      lnTig = tigs.newUnitig(false);

    //  Count what we do

    numT++;

    if (fnTig)            numF++;
    if (fnTig && lnTig)   numB++;
    if (lnTig)            numL++;

    //  Move reads to their new unitig.

    sprintf(fnMsg, "                                      ");
    sprintf(nnMsg, "                          ");
    sprintf(lnMsg, "");

    for (uint32 cc=0, tt=0; tt<tig->ufpath.size(); tt++) {
      ufNode  &read = tig->ufpath[tt];

      if        (read.ident == fn) { 
        sprintf(fnMsg, "first read %9u to tig %7u --", read.ident, fnTig->id());
        fnTig->addRead(read, -read.position.min(), false);

      } else if (read.ident == ln) {
        sprintf(lnMsg, "-- last read %9u to tig %7u", read.ident, lnTig->id());
        lnTig->addRead(read, -read.position.min(), false);

      } else {
        if (nnOff == INT32_MAX) {
          sprintf(nnMsg, "other reads to tig %7u", nnTig->id());
          nnOff = read.position.min();
        }
        nnTig->addRead(read, -nnOff, false);
      }
    }

    writeLog("dropDeadEnds()-- tig %7u --> %s %s %s\n", tig->id(), fnMsg, nnMsg, lnMsg);

    if (fnTig)  fnTig->cleanUp();   //  Probably not neeeded, but cheap.
    if (lnTig)  lnTig->cleanUp();   //  Probably not neeeded, but cheap.
    if (nnTig)  nnTig->cleanUp();   //  Most likely needed.

    //  Old tig is now junk.

    delete  tigs[ti];
    tigs[ti] = NULL;
  }

  writeStatus("dropDeadEnds()-- Modified %u tigs.  Dropped %u first and %u last reads, %u tig%s had both reads dropped.\n",
              numT, numF, numL, numB, (numB == 1) ? "" : "s");
}
