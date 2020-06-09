
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "AS_BAT_SplitTig.H"


uint32
splitTigByRegions(TigVector                &tigs,
                  Unitig                   *tig,
                  vector<breakPointCoords> &BP,
                  Unitig                  **newTigs,
                  int32                    *lowCoord,
                  uint32                   *nRepeat,
                  uint32                   *nUnique,
                  bool                      doMove) {

  if (doMove == true) {
    memset(newTigs,  0, sizeof(Unitig *) * BP.size());
    memset(lowCoord, 0, sizeof(int32)    * BP.size());
  } else {
    memset(nRepeat,  0, sizeof(uint32)   * BP.size());
    memset(nUnique,  0, sizeof(uint32)   * BP.size());
  }

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     &frg    = tig->ufpath[fi];
    int32       frgbgn = min(frg.position.bgn, frg.position.end);
    int32       frgend = max(frg.position.bgn, frg.position.end);

    //  Search for the region that matches the read.  BP's are sorted in increasing order.  It
    //  probably doesn't matter, but makes the logging a little easier to read.

    uint32      rid = UINT32_MAX;
    bool        rpt = false;

    //fprintf(stderr, "Searching for placement for read %u at %d-%d\n", frg.ident, frgbgn, frgend);

    for (uint32 ii=0; ii<BP.size(); ii++) {
      int32   rgnbgn = BP[ii]._bgn;
      int32   rgnend = BP[ii]._end;
      bool    repeat = BP[ii]._rpt;

      //  For repeats, the read must be contained fully.

      if ((repeat == true) && (rgnbgn <= frgbgn) && (frgend <= rgnend)) {
        rid = ii;
        rpt = true;
        break;
      }

      //  For non-repeat, the read just needs to intersect.

      if ((repeat == false) && (rgnbgn < frgend) && (frgbgn < rgnend)) {
        rid = ii;
        rpt = false;
        break;
      }
    }

    if (rid == UINT32_MAX) {
      fprintf(stderr, "Failed to place read %u at %d-%d with %lu BPs:\n", frg.ident, frgbgn, frgend, BP.size());
      for (uint32 ii=0; ii<BP.size(); ii++)
        fprintf(stderr, "BP[%3u] at %8u-%8u repeat %u\n", ii, BP[ii]._bgn, BP[ii]._end, BP[ii]._rpt);
      flushLog();
    }
    assert(rid != UINT32_MAX);  //  We searched all the BP's, the read had better be placed!

    //  If moving reads, move the read!

    if (doMove) {
      if (newTigs[rid] == NULL) {
        lowCoord[rid] = frgbgn;

        newTigs[rid]  = tigs.newUnitig(true);
        newTigs[rid]->_isBubble = tig->_isBubble;

        if (nRepeat[rid] > nUnique[rid])
          newTigs[rid]->_isRepeat = true;
      }

      newTigs[rid]->addRead(frg, -lowCoord[rid], false);
    }

    //  Else, we're not moving, just count how many reads came from repeats or uniques.

    else {
      if (rpt)
        nRepeat[rid]++;
      else
        nUnique[rid]++;
    }
  }

  //  Return the number of tigs created.

  uint32  nTigsCreated = 0;

  for (uint32 ii=0; ii<BP.size(); ii++)
    if (nRepeat[ii] + nUnique[ii] > 0)
      nTigsCreated++;

  return(nTigsCreated);
}



