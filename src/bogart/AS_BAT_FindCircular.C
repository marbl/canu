
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

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_FindCircular.H"


class readCheck {
public:
  void         init(ufNode *node) {
    read    = node;

    ovlLen  = 0;
    ovl     = OC->getOverlaps(read->ident, ovlLen);
  };

  ufNode      *read;

  uint32       ovlLen;
  BAToverlap  *ovl;
};


//  Checks if all the overlaps implied by circularizing the tig actually
//  exist.
//
//  This checks evry read in the path
//      rdFid .. the end of the tig --circular-> the other end .. rdTid
//  and ensures that the tested read has an overlap to all later reads.
//
//  [OLD COMMENT - names are incorrect, but the idea is the same]
//
//  Decide if the external edge in externalSco can plausibly circularize the tig.
//  For this to happen, the external egdge read must overlap all the reads in the path
//  to rdA.
//
//         =======================         <- another copy of the extern read
//          -----------------------        <- another copy of the last read
//            1--------------------        <- the first read
//                 2---------------------  <- rdA
//                    3-------------------------
//                       . . . . . . . . . . . . . . .
//                            8==============================  <- extern
//                             9------------------------------ <- the last read
//                                -------------    ------
// 
//  Here, the extern read is causing rdA to look confused on it's low end.
//  In order for this to circularize, we need to check that all the
//  (non-contained) reads from the extern read to the last read have overlaps
//  to all the reads between itself and rdA (going around the circle).
//
//  This will correcly not find a circle if the extern read is too far from
//  the end, or if rdA is too far from the beginning.
//

bool
isCircularizingEdge(Unitig   *tig,
                    uint32    rdFid,     //  read From
                    uint32    rdTid) {   //  read To
  bool   isCircular = true;

  if ((rdTid == 0) || (tig->inUnitig(rdTid) != tig->id()) ||   //  Quick sanity checks:
      (rdFid == 0) || (tig->inUnitig(rdFid) != tig->id()))     //    the reads must exist.
    return(false);                                             //    the reads must be in the supplied tig.

  //  Based on the position of the reads in the tig, decide if we are
  //  checking overlaps from end-to-start or from start-to-end.  Then build a
  //  list of the reads we need to check.
  //
  //  forward == true if we iterate:
  //    from the F read to the high end of the tig
  //                    across the circle
  //                    to the low end of the tig
  //                    and up to the T read.
  //
  //  forward == false if we iterate:
  //    from the F read to the low end of the tig
  //                    across the circle
  //                    to the high end of the tig
  //                    and down to the T read.

  uint32      rdTidx  = tig->ufpathIdx(rdTid);
  uint32      rdFidx  = tig->ufpathIdx(rdFid);
  bool        forward = false;

  uint32      rcMax   = 0;
  uint32      rcLen   = 0;
  readCheck  *rc      = nullptr;

  //  -- rdT is at the start of the tig, rdF is at the end.
  if (rdTidx < rdFidx) {
    rcMax   = tig->ufpath.size() - rdFidx + 1 + rdTidx + 1;
    rc      = new readCheck [rcMax];
    forward = true;

    for (uint32 rd=rdFidx; rd < tig->ufpath.size(); rd++)
      if (OG->isContained(tig->ufpath[rd].ident) == false)
        rc[rcLen++].init(&tig->ufpath[rd]);

    for (uint32 rd=0; rd <= rdTidx; rd++)
      if (OG->isContained(tig->ufpath[rd].ident) == false)
        rc[rcLen++].init(&tig->ufpath[rd]);
  }

  //  -- rdF is at the start of the tig, rdT is at the end.
  else {
    rcMax   = tig->ufpath.size() - rdTidx + 1 + rdFidx + 1;
    rc      = new readCheck [rcMax];
    forward = false;

    for (uint32 rd=rdFidx+1; rd-- > 0; )
      if (OG->isContained(tig->ufpath[rd].ident) == false)
        rc[rcLen++].init(&tig->ufpath[rd]);

    for (uint32 rd=tig->ufpath.size(); rd-- > rdTidx; )
      if (OG->isContained(tig->ufpath[rd].ident) == false)
        rc[rcLen++].init(&tig->ufpath[rd]);
  }

  assert(rcLen < rcMax);

  //  Iterate over all the reads, checking for an overlap to all the later
  //  reads in the list.
  //
  //  Once we get a pair, decide which read ends the overlap should be on.
  //
  //  In the 'forward' case, we start with rdA being the read near the high
  //  end of the tig (passed in as rdF) and check for overlaps to each rdB
  //  'after' it in the circular tig.
  //
  //    rc[2]  ----------
  //    rc[3]     ------------ <- B read
  //                                A read -> --------->       rc[0]
  //                                             <-----------  rc[1]
  //
  //  In the 'reverse' case, we start with rdA being the read near the low
  //  end of the tig (passed in as rdF) and check for overlaps to each rdB
  //  'before; it in the circular tig.
  //
  //    rc[1]  --------->
  //    rc[0]     -----------> <- A read
  //                                B read -> ----------       rc[3]
  //                                             ------------  rc[2]
  //
  //  In both cases, we check for an overlap between rc[i] and rc[j], i<j.
  //  'forward' or 'reverse' is used to decide which end of the read we need
  //  to find the overlap on.  The orientation on the read in the picture is
  //  to help visualizing ovlA3 and ovlB3 below).
  //
  //  For each rc[] pair, we just scan all the overlaps for rdA and check that
  //    the overlap is to the correct rdB read
  //    the overlap is on the correct ends of each read
  //    the oveerlap is of good quality.

  for (uint32 ai=0; (isCircular == true) && (ai<rcLen); ai++) {
    ufNode     *rdA    = rc[ai].read;
    uint32      ovlLen = rc[ai].ovlLen;
    BAToverlap *ovl    = rc[ai].ovl;
    uint32      rdAid  = rdA->ident;

    for (uint32 bi=ai+1; (isCircular == true) && (bi<rcLen); bi++) {
      ufNode   *rdB    = rc[bi].read;
      uint32    rdBid  = rdB->ident;

      bool  ovlA3 = (forward) ? (rdA->isForward()) : (rdA->isReverse());   //  Need overlap on 3' end of A read.
      bool  ovlB3 = (forward) ? (rdB->isReverse()) : (rdB->isForward());
      bool  found = false;

      for (uint32 oo=0; (oo<ovlLen) && (ovl[oo].b_iid <= rdBid) && (found == false); oo++) {
        if ((ovl[oo].b_iid == rdBid) &&
            (ovl[oo].isDovetail() == true) &&
            (ovl[oo].AEndIs3prime() == ovlA3) &&
            (ovl[oo].BEndIs3prime() == ovlB3) &&
            (ovl[oo].erate() <= OG->reportErrorLimit())) {
          found = true;
          writeLog("                    overlap from read %9u %s idx=%-7u to read %9u %s idx=%-7u at %6.4f%% <= %6.4f%%\n",
                   rdAid, (rdA->isForward()) ? "->" : "<-", tig->ufpathIdx(rdAid),
                   rdBid, (rdB->isForward()) ? "->" : "<-", tig->ufpathIdx(rdBid), 100.0 * ovl[oo].erate(), 100.0 * OG->reportErrorLimit());
        }
      }

      //  If no overlap found, we cannot circularize this tig.
      if (found == false) {
        writeLog("                    overlap from read %9u %s idx=%-7u to read %9u %s idx=%-7u not found\n",
                 rdAid, (rdA->isForward()) ? "->" : "<-", tig->ufpathIdx(rdAid),
                 rdBid, (rdB->isForward()) ? "->" : "<-", tig->ufpathIdx(rdBid));
        isCircular = false;
      }
    }
  }

  return(isCircular);
}



void
findCircularContigs(TigVector &tigs,
                    const char *prefix) {
  FILE *CIRC = AS_UTL_openOutputFile(prefix, '.', "circular");

  fprintf(CIRC, "          ----------------------------------------------------------------------  ----------------------------------------------------------------------\n");
  fprintf(CIRC, "              first                         overlap                                    last                         overlap\n");
  fprintf(CIRC, "   tigID     readID        bgn-end           readID        bgn-end        length     readID        bgn-end           readID        bgn-end        length\n");
  fprintf(CIRC, "--------  --------- ---------- ---------- --------- ---------- ---------- ------  --------- ---------- ---------- --------- ---------- ---------- ------\n");

  for (uint32 ti=1; ti<tigs.size(); ti++) {
    Unitig  *tig = tigs[ti];

    if ((tig == NULL) ||
        (tig->_isUnassembled == true))
      continue;

    //  Grab the first and last reads in the tig, then find the edge that
    //  points out of the tig.
  
    ufNode          *fRead = tig->firstRead();
    ufNode          *lRead = tig->lastRead();

    uint32      circularLength = 0;
    uint32      ovlLen = 0;
    BAToverlap *ovl    = OC->getOverlaps(lRead->ident, ovlLen);

    for (uint32 oo=0; oo<ovlLen; oo++) {
      if ((ovl[oo].b_iid == fRead->ident) &&
          (((lRead->position.isForward() ==  true) && (ovl[oo].AEndIs3prime() ==  true)) ||
           ((lRead->position.isForward() == false) && (ovl[oo].AEndIs5prime() ==  true))))

        //  Circular!
        circularLength = RI->overlapLength(lRead->ident, fRead->ident, ovl[oo].a_hang, ovl[oo].b_hang);
    }

    if (circularLength == 0) 
       continue;

    tig->_isCircular     = true;
    tig->_circularLength = circularLength;

    ufNode *foRead = &tig->ufpath[ tig->ufpathIdx(lRead->ident) ];
    ufNode *loRead = &tig->ufpath[ tig->ufpathIdx(fRead->ident) ];

    fprintf(CIRC, "%8u  %9u %10u-%-10u %9u %10u-%-10u %6u  %9u %10u-%-10u %9u %10u-%-10u %6u\n",
            ti,
            fRead->ident, fRead->position.bgn, fRead->position.end, foRead->ident, foRead->position.bgn, foRead->position.end, tig->_circularLength,
            lRead->ident, lRead->position.bgn, lRead->position.end, loRead->ident, loRead->position.bgn, loRead->position.end, tig->_circularLength);
  }

  AS_UTL_closeFile(CIRC, prefix, '.', "circular");
}
