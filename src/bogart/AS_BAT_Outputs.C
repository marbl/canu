
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
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_Outputs.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-MAR-31
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUN-05
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-04
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-MAR-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "tgStore.H"


void
unitigToTig(tgTig       *tig,
            uint32       tigid,
            Unitig      *utg) {

  //  Initialize the output tig.

  tig->clear();

  tig->_tigID           = tigid;
  utg->_tigID           = tigid;

  tig->_coverageStat    = 1.0;  //  Default to just barely unique
  tig->_microhetProb    = 1.0;  //  Default to 100% probability of unique

  //  Set the class.

  if      (utg->_isUnassembled == true)
    tig->_class = tgTig_unassembled;

  else if (utg->_isBubble == true)
    tig->_class = tgTig_bubble;

  else
    tig->_class = tgTig_contig;

  tig->_suggestRepeat   = (utg->_isRepeat   == true);
  tig->_suggestCircular = (utg->_isCircular == true);

  tig->_layoutLen       = utg->getLength();

  //  Transfer reads from the bogart tig to the output tig.

  resizeArray(tig->_children, tig->_childrenLen, tig->_childrenMax, utg->ufpath.size(), resizeArray_doNothing);

  for (uint32 ti=0; ti<utg->ufpath.size(); ti++) {
    ufNode        *frg   = &utg->ufpath[ti];

    tig->addChild()->set(frg->ident,
                         frg->parent, frg->ahang, frg->bhang,
                         frg->position.bgn, frg->position.end);
  }
}



void
writeUnitigsToStore(UnitigVector  &unitigs,
                    char          *fileprefix,
                    char          *tigStorePath,
                    uint32         frg_count_target,
                    bool           isFinal) {
  uint32      utg_count              = 0;
  uint32      frg_count              = 0;
  uint32      prt_count              = 1;
  char        filename[FILENAME_MAX] = {0};

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

  tgStore     *tigStore = new tgStore(tigStorePath);
  tgTig       *tig      = new tgTig;

  for (uint32 tigID=0, ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if ((utg == NULL) || (utg->getNumFrags() == 0))
      continue;

    assert(utg->getLength() > 0);

    //  Convert the bogart tig to a tgTig and save to the store.

    unitigToTig(tig, (isFinal) ? tigID : ti, utg);
    tigID++;

    tigStore->insertTig(tig, false);

    //  Increment the partition if the current one is too large.

    if ((frg_count + utg->getNumFrags() >= frg_count_target) &&
        (frg_count                      >  0)) {
      fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
              prt_count, utg_count, frg_count);

      prt_count++;
      utg_count = 0;
      frg_count = 0;
    }

    //  Note that the tig is included in this partition.

    utg_count += 1;
    frg_count += utg->getNumFrags();

    //  Map the tig to a partition, and log both the tig-to-partition map and the partition-to-read map.

    fprintf(iidm, "bogart "F_U32" -> tig "F_U32" (in partition "F_U32" with "F_U32" frags)\n",
            utg->id(),
            utg->tigID(),
            prt_count,
            utg->getNumFrags());

    for (uint32 fragIdx=0; fragIdx<utg->getNumFrags(); fragIdx++)
      fprintf(part, "%d\t%d\n", prt_count, utg->ufpath[fragIdx].ident);
  }

  fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",   //  Don't forget to log the last partition!
          prt_count, utg_count, frg_count);

  fclose(pari);
  fclose(part);
  fclose(iidm);

  delete    tig;
  delete    tigStore;
}



class  rawEdge_t {
public:
  rawEdge_t(uint32 o, uint32 t, int32 ab, int32 ae, int32 bb, int32 be) {
    oi     = o;
    tigID  = t;

    Abgn   = ab;
    Aend   = ae;

    Bbgn   = bb;
    Bend   = be;
  };

  uint32   oi;
  int32    tigID;

  int32    Abgn;  //  Overlapping read placement.
  int32    Aend;

  int32    Bbgn;  //  Parent placement.
  int32    Bend;

  bool  operator<(rawEdge_t const &that) const {
    if (tigID != that.tigID)
      return(tigID < that.tigID);

    return(Abgn < that.Abgn);
  }
};



void
findUnusedEdges(UnitigVector  &unitigs,
                ufNode        *rdA,         //  Read we're finding edges for
                bool           rdA3p,       //  Overlaps from the 3' end of the read
                set<uint32>    edgeReads,
                FILE          *EF) {

  uint32      rdAid    =  rdA->ident;
  uint32      rdAlen   =  FI->fragmentLength(rdAid);
  bool        rdAfwd   =  rdA->isForward();
  int32       rdAlo    =  (rdAfwd) ? (rdA->position.bgn) : (rdA->position.end);
  int32       rdAhi    =  (rdAfwd) ? (rdA->position.end) : (rdA->position.bgn);
  uint32      rdAtigID =  Unitig::fragIn(rdAid);
  Unitig     *rdAtig   =  unitigs[rdAtigID];

  uint32      ovlLen   = 0;
  BAToverlap *ovl      = OC->getOverlaps(rdA->ident, AS_MAX_ERATE, ovlLen);

  vector<rawEdge_t>   rawEdges;

  //fprintf(stderr, "WORKING ON read rdA=%u 3p=%d\n", rdA->ident, rdA3p);

  //  Over all overlaps for this read, find and report edges to 'edgeReads'.  Though
  //  edgeReads should be just one read per tig end, the code below was originally written
  //  to find all edges to all reads, then pick the longest for each cluster.

  for (uint32 oi=0; oi<ovlLen; oi++) {
    if ((ovl[oi].AisContainer()) ||          //  Not interested in container overlaps.
        (ovl[oi].AisContained()) ||          //  Allow A-is-contained overlaps?  Should be OK, but only really care about dovetails.
        (ovl[oi].AEndIs3prime() != rdA3p))   //  Overlap off the wrong end of A.
      continue;

    uint32    rdBid    = ovl[oi].b_iid;
    uint32    rdBtigID = Unitig::fragIn(rdBid);
    Unitig   *rdBtig   = unitigs[rdBtigID];

    if ((rdBtig == NULL) ||
        (rdBtig->getNumFrags() == 0) ||    //  Not interested in edges to singletons
        (rdBtig->_isUnassembled == true))  //  Or other unassembled crap.  rdA filtered outside here.
      continue;

    if ((rdAtigID != rdBtigID) &&          //  Not to self (circular) and
        (edgeReads.count(rdBid) == 0))     //  not a read we can overlap to.
      continue;

    ufNode   *rdB      = &rdBtig->ufpath[ Unitig::pathPosition(rdBid) ];
    bool      rdBfwd   =  rdB->isForward();
    int32     rdBlo    = (rdBfwd) ? (rdB->position.bgn) : (rdB->position.end);
    int32     rdBhi    = (rdBfwd) ? (rdB->position.end) : (rdB->position.bgn);

    //  Exclude overlaps satisfied in the same tig.

    if ((rdAtigID == rdBtigID) && (rdAlo < rdBhi) && (rdBlo < rdAhi))
      continue;

    //  Exclude overlaps that are higher than expected error.

    ;

    //  Compute the placement of rdA on rdBtig.

    ufNode           placed;
    BestEdgeOverlap  edge(ovl[oi]);

    rdBtig->placeFrag(placed,
                      rdAid,
                      rdA3p,
                      &edge);

    //writeLog("placed tig %u rdA %u %d-%d on tig %u %d-%d from rdB %u %d-%d oi %u\n",
    //         rdAtigID, rdAid, rdAlo, rdAhi, rdBtigID, placed.position.bgn, placed.position.end, rdBid, rdBlo, rdBhi, oi);

    //  Save the overlap.

    rawEdges.push_back(rawEdge_t(oi, rdBtigID, placed.position.min(), placed.position.max(), rdBlo, rdBhi));
  }

  //  We've now got a pile of (unsorted) overlaps to reads in other tigs.  We need to pick one
  //  overlap (the longest?) from each pile and output it.

  sort(rawEdges.begin(), rawEdges.end());

  //  We expect to have a pile of placements that are 'the same', generated by each one of the
  //  overlapping reads in the target tig.  We need to group these placements together and pick
  //  one exemplar overlap to output the edge for.
  //
  //  A complication is caused by large tandem repeats.  We can get two distinct placements that
  //  overlap:
  //
  //     [rrrr][rrrr]
  //     ---------------           (rdA aligning to the first and second repeat)
  //           ----------------    (rdA aligning to only the second repeat)
  //
  //  These are just overlaps, and we don't know that the rest of rdA fails to align.
  //
  //  Overlaps are sorted by the start of rdA on rdBtig.  We'll use the simple and largely unvalidated
  //  heuristic of any placement that starts within 500bp of the last is for the same placement.

  for (uint32 ri=0, rj=0; ri<rawEdges.size(); ri = rj) {
    for (rj=ri+1; ((rj < rawEdges.size()) &&
                   (rawEdges[rj].tigID == rawEdges[ri].tigID) &&
                   (rawEdges[rj-1].Abgn + 500 >= rawEdges[rj].Abgn)); )
      rj++;

    //  Scan overlaps from ri to rj, retain the thickest.

    //fprintf(stderr, "Scan batch from ri=%u to rj=%u\n", ri, rj);

    uint32   rrMax = 0;
    int32    rrIdx = INT32_MAX;

    for (uint32 rr=ri; rr<rj; rr++) {
      int32  olapLen = 0;

      if (rawEdges[rr].Abgn < rawEdges[rr].Bbgn) {
        assert(rawEdges[rr].Bend >= rawEdges[rr].Abgn);
        olapLen = rawEdges[rr].Bend - rawEdges[rr].Abgn;
      } else {
        assert(rawEdges[rr].Aend >= rawEdges[rr].Bbgn);
        olapLen = rawEdges[rr].Aend - rawEdges[rr].Bbgn;
      }

      if (rrMax < olapLen) {
        rrMax = olapLen;
        rrIdx = rr;
      }
    }

    //  Emit the edge.

    uint32 oi = rawEdges[rrIdx].oi;

    uint32    rdBid    = ovl[oi].b_iid;
    uint32    rdBtigID = Unitig::fragIn(rdBid);
    Unitig   *rdBtig   = unitigs[rdBtigID];

    ufNode   *rdB      = &rdBtig->ufpath[ Unitig::pathPosition(rdBid) ];
    bool      rdBfwd   =  rdB->isForward();
    int32     rdBlo    = (rdBfwd) ? (rdB->position.bgn) : (rdB->position.end);
    int32     rdBhi    = (rdBfwd) ? (rdB->position.end) : (rdB->position.bgn);

    char  rdAEnd, rdBEnd;

    if (ovl[oi].isDovetail()) {
      rdAEnd = ovl[oi].AEndIs5prime() ? '5' : '3';
      rdBEnd = ovl[oi].BEndIs5prime() ? '5' : '3';
    } else {
      rdAEnd = ovl[oi].AisContainer() ? 'C' : 'c';
      rdBEnd = ovl[oi].AisContainer() ? 'c' : 'C';
    }

    char  ori = (ovl[oi].flipped) ? '<' : '>';

    fprintf(EF, "tig %7u %c read %8u at %9u %-9u %8d %c %-8d tig %7u %c read %8u at %9u %-9u\n",
            rdAtig->_tigID, rdAtig->type(), rdAid, rdA->position.bgn, rdA->position.end,
            ovl[oi].a_hang, ori, ovl[oi].b_hang,
            rdBtig->_tigID, rdBtig->type(), rdBid, rdB->position.bgn, rdB->position.end);
  }
}




void
writeUnusedEdges(UnitigVector  &unitigs,
                 char          *fileprefix) {
  char        filename[FILENAME_MAX] = {0};

  sprintf(filename, "%s.unused.edges", fileprefix);
  FILE *EF = fopen(filename, "w");
  if (errno)
    fprintf(stderr, "Failed to create unused edge output '%s': %s\n", filename, strerror(errno)), exit(1);

  //  Find reads we're allowed to find edges to.  We can pick either the outer-most non-contained reads,
  //  or just the reads touching the edge of the tig.

  set<uint32>  edgeReads;  //  Reads at the end of the tig
  set<uint32>  nearReads;  //  Reads close to the end of the tig


  //  Find the outer-most non-contained reads in each unitig.

#if 0
  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) ||
        (tig->getNumFrags() == 0) ||
        (tig->_isUnassembled == true))
      continue;

    //  Find reads at the start of the tig

    for (uint32 ct=0, fi=0; (ct < 5) && (fi < tig->ufpath.size()); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      if (OG->isContained(frg->ident) == false) {
        if (ct == 0)
          edgeReads5e.insert(frg->ident);
        else
          nearReads5m.insert(frg->ident);
        ct++;
      }
    }

    //  Find reads at the end of the tig

    for (uint32 ct=0, fi=tig->ufpath.size(); (ct < 5) && (fi-- > 0); ) {
      ufNode  *frg = &tig->ufpath[fi];

      if (OG->isContained(frg->ident) == false) {
        if (ct == 0)
          edgeReads3e.insert(frg->ident);
        else
          nearReads3m.insert(frg->ident);
        ct++;
      }
    }
  }
#endif

  //  Find the reads at the ends of the tig.

#if 1
  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) ||
        (tig->getNumFrags() == 0) ||
        (tig->_isUnassembled == true))
      continue;

    edgeReads.insert(tig->firstRead()->ident);
    edgeReads.insert(tig->lastRead()->ident);
  }
#endif


  //  Step through all the unitigs, find all unused overlaps off the ends of the tig.


  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) ||
        (tig->getNumFrags() == 0) ||
        (tig->_isUnassembled == true))
      continue;

    assert(tig->getLength() > 0);

    //  Find the first/last non-contained reads in the tig.

#if 0
    ufNode  *rd5 = &tig->ufpath.front();
    ufNode  *rd3 = &tig->ufpath.back();

    for (uint32 fi=1; (fi < tig->ufpath.size()) && (OG->isContained(rd5->ident) == true); fi++)
      rd5 = &tig->ufpath[fi];

    for (uint32 fi=tig->ufpath.size()-1; (fi-- > 0) && (OG->isContained(rd3->ident) == true); )
      rd3 = &tig->ufpath[fi];

    //  What to do if either of those reads are contained?  If so (then both will be contained; no
    //  dovetail at all) we've swapped the meaning of 5' and 3'.

    if ((OG->isContained(rd5) == true) || (OG->isContained(rd3) == true)) {
      rd5 = &tig->ufpath.front();
      rd3 = &tig->ufpath.back();
    }
#endif

    //  Find the smallest/largest read position - the two reads that are at the end of the tig.

#if 1
    ufNode  *rd5 = tig->firstRead();
    ufNode  *rd3 = tig->lastRead();
#endif

    //  Finally, we probably should be finding just the reads touching the ends of the unitig, not the
    //  first/last non-contained read.

    findUnusedEdges(unitigs, rd5, rd5->isReverse(), edgeReads, EF);   //  First read, if reverse, find edges off 3' end
    findUnusedEdges(unitigs, rd3, rd3->isForward(), edgeReads, EF);   //  Last  read, if forward, find edges off 3' end
  }

  fclose(EF);
}

