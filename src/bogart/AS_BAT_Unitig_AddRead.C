
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

#include "AS_BAT_Unitig.H"




void
Unitig::addRead(ufNode node, int offset, bool report) {

  node.position.bgn += offset;
  node.position.end += offset;

  assert(node.ident > 0);

  // keep track of the unitig a read is in
  _vector->registerRead(node.ident, _id, ufpath.size());

  // keep track of max position in unitig
  int32 frgEnd = max(node.position.bgn, node.position.end);
  if (frgEnd > _length)
    _length = frgEnd;

  ufpath.push_back(node);

  if ((report) || (node.position.bgn < 0) || (node.position.end < 0)) {
    int32 trulen = RI->readLength(node.ident);
    int32 poslen = (node.position.end > node.position.bgn) ? (node.position.end - node.position.bgn) : (node.position.bgn - node.position.end);

    if (node.contained)
      writeLog("Added read %d (len %d) to unitig %d at %d,%d (idx %lu) (lendiff %d) (contained in %d)\n",
              node.ident, trulen, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              poslen - trulen,
              node.contained);
    else
      writeLog("Added read %d (len %d) to unitig %d at %d,%d (idx %lu) (lendiff %d)\n",
              node.ident, trulen, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              poslen - trulen);

    assert(poslen / trulen < 10);
    assert(trulen / poslen < 10);
  }

  assert(node.position.bgn >= 0);
  assert(node.position.end >= 0);
}



//  Percolate the last read to the correct spot in the list.
#if 0
void
Unitig::bubbleSortLastRead(void) {
  uint32   previd  = ufpath.size() - 2;
  uint32   lastid  = ufpath.size() - 1;

  ufNode   last    = ufpath[lastid];
  uint32   lastbgn = min(last.position.bgn, last.position.end);

  while ((lastid > 0) &&
         (lastbgn < min(ufpath[previd].position.bgn, ufpath[previd].position.end))) {
    ufpath[lastid] = ufpath[previd];

    _ufpathIdx[ufpath[lastid].ident] = lastid;

    lastid--;
    previd--;
  }

  _ufpathIdx[last.ident] = lastid;

  if (lastid < ufpath.size() - 1)
    ufpath[lastid] = last;
}
#endif
