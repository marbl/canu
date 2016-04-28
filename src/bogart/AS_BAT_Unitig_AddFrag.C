
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
 *    src/AS_BAT/AS_BAT_Unitig_AddFrag.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010,2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"




void
Unitig::addFrag(ufNode node, int offset, bool report) {

  node.position.bgn += offset;
  node.position.end += offset;

  assert(node.ident > 0);

  // keep track of the unitig a frag is in
  _inUnitig[node.ident]     = _id;
  _pathPosition[node.ident] = ufpath.size();

  // keep track of max position in unitig
  int32 frgEnd = MAX(node.position.bgn, node.position.end);
  if (frgEnd > _length)
    _length = frgEnd;

  ufpath.push_back(node);

  if ((report) || (node.position.bgn < 0) || (node.position.end < 0)) {
    int32 trulen = FI->fragmentLength(node.ident);
    int32 poslen = (node.position.end > node.position.bgn) ? (node.position.end - node.position.bgn) : (node.position.bgn - node.position.end);

    if (node.contained)
      writeLog("Added frag %d (len %d) to unitig %d at %d,%d (idx %lu) (lendiff %d) (contained in %d)\n",
              node.ident, trulen, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              poslen - trulen,
              node.contained);
    else
      writeLog("Added frag %d (len %d) to unitig %d at %d,%d (idx %lu) (lendiff %d)\n",
              node.ident, trulen, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              poslen - trulen);

    assert(poslen / trulen < 10);
    assert(trulen / poslen < 10);
  }

  assert(node.position.bgn >= 0);
  assert(node.position.end >= 0);
}



//  Percolate the last fragment to the correct spot in the list.
#if 0
void
Unitig::bubbleSortLastFrag(void) {
  uint32   previd  = ufpath.size() - 2;
  uint32   lastid  = ufpath.size() - 1;

  ufNode   last    = ufpath[lastid];
  uint32   lastbgn = MIN(last.position.bgn, last.position.end);

  while ((lastid > 0) &&
         (lastbgn < MIN(ufpath[previd].position.bgn, ufpath[previd].position.end))) {
    ufpath[lastid] = ufpath[previd];

    _pathPosition[ufpath[lastid].ident] = lastid;

    lastid--;
    previd--;
  }

  _pathPosition[last.ident] = lastid;

  if (lastid < ufpath.size() - 1)
    ufpath[lastid] = last;
}
#endif
