
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
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/AS_CNS/RefreshMANode.C
 *    src/AS_CNS/RefreshMANode.c
 *    src/utgcns/libcns/RefreshMANode.C
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2005-MAR-30 to 2008-FEB-13
 *      are Copyright 2005-2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Gennady Denisov from 2005-MAY-09 to 2008-JUN-06
 *      are Copyright 2005-2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-AUG-01
 *      are Copyright 2005-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2005-SEP-29 to 2006-OCT-03
 *      are Copyright 2005-2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-FEB-27 to 2009-SEP-25
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz beginning on 2014-NOV-17
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static char *rcsid = "$Id$";

#include "abAbacus.H"

//  Recall consensus in all columns.
//  Rebuild the list of columns in the multialign.
//  Rebuild the column to position map.

void
abAbacus::refreshMultiAlign(abMultiAlignID  mid,
                            bool            recallBase,         //  (false) If true, recall the base
                            bool            highQuality) {      //  (false) If true, use the high quality base call algorithm

  abMultiAlign *ma = getMultiAlign(mid);

  //fprintf(stderr, "abMultiAlign::refreshMultiAlign()--  legnth %u recallBase=%d highQuality=%d\n",
  //        ma->length(), recallBase, highQuality);

  ma->columns().clear();

  //  The first column MUST be correct, but we can then update everything else.
  abColID  cid = ma->firstColumn();

  for (uint32 index=0; (cid.isValid() == true); index++) {
    abColumn  *column = getColumn(cid);

    if (recallBase == true)
      baseCall(cid, highQuality);

    column->ma_position = index;                      //  Position of the column in the gapped consensus.
    ma->columns().push_back(column->ident());         //  Just a vector of columns, for random access.
    ma->last = cid;                                   //  Wherever it ends, it ends.

    assert(ma->columns()[column->ma_position] == column->ident());

    cid = column->next;
  }

  //  Check column pointers - first/last have no prev/next, and all the other prev/next agree.

  assert(getColumn(ma->firstColumn())->prevID() == abColID());
  assert(getColumn(ma->lastColumn())->nextID()  == abColID());

  abColID   pid = ma->firstColumn();  //  Previous column
  abColumn *pol = getColumn(pid);

  abColID   nid = pol->next;          //  Next column
  abColumn *nol = getColumn(nid);

  while (nid.isValid()) {
    assert(pol->next == nid);  //  We're iterating over this, must be true.
    assert(pid == nol->prev);  //  What we're testing

    pid = nid;
    pol = nol;

    nid = pol->next;          //  Next column
    nol = getColumn(nid);
  }
}
