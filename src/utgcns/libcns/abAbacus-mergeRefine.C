
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
 *    src/AS_CNS/MergeRefine.C
 *    src/AS_CNS/MergeRefine.c
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/MergeRefine.C
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
 *    Sergey Koren from 2008-FEB-27 to 2009-MAY-14
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-MAR-06
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-14
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static char *rcsid = "$Id$";

#include "abAbacus.H"



//  Extends the read represented by column/beadLink into this column.

uint16
abColumn::extendRead(abColumn *column, uint16 beadLink) {

  increaseArray(_beads, _beadsLen, _beadsMax, 1);

  uint32  link = _beadsLen++;

  _beads[link]._unused     = column->_beads[beadLink]._unused;
  _beads[link]._isRead     = column->_beads[beadLink]._isRead;
  _beads[link]._isUnitig   = column->_beads[beadLink]._isUnitig;
  _beads[link]._base       = '-';
  _beads[link]._qual       = 0;

  if (column->_beads[beadLink]._prevOffset == UINT16_MAX) {
    assert(column->_beads[beadLink]._nextOffset != UINT16_MAX);
    assert(_nextColumn == column);

    _beads[link]._prevOffset = UINT16_MAX;
    _beads[link]._thisOffset = link;
    _beads[link]._nextOffset = beadLink;

    column->_beads[beadLink]._prevOffset = link;

  } else {
    assert(column->_beads[beadLink]._prevOffset != UINT16_MAX);
    assert(_prevColumn == column);

    _beads[link]._prevOffset = beadLink;
    _beads[link]._thisOffset = link;
    _beads[link]._nextOffset = UINT16_MAX;

    column->_beads[beadLink]._nextOffset = link;
  }

  return(link);
}



//  Merges the next column into this column, if possible.  Possible if no read has
//  an actual base in both columns.

bool
abColumn::mergeWithNext(bool highQuality) {
  abColumn *lcolumn = this;
  abColumn *rcolumn = next();

  assert(lcolumn != NULL);
  assert(rcolumn != NULL);

  assert(lcolumn->next() == rcolumn);
  assert(rcolumn->prev() == lcolumn);

  //  If both columns have a non-gap (for a single read), we cannot merge.

  for (uint32 ii=0; ii<lcolumn->_beadsLen; ii++) {
    uint32  jj = lcolumn->_beads[ii].nextOffset();

    if ((jj < UINT16_MAX) &&
        (lcolumn->_beads[ii].base() != '-') &&
        (rcolumn->_beads[jj].base() != '-'))
      return(false);
  }

  //  OK to merge.  Merge all the bases from the right column to the current column.  We already
  //  checked that whenever the right column has a base, the left column has a gap, so just march
  //  down the right column and move those bases over!

  for (uint32 rr=0; rr<rcolumn->_beadsLen; rr++) {
    uint32  ll = rcolumn->_beads[rr].prevOffset();

    //  Ignore the gaps.

    if (rcolumn->_beads[rr].base() == '-')
      continue;

    //  Oh, great.  We just found the end of a read.  We need to link in the gap (in lcolumn)
    //  before we can swap.  Correction: we need to ADD a gap (in lcolumn) before we can swap.

    if (ll == UINT16_MAX)
      ll = lcolumn->extendRead(rcolumn, rr);

    //  The simple case, just swap the contents.

    swap(lcolumn->_beads[ll], rcolumn->_beads[rr]);
  }

  //  Check that all bases in the rcolumn are now empty.

  for (uint32 rr=0; rr<rcolumn->_beadsLen; rr++)
    assert(rcolumn->_beads[rr].base() == '-');

  //  Yank the column out of the list and dispose of it properly.  Please don't litter!

  if (rcolumn->_prevColumn)   rcolumn->_prevColumn->_nextColumn = rcolumn->_nextColumn;
  if (rcolumn->_nextColumn)   rcolumn->_nextColumn->_prevColumn = rcolumn->_prevColumn;

  delete rcolumn;

  //  Finally, recall the base for this column.

  baseCall(highQuality);

  return(true);
}




//  Simple sweep through the MultiAlignment columns, looking for columns to merge and removing null
//  columns
//
//  Loop over all columns.  If we do not merge, advance to the next column, otherwise, stay here and
//  merge to the now different next column (mergeCompatible removes the column that gets merged into
//  the current column).
//
//  Note that _firstColumn is never removed.  The second column could be merged into the first,
//  and the second one then removed.
void
abAbacus::mergeColumns(bool highQuality) {
  assert(_firstColumn != NULL);

  abColumn   *column = _firstColumn;

  bool        somethingMerged = false;

  assert(column->prev() == NULL);

  //  If we merge, update the base call, and stay here to try another merge of the now different
  //  next column.  Otherwise, we didn't merge anything, so advance to the next column.

  while (column->next()) {
    if (column->mergeWithNext(highQuality) == true)
      somethingMerged = true;
    else
      column = column->next();
  }

  //  If any merges were performed, refresh.  This updates the column list.

  if (somethingMerged)
    refreshColumns();
}
