
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
    _beads[link]._nextOffset = beadLink;

    column->_beads[beadLink]._prevOffset = link;

  } else {
    assert(column->_beads[beadLink]._prevOffset != UINT16_MAX);
    assert(_prevColumn == column);

    _beads[link]._prevOffset = beadLink;
    _beads[link]._nextOffset = UINT16_MAX;

    column->_beads[beadLink]._nextOffset = link;
  }

  return(link);
}



//  Merges the next column into this column, if possible.  Possible if no read has
//  an actual base in both columns.

bool
abColumn::mergeWithNext(abAbacus *abacus, bool highQuality) {
  abColumn *lcolumn = this;
  abColumn *rcolumn = next();
  abColumn *ncolumn = next()->next();  //  The column after rcolumn.

  assert(lcolumn != NULL);
  assert(rcolumn != NULL);

  assert(lcolumn->next() == rcolumn);
  assert(rcolumn->prev() == lcolumn);

#if 0
  lcolumn->checkLinks();
  rcolumn->checkLinks();
#endif

  //  If both columns have a non-gap (for a single read), we cannot merge.

  for (uint32 ii=0; ii<lcolumn->_beadsLen; ii++) {
    uint32  jj = lcolumn->_beads[ii].nextOffset();

    if ((jj < UINT16_MAX) &&
        (lcolumn->_beads[ii].base() != '-') &&
        (rcolumn->_beads[jj].base() != '-'))
      return(false);
  }

#if 0
  fprintf(stderr, "MERGE columns %d %p <-  %d %p\n",
          lcolumn->position(), lcolumn,
          rcolumn->position(), rcolumn);

  fprintf(stderr, "rcolumn links\n");
  rcolumn->showLinks();

  lcolumn->checkLinks();
  rcolumn->checkLinks();
  ncolumn->checkLinks();
#endif

  //  OK to merge.  Merge all the bases from the right column to the current column.  We already
  //  checked that whenever the right column has a base, the left column has a gap, so just march
  //  down the right column and move those bases over!

  for (uint16 rr=0; rr<rcolumn->_beadsLen; rr++) {
    uint16  ll = rcolumn->_beads[rr].prevOffset();

    //  Ignore the gaps.

    if (rcolumn->_beads[rr].base() == '-')
      continue;

    //  Oh, great.  We just found the end of a read.  We need to link in the gap (in lcolumn)
    //  before we can swap.  Correction: we need to ADD a gap (in lcolumn) before we can swap.

    if (ll == UINT16_MAX) {
      //fprintf(stderr, "EXTEND READ at rr=%d\n", rr);
      ll = lcolumn->extendRead(rcolumn, rr);
    }

    //  The simple case: just swap the contents.

#if 0
    fprintf(stderr, "mergeWithNext()-- swap beads lcolumn %d %c and rcolumn %d %c\n",
            ll, lcolumn->_beads[ll].base(),
            rr, rcolumn->_beads[rr].base());
#endif

#ifdef BASECOUNT
    lcolumn->baseCountDecr(lcolumn->_beads[ll].base());
    rcolumn->baseCountDecr(rcolumn->_beads[rr].base());  //  We don't really care about rcolumn.
#endif

    swap(lcolumn->_beads[ll], rcolumn->_beads[rr]);

#ifdef BASECOUNT
    lcolumn->baseCountIncr(lcolumn->_beads[ll].base());
    rcolumn->baseCountIncr(rcolumn->_beads[rr].base());
#endif

    //  While we're here, update the bead-to-read maps.

    beadID oldb(rcolumn, rr);
    beadID newb(lcolumn, ll);

    map<beadID,uint32>::iterator  fit = abacus->fbeadToRead.find(oldb);  //  Does old bead exist
    map<beadID,uint32>::iterator  lit = abacus->lbeadToRead.find(oldb);  //  in either map?

    if (fit != abacus->fbeadToRead.end()) {
      uint32  rid = fit->second;

      //fprintf(stderr, "mergeWithNext()-- move fbeadToRead from %p/%d to %p/%d for read %d\n",
      //        rcolumn, rr, lcolumn, ll, rid);

      abacus->fbeadToRead.erase(fit);     //  Remove the old bead to read pointer

      abacus->fbeadToRead[newb] = rid;    //  Add a new bead to read pointer
      abacus->readTofBead[rid]  = newb;   //  Update the read to bead pointer
    }

    if (lit != abacus->lbeadToRead.end()) {
      uint32  rid = lit->second;

      //fprintf(stderr, "mergeWithNext()-- move lbeadToRead from %p/%d to %p/%d for read %d\n",
      //        rcolumn, rr, lcolumn, ll, rid);

      abacus->lbeadToRead.erase(lit);

      abacus->lbeadToRead[newb] = rid;
      abacus->readTolBead[rid]  = newb;
    }
  }

  //  The rcolumn should now be full of gaps.  (We could just test that baseCount('-') == depth()

  for (uint32 rr=0; rr<rcolumn->_beadsLen; rr++)
    assert(rcolumn->_beads[rr].base() == '-');

#if 0
  lcolumn->checkLinks();
  rcolumn->checkLinks();
  ncolumn->checkLinks();
#endif

  //  To make checkLinks() work, we need to unlink rcolumn from the column list right now.

  if (rcolumn->_prevColumn)   rcolumn->_prevColumn->_nextColumn = rcolumn->_nextColumn;
  if (rcolumn->_nextColumn)   rcolumn->_nextColumn->_prevColumn = rcolumn->_prevColumn;

  assert(ncolumn == rcolumn->next());
  assert(ncolumn == next());

  //  Before the rcolumn can be removed, we need to unlink it from the bead link list.  If there is a column after rcolumn,
  //  we need to move rcolumn's link pointers to lcolumn (prev) and ncolumn (next).
  //
  //  The actual example (,'s indicate no bases because the read ended):
  //
  //          1234    Column 2 is merged into column 1, and then we delete column 2.
  //    1  -aa-aaa
  //    2  gaaT,,,       Column 1 read 4 has _beads position 3.            (next = 1)
  //    3  gaaT,,,       Column 2 read 4 has _beads position 1.  (prev = 3, next = 1)
  //    4  -aa-aaa       Column 3 read 4 has _beads position 1.  (prev = 1)
  //    5  -aa-aaa     When column 2 is deleted, we're left with a busted link back from col 3 read 4; it should be 3
  //    6  gaa,,,,

  if (ncolumn != NULL) {
    for (uint32 rr=0; rr<rcolumn->_beadsLen; rr++) {
      uint16  bl = rcolumn->_beads[rr].prevOffset();  //  back link from deleted column to lcolumn -> set as ncolumns back link (known as ll above)
      uint16  fl = rcolumn->_beads[rr].nextOffset();  //  forw link from deleted column to ncolumn -> set as lcolumns forw link

      if (bl != UINT16_MAX)   lcolumn->_beads[bl]._nextOffset = fl;
      if (fl != UINT16_MAX)   ncolumn->_beads[fl]._prevOffset = bl;
    }
  }

  lcolumn->checkLinks();
  ncolumn->checkLinks();

  //  Now, finally, we're done.  Remove the old column, recall the base, and do a final check.

  //fprintf(stderr, "mergeWithNext()--  Remove rcolumn %d %p\n", rcolumn->position(), rcolumn);

  delete rcolumn;

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

#if 0
  fprintf(stderr, "mergeColumns()--\n");
  display(stderr);
#endif

  //  If we merge, update the base call, and stay here to try another merge of the now different
  //  next column.  Otherwise, we didn't merge anything, so advance to the next column.

  while (column->next()) {
    if (column->mergeWithNext(this, highQuality) == true)
      somethingMerged = true;
    else
      column = column->next();
  }

  //  If any merges were performed, refresh.  This updates the column list.

  if (somethingMerged)
    refreshColumns();
}
