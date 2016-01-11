
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

//  (Optionally) recall consensus in all columns.
//  Rebuild the list of columns in the multialign.
//  Rebuild the column to position map.

void
abAbacus::refreshColumns(void) {

  //fprintf(stderr, "abAbacus::refreshColumns()--\n");

  //  Given that _firstColumn is a valid column, walk to the start of the column list.

  while (_firstColumn->_prevColumn != NULL)
    _firstColumn = _firstColumn->_prevColumn;

  //  Number the columns, so we can make sure the _columns array has enough space.  Probably not
  //  needed to be done first, but avoids having the resize call in the next loop.

  uint32 cn = 0;

  for (abColumn *column = _firstColumn; column; column = column->next())
    column->_columnPosition = cn++;  //  Position of the column in the gapped consensus.

  //  Fake out resizeArray so it will work on three arrays.

  uint32  cm = _columnsMax;

  resizeArray(_columns,  0, cm, cn+1, resizeArray_doNothing);  cm = _columnsMax;
  resizeArray(_cnsBases, 0, cm, cn+1, resizeArray_doNothing);  cm = _columnsMax;
  resizeArray(_cnsQuals, 0, cm, cn+1, resizeArray_doNothing);  _columnsMax = cm;

  //  Build the list of columns and update consensus and quals while we're there.

  _columnsLen = 0;

  for (abColumn *column = _firstColumn; column; column = column->next()) {
    _columns [_columnsLen] = column;
    _cnsBases[_columnsLen] = column->baseCall();
    _cnsQuals[_columnsLen] = column->baseQual();
    _columnsLen++;
  }

  _cnsBases[_columnsLen] = 0;
  _cnsQuals[_columnsLen] = 0;  //  Not actually zero terminated.

  //for (abColumn *column = _firstColumn; column; column = column->next())
  //  fprintf(stderr, "refreshColumns()--  column %p is at position %d\n",
  //          column, column->position());
}


void
abAbacus::recallBases(bool highQuality) {

  //fprintf(stderr, "abAbacus::recallBases()--  highQuality=%d\n", highQuality);

  //  Given that _firstColumn is a valid column, walk to the start of the column list.
  //  We could use _columns[] instead.

  while (_firstColumn->_prevColumn != NULL)
    _firstColumn = _firstColumn->_prevColumn;

  //  Number the columns, so we can make sure the _columns array has enough space.  Probably not
  //  needed to be done first, but avoids having the resize call in the next loop.

  for (abColumn *column = _firstColumn; column; column = column->next())
    column->baseCall(highQuality);

  //  After calling bases, we need to refresh to copy the bases from each column into
  //  _cnsBases and _cnsQuals.

  refreshColumns();
}
