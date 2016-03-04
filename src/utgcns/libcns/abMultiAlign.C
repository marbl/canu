
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
 *    src/AS_CNS/MultiAlignment_CNS.C
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/MultiAlignment_CNS.C
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
 *    Brian P. Walenz from 2014-NOV-17 to 2015-AUG-11
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-01
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "abAbacus.H"


void
abAbacus::getConsensus(tgTig *tig) {

  //  Resize the bases/quals storage to include a NUL byte.

  resizeArrayPair(tig->_gappedBases, tig->_gappedQuals, 0, tig->_gappedMax, _columnsLen + 1, resizeArray_doNothing);

  //  Copy in the bases.

  memcpy(tig->_gappedBases, _cnsBases, sizeof(char)  * _columnsLen);
  memcpy(tig->_gappedQuals, _cnsQuals, sizeof(uint8) * _columnsLen);

  //  Terminate the strings.

  tig->_gappedBases[_columnsLen] = 0;
  tig->_gappedQuals[_columnsLen] = 0;

  //  And set the length.

  tig->_gappedLen = _columnsLen;
}





//  Crap.  How can I figure out where a read starts in the multialign??
//  Something needs to store the first column and index for the first bead....but
//  the column could be merged and removed.
//
//  Put it in the column?  A list of:
//    uint32  readId - the read id that starts here
//    uint16  link   - the bead for that read
//  Adds one pointer to each column, most pointers are null.  If the pointer is set, and we move bases,
//  we need to update.  On a 32m unitig (pretty big) this adds 256mb + 1mb (100k reads at 8 bytes each)
//
//  Put it in the read?  If we move the first bead in the read, we know we need to update something,
//  but we have no link back to the read from the bead.
//
//  It's so sparse - can we just hash it?  Need to store:
//    readID      -> column/link of first bead
//    column/link -> readID
//  The column probably needs to be a pointer (not the columnid) since the ID changes frequently.
//  Moving beads needs to update this structure, but any structure we make has the same problem.



  //  index < length eliminates any endgaps from the delta list KAR, 09/19/02


uint32
abAbacus::getSequenceDeltas(uint32    sid,
                            int32    *deltas) {
  abSequence         *seq    = getSequence(sid);
  beadID              bid    = readTofBead[sid];
  abColumn           *column = bid.column;
  uint16              link   = bid.link;

  uint32              dl     = 0;
  uint32              bp     = 0;

  while (link != UINT16_MAX) {
    assert(column != NULL);

    char   base = column->_beads[link].base();

    if (base == '-') {
      if (deltas)
        deltas[dl] = bp;
      dl++;
    }

    else
      bp++;

    link   = column->_beads[link].nextOffset();
    column = column->next();
  }

  if (deltas)
    deltas[dl] = 0;
  dl++;

  return(dl);
}



void
abAbacus::getPositions(tgTig *tig) {

  uint32  nd = 0;

  for (uint32 si=0; si<numberOfSequences(); si++)
    nd += getSequenceDeltas(si, NULL);

  resizeArray(tig->_childDeltas, tig->_childDeltasLen, tig->_childDeltasMax, nd, resizeArray_doNothing);

  tig->_childDeltasLen = 0;

  int32  maxPos = 0;

  for (uint32 si=0; si<numberOfSequences(); si++) {
    abSequence *seq   = getSequence(si);
    tgPosition *child = tig->getChild(si);  //  Assumes one-to-one map of seqs to children;

    assert(seq->gkpIdent() == child->ident());

    beadID    fBead = readTofBead[si];
    beadID    lBead = readTolBead[si];

    if (fBead.column == NULL) {
      fprintf(stderr, "WARNING: read %u not in multialignment; position set to 0,0.\n", seq->gkpIdent());
      child->setMinMax(0, 0);
      continue;
    }

    assert(fBead.column != NULL);
    assert(lBead.column != NULL);

    //  Positions are zero-based and inclusive.  The end position gets one added to it to make it true space-based.

    int32 min = fBead.column->position();
    int32 max = lBead.column->position() + 1;

    if (maxPos < min)   maxPos = min;
    if (maxPos < max)   maxPos = max;

    if (seq->isRead() == true) {
      child->setMinMax(min, max);      //  The child remembers its orientation, and sets min/max appropriately.

      //  Grab the deltas for real now, and set the deltaLen to one less than the number returned
      //  (we don't care about the terminating zero!)

      child->_deltaOffset   = tig->_childDeltasLen;
      child->_deltaLen      = getSequenceDeltas(si, tig->_childDeltas + tig->_childDeltasLen);

      tig->_childDeltasLen += child->_deltaLen;

      child->_deltaLen--;

      //fprintf(stderr, "getPositions()-- si %u deltaLen %u offset %u\n",
      //        si.get(), child->_deltaLen, child->_deltaOffset);
    }
  }

  //  Hopefully, maxPos and _gappedLen agree.  If not....?  Note that positions are space-based

  if (maxPos != tig->_gappedLen)
    fprintf(stderr, "WARNING:  maxPos=%d differs from gappedLen=%d.  Huh?\n", maxPos, tig->_gappedLen);

  tig->_layoutLen = tig->_gappedLen;
}





//  If called from unitigConsensus, this needs a rebuild() first.

void
abAbacus::display(FILE *F) {
  int32   pageWidth = 250;

  //  For display, offset the quals to Sanger spec.
  for (uint32 ii=0; ii<numberOfColumns(); ii++)
    _cnsQuals[ii] += '!';

  uint32 numSeqs = numberOfSequences();

  beadID                 *fit       = new beadID [numSeqs];
  uint32                 *fid       = new uint32 [numSeqs];
  char                   *type      = new char   [numSeqs];  //  always 'r' for read
  uint32                 *bgn       = new uint32 [numSeqs];  //  former positions
  uint32                 *end       = new uint32 [numSeqs];

  for (uint32 i=0; i<numSeqs; i++) {
    abSequence  *seq       = getSequence(i);
    abColumn    *bgnColumn = readTofBead[i].column;
    abColumn    *endColumn = readTolBead[i].column;

    if ((bgnColumn == NULL) ||
        (endColumn == NULL)) {
      fid[i]        = 0;
      fit[i].column = NULL;
      fit[i].link   = UINT16_MAX;
      type[i]       = '?';
      bgn[i]        = 0;
      end[i]        = 0;

    } else {
      fid[i]        = seq->gkpIdent();
      fit[i].column = NULL;
      fit[i].link   = UINT16_MAX;
      type[i]       = (seq->isRead() == true) ? 'R' : '?';
      bgn[i]        = bgnColumn->position();
      end[i]        = endColumn->position();
    }
  }

  fprintf(F,"\n==================== abAbacus::display ====================\n");

  for (uint32 window_start=0; window_start < numberOfColumns(); ) {
    fprintf(F,"\n");
    fprintf(F,"%d - %d (length %d)\n", window_start, window_start + pageWidth, numberOfColumns());
    fprintf(F,"%-*.*s <<< consensus\n", pageWidth, pageWidth, _cnsBases + window_start);
    fprintf(F,"%-*.*s <<< quality\n\n", pageWidth, pageWidth, _cnsQuals  + window_start);

    for (uint32 i=0; i<numSeqs; i++) {
      if (fid[i] == 0)
        continue;

      if (end[i] < window_start)
        continue;

      if (window_start + pageWidth < bgn[i])
        continue;


      for (uint32 wi=window_start; wi<window_start + pageWidth; wi++) {

        //  Spaces before the read
        if (wi < bgn[i]) {
          fprintf(F, ".");
        }

        //  Starting in the middle of the sequence.
        else if ((bgn[i] <= wi) &&
                 (wi     <= end[i])) {

          if (fit[i].column == NULL) {
            fit[i] = readTofBead[i];

            while (fit[i].column->position() < wi) {
              fit[i].link   = fit[i].column->_beads[fit[i].link].nextOffset();
              fit[i].column = fit[i].column->next();
            }
          }

          char pc = fit[i].column->_beads[fit[i].link].base();

          if (pc == _cnsBases[wi])
            pc = tolower(pc);
          else
            pc = toupper(pc);

          fprintf(F, "%c", pc);

          fit[i].link   = fit[i].column->_beads[fit[i].link].nextOffset();
          fit[i].column = fit[i].column->next();
        }

        //  Spaces after the read
        else if (end[i] < wi) {
          fprintf(F, ",");
        }

        //  Sequence isn't involved in this window.
        else {
          fprintf(stderr, "bgn=%d end=%d wi=%d\n", bgn[i], end[i], wi);
          assert(0);
        }


        //  If at the end of the window, print the id of the object we just printed.
        //
        if (wi == window_start + pageWidth - 1)
          fprintf(F," <<< %d (%c)\n", fid[i], type[i]);
      }
    }

    window_start += pageWidth;
  }

  //  Unoffset the quals back to integers.

  for (uint32 ii=0; ii<numberOfColumns(); ii++)
    _cnsQuals[ii] -= '!';

  //  Cleanup, goodbye.

  delete [] fit;
  delete [] fid;
  delete [] type;
  delete [] bgn;
  delete [] end;
}
