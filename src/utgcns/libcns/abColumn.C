
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
 *    Brian P. Walenz beginning on 2014-NOV-17
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "abAbacus.H"


void
abColumn::CheckBaseCounts(abAbacus *abacus) {
  uint32 counts[256] = {0};
  uint32 nBeads = 0;

  //  Why skip the last column??
  //if (next == -1)
  //  return;

  abBead *bead = abacus->getBead(callID());

  while (bead->downID().isValid()) {
    bead = abacus->getBead(bead->downID());

    counts[ abacus->getBase(bead->ident()) ]++;

    nBeads++;
  }

  if (counts['A'] != GetColumnBaseCount('A'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u A computed %u != stored %u (%u beads)\n", ident().get(), counts['A'], GetColumnBaseCount('A'), nBeads);

  if (counts['C'] != GetColumnBaseCount('C'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u C computed %u != stored %u (%u beads)\n", ident().get(), counts['C'], GetColumnBaseCount('C'), nBeads);

  if (counts['G'] != GetColumnBaseCount('G'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u G computed %u != stored %u (%u beads)\n", ident().get(), counts['G'], GetColumnBaseCount('G'), nBeads);

  if (counts['T'] != GetColumnBaseCount('T'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u T computed %u != stored %u (%u beads)\n", ident().get(), counts['T'], GetColumnBaseCount('T'), nBeads);

  if (counts['-'] != GetColumnBaseCount('-'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u - computed %u != stored %u (%u beads)\n", ident().get(), counts['-'], GetColumnBaseCount('-'), nBeads);
};
