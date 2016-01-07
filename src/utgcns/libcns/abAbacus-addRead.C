
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
 *    src/utgcns/libcns/abacus-addRead.C
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
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "abAbacus.H"

//  CA8 code for adding a unitig (and expanding it to include all the reads) exists
//  last in b8cc87300a0b5da87513ea1a6c02e8280af30cd0.


abSequence::abSequence(uint32  readID,
                       uint32  length,
                       char   *seq,
                       char   *qlt,
                       uint32  complemented) {
  _iid              = readID;

  _is_read          = true;

  _length           = length;

  _complement       = complemented;

  _bases            = new char  [_length + 1];
  _quals            = new uint8 [_length + 1];

  //  Make a complement table

  char inv[256] = {0};

  inv['a'] = 't';  inv['A'] = 'T';
  inv['c'] = 'g';  inv['C'] = 'G';
  inv['g'] = 'c';  inv['G'] = 'C';
  inv['t'] = 'a';  inv['T'] = 'A';
  inv['n'] = 'n';  inv['N'] = 'N';
  inv['-'] = '-';

  //  Stash the bases/quals

  for (uint32 ii=0; ii<_length; ii++)
    assert((seq[ii] == 'A') ||
           (seq[ii] == 'C') ||
           (seq[ii] == 'G') ||
           (seq[ii] == 'T') ||
           (seq[ii] == 'N'));

  if (complemented == false)
    for (uint32 ii=0, pp=0; ii<_length; ii++, pp++) {
      _bases[pp] = seq[ii];
      _quals[pp] = qlt[ii];
    }

  else
    for (uint32 ii=_length, pp=0; ii-->0; pp++) {
      _bases[pp] = inv[ seq[ii] ];
      _quals[pp] =      qlt[ii];
    }

  _bases[_length] = 0;  //  NUL terminate the strings so we can use them in aligners.
  _quals[_length] = 0;  //  Not actually a string, the 0 is just another QV=0 entry.
};




void
abAbacus::addRead(gkStore *gkpStore,
                  uint32   readID,
                  uint32   askip, uint32 bskip,
                  bool     complemented,
                  map<uint32, gkRead *>     *inPackageRead,
                  map<uint32, gkReadData *> *inPackageReadData) {

  //  Grab the read.  If there is no package, load the read from the store.  Otherwise, load the
  //  read from the package.  This REQUIRES that the package be in-sync with the unitig.  We fail
  //  otherwise.  Hey, it's used for debugging only...

  gkRead      *read     = NULL;
  gkReadData  *readData = NULL;

  if (inPackageRead == NULL) {
    read     = gkpStore->gkStore_getRead(readID);
    readData = new gkReadData;

    gkpStore->gkStore_loadReadData(read, readData);
  }

  else {
    read     = (*inPackageRead)[readID];
    readData = (*inPackageReadData)[readID];
  }

  assert(read     != NULL);
  assert(readData != NULL);

  //  Grab seq/qlt from the read, offset to the proper begin and length.

  uint32  seqLen = read->gkRead_sequenceLength() - askip - bskip;
  char   *seq    = readData->gkReadData_getSequence()  + ((complemented == false) ? askip : bskip);
  char   *qlt    = readData->gkReadData_getQualities() + ((complemented == false) ? askip : bskip);

  //  Tell abacus about it.  We could pre-allocate _sequences (in the constructor) but this is
  //  relatively painless and makes life easier outside here.

  increaseArray(_sequences, _sequencesLen, _sequencesMax, 1);

  _sequences[_sequencesLen++] = new abSequence(readID, seqLen, seq, qlt, complemented);

  delete readData;
}


