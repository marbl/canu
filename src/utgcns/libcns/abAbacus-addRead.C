
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

abSeqID
abAbacus::addRead(gkStore *gkpStore,
                  uint32   readID,
                  uint32   askip, uint32 bskip,
                  bool     complemented,
                  map<uint32, gkRead *>     *inPackageRead,
                  map<uint32, gkReadData *> *inPackageReadData) {

  //  Grab the read

  gkRead      *read     = NULL;
  gkReadData  *readData = NULL;

  //  If no package, load the read from the store.  Otherwise, load the read from the package.  This
  //  REQUIRES that the package be in-sync with the unitig.  We fail otherwise.  Hey, it's used for
  //  debugging only...

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

  uint32  seqLen = read->gkRead_sequenceLength() - askip - bskip;

  //  Tell abacus about it.

  increaseArray(_sequences, _sequencesLen, _sequencesMax, 1);

  abSequence   *ns = _sequences + _sequencesLen;

  ns->initialize(readID, _sequencesLen++, seqLen, complemented, _basesLen, _beadsLen);

  //  Make a complement table

  char inv[256] = {0};

  inv['a'] = 't';  inv['A'] = 'T';
  inv['c'] = 'g';  inv['C'] = 'G';
  inv['g'] = 'c';  inv['G'] = 'C';
  inv['t'] = 'a';  inv['T'] = 'A';
  inv['n'] = 'n';  inv['N'] = 'N';
  inv['-'] = '-';

  //  Stash the bases/quals

  {
    char  *seq = readData->gkReadData_getSequence()  + ((complemented == false) ? askip : bskip);
    char  *qlt = readData->gkReadData_getQualities() + ((complemented == false) ? askip : bskip);

    while (_basesMax <= _basesLen + seqLen + 1)
      resizeArrayPair(_bases, _quals, _basesLen, _basesMax, 2 * _basesMax);

    for (uint32 ii=0; ii<seqLen; ii++)
      assert((seq[ii] == 'A') ||
             (seq[ii] == 'C') ||
             (seq[ii] == 'G') ||
             (seq[ii] == 'T') ||
             (seq[ii] == 'N'));

    if (complemented == false)
      for (uint32 ii=0, pp=_basesLen; ii<seqLen; ii++, pp++, _basesLen++) {
        _bases[pp] = seq[ii];
        _quals[pp] = qlt[ii];

        assert(CNS_MIN_QV <= _quals[pp]);
        assert(_quals[pp] <= CNS_MAX_QV);
      }

    else
      for (uint32 ii=seqLen, pp=_basesLen; ii-->0; pp++, _basesLen++) {
        _bases[pp] = inv[ seq[ii] ];
        _quals[pp] =      qlt[ii];

        assert(CNS_MIN_QV <= _quals[pp]);
        assert(_quals[pp] <= CNS_MAX_QV);
      }

    _bases[_basesLen] = 0;  //  NUL terminate the strings so we can use them in aligners
    _quals[_basesLen] = 0;
    _basesLen++;
  }

  delete readData;

  //  Make beads for each base, set the pointer to the first bead

  {
    increaseArray(_beads, _beadsLen, _beadsMax, seqLen);

    uint32  firstBead = _beadsLen;

    for (uint32 bp=0; bp<ns->length(); bp++, _beadsLen++) {
      _beads[_beadsLen].boffset.set(_beadsLen);                   //  Offset into the beads array
      _beads[_beadsLen].soffset.set(ns->firstBase().get() + bp);  //  Offset into the sequence array
      _beads[_beadsLen].foffset = bp;                             //  Offset into the read itself

      //  Check that nothing odd happened with the ident.

      assert(_beads[_beadsLen].ident().get() == ns->firstBead().get() + bp);

      //  Set previous/next bead appropriately.

      if (_beads[_beadsLen].foffset == 0)
        _beads[_beadsLen].prev = abBeadID();
      else
        _beads[_beadsLen].prev.set(_beads[_beadsLen].ident().get() - 1);

      if (_beads[_beadsLen].foffset == ns->length() - 1)
        _beads[_beadsLen].next = abBeadID();
      else
        _beads[_beadsLen].next.set(_beads[_beadsLen].ident().get() + 1);

      _beads[_beadsLen].up           = abBeadID();   //  No up bead yet.
      _beads[_beadsLen].down         = abBeadID();   //  No down bead yet.

      _beads[_beadsLen].frag_index   = ns->ident();  //  Bead is for this read idx.
      _beads[_beadsLen].column_index = abColID();    //  Isn't in a column yet.
    }
  }

  assert(_beads[_beadsLen-1].ident() == ns->lastBead());

  //  Return the (internal) index we saved this read at.

#if 0
  fprintf(stderr, "read %d firstBead %d lastBead %d _basesLen %u\n",
          ns->ident().get(),
          ns->firstBead().get(),
          ns->lastBead().get(),
          _basesLen);
#endif

  return(ns->ident());
}


