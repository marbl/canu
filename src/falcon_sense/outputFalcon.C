
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-APR-20 to 2015-MAY-20
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "outputFalcon.H"

#include "AS_UTL_reverseComplement.H"


//  The falcon consensus format:
//
//  name sequence
//  read sequence
//  read sequence
//  read sequence
//  read sequence
//  + +            #  generate consensus for the 'name' sequence using 'read' sequences
//  ...
//  - -            #  To end processing
//


void
outputFalcon(gkStore      *gkpStore,
             tgTig        *tig,
             bool          trimToAlign,
             FILE         *F,
             gkReadData   *readData) {

  gkpStore->gkStore_loadReadData(tig->tigID(), readData);

  fprintf(F, "read"F_U32" %s\n", tig->tigID(), readData->gkReadData_getSequence());

  for (uint32 cc=0; cc<tig->numberOfChildren(); cc++) {
    tgPosition  *child = tig->getChild(cc);

    gkpStore->gkStore_loadReadData(child->ident(), readData);

    if (child->isReverse())
      reverseComplementSequence(readData->gkReadData_getSequence(),
                                readData->gkReadData_getRead()->gkRead_sequenceLength());

    //  For debugging/testing, skip one orientation of overlap.
    //
    //if (child->isReverse() == false)
    //  continue;
    //if (child->isReverse() == true)
    //  continue;

    //  Trim the read to the aligned bit
    char   *seq = readData->gkReadData_getSequence();

    if (trimToAlign) {
      seq += child->_askip;
      seq[ readData->gkReadData_getRead()->gkRead_sequenceLength() - child->_askip - child->_bskip ] = 0;
    }

    fprintf(F, "data"F_U32" %s\n", tig->getChild(cc)->ident(), seq);
  }

  fprintf(F, "+ +\n");
}

