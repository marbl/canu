
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
 *    src/stores/gkStore.C
 *    src/stores/gkStoreInfo.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-26 to 2015-AUG-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sqStore.H"


sqStoreInfo::sqStoreInfo() {
  _sqMagic            = SQ_MAGIC;
  _sqVersion          = SQ_VERSION;

  _sqLibrarySize      = sizeof(sqLibrary);
  _sqReadSize         = sizeof(sqRead);
  _sqMaxLibrariesBits = AS_MAX_LIBRARIES_BITS;
  _sqLibraryNameSize  = LIBRARY_NAME_SIZE;
  _sqMaxReadBits      = AS_MAX_READS_BITS;
  _sqMaxReadLenBits   = AS_MAX_READLEN_BITS;

  _numLibraries       = 0;
  _numReads           = 0;
  _numBlobs           = 0;

  _numRawReads        = 0;
  _numCorrectedReads  = 0;
  _numTrimmedReads    = 0;

  _numRawBases        = 0;
  _numCorrectedBases  = 0;
  _numTrimmedBases    = 0;
}



sqStoreInfo::~sqStoreInfo() {
}



bool
sqStoreInfo::checkInfo(void) {
  uint32  failed = 0;

  if (_sqMagic            != SQ_MAGIC)
    failed += fprintf(stderr, "ERROR:  sqMagic in store = 0x%016" F_X64P ", differs from executable = 0x%016" F_X64P "\n",
                      _sqMagic, SQ_MAGIC);

  if (_sqVersion          != SQ_VERSION)
    failed += fprintf(stderr, "ERROR:  sqVersion in store = 0x%016" F_X64P ", differs from executable = 0x%016" F_X64P "\n",
                      _sqVersion, SQ_VERSION);

  if (_sqLibrarySize      != sizeof(sqLibrary))
    failed += fprintf(stderr, "ERROR:  sqLibrary size in store = " F_U32 ", differs from executable = " F_SIZE_T "\n",
                      _sqLibrarySize, sizeof(sqLibrary));

  if (_sqReadSize         != sizeof(sqRead))
    failed += fprintf(stderr, "ERROR:  sqRead size in store = " F_U32 ", differs from executable = " F_SIZE_T "\n",
                      _sqReadSize, sizeof(sqRead));

  if (_sqMaxLibrariesBits != AS_MAX_LIBRARIES_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_LIBRARIES_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _sqMaxLibrariesBits, AS_MAX_LIBRARIES_BITS);

  if (_sqLibraryNameSize  != LIBRARY_NAME_SIZE)
    failed += fprintf(stderr, "ERROR:  LIBRARY_NAME_SIZE in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _sqLibraryNameSize, LIBRARY_NAME_SIZE);

  if (_sqMaxReadBits      != AS_MAX_READS_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_READS_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _sqMaxReadBits, AS_MAX_READS_BITS);

  if (_sqMaxReadLenBits   != AS_MAX_READLEN_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_READLEN_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _sqMaxReadLenBits, AS_MAX_READLEN_BITS);

  return(failed == 0);
}



void
sqStoreInfo::recountReads(sqRead *reads) {

  _numRawReads = _numCorrectedReads = _numTrimmedReads = 0;
  _numRawBases = _numCorrectedBases = _numTrimmedBases = 0;

  for (uint32 ii=0; ii<_numReads + 1; ii++) {
    uint32 rr = reads[ii].sqRead_sequenceLength(sqRead_raw);
    uint32 rc = reads[ii].sqRead_sequenceLength(sqRead_corrected);
    uint32 rt = reads[ii].sqRead_sequenceLength(sqRead_trimmed);

    if (rr > 0) {
      _numRawReads++;
      _numRawBases += rr;
    }

    if (rc > 0) {
      _numCorrectedReads++;
      _numCorrectedBases += rc;
    }

    if (rt > 0) {
      _numTrimmedReads++;
      _numTrimmedBases += rt;
    }
  }
}



void
sqStoreInfo::setLastBlob(sqStoreBlobWriter *writer) {
  if (writer)
    _numBlobs = writer->writtenBlob();
}



void
sqStoreInfo::writeInfoAsText(FILE *F) {
  fprintf(F, "sqMagic            = 0x" F_X64 "\n", _sqMagic);
  fprintf(F, "sqVersion          = 0x" F_X64 "\n", _sqVersion);
  fprintf(F, "\n");
  fprintf(F, "sqLibrarySize      = " F_U32 "\n", _sqLibrarySize);
  fprintf(F, "sqReadSize         = " F_U32 "\n", _sqReadSize);
  fprintf(F, "sqMaxLibrariesBits = " F_U32 "\n", _sqMaxLibrariesBits);
  fprintf(F, "sqLibraryNameSize  = " F_U32 "\n", _sqLibraryNameSize);
  fprintf(F, "sqMaxReadBits      = " F_U32 "\n", _sqMaxReadBits);
  fprintf(F, "sqMaxReadLenBits   = " F_U32 "\n", _sqMaxReadLenBits);
  fprintf(F, "\n");
  fprintf(F, "numLibraries       = " F_U32 "\n", _numLibraries);
  fprintf(F, "numReads           = " F_U32 "\n", _numReads);
  fprintf(F, "numBlobs           = " F_U32 "\n", _numBlobs);
  fprintf(F, "\n");
  fprintf(F, "numRawReads        = " F_U32 "\n", _numRawReads);
  fprintf(F, "numCorrectedReads  = " F_U32 "\n", _numCorrectedReads);
  fprintf(F, "numTrimmedReads    = " F_U32 "\n", _numTrimmedReads);
  fprintf(F, "\n");
  fprintf(F, "numRawBases        = " F_U64 "\n", _numRawBases);
  fprintf(F, "numCorrectedBases  = " F_U64 "\n", _numCorrectedBases);
  fprintf(F, "numTrimmedBases    = " F_U64 "\n", _numTrimmedBases);
}
