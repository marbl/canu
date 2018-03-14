
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
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-MAR-13
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "gkStore.H"


gkStoreInfo::gkStoreInfo() {
  _gkMagic            = GK_MAGIC;
  _gkVersion          = GK_VERSION;

  _gkLibrarySize      = sizeof(gkLibrary);
  _gkReadSize         = sizeof(gkRead);
  _gkMaxLibrariesBits = AS_MAX_LIBRARIES_BITS;
  _gkLibraryNameSize  = LIBRARY_NAME_SIZE;
  _gkMaxReadBits      = AS_MAX_READS_BITS;
  _gkMaxReadLenBits   = AS_MAX_READLEN_BITS;

  _numLibraries       = 0;
  _numReads           = 0;

  _numRawReads        = 0;
  _numCorrectedReads  = 0;
  _numTrimmedReads    = 0;

  _numRawBases        = 0;
  _numCorrectedBases  = 0;
  _numTrimmedBases    = 0;
}



gkStoreInfo::~gkStoreInfo() {
}



bool
gkStoreInfo::checkInfo(void) {
  uint32  failed = 0;

  if (_gkMagic            != GK_MAGIC)
    failed += fprintf(stderr, "ERROR:  gkMagic in store = 0x%016" F_X64P ", differs from executable = 0x%016" F_X64P "\n",
                      _gkMagic, GK_MAGIC);

  if (_gkVersion          != GK_VERSION)
    failed += fprintf(stderr, "ERROR:  gkVersion in store = 0x%016" F_X64P ", differs from executable = 0x%016" F_X64P "\n",
                      _gkVersion, GK_VERSION);

  if (_gkLibrarySize      != sizeof(gkLibrary))
    failed += fprintf(stderr, "ERROR:  gkLibrary size in store = " F_U32 ", differs from executable = " F_SIZE_T "\n",
                      _gkLibrarySize, sizeof(gkLibrary));

  if (_gkReadSize         != sizeof(gkRead))
    failed += fprintf(stderr, "ERROR:  gkRead size in store = " F_U32 ", differs from executable = " F_SIZE_T "\n",
                      _gkReadSize, sizeof(gkRead));

  if (_gkMaxLibrariesBits != AS_MAX_LIBRARIES_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_LIBRARIES_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _gkMaxLibrariesBits, AS_MAX_LIBRARIES_BITS);

  if (_gkLibraryNameSize  != LIBRARY_NAME_SIZE)
    failed += fprintf(stderr, "ERROR:  LIBRARY_NAME_SIZE in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _gkLibraryNameSize, LIBRARY_NAME_SIZE);

  if (_gkMaxReadBits      != AS_MAX_READS_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_READS_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _gkMaxReadBits, AS_MAX_READS_BITS);

  if (_gkMaxReadLenBits   != AS_MAX_READLEN_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_READLEN_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _gkMaxReadLenBits, AS_MAX_READLEN_BITS);

  return(failed == 0);
}



void
gkStoreInfo::recountReads(gkRead *reads) {

  _numRawReads = _numCorrectedReads = _numTrimmedReads = 0;
  _numRawBases = _numCorrectedBases = _numTrimmedBases = 0;

  for (uint32 ii=0; ii<_numReads + 1; ii++) {
    if (reads[ii].gkRead_rawLength() > 0) {
      _numRawReads++;
      _numRawBases += reads[ii].gkRead_rawLength();
    }

    if (reads[ii].gkRead_correctedLength() > 0) {
      _numCorrectedReads++;
      _numCorrectedBases += reads[ii].gkRead_correctedLength();
    }

    if (reads[ii].gkRead_clearBgn() < reads[ii].gkRead_clearEnd()) {
      _numTrimmedReads++;
      _numTrimmedBases += reads[ii].gkRead_clearEnd() - reads[ii].gkRead_clearBgn();
    }
  }
}



void
gkStoreInfo::writeInfoAsText(FILE *F) {
  fprintf(F, "gkMagic            = 0x" F_X64 "\n", _gkMagic);
  fprintf(F, "gkVersion          = 0x" F_X64 "\n", _gkVersion);
  fprintf(F, "\n");
  fprintf(F, "gkLibrarySize      = " F_U32 "\n", _gkLibrarySize);
  fprintf(F, "gkReadSize         = " F_U32 "\n", _gkReadSize);
  fprintf(F, "gkMaxLibrariesBits = " F_U32 "\n", _gkMaxLibrariesBits);
  fprintf(F, "gkLibraryNameSize  = " F_U32 "\n", _gkLibraryNameSize);
  fprintf(F, "gkMaxReadBits      = " F_U32 "\n", _gkMaxReadBits);
  fprintf(F, "gkMaxReadLenBits   = " F_U32 "\n", _gkMaxReadLenBits);
  fprintf(F, "\n");
  fprintf(F, "numLibraries       = " F_U32 "\n", _numLibraries);
  fprintf(F, "numReads           = " F_U32 "\n", _numReads);
  fprintf(F, "\n");
  fprintf(F, "numRawReads        = " F_U32 "\n", _numRawReads);
  fprintf(F, "numCorrectedReads  = " F_U32 "\n", _numCorrectedReads);
  fprintf(F, "numTrimmedReads    = " F_U32 "\n", _numTrimmedReads);
  fprintf(F, "\n");
  fprintf(F, "numRawBases        = " F_U64 "\n", _numRawBases);
  fprintf(F, "numCorrectedBases  = " F_U64 "\n", _numCorrectedBases);
  fprintf(F, "numTrimmedBases    = " F_U64 "\n", _numTrimmedBases);
}
