
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

  for (uint32 ii=0; ii<sqRead_allset; ii++) {
    _reads[ii] = 0;
    _bases[ii] = 0;
  }
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

  if (failed) {
    fprintf(stderr, "\n");
    fprintf(stderr, "sqMagic            = 0x" F_X64 "\n", _sqMagic);
    fprintf(stderr, "sqVersion          = 0x" F_X64 "\n", _sqVersion);
    fprintf(stderr, "\n");
    fprintf(stderr, "sqLibrarySize      = " F_U32 "\n", _sqLibrarySize);
    fprintf(stderr, "sqReadSize         = " F_U32 "\n", _sqReadSize);
    fprintf(stderr, "sqMaxLibrariesBits = " F_U32 "\n", _sqMaxLibrariesBits);
    fprintf(stderr, "sqLibraryNameSize  = " F_U32 "\n", _sqLibraryNameSize);
    fprintf(stderr, "sqMaxReadBits      = " F_U32 "\n", _sqMaxReadBits);
    fprintf(stderr, "sqMaxReadLenBits   = " F_U32 "\n", _sqMaxReadLenBits);
    fprintf(stderr, "\n");
    fprintf(stderr, "numLibraries       = " F_U32 "\n", _numLibraries);
    fprintf(stderr, "numReads           = " F_U32 "\n", _numReads);
    fprintf(stderr, "\n");
  }

  return(failed == 0);
}



bool
sqStoreInfo::examineRead(sqReadSeq     *seq,
                         sqRead_which   w) {
  sqRead_which  t = w | sqRead_trimmed;

  uint32      length  = (seq->sqReadSeq_valid()   == false) ? 0 : seq->sqReadSeq_length();
  uint32      trimmed = (seq->sqReadSeq_trimmed() == false) ? 0 : seq->sqReadSeq_clearEnd() - seq->sqReadSeq_clearBgn();

  bool        exists  = false;

  if (seq->sqReadSeq_ignoreU() == true)   { length  = 0; }
  if (seq->sqReadSeq_ignoreT() == true)   { trimmed = 0; }

  assert(trimmed <= length);

  if (seq->sqReadSeq_ignoreU())
    assert(seq->sqReadSeq_ignoreT());  //  If untrimmed ignored, trimmed must be ignored too.

  if (length > 0) {
    exists     = true;
    _reads[w] += 1;
    _bases[w] += length;
  }

  if (trimmed > 0) {
    exists     = true;
    _reads[t] += 1;
    _bases[t] += trimmed;
  }

  return(exists);
}


void
sqStoreInfo::update(sqReadSeq  *rawU,
                    sqReadSeq  *rawC,
                    sqReadSeq  *corU,
                    sqReadSeq  *corC) {

  for (sqRead_which ii=0; ii<sqRead_allset; ii++) {
    _reads[ii] = 0;
    _bases[ii] = 0;
  }

  if (rawU)
    for (uint32 ii=0; ii<_numReads + 1; ii++)
      examineRead(&rawU[ii], sqRead_raw);

  if (rawC)
    for (uint32 ii=0; ii<_numReads + 1; ii++)
      examineRead(&rawC[ii], sqRead_raw | sqRead_compressed);

  if (corU)
    for (uint32 ii=0; ii<_numReads + 1; ii++)
      examineRead(&corU[ii], sqRead_corrected);

  if (corC)
    for (uint32 ii=0; ii<_numReads + 1; ii++)
      examineRead(&corC[ii], sqRead_corrected | sqRead_compressed);

  //  For convenience, we store the total number of reads
  //  in the sqRead_unset space.  There's no comparable
  //  meaning for _bases[sqRead_unset] though.

  _reads[sqRead_unset] = _numReads;
  _bases[sqRead_unset] = 0;
}



void
sqStoreInfo::writeInfoAsText(FILE *F) {
  sqRead_which  raw = sqRead_raw;
  sqRead_which  cor = sqRead_corrected;
  sqRead_which  cmp = sqRead_compressed;
  sqRead_which  tri = sqRead_trimmed;
  sqRead_which  x   = sqRead_unset;

  fprintf(F, "     Reads        Bases Read Type\n");
  fprintf(F, "---------- ------------ ----------------------------------------\n");

  ;                      fprintf(F, "%10" F_U32P "            - total-reads\n", _numReads);

  x = raw;               fprintf(F, "%10" F_U64P " %12" F_U64P " %s\n", _reads[x], _bases[x], toString(x));
  x = raw | tri;         fprintf(F, "%10" F_U64P " %12" F_U64P " %s\n", _reads[x], _bases[x], toString(x));
  x = raw | cmp;         fprintf(F, "%10" F_U64P " %12" F_U64P " %s\n", _reads[x], _bases[x], toString(x));
  x = raw | cmp | tri;   fprintf(F, "%10" F_U64P " %12" F_U64P " %s\n", _reads[x], _bases[x], toString(x));

  x = cor;               fprintf(F, "%10" F_U64P " %12" F_U64P " %s\n", _reads[x], _bases[x], toString(x));
  x = cor | tri;         fprintf(F, "%10" F_U64P " %12" F_U64P " %s\n", _reads[x], _bases[x], toString(x));
  x = cor | cmp;         fprintf(F, "%10" F_U64P " %12" F_U64P " %s\n", _reads[x], _bases[x], toString(x));
  x = cor | cmp | tri;   fprintf(F, "%10" F_U64P " %12" F_U64P " %s\n", _reads[x], _bases[x], toString(x));
}
