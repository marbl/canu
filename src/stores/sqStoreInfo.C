
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "sqStore.H"
#include "files.H"


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

  for (uint32 ii=0; ii<sqRead_largest; ii++) {
    _reads[ii] = 0;
    _bases[ii] = 0;
  }

  checkInfo();
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
    fprintf(stderr, "numBlobs           = " F_U32 "\n", _numBlobs);
    fprintf(stderr, "\n");
  }

  return(failed == 0);
}



bool
sqStoreInfo::examineRead(uint32         ii,
                         sqReadSeq     *seq,
                         sqRead_which   w) {
  sqRead_which  t = w | sqRead_trimmed;

  uint32      length  = (seq->sqReadSeq_valid()   == false) ? 0 : seq->sqReadSeq_length();
  uint32      trimmed = (seq->sqReadSeq_trimmed() == false) ? 0 : seq->sqReadSeq_clearEnd() - seq->sqReadSeq_clearBgn();

  bool        exists  = false;

  if (seq->sqReadSeq_ignoreU() == true)   { length  = 0; }
  if (seq->sqReadSeq_ignoreT() == true)   { trimmed = 0; }

  if (trimmed > length) {
    fprintf(stderr, "sqStoreInfo::examineRead()-- %s read %u length %u (%s) < trimmed %u - %u = %u (%s)\n",
            toString(w),
            ii,
            length,  (seq->sqReadSeq_ignoreU() == true) ? "ignored" : "used",
            seq->sqReadSeq_clearEnd(),
            seq->sqReadSeq_clearBgn(),
            trimmed, (seq->sqReadSeq_ignoreT() == true) ? "ignored" : "used");
  }
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

  for (sqRead_which ii=0; ii<sqRead_largest; ii++) {
    _reads[ii] = 0;
    _bases[ii] = 0;
  }

  if (rawU)
    for (uint32 ii=1; ii<_numReads + 1; ii++)
      examineRead(ii, &rawU[ii], sqRead_raw);

  if (rawC)
    for (uint32 ii=1; ii<_numReads + 1; ii++)
      examineRead(ii, &rawC[ii], sqRead_raw | sqRead_compressed);

  if (corU)
    for (uint32 ii=1; ii<_numReads + 1; ii++)
      examineRead(ii, &corU[ii], sqRead_corrected);

  if (corC)
    for (uint32 ii=1; ii<_numReads + 1; ii++)
      examineRead(ii, &corC[ii], sqRead_corrected | sqRead_compressed);

  //  For convenience, we store the total number of reads
  //  in the sqRead_unset space.  There's no comparable
  //  meaning for _bases[sqRead_unset] though.

  _reads[sqRead_unset] = _numReads;
  _bases[sqRead_unset] = 0;
}



bool
sqStoreInfo::readInfo8(char *storeName) {
  uint64  magic;
  uint64  version;

  FILE  *I = AS_UTL_openInputFile(storeName, '/', "info");

  loadFromFile(&_sqMagic,            "sqInfo_sqMagic", I);
  loadFromFile(&_sqVersion,          "sqInfo_sqVersion", I);

  if (_sqMagic != SQ_MAGIC8) {   //  The magic number changed at version 9, to fix a type,
    AS_UTL_closeFile(I);         //  and to get rid of GKP.
    return(false);
  }

  if (_sqVersion < 8) {
    fprintf(stderr, "sqStore_loadInfo()-- Error: Unsupported historical version %lu detected.\n", _sqVersion);
    AS_UTL_closeFile(I);
    exit(1);
  }

  if (_sqVersion == 8) {
    fprintf(stderr, "sqStore_loadInfo()-- Warning: obsolete version 8 detected; adjusting.\n");

    loadFromFile(&_sqLibrarySize,      "sqInfo_sqLibrarySize", I);
    loadFromFile(&_sqReadSize,         "sqInfo_sqReadSize", I);
    loadFromFile(&_sqMaxLibrariesBits, "sqInfo_sqMaxLibrariesBits", I);
    loadFromFile(&_sqLibraryNameSize,  "sqInfo_sqLibraryNameSize", I);
    loadFromFile(&_sqMaxReadBits,      "sqInfo_sqMaxReadBits", I);
    loadFromFile(&_sqMaxReadLenBits,   "sqInfo_sqMaxReadLenBits", I);

    loadFromFile(&_numLibraries,     "sqInfo_numLibraries", I);
    loadFromFile(&_numReads,         "sqInfo_numReads", I);

    loadFromFile(_reads, "sqInfo_reads", sqRead_largest, I);
    loadFromFile(_bases, "sqInfo_bases", sqRead_largest, I);

    //  Figure out how many blobs files there are in the store.

    {
      char    bName[FILENAME_MAX + 1] = {0};
      uint32  bNum = 0;

      makeBlobName(storeName, bNum, bName);

      while (fileExists(bName) == true)
        makeBlobName(storeName, ++bNum, bName);

      _numBlobs  = bNum;

      fprintf(stderr, "sqStore_loadInfo()-- Warning: found numBlobs %u\n", _numBlobs);
    }

    //  "Upgrade" this to version 9 data.

    _sqMagic   = SQ_MAGIC;
    _sqVersion = 9;

    AS_UTL_closeFile(I);
    return(true);
  }

  if (_sqVersion > 8) {
    AS_UTL_closeFile(I);
    return(false);
  }

  assert(0);
  return(false);
}



void
sqStoreInfo::readInfo(char *storePath) {

  //  If no file, don't read it.

  if (fileExists(storePath, '/', "info") == false)
    return;

  //  Try to read a v8 info file.  If that fails, read the current version.

  if (readInfo8(storePath) == false) {
    readBuffer  *B = new readBuffer(storePath, '/', "info");

    B->readIFFobject("MAGC", _sqMagic);
    B->readIFFobject("VERS", _sqVersion);

    B->readIFFobject("LSIZ", _sqLibrarySize);
    B->readIFFobject("RSIZ", _sqReadSize);
    B->readIFFobject("MLB ", _sqMaxLibrariesBits);
    B->readIFFobject("LNS ", _sqLibraryNameSize);
    B->readIFFobject("MRB ", _sqMaxReadBits);
    B->readIFFobject("MRLB", _sqMaxReadLenBits);

    B->readIFFobject("NLIB", _numLibraries);
    B->readIFFobject("NREA", _numReads);
    B->readIFFobject("NBLO", _numBlobs);

    B->readIFFarray("READ", _reads, sqRead_largest);
    B->readIFFarray("BASE", _bases, sqRead_largest);

    delete B;
  }

  //  Check that the info is compatible with our code.

  if (checkInfo() == false) {
    fprintf(stderr, "\n");
    fprintf(stderr, "ERROR:  Can't open store '%s': parameters in sqStore.H and sqRead.H are incompatible with the store.\n", storePath);
    exit(1);
  }
}


void
sqStoreInfo::writeInfo(char *storePath) {
  writeBuffer   *B = new writeBuffer(storePath, '/', "info", "w", 16384);

  B->writeIFFobject("MAGC", _sqMagic);
  B->writeIFFobject("VERS", _sqVersion);

  B->writeIFFobject("LSIZ", _sqLibrarySize);
  B->writeIFFobject("RSIZ", _sqReadSize);
  B->writeIFFobject("MLB ", _sqMaxLibrariesBits);
  B->writeIFFobject("LNS ", _sqLibraryNameSize);
  B->writeIFFobject("MRB ", _sqMaxReadBits);
  B->writeIFFobject("MRLB", _sqMaxReadLenBits);

  B->writeIFFobject("NLIB", _numLibraries);
  B->writeIFFobject("NREA", _numReads);
  B->writeIFFobject("NBLO", _numBlobs);

  B->writeIFFarray("READ", _reads, sqRead_largest);
  B->writeIFFarray("BASE", _bases, sqRead_largest);

  delete B;
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
