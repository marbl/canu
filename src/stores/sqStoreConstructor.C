
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
#include "ovOverlap.H"

#include "files.H"



void
sqStore::sqStore_loadMetadata(void) {
  char    name[FILENAME_MAX+1];

  _librariesAlloc = _info.sqInfo_lastLibraryID() + 1;
  _readsAlloc     = _info.sqInfo_lastReadID()    + 1;

  _libraries  = new sqLibrary  [_librariesAlloc];
  _meta       = new sqReadMeta [_readsAlloc];
  _rawU       = NULL;
  _rawC       = NULL;
  _corU       = NULL;
  _corC       = NULL;

  AS_UTL_loadFile(_storePath, '/', "libraries", _libraries, _librariesAlloc);
  AS_UTL_loadFile(_storePath, '/', "reads",     _meta,      _readsAlloc);

  //  If the user hasn't set a default version (by calling
  //  sqRead_setDefaultVersion() before a sqStore object is constructed),
  //  pick the plausible most recent version.

  if (sqRead_defaultVersion == sqRead_unset) {
    if (sqStore_getNumReads(sqRead_raw)                        > 0)   sqRead_defaultVersion = sqRead_raw;
    if (sqStore_getNumReads(sqRead_raw | sqRead_trimmed)       > 0)   sqRead_defaultVersion = sqRead_raw       | sqRead_trimmed;
    if (sqStore_getNumReads(sqRead_corrected)                  > 0)   sqRead_defaultVersion = sqRead_corrected;
    if (sqStore_getNumReads(sqRead_corrected | sqRead_trimmed) > 0)   sqRead_defaultVersion = sqRead_corrected | sqRead_trimmed;
  }

  //  Otherwise, if neither raw or corrected reads were specified, decide
  //  which one is the most recent.  This is for when the user only asks
  //  for 'trimmed' reads, not caring if they are raw or corrected.

  else if (((sqRead_defaultVersion & sqRead_raw)       == sqRead_unset) &&
           ((sqRead_defaultVersion & sqRead_corrected) == sqRead_unset)) {
    if      (sqStore_getNumReads(sqRead_corrected) > 0)   sqRead_defaultVersion |= sqRead_corrected;
    else if (sqStore_getNumReads(sqRead_raw)       > 0)   sqRead_defaultVersion |= sqRead_raw;
  }

  //  The store itself can now insist that 'compressed' reads be used by default.

  if (((sqRead_defaultVersion & sqRead_normal) == sqRead_unset) &&
      (fileExists(sqStore_path(), '/', "homopolymerCompression") == true)) {
    sqRead_defaultVersion &= ~sqRead_normal;
    sqRead_defaultVersion |=  sqRead_compressed;
  }

  //  Default version MUST be set now.

  assert(sqRead_defaultVersion != sqRead_unset);

  //  We can maybe, eventually, be clever and load these on-demand.
  //  For now, load all the metadata.

  AS_UTL_loadFile(_storePath, '/', "reads-rawu", _rawU = new sqReadSeq[_readsAlloc], _readsAlloc);
  AS_UTL_loadFile(_storePath, '/', "reads-rawc", _rawC = new sqReadSeq[_readsAlloc], _readsAlloc);
  AS_UTL_loadFile(_storePath, '/', "reads-coru", _corU = new sqReadSeq[_readsAlloc], _readsAlloc);
  AS_UTL_loadFile(_storePath, '/', "reads-corc", _corC = new sqReadSeq[_readsAlloc], _readsAlloc);
}






sqStore::sqStore(char const    *storePath_,
                 sqStore_mode   mode_) {
  char    nameL[FILENAME_MAX+1];
  char    nameR[FILENAME_MAX+1];
  char    nameB[FILENAME_MAX+1];

  //  Set the sqStore pointer in an overlap.

  ovOverlap::sqStoreAttach(this);

  //  Clear ourself, to make valgrind happier.

  memset(_storePath, 0, sizeof(char) * (FILENAME_MAX + 1));

  _mode                   = mode_;

  _librariesAlloc         = 0;
  _libraries              = NULL;

  _readsAlloc             = 0;
  _meta                   = NULL;
  _rawU                   = NULL;
  _rawC                   = NULL;
  _corU                   = NULL;
  _corC                   = NULL;

  _blobReader             = NULL;
  _blobWriter             = NULL;

  //  Save the path and name.

  if (storePath_)   strncpy(_storePath, storePath_, FILENAME_MAX);   //  storePath must always exist though.

  //  Load the info file, if it exists.  And if not, do nothing.

  _info.readInfo(_storePath);

  //
  //  CREATE - allocate some memory for saving libraries and reads, and create a file to dump the data into.
  //

  if (_mode == sqStore_create) {
    if (directoryExists(_storePath) == true)
      fprintf(stderr, "ERROR:  Can't create store '%s': store already exists.\n", _storePath), exit(1);

    AS_UTL_mkdir(_storePath);

    _librariesAlloc = 32;           //  _libraries and
    _readsAlloc     = 32768;        //  _reads MUST be preallocated.

    _libraries      = new sqLibrary  [_librariesAlloc];
    _meta           = new sqReadMeta [_readsAlloc];
    _rawU           = new sqReadSeq  [_readsAlloc];
    _rawC           = new sqReadSeq  [_readsAlloc];
    _corU           = new sqReadSeq  [_readsAlloc];
    _corC           = new sqReadSeq  [_readsAlloc];

    for (uint32 ii=0; ii<_librariesAlloc; ii++)
      _libraries[ii].sqLibrary_initialize();

    for (uint32 ii=0; ii<_librariesAlloc; ii++) {
      _meta[ii].sqReadMeta_initialize();
      _rawU[ii].sqReadSeq_initialize();
      _rawC[ii].sqReadSeq_initialize();
      _corU[ii].sqReadSeq_initialize();
      _corC[ii].sqReadSeq_initialize();
    }
    
    _blobWriter     = new sqStoreBlobWriter(_storePath, &_info);

    return;
  }

  //
  //  Not creating, so the store MUST exist.
  //  Load metadata and check it is compatible.
  //  Make a writer if we're extending, then make a reader.
  //

  if (directoryExists(_storePath) == false)
    fprintf(stderr, "sqStore()--  failed to open '%s' for read-only access: store doesn't exist.\n", _storePath), exit(1);

  sqStore_loadMetadata();

  if (_mode == sqStore_extend)
    _blobWriter = new sqStoreBlobWriter(_storePath, &_info);

  _blobReader = new sqStoreBlobReader(_storePath);
}






sqStore::~sqStore() {
  char    No[FILENAME_MAX+1];
  char    Nn[FILENAME_MAX+1];
  uint32  V = 1;

  //  Save original metadata.

  if (_mode == sqStore_extend) {
    snprintf(No, FILENAME_MAX, "%s/version.%03" F_U32P, _storePath, V);
    while (directoryExists(No) == true) {
      V++;
      snprintf(No, FILENAME_MAX, "%s/version.%03" F_U32P, _storePath, V);
    }

    AS_UTL_mkdir(No);

    snprintf(No, FILENAME_MAX, "%s/libraries", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/libraries", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads-rawu", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads-rawu", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads-rawc", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads-rawc", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads-coru", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads-coru", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads-corc", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads-corc", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/info", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/info", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/info.txt", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/info.txt", _storePath, V);
    AS_UTL_rename(No, Nn);
  }

  //  Recount.  This is significantly easier than trying to adjust stats
  //  as reads are added - especially for trimming.

  if ((_mode == sqStore_create) ||
      (_mode == sqStore_extend))
    _info.update(_rawU, _rawC, _corU, _corC);

  //  Write updated metadata.

  if ((_mode == sqStore_create) ||
      (_mode == sqStore_extend)) {
    AS_UTL_saveFile(_storePath, '/', "libraries",   _libraries, sqStore_lastLibraryID() + 1);
    AS_UTL_saveFile(_storePath, '/', "reads",       _meta,      sqStore_lastReadID()    + 1);
    AS_UTL_saveFile(_storePath, '/', "reads-rawu",  _rawU,      sqStore_lastReadID()    + 1);
    AS_UTL_saveFile(_storePath, '/', "reads-rawc",  _rawC,      sqStore_lastReadID()    + 1);
    AS_UTL_saveFile(_storePath, '/', "reads-coru",  _corU,      sqStore_lastReadID()    + 1);
    AS_UTL_saveFile(_storePath, '/', "reads-corc",  _corC,      sqStore_lastReadID()    + 1);

    _info.writeInfo(_storePath);

    FILE *F = AS_UTL_openOutputFile(_storePath, '/', "info.txt");   //  Used by Canu/Gatekeeper.pm
    _info.writeInfoAsText(F);                                       //  Do not remove!
    AS_UTL_closeFile(F, _storePath, '/', "info.txt");
  }

  //  Clean up.

  delete [] _libraries;
  delete [] _meta;
  delete [] _rawU;
  delete [] _rawC;
  delete [] _corU;
  delete [] _corC;

  delete    _blobWriter;
  delete    _blobReader;
};



void
sqStore::sqStore_delete(void) {
  char path[FILENAME_MAX+1];

  snprintf(path, FILENAME_MAX, "%s/info",       _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/libraries",  _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/reads",      _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/reads-rawu", _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/reads-rawc", _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/reads-coru", _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/reads-corc", _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/blobs",      _storePath);  AS_UTL_unlink(path);

  AS_UTL_rmdir(_storePath);
}
