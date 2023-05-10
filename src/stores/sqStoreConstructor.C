
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


//  Returns the last version available.
//             0 - no historical versions exist, but store is valid
//    UINT32_MAX - no store exists at that location
//
uint32
sqStore::sqStore_lastVersion(char const *storePath) {
  char    versionName[FILENAME_MAX+1];
  uint32  version = 0;

  if (directoryExists(storePath) == false)
    return(UINT32_MAX);

  snprintf(versionName, FILENAME_MAX, "%s/version.%03u/info", storePath, ++version);

  while (fileExists(versionName) == true)
    snprintf(versionName, FILENAME_MAX, "%s/version.%03u/info", storePath, ++version);

  return(version - 1);
}



//  Deletes versions higher than 'version' and
//  makes 'version' be the current.
void
sqStore::sqStore_revertVersion(char const *storePath, uint32 version) {
  char    versionPath[FILENAME_MAX+1];

  if (version == 0)
    return;

  snprintf(versionPath, FILENAME_MAX, "%s/version.%03u", storePath, version);

  if (fileExists(versionPath, '/', "info") == false) {
    fprintf(stderr, "Cannot revert to version %u; version not found.\n", version);
    return;
  }

  //  Unlink the current data.

  fprintf(stderr, "Removing current data.\n");

  AS_UTL_unlink(storePath, '/', "info");
  AS_UTL_unlink(storePath, '/', "info.txt");
  AS_UTL_unlink(storePath, '/', "libraries");
  AS_UTL_unlink(storePath, '/', "reads");
  AS_UTL_unlink(storePath, '/', "reads-corc");
  AS_UTL_unlink(storePath, '/', "reads-coru");
  AS_UTL_unlink(storePath, '/', "reads-rawc");
  AS_UTL_unlink(storePath, '/', "reads-rawu");

  //  Move the precious version to current.

  fprintf(stderr, "Restoring version %u.\n", version);

  AS_UTL_rename(versionPath, '/', "info",       storePath, '/', "info");
  AS_UTL_rename(versionPath, '/', "info.txt",   storePath, '/', "info.txt");
  AS_UTL_rename(versionPath, '/', "libraries",  storePath, '/', "libraries");
  AS_UTL_rename(versionPath, '/', "reads",      storePath, '/', "reads");
  AS_UTL_rename(versionPath, '/', "reads-corc", storePath, '/', "reads-corc");
  AS_UTL_rename(versionPath, '/', "reads-coru", storePath, '/', "reads-coru");
  AS_UTL_rename(versionPath, '/', "reads-rawc", storePath, '/', "reads-rawc");
  AS_UTL_rename(versionPath, '/', "reads-rawu", storePath, '/', "reads-rawu");

  AS_UTL_rmdir(versionPath);

  //  Remove all the now-obsolete versions.

  while (1) {
    snprintf(versionPath, FILENAME_MAX, "%s/version.%03u", storePath, ++version);

    if (fileExists(versionPath, '/', "info") == false)
      break;

    fprintf(stderr, "Removing version %u.\n", version);

    AS_UTL_unlink(versionPath, '/', "info");
    AS_UTL_unlink(versionPath, '/', "info.txt");
    AS_UTL_unlink(versionPath, '/', "libraries");
    AS_UTL_unlink(versionPath, '/', "reads");
    AS_UTL_unlink(versionPath, '/', "reads-corc");
    AS_UTL_unlink(versionPath, '/', "reads-coru");
    AS_UTL_unlink(versionPath, '/', "reads-rawc");
    AS_UTL_unlink(versionPath, '/', "reads-rawu");

    AS_UTL_rmdir(versionPath);
  }

  fprintf(stderr, "Done.\n");
}



void
sqStore::sqStore_loadMetadata(void) {

  _librariesAlloc = _info.sqInfo_lastLibraryID() + 1;
  _readsAlloc     = _info.sqInfo_lastReadID()    + 1;

  _libraries  = new sqLibrary  [_librariesAlloc];
  _meta       = new sqReadMeta [_readsAlloc];
  _rawU       = NULL;
  _rawC       = NULL;
  _corU       = NULL;
  _corC       = NULL;

  AS_UTL_loadFile(_metaPath, '/', "libraries", _libraries, _librariesAlloc);
  AS_UTL_loadFile(_metaPath, '/', "reads",     _meta,      _readsAlloc);

  //  If the user hasn't set a default version, or has only requested
  //  uncompressed reads, pick the plausible most recent version then add
  //  that to the default.
  //
  //  Otherwise, the user HAS requested a default.  If this isn't requesting
  //  'raw' or 'corrected' (so, e.g., only requesting the latest 'trimmed'
  //  reads), decide which of 'raw' or 'corrected' is the most recent.

  if ((sqRead_defaultVersion == sqRead_normal) ||
      (sqRead_defaultVersion == sqRead_unset)) {
    sqRead_which  mr;   //  "Most Recent"

    if (sqStore_getNumReads(sqRead_raw)                        > 0)   mr = sqRead_raw;
    if (sqStore_getNumReads(sqRead_raw | sqRead_trimmed)       > 0)   mr = sqRead_raw       | sqRead_trimmed;
    if (sqStore_getNumReads(sqRead_corrected)                  > 0)   mr = sqRead_corrected;
    if (sqStore_getNumReads(sqRead_corrected | sqRead_trimmed) > 0)   mr = sqRead_corrected | sqRead_trimmed;

    sqRead_defaultVersion |= mr;
  }

  else if (((sqRead_defaultVersion & sqRead_raw)       == sqRead_unset) &&
           ((sqRead_defaultVersion & sqRead_corrected) == sqRead_unset)) {
    if      (sqStore_getNumReads(sqRead_corrected) > 0)   sqRead_defaultVersion |= sqRead_corrected;
    else if (sqStore_getNumReads(sqRead_raw)       > 0)   sqRead_defaultVersion |= sqRead_raw;
  }

  //  If the store is set to return compressed reads by default and the user
  //  didn't explicitly say to use uncompressed reads, enable compression.

  if (((sqRead_defaultVersion & sqRead_normal) == sqRead_unset) &&
      (fileExists(sqStore_path(), '/', "homopolymerCompression") == true)) {
    sqRead_defaultVersion &= ~sqRead_normal;
    sqRead_defaultVersion |=  sqRead_compressed;
  }

  //  The default version MUST be set now:
  //   - either raw or corrected is set, but not both.
  //   - at most one of normal or compressed is set (but if neither is set, it's treated as 'normal').

  bool  validMode = true;

  if (sqRead_defaultIsNot(sqRead_raw) &&
      sqRead_defaultIsNot(sqRead_corrected)) {
    fprintf(stderr, "sqStore_loadMetadata()-- At least one of 'raw' or 'corrected' must be specified.\n");
    validMode = false;
  }

  if (sqRead_defaultIs(sqRead_raw) &&
      sqRead_defaultIs(sqRead_corrected)) {
    fprintf(stderr, "sqStore_loadMetadata()-- Cannot read both 'raw' and 'corrected' reads at the same time.\n");
    validMode = false;
  }

  if (sqRead_defaultIs(sqRead_normal) &&
      sqRead_defaultIs(sqRead_compressed)) {
    fprintf(stderr, "sqStore_loadMetadata()-- Cannot read both 'normal' and 'compressed' reads at the same time.\n");
    validMode = false;
  }

  assert(validMode == true);

  //  We can maybe, eventually, be clever and load these on-demand.
  //  For now, load all the metadata.

  AS_UTL_loadFile(_metaPath, '/', "reads-rawu", _rawU = new sqReadSeq[_readsAlloc], _readsAlloc);
  AS_UTL_loadFile(_metaPath, '/', "reads-rawc", _rawC = new sqReadSeq[_readsAlloc], _readsAlloc);
  AS_UTL_loadFile(_metaPath, '/', "reads-coru", _corU = new sqReadSeq[_readsAlloc], _readsAlloc);
  AS_UTL_loadFile(_metaPath, '/', "reads-corc", _corC = new sqReadSeq[_readsAlloc], _readsAlloc);
}






sqStore::sqStore(char const    *storePath_,
                 sqStore_mode   mode_,
                 uint32         version_) {

  //  Set the sqStore pointer in an overlap.

  ovOverlap::sqStoreAttach(this);

  //  Clear ourself, to make valgrind happier.

  assert(storePath_ != nullptr);

  memset(_storePath, 0, sizeof(char) * (FILENAME_MAX + 1));
  memset(_metaPath,  0, sizeof(char) * (FILENAME_MAX + 1));

  _mode                   = mode_;
  _version                = version_;

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

  //  Save the path to the store and make a metadata path name.

  strncpy(_storePath, storePath_, FILENAME_MAX);

  if (_version > 0)
    snprintf(_metaPath, FILENAME_MAX, "%s/version.%03u", sqStore_path(), _version);
  else
    snprintf(_metaPath, FILENAME_MAX, "%s", sqStore_path());

  //  Load the info file, if it exists.  And if not, do nothing.

  _info.readInfo(_metaPath);

  //
  //  CREATE - allocate some memory for saving libraries and reads, and create a file to dump the data into.
  //

  if (_mode == sqStore_create) {
    if (directoryExists(sqStore_path()) == true)
      fprintf(stderr, "ERROR:  Can't create store '%s': store already exists.\n", sqStore_path()), exit(1);

    AS_UTL_mkdir(sqStore_path());

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
    
    _blobWriter     = new sqStoreBlobWriter(sqStore_path(), &_info);

    return;
  }

  //
  //  READING or EXTENDING.
  //

  //  If the store directory doesn't exist, just blow up immediately.

  if (directoryExists(sqStore_path()) == false)
    fprintf(stderr, "sqStore()--  failed to open '%s' for read-only access: store doesn't exist.\n", sqStore_path()), exit(1);

  //  But if the metadata directory doesn't exist, silently return.  This
  //  allows sqStoreDumpMetaData to step through the versions and ask if
  //  there are reads present.  Any other client that wants to use metadata
  //  versions will need to do a similar check.

  if (directoryExists(_metaPath) == false)
    return;

  //  Now that we know the store and metadata directories exist, load the
  //  metadata, make a writer if we're extending, and make a reader.

  sqStore_loadMetadata();

  if (_mode == sqStore_extend)
    _blobWriter = new sqStoreBlobWriter(sqStore_path(), &_info);

  _blobReader = new sqStoreBlobReader(sqStore_path());
}






sqStore::~sqStore() {
  char    No[FILENAME_MAX+1];
  char    Nn[FILENAME_MAX+1];
  uint32  V = 1;

  //  Save original metadata.

  if (_mode == sqStore_extend) {
    snprintf(No, FILENAME_MAX, "%s/version.%03" F_U32P, sqStore_path(), V);
    while (directoryExists(No) == true) {
      V++;
      snprintf(No, FILENAME_MAX, "%s/version.%03" F_U32P, sqStore_path(), V);
    }

    AS_UTL_mkdir(No);

    snprintf(No, FILENAME_MAX, "%s/libraries", sqStore_path());
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/libraries", sqStore_path(), V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads", sqStore_path());
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads", sqStore_path(), V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads-rawu", sqStore_path());
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads-rawu", sqStore_path(), V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads-rawc", sqStore_path());
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads-rawc", sqStore_path(), V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads-coru", sqStore_path());
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads-coru", sqStore_path(), V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads-corc", sqStore_path());
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads-corc", sqStore_path(), V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/info", sqStore_path());
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/info", sqStore_path(), V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/info.txt", sqStore_path());
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/info.txt", sqStore_path(), V);
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
    AS_UTL_saveFile(sqStore_path(), '/', "libraries",   _libraries, sqStore_lastLibraryID() + 1);
    AS_UTL_saveFile(sqStore_path(), '/', "reads",       _meta,      sqStore_lastReadID()    + 1);
    AS_UTL_saveFile(sqStore_path(), '/', "reads-rawu",  _rawU,      sqStore_lastReadID()    + 1);
    AS_UTL_saveFile(sqStore_path(), '/', "reads-rawc",  _rawC,      sqStore_lastReadID()    + 1);
    AS_UTL_saveFile(sqStore_path(), '/', "reads-coru",  _corU,      sqStore_lastReadID()    + 1);
    AS_UTL_saveFile(sqStore_path(), '/', "reads-corc",  _corC,      sqStore_lastReadID()    + 1);

    _info.writeInfo(sqStore_path());

    FILE *F = AS_UTL_openOutputFile(sqStore_path(), '/', "info.txt");   //  Used by Canu/Gatekeeper.pm
    _info.writeInfoAsText(F);                                       //  Do not remove!
    AS_UTL_closeFile(F, sqStore_path(), '/', "info.txt");
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
