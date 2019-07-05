
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

#include "files.H"


sqRead_which    sqRead_defaultVersion = sqRead_unset;





//  Fetch the blob data from a readBuffer.  Do NOT position the buffer,
//  since this is used both for loading data from a store (random access)
//  and the correction/consensus 'package' files (sequential access).
void
sqRead::sqRead_fetchBlob(readBuffer *B) {

  B->readIFFchunk(_blobName, _blob, _blobLen, _blobMax);

  if (strncmp(_blobName, "BLOB", 4) != 0)
    fprintf(stderr, "Index error in read " F_U32 " mSegm " F_U64 " mByte " F_U64 " expected BLOB, got %02x %02x %02x %02x '%c%c%c%c'\n",
            _meta->sqRead_readID(),
            _meta->sqRead_mSegm(), _meta->sqRead_mByte(),
            _blobName[0], _blobName[1], _blobName[2], _blobName[3],
            _blobName[0], _blobName[1], _blobName[2], _blobName[3]), exit(1);
}


//  Return a readBuffer, correctly positioned, to load data for read 'readID'.
readBuffer *
sqStore::sqStore_getReadBuffer(uint32 readID) {
  readBuffer *buffer = _blobReader->getBuffer(_meta[readID]);

  buffer->seek(_meta[readID].sqRead_mByte());

  return(buffer);
}



//  Set pointers to the metadata, forget whatever sequence we're
//  remembering, and (optionally) load bases from the blob.
//
sqRead *
sqStore::sqStore_getRead(uint32 readID, sqRead *read) {

  read->_meta     =           (_meta + readID);
  read->_rawU     = (_rawU) ? (_rawU + readID) : (NULL);
  read->_rawC     = (_rawC) ? (_rawC + readID) : (NULL);
  read->_corU     = (_corU) ? (_corU + readID) : (NULL);
  read->_corC     = (_corC) ? (_corC + readID) : (NULL);

  read->_library  = sqStore_getLibrary(read->_meta->sqRead_libraryID());

  read->_retFlags = 0;

  if (true) {
    read->sqRead_fetchBlob(sqStore_getReadBuffer(readID));
    read->sqRead_decodeBlob();
  }

  return(read);
}



//  Load read metadata and data from a stream.
//
void
sqStore::sqStore_loadReadFromBuffer(readBuffer *B, sqRead *read) {

  //  Load the read and sequence metadata.

  if (read->_metaA == NULL) {
    read->_metaA = new sqReadMeta [1];
    read->_rseqA = new sqReadSeq  [4];

    read->_meta = read->_metaA;
    read->_rawU = read->_rseqA + 0;
    read->_rawC = read->_rseqA + 1;
    read->_corU = read->_rseqA + 2;
    read->_corC = read->_rseqA + 3;
  }

  B->read(read->_meta, sizeof(sqReadMeta));
  B->read(read->_rawU, sizeof(sqReadSeq));
  B->read(read->_rawC, sizeof(sqReadSeq));
  B->read(read->_corU, sizeof(sqReadSeq));
  B->read(read->_corC, sizeof(sqReadSeq));

  read->_library = NULL;

  //  Load the read sequence data.

  read->sqRead_fetchBlob(B);
  read->sqRead_decodeBlob();
}



//  Dump the read metadata and read data to a stream.
//
void
sqStore::sqStore_saveReadToBuffer(writeBuffer *B, uint32 id, sqRead *rd, sqReadDataWriter *wr) {
  sqReadSeq   emptySeq;

  //  Write the read metadata.

  B->write(&_meta[id], sizeof(sqReadMeta));

  //  Write the sequence metadata, or an empty record if no metadata exists.

  B->write((_rawU) ? (&_rawU[id]) : (&emptySeq), sizeof(sqReadSeq));
  B->write((_rawC) ? (&_rawC[id]) : (&emptySeq), sizeof(sqReadSeq));
  B->write((_corU) ? (&_corU[id]) : (&emptySeq), sizeof(sqReadSeq));
  B->write((_corC) ? (&_corC[id]) : (&emptySeq), sizeof(sqReadSeq));

  //  Load (or reload) the sequence data, then write it out.

  rd->sqRead_fetchBlob(sqStore_getReadBuffer(id));
  rd->sqRead_decodeBlob();

  wr->sqReadDataWriter_importData(rd);
  wr->sqReadDataWriter_writeBlob(B);
}



sqLibrary *
sqStore::sqStore_addEmptyLibrary(char const *name) {

  assert(_info.sqInfo_lastLibraryID() <= _librariesAlloc);

  //  Just like with reads below, there is no _libraries[0] element.

  _info.sqInfo_addLibrary();

  increaseArray(_libraries, _info.sqInfo_lastLibraryID(), _librariesAlloc, 128);

  //  Initialize the new library.

  _libraries[_info.sqInfo_lastLibraryID()] = sqLibrary();
  _libraries[_info.sqInfo_lastLibraryID()]._libraryID = _info.sqInfo_lastLibraryID();

  //  Bullet proof the library name - so we can make files with this prefix.

  char   *libname    = _libraries[_info.sqInfo_lastLibraryID()]._libraryName;
  uint32  libnamepos = 0;

  memset(libname, 0, sizeof(char) * LIBRARY_NAME_SIZE);

  for (char const *orig=name; *orig; orig++) {
    if        (*orig == '/') {
      libname[libnamepos++] = '_';

    } else if (isspace(*orig) == 0) {
      libname[libnamepos++] = *orig;

    } else {
      libname[libnamepos++] = '_';
    }

    if (libnamepos >= LIBRARY_NAME_SIZE) {
      libname[LIBRARY_NAME_SIZE-1] = 0;
      break;
    }
  }

  return(_libraries + _info.sqInfo_lastLibraryID());
}




sqReadDataWriter *
sqStore::sqStore_addEmptyRead(sqLibrary *lib) {

  assert(_info.sqInfo_lastReadID() < _readsAlloc);
  assert(_mode != sqStore_readOnly);

  //  We reserve the zeroth read for "null".  This is easy to accomplish
  //  here, just pre-increment the number of reads.  However, we need to be sure
  //  to iterate up to and including _info.sqInfo_lastReadID().

  _info.sqInfo_addRead();

  if (_readsAlloc <= _info.sqInfo_lastReadID()) {
    uint32  newMax = _readsAlloc + _info.sqInfo_lastReadID() / 2;

    setArraySize(_meta, _info.sqInfo_lastReadID(), _readsAlloc, newMax);
    setArraySize(_rawU, _info.sqInfo_lastReadID(), _readsAlloc, newMax);
    setArraySize(_rawC, _info.sqInfo_lastReadID(), _readsAlloc, newMax);
    setArraySize(_corU, _info.sqInfo_lastReadID(), _readsAlloc, newMax);
    setArraySize(_corC, _info.sqInfo_lastReadID(), _readsAlloc, newMax);
  }

  //  Initialize the new read.

  uint32  rID = _info.sqInfo_lastReadID();
  uint32  lID = lib->sqLibrary_libraryID();

  _meta[rID] = sqReadMeta(rID, lID);

  //  With the read set up, set pointers in the readData.  Whatever data is in there can stay.

  return(new sqReadDataWriter(&_meta[rID],
                              &_rawU[rID],
                              &_rawC[rID],
                              &_corU[rID],
                              &_corC[rID]));
}
