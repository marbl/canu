
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
bool
sqStore::sqStore_loadReadFromBuffer(readBuffer *B, sqRead *read) {

  //  If no buffer, or it's at the end, stop.

  if ((B == NULL) ||
      (B->eof() == true))
    return(false);

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

  return(true);
}



//  Dump the read metadata and read data to a stream.
//    rd must be allocated.  it is overwritten with read 'id's data.
//    wr must be allocated and uninitialized.
//
void
sqStore::sqStore_saveReadToBuffer(writeBuffer *B, uint32 id, sqRead *rd, sqReadDataWriter *wr) {
  sqReadSeq   emptySeq;   emptySeq.sqReadSeq_initialize();

  //  Write the read metadata.

  B->write(&_meta[id], sizeof(sqReadMeta));

  //  Write the sequence metadata, or an empty record if no metadata exists.

  B->write((_rawU) ? (&_rawU[id]) : (&emptySeq), sizeof(sqReadSeq));
  B->write((_rawC) ? (&_rawC[id]) : (&emptySeq), sizeof(sqReadSeq));
  B->write((_corU) ? (&_corU[id]) : (&emptySeq), sizeof(sqReadSeq));
  B->write((_corC) ? (&_corC[id]) : (&emptySeq), sizeof(sqReadSeq));

  //  Load (or reload) the sequence data, then write it out.

  sqStore_getRead(id, rd);

  rd->sqRead_fetchBlob(sqStore_getReadBuffer(id));
  rd->sqRead_decodeBlob();

  wr->sqReadDataWriter_importData(rd);
  wr->sqReadDataWriter_writeBlob(B);
}



sqLibrary *
sqStore::sqStore_addEmptyLibrary(char const *name, sqLibrary_tech techType) {

  assert(_info.sqInfo_lastLibraryID() <= _librariesAlloc);

  //  Just like with reads below, there is no _libraries[0] element, so we
  //  pre-increment the library number then allocate a new library.

  _info.sqInfo_addLibrary();

  uint32   libIdx = _info.sqInfo_lastLibraryID();

  increaseArray(_libraries, libIdx, _librariesAlloc, 128);

  //  Initialize the new library.

  _libraries[libIdx].sqLibrary_initialize(name, libIdx, techType);

  //  And return the freshly created library.

  return(_libraries + libIdx);
}



//  Allocate and return a new sqReadDataWriter that can be used to add a new
//  read to the store.  This function does NOT actually add the read to the
//  store; the sqReadDataWriter is used to collect all the info about the
//  read (read ID, read name, bases, quals, trim points) and then that object
//  is added to the store after all info is added.
//
//  Because the read isn't added until later, two consecutive calls to
//  createEmptyRead() will result in both sqReadDataWriter objects referring
//  to the same sqRead.
//
//  (The reason for this annoyance is so that sqStoreCreate can test that the
//  _homopoly_compressed_ length is big enough, and that length is only
//  computed by sqReadDataWriter::setRawBases().  Thus, we need to populate
//  the sqReadDataWriter object, test it, then discard it if too short.)
//
sqReadDataWriter *
sqStore::sqStore_createEmptyRead(sqLibrary *lib, const char *name) {

  assert(_info.sqInfo_lastReadID() < _readsAlloc);
  assert(_mode != sqStore_readOnly);

  //  The zeroth read is "null".  To get the next valid read ID, just
  //  add one to the current last read ID.

  uint32  rID = _info.sqInfo_lastReadID() + 1;
  uint32  lID = lib->sqLibrary_libraryID();

  //  Grow the metadata arrays if they're too small.

  if (_readsAlloc <= rID) {
    uint32  newMax = _readsAlloc + rID / 2;

    setArraySize(_meta, rID, _readsAlloc, newMax);
    setArraySize(_rawU, rID, _readsAlloc, newMax);
    setArraySize(_rawC, rID, _readsAlloc, newMax);
    setArraySize(_corU, rID, _readsAlloc, newMax);
    setArraySize(_corC, rID, _readsAlloc, newMax);
  }

  //  Initialize the metadata.

  _meta[rID].sqReadMeta_initialize(rID, lID);
  _rawU[rID].sqReadSeq_initialize();
  _rawC[rID].sqReadSeq_initialize();
  _corU[rID].sqReadSeq_initialize();
  _corC[rID].sqReadSeq_initialize();

  //  Make a new writer object, and initialize what we can.

  sqReadDataWriter  *rdw = new sqReadDataWriter(&_meta[rID],
                                                &_rawU[rID],
                                                &_rawC[rID],
                                                &_corU[rID],
                                                &_corC[rID]);

  rdw->sqReadDataWriter_setName(name);

  return(rdw);
}



//  Add a fully initialized read in sqReadDataWriter to the store.
//
//  This is inherently dangerous.  It assumes that the supplied sRDW is
//  actually for the next read.
//
//    If we're adding a new read - the ID _must_ be one more than the number
//    of reads in the store - increment the number of reads in the metadata.
//
//    If we're trying to add/modify a read that doesn't exist, blow up.
//
//    Otherwise, we're modifying an existing read, and need to only write the
//    read data.
//
void
sqStore::sqStore_addRead(sqReadDataWriter *rdw) {

  if        (rdw->_meta->sqRead_readID() == _info.sqInfo_lastReadID() + 1)
    _info.sqInfo_addRead();

  if (rdw->_meta->sqRead_readID() > _info.sqInfo_lastReadID()) {
    fprintf(stderr, "ERROR:  Attempt to add/modify read %u in a store with only %u reads.\n",
            rdw->_meta->sqRead_readID(), _info.sqInfo_lastReadID());
    assert(0);
  }

  _blobWriter->writeData(rdw);
}



void
sqStore::sqStore_setIgnored(uint32       id,
                            bool         untrimmed,
                            bool         trimmed,
                            sqRead_which w) {

  if (untrimmed) {
    sqStore_getReadSeq(id, w & ~sqRead_compressed)->sqReadSeq_setIgnoreU();
    sqStore_getReadSeq(id, w |  sqRead_compressed)->sqReadSeq_setIgnoreU();
  }

  if (trimmed) {
    sqStore_getReadSeq(id, w & ~sqRead_compressed)->sqReadSeq_setIgnoreT();
    sqStore_getReadSeq(id, w |  sqRead_compressed)->sqReadSeq_setIgnoreT();
  }
}



void
sqStore::sqStore_setClearRange(uint32 id,
                               uint32 bgn, uint32 end, bool bogus,
                               sqRead_which w) {
  sqRead_which   norm = w & ~sqRead_compressed;
  sqRead_which   comp = w |  sqRead_compressed;

  //  If we are bogus, just mark the normal and compressed sequences
  //  for ignore.

  if (bogus == true) {
    sqStore_setIgnored(id, false, true, w);
    return;
  }

  //  Grab the uncompressed sequence, build a map between that and the compressed
  //  sequence, then use the map to set clear ranges for both
  //  the normal and compressed.

  uint32    nlen  = sqStore_getReadLength(id, norm);

  sqRead   *read  = sqStore_getRead(id, new sqRead());
  uint32   *ntoc  = new uint32 [ nlen + 1 ];

  uint32    clen  = homopolyCompress(read->sqRead_sequence(norm), nlen, NULL, ntoc);

  assert(clen == sqStore_getReadLength(id, sqRead_corrected | sqRead_compressed));

  uint32    nbgn=0, nend=0;   //  Clear range in normal sequence
  uint32    cbgn=0, cend=0;   //  Clear range in compressed sequence

  //  Short clear ranges should be handled outside this function - by
  //  flagging the read to be ignored.  'end' == 0 when compressed is true
  //  causes an overflow in ntoc[nend].  Rather than clutter up the code with
  //  more special cases, we just require non-zero clear ranges.
  //
  assert(bgn < end);

  //  If we've got clear ranges for the normal version, we can directly
  //  map to the compressed version.
  //
  if ((w & sqRead_compressed) == sqRead_unset) {
    assert(end <= nlen);

    nbgn = bgn;
    nend = end;

    cbgn = ntoc[nbgn];
    cend = ntoc[nend];
  }

  //  But if we've got clear ranges for the compressed version, we need
  //  to invert the map before we can find the corresponding coordinates.
  //
  else {
    assert(end <= clen);

    cbgn = bgn;
    cend = end;

    nbgn = 0;
    nend = nlen;

    while ((nbgn < nend) && (ntoc[nbgn] < cbgn))
      nbgn++;

    while ((nbgn < nend) && (cend < ntoc[nend]))
      nend--;
  }

  //  Probably need real bounds checking.

  assert(nbgn <= nend);
  assert(nend <= nlen);

  assert(cbgn <= cend);
  assert(cend <= clen);

  sqStore_getReadSeq(id, norm)->sqReadSeq_setClearRange(nbgn, nend);
  sqStore_getReadSeq(id, comp)->sqReadSeq_setClearRange(cbgn, cend);

  fprintf(stderr, "id %5u length %7u %7u clear  norm %8u-%-8u  compressed %8u-%-8u\n",
          id, nlen, clen, nbgn, nend, cbgn, cend);

  delete [] ntoc;
  delete    read;
}
