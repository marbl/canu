
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
 *    src/stores/gkStoreBlobWriter.H
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

#include "sqStore.H"
#include "files.H"
#include "objectStore.H"




//  Scan the 'storePath' for the next free blob number, and set up
//  to write new data there.
//
sqStoreBlobWriter::sqStoreBlobWriter(const char *storePath, sqStoreInfo *info) {

  memset(_storePath, 0, sizeof(char) * FILENAME_MAX);    //  Clear the path.
  memset(_blobName,  0, sizeof(char) * FILENAME_MAX);    //  Clear the name.

  strncpy(_storePath, storePath, FILENAME_MAX);          //  Copy path to our path.

  _info       =  info;                                   //  Remember the info!
  _blobNumber = _info->_numBlobs++;                      //  Set up for the next avail blob.

  makeBlobName(_storePath, _blobNumber, _blobName);      //  Construct the name of the blob file.

  _buffer     = new writeBuffer(_blobName, "w");         //  And open it.
}


sqStoreBlobWriter::~sqStoreBlobWriter() {
  delete _buffer;
}



//  Add the data in sqReadDataWriter to the current blob file, starting
//  a new blob if the current one is too big.
//
void
sqStoreBlobWriter::writeData(sqReadDataWriter *rdw) {

  if (_buffer->tell() > AS_BLOBFILE_MAX_SIZE) {
    delete _buffer;

    makeBlobName(_storePath, ++_blobNumber, _blobName);

    _buffer = new writeBuffer(_blobName, "w");
  }

  //  Save the current position in the blob file in the sqStore
  //  metadata, then tell the rdw to dump data.

  rdw->_meta->sqRead_setPosition(_blobNumber, _buffer->tell());
  rdw->sqReadDataWriter_writeBlob(_buffer);
}



sqStoreBlobReader::sqStoreBlobReader(const char *storePath) {

  memset(_storePath, 0, sizeof(char) * FILENAME_MAX);
  memset(_blobName,  0, sizeof(char) * FILENAME_MAX);

  strncpy(_storePath, storePath, FILENAME_MAX);

  _buffersMax = 0;
  _buffers    = NULL;

  resizeArray(_buffers, _buffersMax, _buffersMax, 128, resizeArray_copyData | resizeArray_clearNew);
}



sqStoreBlobReader::~sqStoreBlobReader() {
  for (uint32 ii=0; ii<_buffersMax; ii++)
    delete _buffers[ii];
  delete [] _buffers;
}



readBuffer *
sqStoreBlobReader::getBuffer(sqReadMeta *meta) {
  uint32  file = meta->sqRead_mSegm();
  uint32  posn = meta->sqRead_mByte();

  while (_buffersMax <= file)
    resizeArray(_buffers, _buffersMax, _buffersMax, _buffersMax * 2, resizeArray_copyData | resizeArray_clearNew);

  if (_buffers[file] == NULL) {
    makeBlobName(_storePath, file, _blobName);

#pragma omp critical
    //  Fetch from object store, if needed and possible.
    fetchFromObjectStore(_blobName);

    _buffers[file] = new readBuffer(_blobName, 1024 * 1024);
  }

  _buffers[file]->seek(posn);

  return(_buffers[file]);
}

