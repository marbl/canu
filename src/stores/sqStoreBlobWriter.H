
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

#ifndef GKSTOREBLOBWRITER_H
#define GKSTOREBLOBWRITER_H


class sqStoreBlobWriter {
public:
  sqStoreBlobWriter(const char *storePath) {

    //  Initialize us.

    strncpy(_storePath, storePath, FILENAME_MAX);

    _bufferCount = 0;
    _buffer      = NULL;

    //  Find the first available file.

    makeNextName(0);  //  Don't increment _bufferCount.
    while (AS_UTL_fileExists(_blobName) == true)
      makeNextName();

    //  Make a new write buffer;

    _buffer = new writeBuffer(_blobName, "w");
  };

  ~sqStoreBlobWriter() {
    delete _buffer;
  };


  void           makeNextName(uint32 next=1) {
    _bufferCount += next;
    snprintf(_blobName, FILENAME_MAX, "%s/blobs.%04" F_U32P , _storePath, _bufferCount);
  };

  void           writeData(uint8 *data, uint64 dataLen) {

    if (_buffer->tell() > AS_BLOBFILE_MAX_SIZE) {
      delete _buffer;

      makeNextName();

      _buffer = new writeBuffer(_blobName, "w");
    }

    _writtenBC = _bufferCount;
    _writtenBP = _buffer->tell();

    _buffer->write(data, dataLen);
  };

  uint32         writtenIndex(void)    { return(_writtenBC); };
  uint64         writtenPosition(void) { return(_writtenBP); };

private:
  char          _storePath[FILENAME_MAX+1];        //  Path to the seqStore.
  char          _blobName[FILENAME_MAX+1];         //  A temporary to make life easier.

  uint32        _writtenBC;                        //  The position before the
  uint64        _writtenBP;                        //  last writeData().

  uint32        _bufferCount;
  writeBuffer  *_buffer;
};



#endif  //  GKSTOREBLOBWRITER_H
