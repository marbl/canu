
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
  sqStoreBlobWriter(const char *storePath, uint32 blobNumber) {

    //  Initialize us.

    strncpy(_storePath, storePath, FILENAME_MAX);

    _writtenBC   = blobNumber;
    _writtenBP   = 0;

    _bufferCount = blobNumber;
    _buffer      = NULL;

    //  Make a filename.

    makeName();

    //  Fail if it exists already.

    if (fileExists(_blobName) == true)
      fprintf(stderr, "sqStoreBlobWriter()-- blob file '%s' already exists.\n", _blobName), exit(1);

    //  Make a new write buffer;

    _buffer = new writeBuffer(_blobName, "w");
  };

  ~sqStoreBlobWriter() {
    delete _buffer;
  };


  void           makeName(void) {
    snprintf(_blobName, FILENAME_MAX, "%s/blobs.%04" F_U32P , _storePath, _bufferCount);
  };

  void           makeNextName(void) {
    _bufferCount++;
    makeName();
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

  uint32         writtenIndex(void)    { return(_writtenBC);       };
  uint64         writtenPosition(void) { return(_writtenBP);       };
  uint32         writtenBlob(void)     { return(_bufferCount + 1); };

private:
  char          _storePath[FILENAME_MAX+1];        //  Path to the seqStore.
  char          _blobName[FILENAME_MAX+1];         //  A temporary to make life easier.

  uint32        _writtenBC;                        //  The position before the
  uint64        _writtenBP;                        //  last writeData().

  uint32        _bufferCount;
  writeBuffer  *_buffer;
};



#endif  //  GKSTOREBLOBWRITER_H
