
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
 *    src/utility/readBuffer.H
 *    src/utility/writeBuffer.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2019-FEB-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef FILES_BUFFERED_H
#define FILES_BUFFERED_H

//  Do not include directly.  Use 'files.H' instead.

class readBuffer {
public:
  readBuffer(const char *prefix, char separator, const char *suffix,
             uint64      bufferMax = 32 * 1024);
  readBuffer(const char *filename,
             uint64      bufferMax = 32 * 1024);
  readBuffer(FILE *F,
             uint64      bufferMax = 32 * 1024);
  ~readBuffer();

private:
  void                 initialize(const char *filename,
                                  uint64      bufferMax);
public:
  //  Returns true if the buffer is at the end.
  //
  //  Note that if the last read was successful, but we're now at the end of
  //  file, this is true.
  bool                 eof(void) { return(_eof); };

  char                 peek(void);

  char                 read(void);
  char                 readuntil(char stop);
  uint64               read(void *buf, uint64 len);
  uint64               read(void *buf, uint64 maxlen, char stop);

  bool                 peekIFFchunk(char *name, uint32 &dataLen);


  //  The next four will read a specific chunk into a pre-allocated buffer.
  template<typename OBJ>
  bool                 readIFFarray(char const *name, OBJ *array, uint32 nObjects) {
    return(readIFFchunk(name, array, sizeof(OBJ) * nObjects));
  };

  template<typename OBJ>
  bool                 readIFFobject(char const *name, OBJ &object) {
    return(readIFFchunk(name, &object, sizeof(OBJ)));
  };

  template<typename OBJ>
  bool                 readIFFobject(char const *name, OBJ *object) {
    return(readIFFchunk(name, object, sizeof(OBJ)));
  };

  bool                 readIFFchunk(char const *name,
                                    void       *data,
                                    uint32      dataLen);

  //  Read a generic chunk into re-allocatable memory.
  void                 readIFFchunk(char   *name,
                                    uint8 *&data,
                                    uint32 &dataLen,
                                    uint32 &dataMax);

  void                 skipAhead(char stop, bool after=false);
  uint64               copyUntil(char stop, char *dest, uint64 destLen);

  void                 seek(uint64 pos, uint64 extra=0);
  uint64               tell(void) { return(_filePos); };

  const char          *filename(void) { return(_filename); };

private:
  void                 fillBuffer(uint64 extra=0);
  void                 init(int fileptr, const char *filename, uint64 bufferMax);

  char                _filename[FILENAME_MAX+1];

  int                 _file;        //
  uint64              _filePos;     //  Position in the file we're at.

  bool                _stdin;

  bool                _eof;

  bool                _ignoreCR;    //  Ignore blasted DOS CR letters in read() and readuntil().

  uint64              _bufferBgn;   //  File position where this buffer is from.
  uint64              _bufferLen;   //  Length of the valid data in the buffer.

  uint64              _bufferPos;   //  Position in the buffer we're at.

  uint64              _bufferMax;   //  Size of _buffer allocation.
  char               *_buffer;      //  Data!
};



class writeBuffer {
public:
  writeBuffer(const char *prefix, char separator, const char *suffix,
              const char *filemode,
              uint64      bufferMax = 1024 * 1024);
  writeBuffer(const char *filename,
              const char *filemode,
              uint64      bufferMax = 1024 * 1024);
  ~writeBuffer();

private:
  void                 initialize(const char *filename,
                                  const char *filemode,
                                  uint64      bufferMax);

public:
  const char          *filename(void) { return(_filename); };
  uint64               tell(void)     { return(_filePos);  };

  void                 write(void const *data, uint64 length);
  void                 flush(void);

  //  If deata length is zero, this chunk is the start of a recursive record.
  //  The chunk header is written (NAME<length>) but length is set to zero.
  //  The position of the length field is is pushed onto an internal stack.
  //  When the chunk is closed, the length field is updated with
  //  the size of blocks added to it.
  //
  //  To get an empty chunk:
  //    writeIFFchunk("NAME", NULL, 0);
  //    closeIFFchunk();
  //
private:
  void   appendIFFdata(void *data, uint32 dataLength);

public:
  template<typename OBJ>
  void   writeIFFarray(char *name, OBJ *array, uint32 nObjects) {
    writeIFFchunk(name, array, sizeof(OBJ) * nObjects);
  };

  template<typename OBJ>
  void   writeIFFobject(char *name, OBJ &object) {
    writeIFFchunk(name, &object, sizeof(OBJ));
  };

  template<typename OBJ>
  void   writeIFFobject(char *name, OBJ *object) {
    writeIFFchunk(name, object, sizeof(OBJ));
  };

  void   writeIFFchunk(char *name, void *data=NULL, uint32 dataLength=0);

  void   closeIFFchunk(char *name=NULL);

private:
  void                 open(void);
  void                 writeToDisk(void const *data, uint64 length);

  char                _filename[FILENAME_MAX+1];
  char                _filemode[17];

  FILE               *_file;
  uint64              _filePos;

  uint64              _bufferLen;
  uint64              _bufferMax;
  char               *_buffer;

  uint64              _chunkBufferLen;         //  For building up recursive chunks,
  uint64              _chunkBufferMax;         //  another buffer of data.
  uint8              *_chunkBuffer;

  uint32              _chunkStartsLen;   //  A stack of start positions.
  uint32              _chunkStartsMax;
  uint64             *_chunkStarts;
  uint64             *_chunkSizes;
};



#endif  //  FILES_BUFFERED_H
