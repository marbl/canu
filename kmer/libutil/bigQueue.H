#ifndef BIGQUEUE_H
#define BIGQUEUE_H

#include "util++.H"

//  A disk-backed list of user-defined objects.
//
//  At creation time, you can opt to have it sorted, using a
//  user-defined function.

//  An list based on a variable length object (let alone a sort!)
//  must use some form of dereferencing scheme.  So, if you want to
//  use variable length records, you have to use pointers, and supply
//  functions to do everything (compare, read, write).
//
//  On the otherhand, it would be quite more convenient (to use) if we
//  used objects (would need copy, compare, read, write).
//
//  1) Restrict to void*, fixed block size, functions for compare,
//  destroy.  read, write and copy done with fread(), fwrite() and
//  memcpy().
//
//  2) Restrict to void*, functions for compare, read, write and
//  destroy.  I allocate an array of pointers.  Assume shallow copies
//  are ok (qsort will be used).  On construct, we need to know the
//  size of the data so we know how many objects to buffer before
//  sorting and writing.  It's also possible to use fread() and
//  fwrite().
//
//  3) Restrict to objects, operators for copy, compare, read, write,
//  default construct, destroy.  I allocate an array of objects.
//
//  1 is the easiest to write, 2 and 3 are conceptually the same.  1
//  cannot write out deep data (pointer to string).  2 is a trivial
//  extenstion to 1, and fixes that.  3 is the correct version, but I
//  don't want to deal with streams io.  So, 2 it is.
//

class bigQueue {
public:
  //  Initialize the bigQueue for anonymous storage, with an
  //  option to later save the array.
  //
  bigQueue(bool   (*readfcn)(FILE *, void *),
           bool   (*writfcn)(FILE *, void *),
           void   (*killfcn)(void *),
           uint32   objectSize,
           char    *tmpPath) {
    _initialize(0L, readfcn, writfcn, killfcn, objectSize, 0, tmpPath, 0L);
  };


  //  Initialize the bigQueue with a file of objects, presumabely from
  //  a previous invocation of bigQueue.
  //
  bigQueue(bool   (*readfcn)(FILE *, void *),
           bool   (*writfcn)(FILE *, void *),
           void   (*killfcn)(void *),
           uint32   objectSize,
           char    *tmpPath,
           char    *filename) {
    _initialize(0L, readfcn, writfcn, killfcn, objectSize, 0, tmpPath, filename);
  };


  //  Initialize the bigQueue for sorting.
  //
  bigQueue(int    (*sortfcn)(const void *a, const void *b),
           bool   (*readfcn)(FILE *, void *),
           bool   (*writfcn)(FILE *, void *),
           void   (*killfcn)(void *),
           uint32   objectSize,
           uint32   memoryToUse,
           char    *tmpPath) {
    _initialize(sortfcn, readfcn, writfcn, killfcn, objectSize, memoryToUse, tmpPath, 0L);
  };

private:
  void      _initialize(int    (*sortfcn)(const void *a, const void *b),
                        bool   (*readfcn)(FILE *f, void *a),
                        bool   (*writfcn)(FILE *f, void *a),
                        void   (*killfcn)(void *),
                        uint32   objectSize,
                        uint32   memoryToUse,
                        char    *tmppath,
                        char    *filename);

public:
  ~bigQueue();

  //  Add elements to the end of the array.
  void    add(void *);

  //  We are designed for streaming access.
  bool    next(void);
  void   *get(void);

  //  Rewind to the start.  Sortable must be sorted.
  void    rewind(void);

  //  Save the anonymous array into a real file.
  void    save(char *filepath);

  //  Sort the sortable.  Flush the flushable.
  void    sort(void);
  void    flush(void);

private:
  void       sortAndWriteBuffer(void);
  void       clearBuffer(void);
  void       mergeTemporaryFiles(void);

  char     *_saveFile;
  char     *_tmpPath;

  int     (*_sortFunction)(const void *a, const void *b);
  bool    (*_writFunction)(FILE *f, void *a);
  bool    (*_readFunction)(FILE *f, void *a);
  void    (*_killFunction)(void *a);

  uint32    _objectSize;
  uint32    _memoryToUse;

  uint32    _maxOpenFiles;
  uint32    _numTemporaryFiles;
  uint32    _numMergeFiles;

  //  _temporaryFiles is all the opened output files.  If we aren't
  //  sorting, then only the first one is opened.
  //
  //  _inputFile is a dup of the first temporary file.  If we are
  //  sorting, and you start reading before you sort, then you'll get
  //  a very short read.
  //
  FILE    **_temporaryFiles;
  FILE     *_inputFile;

  //  Stores things read back from disk, for return to the user.
  //  Currently just one, but should be extended to many.
  //
  void     *_thingBuffer;

  uint32    _bufferMax;
  uint32    _bufferLen;
  void    **_buffer;
};



#endif  //  BIGQUEUE_H
