#include "bigQueue.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

//  Kaz Kylheku <kaz@ashi.footprints.net> library.
#include "kazlib/dict.h"
#include "kazlib/except.h"
#include "kazlib/hash.h"
#include "kazlib/list.h"
#include "kazlib/sfx.h"

//  qsort and kazlib are incombatible.  qsort passes a pointer to the data, kaz lib passes
//  the data (which it assumes is a pointer to begin with).


void
bigQueue::_initialize(int    (*sortfcn)(const void *a, const void *b),
                      bool   (*readfcn)(FILE *f, void *a),
                      bool   (*writfcn)(FILE *f, void *a),
                      void   (*killfcn)(void *a),
                      uint32   objectSize,
                      uint32   memoryToUse,
                      char    *tmppath,
                      char    *filename) {
  _saveFile     = 0L;
  _tmpPath      = 0L;

  if (filename) {
    _saveFile = new char [strlen(filename) + 1];
    strcpy(_saveFile, filename);
  }
  if (tmppath) {
    _tmpPath = new char [strlen(tmppath) + 1];
    strcpy(_tmpPath, tmppath);
  }

  _sortFunction      = sortfcn;
  _writFunction      = writfcn;
  _readFunction      = readfcn;
  _killFunction      = killfcn;

  _objectSize        = objectSize;
  _memoryToUse       = memoryToUse;

  _maxOpenFiles      = getdtablesize() - 8;
  _numTemporaryFiles = 0;
  _numMergeFiles     = 0;

  _temporaryFiles = new FILE* [_maxOpenFiles];

  for (uint32 i=0; i<_maxOpenFiles; i++)
    _temporaryFiles[i] = 0L;

  //  Open the first temporary file for writing.
  //
  _temporaryFiles[_numTemporaryFiles++] = makeTempFile(_tmpPath);

  //  XXX: It would be rather convenient if we could get another file
  //  handle given an existing handle (no, dup(2) doesn't do that).
  //  In particular, we want two file pointers, one for read, one for
  //  write.
  //
  //_inputFile  = fdopen(dup(fileno(_temporaryFiles[0])), "w+");

  _thingBuffer = new uint64 [_objectSize / 8 + 1];

  _bufferMax = 0;
  _bufferLen = 0;
  _buffer    = 0L;

  if (_sortFunction) {
    _bufferMax = (uint64)memoryToUse * 1024 * 1024 / ((uint64)sizeof(void *) + objectSize);
    _bufferLen = 0;
    _buffer    = new void* [_bufferMax];
  }
}



bigQueue::~bigQueue() {
  delete [] _saveFile;
  delete [] _tmpPath;

  for (uint32 i=0; i<_numTemporaryFiles; i++)
    fclose(_temporaryFiles[i]);

  delete [] _temporaryFiles;

  //fclose(_inputFile);

  clearBuffer();
}






//  Add elements to the end of the array.
void    
bigQueue::add(void *thing) {

  if (_buffer == 0L) {
    if (_writFunction)
      (*_writFunction)(_temporaryFiles[_numTemporaryFiles-1], thing);
    else
      fwrite(thing, _objectSize, 1, _temporaryFiles[_numTemporaryFiles-1]);
  } else {

    //  No space in the buffer?  Sort it, write it out and make a new
    //  one.
    //
    if (_bufferLen >= _bufferMax) {
      sortAndWriteBuffer();

      if (_numTemporaryFiles+1 >= _maxOpenFiles)
        mergeTemporaryFiles();

      _temporaryFiles[_numTemporaryFiles++] = makeTempFile(_tmpPath);
    }

    _buffer[_bufferLen++] = thing;
  }
}



void
bigQueue::sortAndWriteBuffer(void) {

  if (_bufferLen > 0) {

    //  Sort!
    //
    qsort(_buffer, _bufferLen, sizeof(void *), _sortFunction);

    //  Write!
    //
    if (_writFunction) {
      for (uint32 i=0; i<_bufferLen; i++)
        (*_writFunction)(_temporaryFiles[_numTemporaryFiles-1], _buffer[i]);
    } else {
      for (uint32 i=0; i<_bufferLen; i++)
        fwrite(_buffer[i], _objectSize, 1, _temporaryFiles[_numTemporaryFiles-1]);
    }

    //  Flush and rewind the file!
    //
    fflush(_temporaryFiles[_numTemporaryFiles-1]);
    ::rewind(_temporaryFiles[_numTemporaryFiles-1]);

    clearBuffer();
  }
}


void
bigQueue::clearBuffer(void) {

  if (_killFunction)
    for (uint32 i=0; i<_bufferLen; i++)
      (*_killFunction)(_buffer[i]);
  else
    for (uint32 i=0; i<_bufferLen; i++)
      free(_buffer[i]);

  _bufferLen = 0;
}


void
bigQueue::mergeTemporaryFiles(void) {

  if (_numTemporaryFiles > 1) {
    dict_t    *sorted;
    dnode_t   *nodes = new dnode_t [_maxOpenFiles];

    //  To be efficient, we need to maintain a sorted queue of the head
    //  elements of each temporary file.  A red-black tree would do
    //  nicely, eh?
    //
    sorted = dict_create(DICTCOUNT_T_MAX, _sortFunction);

    //  Grab the first thing off each file, insert it into the dictionary.
    //  The 'key' is our chunk of data, and the 'value' is the file number
    //  it came from.
    //
    for (uint32 i=0; i<_numTemporaryFiles; i++) {
      if (_temporaryFiles[i]) {

        //  Rewind all the temporary files.  XXXX This is probably done
        //  already.
        //
        ::rewind(_temporaryFiles[i]);

        void *thing = malloc(_objectSize);

        if (_readFunction)
          (*_readFunction)(_temporaryFiles[i], thing);
        else
          fread(thing, _objectSize, 1, _temporaryFiles[i]);

        if (feof(_temporaryFiles[i])) {
          fclose(_temporaryFiles[i]);
          _temporaryFiles[i] = 0L;
        } else {
          //  initialize the node with the value
          dnode_init(&nodes[i], (void *)(unsigned long)i);

          //  insert the node into the tree using the key
          dict_insert(sorted, &nodes[i], thing);
        }
      }
    }

    FILE *mergeFile = makeTempFile(_tmpPath);

    //  while there is stuff in the tree

    while (dict_isempty(sorted) == 0) {

      //  pop the head element off, and print it
      dnode_t  *head = dict_first(sorted);

      //  XXX: should be const thing

      void   *thing  = (void *)dnode_getkey(head);
      long    fileid = (long)dnode_get(head);

      if (_writFunction)
        (*_writFunction)(mergeFile, thing);
      else
        fwrite(thing, _objectSize, 1, mergeFile);

      //  delete the node from the tree
      dict_delete(sorted, head);

      //  destroy the thing
      if (_killFunction)
        (*_killFunction)(thing);
      else
        free(thing);

      //  load the next element from the same file that the head was
      //  from (that's stored as the value of the head element)

      thing = malloc(_objectSize);

      if (_readFunction)
        (*_readFunction)(_temporaryFiles[fileid], thing);
      else
        fread(thing, _objectSize, 1, _temporaryFiles[fileid]);

      //  if there was a next element in that file, insert it
      //  into the tree.  if not, close the temporary file.
      //
      if (feof(_temporaryFiles[fileid])) {
        fclose(_temporaryFiles[fileid]);
        _temporaryFiles[fileid] = 0;
        free(thing);
      } else {
        //  initialize the node with the value
        dnode_init(&nodes[fileid], (void *)fileid);

        //  insert the node into the tree using the key
        dict_insert(sorted, &nodes[fileid], thing);
      }
    }

    dict_free(sorted);
    delete [] nodes;

    _numTemporaryFiles = 1;
    _temporaryFiles[0] = mergeFile;
  }

  ::rewind(_temporaryFiles[0]);

#if 0
  fclose(_inputFile);
  errno = 0;
  _inputFile  = fdopen(dup(fileno(_temporaryFiles[0])), "w+");
  if (errno)
    fprintf(stderr, "bigQueue::mergeTemporaryFiles()-- _inputFile = fdopen() failed: %s\n", strerror(errno)), exit(1);

  ::rewind(_inputFile);
#endif
}


bool
bigQueue::next(void) {

  if (_readFunction) {
    //(*_readFunction)(_inputFile, _thingBuffer);
    (*_readFunction)(_temporaryFiles[0], _thingBuffer);
  } else {
    //fread(_thingBuffer, _objectSize, 1, _inputFile);
    fread(_thingBuffer, _objectSize, 1, _temporaryFiles[0]);
  }

#if 0
  if (feof(_inputFile))
    return(false);
#endif

  if (feof(_temporaryFiles[0]))
    return(false);

  return(true);
}


void*
bigQueue::get(void) {
  return(_thingBuffer);
}

void
bigQueue::rewind(void) {
  //::rewind(_inputFile);
  ::rewind(_temporaryFiles[0]);
  next();
}

void
bigQueue::save(char *filepath) {
  fprintf(stderr, "bigQueue::save()-- not implemented.\n");
}

void
bigQueue::sort(void) {
  sortAndWriteBuffer();
  mergeTemporaryFiles();
}

void
bigQueue::flush(void) {
  fflush(_temporaryFiles[_numTemporaryFiles-1]);
}
