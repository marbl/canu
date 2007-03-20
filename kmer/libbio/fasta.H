#ifndef FASTA_H
#define FASTA_H

#ifndef BIO_PLUS_PLUS_H
#error Please use bio++.H instead of fasta.H
#endif

#include "fasta-types.H"
#include "fasta-sequenceondisk.H"
#include "fasta-sequenceincore.H"

class FastAWrapper {
public:
  FastAWrapper(char const *filename, u32bit bufferSize=~u32bitZERO);
  ~FastAWrapper();

  //  Return some accessors for whatever sequence we are currently at.
  //
  FastASequenceInCore *getSequence(void);
  FastASequenceOnDisk *getSequenceOnDisk(void);

  //  Return information about the sequence, only if we have an index.
  //  Otherwise, returns 0.
  //
  u32bit       sequenceLength(IID_t iid);
  u32bit       headerLength(IID_t iid);
  u32bit       currentIID(void);

  //  Sequential access to sequences (also use getSequence*()).
  //
  bool         eof(void) { return(_filebuffer->eof()); };

  //  Random access to sequences.
  //
  bool         find(IID_t  iid);
  bool         find(char  *id);

  //  Returns true if a call to openIndex() will find a valid index to
  //  load.  Returns false if openIndex() will need to build an index.
  //
  bool         isIndexValid(u32bit indextypetoload=FASTA_INDEX_ANY,
                            bool   beVerbose=false);

  //  Opens an existing index, or creates and opens an index if none
  //  exists.
  //
  void         openIndex(u32bit indextypetoload=FASTA_INDEX_ANY);

  //  Tell the wrapper that we are going to be doing mostly random
  //  access, and it should adjust it's buffer size appropriately.
  //
  void         optimizeRandomAccess(void);

  u32bit       getNumberOfSequences(void);
  char const  *getSourceName(void);

  void         printDescription(FILE *out, char *name);

  bool         isRandomAccess(void)  { return(_isRandomAccess); };
  bool         isSqueezed(void)      { return(_theGlobalDesc._squeezedSequences); };

private:
  void             createIndex(u32bit indextype);
  const char      *indexTypeNames(u32bit indextype);
  char            *moveToNextName(char *ptr, u32bit iid);

  _idxfa_desc    *_theSeqs;
  u32bit          _theNamesLen;
  char           *_theNames;

  readBuffer     *_filebuffer;
  char           *_filename;
  char           *_indexname;

  u32bit          _currentSequenceNumber;

  bool            _isStreamInput;      //  Forced to be not random access
  bool            _isRandomAccess;     //  Set if an index is open
  bool            _isRandomAccessOpt;  //  Set if we resized the readBuffer

  _idxfa_global   _theGlobalDesc;
};

#endif  //  FASTA_H

