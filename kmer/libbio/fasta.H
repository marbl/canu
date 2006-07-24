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

  //  Computes the md5 checksum on the current sequence.  Performs an
  //  implicit getSequence().
  //
  //  XXX: Implementation is pretty ugly, lots of copied code from
  //  getSequence().  Maintenance nightmare.  Should use a sequence on
  //  disk, but that only works with squeezed sequences right now.
  //
  //void         computeMD5(md5_s &md5);

  md5_s *getMD5(u32bit iid) {
    return(_theMD5s + iid);
  };

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

  //  Returns true if a call to openIndex() will find a valid index to load.
  //  Returns false if openIndex() will need to build an index.
  //
  bool         isIndexValid(u32bit indextypetoload=FASTA_INDEX_ANY,
                            bool   beVerbose=false);

  //  Opens an existing index, or creates and opens an index if none exists.
  //
  void         openIndex(u32bit indextypetoload=FASTA_INDEX_ANY);

  //  Tell the wrapper that we are going to be doing mostly random access, and
  //  it should adjust it's buffer size appropriately.
  //
  void         optimizeRandomAccess(void);

  //  Various statistics about the file itself
  //
  u32bit       getNumberOfSequences(void);
  char const  *getSourceName(void);

  void         printATADescription(FILE *out, char *name);
  void         printTextDescription(FILE *out);

  bool         isRandomAccess(void)  { return(_isRandomAccess); };
  bool         isSqueezed(void)      { return(_theGlobalDesc._squeezedSequences); };
  bool         isChecksummed(void)   { return(_isChecksummed); };
private:

  //  Creates an index on disk
  void             createIndex(u32bit indextype);

  char            *indexTypeNames(u32bit indextype);

  char            *moveToNextName(char *ptr, u32bit iid);

  md5_s          *_theMD5s;
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
  bool            _isChecksummed;      //  Set if the index has md5's

  _idxfa_global   _theGlobalDesc;
};

#endif  //  FASTA_H

