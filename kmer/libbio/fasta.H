#ifndef FASTA_H
#define FASTA_H

#ifndef BRI_PLUS_PLUS_H
#error Please use bri++.H instead of fasta.H
#endif

//#include "bri++.H"

#ifndef FILE_WRAPPER_H
#define FILE_WRAPPER_H

//  If we don't have RASCAL, we need to define IID_t
//
typedef u32bit IID_t;
typedef u32bit SeqOff_t;

class FileWrapper {
public:
  FileWrapper() {};
  virtual ~FileWrapper() {};

  //  Sequential access to sequences.
  //
  virtual bool        eof(void) =0;

  //  Random access to sequences.
  //
  virtual bool        find(IID_t  iid) =0;  //  throw(Ex)
  virtual bool        find(char  *id) =0;  //   throw(Ex)

  //  Various statistics about the file itself
  //
  virtual u32bit      getNumberOfSequences(void) =0;
  virtual const char *getSourceName(void) =0;
};

#endif  //  FILE_WRAPPER_H



#include "fasta-types.H"
#include "fasta-sequenceondisk.H"
#include "fasta-sequenceincore.H"


//  Version Numbers:
//
//  versionNumber 1 contains the simple index, name or defline.
//
//  versionNumber 2 extends that to contain md5 checksums.  A version 1
//  index is binary compatible with version 2.
//
#define FASTA_MAGICNUMBER   0x4abeadeddadababeULL
#define FASTA_VERSIONNUMBER 2


class FastAWrapper : public FileWrapper {
public:
  FastAWrapper(const char *filename);
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
  bool         isIndexValid(u32bit indextypetoload, const char *indexname);
  bool         isIndexValid(const char *indexname=0L)  { return(isIndexValid(FASTA_INDEX_ANY, indexname)); };
  bool         isIndexValid(u32bit indextypetoload)    { return(isIndexValid(indextypetoload, 0L));        };

  //  Opens an existing index, or creates and opens an index if none exists.
  //
  void         openIndex(u32bit indextypetoload, const char *indexname);
  void         openIndex(const char *indexname=0L)  { openIndex(FASTA_INDEX_ANY, indexname); };
  void         openIndex(u32bit indextypetoload)    { openIndex(indextypetoload, 0L);        };

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

  //  Tests on a disk resident index.
  //    Returns true if there is an index on disk and it is compatible
  //    with what we want.
  //
  bool             isIndexOnDiskCompatible(u32bit indextype, const char *indexname, bool beVerbose=true);

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

  bool            _isStreamInput;    //  Forced to be not random access
  bool            _isRandomAccess;   //  Set if an index is open
  bool            _isChecksummed;    //  Set if the index has md5's

  _idxfa_global   _theGlobalDesc;
};

#endif  //  FASTA_H

