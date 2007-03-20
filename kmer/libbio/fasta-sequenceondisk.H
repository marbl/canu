#ifndef FASTA_SEQUENCEONDISK_H
#define FASTA_SEQUENCEONDISK_H

#include "bio.h"
#include "fasta-types.H"


//  Leave the sequence mostly on disk.  A small buffer is used to
//  cache reads.  The defline is read into core.
//
//  Two operations are currently supported:
//    char  getChar(void)
//    char *getChars(char *block, u32bit position, u32bit length)
//
//  When the end of the sequence is reached, getChar() returns 0.
//
//  Only supported for indexed and squeezed fasta file sources
//


//
//  Three types of sequence
//   -- fully squeezed: use a readBuffer directly, both stream and block copy
//   -- fixed width: stream and block copy
//   -- free form: only stream
//



class FastASequenceOnDisk {
public:
  FastASequenceOnDisk(char const *filename,
                      _idxfa_pos hdrstart, _idxfa_len hdrlen,
                      _idxfa_pos seqstart, _idxfa_len seqlen,
                      IID_t iid,
                      bool isSqueezed,
                      bool isFixedWidth,
                      _idxfa_len sl,
                      _idxfa_len ss);
  FastASequenceOnDisk(IID_t iid,
                      char *hdr, u32bit hdrlen,
                      char *seq, u32bit seqlen);
  ~FastASequenceOnDisk();

  char           *header(void) const         { return(_header);       };
  u32bit          headerLength(void) const   { return(_headerLength); };

  char           *sequence(void) const       { return(_sequence);       };
  u32bit          sequenceLength(void) const { return(_sequenceLength); };

  u32bit          getIID(void) const         { return(_idx); };

  //  XXX:  Probably fails if the sequence ends with whitespace
  //#warning FastASequenceOnDisk might fail if the sequence ends with whitespace!

  //  next() returns true on eof.
  //  get() returns 0 on eof.
  //
  bool            next(void);
  char            get(void);
  char           *getChars(char *block, u32bit position, u32bit length);
private:

  readBuffer    *_readBuffer;

  //  The iid of this guy
  //
  u32bit         _idx;

  //  Information about the sequence on disk
  //
  _idxfa_len     _headerLength;
  _idxfa_len     _sequenceLength;
  _idxfa_pos     _headerStart;
  _idxfa_pos     _sequenceStart;

  //  A copy of the defline.  We get it for nearly free when we make
  //  the first read buffer.
  //
  char          *_header;

  //  If defined, this means we are not readBuffer backed, but
  //  we are faking it.
  //
  char          *_sequence;

  //  The position we are at in the sequence.  Used to decide when we
  //  have seen all the sequence.  If the source is not squeezed, this
  //  is not the same as the file position.
  //
  _idxfa_len     _sequencePosition;

  //  Cache the type of the source
  //    0 - squeezed
  //    1 - fixed width
  //    2 - freeform
  //
  u32bit         _sourceType;

  //  If the source file is not squeezed, we'll do on-the-fly
  //  squeezing.  If the source file is fixed width, we know the
  //  length of each line, and the number of bytes between lines.
  //
  _idxfa_len     _lineLength;
  _idxfa_len     _lineSep;
};



//  Advances to the next non-space letter.
//  Returns false if there is no next non-space letter.
//
inline
bool
FastASequenceOnDisk::next(void) {

  //  Check if we're already at the end
  //
  if (_sequencePosition >= _sequenceLength)
    return(false);

  //  Nope, we know there is another non-space letter out there
  //  somewhere.

  if (_readBuffer) {
    _readBuffer->next();

    if (_sourceType != 0)
      while (whitespaceSymbol[(int)_readBuffer->get()])
        _readBuffer->next();    
  }

  _sequencePosition++;

  return(true);
}



//  Returns 0 on eof (or end of sequence).  Does not advance the
//  file position.
//
inline
char
FastASequenceOnDisk::get(void) {

  if (_sequencePosition >= _sequenceLength)
    return(0);

  if (_readBuffer)
    return(_readBuffer->get());
  else
    return(_sequence[_sequencePosition]);
}



#endif  //  FASTA_SEQUENCEONDISK_H

