
//  Leave the sequence mostly on disk.  A small buffer is used to
//  cache reads.  The defline is read into core.
//
//  Two operations are currently supported:
//    char  getChar(void)
//    char *getChars(char *block, u32bit position, u32bit length)
//
//  When the end of the sequence is reached, getChar() returns 0.

class seqOnDisk {
public:
  seqOnDisk(char const *filename,
            u64bit hdrstart, u32bit hdrlen,
            u64bit seqstart, u32bit seqlen,
            seqIID iid,
            bool isSqueezed,
            bool isFixedWidth,
            u32bit sl,
            u32bit ss);
  seqOnDisk(seqIID iid,
            char *hdr, u32bit hdrlen,
            char *seq, u32bit seqlen);
  ~seqOnDisk();

  char           *header(void) const         { return(_header);       };
  u32bit          headerLength(void) const   { return(_headerLength); };

  char           *sequence(void) const       { return(_sequence);       };
  u32bit          sequenceLength(void) const { return(_sequenceLength); };

  u32bit          getIID(void) const         { return(_idx); };

  //  next() returns true on eof.
  //  get() returns 0 on eof.
  //
  bool            next(void);
  char            get(void);
  char           *getChars(char *block, u32bit position, u32bit length);
private:

  readBuffer    *_readBuffer;

  u32bit         _idx;

  u32bit         _headerLength;
  u32bit         _sequenceLength;
  u64bit         _headerStart;
  u64bit         _sequenceStart;

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
  u32bit         _sequencePosition;

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
  u32bit         _lineLength;
  u32bit         _lineSep;
};


//  Advances to the next non-space letter.
//  Returns false if there is no next non-space letter.
//
inline
bool
seqOnDisk::next(void) {

  //  Check if we're already at the end
  //
  if (_sequencePosition >= _sequenceLength)
    return(false);

  //  Nope, we know there is another non-space letter out there
  //  somewhere.

  if (_readBuffer) {
    _readBuffer->next();

    if (_sourceType != 0)
      while (whitespaceSymbol[_readBuffer->get()])
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
seqOnDisk::get(void) {

  if (_sequencePosition >= _sequenceLength)
    return(0);

  if (_readBuffer)
    return(_readBuffer->get());
  else
    return(_sequence[_sequencePosition]);
}
