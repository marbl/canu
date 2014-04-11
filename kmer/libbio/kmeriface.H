
#if 0

//  Documentation, really.

//  Incomplete too.

class kMerInterface {
  kMerInterface() {};
  virtual ~kMerInterface() {};

  //  Reverse all the words, reverse and complement the bases in
  //  each word, then shift right to align the edge.
  //
  virtual kMerInterface  &reverseComplement(void) = 0;
  virtual void            clear(void);

  //  Construct a mer by shifting bases onto the end:
  //    += shifts onto the right end
  //    -= shifts onto the left end
  //
  virtual void   operator+=(uint64 x) = 0;
  virtual void   operator-=(uint64 x) = 0;

  //  used by merStream at least
  //
  virtual void     mask(bool) = 0;

  //  Return the mer, as a 64-bit integer.  If the mer is more than
  //  32-bases long, then the left-most (the earliest, the start, etc)
  //  bases are used.
  //
  virtual operator uint64 () const = 0;

  //  These are written/read in 5'endian, which isn't the most natural
  //  implementation.  It's done this way to keep the sequence in
  //  order (e.g., the merStreamFile).  Don't change the order.
  //
  //  On the otherhand, the implementation (of write anyway) is
  //  basically the same as merToString().
  //
  //  Takes an optional number of BITS to write, pulled from the
  //  END of the mer.
  //
  virtual void   writeToBitPackedFile(bitPackedFile *BPF, uint32 numBits=0) const = 0;
  virtual void   readFromBitPackedFile(bitPackedFile *BPF, uint32 numBits=0) = 0;

  //  Returns a sub-mer from either the start (left end) or the end
  //  (right end) of the mer.  The sub-mer must be at most 64 bits
  //  long.  Yes, BITS.
  //
  //  The start is difficult, because it can span multiple words.  The
  //  end is always in the first word.
  //
  virtual uint64   startOfMer(uint32 bits) const = 0;
  virtual uint64   endOfMer(uint32 bits) const = 0;

  //  Set 'numbits' bits from (the end of) 'val' at bit position 'pos'
  //  in the mer.  This is wildly low-level, but merylStreamReader
  //  needs it.
  //
  //  The position is measured from the right end.
  //    (0, 8, X) would copy the bits 7 to 0 of X to bits 7 to 0 of the mer.
  //
  //  Argh!  Can't use set/getDecodedValue because that is doing things in the wrong order.
  //
  //  Meryl
  //
  virtual uint64   getWord(uint32 wrd) const        = 0;  // { return(MERWORD(wrd)); };
  virtual void     setWord(uint32 wrd, uint64 val)  = 0;  // { MERWORD(wrd) = val;   };

  //  Show the mer as ascii
  //
  //  Doesn't print the last full word, if it's on the word boundary
  //
  //  We build the string right to left, print any partial word first,
  //  then print whole words until we run out of words to print.
  //
  virtual char    *merToString(char *instr) const = 0;
};


#endif
