
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
 *    src/stores/ovStore.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-OCT-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_OVOVERLAP_H
#define AS_OVOVERLAP_H


//  Error rates are encoded as a 12-bit fixed-point value.  This gives us up to 40.95% error, with
//  0.01% resolution.  Changing the number of bits WILL break the carefully structured
//  ovOverlapDAT.
//
//  The decoded value is a double representing fraction error -- between 0.0000 and 1.0000.
//  The encoded value is an integer type (see the ovsOvelrapDAT below).

#define AS_MAX_EVALUE_BITS         12
#define AS_MAX_EVALUE              ((1 << AS_MAX_EVALUE_BITS) - 1)

#define AS_OVS_decodeEvalue(E)   ((E) / 10000.0)
#define AS_OVS_encodeEvalue(Q)   (((Q) < AS_OVS_decodeEvalue(AS_MAX_EVALUE)) ? (int)(10000.0 * (Q) + 0.5) : AS_MAX_EVALUE)

#define AS_MAX_ERATE               AS_OVS_decodeEvalue(AS_MAX_EVALUE)

//  The old implementation allowed up to 20-bit reads, and used 3 32-bit words.  No alignment was
//  stored.
//
//  The new implementation uses either 3 32-bit words (for EXACTLY 16-bit reads), 2 64-bit words
//  (for up to 21-bit reads) or 6 32-bit words.  It can optionally use 64 bits for storing a pointer
//  to the alignments, though this is not actually implemented.

#define DO_NOT_STORE_ALIGN_PTR


#if AS_MAX_READLEN_BITS < 17

#define        ovOverlapNWORDS  3
#define        ovOverlapWORDSZ  32
typedef uint32 ovOverlapWORD;
#define        F_OV   F_U32
#define        F_OVP  F_U32P

class ovOverlapDAT {
public:
  ovOverlapWORD  ahg5           : AS_MAX_READLEN_BITS;   //  16
  ovOverlapWORD  ahg3           : AS_MAX_READLEN_BITS;   //  16

  ovOverlapWORD  bhg5           : AS_MAX_READLEN_BITS;   //  16
  ovOverlapWORD  bhg3           : AS_MAX_READLEN_BITS;   //  16

  ovOverlapWORD  span           : AS_MAX_READLEN_BITS;   //  16
  ovOverlapWORD  evalue         : AS_MAX_EVALUE_BITS;    //  12
  ovOverlapWORD  flipped        : 1;                     //  1

  ovOverlapWORD  forOBT         : 1;                     //  1
  ovOverlapWORD  forDUP         : 1;                     //  1
  ovOverlapWORD  forUTG         : 1;                     //  1

#ifndef DO_NOT_STORE_ALIGN_PTR
#undef  ovOverlapNWORDS
#define ovOverlapNWORDS 5
  ovOverlapWORD  alignSwapped   : 1;                     //  Our IDs are opposite those in the alignment
  ovOverlapWORD  alignFile      : 19;                    //  Which file of overlap alignments
  ovOverlapWORD  alignPosHi     : 12;                    //  Position in that file (high-order bits)
  ovOverlapWORD  alignPosLo     : 32;                    //  Position in that file (low-order bits)
#endif
};

#elif AS_MAX_READLEN_BITS < 22

#define        ovOverlapNWORDS  2
#define        ovOverlapWORDSZ  64
typedef uint64 ovOverlapWORD;
#define        F_OV   F_U64
#define        F_OVP  F_U64P

class ovOverlapDAT {
public:
  ovOverlapWORD  ahg5           : AS_MAX_READLEN_BITS;   //  17-21
  ovOverlapWORD  ahg3           : AS_MAX_READLEN_BITS;   //  17-21
  ovOverlapWORD  evalue         : AS_MAX_EVALUE_BITS;    //  12
  ovOverlapWORD  flipped        : 1;                     //  1
  ovOverlapWORD  forOBT         : 1;                     //  1
  ovOverlapWORD  forDUP         : 1;                     //  1
  ovOverlapWORD  forUTG         : 1;                     //  1
  ovOverlapWORD  extra1         : 64 - 2 * AS_MAX_READLEN_BITS - AS_MAX_EVALUE_BITS - 1 - 1 - 1 - 1;  //  Between 15 and 7

  ovOverlapWORD  bhg5           : AS_MAX_READLEN_BITS;   //  17-21
  ovOverlapWORD  bhg3           : AS_MAX_READLEN_BITS;   //  17-21
  ovOverlapWORD  span           : AS_MAX_READLEN_BITS;   //  17-21
  ovOverlapWORD  extra2         : 64 - 3 * AS_MAX_READLEN_BITS;  //  Between 13 and 1

#ifndef DO_NOT_STORE_ALIGN_PTR
#undef  ovOverlapNWORDS
#define ovOverlapNWORDS 3
  ovOverlapWORD  alignSwapped   : 1;                     //  Our IDs are opposite those in the alignment
  ovOverlapWORD  alignFile      : 19;                    //  Which file of overlap alignments
  ovOverlapWORD  alignPos       : 44;                    //  Position in that file
#endif
};

#else

#define        ovOverlapNWORDS  6
#define        ovOverlapWORDSZ  32
typedef uint32 ovOverlapWORD;
#define        F_OV   F_U32
#define        F_OVP  F_U32P

class ovOverlapDAT {
public:
  ovOverlapWORD  ahg5;
  ovOverlapWORD  ahg3;
  ovOverlapWORD  bhg5;
  ovOverlapWORD  bhg3;
  ovOverlapWORD  span;

  ovOverlapWORD  evalue         : AS_MAX_EVALUE_BITS;    //  12
  ovOverlapWORD  flipped        : 1;                     //  1
  ovOverlapWORD  forOBT         : 1;                     //  1
  ovOverlapWORD  forDUP         : 1;                     //  1
  ovOverlapWORD  forUTG         : 1;                     //  1
  ovOverlapWORD  extra          : 32 - AS_MAX_EVALUE_BITS - 1 - 1 - 1 - 1;  //  Between 15 and 7

#ifndef DO_NOT_STORE_ALIGN_PTR
#undef  ovOverlapNWORDS
#define ovOverlapNWORDS 8
  ovOverlapWORD  alignSwapped   : 1;                     //  Our IDs are opposite those in the alignment
  ovOverlapWORD  alignFile      : 19;                    //  Which file of overlap alignments
  ovOverlapWORD  alignPosHi     : 12;                    //  Position in that file (high-order bits)
  ovOverlapWORD  alignPosLo     : 32;                    //  Position in that file (low-order bits)
#endif
};

#endif



enum ovOverlapDisplayType {
  ovOverlapAsHangs      = 0,  //  Show a and b hang
  ovOverlapAsCoords     = 1,  //  Show bgn,end for each read
  ovOverlapAsRaw        = 2,  //  Show all four hangs
  ovOverlapAsCompat     = 3,  //  Show in a format more-or-less compatible with CA 8.3
  ovOverlapAsPaf        = 4,  //  Show in a format compatible with miniasm
};




class ovOverlap {
private:
  ovOverlap() {
    g = NULL;
    clear();
  };

public:
  ovOverlap(gkStore *gkp) {
    g = gkp;
    clear();
  };

  ~ovOverlap() {
  };

  static
  ovOverlap  *allocateOverlaps(gkStore *gkp, uint64 num) {
    ovOverlap *r = new ovOverlap [num];

    for (uint32 ii=0; ii<num; ii++)
      r[ii].g = gkp;

    return(r);
  };


  //  Dovetail if any of the following are true:
  //    ahg3 == 0  &&  ahg5 == 0  (a is contained)
  //    ahg3 == 0  &&  bhg5 == 0  (a3' dovetail b5')
  //
  //    bhg3 == 0  &&  bhg5 == 0  (b is contained)
  //    bhg3 == 0  &&  ahg5 == 0  (a5' dovetail b3')
  //
  //  In general, if the 3' hang of either A or B is zero, AND the 5' hang of either A or B is zero.
  //
  uint32     overlapIsDovetail(void) const {
    return(((dat.ovl.ahg5 == 0) || (dat.ovl.bhg5 == 0)) &&
           ((dat.ovl.ahg3 == 0) || (dat.ovl.bhg3 == 0)));
  };


  //  These assume that at most one of ahg5 and bhg5 (or 3') is positive.  If two are positive, then the overlap is partial.
  //
  //  The conversion from a_hang is trivial:
  //    a_hang > 0 ---> ahg5 > 0 (and bhg5 == 0)
  //    a_hang < 0 ---> bhg5 > 0 (and ahg5 == 0)
  //
  //    b_hang > 0 ---> bhg3 > 0 (and ahg3 == 0)
  //    b_hang < 0 ---> ahg3 > 0 (and bhg3 == 0)
  //

  //  Convenience functions.
  int32      a_hang(void) const         { return((int32)dat.ovl.ahg5 - (int32)dat.ovl.bhg5); };
  int32      b_hang(void) const         { return((int32)dat.ovl.bhg3 - (int32)dat.ovl.ahg3); };

  void       a_hang(int32 a)            { dat.ovl.ahg5 = (a < 0) ? 0 : a;  dat.ovl.bhg5 = (a < 0) ? -a : 0; };
  void       b_hang(int32 b)            { dat.ovl.bhg3 = (b < 0) ? 0 : b;  dat.ovl.ahg3 = (b < 0) ? -b : 0; };

  //  These return the actual coordinates on the read.  For reverse B reads, the coordinates are in the reverse-complemented
  //  sequence, and are returned as bgn > end to show this.
  uint32     a_bgn(void) const          { return(dat.ovl.ahg5); };
  uint32     a_end(void) const          { return(g->gkStore_getRead(a_iid)->gkRead_sequenceLength() - dat.ovl.ahg3); };

  uint32     b_bgn(void) const          { return((dat.ovl.flipped) ? (g->gkStore_getRead(b_iid)->gkRead_sequenceLength() - dat.ovl.bhg5) : (dat.ovl.bhg5)); };
  uint32     b_end(void) const          { return((dat.ovl.flipped) ? (dat.ovl.bhg3) : (g->gkStore_getRead(b_iid)->gkRead_sequenceLength() - dat.ovl.bhg3)); };

  uint32     span(void) const           { return(dat.ovl.span); };
  void       span(uint32 s)             { dat.ovl.span = s; };

#if 0
  //  Return an approximate span as the average of the read span aligned.
  uint32     span(void) const {
    if (dat.ovl.span > 0)
      return(dat.ovl.span);
    else {
      uint32 ab = a_bgn(), ae = a_end();
      uint32 bb = b_bgn(), be = b_end();

      if (bb < be)
        return(((ae - ab) + (be - bb)) / 2);
      else
        return(((ae - ab) + (bb - be)) / 2);
    }
  }
#endif

  void       flipped(uint32 f)          { dat.ovl.flipped = f; };
  uint32     flipped(void) const        { return(dat.ovl.flipped == true); };

  void       erate(double e)            { dat.ovl.evalue = AS_OVS_encodeEvalue(e); };
  double     erate(void) const          { return(AS_OVS_decodeEvalue(dat.ovl.evalue)); };

  void       evalue(uint64 e)           { dat.ovl.evalue = e; };
  uint64     evalue(void) const         { return(dat.ovl.evalue); };

  bool       forOBT(void)               { return(dat.ovl.forOBT); };
  bool       forDUP(void)               { return(dat.ovl.forDUP); };
  bool       forUTG(void)               { return(dat.ovl.forUTG); };

  //  These are true only if the overlap is dovetail, which is the usual case, and isn't checked.

  uint32     overlapAEndIs5prime(void) const  { return((dat.ovl.bhg5 > 0) && (dat.ovl.ahg3  > 0));  };
  uint32     overlapAEndIs3prime(void) const  { return((dat.ovl.ahg5 > 0) && (dat.ovl.bhg3  > 0));  };

  uint32     overlapBEndIs5prime(void) const  { return((overlapAEndIs5prime() && (dat.ovl.flipped == true)) ||
                                                       (overlapAEndIs3prime() && (dat.ovl.flipped == false))); };
  uint32     overlapBEndIs3prime(void) const  { return((overlapAEndIs5prime() && (dat.ovl.flipped == false)) ||
                                                       (overlapAEndIs3prime() && (dat.ovl.flipped == true))); };

  uint32     overlapAIsContained(void) const  { return((dat.ovl.ahg5 == 0) && (dat.ovl.ahg3 == 0));  };
  uint32     overlapBIsContainer(void) const  { return((dat.ovl.ahg5 == 0) && (dat.ovl.ahg3 == 0));  };

  uint32     overlapAIsContainer(void) const  { return((dat.ovl.bhg5 == 0) && (dat.ovl.bhg3 == 0));  };
  uint32     overlapBIsContained(void) const  { return((dat.ovl.bhg5 == 0) && (dat.ovl.bhg3 == 0));  };

  //  Test if the overlap is dovetail or partial.

  uint32     overlap5primeIsPartial(void) const { return((dat.ovl.ahg5 > 0) && (dat.ovl.bhg5 > 0)); };
  uint32     overlap3primeIsPartial(void) const { return((dat.ovl.ahg3 > 0) && (dat.ovl.bhg3 > 0)); };

  uint32     overlapIsPartial(void)       const { return(overlap5primeIsPartial() || overlap3primeIsPartial()); };

  char      *toString(char *str, ovOverlapDisplayType type, bool newLine);

  void       swapIDs(ovOverlap const &orig);

  void       clear(void) {
    //g        = NULL;    //  Explicitly DO NOT clear the pointer to gkpStore.

    for (uint32 ii=0; ii<ovOverlapNWORDS; ii++)
      dat.dat[ii] = 0;

    a_iid      = 0;
    b_iid      = 0;
  };

  bool
  operator<(const ovOverlap &that) const {
    if (a_iid      < that.a_iid)       return(true);
    if (a_iid      > that.a_iid)       return(false);
    if (b_iid      < that.b_iid)       return(true);
    if (b_iid      > that.b_iid)       return(false);

    for (uint32 ii=0; ii<ovOverlapNWORDS; ii++) {
      if (dat.dat[ii] < that.dat.dat[ii])  return(true);
      if (dat.dat[ii] > that.dat.dat[ii])  return(false);
    }

    return(false);
  };

public:
  gkStore             *g;

public:
  uint32               a_iid;
  uint32               b_iid;

  union {
    ovOverlapWORD     dat[ovOverlapNWORDS];
    ovOverlapDAT      ovl;
  } dat;
};


#endif  //  AS_OVOVERLAP_H
