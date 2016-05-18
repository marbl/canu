
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-DEC-09 to 2015-JUL-01
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-MAR-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_OVERLAP_H
#define AS_OVERLAP_H

#include "AS_global.H"
#include "gkStore.H"

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
//  The new implementation uses either 5 (for EXACTLY 16-bit reads) or 6 32-bit words.  It uses 2
//  32-bit words for storing a pointer to the alignments.
//
//  Note that the 5-word version must use uint32, and the 6-word version must use uint64.
//  Note that the 5-word version needs to split out the alignPos into two words.

#if AS_MAX_READLEN_BITS < 17

#define        ovOverlapNWORDS  5
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

  ovOverlapWORD  alignSwapped   : 1;                     //  Our IDs are opposite those in the alignment
  ovOverlapWORD  alignFile      : 19;                    //  Which file of overlap alignments
  ovOverlapWORD  alignPosHi     : 12;                    //  Position in that file (high-order bits)
  ovOverlapWORD  alignPosLo     : 32;                    //  Position in that file (low-order bits)
};

#elif AS_MAX_READLEN_BITS < 22

#define        ovOverlapNWORDS  3
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

  ovOverlapWORD  alignSwapped   : 1;                     //  Our IDs are opposite those in the alignment
  ovOverlapWORD  alignFile      : 19;                    //  Which file of overlap alignments
  ovOverlapWORD  alignPos       : 44;                    //  Position in that file
};

#else

#define        ovOverlapNWORDS  8
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

  ovOverlapWORD  alignSwapped   : 1;                     //  Our IDs are opposite those in the alignment
  ovOverlapWORD  alignFile      : 19;                    //  Which file of overlap alignments
  ovOverlapWORD  alignPosHi     : 12;                    //  Position in that file (high-order bits)
  ovOverlapWORD  alignPosLo     : 32;                    //  Position in that file (low-order bits)
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
    dat.dat[0] = 0;
    dat.dat[1] = 0;
    dat.dat[2] = 0;
#if (ovOverlapNWORDS == 5)
    dat.dat[3] = 0;
    dat.dat[4] = 0;
#endif
  };

  bool
  operator<(const ovOverlap &that) const {
    if (a_iid      < that.a_iid)       return(true);
    if (a_iid      > that.a_iid)       return(false);
    if (b_iid      < that.b_iid)       return(true);
    if (b_iid      > that.b_iid)       return(false);
    if (dat.dat[0] < that.dat.dat[0])  return(true);
    if (dat.dat[0] > that.dat.dat[0])  return(false);
    if (dat.dat[1] < that.dat.dat[1])  return(true);
    if (dat.dat[1] > that.dat.dat[1])  return(false);
    if (dat.dat[2] < that.dat.dat[2])  return(true);
    if (dat.dat[2] > that.dat.dat[2])  return(false);
#if (ovOverlapNWORDS == 5)
    if (dat.dat[3] < that.dat.dat[3])  return(true);
    if (dat.dat[3] > that.dat.dat[3])  return(false);
    if (dat.dat[4] < that.dat.dat[4])  return(true);
    if (dat.dat[4] > that.dat.dat[4])  return(false);
#endif
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


//  The default, no flags, is to open for normal overlaps, read only.  Normal overlaps mean they
//  have only the B id, i.e., they are in a fully built store.
//
//  Output of overlapper (input to store building) should be ovFileFullWrite.  The specialized
//  ovFileFullWriteNoCounts is used internally by store creation.
//
enum ovFileType {
  ovFileNormal              = 0,  //  Reading of b_id overlaps (aka store files)
  ovFileNormalWrite         = 1,  //  Writing of b_id overlaps
  ovFileFull                = 2,  //  Reading of a_id+b_id overlaps (aka dump files)
  ovFileFullWrite           = 3,  //  Writing of a_id+b_id overlaps
  ovFileFullWriteNoCounts   = 4   //  Writing of a_id+b_id overlaps, omitting the counts of olaps per read
};




class ovFile {
public:
  ovFile(const char  *name,
         ovFileType   type = ovFileNormal,
         uint32       bufferSize = 1 * 1024 * 1024);
  ~ovFile();

  void    flushOverlaps(void);

  void    writeOverlap(ovOverlap *overlap);
  void    writeOverlaps(ovOverlap *overlaps, uint64 overlapLen);

  bool    readOverlap(ovOverlap *overlap);
  uint64  readOverlaps(ovOverlap *overlaps, uint64 overlapMax);

  void    seekOverlap(off_t overlap);

  //  The size of an overlap record is 1 or 2 IDs + the size of a word times the number of words.
  uint64  recordSize(void) {
    return(sizeof(uint32) * ((_isNormal) ? 1 : 2) + sizeof(ovOverlapWORD) * ovOverlapNWORDS);
  };

private:
  uint32                  _bufferLen;    //  length of valid data in the buffer
  uint32                  _bufferPos;    //  position the read is at in the buffer
  uint32                  _bufferMax;    //  allocated size of the buffer
  uint32                 *_buffer;

  uint32                  _olapsPerReadAlloc;
  uint32                  _olapsPerReadLast;
  uint32                 *_olapsPerRead;

  bool                    _isOutput;     //  if true, we can writeOverlap()
  bool                    _isSeekable;   //  if true, we can seekOverlap()
  bool                    _isNormal;     //  if true, 3 words per overlap, else 4

  compressedFileReader   *_reader;
  compressedFileWriter   *_writer;

  char                    _prefix[FILENAME_MAX];
  FILE                   *_file;
};





class ovStoreInfo {
private:
  uint64    _ovsMagic;
  uint64    _ovsVersion;
  uint64    _UNUSED;
  uint64    _smallestIID;         //  smallest frag iid in the store
  uint64    _largestIID;          //  largest frag iid in the store
  uint64    _numOverlapsTotal;    //  number of overlaps in the store
  uint64    _highestFileIndex;
  uint64    _maxReadLenInBits;    //  length of a fragment

  friend class ovStore;

  friend
  void       writeOverlaps(char       *storePath,
                           ovOverlap *ovls,
                           uint64      ovlsLen,
                           uint32      fileID);

  friend
  bool
  testIndex(char     *storePath,
            bool      doFixes);

  friend
  void
  mergeInfoFiles(char       *storePath,
                 uint32      nPieces);
};


class ovStoreOfft {
private:
  uint32    _a_iid;      //  read ID for this block of overlaps.

  uint32    _fileno;     //  the file that contains this a_iid
  uint32    _offset;     //  offset to the first overlap for this iid
  uint32    _numOlaps;   //  number of overlaps for this iid

  uint64    _overlapID;  //  overlapID for the first overlap in this block.  in memory, this is the id of the next overlap.

  void       clear(void) {
    _a_iid     = 0;
    _fileno    = 0;
    _offset    = 0;
    _numOlaps  = 0;
    _overlapID = 0;
  };

  friend class ovStore;

  friend
  void
  writeOverlaps(char       *storePath,
                ovOverlap *ovls,
                uint64      ovlsLen,
                uint32      fileID);

  friend
  bool
  testIndex(char     *storePath,
            bool      doFixes);

  friend
  void
  mergeInfoFiles(char       *storePath,
                 uint32      nPieces);
};


//  The default here is to open a read only store.
//
enum ovStoreType {
  ovStoreReadOnly  = 0,
  ovStoreWrite     = 1,  //  Open for write, fail if one exists already
  ovStoreOverwrite = 2,  //  Open for write, and obliterate an existing store
};


class ovStore {
private:
  void       ovStore_read(void);
  void       ovStore_write(void);

public:
  ovStore(const char *name, gkStore *gkp, ovStoreType cType=ovStoreReadOnly);
  ~ovStore();

  //  Read the next overlap from the store.  Return value is the number of overlaps read.
  uint32     readOverlap(ovOverlap *overlap);

  //  Return the number of overlaps that would be read.  Basically the same as the next readOverlaps() call.
  uint32     numberOfOverlaps(void);

  //  Read ALL remaining overlaps for the current A_iid.  Return value is the number of overlaps read.
  uint32     readOverlaps(ovOverlap *&overlaps,
                          uint32     &maxOverlaps,
                          bool        restrictToIID=true);

  //  Append ALL remaining overlaps for the current A_iid to the overlaps in ovl.  Return value is
  //  the number of overlaps in ovl that are for A_iid == iid.
  //
  //  It is up to the client to verify that ovl[0] is the same as iid (e.g., that the return value
  //  is not zero); ovlLen is the number of overlaps in ovl, NOT the number of overlaps in ovl that
  //  are the same as iid.
  //
  uint32       readOverlaps(uint32       iid,
                            ovOverlap  *&ovl,
                            uint32      &ovlLen,
                            uint32      &ovlMax);

  void         setRange(uint32 low, uint32 high);
  void         resetRange(void);

  uint64       numOverlapsInRange(void);
  uint32 *     numOverlapsPerFrag(uint32 &firstFrag, uint32 &lastFrag);

  //  The (mostly) private interface for adding overlaps to a store.  Overlaps must be sorted already.

  void         writeOverlap(ovOverlap *olap);
  void         writeOverlap(ovOverlap *overlap, uint32 maxOverlapsThisFile);

  //  Write a block of sorted overlaps to store file 'fileID', saving the info and index into
  //  'fileID.info' and 'fileID.index'

  friend
  void         writeOverlaps(char       *storePath,
                             ovOverlap *ovls,
                             uint64      ovlsLen,
                             uint32      fileID);


  //  Add new evalues for reads between bgnID and endID.  No checking of IDs is done, but the number
  //  of evalues must agree.

  void       addEvalues(uint32 bgnID, uint32 endID, uint16 *evalues, uint64 evaluesLen);

private:
  char               _storePath[FILENAME_MAX];

  bool               _isOutput;

  ovStoreInfo        _info;

  uint32             _firstIIDrequested;
  uint32             _lastIIDrequested;

  FILE              *_offtFile;   //  For writing overlaps, a place to dump ovStoreOfft's.
  ovStoreOfft        _offt;       //  For writing overlaps, the current ovStoreOfft.
  ovStoreOfft        _offm;       //  For writing overlaps, an empty ovStoreOfft, for reads with no overlaps.

  memoryMappedFile  *_evaluesMap;
  uint16            *_evalues;

  uint64             _overlapsThisFile;  //  Count of the number of overlaps written so far
  uint32             _currentFileIndex;
  ovFile            *_bof;

  gkStore           *_gkp;
};


//  This should be part of ovStore, but when it is used, in ovStoreSorter, we don't
//  have a store opened.
void
writeOverlaps(char       *storePath,
              ovOverlap *ovls,
              uint64      ovlsLen,
              uint32      fileID);

bool
testIndex(char     *storePath,
          bool      doFixes);

void
mergeInfoFiles(char       *storePath,
               uint32      nPieces);







//  For store construction.  Probably should be in either ovOverlap or ovStore.

class ovStoreFilter {
public:
  ovStoreFilter(gkStore *gkp_, double maxErate) {
    gkp             = gkp_;

    resetCounters();

    maxID     = gkp->gkStore_getNumReads() + 1;
    maxEvalue = AS_OVS_encodeEvalue(maxErate);

    skipReadOBT     = new char [maxID];
    skipReadDUP     = new char [maxID];

    memset(skipReadOBT, 0, sizeof(char) * maxID);
    memset(skipReadDUP, 0, sizeof(char) * maxID);


    uint32  numSkipOBT = 0;
    uint32  numSkipDUP = 0;

    fprintf(stderr, "Marking fragments to skip overlap based trimming.\n");

    fprintf(stderr, "LIB 1 - dup=%d trim=%d spur=%d chimera=%d subreads=%d\n",
            gkp->gkStore_getLibrary(1)->gkLibrary_removeDuplicateReads(),
            gkp->gkStore_getLibrary(1)->gkLibrary_finalTrim(),
            gkp->gkStore_getLibrary(1)->gkLibrary_removeSpurReads(),
            gkp->gkStore_getLibrary(1)->gkLibrary_removeChimericReads(),
            gkp->gkStore_getLibrary(1)->gkLibrary_checkForSubReads());

    for (uint64 iid=0; iid<maxID; iid++) {
      uint32     Lid = gkp->gkStore_getRead(iid)->gkRead_libraryID();
      gkLibrary *L   = gkp->gkStore_getLibrary(Lid);

      if ((L->gkLibrary_removeDuplicateReads()     == false) &&
          (L->gkLibrary_finalTrim()                == GK_FINALTRIM_NONE) &&
          (L->gkLibrary_removeSpurReads()          == false) &&
          (L->gkLibrary_removeChimericReads()      == false) &&
          (L->gkLibrary_checkForSubReads()         == false)) {
        numSkipOBT++;
        skipReadOBT[iid] = true;
      }

      if (L->gkLibrary_removeDuplicateReads() == false) {
        numSkipDUP++;
        skipReadDUP[iid] = true;
      }
    }

    fprintf(stderr, "Marked "F_U32" reads so skip OBT, "F_U32" reads to skip dedupe.\n", numSkipOBT, numSkipDUP);
  };

  ~ovStoreFilter() {
    delete [] skipReadOBT;
    delete [] skipReadDUP;
  };


  void    filterOverlap(ovOverlap     &foverlap,
                        ovOverlap     &roverlap);

  void    reportFate(void);
  void    resetCounters(void);

public:
  gkStore *gkp;

  uint32   maxID;
  uint32   maxEvalue;

  uint64   saveUTG;
  uint64   saveOBT;
  uint64   saveDUP;

  uint64   skipERATE;

  uint64   skipOBT;        //  OBT not requested for the A read
  uint64   skipOBTbad;     //  Overlap too similiar
  uint64   skipOBTshort;   //  Overlap is too short

  uint64   skipDUP;        //  DUP not requested for the A read
  uint64   skipDUPdiff;    //  Overlap isn't remotely similar
  uint64   skipDUPlib;

  //  Not really stats, but global state for the filter.

  char    *skipReadOBT;
  char    *skipReadDUP;
};



#endif
