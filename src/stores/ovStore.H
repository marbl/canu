

#ifndef AS_OVERLAP_H
#define AS_OVERLAP_H

static const char *rcsid_ovStore_H = "$Id:  $";

#include "AS_global.H"
#include "gkStore.H"

//  Error rates are encoded as a 12-bit fixed-point value.  This gives us up to 40.95% error, with
//  0.01% resolution.  Changing the number of bits WILL break the carefully structured
//  ovsOverlapDAT.
//
//  The decoded value is a double representing fraction error -- between 0.0000 and 1.0000.
//  The encoded value is an integer type (see the ovsOvelrapDAT below).

#define AS_MAX_ERATE_BITS         12
#define AS_MAX_ERATE              ((1 << AS_MAX_ERATE_BITS) - 1)

#define AS_OVS_decodeQuality(E)   ((E) / 10000.0)
#define AS_OVS_encodeQuality(Q)   (((Q) < AS_OVS_decodeQuality(AS_MAX_ERATE)) ? (int)(10000.0 * (Q) + 0.5) : AS_MAX_ERATE)


//  The old implementation allowed up to 20-bit reads, and used 3 32-bit words.  No alignment was
//  stored.
//
//  The new implementation uses either 5 (for EXACTLY 16-bit reads) or 6 32-bit words.  It uses 2
//  32-bit words for storing a pointer to the alignments.
//
//  Note that the 5-word version must use uint32, and the 6-word version must use uint64.
//  Note that the 5-word version needs to split out the alignPos into two words.

#if AS_MAX_READLEN_BITS < 17

#define        ovsOverlapNWORDS  5
typedef uint32 ovsOverlapWORD;
#define        F_OV   F_U32
#define        F_OVP  F_U32P

class ovsOverlapDAT {
public:
  ovsOverlapWORD  ahg5           : AS_MAX_READLEN_BITS;   //  16
  ovsOverlapWORD  ahg3           : AS_MAX_READLEN_BITS;   //  16

  ovsOverlapWORD  bhg5           : AS_MAX_READLEN_BITS;   //  16
  ovsOverlapWORD  bhg3           : AS_MAX_READLEN_BITS;   //  16

  ovsOverlapWORD  span           : AS_MAX_READLEN_BITS;   //  16
  ovsOverlapWORD  erate          : AS_MAX_ERATE_BITS;     //  12
  ovsOverlapWORD  flipped        : 1;                     //  1

  ovsOverlapWORD  forOBT         : 1;                     //  1
  ovsOverlapWORD  forDUP         : 1;                     //  1
  ovsOverlapWORD  forUTG         : 1;                     //  1

  ovsOverlapWORD  alignSwapped   : 1;                     //  Our IDs are opposite those in the alignment
  ovsOverlapWORD  alignFile      : 19;                    //  Which file of overlap alignments
  ovsOverlapWORD  alignPosHi     : 12;                    //  Position in that file (high-order bits)
  ovsOverlapWORD  alignPosLo     : 32;                    //  Position in that file (low-order bits)

#if 0
  void  clear(double quality) {
    dat.ovl.ahg5         = 0;
    dat.ovl.ahg3         = 0;
    dat.ovl.bhg5         = 0;
    dat.ovl.bhg3         = 0;
    dat.ovl.span         = 0;
    dat.ovl.erate        = AS_OVS_encodeQuality(quality);
    dat.ovl.flipped      = false;
    dat.ovl.forOBT       = false;
    dat.ovl.forDUP       = false;
    dat.ovl.forUTG       = false;
    dat.ovl.alignSwapped = false;
    dat.ovl.alignFile    = 0;
    dat.ovl.alignPosHi   = 0;
    dat.ovl.alignPosLo   = 0;
  };
#endif
};

#else

#define        ovsOverlapNWORDS  3
typedef uint64 ovsOverlapWORD;
#define        F_OV   F_U64
#define        F_OVP  F_U64P

class ovsOverlapDAT {
puublic:
  ovsOverlapWORD  ahg5           : AS_MAX_READLEN_BITS;   //  17-21
  ovsOverlapWORD  ahg3           : AS_MAX_READLEN_BITS;   //  17-21
  ovsOverlapWORD  erate          : AS_MAX_ERATE_BITS;     //  12
  ovsOverlapWORD  flipped        : 1;                     //  1
  ovsOverlapWORD  forOBT         : 1;                     //  1
  ovsOverlapWORD  forDUP         : 1;                     //  1
  ovsOverlapWORD  forUTG         : 1;                     //  1
  ovsOverlapWORD  extra1         : 64 - 2 * AS_MAX_READLEN_BITS - AS_MAX_ERATE_BITS - 1 - 1 - 1 - 1;  //  Between 15 and 7

  ovsOverlapWORD  bhg5           : AS_MAX_READLEN_BITS;   //  17-21
  ovsOverlapWORD  bhg3           : AS_MAX_READLEN_BITS;   //  17-21
  ovsOverlapWORD  span           : AS_MAX_READLEN_BITS;   //  17-21
  ovsOverlapWORD  extra2         : 64 - 3 * AS_MAX_READLEN_BITS;                          //  Between 13 and 1

  ovsOverlapWORD  alignSwapped   : 1;                     //  Our IDs are opposite those in the alignment
  ovsOverlapWORD  alignFile      : 19;                    //  Which file of overlap alignments
  ovsOverlapWORD  alignPos       : 44;                    //  Position in that file

#if 0
  void  clear(double quality) {
    dat.ovl.ahg5         = 0;
    dat.ovl.ahg3         = 0;
    dat.ovl.bhg5         = 0;
    dat.ovl.bhg3         = 0;
    dat.ovl.span         = 0;
    dat.ovl.erate        = AS_OVS_encodeQuality(quality);
    dat.ovl.flipped      = false;
    dat.ovl.forOBT       = false;
    dat.ovl.forDUP       = false;
    dat.ovl.forUTG       = false;
    dat.ovl.extra1       = 0;
    dat.ovl.extra2       = 0;
    dat.ovl.alignSwapped = false;
    dat.ovl.alignFile    = 0;
    dat.ovl.alignPosHi   = 0;
    dat.ovl.alignPosLo   = 0;
  };
#endif
};

#endif



enum ovOverlapDisplayType {
  ovOverlapAsHangs      = 0,  //  Show a and b hang
  ovOverlapAsCoords     = 1,  //  Show bgn,end for each read
  ovOverlapAsRaw        = 2,  //  Show all four hangs
  ovOverlapAsCompat     = 3,  //  Show in a format more-or-less compatible with CA 8.3
};




class ovsOverlap {
public:
  ovsOverlap() {
    clear();
  };
  ~ovsOverlap() {
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
  uint32     a_bgn(gkStore *g=NULL) const          { return(dat.ovl.ahg5); };
  uint32     a_end(gkStore *g) const    { return(g->gkStore_getRead(a_iid)->gkRead_sequenceLength() - dat.ovl.ahg3); };

  uint32     b_bgn(gkStore *g) const    { return((dat.ovl.flipped) ? (g->gkStore_getRead(b_iid)->gkRead_sequenceLength() - dat.ovl.bhg5) : (dat.ovl.bhg5)); };
  uint32     b_end(gkStore *g) const    { return((dat.ovl.flipped) ? (dat.ovl.bhg3) : (g->gkStore_getRead(b_iid)->gkRead_sequenceLength() - dat.ovl.bhg3)); };

  uint32     span(void) const           { return(dat.ovl.span); };
  void       span(uint32 s)             { dat.ovl.span = s; };

  void       flipped(uint32 f)          { dat.ovl.flipped = f; };
  uint32     flipped(void) const        { return(dat.ovl.flipped == true); };

  void       erate(double e)            { dat.ovl.erate = AS_OVS_encodeQuality(e); };
  double     erate(void) const          { return(AS_OVS_decodeQuality(dat.ovl.erate)); };

  void       evalue(uint64 e)           { dat.ovl.erate = e; };
  uint64     evalue(void) const         { return(dat.ovl.erate); };

  bool       forOBT(void)               { return(dat.ovl.forOBT); };
  bool       forDUP(void)               { return(dat.ovl.forDUP); };
  bool       forUTG(void)               { return(dat.ovl.forUTG); };


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

  char      *toString(char *str, gkStore *gkp, ovOverlapDisplayType type);

  void       swapIDs(ovsOverlap const &orig);

  void       clear(void) {
    dat.dat[0] = 0;
    dat.dat[1] = 0;
    dat.dat[2] = 0;
#if (ovsOverlapNWORDS == 5)
    dat.dat[3] = 0;
    dat.dat[4] = 0;
#endif
  };

  bool
  operator<(const ovsOverlap &that) const {
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
#if (ovsOverlapNWORDS == 5)
    if (dat.dat[3] < that.dat.dat[3])  return(true);
    if (dat.dat[3] > that.dat.dat[3])  return(false);
    if (dat.dat[4] < that.dat.dat[4])  return(true);
    if (dat.dat[4] > that.dat.dat[4])  return(false);
#endif
    return(false);
  };

private:
public:
  uint32               a_iid;
  uint32               b_iid;

  union {
    ovsOverlapWORD     dat[ovsOverlapNWORDS];
    ovsOverlapDAT      ovl;
  } dat;
};


//  If ovsOverlapWORD is uint32, we have no problems.

//  If ovsOverlapWORD is uint64, and we are internal, then we have a problem with the internal format.
//  We want to store just one uint32 id and 3 uint64's, but our buffer wants to be uint64.



//  The flags are a bit goofy.
//
//  The default, no flags, is to open for normal overlaps, read only.  Normal overlaps mean they
//  have only the B id, i.e., they are in a fully built store.
//
//  To write those overlaps, open with ovFileNormal | ovFileWrite.
//
//  To read/write intermediate overlaps (both ids are present), use ovFileFull and ovFileFull | ovFileWrite
//
//  Stupid enum, can't be combined - ofFileNormal | ofFileWrite gives an error.
//
enum ovFileType {
  ovFileNormal      = 0,  //  Reading of b_id overlaps (aka store files)
  ovFileNormalWrite = 1,  //  Writing of b_id overlaps
  ovFileFull        = 2,  //  Reading of a_id+b_id overlaps (aka dump files)
  ovFileFullWrite   = 3   //  Writing of a_id+b_id overlaps
};




class ovFile {
public:
  ovFile(const char  *name,
         ovFileType   type = ovFileNormal,
         uint32       bufferSize = 1 * 1024 * 1024);
  ~ovFile();

  void    flushOverlaps(void);

  void    writeOverlap(ovsOverlap *overlap);
  void    writeOverlaps(ovsOverlap *overlaps, uint64 overlapLen);

  bool    readOverlap(ovsOverlap *overlap);
  uint64  readOverlaps(ovsOverlap *overlaps, uint64 overlapMax);

  void    seekOverlap(off_t overlap);

  //  The size of an overlap record is 1 or 2 IDs + the size of a word times the number of words.
  uint64  recordSize(void) {
    return(sizeof(uint32) * ((_isNormal) ? 1 : 2) + sizeof(ovsOverlapWORD) * ovsOverlapNWORDS);
  };

private:
  uint32                  _bufferLen;    //  length of valid data in the buffer
  uint32                  _bufferPos;    //  position the read is at in the buffer
  uint32                  _bufferMax;    //  allocated size of the buffer
  uint32                 *_buffer;

  bool                    _isOutput;     //  if true, we can writeOverlap()
  bool                    _isSeekable;   //  if true, we can seekOverlap()
  bool                    _isNormal;     //  if true, 3 words per overlap, else 4

  compressedFileReader   *_reader;
  compressedFileWriter   *_writer;

  FILE                   *_file;
};





class ovStoreInfo {
private:
  uint64    _ovsMagic;
  uint64    _ovsVersion;
  uint64    _numOverlapsPerFile;  //  on create, how big to make each file
  uint64    _smallestIID;         //  smallest frag iid in the store
  uint64    _largestIID;          //  largest frag iid in the store
  uint64    _numOverlapsTotal;    //  number of overlaps in the store
  uint64    _highestFileIndex;
  uint64    _maxReadLenInBits;    //  length of a fragment

  friend class ovStore;

  friend
  void       writeOverlaps(char       *storePath,
                           ovsOverlap *ovls,
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
                ovsOverlap *ovls,
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
  ovStore(const char *name, ovStoreType cType=ovStoreReadOnly);
  ~ovStore();

  //  Read the next overlap from the store.  Return value is the number of overlaps read.
  uint32     readOverlap(ovsOverlap *overlap);

  //  Read ALL remaining overlaps for the current A_iid.  Return value is the number of overlaps read.
  uint32     readOverlaps(ovsOverlap *overlaps,
                          uint32      maxOverlaps,
                          bool        restrictToIID=true);

  //  Append ALL remaining overlaps for the current A_iid, as long as A_iid == iid, to the overlaps
  //  in ovl.  Return value is the number of overlaps in ovl that are for 'iid'.
  //
  //  It is up to the client to verify that ovl[0] is the same as iid; ovlLen is the number of
  //  overlaps in ovl, NOT the number of overlaps in ovl that are the same as iid.
  //
  uint32       readOverlaps(uint32       iid,
                            ovsOverlap *&ovl,
                            uint32      &ovlLen,
                            uint32      &ovlMax);

  void       setRange(uint32 low, uint32 high);
  void       resetRange(void);

  uint64     numOverlapsInRange(void);
  uint32 *   numOverlapsPerFrag(uint32 &firstFrag, uint32 &lastFrag);

  //  The (mostly) private interface for adding overlaps to a store.  Overlaps must be sorted already.

  void       writeOverlap(ovsOverlap *olap);
  void       writeOverlap(ovsOverlap *overlap, uint32 maxOverlapsThisFile);

  //  Write a block of sorted overlaps to store file 'fileID', saving the info and index into
  //  'fileID.info' and 'fileID.index'

  friend
  void       writeOverlaps(char       *storePath,
                           ovsOverlap *ovls,
                           uint64      ovlsLen,
                           uint32      fileID);


  //  Add new erates for reads between bgnID and endID.  No checking of IDs is done, but the number
  //  of erates must agree.

  void       addErates(uint32 bgnID, uint32 endID, uint16 *erates, uint64 eratesLen);

private:
  char               _storePath[FILENAME_MAX];

  bool               _isOutput;

  ovStoreInfo        _info;

  uint32             _firstIIDrequested;
  uint32             _lastIIDrequested;

  //memoryMappedFile  *_offtMMap;   //  For reading overlaps, a mapped ovStoreOfft file.
  //ovStoreOfft       *_offts;      //  For reading overlaps, an array of mapped ovStoreOfft's.
  //uint64             _offtLength; //  For reading overlaps, the number of ovStoreOfft's.

  FILE              *_offtFile;   //  For writing overlaps, a place to dump ovStoreOfft's.
  ovStoreOfft        _offt;       //  For writing overlaps, the current ovStoreOfft.
  ovStoreOfft        _offm;       //  For writing overlaps, an empty ovStoreOfft, for reads with no overlaps.

  memoryMappedFile  *_eratesMap;
  uint16            *_erates;

  uint64             _overlapsThisFile;  //  Count of the number of overlaps written so far
  uint32             _currentFileIndex;
  ovFile            *_bof;
};


//  This should be part of ovStore, but when it is used, in ovStoreSorter, we don't
//  have a store opened.
void
writeOverlaps(char       *storePath,
              ovsOverlap *ovls,
              uint64      ovlsLen,
              uint32      fileID);

bool
testIndex(char     *storePath,
          bool      doFixes);

void
mergeInfoFiles(char       *storePath,
               uint32      nPieces);







//  For store construction.  Probably should be in either ovsOverlap or ovStore.

class ovStoreFilter {
public:
  ovStoreFilter(gkStore *gkp_, double maxErate) {
    gkp             = gkp_;

    resetCounters();

    maxID     = gkp->gkStore_getNumReads() + 1;
    maxEvalue = AS_OVS_encodeQuality(maxErate);

    skipReadOBT     = new char [maxID];
    skipReadDUP     = new char [maxID];

    memset(skipReadOBT, 0, sizeof(char) * maxID);
    memset(skipReadDUP, 0, sizeof(char) * maxID);


    uint32  numSkipOBT;
    uint32  numSkipDUP;

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
          (L->gkLibrary_finalTrim()                == FINALTRIM_NONE) &&
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

    fprintf(stderr, "Marked "F_U64" reads so skip OBT, "F_U64" reads to skip dedupe.\n", numSkipOBT, numSkipDUP);
  };

  ~ovStoreFilter() {
    delete [] skipReadOBT;
    delete [] skipReadDUP;
  };


  void    filterOverlap(ovsOverlap     &foverlap,
                        ovsOverlap     &roverlap);

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
