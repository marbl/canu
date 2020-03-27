
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef AS_OVSTORE_H
#define AS_OVSTORE_H

#include "runtime.H"
#include "files.H"

#include "sqStore.H"

#include "ovOverlap.H"
#include "ovStoreFile.H"
#include "ovStoreHistogram.H"



const uint64 ovStoreVersion         = 4;
const uint64 ovStoreMagic           = 0x53564f3a756e6163;   //  == "canu:OVS - store complete
//const uint64 ovStoreMagicIncomplete = 0x50564f3a756e6163;   //  == "canu:OVP - store under construction

#define  OVSTORE_MEMORY_OVERHEAD     (256 * 1024 * 1024)



class ovStoreInfo {
public:
  ovStoreInfo(uint32 maxID=0) {
    clear(maxID);
  };

  void     clear(uint32 maxID) {
    _ovsMagic      = 0;
    _ovsVersion    = 0;
    _readLenInBits = AS_MAX_READLEN_BITS;
    _bgnID         = UINT32_MAX;
    _endID         = 0;
    _maxID         = maxID;
    _numOlaps      = 0;
  };

  void       load(const char *path, uint32 index=UINT32_MAX, bool temporary=false) {
    char    name[FILENAME_MAX];
    uint32  failed = 0;

    if (temporary == false)
      snprintf(name, FILENAME_MAX, "%s/info", path);
    else
      snprintf(name, FILENAME_MAX, "%s/%04u.info", path, index);

    AS_UTL_loadFile(name, this, 1);

    if (_ovsMagic != ovStoreMagic)
      failed += fprintf(stderr, "ERROR:  directory '%s' is not an ovStore.\n", path);

    if (_ovsVersion != ovStoreVersion)
      failed += fprintf(stderr, "ERROR:  directory '%s' is not a supported ovStore version (store version " F_U64 "; supported version " F_U64 ".\n",
                        path, _ovsVersion, ovStoreVersion);

    if (_readLenInBits != AS_MAX_READLEN_BITS)
      failed += fprintf(stderr, "ERROR:  directory '%s' is not a supported read length (store is " F_U32 " bits, AS_MAX_READLEN_BITS is " F_U32 ").\n",
                        path, _readLenInBits, AS_MAX_READLEN_BITS);

    if (failed)
      exit(1);
  };

  void       save(const char *path, uint32 index=UINT32_MAX, bool temporary=false) {
    char  name[FILENAME_MAX];

    if (temporary == false)
      snprintf(name, FILENAME_MAX, "%s/info", path);
    else
      snprintf(name, FILENAME_MAX, "%s/%04u.info", path, index);

    _ovsMagic   = ovStoreMagic;
    _ovsVersion = ovStoreVersion;

    if (_numOlaps == 0) {
      fprintf(stderr, "WARNING:\n");
      fprintf(stderr, "WARNING:  No overlaps loaded in to the store.\n");
      fprintf(stderr, "WARNING:\n");
    }

    else {
      assert(_bgnID    != UINT32_MAX);
      assert(_endID    != 0);
      assert(_maxID    != 0);
      assert(_numOlaps != 0);
    }

    //fprintf(stderr, "ovStoreInfo::save()-- bgnID=%u endID=%u maxID=%u numOlaps=%lu\n",
    //        _bgnID, _endID, _maxID, _numOlaps);

    AS_UTL_saveFile(name, this, 1);
  };

  uint32     bgnID(void)  { return(_bgnID); };
  uint32     endID(void)  { return(_endID); };
  uint32     maxID(void)  { return(_maxID); };

  void       addOverlaps(uint32 curID, uint32 nOverlaps=1)   {
    _bgnID = min(_bgnID, curID);
    _endID = max(_endID, curID);

    _numOlaps += nOverlaps;
  };
  uint64     numOverlaps(void) {
    return(_numOlaps);
  };

private:
  uint64    _ovsMagic;
  uint64    _ovsVersion;

  uint32    _readLenInBits;

  uint32    _bgnID;               //  First ID with overlaps
  uint32    _endID;               //  Last ID with overlaps
  uint32    _maxID;               //  ID of the last read in the assembly.

  uint64    _numOlaps;            //  number of overlaps in the store
};



class ovStoreOfft {
public:
  ovStoreOfft() {
    _slice     = 0;
    _piece     = 0;
    _offset    = 0;
    _numOlaps  = 0;
    _overlapID = 0;
  };

  void       addOverlap(uint32 slice, uint32 piece, uint32 offset, uint64 overlapID) {
    if (_numOlaps == 0) {      //  If the first overlap to be added,
      _slice     = slice;      //  set all the good info.
      _piece     = piece;
      _offset    = offset;
      _numOlaps  = 0;
      _overlapID = overlapID;
    }

    _numOlaps++;
  }

  uint16    _slice;           //  Which slice are these overlaps in?
  uint16    _piece;           //  Which piece are these overlaps in?
  uint32    _offset;          //  Offset (in overlaps) in the piece file.
  uint32    _numOlaps;        //  number of overlaps for this iid

  uint64    _overlapID;       //  index into erates for this block.
};



//  For sequential construction, there is only a constructor, destructor and writeOverlap().
//  Overlaps must be sorted by a_iid (then b_iid) already.

class ovStoreWriter {
public:
  ovStoreWriter(const char *path, sqStore *seq);
  ~ovStoreWriter();

  void                writeOverlap(ovOverlap *olap);

private:
  char               _storePath[FILENAME_MAX+1];

  ovStoreInfo        _info;
  sqStore           *_seq;

  ovStoreOfft       *_index;

  ovFile            *_bof;
  uint32             _bofSlice;
  uint32             _bofPiece;

  ovStoreHistogram  *_histogram;         //  When constructing a sequential store, collects all the stats from each file
};



//  For parallel construction, usage is much more complicated.  The constructor
//  will write a single file of sorted overlaps, and each file has it's own metadata.
//  After all files are written, the metadata is merged into one file.

class ovStoreSliceWriter {
public:
  ovStoreSliceWriter(const char *path, sqStore *seq, uint32 sliceNum, uint32 numSlices, uint32 numBuckets);
  ~ovStoreSliceWriter();

  uint64       loadBucketSizes(uint64 *bucketSizes);
  void         loadOverlapsFromBucket(uint32 bucket, uint64 expectedLen, ovOverlap *ovls, uint64& ovlsLen);

  void         writeOverlaps(ovOverlap *ovls, uint64 ovlsLen);

  void         mergeInfoFiles(void);
  void         mergeHistogram(void);

  void         removeOverlapSlice(void);
  void         checkSortingIsComplete(void);
  void         removeAllIntermediateFiles(void);

private:
  char               _storePath[FILENAME_MAX+1];

  ovStoreInfo        _info;
  sqStore           *_seq;

  uint32             _sliceNum;
  uint32             _pieceNum;
  uint32             _numSlices;
  uint32             _numBuckets;
};



class ovStore {
public:
  ovStore(const char *name, sqStore *seq);
  ~ovStore();

public:
  void               testStore(bool verbose=true);

public:
  //  Read the next overlap from the store.  Return value is the number of overlaps read.
  uint32             readOverlap(ovOverlap *overlap);

  //  Loads the overlaps for a single read, returning the number of overlaps loaded.
  uint32             loadOverlapsForRead(uint32       id,
                                         ovOverlap  *&ovl,
                                         uint32      &ovlMax);

  //  Try not to use this interface.  It's gross.  Then again, so is the
  //  previous one.  The intent was to load exactly ovlMax overlaps, but the
  //  implementation requires all overlaps for a read to be loaded, so we end
  //  up with fewer than ovlMax overlaps.
  uint32             loadBlockOfOverlaps(ovOverlap *&ovl,
                                         uint32     &ovlMax);

  void               setRange(uint32 bgnID, uint32 endID);

  void               restartIteration(void);    //  UNTESTED, probably needs to seekOverlap() too
  void               endIteration(void);

  uint32             numOverlaps(uint32 readID)   {  return(_index[readID]._numOlaps);  };
  uint64             numOverlapsInRange(void);
  uint32            *numOverlapsPerRead(void);

  //  Add new evalues for reads between bgnID and endID.  No checking of IDs is done, but the number
  //  of evalues must agree.

  void               addEvalues(vector<char const *> &fileList);

  //  Return the statistics associated with this store

  ovStoreHistogram  *getHistogram(void) {
    return(new ovStoreHistogram(_storePath));
  };

public:
  void                dumpMetaData(uint32 bgnID, uint32 endID);

private:
  char               _storePath[FILENAME_MAX+1];

  ovStoreInfo        _info;
  sqStore           *_seq;

  uint32             _bgnID;    //  First ID requested
  uint32             _endID;    //  Last ID requested

  uint32             _curID;    //  Current ID being read
  uint32             _curOlap;  //  Current overlap being read (0 .. N)

  ovStoreOfft       *_index;

  memoryMappedFile  *_evaluesMap;
  uint16            *_evalues;

  ovFile            *_bof;
  uint32             _bofSlice;
  uint32             _bofPiece;
};





//  For store construction.  Probably should be in either ovOverlap or ovStore.

class ovStoreFilter {
public:
  ovStoreFilter(sqStore *seq_, double maxErate);
  ~ovStoreFilter();

  void     filterOverlap(ovOverlap     &foverlap,
                         ovOverlap     &roverlap);

  void     resetCounters(void);

  uint64   savedUnitigging(void)    { return(saveUTG);      };
  uint64   savedTrimming(void)      { return(saveOBT);      };

  uint64   filteredNoTrim(void)     { return(skipOBT);      };
  uint64   filteredErate(void)      { return(skipERATE);    };
  uint64   filteredFlipped(void)    { return(skipFLIPPED);  };

public:
  sqStore *seq;

  uint32   maxID;
  uint32   maxEvalue;

  uint64   saveUTG;
  uint64   saveOBT;

  uint64   skipOBT;        //  OBT not requested for the A read
  uint64   skipERATE;
  uint64   skipFLIPPED;

private:
  char    *skipReadOBT;    //  State of the filter.
};


#endif  //  AS_OVSTORE_H
