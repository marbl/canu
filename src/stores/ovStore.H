
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

#ifndef AS_OVSTORE_H
#define AS_OVSTORE_H

#include "AS_global.H"
#include "gkStore.H"

#define SNAPPY

#include "ovOverlap.H"
#include "ovStoreFile.H"
#include "ovStoreHistogram.H"



const uint64 ovStoreVersion         = 2;
const uint64 ovStoreMagic           = 0x53564f3a756e6163;   //  == "canu:OVS - store complete
const uint64 ovStoreMagicIncomplete = 0x50564f3a756e6163;   //  == "canu:OVP - store under construction


class ovStoreInfo {
public:
  ovStoreInfo() {
    clear();
  };

  void     clear(void) {
    _ovsMagic         = ovStoreMagicIncomplete;  //  Appropriate for a new store.
    _ovsVersion       = ovStoreVersion;
    _UNUSED           = 0;
    _smallestIID      = UINT64_MAX;
    _largestIID       = 0;
    _numOverlapsTotal = 0;
    _highestFileIndex = 0;
    _maxReadLenInBits = AS_MAX_READLEN_BITS;
  };

  bool       load(const char *path, uint32 index=UINT32_MAX, bool temporary=false) {
    char  name[FILENAME_MAX];

    if (temporary == false)
      snprintf(name, FILENAME_MAX, "%s/info", path);
    else
      snprintf(name, FILENAME_MAX, "%s/%04u.info", path, index);

    if (AS_UTL_fileExists(name, false, false) == false) {
      fprintf(stderr, "ERROR: directory '%s' is not an overlapStore; didn't find file '%s': %s\n",
              path, name, strerror(errno));
      return(false);
    }

    errno = 0;
    FILE *ovsinfo = fopen(name, "r");
    if (errno) {
      fprintf(stderr, "ERROR: directory '%s' is not an overlapStore; failed to open '%s': %s\n",
              path, name, strerror(errno));
      return(false);
    }

    AS_UTL_safeRead(ovsinfo, this, "ovStore::ovStore::info", sizeof(ovStoreInfo), 1);

    fclose(ovsinfo);

    return(true);
  };

  bool       test(const char *path) {
    char  name[FILENAME_MAX];

    snprintf(name, FILENAME_MAX, "%s/info", path);

    if (AS_UTL_fileExists(name, false, false) == false)
      return(false);

    errno = 0;
    FILE *ovsinfo = fopen(name, "r");
    if (errno)
      fprintf(stderr, "ERROR: failed to load '%s'; can't check if this is a valid ovStore: %s\n",
              name, strerror(errno)), exit(1);

    AS_UTL_safeRead(ovsinfo, this, "ovStore::ovStore::info", sizeof(ovStoreInfo), 1);

    return(checkMagic());
  };

  void       save(const char *path, uint32 index=UINT32_MAX, bool temporary=false) {
    char  name[FILENAME_MAX];

    if (temporary == false)
      snprintf(name, FILENAME_MAX, "%s/info", path);
    else
      snprintf(name, FILENAME_MAX, "%s/%04u.info", path, index);

    if (temporary == false) {
      _ovsMagic         = ovStoreMagic;
      _ovsVersion       = ovStoreVersion;
      _highestFileIndex = index;
    } else {
    }

    errno = 0;
    FILE *ovsinfo = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: failed to save '%s': %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(ovsinfo, this, "ovStore::ovStore::saveinfo", sizeof(ovStoreInfo), 1);

    fclose(ovsinfo);
  };

  void       addOverlap(uint32 id, uint32 nOverlaps=1) {
    if (_smallestIID > id) _smallestIID = id;
    if (_largestIID  < id) _largestIID  = id;

    _numOverlapsTotal += nOverlaps;
  };

  bool       checkIncomplete(void)    { return(_ovsMagic         == ovStoreMagicIncomplete);  };
  bool       checkMagic(void)         { return(_ovsMagic         == ovStoreMagic);            };
  bool       checkVersion(void)       { return(_ovsVersion       == ovStoreVersion);          };
  bool       checkSize(void)          { return(_maxReadLenInBits == AS_MAX_READLEN_BITS);     };

  uint32     getVersion(void)         { return((uint32)_ovsVersion);          };
  uint32     getCurrentVersion(void)  { return((uint32)ovStoreVersion);       };
  uint32     getSize(void)            { return((uint32)_maxReadLenInBits);    };

  uint64     numOverlaps(void)        { return(_numOverlapsTotal); };
  uint32     smallestID(void)         { return(_smallestIID);      };
  uint32     largestID(void)          { return(_largestIID);       };

  uint32     lastFileIndex(void)      { return(_highestFileIndex); };

private:
  uint64    _ovsMagic;
  uint64    _ovsVersion;
  uint64    _UNUSED;              //  needed to keep the file layout the same
  uint64    _smallestIID;         //  smallest frag iid in the store
  uint64    _largestIID;          //  largest frag iid in the store
  uint64    _numOverlapsTotal;    //  number of overlaps in the store
  uint64    _highestFileIndex;
  uint64    _maxReadLenInBits;    //  length of a fragment
};



class ovStoreOfft {
public:
  ovStoreOfft() {
    clear();
  };
  ~ovStoreOfft() {
  };

  void       clear(void) {
    _a_iid     = 0;
    _fileno    = 0;
    _offset    = 0;
    _numOlaps  = 0;
    _overlapID = 0;
  };

private:
  uint32    _a_iid;      //  read ID for this block of overlaps.

  uint32    _fileno;     //  the file that contains this a_iid
  uint32    _offset;     //  offset to the first overlap for this iid
  uint32    _numOlaps;   //  number of overlaps for this iid

  uint64    _overlapID;  //  overlapID for the first overlap in this block.  in memory, this is the id of the next overlap.

  friend class ovStore;
  friend class ovStoreWriter;

  friend
  void
  writeOverlaps(gkStore    *gkp,
                char       *storePath,
                ovOverlap  *ovls,
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



class ovStoreWriter {
public:
  ~ovStoreWriter();

  //  For sequential construction, there is only a constructor, destructor and writeOverlap().
  //  Overlaps must be sorted by a_iid (then b_iid) already.

  ovStoreWriter(const char *path, gkStore *gkp);

  void         writeOverlap(ovOverlap *olap);

  //  For parallel construction, usage is much more complicated.  The constructor
  //  will write a single file of sorted overlaps, and each file has it's own metadata.
  //  After all files are written, the metadata is merged into one file.

  ovStoreWriter(const char *path, gkStore *gkp, uint32 fileLimit, uint32 fileID, uint32 jobIdxMax);

  void         writeOverlaps(ovOverlap *ovls, uint64 ovlsLen);

  uint64       loadBucketSizes(uint64 *bucketSizes);
  void         loadOverlapsFromSlice(uint32 slice, uint64 expectedLen, ovOverlap *ovls, uint64& ovlsLen);
  void         removeOverlapSlice(void);

  void         mergeInfoFiles(void);
  void         mergeHistogram(void);

  bool         testIndex(bool doFixes);

  void         checkSortingIsComplete(void);
  void         removeAllIntermediateFiles(void);

private:
  char               _storePath[FILENAME_MAX];

  ovStoreInfo        _info;
  gkStore           *_gkp;

  FILE              *_offtFile;   //  For writing overlaps, a place to dump ovStoreOfft's.
  ovStoreOfft        _offt;       //  For writing overlaps, the current ovStoreOfft.
  ovStoreOfft        _offm;       //  For writing overlaps, an empty ovStoreOfft, for reads with no overlaps.

  memoryMappedFile  *_evaluesMap;
  uint16            *_evalues;

  uint64             _overlapsThisFile;  //  Count of the number of overlaps written so far
  uint64             _overlapsThisFileMax;
  uint32             _currentFileIndex;
  ovFile            *_bof;

  ovStoreHistogram  *_histogram;         //  When constructing a sequential store, collects all the stats from each file

  //  Parallel store support

  uint32             _fileLimit;   //  number of slices used in bucketizing/sorting
  uint32             _fileID;      //  index of the overlap file we're processing
  uint32             _jobIdxMax;   //  total number of overlap files
};



class ovStore {
public:
  ovStore(const char *name, gkStore *gkp);
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

  //  Add new evalues for reads between bgnID and endID.  No checking of IDs is done, but the number
  //  of evalues must agree.

  void       addEvalues(vector<char *> &fileList);
  void       addEvalues(uint32 bgnID, uint32 endID, uint16 *evalues, uint64 evaluesLen);

  //  Return the statistics associated with this store

  ovStoreHistogram  *getHistogram(void) {
    return(new ovStoreHistogram(_gkp, _storePath));
  };

private:
  char               _storePath[FILENAME_MAX];

  ovStoreInfo        _info;
  gkStore           *_gkp;

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
};





//  For store construction.  Probably should be in either ovOverlap or ovStore.

class ovStoreFilter {
public:
  ovStoreFilter(gkStore *gkp_, double maxErate);
  ~ovStoreFilter();

  void     filterOverlap(ovOverlap     &foverlap,
                         ovOverlap     &roverlap);

  void     resetCounters(void);

  uint64   savedUnitigging(void)    { return(saveUTG);      };
  uint64   savedTrimming(void)      { return(saveOBT);      };
  uint64   savedDedupe(void)        { return(saveDUP);      };

  uint64   filteredErate(void)      { return(skipERATE);    };

  uint64   filteredNoTrim(void)     { return(skipOBT);      };
  uint64   filteredBadTrim(void)    { return(skipOBTbad);   };
  uint64   filteredShortTrim(void)  { return(skipOBTshort); };

  uint64   filteredNoDedupe(void)   { return(skipDUP);      };
  uint64   filteredNotDupe(void)    { return(skipDUPdiff);  };
  uint64   filteredDiffLib(void)    { return(skipDUPlib);   };

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

  char    *skipReadOBT;    //  State of the filter.
  char    *skipReadDUP;
};


#endif  //  AS_OVSTORE_H
