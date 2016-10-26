
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
  void       writeOverlaps(gkStore    *gkp,
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

  //  Write a block of sorted overlaps to store file 'fileID', saving the info and index into
  //  'fileID.info' and 'fileID.index'

  friend
  void         writeOverlaps(gkStore    *gkp,
                             char       *storePath,
                             ovOverlap  *ovls,
                             uint64      ovlsLen,
                             uint32      fileID);


  //  Add new evalues for reads between bgnID and endID.  No checking of IDs is done, but the number
  //  of evalues must agree.

  void       addEvalues(uint32 bgnID, uint32 endID, uint16 *evalues, uint64 evaluesLen);

  //  Return the statistics associated with this store

  ovStoreHistogram  *getHistogram(void) {
    return(new ovStoreHistogram(_gkp, _storePath));
  };

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

  ovStoreHistogram  *_histogram;         //  When constructing a sequential store, collects all the stats from each file

  gkStore           *_gkp;
};


//  This should be part of ovStore, but when it is used, in ovStoreSorter, we don't
//  have a store opened.
void
writeOverlaps(gkStore    *gkp,
              char       *storePath,
              ovOverlap  *ovls,
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

    fprintf(stderr, "Marked " F_U32 " reads so skip OBT, " F_U32 " reads to skip dedupe.\n", numSkipOBT, numSkipDUP);
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



#endif  //  AS_OVSTORE_H
