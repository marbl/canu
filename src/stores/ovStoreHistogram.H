
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
 *    Brian P. Walenz beginning on 2016-OCT-25
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_OVSTOREHISTOGRAM_H
#define AS_OVSTOREHISTOGRAM_H

//  Automagically gathers statistics on overlaps as they're written:
//    from overlappers, the number of overlaps per read.
//    in the store, the number of overlaps per (evalue,overlapLength)

#include "AS_global.H"
#include "gkStore.H"

#include "ovStoreFile.H"



#define  N_OVL_SCORE   16   //  Number of overlap scores to save per read

class oSH_ovlSco {
public:
  oSH_ovlSco()  {}
  ~oSH_ovlSco() {};

  uint16    points[N_OVL_SCORE];
  uint16    scores[N_OVL_SCORE];
};



class ovStoreHistogram {
private:
  void    initialize(gkStore *gkp=NULL);
public:
  void    clear(void);

public:
  ovStoreHistogram();                                //  Used when loading data, user must loadData() later
  ovStoreHistogram(char *path);                      //  Used when loading data, calls loadData() for you
  ovStoreHistogram(gkStore *gkp, ovFileType type);   //  Used when writing ovFile
  ~ovStoreHistogram();

#if 0
  double    minErate(void)             {  return(AS_OVS_decodeEvalue(0));           };
  double    maxErate(void)             {  return(AS_OVS_decodeEvalue(_maxEvalue));  };

  uint32    minEvalue(void)            {  return(0);           };
  uint32    maxEvalue(void)            {  return(_maxEvalue);  };
#endif

  uint32    numEvalueBuckets(void)     {  return(AS_MAX_EVALUE + 1);  };
  uint32    numLengthBuckets(void)     {  return(_opelLen);           };

  //uint32    minOverlapLength(void)     {  return(0);                  };
  //uint32    maxOverlapLength(void)     {  return(_maxLength * _bpb);  };

  uint32    evaluePerBucket(void)      {  return(_epb);  };
  uint32    basesPerBucket(void)       {  return(_bpb);  };

  uint32    numOverlaps(uint32 eb, uint32 lb) {
    assert(eb < numEvalueBuckets());
    assert(lb < numLengthBuckets());

    return((_opel[eb] == NULL) ? 0 : _opel[eb][lb]);
  };

  //uint32    numOverlaps(uint32 id);
  //uint32    numOverlaps(uint32 evalue, uint32 length);

  //  In an ovFile, add a single value to the histogram

  void      addOverlap(ovOverlap *overlap);

  //  In an ovStore, load the histogram saved in a file, and add it to our current data.

  void      saveData(char *prefix);
  void      loadData(char *prefix, uint32 maxIID=UINT32_MAX);

  //  Remove data associated with some prefix.

  static
  void      removeData(char *prefix);

  //  Add in the data from histogram 'input' to this histogram

  void      add(ovStoreHistogram *input);

  //  Copies the number of overlaps per read into oprOut.  This array is assumed to sized using
  //  gkpStore to get the number of reads.  The data in the histogram can be shorter, but shouldn't
  //  be longer.  If so, it will fail and exit.
  //
  //  Returns total number of overlaps in this histogram.

  uint64    getOverlapsPerRead(uint32 *oprOut, uint32 oprOutLen);

  //  Dump a gnuplot-friendly data file of the evalues-length.

  void      dumpEvalueLength(FILE *out);

  //  Report the ith overlap score estimate.  The top i overlaps will (usually) have score bigger
  //  than this.

  uint16    overlapScoreEstimate(uint32 id, uint32 i);

private:
  gkStore   *_gkp;

  uint32     _maxOlength;   //  Max overlap length seen - NOT USED
  uint32     _maxEvalue;    //  Max evalue seen         - NOT USED

  //  Overlaps per evalue-length, for overlaps in the store.

  uint32     _epb;      //  Evalues per bucket
  uint32     _bpb;      //  Bases per bucket

  uint32     _opelLen;  //  Length of the data vector for one evalue
  uint32   **_opel;     //  Overlaps per evalue-length

  //  Overlaps per read, for overlaps out of overlapper.

  uint32     _oprLen;   //  Length of opr valid data
  uint32     _oprMax;   //  Last allocated opr
  uint32    *_opr;      //  Overlaps per read

  //  Overlap score for the top overlaps.  Used during correction.
  //  Want to store ~11 values per read, 16 bits each, so 22 bytes.
  //  Human has 14,625,060 reads -> 160,875,660 bytes data.
  //  UINT32_MAX indicates the no more data.

  void          processScores(uint32 Aid=UINT32_MAX);

  uint32       _scoresListLen;
  uint32       _scoresListMax;
  uint16      *_scoresList;
  uint32       _scoresListAid;  //  Current ID being stored in _scoresList

  uint32       _scoresBgn;      //  First ID stored in _scores
  uint32       _scoresLen;      //  Number of valid scores
  uint32       _scoresMax;      //  Number of allocated scores
  oSH_ovlSco  *_scores;
};

#endif  //  AS_OVSTOREHISTOGRAM_H
