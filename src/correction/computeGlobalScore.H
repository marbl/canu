
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

#ifndef COMPUTEGLOBALSCORES_H
#define COMPUTEGLOBALSCORES_H

#include "runtime.H"
#include "ovStore.H"

class globalScoreStats {
public:
  globalScoreStats() {
    totalOverlaps = 0;
    lowErate     = 0;
    highErate    = 0;
    tooShort     = 0;
    tooLong      = 0;
    belowCutoff  = 0;
    retained     = 0;

    reads00OlapsFiltered  = 0;
    reads50OlapsFiltered  = 0;
    reads80OlapsFiltered  = 0;
    reads95OlapsFiltered  = 0;
    reads99OlapsFiltered  = 0;
  };

  uint64      totalOverlaps;
  uint64      lowErate;
  uint64      highErate;
  uint64      tooShort;
  uint64      tooLong;
  uint64      belowCutoff;
  uint64      retained;

  uint64      reads00OlapsFiltered;
  uint64      reads50OlapsFiltered;
  uint64      reads80OlapsFiltered;
  uint64      reads95OlapsFiltered;
  uint64      reads99OlapsFiltered;
};



class globalScore {
public:
  globalScore(uint32  minOvlLength_,
              uint32  maxOvlLength_,
              double  minErate_,
              double  maxErate_,
              FILE   *logFile_ = NULL,
              bool    doStats  = false) {
    hist         = NULL;
    histLen      = 0;
    histMax      = 0;

    stats        = NULL;

    minOvlLength = minOvlLength_;
    maxOvlLength = maxOvlLength_;
    minEvalue    = AS_OVS_encodeEvalue(minErate_);
    maxEvalue    = AS_OVS_encodeEvalue(maxErate_);

    logFile      = logFile_;
    stats        = (doStats) ? new globalScoreStats : NULL;
  };

  ~globalScore() {
    delete [] hist;
    delete    stats;
  };

  uint16    compute(uint32             ovlLen,
                    ovOverlap         *ovl,
                    uint32             expectedCoverage,
                    uint32             thresholdsLen,
                    uint16            *thresholds);

  void      estimate(uint32            ovlLen,
                     uint32            expectedCoverage);

  uint64      totalOverlaps(void)           { return(stats->totalOverlaps); };
  uint64      lowErate(void)                { return(stats->lowErate);      };
  uint64      highErate(void)               { return(stats->highErate);     };
  uint64      tooShort(void)                { return(stats->tooShort);      };
  uint64      tooLong(void)                 { return(stats->tooLong);       };
  uint64      belowCutoff(void)             { return(stats->belowCutoff);   };
  uint64      retained(void)                { return(stats->retained);      };

  uint64      reads00OlapsFiltered(void)    { return(stats->reads00OlapsFiltered); };
  uint64      reads50OlapsFiltered(void)    { return(stats->reads50OlapsFiltered); };
  uint64      reads80OlapsFiltered(void)    { return(stats->reads80OlapsFiltered); };
  uint64      reads95OlapsFiltered(void)    { return(stats->reads95OlapsFiltered); };
  uint64      reads99OlapsFiltered(void)    { return(stats->reads99OlapsFiltered); };

private:
  uint64            *hist;
  uint32             histLen;
  uint32             histMax;

  globalScoreStats  *stats;

  uint32             minOvlLength;
  uint32             maxOvlLength;
  uint32             minEvalue;
  uint32             maxEvalue;

  FILE              *logFile;
};


#endif  //  COMPUTEGLOBALSCORES_H
