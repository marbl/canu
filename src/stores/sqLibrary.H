
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
 *    src/stores/gkLibrary.H
 *    src/stores/gkStore.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-26 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SQLIBRARY_H
#define SQLIBRARY_H


//  DO NOT INCLUDE THIS FILE DIRECTLY, include sqStore.H.


//  Per-library options.

//  Read type

#define SQ_READTYPE_GENERIC            0x0000
#define SQ_READTYPE_CONTIG             0x0001
#define SQ_READTYPE_PACBIO_RAW         0x0002
#define SQ_READTYPE_PACBIO_CORRECTED   0x0003
#define SQ_READTYPE_NANOPORE_RAW       0x0004
#define SQ_READTYPE_NANOPORE_CORRECTED 0x0005
#define SQ_READTYPE_PACBIO_HIFI        0x0006

//  Correction algorithm

#define SQ_CORRECTION_NONE             0x0000
#define SQ_CORRECTION_CONSENSUS        0x0001
#define SQ_CORRECTION_MER              0x0002

//  Trimming algorithm

#define SQ_FINALTRIM_NONE              0x0000
#define SQ_FINALTRIM_LARGEST_COVERED   0x0001  //  largest region covered by good overlaps
#define SQ_FINALTRIM_BEST_EDGE         0x0002  //  largest region covered by best overlaps (broken)



class sqLibrary {
public:
  sqLibrary() {
    memset(_libraryName, 0, sizeof(char) * LIBRARY_NAME_SIZE);
    strncpy(_libraryName, "UNDEFINED", LIBRARY_NAME_SIZE-1);

    _libraryID = UINT32_MAX;

    sqLibrary_clearFeatures();
  };

  ~sqLibrary() {
  };

  void sqLibrary_clearFeatures(void) {
    //  DO NOT change defaults without updating sqLibrary_parsePreset().
    _isNonRandom                = 0;

    _readCorrection             = SQ_CORRECTION_NONE;
    _readType                   = SQ_READTYPE_GENERIC;

    _finalTrim                  = SQ_FINALTRIM_LARGEST_COVERED;
    _removeDuplicateReads       = 1;
    _removeSpurReads            = 1;
    _removeChimericReads        = 1;
    _checkForSubReads           = 1;

    _defaultQV                  = 20;
  };

public:

  char const    *sqLibrary_libraryName(void)           { return(_libraryName);          };
  uint32         sqLibrary_libraryID(void)             { return(_libraryID);            };

  uint32         sqLibrary_isNonRandom(void)           { return(_isNonRandom);          };
  uint32         sqLibrary_readType(void)              { return(_readType);             };
  char const    *sqLibrary_readTypeString(void);
  uint32         sqLibrary_readCorrection(void)        { return(_readCorrection);       };
  char const    *sqLibrary_readCorrectionString(void);
  uint32         sqLibrary_finalTrim(void)             { return(_finalTrim);            };
  char const    *sqLibrary_finalTrimString(void);
  uint32         sqLibrary_removeDuplicateReads(void)  { return(_removeDuplicateReads); };
  uint32         sqLibrary_removeSpurReads(void)       { return(_removeSpurReads);      };
  uint32         sqLibrary_removeChimericReads(void)   { return(_removeChimericReads);  };
  uint32         sqLibrary_checkForSubReads(void)      { return(_checkForSubReads);     };

  uint32         sqLibrary_defaultQV(void)             { return(_defaultQV);            };


  void           sqLibrary_setIsNonRandom(bool f)           { _isNonRandom = f;          };
  void           sqLibrary_setReadType(char *f);
  void           sqLibrary_setReadCorrection(char *t);
  void           sqLibrary_setFinalTrim(char *t);
  void           sqLibrary_setRemoveDuplicateReads(bool f)  { _removeDuplicateReads = f; };
  void           sqLibrary_setRemoveSpurReads(bool f)       { _removeSpurReads = f;      };
  void           sqLibrary_setRemoveChimericReads(bool f)   { _removeChimericReads = f;  };
  void           sqLibrary_setCheckForSubReads(bool f)      { _checkForSubReads = f;     };

  void           sqLibrary_setDefaultQV(double qv)          { _defaultQV = qv;           };

  void           sqLibrary_parsePreset(char *t);


private:
  char           _libraryName[LIBRARY_NAME_SIZE];
  uint32         _libraryID;

  //  If set, reads are from a non-random library and shouldn't count toward coverage stats.
  uint32         _isNonRandom;

  //  What generated these reads?
  uint32         _readType;

  //  Should reads be corrected?  How?
  uint32         _readCorrection;

  //  Should reads be trimmed based on overlaps?  How?
  uint32         _finalTrim;

  //  Should duplicate reads (based on overlaps) be removed?
  uint32         _removeDuplicateReads;

  //  Should spur reads be cleaned up?  How?
  uint32         _removeSpurReads;

  //  Should chimeric reads be cleaned up?  How?
  uint32         _removeChimericReads;

  //  Should PacBio circular sub-reads be cleaned up?  How?
  uint32         _checkForSubReads;

  //  For reads with no QVs, use this.
  uint32         _defaultQV;

  friend class sqStore;
};


#endif  //  SQLIBRARY_H
