
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef MERYLINPUT_H
#define MERYLINPUT_H

#include "meryl.H"
#include "kmers.H"
#include "sequence.H"

#ifdef CANU
#include "sqStore.H"
#endif

class merylOperation;


class merylInput {
public:
  merylInput(merylOperation *o);
  merylInput(const char *n, kmerCountFileReader *s, uint32 fileNum=UINT32_MAX);
  merylInput(const char *n, dnaSeqFile *f);
#ifdef CANU
  merylInput(const char *n, sqStore *s, uint32 segment, uint32 segmentMax);
#endif
  ~merylInput();

  char  *inputType(void) {
    if (_operation)   return("meryl-operation");
    if (_stream)      return("meryl-database");
    if (_sequence)    return("sequence-file");
#ifdef CANU
    if (_store)       return("canu-seqStore");
#endif

    return("invalid-input");
  };

  void   initialize(void);
  void   nextMer(void);

  bool   loadBases(char    *seq,
                   uint64   maxLength,
                   uint64  &seqLength,
                   bool    &endOfSequence);

  bool   isFromOperation(void)    { return(_operation != NULL); };
  bool   isFromDatabase(void)     { return(_stream    != NULL); };
  bool   isFromSequence(void)     { return(_sequence  != NULL); };
#ifdef CANU
  bool   isFromStore(void)        { return(_store     != NULL); };
#else
  bool   isFromStore(void)        { return(false);              };
#endif

  bool   isMultiSet(void)         { return(_isMultiSet);        };

  merylOperation           *_operation;
  kmerCountFileReader      *_stream;
  dnaSeqFile               *_sequence;
#ifdef CANU
  sqStore                  *_store;
#endif

  char                      _name[FILENAME_MAX+1];
  bool                      _isMultiSet;

  //  For _operation and _stream, a copy of the 'active' kmer

  kmer                      _kmer;
  uint64                    _value;
  bool                      _valid;

  //  For _store, the position we're at in the store.

#ifdef CANU
  uint32                    _sqBgn;
  uint32                    _sqEnd;

  sqRead                   *_read;
  sqReadData               *_readData;
  uint32                    _readID;
  uint32                    _readPos;
#endif
};

#endif  //  MERYLINPUT_H
