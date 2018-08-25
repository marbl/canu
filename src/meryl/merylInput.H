
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#ifndef MERYLINPUT_H
#define MERYLINPUT_H

#include "meryl.H"
#include "kmers.H"
#include "sequence.H"

#include "sqStore.H"


class merylOperation;


class merylInput {
public:
  merylInput(merylOperation *o);
  merylInput(const char *n, kmerCountFileReader *s, uint32 fileNum=UINT32_MAX);
  merylInput(const char *n, dnaSeqFile *f);
  merylInput(const char *n, sqStore *s, uint32 segment, uint32 segmentMax);
  ~merylInput();

  char  *inputType(void) {
    if (_operation)   return("meryl-operation");
    if (_stream)      return("meryl-database");
    if (_sequence)    return("sequence-file");
    if (_store)       return("canu-seqStore");

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
  bool   isFromStore(void)        { return(_store     != NULL); };

  merylOperation           *_operation;
  kmerCountFileReader      *_stream;
  dnaSeqFile               *_sequence;
  sqStore                  *_store;

  char                      _name[FILENAME_MAX+1];

  //  For _operation and _stream, a copy of the 'active' kmer

  kmer                      _kmer;
  uint64                    _count;
  bool                      _valid;

  //  For _store, the position we're at in the store.

  uint32                    _sqBgn;
  uint32                    _sqEnd;

  sqRead                   *_read;
  sqReadData               *_readData;
  uint32                    _readID;
  uint32                    _readPos;
};

#endif  //  MERYLINPUT_H
