
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


class merylOperation;


class merylInput {
public:
  merylInput(merylOperation *o);
  merylInput(char *n, kmerCountFileReader *s);
  merylInput(char *n, dnaSeqFile *f);
  ~merylInput();

  void   initialize(void);
  void   nextMer(void);

  bool   isFromOperation(void)    { return(_operation != NULL); };
  bool   isFromDatabase(void)     { return(_stream    != NULL); };
  bool   isFromSequence(void)     { return(_sequence  != NULL); };

  merylOperation        *_operation;
  kmerCountFileReader   *_stream;
  dnaSeqFile            *_sequence;

  char                   _name[FILENAME_MAX];

  kmer                   _kmer;
  uint64                 _count;
  bool                   _valid;
};

#endif  //  MERYLINPUT_H
