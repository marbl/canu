
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
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SEQFILE_H
#define SEQFILE_H

#undef DEBUG

#include "AS_global.H"
#include "readBuffer.H"

//  General flow of the constructors is:
//    Clear all data
//    Open the file
//    Set _filename, _typename
//    Read/build the index structure
//    Position the file to the first read
//    Set _numberOfSequences (IMPORTANT, and subtle)

class seqFile {
protected:
  seqFile(const char *UNUSED(filename)) {};
  seqFile() {};

public:
  virtual ~seqFile() {};

protected:
  virtual seqFile      *openFile(const char *filename) = 0;

public:
  virtual const char   *getSourceName(void)    { return(_filename); };
  virtual const char   *getFileTypeName(void)  { return(_typename); };

  virtual bool          randomAccessSupported(void) { return(_randomAccessSupported); };

  virtual uint32        getNumberOfSequences(void) { return(_numberOfSequences); };

public:
  virtual uint32        find(const char *sequencename) = 0;

  virtual uint32        getSequenceLength(uint32 id) = 0;
  virtual bool          getSequence(uint32 id,
                                    char *&h, uint32 &hLen, uint32 &hMax,
                                    char *&s, uint32 &sLen, uint32 &sMax) = 0;
  virtual bool          getSequence(uint32 iid,
                                    uint32 bgn, uint32 end, char *s) = 0;

protected:
  char                 _filename[FILENAME_MAX];
  char                 _typename[FILENAME_MAX];

  bool                 _randomAccessSupported;

  uint32               _numberOfSequences;

  friend class seqFactory;
};

#endif //  SEQFILE_H
