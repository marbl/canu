
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

#ifndef BED_H
#define BED_H

#include "strings.H"

#include <vector>


class bedRecord {
public:
  bedRecord();
  bedRecord(char *inLine);
  ~bedRecord();

  void    load(char *inLine);
  void    save(FILE *outFile);

public:
  char   *_Aname;    //  The 'chromosome'
  uint32  _Aid;      //  Canu specific.

  int32   _bgn;
  int32   _end;

  char   *_Bname;    //  The 'feature'
  uint32  _Bid;      //  Canu specific.

  uint32  _score;
  bool    _Bfwd;
};


class bedFile {
public:
  bedFile(char *inName);
  ~bedFile();

  bool    loadFile(char *inName);
  bool    saveFile(char *outName);

public:
  std::vector<bedRecord *>    _records;
};




#endif  //  BED_H
