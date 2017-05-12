
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
 *    Brian P. Walenz beginning on 2017-APR-04
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_UTL_BED_H
#define AS_UTL_BED_H

#include "AS_global.H"
#include "splitToWords.H"



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
  bedFile(char *inFile);
  ~bedFile();

  bool    loadFile(char *inFile);
  bool    saveFile(char *outFile);

public:
  vector<bedRecord *>    _records;
};




#endif  //  AS_UTL_BED_H
