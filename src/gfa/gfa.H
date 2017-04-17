
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

#ifndef AS_UTL_GFA_H
#define AS_UTL_GFA_H

#include "AS_global.H"
#include "splitToWords.H"



//  Features assumed to hold only the length, and we don't use it.

class gfaSequence {
public:
  gfaSequence();
  gfaSequence(char *inLine);
  ~gfaSequence();

  void    load(char *inLine);
  void    save(FILE *outFile);

public:
  char   *_name;
  uint32  _id;
  char   *_sequence;
  char   *_features;

  uint32  _length;
};




class gfaLink {
public:
  gfaLink();
  gfaLink(char *inLine);
  ~gfaLink();

  void    load(char *inLine);
  void    save(FILE *outFile);

  void    alignmentLength(int32 &queryLen, int32 &refceLen, int32 &alignLen);

public:
  char   *_Aname;
  uint32  _Aid;      //  Canu specific.
  bool    _Afwd;

  char   *_Bname;
  uint32  _Bid;      //  Canu specific.
  bool    _Bfwd;

  char   *_cigar;

  char   *_features;
};




class gfaFile {
public:
  gfaFile(char *inFile);
  ~gfaFile();

  bool    loadFile(char *inFile);
  bool    saveFile(char *outFile);

public:
  char                  *_header;

  vector<gfaSequence *>  _sequences;
  vector<gfaLink *>      _links;
};




#endif  //  AS_UTL_GFA_H
