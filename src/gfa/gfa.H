
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

#ifndef AS_UTL_GFA_H
#define AS_UTL_GFA_H

#include "runtime.H"
#include "strings.H"



//  Features assumed to hold only the length, and we don't use it.

class gfaSequence {
public:
  gfaSequence();
  gfaSequence(char *inLine);
  gfaSequence(char *name, uint32 id, uint32 len);
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
  gfaLink(char *Aname, uint32 Aid, bool Afwd,
          char *Bname, uint32 Bid, bool Bfwd, char *cigar);
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
  gfaFile();
  gfaFile(char const *inName);
  ~gfaFile();

  bool    loadFile(char const *inName);
  bool    saveFile(char const *outName);

public:
  char                  *_header;

  vector<gfaSequence *>  _sequences;
  vector<gfaLink *>      _links;
};




#endif  //  AS_UTL_GFA_H
