// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Dan Fasulo
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//  20050223 - bpw - Removed RASCAL (Sorry, Dan!  This was the only
//  thing that used it.)

#ifndef ATACREADER_H
#define ATACREADER_H

#include <vector>
#include <map>
#include <string>
#include "AtacMatch.h"
#include "AtacComment.h"
#include "AtacAttribute.h"
#include "AtacDataset.h"
#include "MEMatch.H"

#include "bio++.H"


struct MatchExtenderParameters {
  int     minEndRunLen;
  int     maxMMBlock;
  int     minBlockSep;
  double  minIdentity;
  int     maxNbrSep;
  int     maxNbrPathMM;
  bool    beVerbose;
};



class ATACreader {
public:
  ATACreader(MatchExtenderParameters *p,
             string infile,
             string outfile = "");
  ~ATACreader();


  //  Sets the Fasta files to use.  Each is identified by an (index
  //  file, sequence file) pair.  Returns true if both are valid.
  //
  void setFastaFiles(string f1name, string f2name);
  

  // Reads and echos all Atac data up to the first match in the file.
  //
  bool processPreamble(void);


  //  Reads the next match, and all subsequent matches between
  //  scaffolds with the same IDs.  Retrieves the sequences associated
  //  with the matches and places the MEMatch structures into vectors
  //  indexed by diagonal ID.  These vectors are dynamically allocated
  //  to the proper size and should be deleted by the caller when they
  //  are no longer needed.
  //
  bool getNextMatchBatches(vector<MEMatch *>& fwd_matches,
			   vector<MEMatch *>& rev_matches,
			   FastASequenceInCore *s1,
                           FastASequenceInCore *s2,
                           string& axis1,
                           string& axis2);

  //  Returns true iff all matches have been processed.
  //
  bool endOfMatches(void) const { return(_eof); };

  //  Writes the batch of matches to the output file.  The IDs of the
  //  two axes must be provided.
  //
  void writeMatchBatch(vector<MEMatch *>& matches,
		       const string& axis1,
                       const string& axis2);

  //  Write a tag giving the count of the mismatch blocks.  Since this
  //  is generally the last attribute set, it also writes all of the
  //  global data to the Atac file.
  //
  void writeMismatchCount(unsigned long count);

  //  Returns the file prefix associated with the first or second
  //  sequence collection, if present in the Atac file.  Otherwise,
  //  returns "".
  //
  string filePrefix1(void) const { return(_prefix1); }
  string filePrefix2(void) const { return(_prefix2); }
  
  inline u64bit numRead(void)    const { return _numRead; }
  inline u64bit numWritten(void) const { return _numWritten; }

protected:
  
  MatchExtenderParameters  *_params;

  bool                      _isExtender;
  AtacMatch                 _nextMatch;
  bool                      _eof;
  AtacDataset              *_atacData;

  FastACache              *_s1Reader;
  string                    _prefix1;
  string                    _asmID1;
  string                    _curAxis1;
  FastASequenceInCore      *_curAxis1Seq;
  unsigned int              _asmID1Len;

  FastACache               *_s2Reader;
  string                    _prefix2;
  string                    _asmID2;
  string                    _curAxis2;
  FastASequenceInCore      *_curAxis2Seq;
  unsigned int              _asmID2Len;

  u64bit                    _numRead, _numWritten;
};

#endif  //  ATACREADER_H
