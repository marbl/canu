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

#include <cassert>
#include <iostream>
#include "atacreader.H"

using namespace std;

ATACreader::ATACreader(MatchExtenderParameters *p,
                       string infile,
                       string outfile) {
  _params = p;

  _isExtender = (outfile != "");
  //  _nextMatch is uninitialized
  _eof        = false;

  _atacData = new AtacDataset();
  _atacData->setInputFile(infile);

  if (outfile != "")
    _atacData->setOutputFile(outfile);

  _s1Reader    = 0L;
  _prefix1     = "";
  _asmID1      = "";
  _curAxis1    = "";
  _curAxis1Seq = 0L;
  _asmID1Len   = 0;

  _s2Reader    = 0L;
  _prefix2     = "";
  _asmID2      = "";
  _curAxis2    = "";
  _curAxis2Seq = 0L;
  _asmID2Len   = 0;

  _numRead    = 0;
  _numWritten = 0;
}


ATACreader::~ATACreader() {
  delete _atacData;
  delete _s1Reader;
  delete _s2Reader;
}

  
void 
ATACreader::setFastaFiles(string f1name, string f2name) {
  _s1Reader = new FastACache(f1name.c_str(), 1, false, false);
  _s2Reader = new FastACache(f2name.c_str(), 1, false, false);
}


bool 
ATACreader::processPreamble(void) {

  _atacData->readGlobalData();


  if (_atacData->hasAttribute("assemblyFilePrefix1"))
    _prefix1 = _atacData->getAttributeString("assemblyFilePrefix1");
  if (_atacData->hasAttribute("assemblyFilePrefix2"))
    _prefix2 = _atacData->getAttributeString("assemblyFilePrefix2");


  if (_atacData->hasAttribute("assemblyId1")) {
    _asmID1 = _atacData->getAttributeString("assemblyId1");
    _asmID1Len = _asmID1.size();
  } else {
    return false;
  }

  if (_atacData->hasAttribute("assemblyId2")) {
    _asmID2 = _atacData->getAttributeString("assemblyId2");
    _asmID2Len = _asmID2.size();
  } else {
    return false;
  }


  if (_isExtender) {
    // Append parameters to global info if necessary
    if (_params) {
      _atacData->putAttribute("matchExtenderMinEndRunLen",  _params->minEndRunLen);
      _atacData->putAttribute("matchExtenderMaxMMBlock",    _params->maxMMBlock);
      _atacData->putAttribute("matchExtenderMinBlockSep",   _params->minBlockSep);
      _atacData->putAttribute("matchExtenderMinIdentity",   _params->minIdentity);
      _atacData->putAttribute("matchExtenderMaxNbrSep",     _params->maxNbrSep);
      _atacData->putAttribute("matchExtenderMaxNbrPathMM",  _params->maxNbrPathMM);
    }
  
    // Write global information back out at the beginning of the file
    //
    _atacData->writeGlobalData();
  }


  //  Read the first match
  //
  fprintf(stderr, "Reading the next match!\n");
  _eof = !(_atacData->readNextMatch(_nextMatch));

  return true;
}


bool 
ATACreader::getNextMatchBatches(vector<MEMatch *>& fwd_matches,
                                vector<MEMatch *>& rev_matches,
                                FastASequenceInCore *s1,
                                FastASequenceInCore *s2,
                                string& axis1,
                                string& axis2) {
  MEMatch *m;
  string axis_iid;

  assert(_s1Reader && _s2Reader);

  // Retrieve sequences

  axis1 = _nextMatch.getGenomicAxisId(0);
  axis2 = _nextMatch.getGenomicAxisId(1);

  if (axis1 != _curAxis1) {
    //delete _curAxis1Seq;
    if (_params->beVerbose) 
      cerr << " * Fetching sequence for axis1 = " << axis1 << endl;
    axis_iid = axis1.substr(_asmID1Len + 1);
    _curAxis1 = axis1;
    _curAxis1Seq = _s1Reader->getSequence(atoi(axis_iid.c_str()));
    if (_curAxis1Seq == 0L) {
      cerr << " *   ERROR!" << endl;
    }
    if (_params->beVerbose)
      cerr << " * Done" << endl;
  }
  s1 = _curAxis1Seq;
  
  if (axis2 != _curAxis2) { 
    //delete _curAxis2Seq;
    if (_params->beVerbose) 
      cerr << " * Fetching sequence for axis2 = " << axis2 << endl;
    axis_iid = axis2.substr(_asmID2Len + 1);
    _curAxis2 = axis2;
    _curAxis2Seq = _s2Reader->getSequence(atoi(axis_iid.c_str()));
    if (_curAxis2Seq == 0L) {
      cerr << " *   ERROR!" << endl;
    }
    if (_params->beVerbose)
      cerr << " * Done" << endl;
  }
  s2 = _curAxis2Seq;

  fwd_matches.clear();
  rev_matches.clear();
  
  if (_params->beVerbose)
    cerr << " * Reading matches for these axes ..." << endl;

  do {

    //  If we are a mismatch counter, skip all "r" matches
    //
    if ((!_isExtender) && (_nextMatch.getTypeId() == "r"))
      continue;

    if (_nextMatch.getLength(0) != _nextMatch.getLength(1)) {
      fprintf(stderr, "_nextMatch.getLength(0)=%d != _nextMatch.getLength(1)=%d\n",
              _nextMatch.getLength(0),
              _nextMatch.getLength(1));
    }

    if (_nextMatch.getOrientation(0) != 1) {
      fprintf(stderr, "WARNING:  getOrientation(0) is not 1!\n");
    }

    if ((_nextMatch.getOrientation(1) != 1) &&
        (_nextMatch.getOrientation(1) != -1)) {
      fprintf(stderr, "WARNING:  getOrientation(1) is not 1 or -1!\n");
    }

    _numRead++;

    m = MEMatch::getNewMatch(_nextMatch.getId(),
			     s1, s2, 
			     (SeqOff_t) _nextMatch.getStart(0),
			     (SeqOff_t) _nextMatch.getStart(1),
			     (SeqOff_t) _nextMatch.getLength(0),
			     (_nextMatch.getOrientation(1) == -1));

    if (m->isReversed()) {
      rev_matches.push_back(m);
    }
    else {
      fwd_matches.push_back(m);
    }
  }while (!(_eof = !(_atacData->readNextMatch(_nextMatch))) &&
	  (_nextMatch.getGenomicAxisId(0) == axis1) &&
	  (_nextMatch.getGenomicAxisId(1) == axis2));

  if (_params->beVerbose)
    cerr << " * Done" << endl;

  return true;
}



//  This only for "MatchExtender"
void 
ATACreader::writeMatchBatch(vector<MEMatch *>& matches,
                            const string& axis1,
                            const string& axis2) {
  unsigned int d;
  MEMatch *mm;

  for (d = 0; d < matches.size(); ++d) {
    mm = matches[d];
    if (mm->isDeleted())
      continue;
      
    AtacMatch am("u", mm->id(), "", axis1, (int64) mm->pos1(),
		 (int64) mm->len(), 1, axis2, (int64) mm->pos2(),
		 (int64) mm->len(), 
		 (mm->isReversed() ? -1 : 1));
    
    _atacData->writeNextRecord(&am);
    _numWritten++;
  }
}
