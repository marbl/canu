/**************************************************************************
 * This file is part of A2Amapper.
 * Copyright (c) 2004 Applera Corporation
 * Author: Clark Mobarry
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**************************************************************************/


#ifndef heavyChainsMod_h
#define heavyChainsMod_h

#include <string>
typedef unsigned int u32bit;

struct Match {
  int xlo,ylo;
  int xhi,yhi;

  double selfS; // The intrinsic score of the Match.
  double S;     // The computed score of the Match?

  // We need two times the number of dimensions of scores.  That is
  // one score starting from each corner of the bounding box of the
  // space.
  double neS;
  double swS;
  double nwS;
  double seS;
  int filled;  // Is this the same as selfS ?
  char ori;
};

class StrandPair {
private:
  std::vector<Match> P;

public:
  int old_stra1;
  int old_stra2;
  int VERBOSE;
  int MAXJUMP;  // Default maximum intra-run jump allowed in a good run.
  double MINSCORE; // Default minimum of bp filled in a good run.

  // The following are only known after StrandPair::print().
  double sumlen1;
  double sumlen2;
  double maxlen1;
  double maxlen2;
  double maxScoreFwd; // Maximum forward chain score for the strand pair.
  double maxScoreRev; // Maximum reverse chain score for the strand pair.

  StrandPair(int _verbose, int _maxjump, double _minscore);
  ~StrandPair(void);
  void addHit  (
   char   direction,
   u32bit id1,
   u32bit xlo,
   u32bit xln,
   u32bit id2,
   u32bit ylo,
   u32bit yln,
   u32bit filled);
  void process(void);
  long print(
	     FILE *outF,
	     const int output_in_brians_format,
	     const char * assemblyId1, 
	     const char * assemblyId2,
             long matchid
	     );
  void clear(void);
  size_t size(void) const;
};

class TheStats {
private:
  double sumlen1;
  double sumlen2;
  double sumMaxLen1;
  double sumMaxLen2;
  double sumMaxScoreFwd;
  double sumMaxScoreRev;
public:
  TheStats(std::map<std::string, std::string> & globals);
  ~TheStats(void);
  void add(StrandPair *sp);
  void print(FILE *file) const;
  std::map<std::string,std::string> & _globals;
};

#endif // heavyChainsMod_h
