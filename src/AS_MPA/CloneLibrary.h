
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
 *************************************************************************/
/* $Id: CloneLibrary.h,v 1.1.1.1 2004-04-14 13:51:56 catmandew Exp $ */
#ifndef CLONELIBRARY_H
#define CLONELIBRARY_H

#include <iostream>
#include <cstdio>

#include "MPTypes.h"
#include "Orientation.h"

class CloneLibrary
{
public:
  CloneLibrary() {pUID = 0; pMean = 0; pStddev = 0; pOrient = PAIR_INNIE;}

  CloneLibrary(uint64 uid, double mean, double stddev,
               PairOrientation_e orient = PAIR_INNIE)
    {
      pUID = uid;
      pMean = mean;
      pStddev = stddev;
      pOrient = orient;
    }

  void setFromString(char * line)
    {
      sscanf(line, _LU " %lf %lf", &pUID, &pMean, &pStddev);
      pOrient = PAIR_INNIE;
    }

  void setUID(uint64 uid) {pUID = uid;}
  void setMean(double mean) {pMean = mean;}
  void setStddev(double stddev) {pStddev = stddev;}
  void setOrientation(PairOrientation_e orient) {pOrient = orient;}

  uint64 getUID() const {return pUID;}
  double getMean() const {return pMean;}
  double getStddev() const {return pStddev;}
  PairOrientation_e getOrientation() const {return pOrient;}

  friend ostream & operator<<(ostream & os, const CloneLibrary & cl)
    {
      os << cl.pUID << " " << cl.pMean << " " << cl.pStddev;
      return os;
    }
private:
  uint64  pUID;
  double pMean;
  double pStddev;
  PairOrientation_e pOrient;
};

#endif
