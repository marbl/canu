
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
/* $Id: FragmentPosition.h,v 1.2 2004-09-23 20:25:23 mcschatz Exp $ */
#ifndef FRAGMENTPOSITION_H
#define FRAGMENTPOSITION_H

#include "MPTypes.h"
#include "Orientation.h"

#define COORD_TYPE  int64

class FragmentPosition
{
public:
  FragmentPosition(){}

  FragmentPosition(const FragmentPosition & fp)
    {
      pUID = fp.pUID;
      pFiveP = fp.pFiveP;
      pFragOrient = fp.pFragOrient;
      pChrom = fp.pChrom;
    }
  
  FragmentPosition(ID_TYPE uid,
                   COORD_TYPE fiveP,
                   bool pointsRight,
                   int32 pChrom = -1)
    {
      pUID = uid;
      pFiveP = fiveP;
      pFragOrient = (pointsRight ? SINGLE_A_B : SINGLE_B_A);
    }

  FragmentPosition(ID_TYPE uid,
                   COORD_TYPE fiveP,
                   SingleOrientation_e fragOrient,
                   int32 pChrom = -1)
    {
      pUID = uid;
      pFiveP = fiveP;
      pFragOrient = fragOrient;
      pChrom = -1;
    }

  void setUID(ID_TYPE uid) {pUID = uid;}
  void setFiveP(COORD_TYPE fiveP) {pFiveP = fiveP;}
  void setOrientation(SingleOrientation_e fragOrient)
    {
      pFragOrient = fragOrient;
    }
  void setChromosome(int32 c) {pChrom = c;}
  void set(ID_TYPE uid,
           COORD_TYPE fiveP,
           SingleOrientation_e fragOrient = SINGLE_A_B,
           int32 chrom = -1)
    {
      setUID(uid);
      setFiveP(fiveP);
      setOrientation(fragOrient);
      setChromosome(chrom);
    }
  void set(ID_TYPE uid, COORD_TYPE fiveP, char * fragOrient, int32 chrom = -1)
    {
      if(strcmp(fragOrient, "A_B") == 0)
        set(uid, fiveP, SINGLE_A_B, chrom);
      else
      {
        assert(strcmp(fragOrient, "B_A") == 0);
        set(uid, fiveP, SINGLE_B_A, chrom);
      }
    }

  ID_TYPE getUID() const {return pUID;}
  COORD_TYPE getFiveP() const {return pFiveP;}
  SingleOrientation_e getOrientation() const {return pFragOrient;}
  bool pointsRight() const {return pFragOrient == SINGLE_A_B;}
  bool pointsLeft() const {return pFragOrient == SINGLE_B_A;}
  int32 getChromosome() const {return pChrom;}
  
private:
  ID_TYPE pUID;
  COORD_TYPE pFiveP;
  SingleOrientation_e pFragOrient;
  int32  pChrom;
};

#endif
