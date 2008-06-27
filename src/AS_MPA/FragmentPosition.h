
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
/* $Id: FragmentPosition.h,v 1.6 2008-06-27 06:29:16 brianwalenz Exp $ */
#ifndef FRAGMENTPOSITION_H
#define FRAGMENTPOSITION_H

#include "MPTypes.h"
#include "Orientation.h"

class FragmentPosition
{
public:
  FragmentPosition(){}

  FragmentPosition(const FragmentPosition & fp)
    {
      pUID = fp.pUID;
      pFiveP = fp.pFiveP;
      pFragOrient = fp.pFragOrient;
      pSeqID = fp.pSeqID;
    }

  FragmentPosition(ID_TYPE uid,
                   UNIT_TYPE fiveP,
                   bool pointsRight,
                   ID_TYPE pSeqID = BOGUS_ID)
    {
      pUID = uid;
      pFiveP = fiveP;
      pFragOrient = (pointsRight ? SINGLE_A_B : SINGLE_B_A);
    }

  FragmentPosition(ID_TYPE uid,
                   UNIT_TYPE fiveP,
                   SingleOrientation_e fragOrient,
                   ID_TYPE pSeqID = BOGUS_ID)
    {
      pUID = uid;
      pFiveP = fiveP;
      pFragOrient = fragOrient;
    }

  void setUID(ID_TYPE uid) {pUID = uid;}
  void setFiveP(UNIT_TYPE fiveP) {pFiveP = fiveP;}
  void setOrientation(SingleOrientation_e fragOrient)
    {
      pFragOrient = fragOrient;
    }
  void setSequenceID(ID_TYPE c) {pSeqID = c;}
  void set(ID_TYPE uid,
           UNIT_TYPE fiveP,
           SingleOrientation_e fragOrient = SINGLE_A_B,
           ID_TYPE seqID = BOGUS_ID)
    {
      setUID(uid);
      setFiveP(fiveP);
      setOrientation(fragOrient);
      setSequenceID(seqID);
    }
  void set(ID_TYPE uid,
           UNIT_TYPE fiveP,
           char * fragOrient,
           ID_TYPE seqID = BOGUS_ID)
    {
      if(strcmp(fragOrient, "A_B") == 0)
        set(uid, fiveP, SINGLE_A_B, seqID);
      else
      {
        assert(strcmp(fragOrient, "B_A") == 0);
        set(uid, fiveP, SINGLE_B_A, seqID);
      }
    }

  ID_TYPE getUID() const {return pUID;}
  UNIT_TYPE getFiveP() const {return pFiveP;}
  SingleOrientation_e getOrientation() const {return pFragOrient;}
  bool pointsRight() const {return pFragOrient == SINGLE_A_B;}
  bool pointsLeft() const {return pFragOrient == SINGLE_B_A;}
  ID_TYPE getSequenceID() const {return pSeqID;}

private:
  ID_TYPE pUID;
  UNIT_TYPE pFiveP;
  SingleOrientation_e pFragOrient;
  ID_TYPE pSeqID;
};

#endif
