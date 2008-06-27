
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
/* $Id: AssessedMatePair.h,v 1.7 2008-06-27 06:29:16 brianwalenz Exp $ */
#ifndef ASSESSEDMATEPAIR_H
#define ASSESSEDMATEPAIR_H

#include "MatePair.h"
#include "CloneLibrary.h"

class AssessedMatePair :
  public MatePair
{
public:
  AssessedMatePair() :
    MatePair()
    {
    }

  AssessedMatePair(const FragmentPosition & leftFrag,
                   const FragmentPosition & rightFrag,
                   uint64_t libUID) :
    MatePair(leftFrag, rightFrag, libUID)
    {
      _mpi = MPI_UNKNOWN;
    }

  AssessedMatePair(const FragmentPosition & leftFrag,
                   const FragmentPosition & rightFrag,
                   PairOrientation_e orient,
                   uint64_t libUID) :
    MatePair(leftFrag, rightFrag, orient, libUID)
    {
      _mpi = MPI_UNKNOWN;
    }

  AssessedMatePair(const FragmentPosition & leftFrag,
                   const FragmentPosition & rightFrag,
                   PairOrientation_e orient,
                   uint64_t libUID,
                   int32 tooFar,
                   int32 tooClose) :
    MatePair(leftFrag, rightFrag, orient, libUID)
    {
      assess(tooFar, tooClose);
    }

  void assess(int32 tooFar, int32 tooClose)
    {
      if(getLeftSequenceID() != getRightSequenceID())
      {
        _mpi = MPI_INTERSEQUENCE;
        return;
      }
      switch(getOrientation())
      {
        case PAIR_INNIE:
          if(getRightCoord() - getLeftCoord() > tooFar)
          {
            _mpi = MPI_STRETCHED;
          }
          else if(getRightCoord() - getLeftCoord() < tooClose)
          {
            _mpi = MPI_COMPRESSED;
          }
          else
          {
            _mpi = MPI_SATISFIED;
          }
          break;
        case PAIR_OUTTIE:
          _mpi = MPI_OUTTIE;
          break;
        case PAIR_NORMAL:
          _mpi = MPI_NORMAL;
          break;
        case PAIR_ANTINORMAL:
          _mpi = MPI_ANTINORMAL;
          break;
        default:
          _mpi = MPI_UNKNOWN;
          break;
      }
    }

  void set(PairOrientation_e po,
           ID_TYPE leftUID, ID_TYPE rightUID, ID_TYPE libUID,
           UNIT_TYPE left5, UNIT_TYPE right5,
           int32 leftSeqID, int32 rightSeqID,
           UNIT_TYPE tooFar, UNIT_TYPE tooClose)
    {
      FragmentPosition fpl(leftUID, left5,
                           (po == PAIR_INNIE || po == PAIR_NORMAL),
                           leftSeqID);
      FragmentPosition fpr(rightUID, right5,
                          (po == PAIR_OUTTIE || po == PAIR_NORMAL),
                           rightSeqID);
      setLibUID(libUID);
      setLeftFrag(fpl);
      setRightFrag(fpr);
      setOrientation(po);
      assess(tooFar, tooClose);
    }
  void set(PairOrientation_e po,
           ID_TYPE leftUID, ID_TYPE rightUID, ID_TYPE libUID,
           UNIT_TYPE left5, UNIT_TYPE right5, int32 seqID,
           UNIT_TYPE tooFar, UNIT_TYPE tooClose)
    {
      FragmentPosition fpl(leftUID, left5,
                           (po == PAIR_INNIE || po == PAIR_NORMAL));
      FragmentPosition fpr(rightUID, right5,
                          (po == PAIR_OUTTIE || po == PAIR_NORMAL));
      setLibUID(libUID);
      setLeftFrag(fpl);
      setRightFrag(fpr);
      setSequenceID(seqID);
      setOrientation(po);
      assess(tooFar, tooClose);
    }
  MatePairIndex_e getType() const {return _mpi;}

private:
  MatePairIndex_e _mpi;
};

#endif // #ifndef ASSESSEDMATEPAIR_H
