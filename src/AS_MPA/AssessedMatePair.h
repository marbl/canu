
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
/* $Id: AssessedMatePair.h,v 1.2 2004-09-23 20:25:23 mcschatz Exp $ */
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
                   uint64 libUID) :
    MatePair(leftFrag, rightFrag, libUID)
    {
      _mpi = MPI_UNKNOWN;
    }

  AssessedMatePair(const FragmentPosition & leftFrag,
                   const FragmentPosition & rightFrag,
                   PairOrientation_e orient,
                   uint64 libUID) :
    MatePair(leftFrag, rightFrag, orient, libUID)
    {
      _mpi = MPI_UNKNOWN;
    }
  
  AssessedMatePair(const FragmentPosition & leftFrag,
                   const FragmentPosition & rightFrag,
                   PairOrientation_e orient,
                   uint64 libUID,
                   int32 tooFar,
                   int32 tooClose) :
    MatePair(leftFrag, rightFrag, orient, libUID)
    {
      assess(tooFar, tooClose);
    }

  void assess(int32 tooFar, int32 tooClose)
    {
      if(getLeftChromosome() != getRightChromosome())
      {
        _mpi = MPI_INTERCHROMOSOME;
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
           COORD_TYPE left5, COORD_TYPE right5,
           int32 leftChrom, int32 rightChrom,
           COORD_TYPE tooFar, COORD_TYPE tooClose)
    {
      FragmentPosition fpl(leftUID, left5,
                           (po == PAIR_INNIE || po == PAIR_NORMAL),
                           leftChrom);
      FragmentPosition fpr(rightUID, right5,
                          (po == PAIR_OUTTIE || po == PAIR_NORMAL),
                           rightChrom);
      setLibUID(libUID);
      setLeftFrag(fpl);
      setRightFrag(fpr);
      setOrientation(po);
      assess(tooFar, tooClose);
    }
  void set(PairOrientation_e po,
           ID_TYPE leftUID, ID_TYPE rightUID, ID_TYPE libUID,
           COORD_TYPE left5, COORD_TYPE right5, int32 chrom,
           COORD_TYPE tooFar, COORD_TYPE tooClose)
    {
      FragmentPosition fpl(leftUID, left5,
                           (po == PAIR_INNIE || po == PAIR_NORMAL));
      FragmentPosition fpr(rightUID, right5,
                          (po == PAIR_OUTTIE || po == PAIR_NORMAL));
      setLibUID(libUID);
      setLeftFrag(fpl);
      setRightFrag(fpr);
      setChromosome(chrom);
      setOrientation(po);
      assess(tooFar, tooClose);
    }
  MatePairIndex_e getType() const {return _mpi;}
    
private:
  MatePairIndex_e _mpi;
};

#endif // #ifndef ASSESSEDMATEPAIR_H
