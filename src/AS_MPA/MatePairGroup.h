
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
/* $Id: MatePairGroup.h,v 1.2 2004-09-23 20:25:23 mcschatz Exp $ */
#ifndef MATEPAIRGROUP_H
#define MATEPAIRGROUP_H

#include "MatePair.h"

class MatePairGroup
{
public:
  MatePairGroup()
    {
      // setOrientation(PAIR_INNIE);
      setType(MPI_UNKNOWN);
    }
  MatePairGroup(const MatePair & mp)
    {
      appendMP(mp);
      // setOrientation(mp.getOrientation());
      setType(MPI_UNKNOWN);
    }
  MatePairGroup(const MatePairGroup & mpg)
    {
      for(unsigned int i = 0; i < mpg.getNumMPs(); i++)
        appendMP(mpg.getMP(i));
      // setOrientation(mpg.getOrientation());
      setType(mpg.getMPIndex());
    }
  
  unsigned int getNumMPs() const {return pmps.size();}
  const MatePair & getMP(int index) const {return pmps[index];}
  void appendMP(const MatePair & mp) {pmps.push_back(mp);}
  void reset() {pmps.clear();}

  bool isSatisfied() const {return(pmpi == MPI_SATISFIED);}
  bool isStretched() const {return(pmpi == MPI_STRETCHED);}
  bool isCompressed() const {return(pmpi == MPI_COMPRESSED);}
  bool isOuttie() const {return pmpi == MPI_OUTTIE;}
  bool isNormal() const {return pmpi == MPI_NORMAL;}
  bool isAntinormal() const {return pmpi == MPI_ANTINORMAL;}
  bool isUnsatisfied() const {return(pmpi != MPI_SATISFIED);}
  bool isInversion() const {return(pmpi == MPI_INVERSION);}
  bool isTransposition() const {return(pmpi == MPI_TRANSPOSITION);}
  bool isInnie() const
    {
      return(pmpi == MPI_STRETCHED ||
             pmpi == MPI_COMPRESSED ||
             pmpi == MPI_SATISFIED);
    }
  
  MatePairIndex_e getType() const {return pmpi;}
  MatePairIndex_e getMPIndex() const {return getType();}
  // PairOrientation_e getOrientation() const {return pOrient;}
  
  void setType(MatePairIndex_e mpi) {pmpi = mpi;}
  void setMPIndex(MatePairIndex_e mpi) {setType(mpi);}
  // void setOrientation(PairOrientation_e orient) {pOrient - orient;}
  
protected:
  vector<MatePair> pmps;
  // PairOrientation_e pOrient;
  MatePairIndex_e pmpi;
};

#endif // #ifndef MATEPAIRGROUP_H
