
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
/* $Id: MatePair.h,v 1.2 2004-09-23 20:25:23 mcschatz Exp $ */
#ifndef MATEPAIR_H
#define MATEPAIR_H

#include <iostream>

#include "MPTypes.h"
#include "FragmentPosition.h"
#include "Orientation.h"

class MatePair
{
public:
  MatePair(){}
  MatePair(const FragmentPosition & leftFrag,
           const FragmentPosition & rightFrag,
           ID_TYPE libUID) :
    pLeftFrag(leftFrag), pRightFrag(rightFrag)
    {
      pLibUID = libUID;
      setOrientFromFrags();
    }
  MatePair(const FragmentPosition & leftFrag,
           const FragmentPosition & rightFrag,
           PairOrientation_e orient,
           ID_TYPE libUID) :
    pLeftFrag(leftFrag), pRightFrag(rightFrag)
    {
      pLibUID = libUID;
      pOrient = orient;
      setFragOrients();
    }

  void setLibUID(ID_TYPE libUID) {pLibUID = libUID;}
  void setLeftCoord(COORD_TYPE coord) {pLeftFrag.setFiveP(coord);}
  void setLeftFrag(FragmentPosition lfp)
    {
      pLeftFrag = lfp;
      setOrientFromFrags();
    }
  void setRightCoord(COORD_TYPE coord) {pRightFrag.setFiveP(coord);}
  void setRightFrag(FragmentPosition rfp)
    {
      pRightFrag = rfp;
      setOrientFromFrags();
    }
  void setOrientation(PairOrientation_e po)
    {
      pOrient = po;
      setFragOrients();
    }
  void setChromosome(int32 chrom)
    {
      pLeftFrag.setChromosome(chrom);
      pRightFrag.setChromosome(chrom);
    }

  void setFromInterChromosomeString(char * line)
    {
      ID_TYPE leftUID, rightUID;
      int32 leftChrom, rightChrom;
      COORD_TYPE left5p, right5p;
      char leftO[4], rightO[4];
      
      sscanf(line, F_U64 " %d " F_U64 " %s " F_U64 " %d " F_U64 " %s " F_U64,
             &leftUID, &leftChrom, &left5p, leftO,
             &rightUID, &rightChrom, &right5p, rightO,
             &pLibUID);
      pLeftFrag.set(leftUID, left5p, leftO, leftChrom);
      pRightFrag.set(rightUID, right5p, rightO, rightChrom);
      setOrientFromFrags();
    }
  
  void setFromString(char * line, int32 chrom = -1)
    {
      ID_TYPE leftUID, rightUID;
      COORD_TYPE left5p, right5p;
      char orient;
      sscanf(line, "%c " F_U64 " " F_U64 " " F_U64 " " F_U64 " " F_U64,
             &orient, &leftUID, &rightUID, &pLibUID, &left5p, &right5p);
      pLeftFrag.set(leftUID, left5p);
      pRightFrag.set(rightUID, right5p);
      switch(orient)
      {
        case 'I':
          setOrientation(PAIR_INNIE);
          break;
        case 'N':
          setOrientation(PAIR_NORMAL);
          break;
        case 'A':
          setOrientation(PAIR_ANTINORMAL);
          break;
        case 'O':
          setOrientation(PAIR_OUTTIE);
          break;
        default:
          cerr << "Unknown mate pair orientation: " << orient << endl;
          break;
      }
    }
  
  const FragmentPosition & getLeftFrag() const {return pLeftFrag;}
  ID_TYPE getLeftFragUID() const {return pLeftFrag.getUID();}
  COORD_TYPE getLeftCoord() const {return pLeftFrag.getFiveP();}
  SingleOrientation_e getLeftFragOrientation() const
    {
      return pLeftFrag.getOrientation();
    }
  int32 getLeftChromosome() const {return pLeftFrag.getChromosome();}

  const FragmentPosition & getRightFrag() const {return pRightFrag;}
  ID_TYPE getRightFragUID() const {return pRightFrag.getUID();}
  COORD_TYPE getRightCoord() const {return pRightFrag.getFiveP();}
  SingleOrientation_e getRightFragOrientation() const
    {
      return pRightFrag.getOrientation();
    }
  int32 getRightChromosome() const {return pRightFrag.getChromosome();}
  ID_TYPE getLibUID() const {return pLibUID;}
  PairOrientation_e getOrientation() const {return pOrient;}

  int32 getChromosome() const
    {
      return ((pRightFrag.getChromosome() == pLeftFrag.getChromosome()) ?
              pRightFrag.getChromosome() : -1);
    }
      
  bool isWithinDelta(const MatePair & other, COORD_TYPE delta) const
    {
      return(getLeftCoord() + delta >= other.getLeftCoord() &&
             getLeftCoord() - delta <= other.getLeftCoord() &&
             getRightCoord() + delta >= other.getRightCoord() &&
             getRightCoord() - delta <= other.getRightCoord() &&
             getOrientation() == other.getOrientation());
    }
  bool operator<(const MatePair & other) const
    {
      return (getLeftCoord() < other.getLeftCoord());
    }
  
  bool intersects(COORD_TYPE left, COORD_TYPE right) const
    {
      return(getLeftCoord() < right && left < getRightCoord());
    }
  bool intersects(const MatePair & other) const
    {
      return intersects(other.getLeftCoord(), other.getRightCoord());
    }
  
  bool spans(COORD_TYPE left, COORD_TYPE right) const
    {
      return(getLeftCoord() < left && right < getRightCoord());
    }
  bool spans(const MatePair & other) const
    {
      return spans(other.getLeftCoord(), other.getRightCoord());
    }

  bool getIntersection(COORD_TYPE left, COORD_TYPE right,
                       COORD_TYPE & leftI, COORD_TYPE & rightI) const
    {
      if(!intersects(left, right)) return false;
      leftI = max(left, getLeftCoord());
      rightI = min(right, getRightCoord());
      return true;
    }
  bool getIntersection(const MatePair & other,
                       COORD_TYPE & leftI, COORD_TYPE & rightI) const
    {
      return getIntersection(other.getLeftCoord(),
                             other.getRightCoord(),
                             leftI, rightI);
    }

  void printATA(ostream & os, char * assembly, int chromosome,
                char * pKey, int pID, int id) const
    {
      switch(pOrient)
      {
        case PAIR_INNIE:
          os << "F mi mi";
          break;
        case PAIR_OUTTIE:
          os << "F mo mo";
          break;
        case PAIR_NORMAL:
          os << "F mn mn";
          break;
        case PAIR_ANTINORMAL:
          os << "F ma ma";
          break;
        default:
          return;
          break;
      }
      if(pID > -1)
        os << assembly << "-" << pID << "-" << id
           << " " << pKey << assembly << "-" << pID;
      else
        os << assembly << "-" << id
           << " .";

      os << " " << assembly << ":" << chromosome
         << " " << getLeftCoord()
         << " " << getRightCoord() - getLeftCoord()
         << " ."
         << " ."
         << " ."
         << " ."
         << " > /leftUID=" << pLeftFrag.getUID()
         << " > /rightUID=" << pRightFrag.getUID()
         << " /libUID=" << pLibUID << endl;
    }
  
  friend ostream & operator<<(ostream & os, const MatePair & mp)
    {
      switch(mp.pOrient)
      {
        case PAIR_INNIE:
          os << "I ";
          break;
        case PAIR_OUTTIE:
          os << "O ";
          break;
        case PAIR_NORMAL:
          os << "N ";
          break;
        case PAIR_ANTINORMAL:
          os << "A ";
          break;
        default:
          os << "? ";
          break;
      }
      os << mp.pLeftFrag.getUID() << " " << mp.pRightFrag.getUID() << " "
         << mp.pLibUID << " "
         << mp.pLeftFrag.getFiveP() << " " << mp.pRightFrag.getFiveP();
      return os;
    }

private:
  void setFragOrients()
    {
      switch(pOrient)
      {
        case PAIR_INNIE:
          pLeftFrag.setOrientation(SINGLE_A_B);
          pRightFrag.setOrientation(SINGLE_B_A);
          break;
        case PAIR_OUTTIE:
          pLeftFrag.setOrientation(SINGLE_B_A);
          pRightFrag.setOrientation(SINGLE_A_B);
          break;
        case PAIR_NORMAL:
          pLeftFrag.setOrientation(SINGLE_A_B);
          pRightFrag.setOrientation(SINGLE_A_B);
          break;
        case PAIR_ANTINORMAL:
          pLeftFrag.setOrientation(SINGLE_B_A);
          pRightFrag.setOrientation(SINGLE_B_A);
          break;
        default:
          assert(0);
          break;
      }
    }
  
  void setOrientFromFrags()
    {
      if(pLeftFrag.pointsRight())
      {
        if(pRightFrag.pointsRight())
          pOrient = PAIR_NORMAL;
        else
          pOrient = PAIR_INNIE;
      }
      else
      {
        if(pRightFrag.pointsRight())
          pOrient = PAIR_OUTTIE;
        else
          pOrient = PAIR_ANTINORMAL;
      }
    }
  FragmentPosition pLeftFrag;
  FragmentPosition pRightFrag;
  PairOrientation_e pOrient;
  ID_TYPE pLibUID;
};


#endif
