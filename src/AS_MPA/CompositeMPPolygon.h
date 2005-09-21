
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
/* $Id: CompositeMPPolygon.h,v 1.6 2005-09-21 20:13:07 catmandew Exp $ */
#ifndef COMPOSITEMPPOLYGON_H
#define COMPOSITEMPPOLYGON_H

#include "MatePairPolygon.h"
#include "CompositePolygon.h"
#include "MatePairGroup.h"

template <class UnitType>
class CompositeMPPolygon :
  public CompositePolygon<UnitType>,
  public MatePairGroup
{
public:
  CompositeMPPolygon() :
    CompositePolygon<UnitType>(), MatePairGroup()
    {
    }
  CompositeMPPolygon(const MatePairPolygon<UnitType> & mpp) :
    CompositePolygon<UnitType>(mpp), MatePairGroup(mpp)
    {
    }

  void printLeftBPATA(ostream & os,
                   char * assembly,
                   char * sequenceID,
                   int relativeID,
                   char * parentKey,
                   char * bpKey) const
    {
      os << "F " << bpKey << " " << bpKey << assembly << "-" << relativeID
         << " " << parentKey << assembly << "-" << relativeID
         << " " << assembly << ":" << sequenceID
         << " " << this->getMinX()
         << " " << this->getMaxX() - this->getMinX()
         << " . . . ." << endl;
    }
  
  void printRightBPATA(ostream & os,
                    char * assembly,
                    char * sequenceID,
                    int relativeID,
                    char * parentKey,
                    char * bpKey) const
    {
      os << "F " << bpKey << " " << bpKey << assembly << "-" << relativeID
         << " " << parentKey << assembly << "-" << relativeID
         << " " << assembly << ":" << sequenceID
         << " " << this->getMinY()
         << " " << this->getMaxY() - this->getMinY()
         << " . . . ." << endl;
    }
  
  void printATA(ostream & os,
                char * assembly,
                char * sequenceID,
                int relativeID,
                bool printMatePairs) const
    {
      /*
        There are four components to each unsatisfied mate pair:
        1. interval of problem
        2. bounds on left breakpoint
        3. bounds on right breakpoint
        4. mate pairs that contribute

        For deletions, 2 & 3 don't apply
      */
      char key[3];
      switch(getMPIndex())
      {
        case MPI_STRETCHED:
          sprintf(key, "ii");
          // print left end & maximum insertion length
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMinX()
             << " " << this->getMaxY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBPATA(os, assembly, sequenceID, relativeID, key, "il");
          printRightBPATA(os, assembly, sequenceID, relativeID, key, "ir");
          break;
        case MPI_COMPRESSED:
          sprintf(key, "dd");
          // print left end & maximum deletion length
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMinX()  // leftmost possible breakpoint
             << " " << this->getMaxX() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          os << "F " << key << "s " << key
             << assembly << "s-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMinY()  // smallest possible size
             << " " << this->getMaxY()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << endl;
          break;
        case MPI_NORMAL:
          sprintf(key, "nn");
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMinX()
             << " " << this->getMinY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBPATA(os, assembly, sequenceID, relativeID, key, "nl");
          printRightBPATA(os, assembly, sequenceID, relativeID, key, "nr");
          break;
        case MPI_ANTINORMAL:
          sprintf(key, "aa");
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMaxX()
             << " " << this->getMaxY() - this->getMaxX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBPATA(os, assembly, sequenceID, relativeID, key, "al");
          printRightBPATA(os, assembly, sequenceID, relativeID, key, "ar");
          break;
        case MPI_OUTTIE:
          sprintf(key, "oo");
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMaxX()
             << " " << this->getMinY() - this->getMaxX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBPATA(os, assembly, sequenceID, relativeID, key, "ol");
          printRightBPATA(os, assembly, sequenceID, relativeID, key, "or");
          break;
        case MPI_INVERSION:
          sprintf(key, "vv");
          // print left end & maximum insertion length
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMinX()
             << " " << this->getMaxY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBPATA(os, assembly, sequenceID, relativeID, key, "ol");
          printRightBPATA(os, assembly, sequenceID, relativeID, key, "or");
          break;
        case MPI_TRANSPOSITION:
          sprintf(key, "tt");
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMinX()
             << " " << this->getMaxY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBPATA(os, assembly, sequenceID, relativeID, key, "tl");
          printRightBPATA(os, assembly, sequenceID, relativeID, key, "tr");
          break;
        case MPI_SATISFIED:
          sprintf(key, "ss");
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMinX()
             << " " << this->getMaxY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          break;
        case MPI_UNKNOWN:
        default:
          /*
          sprintf(key, "un");
          os << "F " << key << " " << key
             << assembly << "-" << relativeID
             << " ."
             << " " << assembly << ":" << sequenceID
             << " " << this->getMinX()
             << " " << this->getMaxY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          */
          break;
      }

      // print matepairs
      if(printMatePairs)
      {
        for(unsigned int i = 0; i < getNumMPs(); i++)
          pmps[i].printATA(os, assembly, sequenceID, key, relativeID, i);
      }
    }
  
  void printForGnuplot(ostream & os, CompressedRepresentation_e cr) const
    {
      for(unsigned int i = 0; i < getNumMPs(); i++)
        os << "# " << pmps[i] << endl;

      if(isCompressed())
      {
        switch(cr)
        {
          case CR_COMPATIBLE:
            for(unsigned int i = 0; i < this->pPolys.size(); i++)
            {
              Polygon<UnitType> tempP;
              UnitType minX = this->getMinX();
              UnitType maxX = this->getMaxX();
              UnitType minY = this->getMinX();
              UnitType maxY = this->getMaxX();
              tempP.append(minX, minY);
              tempP.append(minX, maxY);
              tempP.append(maxX, maxY);
              tempP.append(maxX, minY);
              tempP.printForGnuplot(os);
            }
            break;
          case CR_NATIVE:
            for(unsigned int i = 0; i < this->pPolys.size(); i++)
            {
              Polygon<UnitType> tempP;
              UnitType minX = this->getMinX();
              UnitType maxX = this->getMaxX();
              UnitType minY = this->getMinY();
              UnitType maxY = this->getMaxY();
              tempP.append(minX, minY);
              tempP.append(minX, maxY);
              tempP.append(maxX, maxY);
              tempP.append(maxX, minY);
              tempP.printForGnuplot(os);
            }
            break;
        }
      }
      else
      {
        for(unsigned int i = 0; i < this->pPolys.size(); i++)
        {
          this->pPolys[i].printForGnuplot(os);
        }
      }
    }
  void printSummary(ostream & os) const
    {
      if(this->getMPIndex() == MPI_COMPRESSED)
      {
        os << this->getMinX() << "\t"
           << this->getMaxX() - this->getMinX() << "\t"
           << this->getMinY() << "\t"
           << this->getMaxY() << "\t"
           << this->getNumMPs() << endl;
      }
      else
      {
        os << this->getMinX() << "\t"
           << this->getMaxX() - this->getMinX() << "\t"
           << this->getMinY() << "\t"
           << this->getMaxY() - this->getMinY() << "\t"
           << ((this->getMaxY() + this->getMinY()) -
               (this->getMaxX() + this->getMinX())) / 2 << "\t"
           << this->getNumMPs() << endl;
      }
    }

private:
};

#endif //#ifndef COMPOSITEMPPOLYGON_H
