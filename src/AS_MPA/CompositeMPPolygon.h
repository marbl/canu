
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
/* $Id: CompositeMPPolygon.h,v 1.4 2005-03-22 19:48:56 jason_miller Exp $ */
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

  void printLeftBP(ostream & os,
                   char * assembly,
                   int chromosome,
                   int id,
                   char * parentKey,
                   char * bpKey) const
    {
      os << "F " << bpKey << " " << bpKey << assembly << "-" << id
         << " " << parentKey << assembly << "-" << id
         << " " << assembly << ":" << chromosome
         << " " << this->getMinX()
         << " " << this->getMaxX() - this->getMinX()
         << " . . . ." << endl;
    }
  
  void printRightBP(ostream & os,
                    char * assembly,
                    int chromosome,
                    int id,
                    char * parentKey,
                    char * bpKey) const
    {
      os << "F " << bpKey << " " << bpKey << assembly << "-" << id
         << " " << parentKey << assembly << "-" << id
         << " " << assembly << ":" << chromosome
         << " " << this->getMinY()
         << " " << this->getMaxY() - this->getMinY()
         << " . . . ." << endl;
    }
  
  void printATA(ostream & os,
                char * assembly,
                int chromosome,
                int id,
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
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
             << " " << this->getMinX()
             << " " << this->getMaxY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBP(os, assembly, chromosome, id, key, "il");
          printRightBP(os, assembly, chromosome, id, key, "ir");
          break;
        case MPI_COMPRESSED:
          sprintf(key, "dd");
          // print left end & maximum deletion length
          os << "F " << key << " " << key
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
             << " " << this->getMinX()  // leftmost possible breakpoint
             << " " << this->getMaxX() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          os << "F " << key << "s " << key
             << assembly << "s-" << id
             << " ."
             << " " << assembly << ":" << chromosome
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
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
             << " " << this->getMinX()
             << " " << this->getMinY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBP(os, assembly, chromosome, id, key, "nl");
          printRightBP(os, assembly, chromosome, id, key, "nr");
          break;
        case MPI_ANTINORMAL:
          sprintf(key, "aa");
          os << "F " << key << " " << key
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
             << " " << this->getMaxX()
             << " " << this->getMaxY() - this->getMaxX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBP(os, assembly, chromosome, id, key, "al");
          printRightBP(os, assembly, chromosome, id, key, "ar");
          break;
        case MPI_OUTTIE:
          sprintf(key, "oo");
          os << "F " << key << " " << key
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
             << " " << this->getMaxX()
             << " " << this->getMinY() - this->getMaxX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBP(os, assembly, chromosome, id, key, "ol");
          printRightBP(os, assembly, chromosome, id, key, "or");
          break;
        case MPI_INVERSION:
          sprintf(key, "vv");
          // print left end & maximum insertion length
          os << "F " << key << " " << key
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
             << " " << this->getMinX()
             << " " << this->getMaxY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBP(os, assembly, chromosome, id, key, "ol");
          printRightBP(os, assembly, chromosome, id, key, "or");
          break;
        case MPI_TRANSPOSITION:
          sprintf(key, "tt");
          os << "F " << key << " " << key
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
             << " " << this->getMinX()
             << " " << this->getMaxY() - this->getMinX()
             << " -1"
             << " 0"
             << " ."
             << " ."
             << " > /weight=" << getNumMPs() << endl;
          printLeftBP(os, assembly, chromosome, id, key, "tl");
          printRightBP(os, assembly, chromosome, id, key, "tr");
          break;
        case MPI_SATISFIED:
          sprintf(key, "ss");
          os << "F " << key << " " << key
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
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
             << assembly << "-" << id
             << " ."
             << " " << assembly << ":" << chromosome
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
          pmps[i].printATA(os, assembly, chromosome, key, id, i);
      }
    }
  void printForGnuplot(ostream & os) const
    {
      for(unsigned int i = 0; i < getNumMPs(); i++)
        os << "# " << pmps[i] << endl;

      if(isCompressed())
      {
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
      }
      else
      {
        for(unsigned int i = 0; i < this->pPolys.size(); i++)
        {
          this->pPolys[i].printForGnuplot(os);
        }
      }
    }
private:
};

#endif //#ifndef COMPOSITEMPPOLYGON_H
