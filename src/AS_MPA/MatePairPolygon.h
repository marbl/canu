
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
/* $Id: MatePairPolygon.h,v 1.1.1.1 2004-04-14 13:52:04 catmandew Exp $ */
#ifndef MATEPAIRPOLYGON_H
#define MATEPAIRPOLYGON_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>

#include "Polygon.h"
#include "CloneLibrary.h"
#include "MatePairGroup.h"

template <class UnitType>
class MatePairPolygon :
  public Polygon<UnitType>,
  public MatePairGroup
{
public:
  MatePairPolygon()
    {
    }
  MatePairPolygon(const vector<Point<UnitType> > & pts, const MatePair & mp) :
    Polygon<UnitType>(pts), MatePairGroup(mp)
    {
    }
  MatePairPolygon(const MatePairPolygon<UnitType> & other) :
    Polygon<UnitType>(other), MatePairGroup(other)
    {
    }
  MatePairPolygon(const MatePair & mp, const CloneLibrary & cl, double n) :
    MatePairGroup(mp)
    {
      Point<UnitType> p;
      UnitType large = (UnitType) (cl.getMean() + n * cl.getStddev());
      UnitType small = (UnitType) (cl.getMean() - n * cl.getStddev());
      switch(mp.getOrientation())
      {
        case PAIR_INNIE:
          if(mp.getRightCoord() > mp.getLeftCoord() + large)
          {
            // stretched
          /* Four points:p
             0. left5p, right5p - (mean + n * stddev)
             1. left5p, right5p - (mean - n * stddev)
             2. left5p + (mean - n * stddev), right5p
             3. left5p + (mean + n * stddev), right5p
             
             |         2     3
             v
                 
             -  1


             -  0
                --->   (     )
                left frag
           */
            p.setX(mp.getLeftCoord());
            p.setY(mp.getRightCoord() - large);
            append(p);
            p.setY(mp.getRightCoord() - small);
            append(p);
            p.setX(mp.getLeftCoord() + small);
            p.setY(mp.getRightCoord());
            append(p);
            p.setX(mp.getLeftCoord() + large);
            append(p);
            setMPIndex(MPI_STRETCHED);
          }
          else if(mp.getRightCoord() < mp.getLeftCoord() + small)
          {
            // compressed
          /* Four points:p
           */
            /*
            p.setX(mp.getLeftCoord());
            p.setY(small);
            append(p);
            p.setY(large);
            append(p);
            p.setX(mp.getRightCoord());
            append(p);
            p.setY(small);
            append(p);
            */
            p.setX(mp.getLeftCoord());
            p.setY(mp.getLeftCoord() + small - mp.getRightCoord());
            append(p);
            p.setY(mp.getLeftCoord() + large - mp.getRightCoord());
            append(p);
            p.setX(mp.getRightCoord());
            append(p);
            p.setY(mp.getLeftCoord() + small - mp.getRightCoord());
            append(p);
            setMPIndex(MPI_COMPRESSED);
          }
          else
          {
            // satisfied
            p.setX(mp.getLeftCoord());
            p.setY(mp.getLeftCoord());
            append(p);
            p.setY(mp.getRightCoord());
            append(p);
            p.setX(mp.getRightCoord());
            append(p);
            p.setY(mp.getLeftCoord());
            append(p);
            setMPIndex(MPI_SATISFIED);
          }
          break;
        case PAIR_NORMAL:
          /* Four points:p
             0. left5p, right5p + (mean - n * stddev)
             1. left5p, right5p + (mean + n * stddev)
             2. left5p + mean + n * stddev, right5p
             3. left5p + mean - n * stddev, right5p
             
            -   1


            -   0

            ^
            |          3     2
                -->    (     )
                left frag
          */
          p.setX(mp.getLeftCoord());
          p.setY(mp.getRightCoord() + small);
          append(p);
          p.setY(mp.getRightCoord() + large);
          append(p);
          p.setX(mp.getLeftCoord() + large);
          p.setY(mp.getRightCoord());
          append(p);
          p.setX(mp.getLeftCoord() + small);
          append(p);
          setMPIndex(MPI_NORMAL);
          break;
        case PAIR_ANTINORMAL:
          /* Four points:p
             0. left5p, right5p + (mean - n * stddev)
             1. left5p, right5p + (mean + n * stddev)
             2. left5p + mean + n * stddev, right5p
             3. left5p + mean - n * stddev, right5p
             
            |  0     1
            v

            -              2


            -              3
               (     )   <--
                         left frag
          */
          p.setX(mp.getLeftCoord() - large);
          p.setY(mp.getRightCoord());
          append(p);
          p.setX(mp.getLeftCoord() - small);
          append(p);
          p.setX(mp.getLeftCoord());
          p.setY(mp.getRightCoord() - small);
          append(p);
          p.setY(mp.getRightCoord() - large);
          append(p);
          setMPIndex(MPI_ANTINORMAL);
          break;
        case PAIR_OUTTIE:
          /* Four points:p
             0. left5p - (mean + n * stddev), right5p
             1. left5p, right5p + (mean + n * stddev)
             2. left5p, right5p + (mean - n * stddev)
             3. left5p - (mean - n * stddev), right5p
             
            -             1


            -             2

            ^
            |   0    3
                (    )   <--
                         left frag
          */
          p.setX(mp.getLeftCoord() - large);
          p.setY(mp.getRightCoord());
          append(p);
          p.setX(mp.getLeftCoord());
          p.setY(mp.getRightCoord() + large);
          append(p);
          p.setY(mp.getRightCoord() + small);
          append(p);
          p.setX(mp.getLeftCoord() - small);
          p.setY(mp.getRightCoord());
          append(p);
          setMPIndex(MPI_OUTTIE);
          break;
        default:
          cerr << "Unknown mate pair orientation: "
               << mp.getOrientation() << endl;
          assert(0);
          break;
      }
    }

  void printForGnuplot(ostream & os) const
    {
      if(ppts.size() == 0) return;

      for(int i = 0; i < getNumMPs(); i++)
        os << "# " << pmps[i] << endl;

      if(isCompressed())
      {
        UnitType minX = getMinX();
        UnitType maxX = getMaxX();
        os << minX << " " << minX << endl;
        os << minX << " " << maxX << endl;
        os << maxX << " " << maxX << endl;
        os << maxX << " " << minX << endl;
        os << minX << " " << minX << endl << endl;
      }
      else
      {
        for(unsigned int i = 0; i < ppts.size(); i++)
          os << ppts[i].getX() << " " << ppts[i].getY() << endl;
        os << ppts[0].getX() << " " << ppts[0].getY() << endl << endl;
      }
    }
  
  void print(ostream & os) const
    {
      if(ppts.size() == 0) return;

      for(int i = 0; i < getNumMPs(); i++)
        os << "# " << pmps[i] << endl;
      
      for(unsigned int i = 0; i < poly.ppts.size(); i++)
        os << " " << poly.ppts[i];
      os << "\n" << poly.getNumMPs() << " mate pairs:\n";
      for(int i = 0; i < poly.getNumMPs(); i++)
        os << poly.pmps[i] << endl;
      return os;
    }
  
private:
};


#endif // MATEPAIRPOLYGON_H
