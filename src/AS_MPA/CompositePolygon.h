
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
/* $Id: CompositePolygon.h,v 1.5 2008-06-27 06:29:16 brianwalenz Exp $ */
#ifndef COMPOSITEPOLYGON_H
#define COMPOSITEPOLYGON_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "Polygon.h"
#include "Interval.h"

template <class UnitType>
class CompositePolygon
{
public:
  CompositePolygon()
    {
    }
  CompositePolygon(const Polygon<UnitType> & poly)
    {
      pPolys.push_back(poly);
    }

  unsigned int size() const {return pPolys.size();}

  const Polygon<UnitType> & operator[](int index) const
    {
      return pPolys[index];
    }
  Polygon<UnitType> & operator[](int index)
    {
      return pPolys[index];
    }

  void getMBR(UnitType & minX, UnitType & maxX,
              UnitType & minY, UnitType & maxY) const
    {
      UnitType nx, xx, ny, xy;
      if(size() <= 0)
      {
        minX = maxX = minY = maxY = -1;
      }
      pPolys[0].getMBR(minX, maxX, minY, maxY);
      for(unsigned int i = 1; i < size(); i++)
      {
        pPolys[i].getMBR(nx, xx, ny, xy);
        minX = (minX < nx) ? minX : nx;
        maxX = (maxX < xx) ? xx : maxX;
        minY = (minY < ny) ? minY : ny;
        maxY = (maxY < xy) ? xy : maxY;
      }
    }

  void append(const CompositePolygon<UnitType> & cp)
    {
      for(unsigned int i = 0; i < cp.size(); i++)
        append(cp[i]);
    }
  void append(const vector<Point<UnitType> > & points)
    {
      Polygon<UnitType> poly(points);
      append(poly);
    }

  void append(const Polygon<UnitType> & poly)
    {
      pPolys.push_back(poly);
    }
  void reset()
    {
      pPolys.clear();
    }

  void setPolygons(const vector<Polygon<UnitType> > & polys)
    {
      pPolys.clear();
      pPolys = polys;
    }
  void setPolygons(const CompositePolygon<UnitType> & cp)
    {
      setPolygons(cp.pPolys);
    }

  void printForGnuplot(ostream & os) const
    {
      for(unsigned int i = 0; i < size(); i++)
        pPolys[i].printForGnuplot(os);
    }
  void print(ostream & os) const
    {
      for(unsigned int i = 0; i < size(); i++)
        pPolys[i].print(os);
    }
  void rotateByRadians(double radians)
    {
      for(unsigned int i = 0; i < size(); i++)
        pPolys[i].rotateByRadians(radians);
    }
  void rotateByDegrees(double degrees)
    {
      rotateByRadians(degrees * M_PI / 180.);
    }

  UnitType getMinX() const
    {
      if(pPolys.size() == 0) return -1;
      UnitType retVal = pPolys[0].getMinX();
      for(unsigned int i = 1; i < pPolys.size(); i++)
      {
        UnitType thisVal = pPolys[i].getMinX();
        retVal = (retVal < thisVal) ? retVal : thisVal;
      }
      return retVal;
    }
  UnitType getLeftMostX() const
    {
      return getMinX();
    }
  UnitType getMaxX() const
    {
      if(pPolys.size() == 0) return -1;
      UnitType retVal = pPolys[0].getMaxX();
      for(unsigned int i = 1; i < pPolys.size(); i++)
      {
        UnitType thisVal = pPolys[i].getMaxX();
        retVal = (retVal > thisVal) ? retVal : thisVal;
      }
      return retVal;
    }
  UnitType getMinY() const
    {
      if(pPolys.size() == 0) return -1;
      UnitType retVal = pPolys[0].getMinY();
      for(unsigned int i = 1; i < pPolys.size(); i++)
      {
        UnitType thisVal = pPolys[i].getMinY();
        retVal = (retVal < thisVal) ? retVal : thisVal;
      }
      return retVal;
    }
  UnitType getMaxY() const
    {
      if(pPolys.size() == 0) return -1;
      UnitType retVal = pPolys[0].getMaxY();
      for(unsigned int i = 1; i < pPolys.size(); i++)
      {
        UnitType thisVal = pPolys[i].getMaxY();
        retVal = (retVal > thisVal) ? retVal : thisVal;
      }
      return retVal;
    }
  friend ostream & operator<<(ostream & os,
                              const CompositePolygon<UnitType> & cpoly)
    {
      os << cpoly.size() << " polygons:\n";
      for(unsigned int i = 0; i < cpoly.size(); i++)
        os << cpoly[i] << endl;
      return os;
    }
protected:
  vector<Polygon<UnitType> > pPolys;
};

#endif // #ifndef COMPOSITEPOLYGON_H
