
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
/* $Id: Polygon.h,v 1.5 2006-10-08 08:47:39 brianwalenz Exp $ */
#ifndef POLYGON_H
#define POLYGON_H

#include <iostream>
#include <vector>
#include <cassert>
#include <inttypes.h>
#include <stdint.h>

#include "Point.h"
#include "Interval.h"

#define MAX_POLYGON_VERTICES 100

template <class UnitType>
class Polygon
{
public:
  Polygon() {pSortX = true; pRangesSet = false;}
  Polygon(const vector<Point<UnitType> > & pts)
    {
      typename vector<Point<UnitType> >::const_iterator iter;
      for(iter = pts.begin(); iter != pts.end(); pts++)
        ppts.push_back(*iter);
      pSortX = true;
      pRangesSet = false;
    }
  Polygon(const Polygon<UnitType> & other)
    {
      for(unsigned int i = 0; i < other.ppts.size(); i++)
        ppts.push_back(other.ppts[i]);
      pSortX = true;
      pRangesSet = false;
    }
  
  unsigned int size() const
    {
      return ppts.size();
    }
  const Point<UnitType> & operator[](int index) const
    {
      return ppts[index];
    }
  Point<UnitType> & operator[](int index)
    {
      return ppts[index];
    }
  unsigned int cwIndex(int i) const
    {
      return((i + 1) % ppts.size());
    }
  unsigned int ccwIndex(int i) const
    {
      return((i - 1 < 0) ? (ppts.size() - 1) : (i - 1));
    }

  unsigned int getLeftMostIndex() const
    {
      if(pRangesSet) return pLeftMostIndex;
      
      int leftMostIndex = 0;
      for(unsigned int i = 1; i < ppts.size(); i++)
      {
        if(ppts[i].getX() < ppts[leftMostIndex].getX() ||
           (ppts[i].getX() == ppts[leftMostIndex].getX() &&
            ppts[i].getY() < ppts[leftMostIndex].getY()))
          leftMostIndex = i;
        return leftMostIndex;
      }
    }

  UnitType getLeftMostX() const
    {
      if(pRangesSet) return pMinX;
      UnitType retVal = ppts[0].getX();
      for(unsigned int i = 0; i < ppts.size(); i++)
        retVal = (retVal > ppts[i].getX()) ? ppts[i].getX() : retVal;
      return retVal;
    }
  UnitType getMinX() const
    {
      return getLeftMostX();
    }
  UnitType getMaxX() const
    {
      if(pRangesSet) return pMaxX;
      UnitType retVal = ppts[0].getX();
      for(unsigned int i = 0; i < ppts.size(); i++)
        retVal = (retVal < ppts[i].getX()) ? ppts[i].getX() : retVal;
      return retVal;
    }
  UnitType getMinY() const
    {
      if(pRangesSet) return pMinY;
      UnitType retVal = ppts[0].getY();
      for(unsigned int i = 0; i < ppts.size(); i++)
        retVal = (retVal > ppts[i].getY()) ? ppts[i].getY() : retVal;
      return retVal;
    }
  UnitType getMaxY() const
    {
      if(pRangesSet) return pMaxY;
      UnitType retVal = ppts[0].getY();
      for(unsigned int i = 0; i < ppts.size(); i++)
        retVal = (retVal < ppts[i].getY()) ? ppts[i].getY() : retVal;
      return retVal;
    }
  void getMBR(UnitType & minX, UnitType & maxX,
              UnitType & minY, UnitType & maxY) const
    {
      if(pRangesSet)
      {
        minX = pMinX;
        maxX = pMaxX;
        minY = pMinY;
        maxY = pMaxY;
        return;
      }
      minX = maxX = ppts[0].getX();
      minY = maxY = ppts[0].getY();
      for(unsigned int i = 1; i < ppts.size(); i++)
      {
        minX = (minX < ppts[i].getX()) ? minX : ppts[i].getX();
        maxX = (maxX > ppts[i].getX()) ? maxX : ppts[i].getX();
        minY = (minY < ppts[i].getY()) ? minY : ppts[i].getY();
        maxY = (maxY > ppts[i].getY()) ? maxY : ppts[i].getY();
      }
    }
  
  int getLeftMostIndex()
    {
      setRanges();
      return pLeftMostIndex;
    }
  UnitType getLeftMostX()
    {
      setRanges();
      return pMinX;
    }
  UnitType getMinX()
    {
      return getLeftMostX();
    }
  UnitType getMaxX()
    {
      setRanges();
      return pMaxX;
    }
  UnitType getMinY()
    {
      setRanges();
      return pMinY;
    }
  UnitType getMaxY()
    {
      setRanges();
      return pMaxY;
    }
  void getMBR(UnitType & minX, UnitType & maxX,
              UnitType & minY, UnitType & maxY)
    {
      setRanges();
      minX = pMinX;
      maxX = pMaxX;
      minY = pMinY;
      maxY = pMaxY;
    }
  
  double getSegmentAngle(int i) const
    {
      Point<double> p1(ppts[i].getX(), ppts[i].getY());
      Point<double> p2(ppts[cwIndex(i)].getX(), ppts[cwIndex(i)].getY());
      p2.subtract(p1);
      p2.normalize();
      
      return M_PI + ((p2.getX() >= 0) ? -acos(-p2.getY()) : acos(-p2.getY()));
    }
  
  void getSegmentAngles(vector<double> & angles) const
    {
      for(unsigned int i = 0; i < ppts.size(); i++)
      {
        double angle = getSegmentAngle(i);
        angles.push_back(angle);
      }
    }
  
  bool append(const Point<UnitType> & point)
    {
      // if it's the same point as the last one, just skip it
      if(ppts.size() > 0 &&
         (point == ppts[0] || point == ppts[ppts.size()-1]))
      {
        cerr << "Attempted to append duplicate vertex to polygon (0)!\n";
        return false;
      }
      
      if(ppts.size() > 1 &&
         (point == ppts[1] || point == ppts[ppts.size()-2]))
      {
        cerr << "Attempted to append duplicate vertex to polygon (1)!\n";
        return false;
      }
      
      if(!isConvexWithPoint(point))
      {
        cerr << "Polygon not convex with new point!\n";
        cerr << "polygon: " << *this << endl;
        cerr << "omitting new point: " << point << endl;
        return false;
      }
      else
        ppts.push_back(point);
      pRangesSet = false;
      setRanges();
      return true;
    }
  
  bool append(UnitType x, UnitType y)
    {
      Point<UnitType> p(x,y);
      return append(p);
    }
  
  bool isConvexWithPoint(UnitType x, UnitType y) const
    {
      Point<UnitType> p(x,y);
      return isConvexWithPoint(p);
    }
  bool isConvexWithPoint(const Point<UnitType> & point) const
    {
      if(ppts.size() < 2) return true;
      
      // point must be to the right (clockwise) from line
      // between previous two points
      if(point.getLineSide(ppts[ppts.size()-2], ppts[ppts.size()-1]) ==
         TP_COUNTERCW)
        return false;
      
      // and 2nd point must be to the right from line
      // between point and ppts[0]
      if(ppts.size() > 2 && ppts[1].getLineSide(point, ppts[0]) ==
         TP_COUNTERCW)
        return false;
      
      return true;
    }

  void setPoints(const vector<Point<UnitType> > & points)
    {
      ppts.clear();
      for(unsigned int i = 0; i < points.size(); i++)
      {
        if(!append(points[i]))
          return false;
      }
      return true;
    }
  void setPoints(const Polygon<UnitType> & other)
    {
      setPoints(other.ppts);
    }
  void reset()
    {
      pSortX = true;
      pRangesSet = false;
      ppts.clear();
    }
  
  void printForGnuplot(ostream & os) const
    {
      if(ppts.size() == 0) return;
      
      for(unsigned int i = 0; i < ppts.size(); i++)
      {
        os << ppts[i].getX() << " " << ppts[i].getY() << endl;
      }
      os << ppts[0].getX() << " " << ppts[0].getY() << endl << endl;
    }
  
  void print(ostream & os) const
    {
      if(ppts.size() == 0) return;
      
      for(unsigned int i = 0; i < ppts.size(); i++)
      {
        os << ppts[i].getX() << "," << ppts[i].getY() << " ";
      }
      os << endl;
    }
  
  void rotateByRadians(double radians)
    {
      for(unsigned int i = 0; i < ppts.size(); i++)
        ppts[i].rotateByRadians(radians);
      pRangesSet = false;
    }
  
  void rotateByDegrees(double degrees)
    {
      rotateByRadians(degrees * M_PI / 180.);
    }

#ifdef NEVER
  bool intersects(const Polygon<UnitType> & other) const
    {
      PolygonIntersector<UnitType> intersector;
      return intersector.polygonsIntersect(*this, other);
    }
  bool intersects(const Interval<int, UnitType> & xi,
                  const Interval<int, UnitType> & yi) const
    {
      Polygon<UnitType> other;
      other.append(xi.getMin(), yi.getMin());
      other.append(xi.getMin(), yi.getMax());
      other.append(xi.getMax(), yi.getMax());
      other.append(xi.getMax(), yi.getMin());
      return intersects(other);
    }
#endif
  
  friend ostream & operator<<(ostream & os, const Polygon<UnitType> & poly)
    {
      os << poly.ppts.size();
      for(unsigned int i = 0; i < poly.ppts.size(); i++)
        os << " " << poly.ppts[i];
      return os;
    }
  
  void setXSort() {pSortX = true;}
  void setYSort() {pSortX = false;}

  bool operator<(Polygon<UnitType> & other)
    {
      setRanges();
      other.setRanges();
      if(pSortX)
      {
        if(pMinX < other.pMinX) return -1;
        else if(pMinX > other.pMinX) return 1;
        else return 0;
      }
      else
      {
        if(pMinY < other.pMinY) return -1;
        else if(pMinY > other.pMinY) return 1;
        else return 0;
      }
    }
  
  
protected:
  void setRanges()
    {
      if(!pRangesSet)
      {
        pLeftMostIndex = 0;
        pMinX = pMaxX = ppts[0].getX();
        pMinY = pMaxY = ppts[0].getY();
        for(unsigned int i = 1; i < ppts.size(); i++)
        {
          if(ppts[i].getX() < ppts[pLeftMostIndex].getX() ||
             (ppts[i].getX() == ppts[pLeftMostIndex].getX() &&
              ppts[i].getY() < ppts[pLeftMostIndex].getY()))
            pLeftMostIndex = i;
          pMinX = (pMinX < ppts[i].getX()) ? pMinX : ppts[i].getX();
          pMaxX = (pMaxX > ppts[i].getX()) ? pMaxX : ppts[i].getX();
          pMinY = (pMinY < ppts[i].getY()) ? pMinY : ppts[i].getY();
          pMaxY = (pMaxY > ppts[i].getY()) ? pMaxY : ppts[i].getY();
        }
        pRangesSet = true;
      }
    }
  
  // points must be in clockwise order
  vector<Point<UnitType> > ppts;
  bool pSortX;
  bool pRangesSet;
  UnitType pMinX, pMaxX, pMinY, pMaxY;
  int pLeftMostIndex;
};


#endif

