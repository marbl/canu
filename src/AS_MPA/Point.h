
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
/* $Id: Point.h,v 1.3 2005-03-22 19:05:28 jason_miller Exp $ */
#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <iostream>
#include <list>
#include <vector>

using namespace std;

#include "MPTypes.h"

typedef enum
{
  TP_COUNTERCW = -1,
  TP_COLINEAR,
  TP_CLOCKWISE
} ThreePointOrientation;

template <class UnitType>
class Point
{
public:
  Point()
    {
      px = py = 0;
    }
  
  Point(const UnitType & x, const UnitType & y)
    {
      px = x;
      py = y;
    }
  Point(const Point<UnitType> & point)
    {
      px = point.px;
      py = point.py;
    }

  UnitType getX() const {return px;}
  UnitType getY() const {return py;}

  void setX(const UnitType & x) {px = x;}
  void setY(const UnitType & y) {py = y;}
  void set(const UnitType & x, const UnitType & y) {px = x; py = y;}
  void set(const Point<UnitType> & p)
    {
      set(p.px, p.py);
    }

  double distanceFrom(const Point<UnitType> & point) const
    {
      return sqrt(((double) px - point.px) * ((double) px - point.px) +
                  ((double) py - point.py) * ((double) py - point.py));
    }
  
  double getDistanceFromOrigin() const
    {
      return sqrt(((double) px * px) + ((double) py * py));
    }
  
  void rotateByRadians(double radians, int64 x, int64 y)
    {
      double r = sqrt(((double) x) * x + ((double) y) * y);
      double oldTheta = (x != 0) ? atan(((double) y) / x) :
        ((y > 0) ? M_PI_2 : -M_PI_2);
      
      oldTheta += (x < 0) ? M_PI : 0;
      px = (UnitType) (.5 + r * cos(oldTheta - radians));
      py = (UnitType) (.5 + r * sin(oldTheta - radians));
    }
  void rotateByRadians(double radians, int x, int y)
    {
      double r = sqrt(((double) x) * x + ((double) y) * y);
      double oldTheta = (x != 0) ? atan(((double) y) / x) :
        ((y > 0) ? M_PI_2 : -M_PI_2);
      
      oldTheta += (x < 0) ? M_PI : 0;
      px = (UnitType) (.5 + r * cos(oldTheta - radians));
      py = (UnitType) (.5 + r * sin(oldTheta - radians));
    }
  void rotateByRadians(double radians, double x, double y)
    {
      double r = sqrt(x * x + y * y);
      double oldTheta = (x != 0) ? atan(y / x) :
        ((y > 0) ? M_PI_2 : -M_PI_2);
      
      oldTheta += (x < 0) ? M_PI : 0;
      px = (UnitType) (r * cos(oldTheta - radians));
      py = (UnitType) (r * sin(oldTheta - radians));
    }
  void rotateByRadians(double radians)
    {
      rotateByRadians(radians, px, py);
    }
  void rotateByDegrees(double degrees)
    {
      rotateByRadians(degrees * M_PI / 180., px, py);
    }

  void add(const Point<UnitType> & point)
    {
      px += point.px;
      py += point.py;
    }
  void subtract(const Point<UnitType> & point)
    {
      px -= point.px;
      py -= point.py;
    }
  void scale(double scaler)
    {
      px *= scaler;
      py *= scaler;
    }
  void normalize()
    {
      double dfo = getDistanceFromOrigin();
      px /= dfo;
      py /= dfo;
    }

  bool isOnSegment(const Point<UnitType> & p1,
                   const Point<UnitType> & p2) const
    {
      return(isOnLine(p1, p2) &&
             ((px >= p1.px && px <= p2.px) || (px <= p1.px && px >= p2.px)) &&
             ((py >= p1.py && py <= p2.py) || (py <= p1.py && py >= p2.py)));
    }

  bool isOnLine(const Point<UnitType> & p1,
                const Point<UnitType> & p2) const
    {
      // see if slope between p1 & p2 equals slope between this & p1
      if(p2.px - p1.px == 0)
      {
        // line is vertical
        if(px - p1.px == 0)
        {
          // px is on line
          return true;
        }
        else
        {
          // px is not on line
          return false;
        }
      }
      
      // line is not vertical
      if(px - p2.px == 0)
      {
        // x is 'vertical' with end of segment
        if(py - p2.py == 0)
        {
          // x,y coincides with end of segment
          return true;
        }
        else
        {
          // x,y is off line
          return false;
        }
      }

      if((p2.py - p1.py) / (p2.px - p1.px) == (p2.py - py) / (p2.px - px))
        return true;
      return false;
    }
  // determines which side of p1p2 is this point
  ThreePointOrientation getLineSide(const Point<UnitType> & p1,
                                    UnitType p2x, UnitType p2y) const
    {
      double posOrNeg =
        ((double)p2x - p1.px) * ((double)p2y - py) -
        ((double)p2x - px) * ((double)p2y - p1.py);
      return((posOrNeg > 1e-14) ? TP_CLOCKWISE :
             ((posOrNeg < -1e-14) ? TP_COUNTERCW : TP_COLINEAR));
    }

  ThreePointOrientation getLineSide(const Point<UnitType> & p1,
                                    const Point<UnitType> & p2) const
    {
      return getLineSide(p1, p2.getX(), p2.getY());
    }

  bool operator!=(const Point<UnitType> & other) const
    {
      return(!(*this==other));
    }
  bool operator==(const Point<UnitType> & other) const
    {
      return(other.px == px && other.py == py);
    }
  bool operator<(const Point<UnitType> & other) const
    {
      if(px < other.px) return true;
      if(px > other.px) return false;
      if(py < other.py) return true;
      if(py > other.py) return false;
      return false;
    }
  friend ostream & operator<<(ostream & os, const Point<UnitType> & point)
    {
      os << point.px << "," << point.py;
      return os;
    }
  
protected:
  UnitType px;
  UnitType py;
};

#endif
