
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
/* $Id: Trapezoid.h,v 1.5 2005-07-13 14:47:55 brianwalenz Exp $ */
#ifndef TILTEDTRAPEZOID_H
#define TILTEDTRAPEZOID_H

#include <cmath>
#include <algorithm>

#include "Point.h"
#include "Interval.h"
#include "Quadrilateral.h"

template <class IDType, class UnitType>
class TiltedTrapezoid : public Quadrilateral
{
public:
  TiltedTrapezoid(IDType id = 0) :
    Quadrilateral(id)
    {
    }
  TiltedTrapezoid(const TiltedTrapezoid<IDType, UnitType> & tt)
    Quadrilateral(tt)
    {
      pUpper = tt.pUpper;
      pPosSlope = tt.pPosSlope;
    }

  TiltedTrapezoid(Point<UnitType> const ps[4], IDType id = 0) :
    Quadrilateral(ps, id)
    {
      CheckAndSetTTDetails();
    }
  
  TiltedTrapezoid(const UnitType x0, const UnitType y0,
                  const UnitType x1, const UnitType y1,
                  const UnitType x2, const UnitType y2,
                  const UnitType x3, const UnitType y3,
                  IDType id = 0) :
    Quadrilateral(x0, y0, x1, y1, x2, y2, x3, y3, id)
    {
      CheckAndSetTTDetails();
    }
  
  TiltedTrapezoid(const Point<UnitType> & p0, const Point<UnitType> & p1,
                  const Point<UnitType> & p2, const Point<UnitType> & p3,
                  IDType id = 0) :
    Quadrilateral(p0, p1, p2, p3, id)
    {
      CheckAndSetTTDetails();
    }
  
  void setPoint(const UnitType & x, const UnitType & y, int pointNumber)
    {
      Quadrilateral::setPoint(x, y, pointNumber);
      CheckAndSetTTDetails();
    }
  void setPoint(const Point<UnitType> & p, int pointNumber)
    {
      Quadrilateral::setPoint(p, pointNumber);
      CheckAndSetTTDetails();
    }

  const Interval<IDType, UnitType> & getXIntercepts() const
    {
      return xIntercepts;
    }
  
  bool xInterceptsIntersect(const TiltedTrapezoid<IDType, UnitType> && tt)
    {
      xIntercepts.intersects(tt.getXIntercepts());
    }
  bool projectionsIntersect(const TiltedTrapezoid<IDType, UnitType> && tt)
    {
      return(xProjectionsIntersect(tt) && yProjectionsIntersect(tt) &&
             xInterceptsIntersect(tt));
    }
  bool intersects(const TiltedTrapezoid<IDType, UnitType> & tt)
    {
      assert(tt.pUpper = pUpper);
      assert(tt.pPosSlope = pPosSlope);

      return(projectionsIntersect(tt));
    }

  bool getIntersection(const TiltedTrapezoid<IDType, UnitType> & tt,
                       TiltedTrapezoid<IDType, UnitType> & rettt)
    {
      if(!intersects(tt)) return false;

      xProjection.getIntersection(tt.getXProjection(), rettt.xProjection);
      yProjection.getIntersection(tt.getYProjection(), rettt.yProjection);
      xIntercepts.getIntersection(tt.getXIntercepts(), rettt.xIntercepts);
      rettt.pUpper = pUpper;
      rettt.pPosSlope = pPosSlope;

      // trick is to set the points
      if(pUpper)
      {
        if(pPosSlope)
        {
          rettt.points[0].setX(rettt.xProjection.getMin());
          rettt.points[0].setY(rettt.yProjection.getMin());
          
          rettt.points[1].setX(rettt.xProjection.getMin());
          rettt.points[1].setY(rettt.xProjection.getMin() +
                               rettt.xIntercepts.getMin());
          
          rettt.points[2].setX(rettt.xIntercepts.getMin() -
                               rettt.yProjection.getMax());
          rettt.points[2].setY(rettt.yProjection.getMax());
          
          rettt.points[3].setX(rettt.xProjection.getMax());
          rettt.points[3].setY(rettt.yProjection.getMax());
        }
        else
        {
          rettt.points[0].setX(rettt.xProjection.getMin());
          rettt.points[0].setY(rettt.yProjection.getMax());

          rettt.points[1].setX(rettt.yProjection.getMax() +
                               rettt.xIntercepts.getMax());
          rettt.points[1].setY(rettt.yProjection.getMax());
          
          rettt.points[2].setX(rettt.xProjection.getMax());
          rettt.points[2].setY(rettt.xProjection.getMax() -
                               rettt.xIntercepts.getMax());

          rettt.points[3].setX(rettt.xProjection.getMax());
          rettt.points[3].setY(rettt.yProjection.getMin());
        }
      }
      else
      {
        if(pPosSlope)
        {
          rettt.points[0].setX(rettt.xProjection.getMin());
          rettt.points[0].setY(rettt.yProjection.getMin());

          rettt.points[1].setX(rettt.xProjection.getMax());
          rettt.points[1].setY(rettt.yProjection.getMax());
          
          rettt.points[2].setX(rettt.xProjection.getMax());
          rettt.points[2].setY(rettt.xProjection.getMax() +
                               rettt.xIntercepts.getMax());

          rettt.points[3].setX(rettt.yProjection.getMin() -
                               rettt.xIntercepts.getMax());
          rettt.points[3].setY(rettt.yProjection.getMin());
        }
        else
        {
          rettt.points[0].setX(rettt.xProjection.getMin());
          rettt.points[0].setY(rettt.xProjection.getMin() -
                               rettt.xIntercepts.getMin());

          rettt.points[1].setX(rettt.xProjection.getMin());
          rettt.points[1].setY(rettt.yProjection.getMax());
          
          rettt.points[2].setX(rettt.xProjection.getMax());
          rettt.points[2].setY(rettt.yProjection.getMin());

          rettt.points[3].setX(rettt.yProjection.getMin() +
                               rettt.xIntercepts.getMin());
          rettt.points[3].setY(rettt.yProjection.getMin());
        }
      }
    }
private:
  void CheckAndSetTTDetails()
    {
      // assume points are ordered clockwise starting with
      // assume 45 degree angle
      // leftmost, lower point
      if(points[0].getY() == points[1].getY())
      {
        /*
          point 1 is to the right of point 0
          upper, negative slope

          0      1


                      2


                      3      
        */
        assert(points[0].getX() < points[1].getX());
        pUpper = true;
        pPosSlope = false;
        xIntercepts.setMin(points[0].getY() + points[0].getX());
        xIntercepts.setMax(points[1].getY() + points[1].getX());
      }
      else if(points[0].getX() != points[1].getX())
      {
        /*
          point 1 is up and to the right of point 0
          lower, positive slope


                     1


                     2

          
          0     3
          
        */
        pUpper = false;
        pPosSlope = true;
        xIntercepts.setMin(points[0].getY() - points[0].getX());
        xIntercepts.setMax(points[3].getY() - points[3].getX());
      }
      else if(points[3].getY() < points[0].getY())
      {
        /*
          point 3 is down and to the right of point 0
          lower, negative slope


          1


          0

          
                 3     2
          
        */
        assert(points[3].getX() >= points[0].getX());
        pUpper = false;
        pPosSlope = false;
        xIntercepts.setMin(points[0].getY() + points[0].getX());
        xIntercepts.setMax(points[1].getY() + points[1].getX());
      }
      else
      {
        /*
          point 3 is up and to the right of point 0
          upper, positive slope


                 2     3

                 
          1


          0
        */
        pUpper = true;
        pPosSlope = true;
        xIntercepts.setMin(points[1].getY() - points[1].getX());
        xIntercepts.setMax(points[0].getY() - points[0].getX());
      }
    }
  
  bool pUpper;
  bool pPosSlope;
  Interval<IDType, UnitType> xIntercepts;
};


#endif
