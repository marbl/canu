
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
/* $Id: Quadrilateral.h,v 1.5 2005-07-13 14:47:55 brianwalenz Exp $ */
#ifndef QUADRILATERAL_H
#define QUADRILATERAL_H

#include <cmath>
#include <algorithm>

#include "Point.h"
#include "Interval.h"

template <class IDType, class UnitType>
class Quadrilateral
{
public:

  /*
   */
  Quadrilateral(IDType id = 0) :
    pid(id)
    {
      Point<UnitType> p(0,0);
      reset(p);
      xProjection.setID(pid);
      yProjection.setID(pid);
    }
  Quadrilateral(const Quadrilateral<IDType, UnitType> & q)
    {
      pid = q.pid;
      for(int i = 0; i < 4; i++)
        points[i] = q.points[i];
      xProjection = q.xProjection;
      yProjection = q.yProjection;
      xProjection.setID(pid);
      yProjection.setID(pid);
    }

  Quadrilateral(Point<UnitType> const ps[4], IDType id = 0) :
    pid(id)
    {
      points[0] = this->p0;
      initializeXProjection(ps[0]);
      initializeYProjection(ps[0]);

      for(int i = 1; i < 4; i++)
        setPoint(ps[i], i);
      xProjection.setID(pid);
      yProjection.setID(pid);
    }

  Quadrilateral(const UnitType x0, const UnitType y0,
                const UnitType x1, const UnitType y1,
                const UnitType x2, const UnitType y2,
                const UnitType x3, const UnitType y3,
                IDType id = 0) :
    pid(id)
    {
      points[0].setX(x0);
      points[0].setY(y0);
      initializeXProjection(points[0]);
      initializeYProjection(points[0]);

      setPoint(x1, y1, 1);
      setPoint(x2, y2, 2);
      setPoint(x3, y3, 3);
    }
  Quadrilateral(const Point<UnitType> & p0, const Point<UnitType> & p1,
                const Point<UnitType> & p2, const Point<UnitType> & p3,
                IDType id = 0) :
    pid(id)
    {
      points[0] = p0;
      initializeXProjection(p0);
      initializeYProjection(p0);

      setPoint(p1, 1);
      setPoint(p2, 2);
      setPoint(p3, 3);
    }

  void reset(const Point<UnitType> & p)
    {
      points[0] = p;
      initializeXProjection(p);
      initializeYProjection(p);
    }

  void reset(const UnitType & x, const UnitType & y)
    {
      points[0].setX(x);
      points[0].setY(y);
      initializeXProjection(points[0]);
      initializeYProjection(points[0]);
    }

  void setPoint(const UnitType & x, const UnitType & y, int pointNumber)
    {
      points[pointNumber].setX(x);
      points[pointNumber].setY(y);
      updateXProjection(points[pointNumber]);
      updateYProjection(points[pointNumber]);
    }
  void setPoint(const Point<UnitType> & p, int pointNumber)
    {
      points[pointNumber] = p;
      updateXProjection(p);
      updateYProjection(p);
    }
  
  void setID(IDType newID)
    {
      pid = newID;
      xProjection.setID(pid);
      yProjection.setID(pid);
    }
  IDType getID() const {return pid;}
  
  const Interval<IDType, UnitType> & getXProjection() const
    {
      return xProjection;
    }

  const Interval<IDType, UnitType> & getYProjection() const
    {
      return yProjection;
    }

  bool xProjectionsIntersect(const Quadrilateral<IDType, UnitType> & q) const
    {
      return(xProjection.intersects(q.getXProjection()));
    }
  bool yProjectionsIntersect(const Quadrilateral<IDType, UnitType> & q) const
    {
      return(yProjection.intersects(q.getYProjection()));
    }
  bool projectionsIntersect(const Quadrilateral<IDType, UnitType> & q) const
    {
      return(xProjectionsIntersect(q) & yProjectionsIntersect(q));
    }

  void rotateByRadians(double radians)
    {
      points[0].rotateByRadians(radians);
      initializeXProjection(points[0]);
      initializeYProjection(points[0]);

      for(int i = 1; i < 4; i++)
      {
        points[i].rotateByRadians(radians);
        updateXProjection(points[i]);
        updateYProjection(points[i]);
      }
    }
  
  void rotateByDegrees(double degrees)
    {
      rotateByRadians(degrees * M_PI / 180.);
    }
  
private:
  void initializeXProjection(const Point<UnitType> & p)
    {
      xProjection.setMin(p.getX());
      xProjection.setMax(p.getX());
    }
  void initializeYProjection(const Point<UnitType> & p)
    {
      yProjection.setMin(p.getY());
      yProjection.setMax(p.getY());
    }
  
  void updateXProjection(const Point<UnitType> & p)
    {
      xProjection.setMin(xProjection.getMin(), p.getX());
      xProjection.setMax(xProjection.getMax(), p.getX());
    }

  void updateYProjection(const Point<UnitType> & p)
    {
      yProjection.setMin(yProjection.getMin(), p.getY());
      yProjection.setMax(yProjection.getMax(), p.getY());
    }

protected:
  IDType pid;
  Interval<IDType, UnitType> xProjection;
  Interval<IDType, UnitType> yProjection;
  Point<UnitType> points[4];
};


#endif
