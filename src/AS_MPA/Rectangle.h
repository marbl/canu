
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
/* $Id: Rectangle.h,v 1.5 2005-10-14 20:34:09 catmandew Exp $ */
#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "Interval.h"

template <class IDType, class UnitType>
class Rectangle
{
public:

  Rectangle(const Rectangle & rect)
    {
      pid = rect.pid;
      pintervals[0] = rect.pintervals[0];
      pintervals[1] = rect.pintervals[1];
      pintervals[0].setID(pid);
      pintervals[1].setID(pid);
    }
  
  Rectangle(const UnitType & x1, const UnitType & x2,
            const UnitType & y1, const UnitType & y2,
            IDType id) :
    pid(id)
    {
      pintervals[0].setInterval(x1, x2);
      pintervals[1].setInterval(y1, y2);
      pintervals[0].setID(pid);
      pintervals[1].setID(pid);
    }
  
  Rectangle(const Interval<IDType, UnitType> & xInterval,
            const Interval<IDType, UnitType> & yInterval)
    {
      pid = xInterval.getID();
      pintervals[0].setInterval(xInterval);
      pintervals[1].setInterval(yInterval);
      pintervals[0].setID(pid);
      pintervals[1].setID(pid);
    }

  Rectangle(const Interval<IDType, UnitType> & xInterval,
            const Interval<IDType, UnitType> & yInterval,
            IDType id) :
    pid(id)
    {
      pintervals[0].setInterval(xInterval);
      pintervals[1].setInterval(yInterval);
      pintervals[0].setID(pid);
      pintervals[1].setID(pid);
    }

  void setID(IDType id)
    {
      pid = id;
    }

  IDType getID() const {return pid;}

  Interval<IDType, UnitType> getInterval(int dimension) const
    {
      return pintervals[dimension];
    }
  Interval<IDType, UnitType> getXInterval() const
    {
      return getInterval(0);
    }
  Interval<IDType, UnitType> getYInterval() const
    {
      return getInterval(1);
    }

  void setInterval(const UnitType & a1, const UnitType & a2, int dimension)
    {
      pintervals[dimension].setInterval(a1, a2);
    }
  void setXInterval(const UnitType & a1, const UnitType & a2)
    {
      setInterval(a1, a2, 0);
    }
  void setYInterval(const UnitType & a1, const UnitType & a2)
    {
      setInterval(a1, a2, 1);
    }

  UnitType getXMin() const
    {
      return pintervals[0].getMin();
    }
  UnitType getXMax() const
    {
      return pintervals[0].getMax();
    }
  void setXMin(const UnitType & xmin)
    {
      pintervals[0].setMin(xmin);
    }
  void setXMax(const UnitType & xmax)
    {
      pintervals[0].setMax(xmax);
    }
  
  UnitType getYMin() const
    {
      return pintervals[1].getMin();
    }
  UnitType getYMax() const
    {
      return pintervals[1].getMax();
    }
  void setYMin(const UnitType & ymin)
    {
      pintervals[0].setMin(ymin);
    }
  void setYMax(const UnitType & ymax)
    {
      pintervals[0].setMax(ymax);
    }

  void setInterval(const Interval<IDType, UnitType> & interval, int dimension)
    {
      pintervals[dimension].setInterval(interval);
    }
  void setXInterval(const Interval<IDType, UnitType> & interval)
    {
      setInterval(interval, 0);
    }
  void setYInterval(const Interval<IDType, UnitType> & interval)
    {
      setInterval(interval, 1);
    }

  bool xIntervalIntersects(const UnitType & min, const UnitType & max) const
    {
      return pintervals[0].intersects(min, max);
    }
  bool xIntervalIntersects(const Interval<IDType, UnitType> & x) const
    {
      return pintervals[0].intersects(x);
    }
  
  bool intersects(const Rectangle<IDType, UnitType> & r) const
    {
      return (pintervals[0].intersects(r.getInterval[0]) &&
              pintervals[1].intersects(r.getInterval[1]));
    }

  bool operator<(const Rectangle<IDType, UnitType> & other) const
    {
      return(pintervals[0].getMin() < other.pintervals[0].getMin());
    }
  friend ostream & operator<<(ostream &os,
                              const Rectangle<IDType, UnitType> & r)
    {
      os << r.pintervals[0].getMin() << " " << r.pintervals[1].getMin() << "\t"
         << r.pintervals[0].getMin() << " " << r.pintervals[1].getMax() << "\t"
         << r.pintervals[0].getMax() << " " << r.pintervals[1].getMax() << "\t"
         << r.pintervals[0].getMax() << " " << r.pintervals[1].getMin() << "\t"
         << r.pid;
      return os;
    }
  
private:
  IDType pid;
  Interval<IDType, UnitType> pintervals[2];
};

#endif
