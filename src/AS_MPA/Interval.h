
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
/* $Id: Interval.h,v 1.4 2005-03-22 19:48:56 jason_miller Exp $ */
#ifndef INTERVAL_H
#define INTERVAL_H

#include <iostream>

using namespace std;

//#define DO_CHECKS_INTERVAL

template <class IDType, class UnitType>
class Interval
{
public:

  Interval()
    {
      pid = 0;
      pmin = pmax = 0;
    }

  Interval(const Interval<IDType, UnitType> & interval)
    {
      *this = interval;
    }
  Interval(const Interval<IDType, UnitType> & interval, const IDType & id)
    {
      pmin = interval.getMin();
      pmax = interval.getMax();
      pid = id;
    }
  Interval(const UnitType & min, const UnitType & max, const IDType & id)
    {
      pmin = min;
      pmax = max;
      pid = id;
    }

  void setID(IDType newID) {pid = newID;}
  const IDType getID() const {return pid;}
  
  void setMin(const UnitType & min)
    {
      pmin = min;
    }
  void setMin(const UnitType & min1, const UnitType & min2)
    {
      pmin = (min1 < min2) ? min1 : min2;
    }
  void setMax(const UnitType & max)
    {
      pmax = max;
    }
  void setMax(const UnitType & max1, const UnitType & max2)
    {
      pmax = (max1 < max2) ? max2 : max1;
    }
  const UnitType getMin() const {return pmin;}
  const UnitType getMax() const {return pmax;}

  void setInterval(const Interval<IDType, UnitType> & iv)
    {
      pid = iv.pid;
      pmin = iv.pmin;
      pmax = iv.pmax;
    }
  void setInterval(const UnitType & min, const UnitType & max)
    {
      setMin(min);
      setMax(max);
    }
  
  void setIntervalToUnion(const Interval & other)
    {
      setMin((getMin() < other.getMin()) ? getMin() : other.getMin());
      setMax((getMax() > other.getMax()) ? getMax() : other.getMax());
    }
  void setIntervalToIntersection(const Interval & other)
    {
#ifdef DO_CHECKS_INTERVAL
      if(getMin() > other.getMax() || getMax() < other.getMin())
      {
        cerr << "No interval intersection found!: ("
             << getMin() << ", " << getMax() << ") vs. ("
             << other.getMin() << ", " << other.getMax() << ")\n";
        return;
      }
#endif
      setMin((getMin() > other.getMin()) ? getMin() : other.getMin());
      setMax((getMax() < other.getMax()) ? getMax() : other.getMax());
    }
  void getIntersection(const Interval & other, Interval & newI)
    {
      newI.setID(pid);
      newI.setMin((getMin() > other.getMin()) ? getMin() : other.getMin());
      newI.setMax((getMax() < other.getMax()) ? getMax() : other.getMax());
    }
  
  bool startsBefore(const UnitType & min) const
    {
      return (getMin() < min);
    }
  bool startsBefore(const Interval<IDType, UnitType> & interval) const
    {
      return startsBefore(interval.getMin());
    }
  
  bool operator<(const Interval<IDType, UnitType> & other) const
    {
      return startsBefore(other);
    }

  bool operator==(const Interval<IDType, UnitType> & other) const
    {
      return(min == getMin() && max == getMax());
    }
  
  bool intersects(const UnitType & min, const UnitType & max) const
    {
      return (min < getMax() && getMin() < max);
    }

  bool intersects(const Interval<IDType, UnitType> & interval) const
    {
      return intersects(interval.getMin(), interval.getMax());
    }

  friend ostream & operator<<(ostream &os,
                              const Interval<IDType, UnitType> & interval)
    {
      os << interval.getID() << ": " << interval.getMin() << ", "
         << interval.getMax();
      return os;
    }
  
private:
  IDType pid;
  UnitType pmin;
  UnitType pmax;
};

#endif
