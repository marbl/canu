
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
/* $Id: IntervalClique.h,v 1.2 2004-09-23 20:25:23 mcschatz Exp $ */
#ifndef CLIQUE_H
#define CLIQUE_H

#include <iostream>
#include <list>

#include "Interval.h"

template <class IDType, class UnitType>
class IntervalClique
{
public:
  IntervalClique()
    {
      this->interval.setMin(0);
      this->interval.setMax(0);
    }

  IntervalClique(const list<Interval<IDType, UnitType> > & intervals,
                 UnitType min, UnitType max)
    {
      typename list<Interval<IDType, UnitType> >::const_iterator iter;
      
      for(iter = intervals.begin(); iter != intervals.end(); iter++)
        pids.push_back(iter->getID());
      
      pids.sort();
      pinterval.setMin(min);
      pinterval.setMax(max);
    }
  
  IntervalClique(const list<IDType> ids,
                 const list<Interval<IDType, UnitType> > & interval)
    {
      populateIDList(pids, ids);
      pinterval = interval;
    }

  IntervalClique(IDType id, const Interval<IDType, UnitType> & interval)
    {
      pids.push_back(id);
      pinterval = interval;
    }

  void addID(IDType id)
    {
      typename list<IDType>::iterator iter;
      for(iter = pids.begin(); *iter < id && iter != pids.end(); iter++);
      pids.insert(iter, id);
    }

  int getNumIDs() const
    {
      return pids.size();
    }
  const list<IDType> & getIDs() const
    {
      return pids;
    }
  void getCommonIDs(const list<IDType> & ids, list<IDType> & commonIDs) const
    {
      typename list<IDType>::const_iterator piter;
      typename list<IDType>::const_iterator iter;

      commonIDs.clear();
      for(iter = ids.begin(), piter = pids.begin();
          iter != ids.end() && piter != pids.end();)
      {
        if(*iter == *piter)
        {
          commonIDs.push_back(*iter);
          iter++;
          piter++;
        }
        else if(*iter < *piter)
          iter++;
        else
          piter++;
      }
    }
      
  void getCommonIDs(const IntervalClique<IDType, UnitType> & clique,
                    list<IDType> & commonIDs) const
    {
      getCommonIDs(clique.pids, commonIDs);
    }

  const Interval<IDType, UnitType> & getInterval() const
    {
      return pinterval;
    }
  UnitType getIntervalMin() const
    {
      return pinterval.getMin();
    }
  UnitType getIntervalMax() const
    {
      return pinterval.getMax();
    }

  int getNumMembers() const
    {
      return pids.size();
    }

  bool operator<(const IntervalClique<IDType, UnitType> & other) const
    {
      return(pinterval.getMin() < other.pinterval.getMin());
    }

  friend ostream & operator<<(ostream &os,
                              const IntervalClique<IDType, UnitType> & c)
    {
      os << c.pinterval.getMin() << " , " << c.pinterval.getMax()
         << " : " << c.getNumMembers() << " members: ";
      typename list<IDType>::const_iterator iter;
      for(iter = c.pids.begin(); iter != c.pids.end(); iter++)
        os << "  " << *iter;
      os << endl;
      return os;
    }
  
private:
  /*
  void populateIDListSorted(list<int> & to,
                            const list<Interval<IDType, UnitType> > & from) const
    {
      to.clear();
      
      typename list<Interval<IDType, UnitType> >::iterator iter = from.begin();
      for(last = iter->getID(); iter != from.end(); iter++)
        to.push_back(iter->getID());
    }
  */
  /*
  void populateIDList(list<int> & to,
                      const list<Interval<IDType, UnitType> > & from) const
    {
      int last;
      bool inOrder = true;

      to.clear();
      
      // iterate through from & add to to
      // find out if they're sorted low to high
      typename list<Interval<IDType, UnitType> >::iterator iter = from.begin();
      for(last = iter->getID(); iter != from.end(); iter++)
      {
        to.push_back(iter->getID());
        inOrder = (iter->getID() < last) ? false : inOrder;
        last = iter->getID();
      }
      
      // sort to, if necessary
      if(!inOrder)
        to.sort();
    }
  */
  void populateIDListSorted(list<IDType> & to, const list<IDType> & from) const
    {
      to = from;
    }

  void populateIDList(list<IDType> & to, const list<IDType> & from) const
    {
      IDType last;
      bool inOrder = true;
      typename list<IDType>::const_iterator iter;
      
      // iterate through from & add to to
      // find out if they're sorted low to high
      for(iter = from.begin(), last = *iter;
          iter != from.end();
          iter++)
      {
        to.push_back(*iter);
        inOrder = (*iter < last) ? false : inOrder;
        last = *iter;
      }
      
      // sort to, if necessary
      if(!inOrder)
        to.sort();
    }
  
  list<IDType> pids;
  Interval<IDType, UnitType> pinterval;
};

#endif
