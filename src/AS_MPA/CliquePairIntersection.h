
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
/* $Id: CliquePairIntersection.h,v 1.3 2005-03-22 19:05:10 jason_miller Exp $ */
#ifndef CLIQUEPAIRINTERSECTION_H
#define CLIQUEPAIRINTERSECTION_H

#include <list>

#include "Interval.h"

//#define DEBUG_CPI

template <class IDType, class UnitType>
class CliquePairIntersection
{
public:
  CliquePairIntersection()
    {
      pxIndex = pyIndex = -1;
    }
  
  CliquePairIntersection(int x, int y, list<IDType> cIDs)
    {
      pxIndex = x;
      pyIndex = y;
      pcommonIDs = cIDs;
    }

  void setXIndex(int x) {pxIndex = x;}
  void setYIndex(int y) {pyIndex = y;}
  void setXInterval(Interval<IDType, UnitType> xInterval)
    {
      pxInterval = xInterval;
    }
  void setCommonIDs(const list<IDType> & cIDs)
    {
      pcommonIDs = cIDs;
    }

  void addID(const IDType & id)
    {
      typename list<IDType>::iterator iter = pcommonIDs.begin();
      while(*iter > id && iter != pcommonIDs.end()) iter++;
      this->ends.insert(iter, id);
    }

  int getNumCommonIDs() const
    {
      return pcommonIDs.size();
    }
  const list<IDType> & getCommonIDs() const
    {
      return pcommonIDs;
    }
  const Interval<IDType, UnitType> & getXInterval() const
    {
      return pxInterval;
    }
  UnitType getXIntervalMin() const
    {
      return pxInterval.getMin();
    }
  UnitType getXIntervalMax() const
    {
      return pxInterval.getMax();
    }
  int getXIndex() const
    {
      return pxIndex;
    }
  int getYIndex() const
    {
      return pyIndex;
    }
  void deleteCommonIDs(const list<IDType> & idList)
    {
      typename list<IDType>::const_iterator iter;
      typename list<IDType>::iterator piter;
      int i = 0;
      int pi = 0;
      int numDeleted = 0;
      // commonIDs & idList should be sorted low to high
      for(iter = idList.begin(), piter = pcommonIDs.begin();
          iter != idList.end() && piter != pcommonIDs.end();)
      {
        if(*iter == *piter)
        {
          numDeleted++;
          piter = pcommonIDs.erase(piter);
          iter++;
          i++;
        }
        else if(*iter < *piter)
        {
          iter++;
          i++;
        }
        else
        {
          piter++;
          pi++;
        }
      }
#ifdef DEBUG_CPI
      cerr << numDeleted << " deleted.\n";
#endif
    }

  bool operator==(const CliquePairIntersection<IDType, UnitType> & other) const
    {
      if(pcommonIDs.size() != other.pcommonIDs.size())
        return false;
      
      typename list<IDType>::const_iterator oiter;
      typename list<IDType>::const_iterator iter;
      for(oiter = other.pcommonIDs.begin(), iter = pcommonIDs.begin();
          oiter != other.pcommonIDs.end() && iter != pcommonIDs.end();
          oiter++, iter++)
      {
        if(*iter != *oiter)
          return false;
      }
      return true;
    }
  
  bool operator<(const CliquePairIntersection<IDType, UnitType> & other) const
    {
      typename list<IDType>::const_iterator oiter;
      typename list<IDType>::const_iterator iter;
      for(oiter = other.pcommonIDs.begin(), iter = pcommonIDs.begin();
          oiter != other.pcommonIDs.end() && iter != pcommonIDs.end();
          oiter++, iter++)
      {
        if(*iter < *oiter)
          return true;
        else if(*iter > *oiter)
          return false;
      }
      if(pcommonIDs.size() < other.pcommonIDs.size())
        return true;
      else
        return false;
    }

  friend ostream & operator<<(ostream &os,
                              const CliquePairIntersection<IDType, UnitType> & cpi)
    {
      os << cpi.pxIndex << " , " << cpi.pyIndex;

      typename list<IDType>::const_iterator iter;
      for(iter = cpi.pcommonIDs.begin(); iter != cpi.pcommonIDs.end(); iter++)
        os << " " << *iter;
      os << endl;
      return os;
    }
  
private:
  Interval<IDType, UnitType> pxInterval;
  int pxIndex;
  int pyIndex;
  list<IDType> pcommonIDs;
};

#endif
