
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
/* $Id: TwoDIntervalClique.h,v 1.1.1.1 2004-04-14 13:52:08 catmandew Exp $ */
#ifndef TWODINTERVALCLIQUE_H
#define TWODINTERVALCLIQUE_H

#include <vector>
#include <list>
#include <map>

#include "Quadrilateral.h"
#include "Rectangle.h"
#include "IntervalClique.h"
#include "CliquePairIntersection.h"

template <class IDType, class UnitType>
class TwoDIntervalClique
{
public:
  TwoDIntervalClique(const Quadrilateral<IDType, UnitType> & quad)
    {
      pids.push_back(quad.getID());
      pxinterval = quad.getXProjection();
      pyinterval = quad.getYProjection();
    }
  
  TwoDIntervalClique(const IntervalClique<IDType, UnitType> clique,
                     bool isXClique,
                     const vector<Quadrilateral<IDType, UnitType> > & quads,
                     map<IDType, int> & id2Quad)
    {
      pids = clique.getIDs();
      list<IDType>::iterator iter = pids.begin();

      if(isXClique)
      {
        pxinterval = clique.getInterval();
        
        pyinterval = quads[id2Quad[(*iter)]].getYProjection();
        iter++;
        for(; iter != pids.end(); iter++)
          pyinterval.setIntervalToUnion(quads[id2Quad[(*iter)]].getYProjection());
      }
      else
      {
        pyinterval = clique.getInterval();
      
        pxinterval = quads[id2Quad[(*iter)]].getXProjection();
        iter++;
        for(; iter != pids.end(); iter++)
          pxinterval.setIntervalToUnion(quads[id2Quad[(*iter)]].getXProjection());
      }
    }

  TwoDIntervalClique(const Rectangle<IDType, UnitType> & r)
    {
      pids.push_back(r.getID());
      pxinterval = r.getXInterval();
      pyinterval = r.getYInterval();
    }
  
  TwoDIntervalClique(const CliquePairIntersection<IDType, UnitType> & cpi,
                     const vector<Quadrilateral<IDType, UnitType> > & quads,
                     map<IDType, int> & id2Quad)
    {
      pids = cpi.getCommonIDs();
      // iterate over ids to find x & y intersection intervals
      list<IDType>::iterator iter = pids.begin();

      pxinterval = quads[id2Quad[IDType (*iter)]].getXProjection();
      pyinterval = quads[id2Quad[IDType (*iter)]].getYProjection();
      for(++iter; iter != pids.end(); iter++)
      {
        pxinterval.setIntervalToIntersection(quads[id2Quad[IDType (*iter)]].getXProjection());
        pyinterval.setIntervalToIntersection(quads[id2Quad[IDType (*iter)]].getYProjection());
      }
    }

  unsigned int getNumIDs() const {return pids.size();}
  const list<IDType> & getIDs() const {return pids;}
  
  void getCommonIDs(const list<IDType> & ids, list<IDType> & commonIDs) const
    {
      list<IDType>::const_iterator piter;
      list<IDType>::const_iterator iter;

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

  void addRectangle(const Rectangle<IDType, UnitType> & r)
    {
      pxinterval.setIntervalToUnion(r.getXInterval());
      pyinterval.setIntervalToUnion(r.getYInterval());
      pids.push_back(r.getID());
    }
  void add(const TwoDIntervalClique<IDType, UnitType> & t)
    {
      pxinterval.setIntervalToUnion(t.getXInterval());
      pyinterval.setIntervalToUnion(t.getYInterval());

      list<IDType>::const_iterator iter;
      for(iter = t.pids.begin(); iter != t.pids.end(); iter++)
        pids.push_back(*iter);
    }
  const Interval<IDType, UnitType> & getXInterval() const
    {
      return pxinterval;
    }
  UnitType getXIntervalMin() const
    {
      return pxinterval.getMin();
    }
  UnitType getXIntervalMax() const
    {
      return pxinterval.getMax();
    }
  const Interval<IDType, UnitType> & getYInterval() const
    {
      return pyinterval;
    }
  UnitType getYIntervalMin() const
    {
      return pyinterval.getMin();
    }
  UnitType getYIntervalMax() const
    {
      return pyinterval.getMax();
    }
  bool xIntervalIntersects(const Interval<IDType, UnitType> & xi) const
    {
      return pxinterval.intersects(xi);
    }
  bool yIntervalIntersects(const Interval<IDType, UnitType> & yi) const
    {
      return pyinterval.intersects(yi);
    }
  bool intersects(const Rectangle<IDType, UnitType> & r) const
    {
      return(pxinterval.intersects(r.getXInterval()) &&
             pyinterval.intersects(r.getYInterval()));
    }
  bool intersects(const TwoDIntervalClique<IDType, UnitType> & t) const
    {
      return(pxinterval.intersects(t.getXInterval()) &&
             pyinterval.intersects(t.getYInterval()));
    }
  
  void rotateByRadians(double radians)
    {
      Quadrilateral<IDType, UnitType> q;
      q.reset(pxinterval.getMin(), pyinterval.getMin());
      q.setPoint(pxinterval.getMin(), pyinterval.getMax(), 1);
      q.setPoint(pxinterval.getMax(), pyinterval.getMax(), 2);
      q.setPoint(pxinterval.getMax(), pyinterval.getMin(), 3);
      q.rotateByRadians(radians);
      pxinterval = q.getXProjection();
      pyinterval = q.getYProjection();
    }
  void rotateByDegrees(double degrees)
    {
      rotateByRadians(degrees * M_PI / 180.);
    }
  bool operator<(const TwoDIntervalClique<IDType, UnitType> & other) const
    {
      return(getXIntervalMax() < other.getXIntervalMax());
    }
  

  friend ostream & operator<<(ostream &os,
                              const TwoDIntervalClique<IDType, UnitType> & t)
    {
      os << t.pxinterval.getMin() << " " << t.pyinterval.getMin() << "\t"
         << t.pxinterval.getMin() << " " << t.pyinterval.getMax() << "\t"
         << t.pxinterval.getMax() << " " << t.pyinterval.getMax() << "\t"
         << t.pxinterval.getMax() << " " << t.pyinterval.getMin() << "\t";
      list<IDType>::const_iterator iter;
      for(iter = t.pids.begin(); iter != t.pids.end(); iter++)
        os << *iter << " ";
      return os;
    }
private:
  Interval<IDType, UnitType> pxinterval;
  Interval<IDType, UnitType> pyinterval;
  list<IDType> pids;
};

#endif
