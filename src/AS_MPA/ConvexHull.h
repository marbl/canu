
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
/* $Id: ConvexHull.h,v 1.4 2005-03-22 19:48:56 jason_miller Exp $ */
#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include <iostream>
#include <list>
#include <vector>

#include "Point.h"
#include "Polygon.h"

//#define DEBUG_CONVEXHULL

template <class UnitType>
class CVHPoint : public Point<UnitType>
{
public:
  CVHPoint() {this->px = 0; this->py = 0; this->ppoly = this->pindex = -1;}
  CVHPoint(const UnitType & x, const UnitType & y,
           int poly = -1, int index = -1) : Point<UnitType> (x, y)
    {
      ppoly = poly;
      pindex = index;
    }
  CVHPoint(const Point<UnitType> & point, int poly = -1, int index = -1) : Point<UnitType> (point)
    {
      ppoly = poly;
      pindex = index;
    }
  CVHPoint(const CVHPoint<UnitType> & cvh) : Point<UnitType>(cvh)
    {
      ppoly = cvh.ppoly;
      pindex = cvh.pindex;
    }
  void set(const UnitType & x, const UnitType & y, int id = -1, int index = -1)
    {
      this->px = x;
      this->py = y;
      this->ppoly = id;
      this->pindex = index;
    }
  
  void setPolygon(int id) {ppoly = id;}
  void setIndex(int index) {pindex = index;}
  unsigned int getPolygon() const {return ppoly;}
  unsigned int getIndex() const {return pindex;}
  
  friend ostream & operator<<(ostream & os, const CVHPoint<UnitType> & point)
    {
      os << point.getX() << "," << point.getY()
         << "," << point.ppoly << "," << point.pindex;
      return os;
    }
  
private:
  unsigned int ppoly;
  unsigned int pindex;
};


template <class UnitType>
class ConvexHull
{
public:
  ConvexHull(){}
  ConvexHull(const list<CVHPoint<UnitType> > & points)
    {
      create(points);
    }
  ConvexHull(const list<Point<UnitType> > & points)
    {
      list<CVHPoint<UnitType> > mp = points;
      create(mp);
    }
  ConvexHull(const Polygon<UnitType> & p1,
             const Polygon<UnitType> & p2)
    {
      create(p1, p2);
    }

  void create(const list<CVHPoint<UnitType> > & points)
    {
      bool debugCreate = false;
      
      if(points.size() < 4)
      {
        typename list<CVHPoint<UnitType> >::const_iterator clfiter;
        for(clfiter = points.begin(); clfiter != points.end(); clfiter++)
          ppts.push_back(*clfiter);
        return;
      }
        
      list<CVHPoint<UnitType> > mp = points;
      mp.sort();

      // eliminate duplicates & keep lower polygon index
      typename list<CVHPoint<UnitType> >::iterator lfiter;
      lfiter = mp.begin();
      CVHPoint<UnitType> lastPoint = *lfiter;
      for(lfiter++; lfiter != mp.end(); lfiter++)
      {
        if(lastPoint.getX() == lfiter->getX() &&
           lastPoint.getY() == lfiter->getY())
        {
          if(lastPoint.getPolygon() == lfiter->getPolygon())
            lfiter--;
          mp.erase(lfiter);
          lfiter--;
        }
        lastPoint = *lfiter;
      }
      
      lfiter = mp.begin();
      for(int i = 0; i < 2; i++, lfiter++)
        ppts.push_back(*lfiter);
      for(; lfiter != mp.end(); lfiter++)
      {
        ppts.push_back(*lfiter);
        while(ppts.size() > 2 &&
              ppts[ppts.size()-1].getLineSide(ppts[ppts.size()-3],
                                              ppts[ppts.size()-2]) != TP_CLOCKWISE)
        {
          typename vector<CVHPoint<UnitType> >::iterator vfiter;
          vfiter = ppts.end();
          vfiter--;
          vfiter--;
          ppts.erase(vfiter);
        }
      }

//#ifdef DEBUG_CONVEXHULL
      if(debugCreate)
      {
        typename vector<CVHPoint<UnitType> >::iterator debugiter;
        cerr << "Upper convex hull:\n";
        for(debugiter = ppts.begin(); debugiter != ppts.end(); debugiter++)
          cerr << *debugiter << endl;
      }
//#endif
      
      vector<CVHPoint<UnitType> > lower;
      typename list<CVHPoint<UnitType> >::reverse_iterator lriter;
      lriter = mp.rbegin();
      for(int i = 0; i < 2; i++, lriter++)
        lower.push_back(*lriter);
      for(; lriter != mp.rend(); lriter++)
      {
        lower.push_back(*lriter);
        while(lower.size() > 2 &&
              lower[lower.size()-1].getLineSide(lower[lower.size()-3],
                                                lower[lower.size()-2]) != TP_CLOCKWISE)
        {
          typename vector<CVHPoint<UnitType> >::iterator vfiter;
          vfiter = lower.end();
          vfiter--;
          vfiter--;
          lower.erase(vfiter);
        }
      }

//#ifdef DEBUG_CONVEXHULL
      if(debugCreate)
      {
        typename vector<CVHPoint<UnitType> >::iterator debugiter;
        cerr << "Lower convex hull:\n";
        for(debugiter = lower.begin(); debugiter != lower.end(); debugiter++)
          cerr << *debugiter << endl;
      }
//#endif
      
      // now add lower points excluding first & last
      for(unsigned int i = 1; i < lower.size() - 1; i++)
        ppts.push_back(lower[i]);
      
//#ifdef DEBUG_CONVEXHULL
      if(debugCreate)
      {
        typename vector<CVHPoint<UnitType> >::iterator debugiter;
        cerr << "Convex hull:\n";
        for(debugiter = ppts.begin(); debugiter != ppts.end(); debugiter++)
          cerr << *debugiter << endl;
      }
//#endif
    }

  void create(const Polygon<UnitType> & p1,
              const Polygon<UnitType> & p2)
    {
      list<CVHPoint<UnitType> > points;
      for(unsigned int i = 0; i < p1.size(); i++)
      {
        CVHPoint<UnitType> p(p1[i], 0, i);
        points.push_back(p);
      }
      for(unsigned int i = 0; i < p2.size(); i++)
      {
        CVHPoint<UnitType> p(p2[i], 1, i);
        points.push_back(p);
      }
      create(points);
    }

  void reset() {ppts.clear();}
  unsigned int size() const {return ppts.size();}
  const CVHPoint<UnitType> & operator[](int index) const
    {
      return ppts[index];
    }
  CVHPoint<UnitType> & operator[](int index){return ppts[index];}

  friend ostream & operator<<(ostream & os, const ConvexHull<UnitType> & cvh)
    {
      os << cvh.size();
      for(unsigned int i = 0; i < cvh.ppts.size(); i++)
        os << " " << cvh.ppts[i];
      return os;
    }
    
private:
  vector<CVHPoint<UnitType> > ppts;
};


#endif // #define CONVEXHULL_H
