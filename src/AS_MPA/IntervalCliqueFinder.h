
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
/* $Id: IntervalCliqueFinder.h,v 1.3 2005-03-22 19:05:15 jason_miller Exp $ */
#ifndef INTERVALCLIQUEFINDER_H
#define INTERVALCLIQUEFINDER_H

#include <list>

#include "IntervalClique.h"

//#define DEBUG_ICF
//#define DO_CHECKS_ICF

template <class IDType, class UnitType>
class IntervalCliqueFinder
{
public:
/*
  Find cliques in one dimension given a set of intervals, I

  IntervalClique is: start, end, list of interval IDs
  
  0. Create empty set of cliques
  1. Sort I low to high on start value
  2. Create dynamic list, L, to insertion sort items into & pop out of
     based on end value
  3. track last action (add, delete interval) and coordinate
  4. iterate over elements in I
     {
       if end coordinate of top element in L is less than start
       start coordinate in next element in I
       {
         if last action was insertion into L
         {
           create new clique
           members are set in L
           start is start of last inserted element
           end is end of top element
         }
         pop top element off L
         set last action to deleted-an-interval
       }
       else
       {
         insertion sort interval into L st. next interval in L
           has higher end value or there is no next interval
         set last action to added-an-interval
         set last a-a-i coordinate to interval start
       }
     }
 */
  void findCliques(const list<Interval<IDType, UnitType > > & intervals,
                   list<IntervalClique<IDType, UnitType> > & cliques)
    {
      list<Interval<IDType, UnitType> > myIntervals(intervals);
      typename list<Interval<IDType, UnitType> >::iterator miter;
      list<Interval<IDType, UnitType> > ends;
      typename list<Interval<IDType, UnitType> >::iterator eiter;
      UnitType lastStart;
      bool lastAction; // true = add, false = delete
      
      cliques.clear();
      if(intervals.size() == 0) return;
      
      myIntervals.sort();
      
      miter = myIntervals.begin();
      lastStart = miter->getMin();
      lastAction = true;
      ends.push_back(*miter);
      
      for(miter++; miter != myIntervals.end();)
      {
        
#ifdef DEBUG_ICF
        cerr << "\nAt top of loop\n";
        printCurrentState(*miter, ends);
#endif
        
        if(ends.size() == 0 || miter->getMin() < ends.front().getMax())
        {
          // next starting interval starts before next ending interval ends
          lastStart = miter->getMin();
          lastAction = true;
          
          eiter = ends.begin();
          while(miter->getMax() > eiter->getMax() &&
                eiter != ends.end()) eiter++;

#ifdef DO_CHECKS_ICF
          if(ends.size() != 0 &&
             eiter != ends.end() &&
             miter->getMax() > eiter->getMax())
          {
            cerr << "Error insertion sorting intervals!\n";
          }
#endif
          
          ends.insert(eiter, *miter);
          miter++;

#ifdef DEBUG_ICF
          cerr << "\nAfter inserting new interval into current set\n";
          printCurrentState(*miter, ends);
#endif

        }
        else
        {
          // next ending interval ends at or before next starting interval starts
          if(lastAction == true)
          {
            // new clique: add followed by delete
            IntervalClique<IDType, UnitType> c(ends,
                                               lastStart,
                                               ends.front().getMax());
            cliques.push_back(c);

#ifdef DEBUG_ICF
            cerr << "\nNEW CLIQUE:\n";
            cerr << c << endl;
#endif

          }
          ends.pop_front();
          lastAction = false;

#ifdef DEBUG_ICF
          cerr << "\nAfter popping off of current intervals\n";
          printCurrentState(*miter, ends);
#endif

        }
      }
      
      if(lastAction == true)
      {
        // new clique: add followed by delete
        IntervalClique<IDType, UnitType> c(ends,
                                           lastStart,
                                           ends.front().getMax());
        cliques.push_back(c);
        
#ifdef DEBUG_ICF
        cerr << "\nNEW CLIQUE:\n";
        cerr << c << endl;
#endif
        
      }
    }

  void printCliques(char * comment,
                    const list<IntervalClique<IDType, UnitType> > & cliques,
                    ostream & fout)
    {
      typename list<IntervalClique<IDType, UnitType> >::const_iterator iter;

      fout << comment << endl;

      int i;
      for(i = 0, iter = cliques.begin(); iter != cliques.end(); iter++, i++)
        fout << i << ": " << *iter;
    }

private:

  void printCurrentState(const Interval<IDType, UnitType> & ni,
                         const list<Interval<IDType, UnitType> > & cis)
    {
      cerr << "Next interval: " << ni << endl;
      cerr << cis.size() << " current intervals: " << endl;
      if(cis.size() == 0)
        cerr << "  none.\n";
      else
      {
        typename list<Interval<IDType, UnitType> >::const_iterator iter;
        for(iter = cis.begin(); iter != cis.end(); iter++)
          cerr << "  " << *iter << endl;
      }
    }
};
  
#endif
