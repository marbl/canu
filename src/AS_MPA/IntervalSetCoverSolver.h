
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
/* $Id: IntervalSetCoverSolver.h,v 1.1.1.1 2004-04-14 13:52:02 catmandew Exp $ */
#ifndef INTERVALSETCOVERSOLVER_H
#define INTERVALSETCOVERSOLVER_H

#include <vector>
#include <list>

#include "IntervalCliqueFinder.h"
#include "TwoDIntervalClique.h"

//#define DEBUG_ISCS
//#define FEEDBACK_ISCS
//#define DO_CHECKS_ISCS
#define DELETE_DUPLICATES

template <class IDType, class UnitType>
class IntervalSetCoverSolver
{
public:
  IntervalSetCoverSolver()
    {
    }

  IntervalSetCoverSolver(vector<Rectangle<IDType, UnitType> > & rects)
    {
      PopulateQuads(rects);
      
      vector<Rectangle<IDType, UnitType> >::iterator iter;
      list<Interval<IDType, UnitType > > xIntervals;
      list<Interval<IDType, UnitType > > yIntervals;
      for(iter = rects.begin(); iter != rects.end(); iter++)
      {
        xIntervals.push_back(iter->getXInterval());
        yIntervals.push_back(iter->getYInterval());
      }
#ifdef DEBUG_ISCS
      list<Interval<IDType, UnitType > >::iterator iiter;
      cerr << "X intervals:\n";
      for(iiter = xIntervals.begin(); iiter != xIntervals.end(); iiter++)
        cerr << *iiter << endl;
      cerr << "Y intervals:\n";
      for(iiter = yIntervals.begin(); iiter != yIntervals.end(); iiter++)
        cerr << *iiter << endl;
#endif

      IntervalCliqueFinder<IDType, UnitType> icf;
      list<IntervalClique<IDType, UnitType> > xCliques;
      icf.findCliques(xIntervals, xCliques);
      list<IntervalClique<IDType, UnitType> > yCliques;
      icf.findCliques(yIntervals, yCliques);
      
#ifdef DEBUG_ISCS
      icf.printCliques("Cliques in the X dimension", xCliques, cerr);
      icf.printCliques("Cliques in the Y dimension", yCliques, cerr);
#endif
      
      SetUpCliques(xCliques, yCliques);
    }
  
  IntervalSetCoverSolver(vector<Quadrilateral<IDType, UnitType> > & quads)
    {
      PopulateQuads(quads);
      
      vector<Quadrilateral<IDType, UnitType> >::iterator qiter;
      list<Interval<IDType, UnitType > > xIntervals;
      list<Interval<IDType, UnitType > > yIntervals;
      for(qiter = quads.begin(); qiter != quads.end(); qiter++)
      {
        Interval<IDType, UnitType> xinterval(qiter->getXProjection(),
                                             qiter->getID());
        xIntervals.push_back(xinterval);
        
        Interval<IDType, UnitType> yinterval(qiter->getYProjection(),
                                             qiter->getID());
        yIntervals.push_back(yinterval);
      }

      IntervalCliqueFinder<IDType, UnitType> icf;
      list<IntervalClique<IDType, UnitType> > xCliques;
      icf.findCliques(xIntervals, xCliques);
      list<IntervalClique<IDType, UnitType> > yCliques;
      icf.findCliques(yIntervals, yCliques);
      
#ifdef DEBUG_ISCS
      icf.printCliques("Cliques in the X dimension", xCliques, cerr);
      icf.printCliques("Cliques in the Y dimension", yCliques, cerr);
#endif
      SetUpCliques(xCliques, yCliques);
    }
  
  IntervalSetCoverSolver(vector<Quadrilateral<IDType, UnitType> > & quads,
                         list<IntervalClique<IDType, UnitType> > & xCliques,
                         list<IntervalClique<IDType, UnitType> > & yCliques)
    {
      PopulateQuads(quads);
      SetUpCliques(xCliques, yCliques);
    }

  /*
    clique IDs in pxes and pyTwoDIs are sorted low to high

    find commonIDs between all x & y cliques

    while there are still remaining IDs not assigned to xyCliques,
      find x,y pair with largest number of common IDs (n)
        make it an xyClique
        decrement number of remaining IDs by n
        remove these IDs from all remaining x,y pair cliques
        delete any x,y clique pairs with 0 common IDs
  */
  void solve(vector<TwoDIntervalClique<IDType, UnitType> > & twoDIs)
    {
      int ix, iy;
      int numIDsRemaining = pquads.size();
      int lpNumIDs = 0;
      list<IntervalClique<IDType, UnitType> >::iterator xiter;
      list<TwoDIntervalClique<IDType, UnitType> >::iterator yiter;
      list<TwoDIntervalClique<IDType, UnitType> >::iterator ysiter;
      list<CliquePairIntersection<IDType, UnitType> >::iterator cpiIter;
      list<CliquePairIntersection<IDType, UnitType> >::iterator cpiIter2;

#ifdef DEBUG_ISCS
      cerr << "X cliques sorted left to right (X dimension end)\n";
      for(ix = 0, xiter = pxes.begin();
          xiter != pxes.end();
          xiter++, ix++)
        cerr << ix << ": " << *xiter << endl;
      
      cerr << "\nY cliques sorted left to right (X dimension end)\n";
      for(iy = 0, yiter = pyTwoDIs.begin();
          yiter != pyTwoDIs.end();
          yiter++, iy++)
        cerr << iy << ": " << *yiter << endl;
#endif
      
#if defined(FEEDBACK_ISCS) || defined(DEBUG_ISCS)
      cerr << "Finding common IDs in x, y clique pairs\n";
#endif
      
      // find extent of each y clique in x dimension
      // and sort left to right to make x/y comparison faster
      
      // find common ids between x cliques and y cliques
      // keep track of the pair with the most IDs in common
      list<CliquePairIntersection<IDType, UnitType> > cpis;
      ysiter = pyTwoDIs.begin();
      for(ix = 0, xiter = pxes.begin();
          xiter != pxes.end();
          xiter++, ix++)
      {
        for(iy = 0, yiter = ysiter;
            yiter != pyTwoDIs.end();
            yiter++, iy++)
        {
          CliquePairIntersection<IDType, UnitType> cpi;
          list<IDType> commonIDs;

          if(xiter->getIntervalMin() > yiter->getXIntervalMax() ||
             xiter->getIntervalMax() < yiter->getXIntervalMin())
            continue;

          // step through x IDs & y IDs to find common ones
          // by definition, they are already sorted
          xiter->getCommonIDs(yiter->getIDs(), commonIDs);
          if(commonIDs.size() != 0)
          {
#ifdef DEBUG_ISCS
            cerr << ix << ":" << iy << ": ";
            list<IDType>::iterator tempI;
            for(tempI = commonIDs.begin();
                tempI != commonIDs.end();
                tempI++)
              cerr << *tempI << " ";
            cerr << endl;
#endif
            cpi.setCommonIDs(commonIDs);
            cpi.setXIndex(ix);
            cpi.setYIndex(iy);
            cpi.setXInterval(yiter->getXInterval());
            cpis.push_back(cpi);
          }
        }
      }

#ifdef DELETE_DUPLICATES
#if defined(FEEDBACK_ISCS) || defined(DEBUG_ISCS)
      cerr << "\nDeleting duplicate xy pairs\n";
#endif

      cpis.sort();
      cpiIter = cpis.begin();
      cpiIter2 = cpiIter;
      int numDeleted = 0;
      for(cpiIter2++;cpiIter2 != cpis.end();)
      {
        if(*cpiIter == *cpiIter2)
        {
          cpiIter2 = cpis.erase(cpiIter2);
          numDeleted++;
        }
        else
        {
          cpiIter++;
          cpiIter2++;
        }
      }

#if defined(FEEDBACK_ISCS) || defined(DEBUG_ISCS)
      cerr << "\nDeleted " << numDeleted << " duplicates\n";
#endif
#endif
      
#if defined(FEEDBACK_ISCS) || defined(DEBUG_ISCS)
      cerr << "\nPerforming greedy selection\n";
#endif

      // find first, largest xy pair
      CliquePairIntersection<IDType, UnitType> largestCPI;
      for(lpNumIDs = 0, cpiIter2 = cpis.begin();
          cpiIter2 != cpis.end();
          cpiIter2++)
      {
        if(cpiIter2->getNumCommonIDs() > lpNumIDs)
        {
          largestCPI = *cpiIter2;
          lpNumIDs = cpiIter2->getNumCommonIDs();
        }
      }
      
      // now greedily select xy pairs
      list<IDType>::const_iterator tdiiter;
      do
      {
        /*
          already have the index of the next xy pair
          create new 2d interval clique from it
          
          decrement numIDsRemaining by the number of commonIDs
          if(numIDsRemaining > 0)
            go through all xy pairs & delete commonIDs, and
            select next xy pair with the most commonIDs
            
          find the xy pair with the most commonIDs
          make it the next xy clique
        */
        TwoDIntervalClique<IDType, UnitType> twoDI(largestCPI,
                                                   pquads, id2Quad);
        numIDsRemaining -= lpNumIDs;

#ifdef DEBUG_ISCS
        cerr << largestCPI << endl;
        cerr << numIDsRemaining << " remaining\n";
#endif
        
        twoDIs.push_back(twoDI);

        // update the list of cliqued IDs
        const list<IDType> tdiids = twoDI.getIDs();
        for(tdiiter = tdiids.begin(); tdiiter != tdiids.end(); tdiiter++)
          idCliqued[*tdiiter] = true;
        
        // update the xy pairs & find the next largest one
        lpNumIDs = 0;
        for(cpiIter = cpis.begin(); cpiIter != cpis.end();)
        {
          if(twoDI.getXIntervalMax() > cpiIter->getXIntervalMin() &&
             twoDI.getXIntervalMin() < cpiIter->getXIntervalMax())
          {
#ifdef DEBUG_ISCS
            cerr << "Deleting common IDs in xy pair "
                 << *tdiiter << ", "
                 << cpiIter->getNumCommonIDs() << endl;
#endif
            cpiIter->deleteCommonIDs(tdiids);

            if(cpiIter->getNumCommonIDs() < 2)
            {
              cpiIter = cpis.erase(cpiIter);
              continue;
            }
          }

          if(cpiIter->getNumCommonIDs() > lpNumIDs)
          {
            largestCPI = *cpiIter;
            lpNumIDs = cpiIter->getNumCommonIDs();
          }
          cpiIter++;
        }
      } while(numIDsRemaining > 0 && cpis.size() > 0 && lpNumIDs > 0);

      // create cliques for uncliqued quadrilaterals
      vector<Quadrilateral<IDType, UnitType> >::iterator qiter;
      for(qiter = pquads.begin(); qiter != pquads.end(); qiter++)
      {
        if(!idCliqued[qiter->getID()])
        {
          TwoDIntervalClique<IDType, UnitType> twoDI(*qiter);
          twoDIs.push_back(twoDI);
        }
      }
    }
  
private:
  void PopulateQuads(const vector<Rectangle<IDType, UnitType> > & rects)
    {
      for(unsigned int i = 0; i < rects.size(); i++)
      {
        Quadrilateral<IDType, UnitType> q(rects[i].getXMin(),
                                          rects[i].getYMin(),
                                          rects[i].getXMin(),
                                          rects[i].getYMax(),
                                          rects[i].getXMax(),
                                          rects[i].getYMax(),
                                          rects[i].getXMax(),
                                          rects[i].getYMin(),
                                          rects[i].getID());
        pquads.push_back(q);
        id2Quad[pquads[i].getID()] = i;
        idCliqued[pquads[i].getID()] = false;
      }
    }
  void PopulateQuads(const vector<Quadrilateral<IDType, UnitType> > & quads)
    {
      pquads = quads;
      for(unsigned int i = 0; i < pquads.size(); i++)
      {
        id2Quad[pquads[i].getID()] = i;
        idCliqued[pquads[i].getID()] = false;
      }
    }

  void SetUpCliques(const list<IntervalClique<IDType, UnitType> > & xCliques,
                    const list<IntervalClique<IDType, UnitType> > & yCliques)
    {
      list<IntervalClique<IDType, UnitType> >::const_iterator iter;

      pxes = xCliques;
      pxes.sort();
                           
      /*
      // convert x Cliques to twoDIs
      for(iter = xCliques.begin(); iter != xCliques.end(); iter++)
      {
        if(iter->getNumIDs() > 1)
        {
          TwoDIntervalClique<IDType, UnitType> twoDI(*iter, true,
                                                     pquads, id2Quad);
          pxTwoDIs.push_back(twoDI);
        }
      }
      // sort on x interval max
      pxTwoDIs.sort();
      */
      
      // convert y Cliques to twoDIs
      for(iter = yCliques.begin(); iter != yCliques.end(); iter++)
      {
        if(iter->getNumIDs() > 1)
        {
          TwoDIntervalClique<IDType, UnitType> twoDI(*iter, false,
                                                     pquads, id2Quad);
          pyTwoDIs.push_back(twoDI);
        }
      }
      // sort on x interval max
      pyTwoDIs.sort();
    }
      
  vector<Quadrilateral<IDType, UnitType> > pquads;
  map<IDType, int> id2Quad;
  map<IDType, bool> idCliqued;
  // list<TwoDIntervalClique<IDType, UnitType> > pxes;
  list<IntervalClique<IDType, UnitType> > pxes;
  list<TwoDIntervalClique<IDType, UnitType> > pyTwoDIs;
};

#endif
