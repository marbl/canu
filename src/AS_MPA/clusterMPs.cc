
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
/* $Id: clusterMPs.cc,v 1.1.1.1 2004-04-14 13:51:58 catmandew Exp $ */
#include <iostream>
#include <vector>
#include <list>
#include <map>

#include "MPTypes.h"
#include "Rectangle.h"
#include "TwoDIntervalClique.h"
#include "IntervalSetCoverSolver.h"
#include "CompositeMPPolygon.h"
#include "PolygonIntersector.h"
#include "clusterMPs.h"

//#define DEBUG_CLUSTERMPS
//#define DEBUG_CLUSTERMPS_BIGTIME

void ReadCloneLibs(map<ID_TYPE, CloneLibrary> & libs, ifstream & fin)
{
  char line[4096];
  CloneLibrary cl;
  while(fin.getline(line, 4095))
  {
    cl.setFromString(line);
    libs[cl.getUID()] = cl;
  }
}


void ReadMatePairs(vector<MatePair> & mps, ifstream & fin)
{
  char line[4096];
  MatePair mp;
  
  while(fin.getline(line, 4095))
  {
    if(line[0] == 'I' || line[0] == 'A' || line[0] == 'N' || line[0] == 'O')
    {
      mp.setFromString(line);
      mps.push_back(mp);
    }
  }
}


void MakeUnionsOfIntersectingRectangles(list<TwoDIntervalClique<int, UNIT_TYPE> > & tdics,
                                           list<Rectangle<int, UNIT_TYPE> > & rects)
{
  // for each un-grouped rectangle
  list<Rectangle<int, UNIT_TYPE> >::iterator mbrIter1;
  for(mbrIter1 = rects.begin(); mbrIter1 != rects.end();)
  {
    // seed a new twoDI & delete the rectangle
    TwoDIntervalClique<int, UNIT_TYPE> tdic(*mbrIter1);
    mbrIter1 = rects.erase(mbrIter1);
    
    // compare with all other remaining rectangles
    list<Rectangle<int, UNIT_TYPE> >::iterator mbrIter2;
    mbrIter2 = mbrIter1;
    for(mbrIter2++; mbrIter2 != rects.end();)
    {
      // if intersection, add to the new tdic & delete the rectangle
      if(tdic.intersects(*mbrIter2))
      {
        tdic.addRectangle(*mbrIter2);
        mbrIter2 = rects.erase(mbrIter2);
      }
      else
        mbrIter2++;
    }
        
    // compare with all other tdics
    list<TwoDIntervalClique<int, UNIT_TYPE> >::iterator titer1;
    for(titer1 = tdics.begin(); titer1 != tdics.end();)
    {
      // if intersection, add to the new tdic & delete the old tdic
      if(tdic.intersects(*titer1))
      {
        tdic.add(*titer1);
        titer1 = tdics.erase(titer1);
      }
      else
        titer1++;
    }
    
    // add the new tdic to the set of tdics
    tdics.push_back(tdic);
  }
}

/*
  incoming - vector of mate pair polygons
  outgoing - vector of clustered mate pair polygons
             vector of problematic mate pair polygons
             vector of solitary mate pair polygons
 */

void ClusterMPPs(const vector<CompositeMPPolygon<UNIT_TYPE> > & mpps,
                 vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                 vector<CompositeMPPolygon<UNIT_TYPE> > & pmpps,
                 vector<CompositeMPPolygon<UNIT_TYPE> > & fmpps,
                 MatePairIndex_e mpii, unsigned int filterThresh)
{
  // populate a list of minimum bounding rectangles - one for each mp
  vector<Rectangle<int, UNIT_TYPE> > mpMBRsFixed;
  list<Rectangle<int, UNIT_TYPE> > mpMBRsSorted;
  for(unsigned int j = 0; j < mpps.size(); j++)
  {
    UNIT_TYPE minX, maxX, minY, maxY;
    mpps[j].getMBR(minX, maxX, minY, maxY);
    Rectangle<int, UNIT_TYPE> r(minX, maxX, minY, maxY, j);
    mpMBRsSorted.push_back(r);
    mpMBRsFixed.push_back(r);
  }
  mpMBRsSorted.sort();

  // find unions of intersecting MBRs
  list<TwoDIntervalClique<int, UNIT_TYPE> > tdics;
  MakeUnionsOfIntersectingRectangles(tdics, mpMBRsSorted);
  
#ifdef DEBUG_CLUSTERMPS
  cerr << "Looking for twoD interval cliques\n";
#endif
  
  // now, within each twoDI, find maximal cliques
  vector<TwoDIntervalClique<int, UNIT_TYPE> > mtdics;
  list<TwoDIntervalClique<int, UNIT_TYPE> >::iterator tdicsIter;
  int ti = 0;
  for(tdicsIter = tdics.begin(); tdicsIter != tdics.end(); tdicsIter++, ti++)
  {
    list<int> ids = tdicsIter->getIDs();
    list<int>::iterator idIter;
    
    if(tdicsIter->getNumIDs() < filterThresh)
    {
      // append associated polygon to vector of solitary polygons
      for(idIter = ids.begin(); idIter != ids.end(); idIter++)
        fmpps.push_back(mpps[*idIter]);
      continue;
    }

    // make a rectangle out of the two-D interval clique
    vector<Rectangle<int, UNIT_TYPE> > uoirects;
    for(idIter = ids.begin(); idIter != ids.end(); idIter++)
      uoirects.push_back(mpMBRsFixed[*idIter]);
    
#ifdef DEBUG_CLUSTERMPS_BIGTIME
    cerr << "Rectangle group " << ti << endl;
    vector<Rectangle<int, UNIT_TYPE> >::iterator uoirIter;
    for(uoirIter = uoirects.begin(); uoirIter != uoirects.end(); uoirIter++)
      cerr << *uoirIter << endl;
#endif
    
    // run the interval set cover solver on this union of intersecting rects
    IntervalSetCoverSolver<int, UNIT_TYPE> iscs(uoirects);
    mtdics.clear();
    iscs.solve(mtdics);
    
    // make sure the polygons intersect the set they're assigned to
    vector<TwoDIntervalClique<int, UNIT_TYPE> >::iterator mtdicsIter;
    PolygonIntersector<UNIT_TYPE> pi;
    for(mtdicsIter = mtdics.begin(); mtdicsIter != mtdics.end(); mtdicsIter++)
    {
      ids = mtdicsIter->getIDs();
      idIter = ids.begin();

      if(mtdicsIter->getNumIDs() < filterThresh)
      {
        for(idIter = ids.begin(); idIter != ids.end(); idIter++)
          fmpps.push_back(mpps[*idIter]);
        continue;
      }

      // see if the polygons intersect
#ifdef DEBUG_CLUSTERMPS
      cerr << "Possibly " << mtdicsIter->getNumIDs() << " intersections\n";
#endif
      bool allIntersect = true;
      CompositeMPPolygon<UNIT_TYPE> polyI(mpps[*idIter]);
      vector<CompositeMPPolygon<UNIT_TYPE> > bps;
      for(++idIter; idIter != ids.end(); idIter++)
      {
        if(!pi.polygonsIntersect(polyI, mpps[*idIter]))
        {
          bool intersectionFound = false;
          // see if it intersects a bp
          for(unsigned int bpi = 0; bpi < bps.size(); bpi++)
          {
            if(pi.polygonsIntersect(bps[bpi], mpps[*idIter]))
            {
              CompositeMPPolygon<UNIT_TYPE> mpp(bps[bpi]);
              pi.intersectPolygons(mpp, mpps[*idIter], bps[bpi]);
              bps[bpi].appendMP(mpps[*idIter].getMP(0));
              intersectionFound = true;
              break;
            }
          }
          if(!intersectionFound)
            bps.push_back(mpps[*idIter]);
          allIntersect = false;
          // break;
        }
        else
        {
          CompositeMPPolygon<UNIT_TYPE> mpp(polyI);
          pi.intersectPolygons(mpp, mpps[*idIter], polyI);
        }
      }

      // if there is a list of non-intersecting polygon ids, try to
      // intersect these
      if(allIntersect)
      {
        polyI.setType((MatePairIndex_e) mpii);
        idIter = ids.begin();
        for(++idIter; idIter != ids.end(); idIter++)
          polyI.appendMP(mpps[*idIter].getMP(0));
        cmpps.push_back(polyI);
      }
      else
      {
        for(unsigned int bpi = 0; bpi < bps.size(); bpi++)
        {
          if(bps[bpi].getNumMPs() > 1)
            cmpps.push_back(bps[bpi]);
          else
            pmpps.push_back(bps[bpi]);
        }
      }
    }
  }
}
  

void DetectTranspositions(const vector<CompositeMPPolygon<UNIT_TYPE> > & compressed,
                          const vector<CompositeMPPolygon<UNIT_TYPE> > & stretched,
                          const vector<CompositeMPPolygon<UNIT_TYPE> > & outties,
                          vector<CompositeMPPolygon<UNIT_TYPE> > & trans,
                          map<ID_TYPE, CloneLibrary> & libs,
                          double numStddevs,
                          unsigned int filterThresh)
{
  /*
    Each outtie implies 3 (adjacent) or 4 (separated) breakpoints - assume 4.
    Process:
      1. Intersect outties & keep only confirmed ones (>=2 mate pairs).
      2. Use stretched to refine left- and right-most breakpoints and
         identify & refine inner-left & inner-right breakpoints
         stretched must have either
         A. right-most read in transposition interval and
            left-most read to the left of transposition and
            left breakpoints intervals intersect
            stretched in such a way that when satisfied, the right-most read
            will still be in the transposition interval OR
         B. left-most read in transposition interval and
            right-most read to the right of the transposition and
            right breakpoint intervals intersect
            stretched in such a way that when satisfied, the left-most read
            will still be in the transposition interval
      3. Use compressed to refine left- and right-most breakpoints
         Similar requirements as stretched - one read must be outside the
         transposition and one must be inside, and the latter must still
         be inside when 'made' satisfied
      
    Or Rule: any single outtie must be corroborated by at least
      one compressed or stretched to the left & at least one
      to the right, given the above descriptions
   */
  trans.clear();
  for(unsigned int i = 0; i < outties.size(); i++)
  {
    CompositeMPPolygon<UNIT_TYPE> mpp(outties[i]);
    
    {
      int numLeftStretched = 0;
      int numRightStretched = 0;
      int numLeftCompressed = 0;
      int numRightCompressed = 0;
      

      /*

         A           B              C              D             E
      -------|--------------|--------------|----------------|----------

   O                -->                        <--
   S1  -->       <--
   S2                                              -->          <--
   C1  -->                                       <--
   C2                -->                                          <--

         A           D                C              B           E
      -------|----------------|--------------|--------------|----------
   O              <--                              -->
   S1  -->                                      <--
   S2                 -->                                       <--
   C1  -->         <--
   C2                                               -->           <--

       */

      // look at all the stretched mate pairs 
      for(unsigned int j = 0; j < stretched.size(); j++)
      {
        double lbound = libs[stretched[j].getMP(0).getLibUID()].getMean() -
          numStddevs * libs[stretched[j].getMP(0).getLibUID()].getStddev();
        double hbound = libs[stretched[j].getMP(0).getLibUID()].getMean() +
          numStddevs * libs[stretched[j].getMP(0).getLibUID()].getStddev();

        if(// left breakpoint intervals intersect
           mpp.getMinX() < stretched[j].getMaxX() &&
           mpp.getMaxX() > stretched[j].getMinX() &&
           // right stretched read is between left outtie read and
           // rightmost end of right breakpoint interval of outtie
           mpp.getMaxX() < stretched[j].getMaxY() &&
           mpp.getMaxY() > stretched[j].getMaxY())
        {
          numLeftStretched++;
          if(numLeftStretched == 1)
            for(unsigned int k = 0; k < stretched[j].getNumMPs(); k++)
              mpp.appendMP(stretched[j].getMP(k));
        }
        else if(// right breakpoint intervals intersect
          mpp.getMaxY() > stretched[j].getMinY() &&
          mpp.getMinY() < stretched[j].getMaxY() &&
          // left stretched read is between right outtie read and
          // leftmost end of left breakpoint interval of outtie
          mpp.getMinY() > stretched[j].getMinX() &&
          mpp.getMinX() < stretched[j].getMinX())
        {
          numRightStretched++;
          if(numRightStretched == 1)
            for(unsigned int k = 0; k < stretched[j].getNumMPs(); k++)
              mpp.appendMP(stretched[j].getMP(k));
        }
      }
    /*
      look for a compressed mate pair off the left or right end
      for the left end, the interval between frags in the compressed matepair
      must intersect the left breakpoint interval of the outtie AND
      the interval where the right fragment of the compressed pair should be
      must be at least partially to the right of the left fragment of the
      outtie pair AND if the compressed pair were made satisfied, the read
      in the transposition before would have to be in it after, too.
     */
    mpp.setType(MPI_TRANSPOSITION);
    for(unsigned int j = 0; j < compressed.size(); j++)
    {
      if(// compressed breakpoint interval intersects left outtie bpt intervavl
        mpp.getMaxX() > compressed[j].getMinX() &&
        mpp.getMinX() < compressed[j].getMaxX() &&
        // right compressed read will end up to right of left outtie read
        mpp.getMaxX() < compressed[j].getMaxX() + compressed[j].getMaxY() &&
        // right compressed read will end up to left of right outtie breakpoint
        mpp.getMaxY() > compressed[j].getMaxX() + compressed[j].getMinY())
      {
        numLeftCompressed++;
        for(unsigned int k = 0; k < compressed[j].getNumMPs(); k++)
          mpp.appendMP(compressed[j].getMP(k));
      }
      else if(// compressed bpt interval intersects right outtie bpt interval
        mpp.getMaxY() > compressed[j].getMinX() &&
        mpp.getMinY() < compressed[j].getMaxX() &&
        // left compressed read will end up to left of right outtie read
        mpp.getMinY() > compressed[j].getMinX() - compressed[j].getMaxY() &&
        // left compressed read will end up to right of left outtie read
        mpp.getMaxX() < compressed[j].getMinX() - compressed[j].getMinY())
      {
        numRightCompressed++;
        for(unsigned int k = 0; k < compressed[j].getNumMPs(); k++)
          mpp.appendMP(compressed[j].getMP(k));
      }
    }

      // require sufficient evidence (off right & left end of outtie)
#ifdef NEVER
      if(numLeftCompressed + numLeftStretched > 0 &&
         numRightCompressed + numRightStretched > 0)
#else
      if(mpp.getNumMPs() >= filterThresh)
        trans.push_back(mpp);
#endif
    }
  }
}

/*
  Recursive function to find index of first mate pair in interval
  (starting at position x)
 */
unsigned int FindFirstMatePair(UNIT_TYPE x,
                               const vector<MatePair> & smpsv,
                               int level, int index, int maxIndex)
{
  if(index <= 0)
    return 0;
  if(index >= maxIndex - 1)
    return (maxIndex - 1);
  if(smpsv[index].getLeftCoord() < x && smpsv[index+1].getLeftCoord() >= x)
    return index+1;
  
  if(smpsv[index].getLeftCoord() == x)
  {
    for(; index >= 0; index--)
      if(smpsv[index].getLeftCoord() < x) return index+1;
    return 0;
  }
  
  level = (level * 2 > maxIndex) ? maxIndex : (level * 2);
  if(smpsv[index].getLeftCoord() > x)
  {
    return FindFirstMatePair(x, smpsv, level,
                             index - maxIndex / level, maxIndex);
  }
  else
    return FindFirstMatePair(x, smpsv, level,
                             index + maxIndex / level, maxIndex);
}


bool RefineIntervalWithSatisfied(CompositeMPPolygon<UNIT_TYPE> & mpp,
                                 UNIT_TYPE & minLeft,
                                 UNIT_TYPE & maxLeft,
                                 UNIT_TYPE & minRight,
                                 UNIT_TYPE & maxRight,
                                 const vector<MatePair> & smpsv)
{
  UNIT_TYPE nl, xl, nr, xr;
  unsigned int si;

  nl = minLeft; xl = maxLeft; nr = minRight; xr = maxRight;
  
  // loop through satisfied mate pairs with a fragment in the interval
  /*
    there are 11 relevant combinations of left/right positions that
    can refine or invalidate an interval

          minLeft   maxLeft               minRight   maxRight
    ----------|-------|----------------------|---------|------------
        1         2               3               4         5
    Define intervals:
        1=to the left of the minLeft
        2=between minLeft and maxLeft
        3=between maxLeft and minRight
        4=between minRight and maxRight
        5=to the right of maxRight
        
    Define interval pairs:
      interval of leftFrag, interval of rightFrag
      1,1: irrelevant
      1,2: set minLeft to thisRight + 1
           to refine the left boundary to be to the right of the matepair
      1,3: invalidates the interval
      1,4: set maxRight to thisRight - 1
           to refine the right boundary to be to the left of the right frag
      1,5: irrelevant
      2,2: set minLeft to thisRight + 1 OR set maxLeft to thisLeft - 1
      2,3: set maxLeft to thisLeft - 1
      2,4: adjust so both fragments are within interval or outside it
      2,5: set minLeft to thisLeft + 1
      3,3: okay
      3,4: set minRight to thisRight + 1
      3,5: invalidates interval
      4,4: set minRight to thisRight + 1 OR set maxRight to thisLeft - 1
      4,5: set maxRight to thisLeft - 1
      5,5: irrelevant

      NOTE: given that intervals 2, 3, and 4 may overlap, these
        alterntatives are not mutually exclusive. However, for simplicity
        they will be treated as if they are, in the order listed above
  */
  si = FindFirstMatePair(minLeft, smpsv, 2, smpsv.size() / 2, smpsv.size());
  assert(si == 0 || si == smpsv.size() - 1 ||
         (smpsv[si-1].getLeftCoord() < minLeft &&
          smpsv[si].getLeftCoord() <= minLeft));
  for(;;si++)
  {
    UNIT_TYPE thisLeft, thisRight;
    if(smpsv[si].getLeftCoord() < smpsv[si].getRightCoord())
    {
      thisLeft = smpsv[si].getLeftCoord();
      thisRight = smpsv[si].getRightCoord();
    }
    else
    {
      thisLeft = smpsv[si].getRightCoord();
      thisRight = smpsv[si].getLeftCoord();
    }

    if(thisLeft < minLeft)
    {
      if(thisRight < minLeft || thisRight > maxRight)
        continue;

      if(thisRight <= maxLeft)
      {
        // 1,2: set minLeft to thisRight + 1
        minLeft = thisRight + 1;
      }
      else if(thisRight <= minRight)
      {
        // 1,3: invalidates the interval - but let it go for now
        continue;
      }
      else
      {
        // 1,4: set maxRight to thisRight - 1
        assert(thisRight <= maxRight);
        maxRight = thisRight - 1;
      }
    }
    else if(thisLeft < maxLeft)
    {
      if(thisRight <= maxLeft)
      {
        // 2,2: set minLeft to thisRight + 1 OR set maxLeft to thisLeft - 1
        // go with whichever refines the least
        if(thisRight - minLeft < maxLeft - thisLeft)
        {
          // adjust so matepair is to the left of the interval
          minLeft = thisRight + 1;
        }
        else
        {
          // adjust so matepair is within the interval
          maxLeft = thisLeft - 1;
        }
      }
      else if(thisRight <= minRight)
      {
        // 2,3: set maxLeft to thisLeft - 1
        maxLeft = thisLeft - 1;
      }
      else if(thisRight <= maxRight)
      {
        // 2,4: adjust so both fragments are within interval or outside it
        // go with whichever refines the least
        if((thisLeft - minLeft) + (maxRight - thisRight) <
           (maxLeft - thisLeft) + (thisRight - minRight))
        {
          // adjust so the mate pair straddles the interval
          minLeft = thisLeft + 1;
          maxRight = thisRight - 1;
        }
        else
        {
          // adjust so the mate pair is within the interval
          maxLeft = thisLeft - 1;
          minRight = thisRight + 1;
        }
      }
      else
      {
        // 2,5: set minLeft to thisLeft + 1
        assert(thisRight > maxRight);
        minLeft = thisLeft + 1;
      }
    }
    else if(thisLeft < minRight)
    {
      if(thisRight <= minRight)
        continue;

      if(thisRight <= maxRight)
      {
        // 3,4: set minRight to thisRight + 1
        minRight = thisRight + 1;
      }
      else
      {
        // 3,5: invalidates interval - but let it go for now
        assert(thisRight > maxRight);
        continue;
      }
    }
    else if(thisLeft < maxRight)
    {
      if(thisRight <= maxRight)
      {
        // 4,4: set minRight to thisRight + 1 OR set maxRight to thisLeft - 1
        // go with whichever refines the least
        if(thisRight - minRight < maxRight - thisLeft)
        {
          // adjust so matepair is within interval
          minRight = thisRight + 1;
        }
        else
        {
          // adjust so matepair is to the right of interval
          maxRight = thisLeft - 1;
        }
      }
      else
      {
        // 4,5: set maxRight to thisLeft - 1
        assert(thisRight > maxRight);
        maxRight = thisLeft - 1;
      }
    }
    else
    {
      // 5,5: irrelevant
      assert(thisLeft >= maxRight);
      break;
    }
  }
  
  if(minLeft > maxLeft || minRight > maxRight)
  {
    cerr << "Interval refinement invalidated interval:\n";
    cerr << " poly: " << mpp << endl;
    cerr << "  Initial ranges: "
         << nl << "," << xl << " - "
         << nr << "," << xr << endl;
    cerr << "  Invalid ranges: " << minLeft << "," << maxLeft << " - "
         << minRight << "," << maxRight << endl;
    return false;
  }

  if(minLeft != nl || maxLeft != xl || minRight != nr || maxRight != xr)
  {
    cerr << "Interval refinement refined interval:\n";
    cerr << "  Initial ranges: "
         << nl << "," << xl << " - "
         << nr << "," << xr << endl;
    cerr << "  Refined ranges: "
         << minLeft << "," << maxLeft << " - "
         << minRight << "," << maxRight << endl;

    Polygon<UNIT_TYPE> refiner;
    PolygonIntersector<UNIT_TYPE> pi;
    refiner.append(minLeft, minRight);
    refiner.append(minLeft, maxRight);
    refiner.append(maxLeft, maxRight);
    refiner.append(maxLeft, minRight);
    if(!pi.polygonsIntersect(mpp, refiner))
    {
      cerr << "Refined interval doesn't intersect original!\n";
      cerr << "  original polygon: " << mpp << endl;
      cerr << "        refinement: " << refiner << endl;
      return false;
    }
    
    CompositeMPPolygon<UNIT_TYPE> tempMPP;
    pi.intersectPolygons(mpp, refiner, tempMPP);
    mpp.append(tempMPP);
  }
  return true;
}

/*
  satisfied mate pair defines 5 regions
  if satisfied left,right intersects CP's left, right,
  at least one intersection of 5 regions with CP must be non-empty
  
    |      |       |
    |  1   |  bad  |  2
  | |------+-------+-------------
  v |      |       |
    | bad  |   3   |  bad
  ^ |      |       |
  | |------+-------+-------------
    |      |       |
    |  4   |  bad  |  5 (bad)
    |      |       |
    -----------------------------
           -->   <--
    1. unsatisfied mate pair straddles satisfied mate pair
    2. unsatisfied mate pair is to the right of satisfied mate pair
    3. satisfied mate pair straddles unsatisfied mate pair
    4. unsatisfied mate pair is to the left of satisfied mate pair
    5. left breakpoint is to the right of right breakpoint...
 */
bool RefineCPWithSatisfied(CompositeMPPolygon<UNIT_TYPE> & cmpp,
                           UNIT_TYPE & minLeft,
                           UNIT_TYPE & maxLeft,
                           UNIT_TYPE & minRight,
                           UNIT_TYPE & maxRight,
                           const vector<MatePair> & smpsv)
{
  unsigned int si;

  if(smpsv.size() == 0) return true;
  
  si = FindFirstMatePair(minLeft, smpsv, 2, smpsv.size() / 2, smpsv.size());
  assert(si == 0 || si == smpsv.size() - 1 ||
         (smpsv[si-1].getLeftCoord() < minLeft &&
          smpsv[si].getLeftCoord() >= minLeft));
  
  PolygonIntersector<UNIT_TYPE> pi;
  for(;si < smpsv.size();si++)
  {
    if(smpsv[si].getLeftCoord() >= maxRight &&
       smpsv[si].getRightCoord() >= maxRight)
      break;
    if(smpsv[si].getLeftCoord() <= minLeft &&
       smpsv[si].getRightCoord() <= minLeft)
      continue;

    UNIT_TYPE left, right;
    if(smpsv[si].getLeftCoord() < smpsv[si].getRightCoord())
    {
      // satisfied in as-is form
      left = smpsv[si].getLeftCoord();
      right = smpsv[si].getRightCoord();
    }
    else
    {
      // satisfied in outtie form
      right = smpsv[si].getLeftCoord();
      left = smpsv[si].getRightCoord();
      
      // check that we haven't already processed this one in innie form
      if(left >= minLeft)
        continue;
    }
    
    // define 5 polygons to intersect
    CompositePolygon<UNIT_TYPE> refiner;
    Polygon<UNIT_TYPE> sp;

    // center
    sp.reset();
    sp.append(left, left);
    sp.append(left, right);
    sp.append(right, right);
    sp.append(right, left);
    refiner.append(sp);

#define NUDGE_BP  1
    
    // upper left
    if(maxRight > right && minLeft < left)
    {
      sp.reset();
      sp.append(minLeft - NUDGE_BP, right);    // lower left
      sp.append(minLeft - NUDGE_BP, maxRight + NUDGE_BP); // upper left
      sp.append(left, maxRight - NUDGE_BP);    // upper right
      sp.append(left, right);       // lower right
      refiner.append(sp);
    }

    // upper right
    if(maxRight > right && maxLeft > right)
    {
      sp.reset();
      sp.append(right, right);
      sp.append(right, maxRight + NUDGE_BP);
      sp.append(maxLeft + NUDGE_BP, maxRight + NUDGE_BP);
      sp.append(maxLeft + NUDGE_BP, right);
      refiner.append(sp);
    }

#ifdef NEVER
    // lower right - shouldn't need this one?
    if(minRight < left && maxLeft > right)
    {
      sp.reset();
      sp.append(right, minRight - NUDGE_BP);
      sp.append(right, left);
      sp.append(maxLeft + NUDGE_BP, left);
      sp.append(maxLeft + NUDGE_BP, minRight - NUDGE_BP);
      refiner.append(sp);
    }
#endif

    // lower left
    if(minRight < left && minLeft < left)
    {
      sp.reset();
      sp.append(minLeft - NUDGE_BP, minRight - NUDGE_BP);
      sp.append(minLeft - NUDGE_BP, left);
      sp.append(left, left);
      sp.append(left, minRight - NUDGE_BP);
      refiner.append(sp);
    }

    if(!pi.polygonsIntersect(cmpp, refiner))
    {
#ifdef DEBUG_CLUSTERMPS
      cerr << "Satisfied matepair invalidates bad matepair polygon\n";
      cerr << "  original polygon: " << cmpp << endl;
      cerr << "        refinement: " << refiner << endl;
#endif
      return false;
    }

    CompositeMPPolygon<UNIT_TYPE> tempMPP;
    pi.intersectPolygons(cmpp, refiner, tempMPP);
    cmpp.setPolygons(tempMPP);
    minLeft = cmpp.getMinX();
    maxRight = cmpp.getMaxY();
  }
  return true;
}


void RefineInversions(vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                     const vector<MatePair> & smpsv, // satisfied matepairs
                     unsigned int filterThresh)
{
  /*
    filter out inversions that are just normal or antinormal
    OR don't have the right properties:
    
    |--------------|------------------------|------------------|
             -->     <--            -->       <--
             rln     lla            rrn       lra
        rln = rightmost left normal coordinate
        lla = leftmost left antinormal coordinate
        rrn = rightmost right normal
        lra = leftmost right antinormal

        rln <= lla (bounds the left breakpoint of the inversion)
        rrn <= lra (bounds the right breakpoint of the inversion)
        lla may be to the right of rrn...

    left breakpoint must be to the left of the right breakpoint
    so, intersect the mpp with rectangle:
      (minY, minX), (minY, maxX), (maxY, maxX), (maxY, minX)
  */
  int numJustOneKind = 0;
  int numSatisfiedInvalidates = 0;
  unsigned int numTotal = cmpps.size();
  vector<CompositeMPPolygon<UNIT_TYPE> >::iterator mppIter;
  for(mppIter = cmpps.begin(); mppIter != cmpps.end();)
  {
    unsigned int numNormals = 0;
    unsigned int numAntinormals = 0;
    UNIT_TYPE rln, lla, lra, rrn;
    for(unsigned int i = 0; i < mppIter->getNumMPs(); i++)
    {
      const MatePair & mp = mppIter->getMP(i);
      if(mp.getOrientation() == PAIR_NORMAL)
      {
        if(numNormals == 0)
        {
          rln = mp.getLeftCoord();
          rrn = mp.getRightCoord();
        }
        else
        {
          rln = (rln < mp.getLeftCoord()) ? mp.getLeftCoord() : rln;
          rrn = (rrn < mp.getRightCoord()) ? mp.getRightCoord() : rrn;
        }
        numNormals++;
      }
      else
      {
        if(numAntinormals == 0)
        {
          lla = mp.getLeftCoord();
          lra = mp.getRightCoord();
        }
        else
        {
          lla = (lla > mp.getLeftCoord()) ? mp.getLeftCoord() : lla;
          lra = (lra > mp.getRightCoord()) ? mp.getRightCoord() : lra;
        }
        numAntinormals++;
      }
    }
#ifdef NEVER
    if(numNormals == 0 || numAntinormals == 0 ||
       numNormals + numAntinormals < filterThresh)
#else
    if(numNormals == 0)
    {
      rln = lla; rrn = lra;
    }
    if(numAntinormals == 0)
    {
      lla = rln; lra = rrn;
    }
    if(numNormals + numAntinormals < filterThresh)
#endif
    {
#ifdef DEBUG_CLUSTERMPS
      cerr << "Found non-inversion (one kind of unsatisfied):\n";
      cerr << " " << numNormals << " normals, "
           << numAntinormals << " antinormals\n";
#endif
      mppIter = cmpps.erase(mppIter);
      numJustOneKind++;
    }
    else
    {
      UNIT_TYPE xl = (lla < rrn) ? lla : rrn;
      UNIT_TYPE nr = (lla < rrn) ? rrn : lla;
      if(RefineCPWithSatisfied(*mppIter, rln, xl, nr, lra, smpsv))
      {
        mppIter++;
      }
      else
      {
#ifdef DEBUG_CLUSTERMPS
        cerr << "Incompatible original & refined inversions\n";
#endif
        mppIter = cmpps.erase(mppIter);
        numSatisfiedInvalidates++;
      }
    }
  }
  cerr << "Deleted "
       << numJustOneKind + numSatisfiedInvalidates
       << " of " << numTotal << " potential inversions.\n";
  cerr << "  " << numJustOneKind << " insufficient total or combination of types.\n";
  cerr << "  " << numSatisfiedInvalidates << " invalidated by satisfied matepairs.\n";
  cerr << "Kept " << cmpps.size() << " inversions.\n";
}


void RefineStretched(vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                     const vector<MatePair> & smpsv) // satisfied matepairs
{
  /*
    Stretched instances should not have one fragment inside the interval
    and the other fragment outside
  */
  int i;
  int numErased = 0;
  unsigned int numTotal = cmpps.size();
  vector<CompositeMPPolygon<UNIT_TYPE> >::iterator mppIter;
  for(i = 0, mppIter = cmpps.begin(); mppIter != cmpps.end();i++)
  {
    UNIT_TYPE minLeft = mppIter->getMinX();
    UNIT_TYPE maxLeft = mppIter->getMaxX();
    UNIT_TYPE minRight = mppIter->getMinY();
    UNIT_TYPE maxRight = mppIter->getMaxY();

    if(RefineCPWithSatisfied(*mppIter, minLeft, maxLeft,
                             minRight, maxRight, smpsv))
    {
      mppIter++;
    }
    else
    {
#ifdef DEBUG_CLUSTERMPS
      cerr << "Stretched interval disagrees with satisfied matepairs.\n";
#endif
      mppIter = cmpps.erase(mppIter);
      numErased++;
    }
  }
  cerr << "Deleted " << numErased << " of "
       << numTotal << " potential stretched matepairs invalidated by satisfied matepairs.\n";
}


/*
  if an mpp has >= minFilterCount matepairs & they're all in the
  same library, assume it's bad tracking in the lab
*/
void FilterMislabelledLibs(vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                           vector<CompositeMPPolygon<UNIT_TYPE> > & mmpps,
                           unsigned int minFilterCount)
{
  vector<CompositeMPPolygon<UNIT_TYPE> >::iterator mppIter;
  for(mppIter = cmpps.begin(); mppIter != cmpps.end();)
  {
    if(mppIter->getNumMPs() < minFilterCount)
    {
      mppIter++;
      continue;
    }

    uint64 libUID = mppIter->getMP(0).getLibUID();
    bool badMP = true;
    for(unsigned int j = 1; j < mppIter->getNumMPs(); j++)
    {
      if(mppIter->getMP(j).getLibUID() != libUID)
      {
        badMP = false;
        break;
      }
    }

    if(badMP)
    {
      mmpps.push_back((CompositeMPPolygon<UNIT_TYPE>) (*mppIter));
      mppIter = cmpps.erase(mppIter);
    }
    else
    {
      mppIter++;
    }
  }
}
