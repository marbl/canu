
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
/* $Id: PolygonIntersector.h,v 1.5 2008-06-27 06:29:16 brianwalenz Exp $ */
#ifndef POLYGONINTERSECTOR_H
#define POLYGONINTERSECTOR_H

#include <iostream>
#include <vector>

#include "Point.h"
#include "Polygon.h"
#include "ConvexHull.h"
#include "CompositePolygon.h"

//#define DEBUG_POLYGONINTERSECTOR

typedef enum
{
  PI_NOINTERSECTION,
  PI_1CONTAINS2,
  PI_2CONTAINS1,
  PI_INTERSECTION,
  PI_NUM_INTERSECTIONS
} PolygonIntersectionType;

template <class UnitType>
class PolygonIntersector
{
public:
  PolygonIntersector(){ppit = PI_NOINTERSECTION; pNumBridges = 0;}

  bool polygonsIntersect(const Polygon<UnitType> & poly0,
                         const Polygon<UnitType> & poly1)
    {
      reset();

      if(poly0.size() < 3 || poly1.size() < 3)
      {
        return PI_NOINTERSECTION;
      }

      // compute the convex hull of the two polygons
      cvh.create(poly0, poly1);

#ifdef DEBUG_POLYGONINTERSECTOR
      cerr << "Intersecting polygons:\n";
      cerr << "poly0: " << poly0 << endl;
      cerr << "poly1: " << poly1 << endl;
      cerr << "Convex hull polyPoints:\n";
      cerr << cvh << endl;
#endif

      // identify vertices upstream of intersection points
      processBridges(poly0, poly1, cvh);
      return(ppit != PI_NOINTERSECTION);
    }

  bool polygonsIntersect(const CompositePolygon<UnitType> & polys0,
                         const Polygon<UnitType> & poly1)
    {
      for(unsigned int i = 0; i < polys0.size(); i++)
        if(polygonsIntersect(polys0[i], poly1) == true) return true;
      return false;
    }
  bool polygonsIntersect(const CompositePolygon<UnitType> & polys0,
                         const CompositePolygon<UnitType> & polys1)
    {
      for(unsigned int i = 0; i < polys1.size(); i++)
        if(polygonsIntersect(polys0, polys1[i]) == true) return true;
      return false;
    }

  PolygonIntersectionType intersectPolygons(const CompositePolygon<UnitType> & cp0,
                                            const CompositePolygon<UnitType> & cp1,
                                            CompositePolygon<UnitType> & cpOut)
    {
      PolygonIntersectionType pit;
      cpOut.reset();
      // loop over cp0 polygons
      for(unsigned int i = 0; i < cp0.size(); i++)
      {
        for(unsigned int j = 0; j < cp1.size(); j++)
        {
          Polygon<UnitType> pInt;
          pit = intersectPolygons(cp0[i], cp1[j], pInt);
          if(pit != PI_NOINTERSECTION)
          {
            cpOut.append(pInt);
          }
        }
      }
      if(cpOut.size() == 0)
        return PI_NOINTERSECTION;
      return PI_INTERSECTION;
    }

  PolygonIntersectionType intersectPolygons(const Polygon<UnitType> & poly0,
                                            const Polygon<UnitType> & poly1,
                                            Polygon<UnitType> & pInt)
    {
      pInt.reset();
      if(polygonsIntersect(poly0, poly1))
      {
        switch(ppit)
        {
          case PI_1CONTAINS2:
            pInt = poly1;
            break;
          case PI_2CONTAINS1:
            pInt = poly0;
            break;
          case PI_INTERSECTION:
            buildIntersection(poly0, poly1, cvh, pInt);
            break;
          default:
            assert(ppit == PI_NOINTERSECTION);
            break;
        }
      }
      return ppit;
    }

private:
  void reset()
    {
      pKeyVs[0].clear();
      pKeyVs[1].clear();
      ppit = PI_NOINTERSECTION;
      cvh.reset();
      pNumBridges = 0;
    }

  void getIntersectionPoint(const Polygon<UnitType> & cwPoly,
                            const Polygon<UnitType> & ccwPoly,
                            int cwIndex, int ccwIndex,
                            Point<UnitType> & p)
    {
      int cwIndex2 = cwPoly.cwIndex(cwIndex);
      int ccwIndex2 = ccwPoly.ccwIndex(ccwIndex);

      // if segments are colinear, return preceding ccw point
      if(cwPoly[cwIndex].getLineSide(cwPoly[cwIndex2], ccwPoly[ccwIndex]) ==
         TP_COLINEAR)
      {
        /*
          assert(cwPoly[cwIndex].getLineSide(ccwPoly[ccwIndex2],
          ccwPoly[ccwIndex]) == TP_COLINEAR);
        */
        p = ccwPoly[ccwIndex2];
      }
      else
      {
        // return intersection between two line segments
        double cwM = cwPoly[cwIndex].getX() - cwPoly[cwIndex2].getX();
        double ccwM = ccwPoly[ccwIndex].getX() - ccwPoly[ccwIndex2].getX();
        if(cwM == 0)
        {
          // cw segment is vertical
          if(ccwM == 0)
          {
            // it's a point & TP_COLINEAR check didn't trigger
            cerr << "WARNING: Parallel vertical lines in convex hull bridge!\n";
            cerr << "     cwpoly: " << cwPoly << endl;
            cerr << "   hull pts: " << cwPoly[cwIndex]
                 << " & " << cwPoly[cwIndex2] << endl;
            cerr << "    ccwpoly: " << ccwPoly << endl;
            cerr << "   hull pts: " << ccwPoly[ccwIndex]
                 << " & " << ccwPoly[ccwIndex2] << endl;
            p.set((UnitType) (cwPoly[cwIndex].getX() +
                              ccwPoly[ccwIndex].getX()) / 2,
                  (UnitType) (cwPoly[cwIndex].getY() +
                              cwPoly[cwIndex2].getY() +
                              ccwPoly[ccwIndex].getY() +
                              ccwPoly[ccwIndex2].getY()) / 4);
            cerr << " setting intersection to : " << p << endl;
          }
          else
          {
            ccwM = (ccwPoly[ccwIndex].getY() -
                    ccwPoly[ccwIndex2].getY()) / ccwM;
            double ccwB = ccwPoly[ccwIndex].getY() -
              ccwM * ccwPoly[ccwIndex].getX();
            p.set((UnitType) cwPoly[cwIndex].getX(),
                  (UnitType) (ccwM * cwPoly[cwIndex].getX() + ccwB));
          }
        }
        else if(ccwM == 0)
        {
          // ccw segment is vertical
          assert(cwM != 0);
          cwM = (cwPoly[cwIndex].getY() - cwPoly[cwIndex2].getY()) / cwM;
          double cwB = cwPoly[cwIndex].getY() -
            cwM * cwPoly[cwIndex].getX();
          p.set((UnitType) ccwPoly[ccwIndex].getX(),
                (UnitType) (cwM * ccwPoly[ccwIndex].getX() + cwB));
        }
        else
        {
          cwM = (cwPoly[cwIndex].getY() - cwPoly[cwIndex2].getY()) / cwM;
          double cwB = cwPoly[cwIndex].getY() -
            cwM * cwPoly[cwIndex].getX();
          ccwM = (ccwPoly[ccwIndex].getY() - ccwPoly[ccwIndex2].getY()) / ccwM;
          double ccwB = ccwPoly[ccwIndex].getY() -
            ccwM * ccwPoly[ccwIndex].getX();
          if(cwM == ccwM)
          {
            cerr << "ERROR finding intersection between line segments!\n";
            cerr << "They must be parallel!\n";
            cerr << "cw segment: " << cwPoly[cwIndex]
                 << " to " << cwPoly[cwIndex2] << endl;
            cerr << "ccw segment: " << ccwPoly[ccwIndex]
                 << " to " << ccwPoly[ccwIndex2] << endl;
            cerr << "cw polygon: " << cwPoly << endl;
            cerr << "ccw polygon: " << ccwPoly << endl;
            assert(0);
          }
          double x = (ccwB - cwB)/(cwM - ccwM);
          if((UnitType) x == x)
          {
            p.setX((UnitType) x);
            p.setY((UnitType) (x * cwM + cwB));
          }
          else
          {
            p.setX((UnitType) (x + .5));
            p.setY((UnitType) (x * cwM + cwB + .5));
          }
        }
      }
    }

  void buildIntersection(const Polygon<UnitType> & poly0,
                         const Polygon<UnitType> & poly1,
                         const ConvexHull<UnitType> & cvh,
                         Polygon<UnitType> & pInt)
    {
      Point<UnitType> p;

      /*
      // sanity check that the first clockwise polygon should be poly0
      if(pKeyVs[0][0].getPolygon() != 0)
      {
      cerr << "Out of order problem with polygons!\n";
      cerr << "poly0: " << poly0 << endl;
      cerr << "poly1: " << poly1 << endl;
      cerr << "cvh: " << cvh << endl;
      assert(0);
      }
      */

#ifdef DEBUG_POLYGONINTERSECTOR
      cerr << "Key points:";
      for(unsigned int i = 0; i < pKeyVs[0].size(); i++)
        cerr << " " << pKeyVs[0][i] << ", " << pKeyVs[1][i];
      cerr << endl;
#endif

      const Polygon<UnitType> * polys[2];
      polys[0] = &poly0;
      polys[1] = &poly1;

      // pKeyVs[0] are clockwise polygon/indices
      // pKeyVs[1] are counter-clockwise polygon/indices

      pInt.reset();

      // add the first point
      // next polygon point is next key point
      getIntersectionPoint(*(polys[pKeyVs[0][0].getPolygon()]),
                           *(polys[pKeyVs[1][0].getPolygon()]),
                           pKeyVs[0][0].getIndex(),
                           pKeyVs[1][0].getIndex(),
                           p);
#ifdef DEBUG_POLYGONINTERSECTOR
      cerr << "Appending point to intersection: " << p << endl;
#endif
      if(!pInt.append(p))
      {
        cerr << "Failed to append point while building polygon intersection\n";
        cerr << *polys[0] << endl;
        cerr << *polys[1] << endl;
        cerr << cvh << endl;
        return;
      }

      unsigned int cpn = pKeyVs[0][0].getPolygon();
      unsigned int cpi = polys[cpn]->cwIndex(pKeyVs[0][0].getIndex());
      unsigned int cki = 1;

      // loop until we've 'looped'
      while(true)
      {
#ifdef DEBUG_POLYGONINTERSECTOR
        cerr << "cpn = " << cpn << endl;
        cerr << "cpi = " << cpi << endl;
        cerr << "cki = " << cki << endl;
#endif
        // stopping criterion - when next clockwise polygon point
        // is the first counter-clockwise key point
        if(cpn == pKeyVs[1][0].getPolygon() &&
           cpi == pKeyVs[1][0].getIndex())
          break;

        if(cki >= pKeyVs[0].size() || cpi != pKeyVs[1][cki].getIndex())
        {
          // no more key points or next polygon point is not next key point
#ifdef DEBUG_POLYGONINTERSECTOR
          cerr << "Appending point to intersection: polygon "
               << cpn << ", index " << cpi << " point: "
               << (*(polys[pKeyVs[0][cki].getPolygon()]))[cpi] << endl;
#endif
          if((*(polys[cpn]))[cpi] == pInt[0])
            break;
          if(!pInt.append((*(polys[cpn]))[cpi]))
          {
            cerr << "Failed to append point while building polygon intersection\n";
            cerr << *polys[0] << endl;
            cerr << *polys[1] << endl;
            cerr << cvh << endl;
            return;
          }
        }
        else
        {
          // next polygon point is next key point
          getIntersectionPoint(*(polys[pKeyVs[0][cki].getPolygon()]),
                               *(polys[pKeyVs[1][cki].getPolygon()]),
                               pKeyVs[0][cki].getIndex(),
                               pKeyVs[1][cki].getIndex(),
                               p);
#ifdef DEBUG_POLYGONINTERSECTOR
          cerr << "Appending point to intersection: " << p << endl;
          cerr << (*(polys[pKeyVs[0][cki].getPolygon()]))[pKeyVs[0][cki].getIndex()] << endl;
          cerr << (*(polys[pKeyVs[1][cki].getPolygon()]))[pKeyVs[1][cki].getIndex()] << endl;
          cerr << pKeyVs[0][cki].getIndex() << endl;
          cerr << pKeyVs[1][cki].getIndex() << endl;
#endif
          if(pInt[pInt.size()-1] != p)
          {
            if(!pInt.append(p))
            {
              cerr << "Failed to append point while building polygon intersection\n";
              cerr << *polys[0] << endl;
              cerr << *polys[1] << endl;
              cerr << cvh << endl;
              return;
            }
          }

#ifdef DEBUG_POLYGONINTERSECTOR
          cerr << "cpn was " << cpn << endl;
          cerr << "cpi was " << cpi << endl;
          cerr << "cki was " << cki << endl;
#endif
          cpn = 1 - cpn;
          cpi = pKeyVs[0][cki].getIndex();
          cki++;
#ifdef DEBUG_POLYGONINTERSECTOR
          cerr << "cpn is " << cpn << endl;
          cerr << "cpi is " << cpi << endl;
          cerr << "cki is " << cki << endl;
#endif
        }
        cpi = (*(polys[cpn])).cwIndex(cpi);
#ifdef DEBUG_POLYGONINTERSECTOR
        cerr << "cpi updated to " << cpi << endl;
#endif
      }
    }

  void getKeyVertices(const Polygon<UnitType> & cwPoly,
                      int cwIndex, int & cwKey,
                      const Polygon<UnitType> & ccwPoly,
                      int ccwIndex, int & ccwKey)
    {
      bool done = false;
      cwKey = cwIndex;
      ccwKey = ccwIndex;
      while(!done)
      {
        done = true;
        while(ccwPoly[ccwPoly.ccwIndex(ccwKey)].
              getLineSide(cwPoly[cwKey], cwPoly[cwPoly.cwIndex(cwKey)]) ==
              TP_COUNTERCW)
        {
          ccwKey = ccwPoly.ccwIndex(ccwKey);
          if(ccwKey == ccwIndex) return;
          done = false;
        }

        while(cwPoly[cwPoly.cwIndex(cwKey)].
              getLineSide(ccwPoly[ccwKey],
                          ccwPoly[ccwPoly.ccwIndex(ccwKey)]) ==
              TP_CLOCKWISE)
        {
          cwKey = cwPoly.cwIndex(cwKey);
          if(cwKey == cwIndex) return;
          done = false;
        }
      }
      ppit = PI_INTERSECTION;
    }

  /* handle special case if hull contains a gap bridge
     detects & adjust for this case:

           ccwPoly /   \     cwPoly
           --------     --------------
           |      |     |            |
       ccwIndex   |     |         cwIndex
              we want these
  */
  void adjustIfContainedGapBridge(const Polygon<UnitType> & cwPoly,
                                  int & cwi,
                                  const Polygon<UnitType> & ccwPoly,
                                  int & ccwi)
    {
      while(true)
      {
        ThreePointOrientation tpoCW2CCW =
          cwPoly[cwi].getLineSide(cwPoly[cwPoly.cwIndex(cwi)], ccwPoly[ccwi]);
        ThreePointOrientation tpoCCW2CW =
          ccwPoly[ccwi].getLineSide(ccwPoly[ccwPoly.ccwIndex(ccwi)], cwPoly[cwi]);

        double distCW2CW = cwPoly[cwi].distanceFrom(cwPoly[cwPoly.cwIndex(cwi)]);
        double distCW2CCW = cwPoly[cwi].distanceFrom(ccwPoly[ccwPoly.ccwIndex(ccwi)]);
        if(tpoCW2CCW != TP_COLINEAR ||
           tpoCCW2CW != TP_COLINEAR ||
           distCW2CW > distCW2CCW) break;

#ifdef DEBUG_POLYGONINTERSECTOR
        cerr << "Adjusting for contained gap bridge\n";
#endif
        /*
        // make sure that an endpoint of one segment isn't on the other segment
        bool onSeg1 = cwPoly[cwPoly.cwIndex(cwi)].isOnSegment(ccwPoly[ccwi],
                                                              ccwPoly[ccwPoly.ccwIndex(ccwi)]);
        bool onSeg2 = ccwPoly[ccwPoly.ccwIndex(ccwi)].isOnSegment(cwPoly[cwi],
                                                                  cwPoly[cwPoly.cwIndex(cwi)]);
        if(onSeg1 || onSeg2) break;
        */
        cwi = cwPoly.cwIndex(cwi);
        ccwi = ccwPoly.ccwIndex(ccwi);
      }
      /*
      while(cwPoly[cwPoly.cwIndex(cwi)].
            getLineSide(ccwPoly[ccwi],
                        ccwPoly[ccwPoly.ccwIndex(ccwi)]) == TP_COLINEAR &&
            cwPoly[cwi].distanceFrom(cwPoly[cwPoly.cwIndex(cwi)]) <=
            cwPoly[cwi].distanceFrom(ccwPoly[ccwPoly.ccwIndex(ccwi)]))
      {
#ifdef DEBUG_POLYGONINTERSECTOR
        cerr << "Adjusting for contained gap bridge\n";
#endif
        // make sure that an endpoint of one segment isn't on the other segment
        if(cwPoly[cwPoly.cwIndex(cwi)].isOnSegment(ccwPoly[ccwi],
                                   ccwPoly[ccwPoly.ccwIndex(ccwi)]) ||
           ccwPoly[ccwPoly.ccwIndex(ccwi)].isOnSegment(cwPoly[cwi],
                                   cwPoly[cwPoly.cwIndex(cwi)]))
          return;
        cwi = cwPoly.cwIndex(cwi);
        ccwi = ccwPoly.ccwIndex(ccwi);
      }
      */
    }

  /* handle special case if hull contains an overlap bridge
     detects & adjust for this case:

            ccwPoly \     /   cwPoly
           ----------=====-------------
           |        |    |            |
       ccwIndex     |    |        cwIndex
                we want these
    Note: this means that ccwPoly must become cwPoly & vice versa
  */
  bool checkForContainedOvlBridge(const Polygon<UnitType> & cwPoly,
                                  int & cwi,
                                  const Polygon<UnitType> & ccwPoly,
                                  int & ccwi)
    {
      bool found = false;
      if(cwPoly[cwPoly.cwIndex(cwi)].
         getLineSide(ccwPoly[ccwi],
                     ccwPoly[cwPoly.ccwIndex(ccwi)]) == TP_COLINEAR &&
         cwPoly[cwi].distanceFrom(ccwPoly[ccwPoly.ccwIndex(ccwi)]) <=
         cwPoly[cwi].distanceFrom(cwPoly[cwPoly.cwIndex(cwi)]))
      {
#ifdef DEBUG_POLYGONINTERSECTOR
        cerr << "Adjusting for contained overlapping bridge\n";
#endif
        cwi = cwPoly.cwIndex(cwi);
        ccwi = ccwPoly.ccwIndex(ccwi);
        found = true;
      }
      return found;
    }


  void processBridge(const Polygon<UnitType> & poly0,
                     const Polygon<UnitType> & poly1,
                     const ConvexHull<UnitType> & cvh,
                     bool poly0IsClockwise,
                     int cwIndex, int ccwIndex)
    {
      CVHPoint<UnitType> chp;
      int cwi, ccwi;

      pNumBridges++;

      if(poly0IsClockwise)
      {
        adjustIfContainedGapBridge(poly0, cwIndex, poly1, ccwIndex);
#ifdef NEVER
        if(checkForContainedOvlBridge(poly0, cwIndex, poly1, ccwIndex))
        {
          cwi = ccwIndex;
          ccwi = cwIndex;
          poly0IsClockwise = (poly0IsClockwise == true) ? false: true;
          ppit = PI_INTERSECTION;
        }
        else
#endif
        {
          getKeyVertices(poly0, cwIndex, cwi, poly1, ccwIndex, ccwi);
        }
      }
      else
      {
        adjustIfContainedGapBridge(poly1, cwIndex, poly0, ccwIndex);
#ifdef NEVER
        if(checkForContainedOvlBridge(poly1, cwIndex, poly0, ccwIndex))
        {
          cwi = ccwIndex;
          ccwi = cwIndex;
          poly0IsClockwise = (poly0IsClockwise == true) ? false: true;
          ppit = PI_INTERSECTION;
        }
        else
#endif
        {
          getKeyVertices(poly1, cwIndex, cwi, poly0, ccwIndex, ccwi);
        }
      }

      if(ppit != PI_NOINTERSECTION)
      {
        chp.set(0, 0, ((poly0IsClockwise) ? 0 : 1), cwi);
        pKeyVs[0].push_back(chp);
        chp.set(0, 0, ((poly0IsClockwise) ? 1 : 0), ccwi);
        pKeyVs[1].push_back(chp);
      }
    }

  void processBridges(const Polygon<UnitType> & poly0,
                      const Polygon<UnitType> & poly1,
                      const ConvexHull<UnitType> & cvh)
    {
      unsigned int i, j;

      // identify vertices upstream of intersection points based on bridges
      for(i = 0, j = 1; i < cvh.size(); i++, j = (j + 1) % cvh.size())
      {
        if(cvh[i].getPolygon() != cvh[j].getPolygon())
          processBridge(poly0, poly1, cvh,
                        (cvh[i].getPolygon() == 0) ? true : false,
                        cvh[i].getIndex(), cvh[j].getIndex());
      }

      if(pNumBridges == 0)
        ppit = (cvh[0].getPolygon() == 0) ? PI_1CONTAINS2 : PI_2CONTAINS1;
      else if(pKeyVs[0].size() == 0)
        ppit = PI_NOINTERSECTION;
      else
        ppit = PI_INTERSECTION;
    }

  friend ostream & operator<<(ostream & os,
                              const PolygonIntersector<UnitType> & pi)
    {
      switch(pi.ppit)
      {
        case PI_NOINTERSECTION:
          os << "No intersection: ";
          break;
        case PI_1CONTAINS2:
          os << "1 contains 2: ";
          break;
        case PI_2CONTAINS1:
          os << "2 contains 1: ";
          break;
        case PI_INTERSECTION:
          os << "intersection: ";
          break;
        default:
          os << "ERROR: ";
          break;
      }
      os << pi.pKeyVs[0].size();
      for(unsigned int i = 0; i < pi.pKeyVs[0].size(); i++)
        os << " " << pi.pKeyVs[0][i] << ":" << pi.pKeyVs[1][i];
      os << endl;
      return os;
    }

  vector<CVHPoint<UnitType> > pKeyVs[2]; // 0=clockwise, 1=counter-clockwise
  ConvexHull<UnitType> cvh;
  PolygonIntersectionType ppit;
  int pNumBridges;
};

#endif
