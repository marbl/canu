
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
/* $Id: testConvexHull.cc,v 1.2 2004-09-23 20:25:24 mcschatz Exp $ */
#include <iostream>

using namespace std;

#include "MPTypes.h"
#include "Point.h"
#include "Polygon.h"
#include "ConvexHull.h"
#include "PolygonIntersector.h"

bool RunTest(Polygon<UNIT_TYPE> & p1, Polygon<UNIT_TYPE> & p2,
             const list<CVHPoint<UNIT_TYPE> > & answer,
             int testNum,
             bool notSuperImposed,
             char * message)
{
  cerr << "\nStarting test " << testNum << ": " << message << endl << endl;

  list<CVHPoint<UNIT_TYPE> > input;

  for(unsigned int i = 0; i < p1.size(); i++)
  {
    CVHPoint<UNIT_TYPE> p(p1[i], 0, i);
    input.push_back(p);
  }
  for(unsigned int i = 0; i < p2.size(); i++)
  {
    CVHPoint<UNIT_TYPE> p(p2[i], 1, i);
    input.push_back(p);
  }
  ConvexHull<UNIT_TYPE> cvh(input);
  ConvexHull<UNIT_TYPE> pcvh(p1, p2);
  
  cout << "  Polygon 0: " << p1 << endl;
  cout << "  Polygon 1: " << p2 << endl;
  cout << "  CVHullList: " << cvh << endl;
  cout << "  CVHullPoly: " << pcvh << endl;
  cout << "  Answer:     " << answer.size();
  list<CVHPoint<UNIT_TYPE> >::const_iterator aiter;
  for(aiter = answer.begin(); aiter != answer.end(); aiter++)
    cout << " " << *aiter;
  cout << endl;

  PolygonIntersector<UNIT_TYPE> pi;
  if(pi.polygonsIntersect(p1, p2))
  {
    Polygon<UNIT_TYPE> newp;
    cout << "Polygons intersect\n";
    pi.intersectPolygons(p1, p2, newp);
    cout << "Intersection:\n";
    cout << newp << endl;
  }
  else
    cout << "Polygons do not intersect\n";
  
  return true;
}

int main(int argc, char ** argv)
{
  Polygon<UNIT_TYPE> p1;
  Polygon<UNIT_TYPE> p2;
  int testNum = 0;
  list<CVHPoint<UNIT_TYPE> > answer;
  CVHPoint<UNIT_TYPE> cvhp;

  /*
    p1 is lower left, p2 is upper right
   */
  p1.reset();
  p1.append(5,5);
  p1.append(5,10);
  p1.append(10,10);
  p1.append(10,5);
  p2.reset();
  p2.append(15,15);
  p2.append(15,30);
  p2.append(30,30);
  p2.append(30,15);
  
  answer.clear();
  cvhp.set(5,5,0,0);
  answer.push_back(cvhp);
  cvhp.set(5,10,0,1);
  answer.push_back(cvhp);
  cvhp.set(15,30,1,1);
  answer.push_back(cvhp);
  cvhp.set(30,30,1,2);
  answer.push_back(cvhp);
  cvhp.set(30,15,1,3);
  answer.push_back(cvhp);
  cvhp.set(10,5,0,3);
  answer.push_back(cvhp);
  if(RunTest(p1, p2, answer, ++testNum, true,
             "similar, non-overlapping rectangles") == false)
    return 1;
  
  answer.clear();
  cvhp.set(5,5,0,0);
  answer.push_back(cvhp);
  cvhp.set(5,10,1,1);
  answer.push_back(cvhp);
  cvhp.set(10,10,1,2);
  answer.push_back(cvhp);
  cvhp.set(10,5,0,3);
  answer.push_back(cvhp);
  if(RunTest(p1, p1, answer, ++testNum, false,
             "super-imposed rectangles") == false)
    return 1;

  // make p2 contain p1
  answer.clear();
  cvhp.set(4,4,1,0);
  answer.push_back(cvhp);
  cvhp.set(4,11,1,1);
  answer.push_back(cvhp);
  cvhp.set(11,11,1,2);
  answer.push_back(cvhp);
  cvhp.set(11,4,1,3);
  answer.push_back(cvhp);

  p2.reset();
  p2.append(4,4);
  p2.append(4,11);
  p2.append(11,11);
  p2.append(11,4);
  if(RunTest(p1, p2, answer, ++testNum, true,
             "contained rectangles") == false)
    return 1;

  // make p2 intersect p1
  p2.reset();
  p2.append(7,7);
  p2.append(7,17);
  p2.append(17,17);
  p2.append(17,7);
  answer.clear();
  cvhp.set(5,5,0,0);
  answer.push_back(cvhp);
  cvhp.set(5,10,0,1);
  answer.push_back(cvhp);
  cvhp.set(7,17,1,1);
  answer.push_back(cvhp);
  cvhp.set(17,17,1,2);
  answer.push_back(cvhp);
  cvhp.set(17,7,1,3);
  answer.push_back(cvhp);
  cvhp.set(10,5,0,3);
  answer.push_back(cvhp);
  if(RunTest(p1, p2, answer, ++testNum, true,
             "intersecting rectangles") == false)
    return 1;

  // p1 has segment within p2 segment
  p1.reset();
  p1.append(1,7);
  p1.append(6,3);
  p1.append(4,2);
  p1.append(3,3);
  p2.reset();
  p2.append(2,4);
  p2.append(5,5);
  p2.append(5,1);
  answer.clear();
  cvhp.set(1,7,0,0);
  answer.push_back(cvhp);
  cvhp.set(5,5,1,1);
  answer.push_back(cvhp);
  cvhp.set(6,3,0,1);
  answer.push_back(cvhp);
  cvhp.set(5,1,1,2);
  answer.push_back(cvhp);
  cvhp.set(2,4,1,0);
  answer.push_back(cvhp);
  if(RunTest(p1, p2, answer, ++testNum, true,
             "left has contained segment") == false)
    return 1;

  // try larger polygons
  p1.reset();
  p1.append(2,-4);
  p1.append(3,-2);
  p1.append(5,-1);
  p1.append(7,-2);
  p1.append(8,-4);
  p1.append(7,-6);
  p1.append(5,-7);
  p1.append(3,-6);
  p2.reset();
  p2.append(-4,0);
  p2.append(-3,2);
  p2.append(-1,3);
  p2.append(2,2);
  p2.append(4,-1);
  p2.append(4,-4);
  p2.append(1,-5);
  p2.append(-3,-3);
  answer.clear();
  cvhp.set(-4,0,0,0);
  answer.push_back(cvhp);
  cvhp.set(-3,2,0,1);
  answer.push_back(cvhp);
  cvhp.set(-1,3,0,2);
  answer.push_back(cvhp);
  cvhp.set(2,2,0,3);
  answer.push_back(cvhp);
  cvhp.set(7,-2,1,3);
  answer.push_back(cvhp);
  cvhp.set(8,-4,1,4);
  answer.push_back(cvhp);
  cvhp.set(7,-6,1,5);
  answer.push_back(cvhp);
  cvhp.set(5,-7,1,6);
  answer.push_back(cvhp);
  cvhp.set(-3,-3,0,7);
  answer.push_back(cvhp);
  if(RunTest(p2, p1, answer, ++testNum, true,
             "intersecting polygons in all quadrants") == false)
    return 1;

  p1.reset();
  p1.append(2,3);
  p1.append(2,7);
  p1.append(7,2);
  p1.append(3,2);
  p2.reset();
  p2.append(4,3);
  p2.append(4,5);
  p2.append(8,1);
  p2.append(6,1);
  answer.clear();
  cvhp.set(2,3,0,0);
  answer.push_back(cvhp);
  cvhp.set(2,7,0,1);
  answer.push_back(cvhp);
  cvhp.set(8,1,1,2);
  answer.push_back(cvhp);
  cvhp.set(6,1,1,3);
  answer.push_back(cvhp);
  cvhp.set(3,2,0,3);
  answer.push_back(cvhp);
  if(RunTest(p1, p2, answer, ++testNum, true,
             "Mate pair simulation") == false)
    return 1;
  
  p1.reset();
  p1.append(2,2);
  p1.append(5,4);
  p1.append(7,4);
  p1.append(4,2);
  p2.reset();
  p2.append(2,4);
  p2.append(4,4);
  p2.append(5,2);
  p2.append(3,2);
  answer.clear();
  cvhp.set(2,2,0,0);
  answer.push_back(cvhp);
  cvhp.set(2,4,1,0);
  answer.push_back(cvhp);
  cvhp.set(7,4,0,2);
  answer.push_back(cvhp);
  cvhp.set(5,2,1,2);
  answer.push_back(cvhp);
  if(RunTest(p1, p2, answer, ++testNum, true,
             "Partial segment overlap") == false)
    return 1;
  
  p1.reset();
  p1.append(2,2);
  p1.append(2,7);
  p1.append(7,7);
  p1.append(7,2);
  p2.reset();
  p2.append(2,2);
  p2.append(2,7);
  p2.append(7,7);
  p2.append(8,4);
  p2.append(7,2);
  answer.clear();
  cvhp.set(2,2,0,0);
  answer.push_back(cvhp);
  cvhp.set(2,7,1,1);
  answer.push_back(cvhp);
  cvhp.set(7,7,1,2);
  answer.push_back(cvhp);
  cvhp.set(8,4,1,3);
  answer.push_back(cvhp);
  cvhp.set(7,2,0,3);
  answer.push_back(cvhp);
  if(RunTest(p1, p2, answer, ++testNum, true,
             "Almost super-imposed") == false)
    return 1;
  
  return 0;
}
