
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
/* $Id: findTTIs.cc,v 1.4 2005-03-22 19:48:58 jason_miller Exp $ */
#include <cstdio>  // for sscanf
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

using namespace std;

#include "Interval.h"
#include "Rectangle.h"
#include "TwoDIntervalClique.h"
#include "IntervalSetCoverSolver.h"

//#define DEBUG_FMPS
//#define FEEDBACK_FMPS

#ifdef NEVER
//#define USE_INT_IDS
#ifdef USE_INT_IDS
#define ID_TYPE          int
#define ID_SCAN_FORMAT   "%d"
#else
#define ID_TYPE          cds_uint64
#define ID_SCAN_FORMAT   F_U64
#endif

#define USE_INT_UNITS
#ifdef USE_INT_UNITS
#define UNIT_TYPE        int
#define UNIT_SCAN_FORMAT "%d"
#else
#define UNIT_TYPE        double
#define UNIT_SCAN_FORMAT "%lf"
#endif
#endif

void ReadQuadrilaterals(vector<Quadrilateral<ID_TYPE, UNIT_TYPE> > & quads,
                        istream & fin)
{
  char line[4096];

  while(fin.getline(line, 4095))
  {
    ID_TYPE id;
    UNIT_TYPE x[4], y[4];
    if(sscanf(line,
              UNIT_SCAN_FORMAT " "
              UNIT_SCAN_FORMAT " "
              UNIT_SCAN_FORMAT " "
              UNIT_SCAN_FORMAT " "
              UNIT_SCAN_FORMAT " "
              UNIT_SCAN_FORMAT " "
              UNIT_SCAN_FORMAT " "
              UNIT_SCAN_FORMAT " "
              ID_SCAN_FORMAT,
              &(x[0]), &(y[0]),
              &(x[1]), &(y[1]),
              &(x[2]), &(y[2]),
              &(x[3]), &(y[3]),
              &id) == 9)
    {
      Quadrilateral<ID_TYPE, UNIT_TYPE> q;
      q.reset(x[0], y[0]);
      q.setID(id);
      for(int i = 1; i < 4; i++)
        q.setPoint(x[i], y[i], i);
      quads.push_back(q);
    }
    else
    {
      cerr << "Line doesn't match format:\n"
           << "x[0] y[0] x[1] y[1] x[2] y[2] x[3] y[3] id\n"
           << line << endl;
    }
  }
}


int main(int argc, char ** argv)
{
  vector<Quadrilateral<ID_TYPE, UNIT_TYPE> > quads;
  vector<Quadrilateral<ID_TYPE, UNIT_TYPE> >::iterator qiter;
  double degrees = 0.0;
  list<TwoDIntervalClique<ID_TYPE, UNIT_TYPE> > twoDIs;
  list<Rectangle<ID_TYPE, UNIT_TYPE> >::iterator iter1;
  list<Rectangle<ID_TYPE, UNIT_TYPE> >::iterator iter2;
  list<TwoDIntervalClique<ID_TYPE, UNIT_TYPE> >::iterator titer1;
  int i;

  ReadQuadrilaterals(quads, cin);
  cerr << "Read " << quads.size() << " quadrilaterals\n";
  if(argc > 1)
  {
    degrees = strtod(argv[1], NULL);
    cerr << "Rotating " << degrees << " degrees before clustering\n";

    if(degrees != 0.0)
    {
      for(qiter = quads.begin(); qiter != quads.end(); qiter++)
        qiter->rotateByDegrees(degrees);
    }
  }
  
  map<ID_TYPE, int> id2Index;
  list<Rectangle<ID_TYPE, UNIT_TYPE> > rects1;
  vector<Rectangle<ID_TYPE, UNIT_TYPE> > rects2;
  for(i = 0, qiter = quads.begin(); qiter != quads.end(); qiter++, i++)
  {
    Rectangle<ID_TYPE, UNIT_TYPE> rect(qiter->getXProjection(),
                                       qiter->getYProjection());
    rects1.push_back(rect);
    rects2.push_back(rect);
    id2Index[rect.getID()] = i;
    
#ifdef DEBUG_FMPS
    cerr << rects2[i] << endl;
#endif
    
  }

#ifdef FEEDBACK_FMPS
  cerr << "Grouping rectangles.\n";
#endif

  rects1.sort();
  // for each un-grouped rectangle
  for(iter1 = rects1.begin(); iter1 != rects1.end();)
  {
    // seed a new twoDI & delete the rectangle
    TwoDIntervalClique<ID_TYPE, UNIT_TYPE> tdi(*iter1);
    iter1 = rects1.erase(iter1);

    // compare with all other remaining rectangles
    iter2 = iter1;
    for(++iter2; iter2 != rects1.end();)
    {
      // if intersection, add to the new tdi & delete the rectangle
      if(tdi.intersects(*iter2))
      {
        tdi.addRectangle(*iter2);
        iter2 = rects1.erase(iter2);
      }
      else
        iter2++;
    }

    // compare with all other tdis
    for(titer1 = twoDIs.begin(); titer1 != twoDIs.end();)
    {
      // if intersection, add to the new tdi & delete the old tdi
      if(tdi.intersects(*titer1))
      {
        tdi.add(*titer1);
        titer1 = twoDIs.erase(titer1);
      }
      else
        titer1++;
    }

    // add the new tdi to the set of tdis
    twoDIs.push_back(tdi);
  }

#ifdef FEEDBACK_FMPS
  cerr << "Reduced to " << twoDIs.size() << " smaller problems\n";
#endif

#ifdef DEBUG_FMPS
  for(titer1 = twoDIs.begin(); titer1 != twoDIs.end(); titer1++)
  {
      cerr << *titer1 << endl;
  }
#endif
  
  // now, within each twoDI, find maximal cliques
  vector<TwoDIntervalClique<ID_TYPE, UNIT_TYPE> > ctdis;
  int ti;
  for(ti = 0, titer1 = twoDIs.begin();
      titer1 != twoDIs.end();
      titer1++, ti++)
  {
    if(titer1->getNumIDs() == 1)
    {
      cout << *titer1 << endl;
    }
    else
    {
      vector<Rectangle<ID_TYPE, UNIT_TYPE> > rects3;
      list<ID_TYPE> ids = titer1->getIDs();
      list<ID_TYPE>::iterator liter;
      
      for(liter = ids.begin(); liter != ids.end(); liter++)
        rects3.push_back(rects2[id2Index[*liter]]);

#ifdef DEBUG_FMPS
      cerr << "Rectangle group " << ti << endl;
      vector<Rectangle<ID_TYPE, UNIT_TYPE> >::iterator r3iter;
      for(r3iter = rects3.begin(); r3iter != rects3.end(); r3iter++)
        cerr << *r3iter << endl;
#endif
      
      IntervalSetCoverSolver<ID_TYPE, UNIT_TYPE> iscs(rects3);
      ctdis.clear();
      iscs.solve(ctdis);

      vector<TwoDIntervalClique<ID_TYPE, UNIT_TYPE> >::iterator titer2;
      for(titer2 = ctdis.begin(); titer2 != ctdis.end(); titer2++)
      {
        if(degrees != 0.0)
          titer2->rotateByDegrees(-degrees);
        cout << *titer2 << endl;
      }
    }
  }
  return 0;
}
