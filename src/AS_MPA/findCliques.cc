
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
/* $Id: findCliques.cc,v 1.4 2005-03-22 19:48:58 jason_miller Exp $ */

#include <cstdio>  // for sscanf
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

using namespace std;

#include "Quadrilateral.h"
#include "cds.h"

#include "Interval.h"
#include "Rectangle.h"
#include "IntervalClique.h"
#include "IntervalCliqueFinder.h"
#include "TwoDIntervalClique.h"
#include "IntervalSetCoverSolver.h"

//#define DEBUG_FC
#define FEEDBACK_FC

/*
  Read in quadrilaterals
  
  Find cliques in x dimension
  Find cliques in y dimension

  Find best x&y clique intersection
 */

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

  ReadQuadrilaterals(quads, cin);
  cerr << "Read " << quads.size() << " quadrilaterals\n";
  if(argc > 1)
  {
    degrees = strtod(argv[1], NULL);
    cerr << "Rotating " << degrees << " degrees before clustering\n";
  }

  if(degrees != 0.0)
  {
    for(qiter = quads.begin(); qiter != quads.end(); qiter++)
      qiter->rotateByDegrees(degrees);
  }

#if defined(FEEDBACK_FC) || defined(DEBUG_FC)
  cerr << "Performing 2D greedy set cover algorithm\n";
#endif
  
  IntervalSetCoverSolver<ID_TYPE, UNIT_TYPE> iscs(quads);
  vector<TwoDIntervalClique<ID_TYPE, UNIT_TYPE> > twoDIs;
  iscs.solve(twoDIs);

#if defined(FEEDBACK_FC) || defined(DEBUG_FC)
  cerr << "Printing 2D cliques\n";
#endif

  vector<TwoDIntervalClique<ID_TYPE, UNIT_TYPE> >::iterator twoDIter;
  for(twoDIter = twoDIs.begin(); twoDIter != twoDIs.end(); twoDIter++)
  {
    if(degrees != 0.0)
      twoDIter->rotateByDegrees(-degrees);
    cout << *twoDIter;
  }
  
  return 0;
}
