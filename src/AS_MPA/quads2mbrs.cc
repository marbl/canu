
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
/* $Id: quads2mbrs.cc,v 1.2 2004-09-23 20:25:24 mcschatz Exp $ */
#include <cstdio>
#include <iostream>
#include <string>

using namespace std;

#include "MPTypes.h"
#include "Point.h"
#include "Quadrilateral.h"


int main(int argc, char ** argv)
{
  Quadrilateral<int, int> q;
  double degrees = 0.;

  if(argc > 1)
    degrees = strtod(argv[1], NULL);
  
  char line[4096];
  while(cin.getline(line, 4095))
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

      if(degrees != 0.)
        q.rotateByDegrees(degrees);

      Interval<ID_TYPE, UNIT_TYPE> xProjection = q.getXProjection();
      Interval<ID_TYPE, UNIT_TYPE> yProjection = q.getYProjection();
      
      cout << xProjection.getMin() << "\t" << yProjection.getMax() << endl;
      cout << xProjection.getMax() << "\t" << yProjection.getMax() << endl;
      cout << xProjection.getMax() << "\t" << yProjection.getMin() << endl;
      cout << xProjection.getMin() << "\t" << yProjection.getMin() << endl;
      cout << xProjection.getMin() << "\t" << yProjection.getMax() << endl;
      cout << endl;
    }
  }
  return 0;
}
