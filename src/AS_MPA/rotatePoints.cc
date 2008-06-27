
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
/* $Id: rotatePoints.cc,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $ */
#include <cstdio>  // for sscanf
#include <iostream>
#include <string>
#include "Point.h"


int main(int argc, char ** argv)
{
  if(argc != 2)
  {
    cerr << "Usage: " << argv[0] << " degrees-to-rotate < infile > outfile\n";
    exit(1);
  }

  Point<int> p;
  int x, y;
  char line1[4096];
  double degrees = strtod(argv[1], NULL);
  while(cin.getline(line1, 4095))
  {
    if(line1[0] == '\0')
      cout << endl;
    else
    {
      sscanf(line1, "%d %d", &x, &y);
      p.setX((int) x);
      p.setY((int) y);
      p.rotateByDegrees(degrees);
      cout << p.getX() << " " << p.getY() << endl;
    }
  }
  return 0;
}
