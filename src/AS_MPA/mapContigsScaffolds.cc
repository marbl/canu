
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
/* $Id: mapContigsScaffolds.cc,v 1.1.1.1 2004-04-14 13:52:03 catmandew Exp $ */
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char ** argv)
{
  if(argc != 3)
  {
    cerr << "Usage: " << argv[0] << "  allGaps.txt   all.out\n";
    cerr << "contigs written to stderr, scaffolds to stdout\n";
    exit(-1);
    
  }

  ifstream fi1(argv[1], ios::in);
  if(!fi1.good())
  {
    cerr << "File " << argv[1] << " ain't workin'\n";
    exit(-2);
  }
  ifstream fi2(argv[2], ios::in);
  if(!fi2.good())
  {
    cerr << "File " << argv[2] << " ain't workin'\n";
    exit(-2);
  }

  // read a line in from all.out
  char line2[2048];
  unsigned long sc, sl, sr, slength;
  fi2.getline(line2, 2047);
  sscanf(line2, "%*s %lu %lu %lu", &sc, &sl, &sr);
  slength = sr-sl;

  int scaffoldID = 0;
  int contigID = 0;
  char line1[2048];

  // last gap
  unsigned long lgc, lgl, lgr, lglength;
  fi1.getline(line1, 2047);
  sscanf(line1, "%lu %lu %lu", &lgc, &lgl, &lglength);
  lgr = lgl + lglength;

  unsigned long cc = lgc;
  unsigned long cl = 0;
  unsigned long cr = lgl;
  unsigned long clength = cr - cl;
  
  cout << sc << " " << scaffoldID << " " << sl << " " << slength << endl;
  if(cl > 0)
  {
    cerr << cc << " " << scaffoldID << " " << contigID << " "
         << cl << " " << clength << endl;
  }
  contigID = -1;

  // more gaps than scaffolds
  while(fi1.getline(line1, 2047))
  {
    // get current gap
    bool printIt = false;
    unsigned long gc, gl, gr, glength;
    sscanf(line1, "%lu %lu %lu", &gc, &gl, &glength);
    gr = gl + glength;

    // compute contig
    if(lglength < 20 && gc==lgc)
    {
      cr = gl;
    }
    else
    {
      cc = gc;
      cl = lgr;
      cr = gl;
      printIt = true;
    }
    clength = cr - cl;

    if(clength < 150)
      printIt = false;
    
    if(gc > sc || (cc == sc && cl > sr))
    {
      fi2.getline(line2, 2047);
      int nsc;
      sscanf(line2, "%*s %lu %lu %lu", &nsc, &sl, &sr);
      slength = sr-sl;
      scaffoldID++;
      cout << nsc << " " << scaffoldID << " " << sl << " " << slength << endl;
      if(nsc > sc)
      {
        cl = 0;
        clength = cr - cl;
      }
      sc = nsc;
    }

    if(sl > cl || sr < cr)
    {
      cerr << "Bad situation with gap & scaffold:\n";
      cerr << "Gap: " << gc << " " << gl << " " << gr << endl;
      cerr << "Scf: " << sc << " " << sl << " " << sr << endl;
      cerr << "Ctg: " << cc << " " << cl << " " << cr << endl;
      exit(1);
    }

    if(printIt)
    {
      contigID++;
      cerr << cc << " " << scaffoldID << " " << contigID << " "
           << cl << " " << clength << endl;
    }
    
    lgc = gc;
    lgl = gl;
    lgr = gr;
    lglength = glength;
  }
  fi1.close();
  fi2.close();
  
  return 0;
}
