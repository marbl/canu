
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
/* $Id: pullFragmentLocations.cc,v 1.2 2004-09-23 20:25:24 mcschatz Exp $ */
#include <cstdio>  // for sscanf
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <map>

using namespace std;


int main(int argc, char ** argv)
{
  if(argc < 3)
  {
    cerr << "Usage: " << argv[0] << " uidFile [intra/inter files]\n";
    exit(0);
  }
  
  ifstream fuids(argv[1], ios::in);
  if(!fuids.good())
    exit(-1);
  map<unsigned long, int> uids;
  char line[2048];
  while(fuids.getline(line, 2047))
  {
    unsigned long uid;
    sscanf(line, "%lu", &uid);
    uids[uid] = 0;
  }
  fuids.close();

  cerr << uids.size() << " uids in table\n";

  for(int i=2; i < argc; i++)
  {
    
    ifstream fin(argv[i], ios::in);
    if(!fin.good())
      exit(-2);

    while(fin.getline(line, 2047))
    {
      int chrom;
      unsigned long uid;
      int fivep, threep;
      sscanf(line, "%lu %d %d %d", &uid, &chrom, &fivep, &threep);
      if(uids.find(uid) != uids.end())
      {
        uids[uid] = uids[uid] + 1;
        if(fivep < threep)
        {
          cout << chrom << " " << fivep << " " << threep << " " << uid << endl;
        }
        else
        {
          cout << chrom << " " << threep << " " << fivep << " " << uid << endl;
        }
      }
    }
    fin.close();
  }

  map<unsigned long, int>::iterator iter;
  for(iter=uids.begin(); iter != uids.end(); iter++)
  {
    if(iter->second == 0)
      cerr << iter->first << endl;
  }
  return 0;
}
