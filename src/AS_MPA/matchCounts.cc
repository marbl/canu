
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
/* $Id: matchCounts.cc,v 1.1.1.1 2004-04-14 13:52:03 catmandew Exp $ */
#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;

#define HUMAN

#ifdef HUMAN
#define NUM_MATED_FRAGS 21525272
#else
#define NUM_MATED_FRAGS  2309852
#endif

int main(int argc, char ** argv)
{
  char line[1024];
  int counts[NUM_MATED_FRAGS];
  int id;

  memset(counts, 0, NUM_MATED_FRAGS * sizeof(int));
  
  while(cin.getline(line, 1023))
    counts[atoi(line)]++;

  cout << "Matches per fragment\n";
  for(int i = 0; i < NUM_MATED_FRAGS; i++)
    cout << counts[i] << endl;
  
  return 0;
}
