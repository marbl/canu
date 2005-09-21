
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
/* $Id: processInterCG.cc,v 1.6 2005-09-21 20:13:07 catmandew Exp $ */
#include <cstdio>  // for sscanf
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <unistd.h> /* man 3 getopt */

using namespace std;

#include "Polygon.h"
#include "MPTypes.h"
#include "Interval.h"
#include "PolygonIntersector.h"
#include "TwoDIntervalClique.h"
#include "IntervalSetCoverSolver.h"
#include "MatePairPolygon.h"
#include "CompositeMPPolygon.h"
#include "CloneLibrary.h"
#include "clusterMPs.h"

#include <unistd.h> /* for getopt */

//#define USE_SATISFIEDS

//#define OUTPUT_ATAC_FORMAT

#define SEQ_OFFSET 500000000

#ifndef LINE_SIZE
#define LINE_SIZE 4096
#endif

char * MatePairLabel[MPI_NUM_INDICES] =
{
  "stretched",
  "compressed",
  "outtie",
  "normal",
  "antinormal",
  "inversion",
  "transposition",
  "satisfied",
  "interSequence",
  "unknown"
};


/*
  PURPOSE:
    Process inter-sequence mate pairs - mate pairs for which each fragment
    was mapped to a different sequence - to identify mis-assemblies and
    perhaps re-sequenced plates or collapsed repeats

  PLAN:
    Read in clone library length estimates
    
    Read in all intra-sequence mate pairs for a given sequence/assembly
      store each mate pair twice (by left and by right frag coord)
      sort by 'left' coorinate
      
    Read in all inter-sequence mate pairs for a given sequence/assembly
      Maintain separate list of mates for each other sequence
        optionally omit lower numbered sequences to avoid duplication
      Sort each list by coordinate on other sequence

    c = this sequence
    For each other sequence, o
      for each mate, m, on that sequence
        create a list, l, of agreeing mate pairs
        determine relative coordinate range & orientation of m
          on c based on m's mate's coordinate & orientation
        for each subsequent mate, n, on o until none for 100kb
          if it agrees with m, add it to l
      if l has 1 member, forget about it
        otherwise print interval on c and o and list of mate pairs & coords
 */

void Usage(char * progname, char * message)
{
  if(message != NULL)
    cerr << endl << message << endl;
  cerr << "Usage: " << progname << " -l lib  -e filename  [-n #]  [-f #]  [-a]  [-g]\n";
  cerr << "\t-l lib        name of clone library file\n";
  cerr << "\t-e filename   name of inter-sequence matepairs file\n";
  cerr << "\t                filename must have form:\n";
  cerr << "\t                assemblyName_#_inter.txt\n";
  cerr << "\t                where # is the scaffold or sequence number\n";
  cerr << "\t-n #          number of stddevs from mean that is excessive\n\n";
  cerr << "\t                default is " << STDDEVS_THRESHOLD << endl;
  cerr << "\t-f #          filter out mate pair sets with fewer members\n";
  cerr << "\t                default is " << CONFIRMATION_THRESHOLD << endl;
  cerr << "\t-a            generate ATA-fomratted output\n";
  cerr << "\t-g            generate gnuplot output\n";
  
  cerr << endl;
  exit(1);
}


void ReadInterSequenceMPs(vector<MatePair> & mps,
                            ifstream & fin,
                            ID_TYPE filterSeqID)
{
  char line[LINE_SIZE];
  while(fin.getline(line, LINE_SIZE-1))
  {
    ID_TYPE otherSeqID;
    sscanf(line, "%*s %*d %*d %*s %*s " F_MPID, &otherSeqID);
    if(otherSeqID == filterSeqID)
    {
      MatePair mp;
      mp.setFromInterSequenceString(line);
      mp.setRightCoord(mp.getRightCoord() + SEQ_OFFSET);
      mps.push_back(mp);
    }
  }
}


int main(int argc, char ** argv)
{
  char * libFilename = NULL;
  char * ewFilename = NULL;
  char assembly[4096];
  char sequence[4096];
  double numStddevs = STDDEVS_THRESHOLD;
  int filterThresh = CONFIRMATION_THRESHOLD;
  bool printGnuplot = false;
  bool printATA = false;

  {
    int ch, errflg = 0;
    while(!errflg && ((ch = getopt(argc, argv, "l:e:n:f:ag")) != EOF))
    {
      switch(ch)
      {
        case 'l':
          libFilename = optarg;
          break;
        case 'e':
          ewFilename = optarg;
          break;
        case 'n':
          numStddevs = atof(optarg);
          break;
        case 'f':
          filterThresh = atoi(optarg);
          break;
        case 'a':
          printATA = true;
          break;
        case 'g':
          printGnuplot = true;
          break;
        default:
          errflg++;
          break;
      }
    }
    if(libFilename == NULL)
      Usage(argv[0], "Please specify a clone library filename");
    if(ewFilename == NULL)
      Usage(argv[0], "Please specify an inter-sequence mate pair filename");
    if(numStddevs <= 0)
      Usage(argv[0], "Please specify a positive number of std deviations");
    {
      char * ptr;
      strcpy(assembly, ewFilename);
      ptr = index(assembly, (int) '_');
      assert(ptr != NULL);
      ptr[0] = '\0';
      strcpy(sequence, ++ptr);
      ptr = index(sequence, (int) '_');
      assert(ptr != NULL);
      ptr[0] = '\0';
    }
  }

  // read the clone library files
  map<ID_TYPE, CloneLibrary> libs;
  ifstream flib(libFilename, ios::in);
  if(!flib.good())
  {
    cerr << "Failed to open " << libFilename << " for reading\n";
    exit(-1);
  }
  // cerr << "Reading clone library data from file " << libFilename << endl;
  ReadCloneLibs(libs, flib);
  flib.close();

  // open & pre-read input file
  set<ID_TYPE> otherSequences;
  char line[LINE_SIZE];
  ifstream fe(ewFilename, ios::in);
  if(!fe.good())
  {
    cerr << "Failed to open " << ewFilename << " for reading\n";
    exit(-1);
  }
  while(fe.getline(line, LINE_SIZE-1))
  {
    ID_TYPE otherSeqID;
    sscanf(line, "%*s %*d %*d %*s %*s " F_MPID, &otherSeqID);
    otherSequences.insert(otherSeqID);
  }
  fe.close();

  // open output file
  char fname[1024];
  ofstream gnuOS;
  ofstream ataOS;
  ID_TYPE atacCount;
  sscanf(sequence, F_MPID, &atacCount);
  atacCount *= 10;
  if(printGnuplot)
  {
    sprintf(fname, "%s.%s.inter.gp", assembly, sequence);
    gnuOS.open(fname, ios::out);
    if(!gnuOS.good())
    {
      cerr << "Failed to open " << fname << " for writing\n";
      exit(-1);
    }
  }
  if(printATA)
  {
    sprintf(fname, "%s.%s.inter.ata", assembly, sequence);
    ataOS.open(fname, ios::out);
    if(!ataOS.good())
    {
      cerr << "Failed to open " << fname << " for writing\n";
      exit(-1);
    }
    // write ata header lines
    ataOS << "! format ata 1.0\n";
    ataOS << "# numStddevs=" << numStddevs << endl;
  }
  ofstream listOS;
  {
    char tempFN[1024];
    sprintf(tempFN, "%s.%s.inter.breakpoints.txt", assembly, sequence);
    listOS.open(tempFN, ios::out);
    if(!listOS.good())
    {
      cerr << "Failed to open " << tempFN << " for writing\n";
      exit(-1);
    }
  }
  
  // iterate over all other sequences in the set
  set<ID_TYPE>::iterator siter;
  for(siter = otherSequences.begin();
      siter != otherSequences.end();
      siter++)
  {
    ID_TYPE otherSeqID = *siter;
    if(otherSeqID >= STR_TO_UID(sequence, NULL, 10))
      break;

    if(printGnuplot)
      gnuOS << "# " << otherSeqID << endl;
    
    vector<MatePair> mps;
    ifstream fel(ewFilename, ios::in);
    ReadInterSequenceMPs(mps, fel, otherSeqID);
    fel.close();
    
    list<MatePair> mpl[MPI_NUM_INDICES];
    int badLibMatePairCount = 0;
    for(int i = 0; i < mps.size(); i++)
    {
      if(libs[mps[i].getLibUID()].getMean() >=
         (numStddevs + .1) * libs[mps[i].getLibUID()].getStddev())
      {
        MatePairPolygon<UNIT_TYPE> mpp(mps[i],
                                       libs[mps[i].getLibUID()],
                                       numStddevs);
        CompositeMPPolygon<UNIT_TYPE> cmpp(mpp);
        mpl[cmpp.getType()].push_back(mps[i]);
      }
      else
      {
        badLibMatePairCount++;
      }
    }
    // cerr << "Omitted " << badLibMatePairCount << " matepair(s) from bad clone libraries\n";

    // process stretched, outtie, normal, & antinormal
    vector<CompositeMPPolygon<UNIT_TYPE> > mpps[MPI_NUM_INDICES]; // input
    map<ID_TYPE, MatePairPolygon<UNIT_TYPE> > mppsMap[MPI_NUM_INDICES];
    vector<CompositeMPPolygon<UNIT_TYPE> > cmpps[MPI_NUM_INDICES];
    for(int mpii = 0; mpii < MPI_INVERSION; mpii++)
    {
      if(mpii == MPI_COMPRESSED) continue;
      
      if(mpl[mpii].size() > 1)
        mpl[mpii].sort();
      
      list<MatePair>::iterator mpli;
      for(mpli = mpl[mpii].begin(); mpli != mpl[mpii].end(); mpli++)
      {
        // add a composite mate pair polygon to the vector
        MatePairPolygon<UNIT_TYPE> mpp(*mpli,
                                       libs[mpli->getLibUID()],
                                       numStddevs);
        CompositeMPPolygon<UNIT_TYPE> cmpp(mpp);
        cmpp.rotateByDegrees(45);
        mpps[mpii].push_back(cmpp);
        for(int j = 0; j < mpp.size(); j++)
            mpp[j].setY(mpp[j].getY() - SEQ_OFFSET);
        mppsMap[mpii][mpli->getLeftFragUID()] = mpp;
      }
      
      if(mpps[mpii].size() > 0)
      {
        // do the clustering
        vector<CompositeMPPolygon<UNIT_TYPE> > pmpps; // problematic
        vector<CompositeMPPolygon<UNIT_TYPE> > fmpps; // filtered out
        
        ClusterMPPs(mpps[mpii], cmpps[mpii], pmpps, fmpps,
                    (MatePairIndex_e) MPI_NORMAL, filterThresh);
        
        // de-rotate & shift composite polygons
        for(int i = 0; i < mpps[mpii].size(); i++)
        {
          mpps[mpii][i].rotateByDegrees(-45);
          for(int j = 0; j < mpps[mpii][i].size(); j++)
            for(int k = 0; k < mpps[mpii][i][j].size(); k++)
              mpps[mpii][i][j][k].setY(mpps[mpii][i][j][k].getY() -
                                       SEQ_OFFSET);
        }
        
        for(int i = 0; i < cmpps[mpii].size(); i++)
        {
          cmpps[mpii][i].rotateByDegrees(-45);
          for(int j = 0; j < cmpps[mpii][i].size(); j++)
            for(int k = 0; k < cmpps[mpii][i][j].size(); k++)
              cmpps[mpii][i][j][k].setY(cmpps[mpii][i][j][k].getY() -
                                        SEQ_OFFSET);

          if(printGnuplot)
          {
            cmpps[mpii][i].printForGnuplot(gnuOS, CR_NATIVE);
            for(int j = 0; j < cmpps[mpii][i].getNumMPs(); j++)
              mppsMap[mpii][cmpps[mpii][i].getMP(j).getLeftFragUID()].printForGnuplot(gnuOS, CR_NATIVE);
          }
          {
            UNIT_TYPE leftThis, rightThis, leftOther, rightOther;
            UNIT_TYPE val =
              cmpps[mpii][i].getMP(0).getRightCoord() - SEQ_OFFSET;
            leftThis = cmpps[mpii][i].getMP(0).getLeftCoord();
            rightThis = leftThis;
            leftOther = rightOther = val;
            for(int j = 1; j < cmpps[mpii][i].getNumMPs(); j++)
            {
              val = cmpps[mpii][i].getMP(j).getRightCoord() - SEQ_OFFSET;
              leftThis =
                (leftThis > cmpps[mpii][i].getMP(j).getLeftCoord()) ?
                cmpps[mpii][i].getMP(j).getLeftCoord() : leftThis;
              rightThis =
                (rightThis < cmpps[mpii][i].getMP(j).getLeftCoord()) ?
                cmpps[mpii][i].getMP(j).getLeftCoord() : rightThis;
              leftOther = (leftOther > val) ? val : leftOther;
              rightOther = (rightOther < val) ? val : rightOther;
            }

            switch(mpii)
            {
              case MPI_STRETCHED:
                if(printATA)
                  ataOS << "M HL HL" << atacCount++ << " . "
                        << assembly << ":" << sequence << " "
                        << leftThis << " " << rightThis - leftThis << " 1 "
                        << assembly << ":" << otherSeqID << " "
                        << leftOther << " " << rightOther - leftOther << " 1 "
                        << "> /weight=" << cmpps[mpii][i].getNumMPs() << endl;
                listOS << "I\t";
                break;
              case MPI_NORMAL:
                if(printATA)
                  ataOS << "M Hl Hl" << atacCount++ << " . "
                        << assembly << ":" << sequence << " "
                        << leftThis << " " << rightThis - leftThis << " 1 "
                        << assembly << ":" << otherSeqID << " "
                        << leftOther << " " << rightOther - leftOther << " -1 "
                        << "> /weight=" << cmpps[mpii][i].getNumMPs() << endl;
                listOS << "N\t";
                break;
              case MPI_ANTINORMAL:
                if(printATA)
                  ataOS << "M hL hL" << atacCount++ << " . "
                        << assembly << ":" << sequence << " "
                        << leftThis << " " << rightThis - leftThis << " -1 "
                        << assembly << ":" << otherSeqID << " "
                        << leftOther << " " << rightOther - leftOther << " 1 "
                        << "> /weight=" << cmpps[mpii][i].getNumMPs() << endl;
                listOS << "A\t";
                break;
              case MPI_OUTTIE:
                if(printATA)
                  ataOS << "M hl hl" << atacCount++ << " . "
                        << assembly << ":" << sequence << " "
                        << leftThis << " " << rightThis - leftThis << " -1 "
                        << assembly << ":" << otherSeqID << " "
                        << leftOther << " " << rightOther - leftOther << " -1 "
                        << "> /weight=" << cmpps[mpii][i].getNumMPs() << endl;
                listOS << "O\t";
                break;
              default:
                cerr << "Unknown matepair type: " << mpii << endl;
                exit(-100);
                break;
            } // switch
            
            // print basic output
            listOS << leftThis << "\t"
                   << rightThis - leftThis << "\t"
                   << leftOther << "\t"
                   << rightOther - leftOther << "\t"
                   << otherSeqID << "\t"
                   << cmpps[mpii][i].getNumMPs() << endl;
            
          } // faux scope
        } // loop over composite polygons
      } // if there are matepairs
    } // loop over matepair types
  } // loop over other sequences

  if(printATA)
    ataOS.close();
  if(printGnuplot)
    gnuOS.close();
  listOS.close();
  
  return 0;
}
