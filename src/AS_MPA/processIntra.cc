
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
/* $Id: processIntra.cc,v 1.1.1.1 2004-04-14 13:52:05 catmandew Exp $ */
#include <cstdio>  // for sscanf
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <map>

using namespace std;

#include "MPTypes.h"
#include "Interval.h"
#include "Polygon.h"
#include "PolygonIntersector.h"
#include "TwoDIntervalClique.h"
#include "IntervalSetCoverSolver.h"
#include "MatePairPolygon.h"
#include "CompositeMPPolygon.h"
#include "CloneLibrary.h"
#include "clusterMPs.h"

//#define DEBUG_PROCESSMPS
//#define DEBUG_PROCESSMPS_BIGTIME

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
  "interChromosome",
  "unknown"
};

void PrintRawMatePairs(vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                       MatePairIndex_e mpii,
                       char * assembly, int chromosome)
{
  char outFilename[1024];
  sprintf(outFilename, "%s.%03d.%s.raw",
          assembly, chromosome, MatePairLabel[mpii]);
  ofstream myOS(outFilename, ios::out);
  
  vector<CompositeMPPolygon<UNIT_TYPE> >::iterator iter;
  for(iter = cmpps.begin(); iter != cmpps.end(); iter++)
  {
    const MatePair & mp = iter->getMP(0);
    switch(mp.getOrientation())
    {
      case PAIR_INNIE:
        myOS << "I ";
        break;
      case PAIR_OUTTIE:
        myOS << "O ";
        break;
      case PAIR_NORMAL:
        myOS << "N ";
        break;
      case PAIR_ANTINORMAL:
        myOS << "A ";
        break;
      default:
        myOS << "U ";
        break;
    }
    myOS << mp.getLeftFragUID()
         << " " << mp.getRightFragUID()
         << " " << mp.getLibUID()
         << " " << mp.getLeftCoord()
         << " " << mp.getRightCoord() << endl;
  }
  myOS.close();
}


void PrintRawSatisfiedMatePairs(vector<MatePair> & smpsv,
                                char * assembly, int chromosome)
{
  char outFilename[1024];
  sprintf(outFilename, "%s.%03d.%s.raw",
          assembly, chromosome, MatePairLabel[MPI_SATISFIED]);
  ofstream myOS(outFilename, ios::out);

  vector<MatePair>::iterator iter;
  for(iter = smpsv.begin(); iter != smpsv.end(); iter++)
  {
    if(iter->getLeftCoord() < iter->getRightCoord())
      myOS << "I " << iter->getLeftFragUID()
           << " " << iter->getRightFragUID()
           << " " << iter->getLibUID()
           << " " << iter->getLeftCoord()
           << " " << iter->getRightCoord() << endl;
  }
  myOS.close();
}


void PrintBasicOutput(vector<CompositeMPPolygon<UNIT_TYPE> > & printmpps,
                      MatePairIndex_e mpii,
                      vector<CompositeMPPolygon<UNIT_TYPE> > & mpps,
                      map<uint64, int> & mppsMap,
                      char * assembly,
                      int chromosome,
                      double numStddevs,
                      int * localUID,
                      bool printGnuplot)
{
  char outFilename[1024];
  sprintf(outFilename, "%s.%03d.%s.ata",
          assembly, chromosome, MatePairLabel[mpii]);
  ofstream myOS(outFilename, ios::out);
  myOS << "! format ata 1.0\n";
  myOS << "# numStddevs=" << numStddevs << endl;
  
  vector<CompositeMPPolygon<UNIT_TYPE> >::iterator mppIter;
  for(mppIter = printmpps.begin(); mppIter != printmpps.end(); mppIter++)
  {
    mppIter->printATA(myOS, assembly, chromosome, (*localUID)++, true);
  }
  myOS.close();
  
  if(printGnuplot)
  {
    sprintf(outFilename, "%s.%03d.%s.gp",
            assembly, chromosome, MatePairLabel[mpii]);
    myOS.open(outFilename, ios::out);
    for(mppIter = printmpps.begin(); mppIter != printmpps.end(); mppIter++)
    {
      mppIter->printForGnuplot(myOS);
      
      for(unsigned int i = 0; i < mppIter->getNumMPs(); i++)
      {
        mpps[mppsMap[(mppIter->getMP(i)).getLeftFragUID()]].printForGnuplot(myOS);
      }
    }
    myOS.close();
  }
}


void Usage(char * progname, char * message)
{
  if(message != NULL)
    cerr << endl << message << endl;
  cerr << "Usage: " << progname << " [-l lib] [-m mps] [-n #stddevs]\n";
  cerr << "\t-l lib        name of clone library file\n";
  cerr << "\t-m mps        name of mate pair file\n";
  cerr << "\t-n #stddevs   number of stddevs from mean that is excessive\n\n";
  cerr << "\t                default is " << STDDEVS_THRESHOLD << endl;
  cerr << "\t-a assembly   ata assembly name\n";
  cerr << "\t-c chrom      chromosome number (index into fasta file)\n";
  cerr << "\t-g            do not print gnuplot output\n";
  cerr << "\t-f #          filter out mate pair sets with fewer members\n";
  cerr << "\t                default is " << CONFIRMATION_THRESHOLD << endl;
  cerr << "\t-e            process 'elsewheres'";
  cerr << "\n\n";
  exit(1);
}
  

int main(int argc, char ** argv)
{
  char * libFilename = NULL;
  char * mpFilename = NULL;
  char * assembly = NULL;
  int chromosome = -1;
  double numStddevs = STDDEVS_THRESHOLD;
  unsigned int filterThresh = CONFIRMATION_THRESHOLD;
  bool printGnuplot = true;
  unsigned int i;
  int localUID = 1;
  MatePairIndex_e maxMPIndex = MPI_TRANSPOSITION;
  
  {
    int ch, errflg = 0;
    while(!errflg && ((ch = getopt(argc, argv, "l:m:n:a:c:gf:e")) != EOF))
    {
      switch(ch)
      {
        case 'l':
          libFilename = optarg;
          break;
        case 'm':
          mpFilename = optarg;
          break;
        case 'n':
          numStddevs = atof(optarg);
          break;
        case 'a':
          assembly = optarg;
          break;
        case 'c':
          chromosome = atoi(optarg);
          break;
        case 'g':
          printGnuplot = false;
          break;
        case 'f':
          filterThresh = atoi(optarg);
          break;
        case 'e':
          maxMPIndex = MPI_ANTINORMAL;
          break;
        default:
          errflg++;
          break;
      }
    }
    if(libFilename == NULL)
      Usage(argv[0], "Please specify a clone library filename");
    if(mpFilename == NULL)
      Usage(argv[0], "Please specify a mate pair filename");
    if(numStddevs <= 0)
      Usage(argv[0], "Please specify a positive number of std deviations");
    if(assembly == NULL)
      Usage(argv[0], "Please specify an assembly name");
    if(chromosome < 0)
      Usage(argv[0], "Please specify a chromosome number");
  }

  // read the clone library files
  map<ID_TYPE, CloneLibrary> libs;
  ifstream flib(libFilename, ios::in);
  if(!flib.good())
  {
    cerr << "Failed to open " << libFilename << " for reading\n";
    exit(-1);
  }
  cerr << "Reading clone library data from file " << libFilename << endl;
  ReadCloneLibs(libs, flib);
  flib.close();

  map<ID_TYPE, CloneLibrary>::iterator liter;
  for(liter = libs.begin(); liter != libs.end(); liter++)
  {
    CloneLibrary lib = (CloneLibrary) (*liter).second;
#ifdef DEBUG_PROCESSMPS_BIGTIME
    cerr << (*liter).second << endl;
#endif
    if(lib.getMean() < (numStddevs + .1) * lib.getStddev())
    {
      cerr << "Omitting mate pairs from library " << lib.getUID() << endl
           << "numStddevs * libStdev is too close to libMean\n"
           << "(" << numStddevs << " * "
           << lib.getStddev() << " ~>= " << lib.getMean() << ")" << endl;
    }
  }

  vector<MatePair> smpsv;
  vector<CompositeMPPolygon<UNIT_TYPE> > mpps[MPI_NUM_INDICES]; // input
  map<uint64, int> mppsMap[MPI_NUM_INDICES];
  {
    vector<MatePair> mps;
    ifstream fmp(mpFilename, ios::in);
    if(!fmp.good())
    {
      cerr << "Failed to open " << mpFilename << " for reading\n";
      exit(-1);
    }
    cerr << "Reading mate pair data in file " << mpFilename << endl;
    ReadMatePairs(mps, fmp);
    fmp.close();

    cerr << "Separating mate pairs by type\n";
    // separate mate pairs by type
    list<MatePair> mpl[MPI_NUM_INDICES];
    list<MatePair>::iterator mpli;
    int badLibMatePairs = 0;
    for(i = 0; i < mps.size(); i++)
    {
      if(libs[mps[i].getLibUID()].getMean() >=
         (numStddevs + .1) * libs[mps[i].getLibUID()].getStddev())
      {
        MatePairPolygon<UNIT_TYPE> mpp(mps[i],
                                       libs[mps[i].getLibUID()],
                                       numStddevs);
        CompositeMPPolygon<UNIT_TYPE> cmpp(mpp);
        mpl[cmpp.getType()].push_back(mps[i]);
        /*
        if(cmpp.getType() == MPI_SATISFIED)
        {
          MatePair revmp(mps[i]);
          revmp.setLeftFrag(mps[i].getRightFrag());
          revmp.setRightFrag(mps[i].getLeftFrag());
          mpl[MPI_SATISFIED].push_back(revmp);
        }
        */
      }
      else
      {
        badLibMatePairs++;
      }
    }
    cerr << "Omitted " << badLibMatePairs << " matepair(s) from bad clone libraries\n";
    
    cerr << "Sorting mate pairs left to right\n";
    // sort by left coordinate, remove coincident pairs, & populate vectors
    int numCoincident = 0;
    int numKept = 0;
    for(int mpii = 0; mpii < MPI_INVERSION; mpii++)
    {
      cerr << "Working on " << MatePairLabel[mpii] << endl;
      if(mpl[mpii].size() > 1)
        mpl[mpii].sort();

      for(mpli = mpl[mpii].begin(); mpli != mpl[mpii].end();)
      {
        if((mpii == MPI_NORMAL || mpii == MPI_ANTINORMAL) &&
           mpli->getLeftCoord() + COINCIDENT_THRESHOLD >=
           mpli->getRightCoord() &&
           mpli->getLeftCoord() - COINCIDENT_THRESHOLD <=
           mpli->getRightCoord())
        {
          numCoincident++;
          mpli = mpl[mpii].erase(mpli);
          continue;
        }
           
        // add a composite mate pair polygon to the vector
        MatePairPolygon<UNIT_TYPE> mpp(*mpli,
                                       libs[mpli->getLibUID()],
                                       numStddevs);
        CompositeMPPolygon<UNIT_TYPE> cmpp(mpp);

        if(mpii != MPI_COMPRESSED)
          cmpp.rotateByDegrees(45);
        mpps[mpii].push_back(cmpp);

        // add it to the map for later reference
        mppsMap[mpii][mpli->getLeftFragUID()] = mpps[mpii].size() - 1;

        // for detecting inversions
        if(maxMPIndex >= MPI_INVERSION &&
           (cmpp.isNormal() || cmpp.isAntinormal()))
        {
          mppsMap[MPI_INVERSION][mpli->getLeftFragUID()] =
            mpps[MPI_INVERSION].size();
          mpps[MPI_INVERSION].push_back(cmpp);
        }
          
        // for detecting transpositions
        if(maxMPIndex >= MPI_TRANSPOSITION &&
           (cmpp.isOuttie() || cmpp.isStretched() || cmpp.isCompressed()))
        {
          if(!cmpp.isCompressed())
            cmpp.rotateByDegrees(-45);
          mppsMap[MPI_TRANSPOSITION][mpli->getLeftFragUID()] =
            mpps[MPI_TRANSPOSITION].size();
          mpps[MPI_TRANSPOSITION].push_back(cmpp);
        }
        
        numCoincident--;
        numKept++;
        // remove this & all coincident mate pairs
        // NOTE: this won't necessarily remove all coincident matepairs...
        MatePair mp(*mpli);
        while(mpli != mpl[mpii].end() &&
              mpli->getLeftCoord() + COINCIDENT_THRESHOLD >=
              mp.getLeftCoord() &&
              mpli->getLeftCoord() - COINCIDENT_THRESHOLD <=
              mp.getLeftCoord() &&
              mpli->getRightCoord() + COINCIDENT_THRESHOLD >=
              mp.getRightCoord() &&
              mpli->getRightCoord() - COINCIDENT_THRESHOLD <=
              mp.getRightCoord())
        {
          mpli = mpl[mpii].erase(mpli);
          numCoincident++;
        }
      }
      PrintRawMatePairs(mpps[mpii], (MatePairIndex_e) mpii,
                        assembly, chromosome);
    }
    
    // filter out coincident satisfied matepairs
    // NOTE: this won't necessarily remove all coincident matepairs...
    mpl[MPI_SATISFIED].sort();
    for(mpli = mpl[MPI_SATISFIED].begin(); mpli != mpl[MPI_SATISFIED].end();)
    {
      list<MatePair>::iterator mpliKeep = mpli;
      mpli++;
      // delete all (directly) subsequent matepairs coincident with this one
      while(mpli != mpl[MPI_SATISFIED].end() &&
            mpli->getLeftCoord() + COINCIDENT_THRESHOLD >=
            mpliKeep->getLeftCoord() &&
            mpli->getLeftCoord() - COINCIDENT_THRESHOLD <=
            mpliKeep->getLeftCoord() &&
            mpli->getRightCoord() + COINCIDENT_THRESHOLD >=
            mpliKeep->getRightCoord() &&
            mpli->getRightCoord() - COINCIDENT_THRESHOLD <=
            mpliKeep->getRightCoord())
      {
        mpli = mpl[MPI_SATISFIED].erase(mpli);
        numCoincident++;
      }
    }
    cerr << "Deleted " << numCoincident << " coincident mate pairs\n";

    // duplicate each entry in satisfied for later filtering
    int numTotal = mpl[MPI_SATISFIED].size();
    int counter = 0;
    for(counter = 0, mpli = mpl[MPI_SATISFIED].begin();
        counter < numTotal;
        mpli++, counter++)
    {
      // append forward & reverse for 'sorting' purposes
      MatePair revmp(*mpli);
      revmp.setLeftFrag(mpli->getRightFrag());
      revmp.setRightFrag(mpli->getLeftFrag());

      mpl[MPI_SATISFIED].push_back(revmp);
    }
    mpl[MPI_SATISFIED].sort();

    // copy satisfied matepairs into vector form
    for(mpli = mpl[MPI_SATISFIED].begin();
        mpli != mpl[MPI_SATISFIED].end();
        mpli++)
    {
      smpsv.push_back(*mpli);
    }
    PrintRawSatisfiedMatePairs(smpsv, assembly, chromosome);
      
    cerr << numKept << " unsatisfied mate pairs to be processed\n";
  }

  for(int l = 0; i < maxMPIndex; i++)
  {
    cerr << mpps[l].size() << " raw " << MatePairLabel[l] << " on input\n";
  }
  if(maxMPIndex <= MPI_SATISFIED)
  {
    cerr << smpsv.size() / 2 << " raw "
         << MatePairLabel[MPI_SATISFIED] << " on input\n";
  }

#ifdef DEBUG_PROCESSMPS
  cerr << "Looking for unions of intersecting MBRs\n";
#endif
  
  cerr << "Processing mate pair polygons\n";
  
  vector<CompositeMPPolygon<UNIT_TYPE> > cmpps[MPI_NUM_INDICES]; // clustered
  list<Rectangle<int, UNIT_TYPE> > rects1;
  for(int mpii = 0; mpii < maxMPIndex; mpii++)
  {
    cerr << "  working on " << MatePairLabel[mpii] << "...\n";
    rects1.clear();
    if(mpps[mpii].size() > 0)
    {
#ifdef DEBUG_PROCESSMPS
      cerr << "Processing " << MatePairLabel[mpii] << endl;
#endif
  
      // do the clustering
      vector<CompositeMPPolygon<UNIT_TYPE> > pmpps; // problematic
      vector<CompositeMPPolygon<UNIT_TYPE> > fmpps; // filtered out
      vector<CompositeMPPolygon<UNIT_TYPE> > mmpps; // mis-labelled libs

      ClusterMPPs(mpps[mpii], cmpps[mpii], pmpps, fmpps,
                  (MatePairIndex_e) mpii, filterThresh);

      // unrotate
      if(mpii != MPI_COMPRESSED)
      {
        for(unsigned int i = 0; i < mpps[mpii].size(); i++)
          mpps[mpii][i].rotateByDegrees(-45);
        for(unsigned int i = 0; i < cmpps[mpii].size(); i++)
          cmpps[mpii][i].rotateByDegrees(-45);
        for(unsigned int i = 0; i < pmpps.size(); i++)
          pmpps[i].rotateByDegrees(-45);
        for(unsigned int i = 0; i < fmpps.size(); i++)
          fmpps[i].rotateByDegrees(-45);
      }

      // filter out compressed/stretched that are in the wrong library
      if(mpii == MPI_STRETCHED || mpii == MPI_COMPRESSED)
      {
        FilterMislabelledLibs(cmpps[mpii], mmpps, 10);
        if(mmpps.size() > 0)
          PrintBasicOutput(mmpps, (MatePairIndex_e) mpii,
                           mpps[mpii], mppsMap[mpii],
                           assembly, 100*(chromosome+1), numStddevs,
                           &localUID, printGnuplot);
      }
      
      // not all overlapping normal/antinormals are consistent with a inversion
      if(mpii == MPI_INVERSION)
        RefineInversions(cmpps[mpii], smpsv, filterThresh);

      if(mpii == MPI_STRETCHED)
        RefineStretched(cmpps[mpii], smpsv);
      
      if(cmpps[mpii].size() > 0)
      {
        PrintBasicOutput(cmpps[mpii], (MatePairIndex_e) mpii,
                         mpps[mpii], mppsMap[mpii],
                         assembly, chromosome, numStddevs,
                         &localUID, printGnuplot);
      }

      if(pmpps.size() > 0)
      {
#ifdef DEBUG_PROCESSMPS
        cerr << pmpps.size() << " problematic mate pair groups:\n";
        
        vector<CompositeMPPolygon<UNIT_TYPE> >::iterator mppIter;
        for(mppIter = pmpps.begin(); mppIter != pmpps.end(); mppIter++)
        {
          cerr << "# No intersection for this set of mate pairs:\n";
          cerr << "# Category: " << MatePairLabel[mpii] << endl;
          cerr << "# MBR intersection: " << *mppIter << endl << endl;
        }
#endif
        PrintBasicOutput(pmpps, (MatePairIndex_e) mpii,
                         mpps[mpii], mppsMap[mpii],
                         assembly, -(chromosome+1), numStddevs,
                         &localUID, printGnuplot);
      }
      cerr << cmpps[mpii].size() << " " << MatePairLabel[mpii]
           << " matepair sets 'confirmed'.\n";
      cerr << fmpps.size() << " " << MatePairLabel[mpii]
           << " matepairs sets with fewer than "
           << filterThresh << " members.\n";
    }
  }

  DetectTranspositions(mpps[MPI_COMPRESSED],
                       mpps[MPI_STRETCHED],
                       cmpps[MPI_OUTTIE],
                       cmpps[MPI_TRANSPOSITION],
                       libs,
                       numStddevs,
                       filterThresh);
  cerr << cmpps[MPI_TRANSPOSITION].size() << " "
       << MatePairLabel[MPI_TRANSPOSITION]
       << " matepair sets 'confirmed'.\n";
  if(cmpps[MPI_TRANSPOSITION].size() > 0)
  {
    PrintBasicOutput(cmpps[MPI_TRANSPOSITION], MPI_TRANSPOSITION,
                     mpps[MPI_TRANSPOSITION], mppsMap[MPI_TRANSPOSITION],
                     assembly, chromosome, numStddevs,
                     &localUID, printGnuplot);
  }

  return 0;
}
