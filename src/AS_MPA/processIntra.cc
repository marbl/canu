
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
/* $Id: processIntra.cc,v 1.11 2008-03-18 07:02:46 brianwalenz Exp $ */
#include <cstdio>  // for sscanf
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <map>
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
  "interSequence",
  "unknown"
};

void PrintRawMatePairs(vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                       MatePairIndex_e mpii,
                       char * assembly, char * seqID)
{
  char outFilename[1024];
  sprintf(outFilename, "%s.%s.%s.raw",
          assembly, seqID, MatePairLabel[mpii]);
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
                                char * assembly, char * seqID)
{
  char outFilename[1024];
  sprintf(outFilename, "%s.%s.%s.raw",
          assembly, seqID, MatePairLabel[MPI_SATISFIED]);
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

void PrintSingleBPFieldLabels(ofstream & listOS)
{
  listOS << "Breakpoint Interval\t\t"
         << "Deleted Length\t\t"
         << "Weight\n";
  listOS << "Left\t"
         << "Length\t"
         << "Min\t"
         << "Max\t"
         << "Mates\n";
}

void PrintDoubleBPFieldLabels(ofstream & listOS)
{
  listOS << "Breakpoint Intervals\n";
  
  listOS << "Left Interval\t\t"
         << "Right Interval\t\t"
         << "Problem Length\t"
         << "Weight\n";

  listOS << "Left\t"
         << "Length\t"
         << "Left\t"
         << "Length\t"
         << "Estimate\t"
         << "Mates\n";
}

void PrintOutput(vector<CompositeMPPolygon<UNIT_TYPE> > & printmpps,
                 MatePairIndex_e mpii,
                 vector<CompositeMPPolygon<UNIT_TYPE> > & mpps,
                 map<uint64_t, int> & mppsMap,
                 char * assembly,
                 char * seqID,
                 char * status,
                 double numStddevs,
                 int * relativeID,
                 ofstream & listOS,
                 bool printATA,
                 bool printGnuplot)
{
  char label[1024];
  char type[1024];
  bool polymorphic;
  char outFilename[1024];
  ofstream ataOS, gnuOS;

  polymorphic = (strcmp(status, "polymorphic") == 0);
  switch(mpii)
  {
    case MPI_STRETCHED:
      sprintf(label, "insertion");
      sprintf(type, "ins");
      break;
    case MPI_COMPRESSED:
      sprintf(label, "deletion");
      sprintf(type, "del");
      break;
    case MPI_INVERSION:
      sprintf(label, MatePairLabel[mpii]);
      sprintf(type, "inv");
      break;
    case MPI_TRANSPOSITION:
      sprintf(label, MatePairLabel[mpii]);
      sprintf(type, "trn");
      break;
    default:
      sprintf(label, MatePairLabel[mpii]);
      sprintf(type, "unk");
      break;
  }
  
  if(printATA)
  {
    sprintf(outFilename, "%s.%s.%s.%s.ata",
            assembly, seqID, status, label);
    ataOS.open(outFilename, ios::out);
    ataOS << "! format ata 1.0\n";
    ataOS << "# numStddevs=" << numStddevs << endl;
  }
    
  if(printGnuplot)
  {
    sprintf(outFilename, "%s.%s.%s.%s.gp",
            assembly, seqID, status, label);
    gnuOS.open(outFilename, ios::out);
  }

  vector<CompositeMPPolygon<UNIT_TYPE> >::iterator mppIter;
  for(mppIter = printmpps.begin(); mppIter != printmpps.end(); mppIter++)
  {
    if(printATA)
      mppIter->printATA(ataOS, assembly, seqID, (*relativeID)++, true);

    if(printGnuplot)
    {
      CompressedRepresentation_e cr = (mpii == MPI_COMPRESSED ?
                                       CR_NATIVE : CR_COMPATIBLE);
      
      mppIter->printForGnuplot(gnuOS, cr);
      
      for(unsigned int i = 0; i < mppIter->getNumMPs(); i++)
      {
        mpps[mppsMap[(mppIter->getMP(i)).getLeftFragUID()]].printForGnuplot(gnuOS, cr);
      }
    }

    mppIter->printSummary(listOS, seqID, type, polymorphic);
  }
  
  if(printGnuplot)
    gnuOS.close();
  
  if(printATA)
    ataOS.close();
}


void Usage(char * progname, char * message)
{
  if(message != NULL)
    cerr << endl << message << endl;
  cerr << "Usage: " << progname << " [-l lib] [-m mps] [-n #stddevs]\n";
  cerr << "\t-l lib        name of clone library file\n";
  cerr << "\t-m mps        name of mate pair file\n";
  cerr << "\t                filename must have form:\n";
  cerr << "\t                assemblyName_#_intra.txt\n";
  cerr << "\t                where # is the sequence number, such as a\n";
  cerr << "\t                scaffold or chromosome number\n";
  cerr << "\t-n #stddevs   number of stddevs from mean that is excessive\n";
  cerr << "\t                default is " << STDDEVS_THRESHOLD << endl;
  cerr << "\t-a            generate ATA-formatted output\n";
  cerr << "\t-g            generate gnuplot output\n";
  cerr << "\t-r            dump raw mate pairs to files\n";
  cerr << "\t-f #          filter out mate pair sets with fewer members\n";
  cerr << "\t                default is " << CONFIRMATION_THRESHOLD << endl;
  cerr << "\n\n";
  exit(1);
}
  

int main(int argc, char ** argv)
{
  char * libFilename = NULL;
  char * mpFilename = NULL;
  char assembly[4096];
  char seqID[4096];
  double numStddevs = STDDEVS_THRESHOLD;
  unsigned int filterThresh = CONFIRMATION_THRESHOLD;
  bool printATA = false;
  bool printGnuplot = false;
  bool printRaw = false;
  unsigned int i;
  int relativeID = 1;
  
  {
    int ch, errflg = 0;
    // while(!errflg && ((ch = getopt(argc, argv, "l:m:n:a:c:gf:e")) != EOF))
    while(!errflg && ((ch = getopt(argc, argv, "l:m:n:agrf:e")) != EOF))
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
          printATA = true;
          break;
        case 'g':
          printGnuplot = true;
          break;
        case 'r':
          printRaw = true;
          break;
        case 'f':
          filterThresh = atoi(optarg);
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

    {
      char * ptr;
      strcpy(assembly, mpFilename);
      ptr = index(assembly, (int) '_');
      assert(ptr != NULL);
      ptr[0] = '\0';
      strcpy(seqID, ++ptr);
      ptr = index(seqID, (int) '_');
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
    exit(1);
  }
  // cerr << "Reading clone library data from file " << libFilename << endl;
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
  map<uint64_t, int> mppsMap[MPI_NUM_INDICES];
  {
    vector<MatePair> mps;
    ifstream fmp(mpFilename, ios::in);
    if(!fmp.good())
    {
      cerr << "Failed to open " << mpFilename << " for reading\n";
      exit(1);
    }
    // cerr << "Reading mate pair data in file " << mpFilename << endl;
    ReadMatePairs(mps, fmp);
    fmp.close();

    // separate mate pairs by type
    // and count number in each library
    // and identify right-most coordinate
    list<MatePair> mpl[MPI_NUM_INDICES];
    list<MatePair>::iterator mpli;
    int badLibMatePairs = 0;
    uint64_t rightMostCoord = 0;
    for(i = 0; i < mps.size(); i++)
    {
      libs[mps[i].getLibUID()].incrementCount();
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
      rightMostCoord = (rightMostCoord > mps[i].getRightCoord() ?
                        rightMostCoord : mps[i].getRightCoord());
    }
    
    cout << rightMostCoord
         << " is right-most coordinate of any mated fragment\n";
    for(liter = libs.begin(); liter != libs.end(); liter++)
    {
      CloneLibrary lib = (CloneLibrary) (*liter).second;
      cout << lib.getCount()
           << " clones from library "
           << lib.getUID() << "\n";
    }
    cout << badLibMatePairs
         << " mate pairs omitted from bad clone libraries\n";
    
    // cerr << "Sorting mate pairs left to right\n";
    // sort by left coordinate, remove coincident pairs, & populate vectors
    int numCoincident = 0;
    int numKept = 0;
    for(int mpii = 0; mpii < MPI_INVERSION; mpii++)
    {
      // cerr << "Working on " << MatePairLabel[mpii] << endl;
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

        // rotate polygon cw 45deg so it's more orthogonally rectangular
        if(mpii != MPI_COMPRESSED)
          cmpp.rotateByDegrees(45);
        mpps[mpii].push_back(cmpp);

        // add it to the map for later reference
        mppsMap[mpii][mpli->getLeftFragUID()] = mpps[mpii].size() - 1;

        // for detecting inversions
        if(cmpp.isNormal() || cmpp.isAntinormal())
        {
          mppsMap[MPI_INVERSION][mpli->getLeftFragUID()] =
            mpps[MPI_INVERSION].size();
          mpps[MPI_INVERSION].push_back(cmpp);
        }
          
        // for detecting transpositions
        if(cmpp.isOuttie() || cmpp.isStretched() || cmpp.isCompressed())
        {
          // matepairs are not rotated relative to genomic axis
          if(!cmpp.isCompressed())
            cmpp.rotateByDegrees(-45);
          mppsMap[MPI_TRANSPOSITION][mpli->getLeftFragUID()] =
            mpps[MPI_TRANSPOSITION].size();
          mpps[MPI_TRANSPOSITION].push_back(cmpp);
        }
        
        numCoincident--;
        numKept++;
        // remove this & all coincident mate pairs
        // NOTE: this won't necessarily remove all coincident mate pairs...
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
      if(printRaw)
      {
        PrintRawMatePairs(mpps[mpii], (MatePairIndex_e) mpii,
                          assembly, seqID);
      }
    }
    
    // filter out coincident satisfied mate pairs
    // NOTE: this won't necessarily remove all coincident mate pairs...
    mpl[MPI_SATISFIED].sort();
    for(mpli = mpl[MPI_SATISFIED].begin(); mpli != mpl[MPI_SATISFIED].end();)
    {
      list<MatePair>::iterator mpliKeep = mpli;
      mpli++;
      // delete all (directly) subsequent mate pairs coincident with this one
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
    cout << numCoincident
         << " coincident mate pairs deleted\n";

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

    // copy satisfied mate pairs into vector form
    for(mpli = mpl[MPI_SATISFIED].begin();
        mpli != mpl[MPI_SATISFIED].end();
        mpli++)
    {
      smpsv.push_back(*mpli);
    }
    if(printRaw)
      PrintRawSatisfiedMatePairs(smpsv, assembly, seqID);
      
    cout << smpsv.size() / 2 << " raw "
         << MatePairLabel[MPI_SATISFIED] << " mate pairs\n";

    cout << numKept << " unsatisfied mate pairs to be processed\n";
    if(numKept == 0)
      return 0;
  }

  for(int l = 0; l < MPI_INVERSION; l++)
  {
    cout << mpps[l].size() << " raw " << MatePairLabel[l] << " mate pairs\n";
  }
#ifdef DEBUG_PROCESSMPS
  cerr << "Looking for unions of intersecting MBRs\n";
#endif
  
  // cerr << "Processing mate pair polygons\n";

  ofstream listOS;
  {
    char tempFN[1024];
    sprintf(tempFN, "%s.%s.intra.breakpoints.txt", assembly, seqID);
    listOS.open(tempFN, ios::out);
    if(!listOS.good())
    {
      cerr << "Failed to open " << tempFN << " for writing\n";
      exit(1);
    }
  }
  listOS << "SeqID\t"
         << "LeftBpILeftC\t"
         << "Length\t"
         << "RightBpILeftC/MinDel\t"
         << "Length/MaxDel\t"
         << "Type\t"
         << "PolyMorphic\t"
         << "NumMPs\n";
  
  vector<CompositeMPPolygon<UNIT_TYPE> > cmpps[MPI_NUM_INDICES]; // clustered
  list<Rectangle<int, UNIT_TYPE> > rects1;
  for(int mpii = 0; mpii <= MPI_INVERSION; mpii++)
  {
    vector<CompositeMPPolygon<UNIT_TYPE> > probs; // problematic
    vector<CompositeMPPolygon<UNIT_TYPE> > fews; // below-threshold
    vector<CompositeMPPolygon<UNIT_TYPE> > polys; // polymorphic
    vector<CompositeMPPolygon<UNIT_TYPE> > misls; // mis-labeled libs

    // cerr << "  working on " << MatePairLabel[mpii] << "...\n";
    rects1.clear();
    if(mpps[mpii].size() > 0)
    {
#ifdef DEBUG_PROCESSMPS
      cerr << "Processing " << MatePairLabel[mpii] << endl;
#endif

      // non-compressed mate pairs are all rotated cw 45deg
      // relative to genomic axis
      ClusterMPPs(mpps[mpii], cmpps[mpii], probs, fews,
                  (MatePairIndex_e) mpii, filterThresh);

      // rotate back to genomic axis
      if(mpii != MPI_COMPRESSED)
      {
        for(unsigned int i = 0; i < mpps[mpii].size(); i++)
          mpps[mpii][i].rotateByDegrees(-45);
        for(unsigned int i = 0; i < cmpps[mpii].size(); i++)
          cmpps[mpii][i].rotateByDegrees(-45);
        for(unsigned int i = 0; i < probs.size(); i++)
          probs[i].rotateByDegrees(-45);
        for(unsigned int i = 0; i < fews.size(); i++)
          fews[i].rotateByDegrees(-45);
      }

      // filter out compressed/stretched that are in the wrong library
      if(mpii == MPI_STRETCHED || mpii == MPI_COMPRESSED)
      {
        // FilterMislabelledLibs is rotation-independent
        FilterMislabelledLibs(cmpps[mpii], misls, 10);
        if(misls.size() > 0)
        {
          // PrintOutput] assumes rotation to genomic axis
          /*
          PrintOutput(misls, (MatePairIndex_e) mpii,
                      mpps[mpii], mppsMap[mpii],
                      assembly, seqID, "mislabeled", numStddevs,
                      &relativeID, listOS, printATA, printGnuplot);
          */
        }
      }
      
      if(mpii == MPI_STRETCHED ||
         mpii == MPI_COMPRESSED ||
         mpii == MPI_INVERSION)
      {
        if(mpii == MPI_STRETCHED || mpii == MPI_INVERSION)
        {
          // RefineWithSatisfied assumes rotation to genomic axis
          if(mpii == MPI_STRETCHED)
            RefineWithSatisfied(cmpps[mpii], polys, smpsv,
                                (MatePairIndex_e) mpii);
          else
            RefineInversions(cmpps[mpii], polys, smpsv, filterThresh);

          if(polys.size() > 0)
          {
            PrintOutput(polys, (MatePairIndex_e) mpii,
                        mpps[mpii], mppsMap[mpii],
                        assembly, seqID, "polymorphic", numStddevs,
                        &relativeID, listOS, printATA, printGnuplot);
          }
        }

        if(probs.size() > 0)
        {
#ifdef DEBUG_PROCESSMPS
          vector<CompositeMPPolygon<UNIT_TYPE> >::iterator mppIter;
          for(mppIter = probs.begin(); mppIter != probs.end(); mppIter++)
          {
            cerr << "# No intersection for this set of mate pairs:\n";
            cerr << "# Category: " << MatePairLabel[mpii] << endl;
            cerr << "# MBR intersection: " << *mppIter << endl << endl;
          }
#endif
          /*
          {
            PrintOutput(probs, (MatePairIndex_e) mpii,
                        mpps[mpii], mppsMap[mpii],
                        assembly, seqID, "problematic", numStddevs,
                        &relativeID, listOS, printATA, printGnuplot);
          }
          */
        }
        
        if(cmpps[mpii].size() > 0)
        {
          PrintOutput(cmpps[mpii], (MatePairIndex_e) mpii,
                      mpps[mpii], mppsMap[mpii],
                      assembly, seqID, "confirmed", numStddevs,
                      &relativeID, listOS, printATA, printGnuplot);
        }
      }
    }
    cout << probs.size() << " problematic " << MatePairLabel[mpii]
         << " mate pair sets\n";
    cout << cmpps[mpii].size() << " confirmed " << MatePairLabel[mpii]
         << " mate pair sets\n";
    cout << fews.size() << " below-threshold "
         << MatePairLabel[mpii]
         << " mate pair sets (<"
         << filterThresh << ")\n";
  }

  DetectTranspositions(mpps[MPI_COMPRESSED], // not rotated, raw
                       mpps[MPI_STRETCHED], // not rotated, raw
                       cmpps[MPI_OUTTIE], // not rotated, clustered
                       cmpps[MPI_TRANSPOSITION], // not rotated, empty
                       libs,
                       numStddevs,
                       filterThresh);
  cout << cmpps[MPI_TRANSPOSITION].size() << " confirmed "
       << MatePairLabel[MPI_TRANSPOSITION]
       << " mate pair sets\n";
  if(cmpps[MPI_TRANSPOSITION].size() > 0)
  {
    PrintOutput(cmpps[MPI_TRANSPOSITION], MPI_TRANSPOSITION,
                mpps[MPI_TRANSPOSITION], mppsMap[MPI_TRANSPOSITION],
                assembly, seqID, "confirmed", numStddevs,
                &relativeID, listOS, printATA, printGnuplot);

    // identify stretched & compressed double-counted in transpositions
    vector<CompositeMPPolygon<UNIT_TYPE> >::iterator mppIter;
    map<ID_TYPE, int> inTransps;
    for(mppIter = cmpps[MPI_TRANSPOSITION].begin();
        mppIter != cmpps[MPI_TRANSPOSITION].end();
        mppIter++)
    {
      for(unsigned int i = 0; i < mppIter->getNumMPs(); i++)
      {
        if(mppIter->getMP(i).getOrientation() == PAIR_INNIE)
          inTransps[mppIter->getMP(i).getLeftFragUID()] = 0;
      }
    }

    for(int q = 0; q < 2; q++)
    {
      MatePairIndex_e mpii = (q == 0 ? MPI_STRETCHED : MPI_COMPRESSED);
      int doubleCounted = 0;
      for(mppIter = cmpps[mpii].begin();
          mppIter != cmpps[mpii].end();
          mppIter++)
      {
        for(unsigned int i = 0; i < mppIter->getNumMPs(); i++)
        {
          if(inTransps.find(mppIter->getMP(i).getLeftFragUID()) !=
             inTransps.end())
          {
            doubleCounted++;
            break;
          }
        }
      }
      cout << doubleCounted << " confirmed "
           << (mpii == MPI_STRETCHED ? "insertions " : "deletions ")
           << "used to confirm transpositions\n";
    }
  }
  listOS.close();

  return 0;
}
