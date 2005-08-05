
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
/* $Id: processInterCG.cc,v 1.5 2005-08-05 00:56:41 catmandew Exp $ */
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

#define CHROM_OFFSET 500000000

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
  "interChromosome",
  "unknown"
};


/*
  PURPOSE:
    Process inter-chromosome mate pairs - mate pairs for which each fragment
    was mapped to a different chromosome - to identify mis-assemblies and
    perhaps re-sequenced plates or collapsed repeats

  PLAN:
    Read in clone library length estimates
    
    Read in all intra-chromosome mate pairs for a given chromosome/assembly
      store each mate pair twice (by left and by right frag coord)
      sort by 'left' coorinate
      
    Read in all inter-chromosome mate pairs for a given chromosome/assembly
      Maintain separate list of mates for each other chromosome
        optionally omit lower numbered chromosomes to avoid duplication
      Sort each list by coordinate on other chromosome

    c = this chromosome
    For each other chromosome, o
      for each mate, m, on that chromosome
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
  cerr << "Usage: " << progname << " [-l lib] [-e filename] [-n #] [-f #] -g\n";
  cerr << "\t-l lib        name of clone library file\n";
  cerr << "\t-e filename   name of inter-chromosome matepairs file\n";
  cerr << "\t                filename must have form:\n";
  cerr << "\t                assemblyName_#_inter.txt\n";
  cerr << "\t                where # is the scaffold or chromosome number\n";
  /*
  cerr << "\t-a assembly   assembly name (e.g., B33A, VAN)\n";
  cerr << "\t-c chrom      chromosome number (index into fasta file)\n";
  */
  cerr << "\t-n #          number of stddevs from mean that is excessive\n\n";
  cerr << "\t                default is " << STDDEVS_THRESHOLD << endl;
  cerr << "\t-f #          filter out mate pair sets with fewer members\n";
  cerr << "\t                default is " << CONFIRMATION_THRESHOLD << endl;
  cerr << "\t-g            print gnuplot output instead of ata output\n";
  
  cerr << endl;
  exit(1);
}


void ReadInterChromosomeMPs(vector<MatePair> & mps,
                            ifstream & fin,
                            int filterChrom)
{
  char line[LINE_SIZE];
  while(fin.getline(line, LINE_SIZE-1))
  {
    int otherChrom;
    sscanf(line, "%*s %*d %*d %*s %*s %d", &otherChrom);
    if(otherChrom == filterChrom)
    {
      MatePair mp;
      mp.setFromInterChromosomeString(line);
      mp.setRightCoord(mp.getRightCoord() + CHROM_OFFSET);
      mps.push_back(mp);
    }
  }
}


int main(int argc, char ** argv)
{
  char * libFilename = NULL;
  char * ewFilename = NULL;
  // char * assembly = NULL;
  // int chromosome = -1;
  char assembly[4096];
  char chromosome[4096];
  double numStddevs = STDDEVS_THRESHOLD;
  int filterThresh = CONFIRMATION_THRESHOLD;
  bool printForGnuplot = false;

  {
    int ch, errflg = 0;
    while(!errflg && ((ch = getopt(argc, argv, "l:e:a:c:n:f:g")) != EOF))
    {
      switch(ch)
      {
        case 'l':
          libFilename = optarg;
          break;
        case 'e':
          ewFilename = optarg;
          break;
        /*
        case 'a':
          assembly = optarg;
          break;
        case 'c':
          chromosome = atoi(optarg);
          break;
        */
        case 'n':
          numStddevs = atof(optarg);
          break;
        case 'f':
          filterThresh = atoi(optarg);
          break;
        case 'g':
          printForGnuplot = true;
          break;
        default:
          errflg++;
          break;
      }
    }
    if(libFilename == NULL)
      Usage(argv[0], "Please specify a clone library filename");
    if(ewFilename == NULL)
      Usage(argv[0], "Please specify an inter-chromosome mate pair filename");
    if(numStddevs <= 0)
      Usage(argv[0], "Please specify a positive number of std deviations");
    /*
    if(assembly == NULL)
      Usage(argv[0], "Please specify an assembly name");
    if(chromosome < 0)
      Usage(argv[0], "Please specify a chromosome number");
    */
    {
      char * ptr;
      strcpy(assembly, ewFilename);
      ptr = index(assembly, (int) '_');
      assert(ptr != NULL);
      ptr[0] = '\0';
      strcpy(chromosome, ++ptr);
      ptr = index(chromosome, (int) '_');
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
  set<int> otherChromosomes;
  char line[LINE_SIZE];
  ifstream fe(ewFilename, ios::in);
  if(!fe.good())
  {
    cerr << "Failed to open " << ewFilename << " for reading\n";
    exit(-1);
  }
  while(fe.getline(line, LINE_SIZE-1))
  {
    int otherChrom;
    sscanf(line, "%*s %*d %*d %*s %*s %d", &otherChrom);
    otherChromosomes.insert(otherChrom);
  }
  fe.close();

  // open output file
  char fname[1024];
  ofstream fo;
  int atacCount = atoi(chromosome) * 10000;
  if(printForGnuplot)
  {
    sprintf(fname, "%s.%s.interChromosome.gp", assembly, chromosome);
    fo.open(fname, ios::out);
    if(!fo.good())
    {
      cerr << "Failed to open " << fname << " for writing\n";
      exit(-1);
    }
  }
  else
  {
    sprintf(fname, "%s.%s.interChromosome.ata", assembly, chromosome);
    fo.open(fname, ios::out);
    if(!fo.good())
    {
      cerr << "Failed to open " << fname << " for writing\n";
      exit(-1);
    }
    // write ata header lines
    fo << "! format ata 1.0\n";
    fo << "# numStddevs=" << numStddevs << endl;
  }
  
  // iterate over all other chromosomes in the set
  set<int>::iterator siter;
  for(siter = otherChromosomes.begin();
      siter != otherChromosomes.end();
      siter++)
  {
    int otherChrom = *siter;
    if(otherChrom >= atoi(chromosome))
      break;

    if(printForGnuplot)
      fo << "# " << otherChrom << endl;
    
    vector<MatePair> mps;
    ifstream fel(ewFilename, ios::in);
    ReadInterChromosomeMPs(mps, fel, otherChrom);
    fel.close();
    
    list<MatePair> mpl[MPI_NUM_INDICES];
    int badLibMatePairs = 0;
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
        badLibMatePairs++;
      }
    }
    // cerr << "Omitted " << badLibMatePairs << " matepair(s) from bad clone libraries\n";

    // process stretched, outtie, normal, & antinormal
    vector<CompositeMPPolygon<UNIT_TYPE> > mpps[MPI_NUM_INDICES]; // input
    map<uint64, MatePairPolygon<UNIT_TYPE> > mppsMap[MPI_NUM_INDICES];
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
            mpp[j].setY(mpp[j].getY() - CHROM_OFFSET);
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
                                       CHROM_OFFSET);
        }
        
        for(int i = 0; i < cmpps[mpii].size(); i++)
        {
          cmpps[mpii][i].rotateByDegrees(-45);
          for(int j = 0; j < cmpps[mpii][i].size(); j++)
            for(int k = 0; k < cmpps[mpii][i][j].size(); k++)
              cmpps[mpii][i][j][k].setY(cmpps[mpii][i][j][k].getY() -
                                        CHROM_OFFSET);


          if(printForGnuplot)
          {
            cmpps[mpii][i].printForGnuplot(fo);
            for(int j = 0; j < cmpps[mpii][i].getNumMPs(); j++)
              mppsMap[mpii][cmpps[mpii][i].getMP(j).getLeftFragUID()].printForGnuplot(fo);
          }
          else
          {
            int leftThis, rightThis, leftOther, rightOther;
            UNIT_TYPE val =
              cmpps[mpii][i].getMP(0).getRightCoord() - CHROM_OFFSET;
            leftThis = cmpps[mpii][i].getMP(0).getLeftCoord();
            rightThis = leftThis;
            leftOther = rightOther = val;
            for(int j = 1; j < cmpps[mpii][i].getNumMPs(); j++)
            {
              val = cmpps[mpii][i].getMP(j).getRightCoord() - CHROM_OFFSET;
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
                fo << "M HL HL" << atacCount++ << " . "
                   << assembly << ":" << chromosome << " "
                   << leftThis << " " << rightThis - leftThis << " 1 "
                   << assembly << ":" << otherChrom << " "
                   << leftOther << " " << rightOther - leftOther << " 1 "
                   << "> /weight=" << cmpps[mpii][i].getNumMPs() << endl;
                break;
              case MPI_NORMAL:
                fo << "M Hl Hl" << atacCount++ << " . "
                   << assembly << ":" << chromosome << " "
                   << leftThis << " " << rightThis - leftThis << " 1 "
                   << assembly << ":" << otherChrom << " "
                   << leftOther << " " << rightOther - leftOther << " -1 "
                   << "> /weight=" << cmpps[mpii][i].getNumMPs() << endl;
                break;
              case MPI_ANTINORMAL:
                fo << "M hL hL" << atacCount++ << " . "
                   << assembly << ":" << chromosome << " "
                   << leftThis << " " << rightThis - leftThis << " -1 "
                   << assembly << ":" << otherChrom << " "
                   << leftOther << " " << rightOther - leftOther << " 1 "
                   << "> /weight=" << cmpps[mpii][i].getNumMPs() << endl;
                break;
              case MPI_OUTTIE:
                fo << "M hl hl" << atacCount++ << " . "
                   << assembly << ":" << chromosome << " "
                   << leftThis << " " << rightThis - leftThis << " -1 "
                   << assembly << ":" << otherChrom << " "
                   << leftOther << " " << rightOther - leftOther << " -1 "
                   << "> /weight=" << cmpps[mpii][i].getNumMPs() << endl;
                break;
              default:
                cerr << "Unknown matepair type: " << mpii << endl;
                exit(-100);
                break;
            }
          }
        } // loop over composite polygons
      } // if there are matepairs
    } // loop over matepair types
  } // loop over other chromosomes
  fo.close();
  
  return 0;
}
