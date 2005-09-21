
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
/* $Id: compareAssemblyMPs.cc,v 1.5 2005-09-21 20:13:07 catmandew Exp $ */
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
#include "AssessedMatePair.h"

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
  Read in all raw & processed mate pairs from two assemblies & compare

  INPUT:
    Intra-sequence
      Read raw mate pairs
        separate by category. store whole record in hashtables key = leftUID
      Read confirmed unsatisfied mate pairs
        keep separate categories. store in arrays of MatePairGroup
          store pointers to each MPG in hashtable key = leftUID

    Inter-sequence
      Read confirmed mate pairs. store MatePairGroups in array
        store pointers to each MPG in hashtable key = leftUID

  PROCESSING:
    For each confirmed unsatisfied intra-sequence mate pair group
      look at status in other assembly:
        # of agreeing satisfieds & location
        # of confirmed dissatisfied by type & location
        of the rest,
        # of raw satisfieds & location
        # of raw dissatisfied by type & location

    For each confirmed inter-sequence mate pair group
      look at status in other assembly:
        # of agreeing, satisfied, intra-sequence mps & location
        # of confirmed, dissatisfied intra-sequence...
        of the rest
        # of raw satisfied intra
        # of raw dissatisfied intra
        # of raw inter-
        # of confirmed inter-
 */

void Usage(char * progname, char * message)
{
  if(message != NULL)
    cerr << endl << message << endl;

  cerr << "Usage: " << progname << " [-l libName] [-n numStddevs] [-d basePath] [-s species] [-a assemblyName]*\n";
  cerr << "\t-l lib        name of clone library file\n";
  cerr << "\t-n #stddevs   number of stddevs from mean that is excessive\n";
  cerr << "\t                default is " << STDDEVS_THRESHOLD << endl;
  cerr << "\t-d dir          path to assembly data\n";
  cerr << "\t                 default is " << DefaultBaseDir << endl;
  cerr << "\t-s species      species of assemblies\n";
  cerr << "\t                 default is " << DefaultSpecies << endl;
  cerr << "\t-a assembly     assembly name - multiple are required\n";

  cerr << "\n";
  cerr << "So, .ata files for intra-sequence analysis of a given assembly\n";
  cerr << "would be in basePath/species/assembly/intraSequence\n\n";
  exit(1);
}


void ReadIntraMPs(vector<AssessedMatePair> & amps,
                  ID_TYPE seqID,
                  map<ID_TYPE, CloneLibrary> & libs,
                  double numStddevs,
                  ifstream & fin)
{
  char line[4096];
  AssessedMatePair amp;

  while(fin.getline(line, 4095))
  {
    bool skip = false;
    PairOrientation_e pairOrient;
    char orient;
    ID_TYPE leftUID, rightUID, libUID;
    UNIT_TYPE left5, right5;
    sscanf(line, "%c " F_U64 " " F_U64 " " F_U64 " %d %d",
           &orient, &leftUID, &rightUID, &libUID, &left5, &right5);
    switch(orient)
    {
      case 'I':
        pairOrient = PAIR_INNIE;
        break;
      case 'O':
        pairOrient = PAIR_OUTTIE;
        break;
      case 'A':
        pairOrient = PAIR_ANTINORMAL;
        break;
      case 'N':
        pairOrient = PAIR_NORMAL;
        break;
      default:
        skip = true;
        break;
    }
    if(skip) continue;
      
    map<ID_TYPE, CloneLibrary>::iterator iter;
    iter = libs.find(libUID);
    amp.set(pairOrient, leftUID, rightUID, libUID, left5, right5, seqID,
            (UNIT_TYPE) ((*iter).second.getMean() +
                          numStddevs * (*iter).second.getStddev()),
            (UNIT_TYPE) ((*iter).second.getMean() -
                          numStddevs * (*iter).second.getStddev()));
    amps.push_back(amp);
  }
}


void ReadInterMPs(vector<AssessedMatePair> & amps,
                  ifstream & fin)
{
  char line[4096];
  AssessedMatePair amp;

  while(fin.getline(line, 4095))
  {
    char leftOrient[4], rightOrient[4];
    ID_TYPE leftUID, rightUID, libUID;
    ID_TYPE leftSeqID, rightSeqID;
    UNIT_TYPE left5, right5;
    
    sscanf(line, F_U64 " %d %d %s " F_U64 " %d %d %s " F_U64,
           &leftUID, &leftSeqID, &left5, leftOrient,
           &rightUID, &rightSeqID, &right5, leftOrient,
           &libUID);

    PairOrientation_e pairOrient;
    if(strcmp(leftOrient, "A_B") == 0)
    {
      if(strcmp(rightOrient, "A_B") == 0)
      {
        pairOrient = PAIR_NORMAL;
      }
      else
      {
        pairOrient = PAIR_INNIE;
      }
    }
    else
    {
      if(strcmp(rightOrient, "A_B") == 0)
      {
        pairOrient = PAIR_OUTTIE;
      }
      else
      {
        pairOrient = PAIR_ANTINORMAL;
      }
    }
      
    amp.set(pairOrient,
            leftUID, rightUID, libUID,
            left5, right5,
            leftSeqID, rightSeqID,
            0, 0);
    
    amps.push_back(amp);
  }
}


void ReadAssemblyMatePairs(vector<AssessedMatePair> & amps,
                           const char * baseDir,
                           const char * species,
                           const string & assembly,
                           map<ID_TYPE, CloneLibrary> & libs,
                           double numStddevs)
{
  char filename[4096];
  int numMPs = 0;

  // loop over all #.txt files & read into mps
  for(ID_TYPE seqID = 0; seqID < 24; seqID++)
  {
    sprintf(filename, "%s/%s/%s/%s/%03d.txt",
            baseDir, species, assembly.c_str(), IntraDir, seqID);
    ifstream fmp(filename, ios::in);
    
    // cerr << "Opening " << filename << ".\n";
    if(!fmp.good())
      continue;

    cerr << "Reading raw mate pairs from " << filename << ".\n";
    ReadIntraMPs(amps, seqID, libs, numStddevs, fmp);
    fmp.close();
    cerr << "  " << amps.size() - numMPs << " mps\n";
    numMPs = amps.size();
  }

  for(ID_TYPE seqID = 24; ; seqID++)
  {
    sprintf(filename, "%s/%s/%s/%s/unmapped/%03d.txt",
            baseDir, species, assembly.c_str(), IntraDir, seqID);
    ifstream fmp(filename, ios::in);
    
    // cerr << "Opening " << filename << ".\n";
    if(!fmp.good())
      break;

    cerr << "Reading raw mate pairs from " << filename << ".\n";
    ReadIntraMPs(amps, seqID, libs, numStddevs, fmp);
    fmp.close();
    cerr << "  " << amps.size() - numMPs << " mps\n";
    numMPs = amps.size();
  }

  // read inter-sequence matepairs
  for(ID_TYPE seqID = 0; seqID < 24; seqID++)
  {
    sprintf(filename, "%s/%s/%s/%s/%03d.txt",
            baseDir, species, assembly.c_str(), InterDir, seqID);
    ifstream fmp(filename, ios::in);
    
    // cerr << "Opening " << filename << ".\n";
    if(!fmp.good())
      continue;

    cerr << "Reading raw mate pairs from " << filename << ".\n";
    ReadInterMPs(amps, fmp);
    fmp.close();
    cerr << "  " << amps.size() - numMPs << " mps\n";
    numMPs = amps.size();
  }
  for(ID_TYPE seqID = 24; ; seqID++)
  {
    sprintf(filename, "%s/%s/%s/%s/unmapped/%03d.txt",
            baseDir, species, assembly.c_str(), InterDir, seqID);
    ifstream fmp(filename, ios::in);
    
    // cerr << "Opening " << filename << ".\n";
    if(!fmp.good())
      break;

    cerr << "Reading raw mate pairs from " << filename << ".\n";
    ReadInterMPs(amps, fmp);
    fmp.close();
    cerr << "  " << amps.size() - numMPs << " mps\n";
    numMPs = amps.size();
  }
}


void CompareRawUnsatisfiedMPs(vector<AssessedMatePair> & rawMPs,
                              map<ID_TYPE, unsigned int> & rawMap,
                              char * baseDir,
                              char * species,
                              vector<string> & assemblies)
{
  /*
    for each assembly beyond 0
      for each sequence
        for each unsatisfied category:
        (stretched, compressed, outtie, normal, antinormal)
          for each mate pair in raw file
            look up in rawMPs via rawMap & print
   */
  for(unsigned int aIndex = 1; aIndex < assemblies.size(); aIndex++)
  {
    for(int seqID = 0; seqID < 24; seqID++)
    {
      int numUnsatisfied = 0;
      int numDiffSeqID = 0;
      int numUnmapped = 0;
      int numInterSeqID = 0;
      int numAbsent = 0;
      int numSeen = 0;
      for(int mpi = 0; mpi < MPI_INVERSION; mpi++)
      {
        char filename[4096];
        sprintf(filename, "%s/%s/%s/%s/%s.%03d.%s.raw",
                baseDir, species, assemblies[aIndex].c_str(), IntraDir,
                assemblies[aIndex].c_str(), seqID,
                MatePairLabel[mpi]);
        ifstream fin(filename, ios::in);
    
        if(!fin.good())
          continue;

        char line[4096];
        MatePair rmp;
        while(fin.getline(line, 4095))
        {
          rmp.setFromString(line);

          cout << "raw ( " << rmp.getLeftFragUID()
               << " " << rmp.getRightFragUID() << " ): ";

          cout << "( " << assemblies[aIndex] << " "
               << seqID << " "
               << MatePairLabel[mpi] << " "
               << rmp.getLeftCoord() << " "
               << rmp.getRightCoord() - rmp.getLeftCoord() << " ) ";

          // look up in map
          map<ID_TYPE, unsigned int>::iterator iter;
          iter = rawMap.find(rmp.getLeftFragUID());
          if(iter == rawMap.end())
            iter = rawMap.find(rmp.getRightFragUID());
          
          if(iter == rawMap.end())
          {
            cout << "not in " << assemblies[0] << endl;
            numAbsent++;
          }
          else
          {
            AssessedMatePair & mp = rawMPs[(*iter).second];
            
            if(mp.getType() == MPI_INTERSEQUENCE)
            {
              numInterSeqID++;
            }
            else
            {
              numUnsatisfied += (mp.getType() != MPI_SATISFIED) ? 1 : 0;
              if(mp.getSequenceID() < 24)
                numDiffSeqID += (mp.getSequenceID() != seqID) ? 1 : 0;
              else
                numUnmapped++;
            }
            
            if(mp.getSequenceID() == -1)
            {
              cout << "( " << assemblies[0] << " "
                   << mp.getLeftSequenceID() << "_"
                   << mp.getRightSequenceID() << " "
                   << MatePairLabel[mp.getType()] << " "
                   << mp.getLeftCoord() << " "
                   << mp.getRightCoord() - mp.getLeftCoord() << " )"
                   << endl;
            }
            else
            {
              cout << "( " << assemblies[0] << " "
                   << mp.getSequenceID() << " "
                   << MatePairLabel[mp.getType()] << " "
                   << mp.getLeftCoord() << " "
                   << mp.getRightCoord() - mp.getLeftCoord() << " )"
                   << endl;
            }
          }
        }
        cout << MatePairLabel[mpi] << " - "
             << " unsatisfied: " << numUnsatisfied
             << " diffSeqID: " << numDiffSeqID
             << " interSeqID: " << numInterSeqID
             << " unmapped: " << numUnmapped
             << " absent: " << numAbsent << endl;
      }
    }
  }
}


void CompareUnsatisfiedMPGs(vector<AssessedMatePair> & rawMPs,
                            map<ID_TYPE, unsigned int> & rawMap,
                            char * baseDir,
                            char * species,
                            vector<string> & assemblies)
{
  /*
    for each assembly beyond 0
      for each sequence
        for each unsatisfied category:
        (stretched, compressed, inversion, transposition)
          for each mate pair group in ata file
            look up in rawMPs via rawMap & print
   */
  vector<MatePairIndex_e> mpIndices;
  mpIndices.push_back(MPI_STRETCHED);
  mpIndices.push_back(MPI_COMPRESSED);
  mpIndices.push_back(MPI_INVERSION);
  mpIndices.push_back(MPI_TRANSPOSITION);
  for(unsigned int aIndex = 1; aIndex < assemblies.size(); aIndex++)
  {
    for(int seqID = 0; seqID < 24; seqID++)
    {
      for(unsigned int mpiIndex = 0; mpiIndex < mpIndices.size(); mpiIndex++)
      {
        char filename[4096];
        sprintf(filename, "%s/%s/%s/%s/%s.%03d.%s.ata",
                baseDir, species, assemblies[aIndex].c_str(), IntraDir,
                assemblies[aIndex].c_str(), seqID,
                MatePairLabel[mpIndices[mpiIndex]]);
        ifstream fin(filename, ios::in);
    
        // cerr << "Opening " << filename << ".\n";
        if(!fin.good())
          continue;

        char line[4096];
        while(fin.getline(line, 4095))
        {
          // cerr << line << endl;
          // get to next mate pair group
          if(strstr(line, "weight") == NULL) continue;

          char mpgName[50];
          UNIT_TYPE start, length;
          int weight;
          sscanf(line, "%*c %*s %s . %*s %d %d %*d %*d . . > /weight=%d",
                 mpgName, &start, &length, &weight);

          int numUnsatisfied = 0;
          int numDiffSeqID = 0;
          int numUnmapped = 0;
          int numPartUnmapped = 0;
          int numInterSeqID = 0;
          int numAbsent = 0;
          int numSeen = 0;
          while(numSeen < weight && fin.getline(line, 4095))
          {
            if(line[2] == 'm')
            {
              numSeen++;

              // get left uid
              UNIT_TYPE start, length;
              ID_TYPE leftUID, rightUID;
              sscanf(line, "%*c %*s %*s %*s %*s %d %d . . . . > /leftUID=" F_U64 " > /rightUID=" F_U64,
                     &start, &length, &leftUID, &rightUID);

              cout << "( " << leftUID << " " << rightUID << " ): ";

              cout << "( " << assemblies[aIndex] << " "
                   << seqID << " "
                   << MatePairLabel[mpIndices[mpiIndex]] << " "
                   << start << " "
                   << length << " ) ";

              // look up in map
              map<ID_TYPE, unsigned int>::iterator iter;
              iter = rawMap.find(leftUID);
              if(iter == rawMap.end())
                iter = rawMap.find(rightUID);
              
              if(iter == rawMap.end())
              {
                cout << "not in " << assemblies[0] << endl;
                numAbsent++;
              }
              else
              {
                AssessedMatePair & mp = rawMPs[(*iter).second];

                if(mp.getType() == MPI_INTERSEQUENCE)
                {
                  numInterSeqID++;
                }
                else
                {
                  numUnsatisfied += (mp.getType() != MPI_SATISFIED) ? 1 : 0;
                  if(mp.getSequenceID() < 24)
                    numDiffSeqID += (mp.getSequenceID() != seqID) ? 1 : 0;
                  else
                    numUnmapped++;
                }

                if(mp.getSequenceID() == -1)
                {
                  cout << "( " << assemblies[0] << " "
                       << mp.getLeftSequenceID() << "_"
                       << mp.getRightSequenceID() << " "
                       << MatePairLabel[mp.getType()] << " "
                       << mp.getLeftCoord() << " "
                       << mp.getRightCoord() - mp.getLeftCoord() << " )"
                       << endl;
                }
                else
                {
                  cout << "( " << assemblies[0] << " "
                       << mp.getSequenceID() << " "
                       << MatePairLabel[mp.getType()] << " "
                       << mp.getLeftCoord() << " "
                       << mp.getRightCoord() - mp.getLeftCoord() << " )"
                       << endl;
                }
              }
            }
          }
          cout << MatePairLabel[mpIndices[mpiIndex]] << " MGP: "
               << start << " , " << length 
               << " weight: " << weight
               << " unsatisfied: " << numUnsatisfied
               << " diffSeqID: " << numDiffSeqID
               << " interSeqID: " << numInterSeqID
               << " unmapped: " << numUnmapped
               << " absent: " << numAbsent << endl;
        }
      }
    }
  }
}


int main(int argc, char ** argv)
{
  vector<string> assemblies;
  char * baseDir = NULL;
  char * species = NULL;
  double numStddevs = STDDEVS_THRESHOLD;
  char * libFilename = NULL;

  int ch, errflg = 0;
  while(!errflg && ((ch = getopt(argc, argv, "l:d:s:a:n:")) != EOF))
  {
    switch(ch)
    {
      case 'l':
        libFilename = optarg;
        break;
      case 'd':
        baseDir = optarg;
        break;
      case 's':
        species = optarg;
        break;
      case 'a':
        assemblies.push_back(optarg);
        break;
      case 'n':
        numStddevs = atof(optarg);
        break;
      default:
        errflg++;
        break;
    }
  }

  if(assemblies.size() < 2)
    Usage(argv[0], "Please specify two or more assemblies.");
  if(libFilename == NULL)
    Usage(argv[0], "Please specify a library filename.");
  if(numStddevs <= 0)
    Usage(argv[0], "Please specify a positive number of std deviations");
  
  if(baseDir == NULL)
  {
    cerr << "Using default directory: " << DefaultBaseDir << endl;
    baseDir = DefaultBaseDir;
  }
  if(species == NULL)
  {
    cerr << "Using default species: " << DefaultSpecies << endl;
    species = DefaultSpecies;
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

  // read in raw intra-sequence mate pairs of first assembly
  vector<AssessedMatePair> rawMPs;
  ReadAssemblyMatePairs(rawMPs, baseDir, species, assemblies[0],
                        libs, numStddevs);
  cerr << "Assembly " << assemblies[0] << " has "
       << rawMPs.size() << " raw mate pairs.\n";

  // create hashtable of raw mate pairs
  map<ID_TYPE, unsigned int > rawMap;
  for(unsigned int i = 0; i < rawMPs.size(); i++)
  {
    rawMap[rawMPs[i].getLeftFragUID()] = i;
  }

  CompareRawUnsatisfiedMPs(rawMPs, rawMap, baseDir, species, assemblies);
  
  // read in intra-sequence unsatisfied mate pair groups of all assemblies
  // evaluate while reading in
  CompareUnsatisfiedMPGs(rawMPs, rawMap, baseDir, species, assemblies);
  
  return 0;
}
