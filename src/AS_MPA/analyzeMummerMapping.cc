
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
/* $Id: analyzeMummerMapping.cc,v 1.1 2005-11-18 23:27:52 catmandew Exp $ */
#include <cstdio>  // for sscanf
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <unistd.h> /* man 3 getopt */

//#define DEBUG_INSTANCE
//#define VERBOSE

using namespace std;

#include "cds.h"

/*
  Purpose: Analyze & summarize output of show-coords program, which itself
    produces readable output from nucmer/mummer

  Inputs:
    .scaff file from reference assembly
    .scaff file from query assembly
    show-coords output file

  Outputs:
    to stdout (initially)

  Algorithm:
    1. create 2 assembly objects
    2. read .scaff file into each
    3. read show-coords file & keep track of best matches
    4. for each assembly
         for each scaffold in the assembly
           check consistency of consecutive matches between assemblies
           identify unmatched sequence
       
 */

typedef enum
{
  Reference = 0,
  Query,
  NumSequenceTypes
} SequenceType;


typedef enum
{
  Forward = 0,
  Reverse,
  NumDirections
} DirectionType;

class Interval
{
public:
  Interval() {set(0,0,0);}
  Interval(CDS_UID_t uid, CDS_COORD_t begin, CDS_COORD_t end)
    {
      set(uid, begin, end);
    }

  void set(CDS_UID_t uid, CDS_COORD_t begin, CDS_COORD_t end)
    {
      setUID(uid);
      setCoords(begin, end);
    }

  void setCoords(CDS_COORD_t begin, CDS_COORD_t end)
    {
      _begin = begin;
      _end = end;
      _min = (_begin < _end ? _begin : _end);
      _max = (_begin > _end ? _begin : _end);
    }
  
  void setUID(CDS_UID_t uid) {_uid = uid;}
  void setBegin(CDS_COORD_t begin) {setCoords(begin, getEnd());}
  void setEnd(CDS_COORD_t end) {setCoords(getBegin(), end);}

  CDS_UID_t getUID() const {return _uid;}
  CDS_COORD_t getBegin() const {return _begin;}
  CDS_COORD_t getEnd() const {return _end;}
  CDS_COORD_t getMin() const {return _min;}
  CDS_COORD_t getMax() const {return _max;}
  bool isVoid() const
    {
      return (getUID() == 0 && getBegin() == 0 && getEnd() == 0);
    }

  bool isForward() const {return getBegin() <= getEnd();}
  bool isReverse() const {return !isReverse();}
  DirectionType getDirection() const
    {
      return (isForward() ? Forward : Reverse);
    }
  
  bool intersects(Interval interval) const
    {
      if(interval.getUID() == getUID() &&
         interval.getMin() <= getMax() &&
         interval.getMax() >= getMin())
        return true;
      return false;
    }
  
private:
  CDS_UID_t _uid;
  CDS_COORD_t _begin;
  CDS_COORD_t _end;
  CDS_COORD_t _min;
  CDS_COORD_t _max;
};

class Match
{
public:
  Match() {_onSeq[Reference].set(0,0,0); _onSeq[Query].set(0,0,0); _pctID = 0;}
  Match(char * line) {set(line);}

  void set(char * line)
    {
      /* tab-delimited line with fields:
       +  1. begin in ref scaffold
       +  2. end in ref scaffold
       +  3. begin in query scaffold
       +  4. end in query scaffold
          5. match length in ref
          6. match length in query
       +  7. % identity of match
          8. total length of scaffold in ref
          9. total length of scaffold in query
         10. % coverage of scaffold in ref
         11. % coverage of scaffold in query
       + 12. scaffold UID in ref
       + 13. scaffold UID in query
      */
      int b1, b2, e1, e2;
      CDS_UID_t uid1, uid2;
      sscanf(line,
             "%d\t%d\t%d\t%d\t%*d\t%*d\t%f\t%*d\t%*d\t%*f\t%*f\t" F_UID "\t" F_UID,
             &b1, &e1, &b2, &e2, &_pctID, &uid1, &uid2);
      // switch to 0-offset inter-bp based coordinates
      b1--;
      b2--;
      _onSeq[0].set(uid1,b1,e1);
      _onSeq[1].set(uid2,b2,e2);
    }

  bool isVoid() const
    {
      return (_onSeq[0].isVoid() && _onSeq[1].isVoid() && getPctID() == 0);
    }
  
  CDS_UID_t getUID(SequenceType which) const {return _onSeq[which].getUID();}
  CDS_COORD_t getBegin(SequenceType which) const {return _onSeq[which].getBegin();}
  CDS_COORD_t getEnd(SequenceType which) const {return _onSeq[which].getEnd();}
  CDS_COORD_t getMin(SequenceType which) const {return _onSeq[which].getMin();}
  CDS_COORD_t getMax(SequenceType which) const {return _onSeq[which].getMax();}
  float getPctID() const {return _pctID;}
  DirectionType getDirection(SequenceType which) const
    {
      return (getBegin(which) == getMin(which) ? Forward : Reverse);
    }
  bool isFlipped() const
    {
      return (getDirection(Reference) != getDirection(Query));
    }
  DirectionType getDirection() const
    {
      return (isFlipped() ? Reverse : Forward);
    }

  void reversePerspective()
    {
      Interval dummy = _onSeq[0];
      _onSeq[0] = _onSeq[1];
      _onSeq[1] = dummy;
    }
  
  void setUID(SequenceType which, CDS_UID_t uid) {_onSeq[which].setUID(uid);}
  void setBegin(SequenceType which, CDS_COORD_t begin) {_onSeq[which].setBegin(begin);}
  void setEnd(SequenceType which, CDS_COORD_t end) {_onSeq[which].setEnd(end);}
  void setPctID(float pctID) {_pctID = pctID;}

  friend ostream & operator<<(ostream & os, const Match & m)
    {
      for(int i = 0; i < NumSequenceTypes; i++)
      {
        os << m.getUID((SequenceType) i) << "\t"
           << m.getBegin((SequenceType) i) << "\t"
           << m.getEnd((SequenceType) i) << "\t";
      }
      os << m.getPctID();
      return os;
    }
  
private:
  Interval _onSeq[2];
  float _pctID;
};


class Scaffold
{
public:
  Scaffold() {reset();}
  Scaffold(CDS_UID_t uid, CDS_COORD_t seqLength, CDS_COORD_t gappedLength)
    {
      set(uid, seqLength, gappedLength);
    }
  Scaffold(char * line)
    {
      set(line);
    }

  int set(char * line)
    {
      int numContigs;
      reset();
      sscanf(line, ">" F_UID " %d %d %d",
             &_uid, &numContigs, &_seqLength, &_gappedLength);
      return numContigs;
    }
  
  void set(CDS_UID_t uid, CDS_COORD_t seqLength, CDS_COORD_t gappedLength)
    {
      _uid = uid;
      _seqLength = seqLength;
      _gappedLength = gappedLength;
      _contigs.clear();
      _gaps.clear();
    }

  void reset()
    {
      set(0,0,0);
    }
  
  void addContig(char * line)
    {
      // Line is space delimited. 2nd field is orientation (ignored)
      CDS_UID_t uid;
      CDS_COORD_t begin = 0;
      CDS_COORD_t end;
      float gapMean;
      sscanf(line, F_UID " %*c%*c %d %f %*f", &uid, &end, &gapMean);

      begin += (getNumGaps() > 0 ? _gaps[getNumGaps() - 1].getEnd() : 0);
      end += begin;

      Interval newContig(uid, begin, end);
      _contigs.push_back(newContig);

      if(fabsf(gapMean) > 0.00001)
      {
        gapMean = (gapMean < 20 ? 20 : gapMean);
        Interval newGap(0, end, (CDS_COORD_t) (end + gapMean + 0.5));
        _gaps.push_back(newGap);
      }
    }

  void considerMatch(Match match, float minIdentity)
    {
      /* given the current list of matches, compare, & do one of the following:
         discard match
         insert (prepend, insert, append)
         replace 1 or more current matches with this one

         first, iterate to match just preceding this one (in coordinates)
           then, iterate to match just after this one
             then,
               if there are intersecting matches
               and all have lower identity
                 replace
               else if there are no intersecting matches
                 insert
       */
      if(match.getPctID() < minIdentity) return;
      
      if(_matches.size() == 0)
      {
        _matches.push_back(match);
        return;
      }

#ifdef DEBUG_INSTANCE
      if(match.getUID(Reference) == 1099520668465 &&
         match.getUID(Query) == 1099520657064)
      {
        cout << "Before adding\n";
        cout << match;
        cout << endl;
        printMatches(cout);
      }
#endif
      
      // find first match that intersects
      list<Match>::iterator minIter;
      for(minIter = _matches.begin();
          minIter != _matches.end() &&
            minIter->getMax(Reference) <= match.getMin(Reference);
          minIter++);

      if(minIter == _matches.end())
      {
#ifdef DEBUG_INSTANCE
        if(match.getUID(Reference) == 1099520668465 &&
           match.getUID(Query) == 1099520657064)
        {
          cout << "push_back\n";
        }
#endif
        _matches.push_back(match);
        return;
      }

      // tally things to compare later & decide which match(es) to keep
      // first check if there actually is an intersection
      int numIntersecting = 0;
      CDS_COORD_t minIntersectCoord = CDS_COORD_MAX;
      CDS_COORD_t maxIntersectCoord = CDS_COORD_MIN;
      CDS_COORD_t bpCovered = 0;
      float sumBPIdentities = 0;
      if(minIter->getMin(Reference) < match.getMax(Reference) &&
         match.getMin(Reference) < minIter->getMax(Reference))
      {
        numIntersecting = 1;
        minIntersectCoord = minIter->getMin(Reference);
        maxIntersectCoord = minIter->getMax(Reference);
        bpCovered = minIter->getMax(Reference) - minIter->getMin(Reference);
        sumBPIdentities =
          (minIter->getPctID() *
           (minIter->getMax(Reference) - minIter->getMin(Reference)));
      }

      // find first match beyond this one
      list<Match>::iterator maxIter = minIter;
      for(++maxIter;
          maxIter != _matches.end() &&
            maxIter->getMin(Reference) < match.getMax(Reference);
          maxIter++)
      {
        numIntersecting++;
        minIntersectCoord = (maxIter->getMin(Reference) < minIntersectCoord ?
                             maxIter->getMin(Reference) : minIntersectCoord);
        maxIntersectCoord = (maxIter->getMax(Reference) > maxIntersectCoord ?
                             maxIter->getMax(Reference) : maxIntersectCoord);
        bpCovered += maxIter->getMax(Reference) - maxIter->getMin(Reference);
        sumBPIdentities +=
          (maxIter->getPctID() *
           (maxIter->getMax(Reference) - maxIter->getMin(Reference)));
      }
      
      if(maxIter == _matches.begin())
      {
#ifdef DEBUG_INSTANCE
        if(match.getUID(Reference) == 1099520668465 &&
           match.getUID(Query) == 1099520657064)
        {
          cout << "push_front\n";
        }
#endif
        _matches.push_front(match);
        return;
      }

      /* if here, there is an intersection to deal with. Need to
         determine whether this match or the intersecting match(es)
         is preferable: age-old problem.

         Since we expect only high-identity matches, select the match(es)
         that cover the most basepairs

         Obviously there is code here to support other heuristics...
      */
      float meanBPIdentity = sumBPIdentities / bpCovered;
      if(bpCovered > (match.getMax(Reference) - match.getMin(Reference)))
      {
#ifdef DEBUG_INSTANCE
        if(match.getUID(Reference) == 1099520668465 &&
           match.getUID(Query) == 1099520657064)
        {
          cout << "rejected\n";
        }
#endif
        return;
      }

      // if here, replace
      _matches.insert(minIter, match);
      if(numIntersecting > 0)
        _matches.erase(minIter, maxIter);
      
#ifdef DEBUG_INSTANCE
      if(match.getUID(Reference) == 1099520668465 &&
         match.getUID(Query) == 1099520657064)
      {
        cout << "After\n";
        printMatches(cout);
      }
#endif
    }

  CDS_UID_t getUID() const {return _uid;}
  CDS_COORD_t getSeqLength() const {return _seqLength;}
  CDS_COORD_t getGappedLength() const {return _gappedLength;}
  
  int getNumContigs() const {return _contigs.size();}
  int getNumGaps() const {return _gaps.size();}
  int getNumMatches() const {return _matches.size();}

  const Interval getContig(int i) const {return _contigs[i];}
  const Interval getGap(int i) const {return _gaps[i];}

  /*
    simple check routine
    are consecutive matches consistent wrt other scaffold & coords
  */
  void checkMatches1(ostream & os, float ratio, int bps) const
    {
      // assumption that this is the 'reference' sequence, relatively
      Match lastMatch;

      list<Match>::const_iterator iter;
      for(iter = _matches.begin(); iter != _matches.end(); iter++)
      {
        if(!lastMatch.isVoid())
        {
          // check if mapping continues in same scaffold
          if(lastMatch.getUID(Query) != iter->getUID(Query))
          {
            os << "Mapping splits scaffold "
               << lastMatch.getUID(Reference) << ":\n";
            os << lastMatch << endl;
            os << *iter << endl << endl;
          }
          else
          {
            // check for a change in orientation
            if(lastMatch.getDirection(Query) != iter->getDirection(Query) ||
               lastMatch.getDirection(Reference) != iter->getDirection(Reference))
            {
              // change in orientation
              os << "Relative inversion breakpoint in scaffold "
                 << lastMatch.getUID(Reference) 
                 << " between "
                 << lastMatch.getMax(Reference) << " and "
                 << iter->getMin(Reference) << ":\n";
              os << lastMatch << endl;
              os << *iter << endl << endl;
            }
            else
            {
              /*
                NOTE: THIS CODE DOESN'T PROPERLY HANDLE INVERSE MAPPINGS
                OF SCAFFOLDS TO EACH OTHER
              */
              
              // check if spacing is about right
              CDS_COORD_t refDelta =
                 iter->getMin(Reference) - lastMatch.getMax(Reference);
              CDS_COORD_t queryDelta =
                (lastMatch.isFlipped() ?
                 lastMatch.getMin(Query) - iter->getMax(Query) :
                 iter->getMin(Query) - lastMatch.getMax(Query));
              CDS_COORD_t diff = refDelta - queryDelta;
              float diffRatio = (queryDelta == 0 ? 1. : diff / queryDelta);
              
              if(diffRatio > ratio || diff > bps)
              {
                os << "Relative insertion in scaffold "
                   << lastMatch.getUID(Reference)
                   << " of " << diff << " bps between "
                   << lastMatch.getMax(Reference) << " and "
                   << iter->getMin(Reference) << ":\n";
                os << lastMatch << endl;
                os << *iter << endl << endl;
              }
              else if(diffRatio < -ratio || diff < -bps)
              {
                os << "Relative deletion in scaffold "
                   << lastMatch.getUID(Reference)
                   << " of " << -diff << " bps between "
                   << lastMatch.getMax(Reference) << " and "
                   << iter->getMin(Reference) << ":\n";
                os << lastMatch << endl;
                os << *iter << endl << endl;
              }
            }
          }
        }
        lastMatch = *iter;
      }
    }
  
  void printMatches(ostream & os) const
    {
      list<Match>::const_iterator iter;
      for(iter = _matches.begin(); iter != _matches.end(); iter++)
      {
        os << *iter << endl;
      }
    }
  
  void printStructure(ostream & os) const
    {
      os << ">" << getUID() << " "
         << getNumContigs() << " "
         << getSeqLength() << " "
         << getGappedLength() << "\n";
      for(int i = 0; i < getNumContigs(); i++)
      {
        const Interval contig = getContig(i);
        os << contig.getUID() << " "
           << (contig.isForward() ? "BE" : "EB") << " "
           << contig.getBegin() << " "
           << contig.getEnd() << " ";
        if(getNumGaps() > i)
        {
          const Interval gap = getGap(i);
          os << gap.getEnd() - gap.getBegin();
        }
        else
          os << "0";
        os << endl;
      }
    }

  friend ostream & operator<<(ostream & os, const Scaffold & s)
    {
      s.printStructure(os);
      os << endl;
      s.printMatches(os);
      return os;
    }
  
private:
  CDS_UID_t _uid;
  CDS_COORD_t _seqLength;
  CDS_COORD_t _gappedLength;
  vector<Interval> _contigs;
  vector<Interval> _gaps;
  list<Match> _matches;
};


// An assembly is simply map<UID, Scaffold>, so no class needed for that

void checkMatches1(ostream & os,
                   map<CDS_UID_t, Scaffold> &assembly,
                   float ratio, int bps)
{
  map<CDS_UID_t, Scaffold>::iterator iter;
  for(iter = assembly.begin(); iter != assembly.end(); iter++)
  {
    ((Scaffold) iter->second).checkMatches1(os, ratio, bps);
  }
}

void printMatches(ostream & os, map<CDS_UID_t, Scaffold> &assembly)
{
  map<CDS_UID_t, Scaffold>::iterator iter;
  for(iter = assembly.begin(); iter != assembly.end(); iter++)
  {
    ((Scaffold) iter->second).printMatches(os);
  }
}

void printScaffs(ostream & os, map<CDS_UID_t, Scaffold> &assembly)
{
  map<CDS_UID_t, Scaffold>::iterator iter;
  for(iter = assembly.begin(); iter != assembly.end(); iter++)
  {
    ((Scaffold) iter->second).printStructure(os);
  }
}

void readScaffFile(char * filename, map<CDS_UID_t, Scaffold> &assembly)
{
  ifstream fin(filename, ios::in);
  if(!fin.good())
  {
    cerr << "Failed to open " << filename << " for reading\n";
    exit(-1);
  }

  cerr << "Reading file " << filename << endl;
  char line[4096];
  while(fin.getline(line, 4095))
  {
    Scaffold scaff;
    int numContigs = scaff.set(line);
    for(int i = 0; i < numContigs; i++)
    {
      fin.getline(line, 4095);
      scaff.addContig(line);
    }
    assembly[scaff.getUID()] = scaff;
  }
  fin.close();
}

void readShowCoordsFile(char * filename,
                        map<CDS_UID_t, Scaffold> assemblies[2],
                        float minIdentity)
{
  ifstream fin(filename, ios::in);
  if(!fin.good())
  {
    cerr << "Failed to open " << filename << " for reading\n";
    exit(-1);
  }

  cerr << "Reading file " << filename << endl;
  char line[4096];
  while(fin.getline(line, 4095))
  {
    Match match(line);
    assemblies[Reference][match.getUID(Reference)].considerMatch(match, minIdentity);
    match.reversePerspective();
    assemblies[Query][match.getUID(Reference)].considerMatch(match, minIdentity);
  }
  fin.close();
}

void Usage(char * progname, char * message)
{
  if(message != NULL)
    cerr << endl << message << endl << endl;;
  cerr << "Usage: " << progname << " [-i minIdentity] -r ref.scaff  -q query.scaff  -s showCoordsFile\n";
  cerr << "\t-r ref.scaff       .scaff file of reference sequence\n";
  cerr << "\t-q query.scaff     .scaff file of query sequence\n";
  cerr << "\t-s showCoordsFile  output file from show-coords for 2 sequences\n";
  cerr << "\t-i minIdentity     minimum pct identity of matches to consider\n";
  cerr << "\t                     default is 95\n";
  cerr << "\nOutput written to stdout\n";
  cerr << "\n\n";
  exit(1);
}


int main(int argc, char ** argv)
{
  char * scaffFilenames[2];
  char * showCoordsFilename = NULL;
  float minIdentity = 95;
  scaffFilenames[Reference] = scaffFilenames[Query] = NULL;

  // process the command line
  {
    int ch, errflg = 0;
    while(!errflg && ((ch = getopt(argc, argv, "hi:r:q:s:")) != EOF))
    {
      switch(ch)
      {
        case 'h':
          Usage(argv[0], "Instructions");
          break;
        case 'r':
          scaffFilenames[Reference] = optarg;
          break;
        case 'q':
          scaffFilenames[Query] = optarg;
          break;
        case 's':
          showCoordsFilename = optarg;
          break;
        case 'i':
          minIdentity = atof(optarg);
          break;
        default:
          errflg++;
          break;
      }
    }
    if(scaffFilenames[Reference] == NULL)
      Usage(argv[0], "Please specify a reference .scaff file");
    if(scaffFilenames[Query] == NULL)
      Usage(argv[0], "Please specify a query .scaff file");
    if(showCoordsFilename == NULL)
      Usage(argv[0], "Please specify a .show-coords file");
    if(minIdentity < 1)
      Usage(argv[0], "Please specify a reasonable minIdentity");
  }

  // read in .scaff files
  map<CDS_UID_t, Scaffold> assemblies[2];
  for(int i = 0; i < NumSequenceTypes; i++)
  {
    readScaffFile(scaffFilenames[i], assemblies[i]);
#ifdef VERBOSE
    printScaffs(cout, assemblies[i]);
#endif
  }

  // read in the show-coords file
  readShowCoordsFile(showCoordsFilename, assemblies, minIdentity);
  
  for(int i = 0; i < NumSequenceTypes; i++)
  {
#ifdef VERBOSE
    cout << "Matches in "
         << ((SequenceType) i == Reference ? "reference " : "query ")
         << "sequence:\n";
    printMatches(cout, assemblies[i]);
    cout << endl;
#endif
    cout << "Mapping discrepancies from the perspective of the "
         << ((SequenceType) i == Reference ? "reference " : "query ")
         << "sequence:\n";
    checkMatches1(cout, assemblies[i], 0.03, 5);
    cout << endl;
  }
  return 0;
}
