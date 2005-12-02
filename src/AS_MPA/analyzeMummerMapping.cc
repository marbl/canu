
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
/* $Id: analyzeMummerMapping.cc,v 1.2 2005-12-02 22:11:52 catmandew Exp $ */
#include <cstdio>  // for sscanf
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <unistd.h> /* man 3 getopt */

using namespace std;

#include "cds.h"

/****************************************************************************
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
    4. analyze mapping

       for each assembly
         for each scaffold
           check intra-scaffold, inter-match consistency
             to the same scaffold?, spacing, orientation
             if not consistent,
               how are ends of matches associated with gaps?
           is there unmatched (unique) sequence?

       
****************************************************************************/


/****************************************************************************/
// ENUMS
typedef enum
{
  AMM_Reference_e = 0,
  AMM_Query_e,
  AMM_NumSequenceTypes_e
} SequenceType;


typedef enum
{
  AMM_Forward_e = 0,
  AMM_Reverse_e,
  AMM_NumOrientations_e
} OrientationType;
/****************************************************************************/


/****************************************************************************/
// CLASSES
class Interval
{
public:
  Interval() {set(0,0);}
  Interval(CDS_COORD_t begin, CDS_COORD_t end)
    {
      set(begin, end);
    }

  void set(CDS_COORD_t begin, CDS_COORD_t end)
    {
      _begin = begin;
      _end = end;
      _min = (_begin < _end ? _begin : _end);
      _max = (_begin > _end ? _begin : _end);
    }
  
  void setBegin(CDS_COORD_t begin) {set(begin, getEnd());}
  void setEnd(CDS_COORD_t end) {set(getBegin(), end);}

  CDS_COORD_t getBegin() const {return _begin;}
  CDS_COORD_t getEnd() const {return _end;}
  CDS_COORD_t getMin() const {return _min;}
  CDS_COORD_t getMax() const {return _max;}
  CDS_COORD_t getLength() const {return getMax() - getMin();}
  bool isVoid() const
    {
      return (getBegin() == 0 && getEnd() == 0);
    }

  bool isForward() const {return getBegin() <= getEnd();}
  bool isReverse() const {return !isReverse();}
  OrientationType getOrientation() const
    {
      return (isForward() ? AMM_Forward_e : AMM_Reverse_e);
    }

  bool intersects(CDS_COORD_t coord) const
    {
      return (coord < getMax() && coord > getMin());
    }
  bool intersects(const Interval & interval) const
    {
      return (interval.getMin() < getMax() && interval.getMax() > getMin());
    }

  bool spans(const Interval & interval) const
    {
      return (getMin() <= interval.getMin() && getMax() >= interval.getMax());
    }

  bool isSpannedBy(const Interval & interval) const
    {
      return interval.spans(*this);
    }

  bool endsBefore(CDS_COORD_t coord) const
    {
      return (getMax() <= coord);
    }
  bool endsBeforeEnd(const Interval & interval) const
    {
      return endsBefore(interval.getMax());
    }
  bool strictlyComesBefore(const Interval & interval) const
    {
      return endsBefore(interval.getMin());
    }

  bool startsAfter(CDS_COORD_t coord) const
    {
      return (getMin() >= coord);
    }
  bool startsAfterStart(const Interval & interval) const
    {
      return startsAfter(interval.getMin());
    }
  bool strictlyComesAfter(const Interval & interval) const
    {
      return startsAfter(interval.getMax());
    }
  
private:
  CDS_COORD_t _begin;
  CDS_COORD_t _end;
  CDS_COORD_t _min;
  CDS_COORD_t _max;
};


class AssemblyInterval : public Interval
{
public:
  AssemblyInterval() {set(0,0,0);}
  AssemblyInterval(CDS_UID_t uid, CDS_COORD_t begin, CDS_COORD_t end)
    {
      set(uid, begin, end);
    }
  
  void set(CDS_UID_t uid, CDS_COORD_t begin, CDS_COORD_t end)
    {
      setUID(uid);
      Interval::set(begin, end);
    }

  void setUID(CDS_UID_t uid) {_uid = uid;}
  CDS_UID_t getUID() const {return _uid;}

  bool isVoid() const
    {
      return (getUID() == 0 && Interval::isVoid());
    }
  
private:
  CDS_UID_t _uid;
};


class Match
{
public:
  Match()
    {
      _onSeq[AMM_Reference_e].set(0,0,0);
      _onSeq[AMM_Query_e].set(0,0,0);
      _pctID = 0;
    }
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
             "%d\t%d\t%d\t%d\t%*d\t%*d\t%f\t%*d\t%*d\t%*f\t%*f\t"
             F_UID "\t" F_UID,
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
  
  const AssemblyInterval getInterval(SequenceType which) const
    {
      return _onSeq[which];
    }
  CDS_COORD_t getBegin(SequenceType which) const
    {
      return _onSeq[which].getBegin();
    }
  CDS_COORD_t getEnd(SequenceType which) const {return _onSeq[which].getEnd();}
  CDS_COORD_t getMin(SequenceType which) const {return _onSeq[which].getMin();}
  CDS_COORD_t getMax(SequenceType which) const {return _onSeq[which].getMax();}
  CDS_COORD_t getLength(SequenceType which) const
    {
      return _onSeq[which].getLength();
    }
  float getPctID() const {return _pctID;}
  OrientationType getOrientation(SequenceType which) const
    {
      return _onSeq[which].getOrientation();
    }
  bool isFlipped() const
    {
      return (getOrientation(AMM_Reference_e) != getOrientation(AMM_Query_e));
    }
  OrientationType getOrientation() const
    {
      return (isFlipped() ? AMM_Reverse_e : AMM_Forward_e);
    }
  bool intersects(SequenceType thisWhich,
                  const Match & other,
                  SequenceType otherWhich) const
    {
      return this->_onSeq[thisWhich].intersects(other._onSeq[otherWhich]);
    }

  void reversePerspective()
    {
      AssemblyInterval dummy = _onSeq[0];
      _onSeq[0] = _onSeq[1];
      _onSeq[1] = dummy;
    }
  
  void setUID(SequenceType which, CDS_UID_t uid) {_onSeq[which].setUID(uid);}
  void setBegin(SequenceType which, CDS_COORD_t begin)
    {
      _onSeq[which].setBegin(begin);
    }
  void setEnd(SequenceType which, CDS_COORD_t end) {_onSeq[which].setEnd(end);}
  void setPctID(float pctID) {_pctID = pctID;}

  friend ostream & operator<<(ostream & os, const Match & m)
    {
      for(int i = 0; i < AMM_NumSequenceTypes_e; i++)
      {
        os << m.getUID((SequenceType) i) << "\t"
           << m.getBegin((SequenceType) i) << "\t"
           << m.getEnd((SequenceType) i) << "\t";
      }
      os << m.getPctID();
      return os;
    }
  
private:
  AssemblyInterval _onSeq[2];
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
      reset();
      _uid = uid;
      _seqLength = seqLength;
      _gappedLength = gappedLength;
    }

  void reset()
    {
      _uid = _seqLength = _gappedLength = 0;
      _contigs.clear();
      _gaps.clear();
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

      AssemblyInterval newContig(uid, begin, end);
      _contigs.push_back(newContig);

      if(fabsf(gapMean) > 0.00001)
      {
        gapMean = (gapMean < 20 ? 20 : gapMean);
        Interval newGap(end, (CDS_COORD_t) (end + gapMean + 0.5));
        _gaps.push_back(newGap);
      }
    }

  CDS_UID_t getUID() const {return _uid;}
  CDS_COORD_t getSeqLength() const {return _seqLength;}
  CDS_COORD_t getGappedLength() const {return _gappedLength;}
  
  int getNumContigs() const {return _contigs.size();}
  int getNumGaps() const {return _gaps.size();}

  const AssemblyInterval getContig(int i) const {return _contigs[i];}
  const Interval getGap(int i) const {return _gaps[i];}

  bool intersectsGap(CDS_COORD_t begin, CDS_COORD_t end) const
    {
      Interval interval(begin, end);
      return intersectsGap(interval);
    }
  bool intersectsGap(const Interval & interval) const
    {
      for(int i = 0; i < getNumGaps(); i++)
      {
        if(_gaps[i].intersects(interval)) return true;
        if(_gaps[i].strictlyComesBefore(interval)) return false;
      }
    }

  bool isSpannedByGap(CDS_COORD_t begin, CDS_COORD_t end) const
    {
      Interval interval(begin, end);
      return isSpannedByGap(interval);
    }
  bool isSpannedByGap(const Interval & interval) const
    {
      for(int i = 0; i < getNumGaps(); i++)
      {
        if(_gaps[i].spans(interval)) return true;
        if(_gaps[i].startsAfterStart(interval)) return false;
      }
    }

  bool isInGap(CDS_COORD_t coord) const
    {
      for(int i = 0; i < getNumGaps(); i++)
      {
        if(_gaps[i].intersects(coord)) return true;
        if(_gaps[i].startsAfter(coord)) return false;
      }
    }

  bool isAtEndOfContig(const Interval & interval) const
    {
      return isNearEndOfContig(interval, 0);
    }

  bool isNearEndOfContig(const Interval & interval,
                         CDS_COORD_t fudgeFactor) const
    {
      for(int i = 0; i < getNumContigs(); i++)
      {
        if((_contigs[i].getMin() - fudgeFactor <= interval.getMin() &&
            _contigs[i].getMin() + fudgeFactor >= interval.getMin()) ||
           (_contigs[i].getMax() - fudgeFactor <= interval.getMax() &&
            _contigs[i].getMax() + fudgeFactor >= interval.getMax()))
          return true;
        if(_contigs[i].strictlyComesAfter(interval)) return false;
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
        const AssemblyInterval contig = getContig(i);
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
      return os;
    }
  
protected:
  CDS_UID_t _uid;
  CDS_COORD_t _seqLength;
  CDS_COORD_t _gappedLength;
  vector<AssemblyInterval> _contigs;
  vector<Interval> _gaps;
};


class MappedScaffold : public Scaffold
{
public:
  MappedScaffold() {set(0,0,0);}
  
  int set(char * line)
    {
      reset();
      return Scaffold::set(line);
    }
  void set(CDS_UID_t uid, CDS_COORD_t seqLength, CDS_COORD_t gappedLength)
    {
      reset();
      Scaffold::set(uid, seqLength, gappedLength);
    }
  
  void reset()
    {
      Scaffold::reset();
      _matches.clear();
      _mappedScaffoldUID = 0;
      _mappedBPs = _mappedScaffoldBPs = 0;
      _mappedScaffoldOrientation = AMM_Forward_e;
      _hasDeletion =
        _hasInsertion =
        _hasInversion =
        _mapsToMultipleScaffolds = false;
    }

  void addContig(char * line)
    {
      Scaffold::addContig(line);
      finalizeMapping();
    }
  
  void considerMatch(const Match & match, float minIdentity)
    {
      /*
        given the current list of matches, compare, & do one of the
        following:
          discard match,
          insert (prepend, insert, append) match, or
          replace 1 or more current matches with this one
       */
      if(match.getPctID() < minIdentity) return;
      
      if(_matches.size() == 0)
      {
        _matches.push_back(match);
        finalizeMapping();
        return;
      }

      /*
        find the first match in the current list that intersects this
        match. '<=' is appropriate because in inter-base coordinates
        the '=' still means the last bp of the minIter match precedes
        the first bp of this match
      */
      list<Match>::iterator minIter;
      for(minIter = _matches.begin();
          minIter != _matches.end() &&
            minIter->getMax(AMM_Reference_e) <= match.getMin(AMM_Reference_e);
          minIter++);

      /*
        if loop got to the end (i.e., beyond the last match) of the
        current list of matches, this match is beyond the last one, so
        append & return
      */
      if(minIter == _matches.end())
      {
        _matches.push_back(match);
        finalizeMapping();
        return;
      }

      /*
        if there is no intersection, then this match precedes the
        one pointed to by minIter and can be inserted cleanly
      */
      if(!match.intersects(AMM_Reference_e, *minIter, AMM_Reference_e))
      {
        _matches.insert(minIter, match);
        finalizeMapping();
        return;
      }

      /*
        If here, there is at least one intersecting match to compare
        with. Start tallying some values to compare later & decide
        which match(es) to keep or delete based on this match
        intersecting one or more in the current list
      */
      int numIntersecting = 1;
      CDS_COORD_t minIntersectCoord = minIter->getMin(AMM_Reference_e);
      CDS_COORD_t maxIntersectCoord = minIter->getMax(AMM_Reference_e);
      CDS_COORD_t bpCovered = minIter->getLength(AMM_Reference_e);
      float sumBPIdentities =
        minIter->getPctID() * minIter->getLength(AMM_Reference_e);

      // find first match beyond the one pointed to by minIter
      list<Match>::iterator maxIter = minIter;
      for(++maxIter;
          maxIter != _matches.end() &&
            maxIter->getMin(AMM_Reference_e) < match.getMax(AMM_Reference_e);
          maxIter++)
      {
        // for each match, update the variables being kept track of
        numIntersecting++;
        minIntersectCoord =
          (maxIter->getMin(AMM_Reference_e) < minIntersectCoord ?
           maxIter->getMin(AMM_Reference_e) : minIntersectCoord);
        maxIntersectCoord =
          (maxIter->getMax(AMM_Reference_e) > maxIntersectCoord ?
           maxIter->getMax(AMM_Reference_e) : maxIntersectCoord);
        bpCovered += maxIter->getLength(AMM_Reference_e);
        sumBPIdentities +=
          maxIter->getPctID() * maxIter->getLength(AMM_Reference_e);
      }

      /*
        Need to determine whether this match or the intersecting
        match(es) is preferable: age-old problem.

        Since we expect only high-identity matches, select the match(es)
        that cover the most basepairs

        Obviously there is code here to support other heuristics...
      */
      // float meanBPIdentity = sumBPIdentities / bpCovered;
      if(bpCovered > match.getLength(AMM_Reference_e))
      {
        return;
      }

      // Insert new match upstream of minIter
      _matches.insert(minIter, match);
      
      /*
        erase any intersecting matches including minIter up to
        but excluding maxIter
      */
      if(numIntersecting > 0)
        _matches.erase(minIter, maxIter);
      
      finalizeMapping();
    }

  CDS_UID_t getMappedScaffoldUID() const
    {
      return _mappedScaffoldUID;
    }

  CDS_COORD_t getMappedScaffoldBPs() const
    {
      return _mappedScaffoldBPs;
    }

  OrientationType getMappedScaffoldOrientation() const
    {
      return _mappedScaffoldOrientation;
    }

  CDS_COORD_t getMappedBPs() const
    {
      return _mappedBPs;
    }

  int getNumMatches() const {return _matches.size();}

  /*
    simple check routine
    are consecutive matches consistent wrt other scaffold & coords
  */
  void checkConsecutiveMatchConsistency(ostream & os,
                                        bool print,
                                        float ratio,
                                        int bps)
    {
      // assumption that this is the 'reference' sequence, relatively
      Match lastMatch;
      _hasDeletion =
        _hasInsertion =
        _hasInversion =
        _mapsToMultipleScaffolds = false;

      list<Match>::const_iterator iter;
      for(iter = _matches.begin(); iter != _matches.end(); iter++)
      {
        if(!lastMatch.isVoid())
        {
          // check if mapping continues in same scaffold
          if(lastMatch.getUID(AMM_Query_e) != iter->getUID(AMM_Query_e))
          {
            if(print)
            {
              os << "Mapping splits scaffold "
                 << lastMatch.getUID(AMM_Reference_e) << ":\n";
              os << lastMatch << endl;
              os << *iter << endl << endl;
            }
            _mapsToMultipleScaffolds = true;
          }
          else
          {
            // check for a change in orientation
            if(lastMatch.getOrientation(AMM_Query_e) !=
               iter->getOrientation(AMM_Query_e) ||
               lastMatch.getOrientation(AMM_Reference_e) !=
               iter->getOrientation(AMM_Reference_e))
            {
              // change in orientation
              if(print)
              {
                os << "Relative inversion breakpoint in scaffold "
                   << lastMatch.getUID(AMM_Reference_e) 
                   << " between "
                   << lastMatch.getMax(AMM_Reference_e) << " and "
                   << iter->getMin(AMM_Reference_e) << ":\n";
                os << lastMatch << endl;
                os << *iter << endl << endl;
              }
              _hasInversion = true;
            }
            else
            {
              // check if spacing is about right
              CDS_COORD_t refDelta = iter->getMin(AMM_Reference_e) -
                lastMatch.getMax(AMM_Reference_e);
              CDS_COORD_t queryDelta =
                (lastMatch.isFlipped() ?
                 lastMatch.getMin(AMM_Query_e) - iter->getMax(AMM_Query_e) :
                 iter->getMin(AMM_Query_e) - lastMatch.getMax(AMM_Query_e));
              CDS_COORD_t diff = refDelta - queryDelta;
              float diffRatio = (queryDelta == 0 ? 1. : diff / queryDelta);
              
              if(diffRatio > ratio || diff > bps)
              {
                if(print)
                {
                  os << "Relative insertion in scaffold "
                     << lastMatch.getUID(AMM_Reference_e)
                     << " of " << diff << " bps between "
                     << lastMatch.getMax(AMM_Reference_e) << " and "
                     << iter->getMin(AMM_Reference_e) << ":\n";
                  os << lastMatch << endl;
                  os << *iter << endl << endl;
                }
                _hasInsertion = true;
              }
              else if(diffRatio < -ratio || diff < -bps)
              {
                if(print)
                {
                  os << "Relative deletion in scaffold "
                     << lastMatch.getUID(AMM_Reference_e)
                     << " of " << -diff << " bps between "
                     << lastMatch.getMax(AMM_Reference_e) << " and "
                     << iter->getMin(AMM_Reference_e) << ":\n";
                  os << lastMatch << endl;
                  os << *iter << endl << endl;
                }
                _hasDeletion = true;
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

  void printBestScaffoldMapping(ostream & os) const
    {
      float pct = 0;
      if(getSeqLength() > 0)
        pct = 100. * getMappedScaffoldBPs() / (float) getSeqLength();
      
      os << getUID() << "\t"
         << getMappedBPs() << "\t"
         << getSeqLength() << "\t"
         << getGappedLength() << "\t"
         << getMappedScaffoldUID() << "\t"
         << getMappedScaffoldBPs() << "\t"
         << pct << "\t"
         << (getMappedScaffoldOrientation() == AMM_Forward_e ? "f" : "r")
         << endl;
    }

  void printMisMappings(ostream & os) const
    {
      if(!_hasDeletion &&
         !_hasInsertion &&
         !_hasInversion &&
         !_mapsToMultipleScaffolds &&
         getSeqLength() == _mappedBPs)
        return;

      /*
        loop over contigs, since there has to be at least one, and
        check matches to find mis-mappings
      */
      CDS_COORD_t currBP;
      list<Match>::const_iterator iter;
      int c;
      for(currBP = 0, c = 0, iter = _matches.begin(); c < getNumContigs(); c++)
      {
        // loop over matches
        while(iter != _matches.end() &&
              iter->getMax(AMM_Reference_e) <= getContig(c).getMax())
        {
          iter++;
        }
      }
      
      os << getUID() << " has mapping issues\n";
    }
  
private:
  list<Match> _matches;
  CDS_UID_t _mappedScaffoldUID;
  CDS_COORD_t _mappedScaffoldBPs;
  OrientationType _mappedScaffoldOrientation;
  CDS_COORD_t _mappedBPs;
  bool _hasDeletion;
  bool _hasInsertion;
  bool _hasInversion;
  bool _mapsToMultipleScaffolds;
  
  void finalizeMapping()
    {
      if(getNumMatches() == 0)
        return;
      
      // iterate over matches, count basepairs matched to each scaffold
      // scaffold with most matching basepairs is the mapped scaffold
      map<CDS_UID_t, CDS_COORD_t> matchBPs[2];
      _mappedBPs = _mappedScaffoldBPs = 0;

      list<Match>::iterator iter;
      for(iter = _matches.begin(); iter != _matches.end(); iter++)
      {
        CDS_UID_t otherUID = iter->getUID(AMM_Query_e);
        int orient = (int) iter->getOrientation();
        
        _mappedBPs += iter->getLength(AMM_Reference_e);

        matchBPs[orient][otherUID] += iter->getLength(AMM_Reference_e);
        matchBPs[1 - orient][otherUID] += 0;
        
        if(matchBPs[AMM_Forward_e][otherUID] +
           matchBPs[AMM_Reverse_e][otherUID] > _mappedScaffoldBPs)
        {
          _mappedScaffoldUID = otherUID;
          _mappedScaffoldBPs = matchBPs[AMM_Forward_e][otherUID] +
            matchBPs[AMM_Reverse_e][otherUID];
          _mappedScaffoldOrientation = (matchBPs[AMM_Forward_e][otherUID] >=
                                        matchBPs[AMM_Reverse_e][otherUID] ?
                                        AMM_Forward_e : AMM_Reverse_e);
        }
      }

      checkConsecutiveMatchConsistency(cout, false, 0.03, 5);
      if(!_mapsToMultipleScaffolds && _mappedScaffoldBPs != _mappedBPs)
      {
        assert(0);
      }
    }
};


// An assembly is simply map<UID, Scaffold>, so no class needed for that
// An assembly pair is simply an array of two assemblies...
// A mapped assembly pair is simply an array of two map<UID, MappedScaffold>
/****************************************************************************/


/****************************************************************************/
// FUNCTIONS
void printBestScaffoldMappings(ostream & os,
                               const map<CDS_UID_t, MappedScaffold> & assembly)
{
  os << "ThisScaffoldUID\t"
     << "TotalBasesMapped\t"
     << "SequenceLength\t"
     << "TotalLength\t"
     << "MappedToScaffoldUID\t"
     << "BasesMappedToScaffold\t"
     << "PctBasesMappedToScaffold\t"
     << "OrientationOfMapping\n";
  
  map<CDS_UID_t, MappedScaffold>::const_iterator iter;
  for(iter = assembly.begin(); iter != assembly.end(); iter++)
  {
    ((MappedScaffold) iter->second).printBestScaffoldMapping(os);
  }
}


void printMisMappings(ostream & os,
                      const map<CDS_UID_t, MappedScaffold> & assembly)
{
  map<CDS_UID_t, MappedScaffold>::const_iterator iter;
  for(iter = assembly.begin(); iter != assembly.end(); iter++)
  {
    ((MappedScaffold) iter->second).printMisMappings(os);
  }
}


void printMatches(ostream & os,
                  const map<CDS_UID_t, MappedScaffold> & assembly)
{
  map<CDS_UID_t, MappedScaffold>::const_iterator iter;
  for(iter = assembly.begin(); iter != assembly.end(); iter++)
  {
    ((MappedScaffold) iter->second).printMatches(os);
  }
}

void printScaffs(ostream & os, const map<CDS_UID_t, MappedScaffold> & assembly)
{
  map<CDS_UID_t, MappedScaffold>::const_iterator iter;
  for(iter = assembly.begin(); iter != assembly.end(); iter++)
  {
    ((MappedScaffold) iter->second).printStructure(os);
  }
}

void readScaffFile(char * filename, map<CDS_UID_t, MappedScaffold> & assembly)
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
    MappedScaffold scaff;
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
                        map<CDS_UID_t, MappedScaffold> assemblies[2],
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
    assemblies[AMM_Reference_e][match.getUID(AMM_Reference_e)].considerMatch(match,
                                                                             minIdentity);
    match.reversePerspective();
    assemblies[AMM_Query_e][match.getUID(AMM_Reference_e)].considerMatch(match,
                                                             minIdentity);
  }
  fin.close();
}

void Usage(char * progname, char * message)
{
  if(message != NULL)
    cerr << endl << message << endl << endl;;
  cerr << "Usage: " << progname << " [-h] [-v level] [-i minIdentity] -r ref.scaff  -q query.scaff  -s showCoordsFile\n"
       << "\t-h                 print usage and exit\n"
       << "\t-v level           set verbosity level\n"
       << "\t-r ref.scaff       .scaff file of reference sequence\n"
       << "\t-q query.scaff     .scaff file of query sequence\n"
       << "\t-s showCoordsFile  output file from show-coords for 2 sequences\n"
       << "\t-i minIdentity     minimum pct identity of matches to consider\n"
       << "\t                     default is 95\n"
       << "\nOutput written to stdout\n"
       << "\n\n";
  exit(1);
}


int main(int argc, char ** argv)
{
  char * scaffFilenames[2];
  char * showCoordsFilename = NULL;
  float minIdentity = 95;
  scaffFilenames[AMM_Reference_e] = scaffFilenames[AMM_Query_e] = NULL;
  int verboseLevel = 0;

  // process the command line
  {
    int ch, errflg = 0;
    while(!errflg && ((ch = getopt(argc, argv, "hv:i:r:q:s:")) != EOF))
    {
      switch(ch)
      {
        case 'h':
          Usage(argv[0], "Instructions");
          break;
        case 'v':
          verboseLevel = atoi(optarg);
          break;
        case 'r':
          scaffFilenames[AMM_Reference_e] = optarg;
          break;
        case 'q':
          scaffFilenames[AMM_Query_e] = optarg;
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
    if(scaffFilenames[AMM_Reference_e] == NULL)
      Usage(argv[0], "Please specify a reference .scaff file");
    if(scaffFilenames[AMM_Query_e] == NULL)
      Usage(argv[0], "Please specify a query .scaff file");
    if(showCoordsFilename == NULL)
      Usage(argv[0], "Please specify a .show-coords file");
    if(minIdentity < 1)
      Usage(argv[0], "Please specify a reasonable minIdentity");
  }

  // read in .scaff files
  map<CDS_UID_t, MappedScaffold> assemblies[2];
  for(int i = 0; i < AMM_NumSequenceTypes_e; i++)
  {
    readScaffFile(scaffFilenames[i], assemblies[i]);
    if(verboseLevel > 1)
    {
      printScaffs(cout, assemblies[i]);
    }
  }

  // read in the show-coords file
  readShowCoordsFile(showCoordsFilename, assemblies, minIdentity);
  
  for(int i = 0; i < AMM_NumSequenceTypes_e; i++)
  {
    if(verboseLevel > 0)
    {
      cout << "Matches in "
           << ((SequenceType) i == AMM_Reference_e ? "reference " : "query ")
           << "sequence:\n";
      printMatches(cout, assemblies[i]);
      cout << endl;
    }
    
    cout << "Mapping inconsistencies in the "
         << ((SequenceType) i == AMM_Reference_e ? "reference " : "query ")
         << "sequence:\n";
    printMisMappings(cout, assemblies[i]);
    cout << endl;
                     
    /*
    cout << "Best scaffold mapping for each scaffold in the "
         << ((SequenceType) i == AMM_Reference_e ? "reference " : "query ")
         << "sequence:\n";
    printBestScaffoldMappings(cout, assemblies[i]);
    cout << endl;

    cout << "Gaps filled in the "
         << ((SequenceType) i == AMM_Reference_e ? "reference " : "query ")
         << "sequence:\n";
    printGapsFilled(cout, assemblies[i]);

    cout << "Unique mid-contig sequence in the "
         << ((SequenceType) i == AMM_Reference_e ? "reference " : "query ")
         << "sequence:\n";
    printUniqueMidContigSequences(cout, assemblies[i]);
    */
  }
  return 0;
}
/****************************************************************************/
// FIN
