
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
/* $Id: getIntervalIntersections.cc,v 1.1.1.1 2004-04-14 13:52:01 catmandew Exp $ */
#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <string>

using namespace std;

class Interval
{
public:
  Interval(){set(0,0);}
  Interval(unsigned int b, unsigned int e) {set(b,e);}
  Interval(const char * line) {set(line);}

  void setBegin(unsigned int b) {_begin = b;}
  unsigned int getBegin() const {return _begin;}
  void setEnd(unsigned int e) {_end = e;}
  unsigned int getEnd() const {return _end;}
  void set(const char * line)
  {
    sscanf(line, "%u %u", &_begin, &_end);
    _end += _begin;
    strncpy(_line, line, _LineLength-1);
    _line[_LineLength-1] = '\0';
    adjustLine();
  }
  const char * getString() const {return _line;}
  void set(unsigned int b, unsigned int e)
    {
      setBegin(b);
      setEnd(e);
      
      sprintf(_line, "%u,%u", b, e);
    }
  void setFromStartAndLength(unsigned int s, unsigned int l)
    {
      set(s, s+l);
    }
  
  bool intersects(unsigned int b, unsigned int e) const
    {
      return (getBegin() < e && b < getEnd());
    }

  bool intersects(const Interval & i) const
    {
      return intersects(i.getBegin(), i.getEnd());
    }

  bool spans(unsigned int b, unsigned int e) const
    {
      return (getBegin() <= b && e <= getEnd());
    }

  bool spans(const Interval & i) const
    {
      return spans(i.getBegin(), i.getEnd());
    }

  bool operator<(const Interval & other) const
    {
      return(getBegin() < other.getBegin());
    }
  bool startsBefore(const Interval & other) const
    {
      return(getBegin() < other.getBegin());
    }
  bool startsBeyond(const Interval & other) const
    {
      return(getBegin() > other.getEnd());
    }

  friend ostream & operator<<(ostream & os, const Interval & iv)
    {
      os << iv.getString();
/*
      os << iv.getBegin() << "," << iv.getEnd()
         << "(" << iv.getEnd() - iv.getBegin() << ")";
*/
      return os;
    }
private:
  static const int _LineLength=128;
  void adjustLine()
    {
      int j = strlen(_line);
      for(int i = 0; i < j; i++)
        _line[i] = isspace(_line[i]) ? ',' : _line[i];
    }
  unsigned int _begin;
  unsigned int _end;
  char _line[_LineLength];
};


#define MLL  4096

void ReadIntervals(char * filename, list<Interval> & ilist, bool keepLines)
{
  unsigned int start, length;
  char line[MLL];
  Interval iv;

  ilist.clear();
  
  ifstream fi(filename);
  if(!fi.good())
  {
    cerr << "Failed to open file " << filename << endl;
    return;
  }
  while(fi.getline(line, MLL-1))
  {
    if(keepLines)
      iv.set(line);
    else
    {
      sscanf(line, "%u %u", &start, &length);
      iv.setFromStartAndLength(start, length);
    }
    ilist.push_back(iv);
  }
  fi.close();
}

void Usage(char * progname, char * message)
{
  if(message != NULL)
    cerr << message << endl;
  
  cerr << "Usage: " << progname << " <file1> <file2> [-s]\n";
  cerr << "  -s requires file2 intervals to span file1 intervals\n";
  cerr << "  -i requires file2 intervals to not span file1 intervals\n";
  cerr << "  -q requires file2 intervals to not be spanned by file1 intervals\n";
  cerr << "\tfields per line in files: start len\n";
  exit(-1);
}

int main(int argc, char ** argv)
{
  int chr;
  bool mustSpan = false;
  bool mustNotSpan = false;
  bool mustNotBeSpannedBy = false;

  if(argc < 3)
  {
    Usage(argv[0], "Missing parameter(s)");
  }
  for(int i = 3; i < argc; i++)
  {
    if(strncmp(argv[i],"-s",2) == 0) mustSpan = true;
    if(strncmp(argv[i],"-i",2) == 0) mustNotSpan = true;
    if(strncmp(argv[i],"-q",2) == 0) mustNotBeSpannedBy = true;
  }

  if(mustSpan && mustNotSpan)
    Usage(argv[0], "Incompatible flags!");

  list<Interval> ilist;
  list<Interval>::iterator iiter;
  ReadIntervals(argv[1], ilist, true);
  if(ilist.size() == 0)
  {
    cerr << "file 1 has no intervals\n";
    return 0;
  }
  cerr << ilist.size() << " intervals in " << argv[1] << endl;
  ilist.sort();

  list<Interval> nlist;
  list<Interval>::iterator niter;
  ReadIntervals(argv[2], nlist, false);
  if(nlist.size() == 0)
  {
    cerr << "file 2 has no intervals\n";
    // regurgitate lines of file 1 & exit
    for(iiter = ilist.begin(); iiter != ilist.end(); iiter++)
      cout << *iiter << ":\n";
    return 0;
  }
  cerr << nlist.size() << " intervals in " << argv[2] << endl;
  nlist.sort();

  // now have ilist & nlist
  list<Interval>::iterator liter;
  for(liter = niter = nlist.begin(), iiter = ilist.begin();
      iiter != ilist.end();
      iiter++)
  {
    cout << *iiter << ": ";
    while(liter != nlist.end() && iiter->startsBeyond(*liter)) liter++;
    niter = liter;
    while(niter != nlist.end() && !niter->startsBeyond(*iiter))
    {
      if(niter->intersects(*iiter) &&
         (mustSpan == false || niter->spans(*iiter)) &&
         (mustNotSpan == false || !niter->spans(*iiter)) &&
         (mustNotBeSpannedBy == false || !iiter->spans(*niter)))
        cout << " " << *niter;
      niter++;
    }
    cout << endl;
  }
  return 0;
}
