/**************************************************************************
 * This file is part of A2Amapper.
 * Copyright (c) 2004 Applera Corporation
 * Author: Ross Lippert
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
**************************************************************************/


#ifndef _INCLUDE_ATACMATCH
#define _INCLUDE_ATACMATCH

#include <iostream>
#include <string>
#include <vector>
#include "AtacRecord.h"

#ifndef _USING_NAMESPACE_STD_
#define _USING_NAMESPACE_STD_
using namespace std;
#endif

/**
 * AtacMatch
 **/
class AtacMatch : public AtacRecord {
 public:
  AtacMatch(); 
  AtacMatch(string typeId, string id, string parentId, string axisId1, int64 start1, int64 length1, int orient1, string axisId2, int64 start2, int64 length2, int orient2);

  // overridden
  char getClassCode() const;
  void print(ostream& os) const;
  void write(FILE *out) const;

  // newly defined
  int getRowIndex() const;
  void setRowIndex(int rowIndex);

  string getTypeId() const;
  string getId() const;
  string getParentId() const;
  
  void addChild(string childId);
  vector<string> *getChildren();
  bool hasChildren();

  string getGenomicAxisId(int which) const ;
  int64 getStart(int which) const ;
  int64 getLength(int which) const ;
  int getOrientation(int which) const ;

  void setDefline(string line);
  string getDefline() const;

 private:
  int _rowIndex;
  string _typeId;
  string _id;
  string _parentId;
  vector<string> _children;

  string _axisId[2];
  int64 _start[2];
  int64 _length[2];
  int _orient[2];

  string _defline;

};

inline ostream &operator<< (ostream&os, const AtacMatch& t) {
  t.print(os); return (os); 
}

#endif
