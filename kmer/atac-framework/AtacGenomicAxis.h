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


#ifndef _INCLUDE_ATACGENOMICAXIS
#define _INCLUDE_ATACGENOMICAXIS

#include <iostream>
#include <string>
#include <vector>
#include "AtacRecord.h"

#ifndef _USING_NAMESPACE_STD_
#define _USING_NAMESPACE_STD_
using namespace std;
#endif

/**
 * AtacGenomicAxis
 **/
class AtacGenomicAxis : public AtacRecord {
 public:
  AtacGenomicAxis(); 
  AtacGenomicAxis(string typeId, string id, string parentId, int64 parentStart, string sequenceId, int64 start, int64 length, int64 offset, int ordinal, string name, string uid);

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
  int64 getParentStart() const ;
  string getSequenceSourceId() const;
  
  int64 getStart() const ;
  int64 getLength() const ;
  int64 getOffset() const ;
  int getOrdinal() const ;

  string getName() const;
  string getUid() const;

  void setDefline(string line);
  string getDefline() const;

 private:
  int _rowIndex;
  string _typeId;
  string _id;
  string _parentId;
  int64 _parentStart;
  string _sequenceSourceId;

  int64 _start;
  int64 _length;
  int64 _offset;
  int _ordinal;

  string _name;
  string _uid;

  string _defline;

};

inline ostream &operator<< (ostream&os, const AtacGenomicAxis& t) {
  t.print(os); return (os); 
}

#endif
