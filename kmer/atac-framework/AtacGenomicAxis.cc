// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Ross Lippert
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include <assert.h>
#include "AtacGenomicAxis.h"
#include "ctype.h"

// Create an invalid instance
AtacGenomicAxis::AtacGenomicAxis()
{
}

// Create a valid instance
AtacGenomicAxis::AtacGenomicAxis(string typeId, string id, string parentId, int64 parentStart, string sequenceSourceId, int64 start, int64 length, int64 offset, int ordinal, string name, string uid)
{
  _typeId = typeId;
  _id = id;
  _parentId = parentId;
  _parentStart = parentStart;
  _sequenceSourceId = sequenceSourceId;
  _start = start;
  _length = length;
  _offset = offset;
  _ordinal = ordinal;
  _name = name;
  _uid = uid;
}

void AtacGenomicAxis::print(ostream& os) const
{
  os << "G " << " " << emptyToDot(_typeId) << " " << _id << " " << emptyToDot(_parentId) << " " << _parentStart << " " << emptyToDot(_sequenceSourceId) << " " << _start << " " << _length << " " << _offset << " " << _ordinal << " " << emptyToDot(_name) << " " << emptyToDot(_uid);

  if (_defline.size() > 0)
    os << ">" << _defline;

  os << "\n";
} 

char AtacGenomicAxis::getClassCode() const
{
  return ('G');
}

void AtacGenomicAxis::write(FILE *out) const
{
  string typeId = emptyToDot(_typeId).c_str();
  string parentId = emptyToDot(_parentId).c_str();
  string sequenceSourceId = emptyToDot(_sequenceSourceId).c_str();
  string name = emptyToDot(_name).c_str();
  string uid = emptyToDot(_uid).c_str();

  fprintf(out, "G %s %s %s %lld %s %lld %lld %lld %d %s %s", 
            typeId.c_str(), _id.c_str(), parentId.c_str(), 
            _parentStart, sequenceSourceId.c_str(), _start,
            _length, _offset, _ordinal, name.c_str(), uid.c_str());

  if (_defline.size() > 0)
    fprintf(out, " >%s\n", _defline.c_str());
  else
    fprintf(out, "\n");
}

string AtacGenomicAxis::getDefline() const { return _defline; }
void AtacGenomicAxis::setDefline(string defline) { _defline = defline; }

int AtacGenomicAxis::getRowIndex() const { return _rowIndex; }
void AtacGenomicAxis::setRowIndex(int rowIndex) { _rowIndex = rowIndex; }

string AtacGenomicAxis::getTypeId() const { return _typeId; }
string AtacGenomicAxis::getId() const { return _id; }
string AtacGenomicAxis::getParentId() const { return _parentId; }
int64 AtacGenomicAxis::getParentStart() const  { return _parentStart; }

string AtacGenomicAxis::getSequenceSourceId() const  { return _sequenceSourceId; }
int64 AtacGenomicAxis::getStart() const  { return _start; }
int64 AtacGenomicAxis::getLength() const  { return _length; }
int64 AtacGenomicAxis::getOffset() const { return _offset; }
int AtacGenomicAxis::getOrdinal() const { return _ordinal; }
