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
#include "AtacMatch.h"
#include "ctype.h"

// Create an invalid instance
AtacMatch::AtacMatch()
{
}

// Create a valid instance
AtacMatch::AtacMatch(string typeId, string id, string parentId, string axisId1, int64 start1, int64 length1, int orient1, string axisId2, int64 start2, int64 length2, int orient2)
{
  _typeId = typeId;
  _id = id;

  _parentId = parentId;
  if (_parentId == ".") _parentId = "";

  _axisId[0] = axisId1;
  _start[0] = start1;
  _length[0] = length1;
  _orient[0] = orient1;

  _axisId[1] = axisId2;
  _start[1] = start2;
  _length[1] = length2;
  _orient[1] = orient2;
}

void AtacMatch::print(ostream& os) const
{
  os << "M " << " " << emptyToDot(_typeId) << " " 
     << _id << " " << emptyToDot(_parentId) << " " 
     << emptyToDot(_axisId[0]) << " " 
     << _start[0] << " " << _length[0] << " " 
     << _orient[0] << " " << emptyToDot(_axisId[1]) << " " 
     << _start[1] << " " << _length[1] << " " << _orient[1];

  if (_defline.size() > 0)
    os << ">" << _defline;

  os << "\n";
} 

char AtacMatch::getClassCode() const
{
  return ('M');
}


void AtacMatch::write(FILE *out) const
{
    string typeId = emptyToDot(_typeId);
    string parentId = emptyToDot(_parentId);
    string axisId1 = emptyToDot(_axisId[0]);
    string axisId2 = emptyToDot(_axisId[1]);

    fprintf(out, "M %s %s %s %s %lld %lld %d %s %lld %lld %d", typeId.c_str(), _id.c_str(), parentId.c_str(), 
            axisId1.c_str(), _start[0], _length[0], _orient[0], 
            axisId2.c_str(), _start[1], _length[1], _orient[1]);

  if (_defline.size() > 0)
    fprintf(out, " >%s\n", _defline.c_str());
  else
    fprintf(out, "\n");
}

string AtacMatch::getDefline() const { return _defline; }
void AtacMatch::setDefline(string defline) { _defline = defline; }

int AtacMatch::getRowIndex() const { return _rowIndex; }
void AtacMatch::setRowIndex(int rowIndex) { _rowIndex = rowIndex; }

string AtacMatch::getTypeId() const { return _typeId; }
string AtacMatch::getId() const { return _id; }
string AtacMatch::getParentId() const { return _parentId; }

void AtacMatch::addChild(string childId) { _children.push_back(childId); }
vector<string> *AtacMatch::getChildren() { return &_children; }
bool AtacMatch::hasChildren() { return _children.size() > 0; }

string AtacMatch::getGenomicAxisId(int which) const  { return _axisId[which]; }
int64 AtacMatch::getStart(int which) const  { return _start[which]; }
int64 AtacMatch::getLength(int which) const  { return _length[which]; }
int AtacMatch::getOrientation(int which) const { return _orient[which]; }
