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
#include "AtacSequenceSource.h"
#include "ctype.h"

// Create an invalid instance
AtacSequenceSource::AtacSequenceSource()
{
}

// Create a valid instance
AtacSequenceSource:: AtacSequenceSource(string id, string fileName, string fileType, string elemType, int64 width, int trim)
{
  _id = id;
  _fileName = fileName;
  int dotIndex = _fileName.find(".", 0);
  _fileNamePrefix = dotIndex < 0 ? fileName : fileName.substr(0, dotIndex);
  _fileType = fileType;
  _elemType = elemType;
  _width = width;
  _trim = trim;
}

void AtacSequenceSource::print(ostream& os) const
{
  os << "S " << " " << _id << " " << _fileName << " " << _fileType << " " << _elemType << " " << _width << " " << _trim << "\n";

  if (_defline.size() > 0)
    os << ">" << _defline;

  os << "\n";
} 

char AtacSequenceSource::getClassCode() const
{
  return ('S');
}

void AtacSequenceSource::write(FILE *out) const
{
    fprintf(out, "S %s %s %s %s %lld %d", _id.c_str(), _fileName.c_str(), _fileType.c_str(), _elemType.c_str(), _width, _trim); 

  if (_defline.size() > 0)
    fprintf(out, " >%s\n", _defline.c_str());
  else
    fprintf(out, "\n");
}

string AtacSequenceSource::getDefline() const { return _defline; }
void AtacSequenceSource::setDefline(string defline) { _defline = defline; }

int AtacSequenceSource::getRowIndex() const { return _rowIndex; }
void AtacSequenceSource::setRowIndex(int rowIndex) { _rowIndex = rowIndex; }

string AtacSequenceSource::getId() const { return _id; }
string AtacSequenceSource::getFileName() const { return _fileName; }
string AtacSequenceSource::getFileNamePrefix() const { return _fileName; }
string AtacSequenceSource::getFileType() const { return _fileType; }
string AtacSequenceSource::getElementType() const { return _elemType; }

int64 AtacSequenceSource::getWidth() const  { return _width; }
int AtacSequenceSource::getTrim() const  { return _trim; }
