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
#include "AtacRecord.h"
#include "ctype.h"

// Create an invalid instance
AtacRecord::AtacRecord()
{
  _valid = true;
}

AtacRecord::~AtacRecord()
{
}

void AtacRecord::print(ostream& os) const
{
  os << getClassCode() << "\n";
} 

char AtacRecord::getClassCode() const
{
  return ('?');
}

bool AtacRecord::isValid() const
{
  return (_valid);
}

void AtacRecord::invalidate()
{
  _valid = false;
}

bool AtacRecord::isImplicit() const
{
  return (_implicit);
}

void AtacRecord::setImplicit(bool state)
{
  _implicit = state;
}

void AtacRecord::write(FILE *out) const
{
    fprintf(out, "?\n");
}

const string AtacRecord::dotToEmpty(string id) const
{
  if (id == ".")
    return "";
  else
    return id;
}

const string AtacRecord::emptyToDot(string id) const
{
  if (id == "")
    return ".";
  else
    return id;
}

