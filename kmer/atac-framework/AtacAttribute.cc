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
#include "AtacAttribute.h"
#include "ctype.h"

AtacAttribute::AtacAttribute()
{
}

// Create a valid instance
AtacAttribute::AtacAttribute(string key, string value) :
  _key(key),
  _value(value)
{
}

AtacAttribute::AtacAttribute(string key, int value) :
  _key(key)
{
   char valueBuffer[32];
   sprintf(valueBuffer, "%d", value);
   _value = valueBuffer;
}

AtacAttribute::AtacAttribute(string key, int64 value) :
  _key(key)
{
   char valueBuffer[64];
   sprintf(valueBuffer, "%lld", value);
   _value = valueBuffer;
}

AtacAttribute::AtacAttribute(string key, double value) :
  _key(key)
{
   char valueBuffer[32];
   sprintf(valueBuffer, "%g", value);
   _value = valueBuffer;
}

void AtacAttribute::print(ostream& os) const
{
  os << "/" << _key << "=" << _value << "\n";
} 

char AtacAttribute::getClassCode() const
{
  return ('/');
}

void AtacAttribute::write(FILE *out) const
{
  fprintf(out, "/%s=%s\n", _key.c_str(), _value.c_str());
}

string AtacAttribute::getKey() const { return _key; }
string AtacAttribute::getValue() const { return _value; }
int AtacAttribute::getIntValue() const { return atoi(_value.c_str()); }
// int64 AtacAttribute::getInt64Value() const { return atoll(_value.c_str()); }
int64 AtacAttribute::getInt64Value() const { return atol(_value.c_str()); }
double AtacAttribute::getDoubleValue() const
{
  double result;
  sscanf(_value.c_str(), "%g", &result);
  return result;
}
