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


#ifndef _INCLUDE_ATACATTRIBUTE
#define _INCLUDE_ATACATTRIBUTE

#include <iostream>
#include <string>
#include "AtacRecord.h"

#ifndef _USING_NAMESPACE_STD_
#define _USING_NAMESPACE_STD_
using namespace std;
#endif

/**
 * AtacAttribute
 **/
class AtacAttribute : public AtacRecord {
 public:
  AtacAttribute(); 
  AtacAttribute(string key, string value);
  AtacAttribute(string key, int value);
  AtacAttribute(string key, int64 value);
  AtacAttribute(string key, double value);

  // overridden
  char getClassCode() const;
  void print(ostream& os) const;
  void write(FILE *out) const;

  // newly defined
  string getKey() const;
  string getValue() const;
  int getIntValue() const;
  int64 getInt64Value() const;
  double getDoubleValue() const;

 private:
  string _key;
  string _value;

};

inline ostream &operator<< (ostream&os, const AtacAttribute& t) {
  t.print(os); return (os); 
}


#endif
