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


#ifndef _INCLUDE_ATACCOMMENT
#define _INCLUDE_ATACCOMMENT

#include <stdio.h>
#include <iostream>
#include <string>
#include "AtacRecord.h"

#ifndef _USING_NAMESPACE_STD_
#define _USING_NAMESPACE_STD_
using namespace std;
#endif

/**
 * AtacComment
 **/
class AtacComment : public AtacRecord {
 public:
  AtacComment(); 
  AtacComment(string text);

  // overridden
  char getClassCode() const;
  void print(ostream& os) const;
  void write(FILE *out) const;

  // newly defined
  string getText() const;

 private:
  string _text;

};

inline ostream &operator<< (ostream&os, const AtacComment& t) {
  t.print(os); return (os); 
}


#endif
