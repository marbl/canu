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


#ifndef _INCLUDE_ATACRECORD
#define _INCLUDE_ATACRECORD

#include <stdio.h>
#include <iostream>
#include <string>
#include "AtacTypes.h"

#ifndef _USING_NAMESPACE_STD_
#define _USING_NAMESPACE_STD_
using namespace std;
#endif

/**
 * AtacRecord
 **/
class AtacRecord {
 public:
  AtacRecord();
  virtual ~AtacRecord();

  virtual char getClassCode() const;
  virtual void invalidate();
  virtual bool isValid() const;
  virtual void setImplicit(bool state);
  virtual bool isImplicit() const;
  virtual void print(ostream& os) const;
  virtual void write(FILE *out) const;

 protected:
  /* Check id string and convert "." to empty strings */
  const string dotToEmpty(string id) const;

  /* Check id string and convert to empty strings to "." */
  const string emptyToDot(string id) const;

 private:
  bool _valid;
  bool _implicit;

};

inline ostream &operator<< (ostream& os, const AtacRecord& t)
{
  t.print(os); return (os); 
}

#endif
