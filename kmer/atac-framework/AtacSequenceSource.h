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


#ifndef _INCLUDE_ATACSEQUENCESOURCE
#define _INCLUDE_ATACSEQUENCESOURCE

#include <iostream>
#include <string>
#include <vector>
#include "AtacRecord.h"

#ifndef _USING_NAMESPACE_STD_
#define _USING_NAMESPACE_STD_
using namespace std;
#endif

/**
 * AtacSequenceSource
 **/
class AtacSequenceSource : public AtacRecord {
 public:
  AtacSequenceSource(); 
  AtacSequenceSource(string id, string fileName, string fileType, string elemType, int64 width, int trim);

  // overridden
  char getClassCode() const;
  void print(ostream& os) const;
  void write(FILE *out) const;

  // newly defined
  int getRowIndex() const;
  void setRowIndex(int rowIndex);

  string getId() const;
  
  string getFileName() const ;
  string getFileNamePrefix() const ;
  string getFileType() const ;
  string getElementType() const ;
  int64 getWidth() const ;
  int getTrim() const ;

  void setDefline(string line);
  string getDefline() const;

 private:
  int _rowIndex;
  string _id;
  string _fileName;
  string _fileNamePrefix;
  string _fileType;
  string _elemType;
  int64 _width;
  int64 _trim;

  string _defline;

};

inline ostream &operator<< (ostream&os, const AtacSequenceSource& t) {
  t.print(os); return (os); 
}

#endif
