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
#include "AtacComment.h"
#include "ctype.h"

AtacComment::AtacComment()
{
}

// Create a valid instance
AtacComment::AtacComment(string text) :
  _text(text)
{
}

void AtacComment::print(ostream& os) const
{
  os << "# " << _text;
} 

char AtacComment::getClassCode() const
{
  return ('#');
}

void AtacComment::write(FILE *out) const
{
  if(0 != _text.c_str()) {
    fprintf(out, "# %s\n", _text.c_str());
  }
}

string AtacComment::getText() const { return _text; }
