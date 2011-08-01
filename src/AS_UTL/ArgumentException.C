/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005, J. Craig Venter Institute. All rights reserved.
 * Author: Brian Walenz
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
 *************************************************************************/

#include "ArgumentException.h"

static const char* RCSID = "$Id: ArgumentException.C,v 1.3 2011-08-01 16:54:03 mkotelbajcvi Exp $";

ArgumentException::ArgumentException(const char* name, const char* message) throw() 
	: RuntimeException(message)
{
	this->name = (char*)name;
}

const char* ArgumentException::what() const throw()
{
	string str("Argument ");
	
	if (this->name != NULL)
	{
		str += "(name=";
		str += this->name;
		str += ")";
	}
	
	str += "exception";
	
	if (this->message != NULL)
	{
		str += ": ";
		str += this->message;
	}
	else
	{
		str += ".";
	}
	
	return str.c_str();
}
