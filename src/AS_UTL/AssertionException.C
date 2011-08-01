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

#include "AssertionException.h"

static const char* RCSID = "$Id: AssertionException.C,v 1.2 2011-08-01 16:54:03 mkotelbajcvi Exp $";

AssertionException::AssertionException(AssertionType type, const char* message, const char* file, int line, const char* function) throw() 
	: RuntimeException(message)
{
	this->type = type;
	this->file = (char*)file;
	this->line = line;
	this->function = (char*)function;
}

const char* AssertionException::what() const throw()
{
	string str(((this->file != NULL) ? this->file : UNKNOWN_LOCATION));
	
	str += ":";
	str += (this->line > 0) ? StringUtils::toString(this->line) : UNKNOWN_LOCATION;
	str += " ";
	
	if (this->function != NULL)
	{
		str += "[";
		str += this->function;
		str += "] ";
	}
	
	str += "Assert ";
	str += assertionTypeToString(this->type);
	str += " failed";
	
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

const char* AssertionException::assertionTypeToString(AssertionType type)
{
	switch (type)
	{
		case ASSERT_FALSE:
			return "false";
		case ASSERT_TRUE:
			return "true";
		case ASSERT_EQUALS:
			return "equals";
		case ASSERT_NOT_EQUALS:
			return "not equals";
		case ASSERT_NULL:
			return "null";
		case ASSERT_NOT_NULL:
			return "not null";
		case ASSERT_EMPTY:
			return "empty";
		case ASSERT_NOT_EMPTY:
			return "not empty";
		default:
			return NULL;
	}
}
