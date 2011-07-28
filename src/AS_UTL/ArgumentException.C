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

static const char* rcsid_ARGUMENT_EXCEPTION_C = "$Id: ArgumentException.C,v 1.1 2011-07-28 04:23:50 mkotelbajcvi Exp $";

ArgumentException::ArgumentException(char* name, char* message) throw()
{
	this->name = name;
	this->message = (char*)((message != NULL) ? message : "");
}

const char* ArgumentException::what() const throw()
{
	char* buffer = new char[1024];
	
	if (this->name != NULL)
	{
		sprintf(buffer, "Argument exception: %s\n", this->message);
	}
	else
	{
		sprintf(buffer, "Argument (name=%s) exception: %s\n", this->name, this->message);
	}
	
	return (const char*)buffer;
}

ArgumentException::operator const char*()
{
	return this->what();
}
