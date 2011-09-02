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

static const char* rcsid = "$Id: RuntimeException.C,v 1.7 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "RuntimeException.h"

RuntimeException::RuntimeException(string message, RuntimeException* cause) throw()
{
	this->message = message;
	this->cause = cause;
	
	this->stackTrace = new StackTrace();
}

RuntimeException::~RuntimeException() throw()
{
}

const char* RuntimeException::what() const throw()
{
	string buffer;
	
	return this->toString(buffer).c_str();
}

string RuntimeException::toString(string& buffer, uint32 depth) const throw()
{
	if (depth < DEFAULT_CAUSE_DEPTH)
	{
		if (depth > 0)
		{
			buffer += "\nCaused by: ";
		}
		
		buffer += this->message;
		
		if (this->stackTrace != NULL)
		{
			for (size_t a = 0; a < this->stackTrace->getLines().size(); a++)
			{
				buffer += "\n\t";
				buffer += this->stackTrace->getLines()[a];
			}
		}
		
		if (this->cause != NULL)
		{
			this->cause->toString(buffer, ++depth);
		}
	}
	else
	{
		buffer += "\n...";
	}
	
	return buffer;
}
