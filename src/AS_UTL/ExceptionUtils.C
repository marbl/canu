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

#include "ExceptionUtils.h"

static const char* RCSID = "$Id: ExceptionUtils.C,v 1.2 2011-08-04 17:08:32 brianwalenz Exp $";

StackTrace& ExceptionUtils::getStackTrace(const char* caller, size depth)
{
	StackTrace* stackTrace = new StackTrace();
	stackTrace->depth = 0;
	stackTrace->lines = NULL;
	
	void** buffer = new void*[depth];
	size actualDepth = backtrace(buffer, depth);
	char** lines = backtrace_symbols(buffer, actualDepth);
	
	int callerIndex = -1;
	
	for (size a = 0; a < actualDepth; a++)
	{
		if (callerIndex != -1)
		{
			stackTrace->lines[a - callerIndex - 1] = new char[strlen(lines[a]) + 1];
			strcpy(stackTrace->lines[a - callerIndex - 1], lines[a]);
		}
		else if (string(lines[a]).rfind(caller) != string::npos)
		{
			callerIndex = a;
			
			stackTrace->depth = actualDepth - callerIndex - 1;
			stackTrace->lines = new char*[stackTrace->depth];
		}
	}
	
	return *stackTrace;
}
