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

static const char* rcsid = "$Id: ExceptionUtils.C,v 1.5 2011-08-31 06:49:27 mkotelbajcvi Exp $";

#include "ExceptionUtils.h"

StackTrace ExceptionUtils::getStackTrace(StackTrace& stackTrace, const char* caller, size_t depth)
{
	void** buffer = new void*[depth];
	size_t actualDepth = backtrace(buffer, depth);
	char** lines = backtrace_symbols(buffer, actualDepth);
	string line;
	int callerIndex = -1;
	
	for (size_t a = 0; a < actualDepth; a++)
	{
		line = string(lines[a]);
		
		if (callerIndex != -1)
		{
			stackTrace.lines.push_back(line);
		}
		else if (line.rfind(caller) != string::npos)
		{
			callerIndex = a;
		}
	}
	
	return stackTrace;
}
