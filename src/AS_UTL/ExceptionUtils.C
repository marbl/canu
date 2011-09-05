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

static const char* rcsid = "$Id: ExceptionUtils.C,v 1.7 2011-09-05 16:49:45 mkotelbajcvi Exp $";

#include "ExceptionUtils.h"

using namespace Utility;

void ExceptionUtils::printStackTrace(FILE* stream, const string lineDelimiter, const string indent, const string caller, size_t depth)
{
	string str;
	
	fputs((getStackTrace(str, lineDelimiter, indent, caller, depth) + lineDelimiter).c_str(), stream);
}

string& ExceptionUtils::getStackTrace(string& buffer, const string lineDelimiter, const string indent, const string caller, size_t depth)
{
	StackTrace stackTrace;
	
	return getStackTrace(stackTrace, caller, depth).toString(buffer, lineDelimiter, indent);
}
		
StackTrace& ExceptionUtils::getStackTrace(StackTrace& stackTrace, const string caller, size_t depth)
{
#if defined(__GLIBC__)
	void** buffer = new void*[depth];
	
	size_t actualDepth = backtrace(buffer, depth);
	char** lines = backtrace_symbols(buffer, actualDepth);
	
	string line;
	bool inCaller = false, pastCaller = false;
	
	for (size_t a = 0; a < actualDepth; a++)
	{
		line = string(lines[a]);
		
		if (pastCaller)
		{
			ReflectionUtils::demangleStackEntry(line);
			
			stackTrace.getLines().push_back(line);
		}
		else if (inCaller)
		{
			if (line.rfind(caller) == string::npos)
			{
				inCaller = false;
				pastCaller = true;
				
				ReflectionUtils::demangleStackEntry(line);
				
				stackTrace.getLines().push_back(line);
			}
		}
		else if (line.rfind(caller) != string::npos)
		{
			inCaller = true;
		}
	}
#endif
	return stackTrace;
}
