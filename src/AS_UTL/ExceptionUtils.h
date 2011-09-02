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

#ifndef EXCEPTIONUTILS_H
#define EXCEPTIONUTILS_H

static const char* rcsid_EXCEPTIONUTILS_H = "$Id: ExceptionUtils.h,v 1.8 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <exception>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "StringUtils.h"

#if defined(__GLIBC__)
#include <execinfo.h>
#else
static int backtrace(void** buffer, int size) { return 0; }
static char** backtrace_symbols(void* const* buffer, int size) { return NULL; }
#endif

namespace Utility
{
	static const char* DEFAULT_STACK_TRACE_LINE_DELIMITER = NEWLINE_STR;
	static const char* DEFAULT_STACK_TRACE_INDENT = "";
	static const char* DEFAULT_STACK_TRACE_CALLER = "ExceptionUtils";
	static const size_t DEFAULT_STACK_TRACE_DEPTH = 10;
	
	class StackTrace
	{
	public:
		StackTrace(size_t reserveDepth = DEFAULT_STACK_TRACE_DEPTH)
		{
			this->lines.reserve(reserveDepth);
		}
		
		string toString()
		{
			string str;
			
			return StringUtils::join(NEWLINE_STR, str, this->lines.begin(), this->lines.end());
		}
		
		vector<string>& getLines()
		{
			return this->lines;
		}
		
	protected:
		vector<string> lines;
	};
	
	class ExceptionUtils
	{
	public:
		__inline__ static void printStackTrace(FILE* stream = stderr, const char* lineDelimiter = DEFAULT_STACK_TRACE_LINE_DELIMITER, 
			const char* indent = DEFAULT_STACK_TRACE_INDENT, const char* caller = DEFAULT_STACK_TRACE_CALLER, size_t depth = DEFAULT_STACK_TRACE_DEPTH)
		{
			string str;
			
			fputs((getStackTrace(str, lineDelimiter, indent, caller, depth) + lineDelimiter).c_str(), stream);
		}
		
		__inline__ static string getStackTrace(string& buffer, const char* lineDelimiter = DEFAULT_STACK_TRACE_LINE_DELIMITER, 
			const char* indent = DEFAULT_STACK_TRACE_INDENT, const char* caller = DEFAULT_STACK_TRACE_CALLER, size_t depth = DEFAULT_STACK_TRACE_DEPTH)
		{
			StackTrace stackTrace;
			getStackTrace(stackTrace, caller, depth);
			
			return indent + StringUtils::join((string(lineDelimiter) + indent).c_str(), buffer, stackTrace.getLines().begin(), stackTrace.getLines().end());
		}
		
		__inline__ static StackTrace getStackTrace(StackTrace& stackTrace, const char* caller = DEFAULT_STACK_TRACE_CALLER, 
			size_t depth = DEFAULT_STACK_TRACE_DEPTH)
		{
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
					stackTrace.getLines().push_back(line);
				}
				else if (inCaller)
				{
					if (line.rfind(caller) == string::npos)
					{
						inCaller = false;
						pastCaller = true;
						
						stackTrace.getLines().push_back(line);
					}
				}
				else if (line.rfind(caller) != string::npos)
				{
					inCaller = true;
				}
			}
			
			return stackTrace;
		}
		
	private:
		ExceptionUtils()
		{
		}
	};
}

#endif
