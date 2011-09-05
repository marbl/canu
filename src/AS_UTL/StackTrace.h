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

#ifndef STACKTRACE_H
#define STACKTRACE_H

static const char* rcsid_STACKTRACE_H = "$Id: StackTrace.h,v 1.1 2011-09-05 16:49:45 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "StringUtils.h"

namespace Utility
{
	static const size_t DEFAULT_STACK_TRACE_DEPTH = 10;
	static const char* DEFAULT_STACK_TRACE_LINE_DELIMITER = NEWLINE_STR;
	static const char* DEFAULT_STACK_TRACE_INDENT = "";
	
	class StackTrace
	{
	public:
		StackTrace(size_t reserveDepth = DEFAULT_STACK_TRACE_DEPTH);
		
		string& toString(string& buffer, const string lineDelimiter = DEFAULT_STACK_TRACE_LINE_DELIMITER, 
			const string indent = DEFAULT_STACK_TRACE_INDENT);
		
		vector<string>& getLines()
		{
			return this->lines;
		}
		
	protected:
		vector<string> lines;
	};
}

#endif
