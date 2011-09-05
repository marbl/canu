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

static const char* rcsid = "$Id: ReflectionUtils.C,v 1.1 2011-09-05 16:49:45 mkotelbajcvi Exp $";

#include "ReflectionUtils.h"

using namespace Utility;

bool ReflectionUtils::demangleStackEntry(string& entry, bool includePath)
{
	size_t start = entry.find(STACK_ENTRY_MANGLED_NAME_PREFIX), 
		end = entry.rfind(STACK_ENTRY_MANGLED_NAME_SUFFIX);
	
	if ((start == string::npos) || (end == string::npos))
	{
		return false;
	}
	
	string entryPath = entry.substr(0, start), 
		entryName = entry.substr(start + 1, end - start - 1);
	
	bool result = demangle(entryName);
	
	entry.clear();
	
	if (entryName.find(STACK_ENTRY_OBJECT_DELIMITER) == -1)
	{
		if (includePath)
		{
			entry += entryPath;
		}
		else
		{
			size_t lastPathDelimiter = entryPath.rfind(PATH_DELIMITER);
			
			if (lastPathDelimiter != string::npos)
			{
				entry += entryPath.substr(lastPathDelimiter + 1, entryPath.size() - lastPathDelimiter - 1);
			}
		}
		
		entry += STACK_ENTRY_OBJECT_DELIMITER;
	}
	
	entry += entryName;
		
	return result;
}

bool ReflectionUtils::demangle(string& name)
{
	if ((name == STACK_ENTRY_MAIN_NAME) || (name == STACK_ENTRY_LIBC_START_MAIN_NAME) || (name == STACK_ENTRY_GXX_PERSONALITY_NAME))
	{
		name += STACK_ENTRY_FUNCTION_SUFFIX;
		
		return true;
	}
	
#if defined(__GLIBCXX__)
	size_t bufferSize = INITIAL_DEMANGLING_NAME_BUFFER_SIZE;
	char* buffer = (char*)safe_malloc(bufferSize);
	int status;
	
	buffer = __cxa_demangle(name.c_str(), buffer, &bufferSize, &status);
	
	if (status == 0)
	{
		name = buffer;
	}
	
	safe_free(buffer);
	
	return (status == 0);
#else
	return false;
#endif
}
