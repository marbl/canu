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

#include "StringUtils.h"

static const char *RCSID = "$Id: StringUtils.C,v 1.3 2011-08-01 16:54:03 mkotelbajcvi Exp $";

const char* StringUtils::trim(const char* str, size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	const char* trimmedStr = trim(str, num, getArgs(num, (const char**)safe_calloc(num, sizeof(const char*)), argsList));
	
	va_end(argsList);
		
	return trimmedStr;
}

const char* StringUtils::trim(const char* str, size num, const char** toTrim)
{
	return trimEnd(trimStart(str, num, toTrim), num, toTrim);
}

const char* StringUtils::trimStart(const char* str, size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	const char* trimmedStr = trimStart(str, num, getArgs(num, (const char**)safe_calloc(num, sizeof(const char*)), argsList));
	
	va_end(argsList);
	
	return trimmedStr;
}

const char* StringUtils::trimStart(const char* str, size num, const char** toTrim)
{
	if (!isEmpty(str))
	{
		string strObj(str);
		
		for (size a = 0; a < num; a++)
		{
			const char* toTrimItem = toTrim[a];
			
			size toTrimItemLength = strlen(toTrimItem);
			
			while ((strObj.length() >= toTrimItemLength) && 
					areEqual(strObj.substr(0, toTrimItemLength).c_str(), toTrimItem))
			{
				strObj = strObj.erase(0, toTrimItemLength);
			}
		}
		
		return strObj.c_str();
	}
	
	return NULL;
}

const char* StringUtils::trimEnd(const char* str, size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	return trimEnd(str, num, getArgs(num, (const char**)safe_calloc(num, sizeof(const char*)), argsList));
}

const char* StringUtils::trimEnd(const char* str, size num, const char** toTrim)
{
	if (!isEmpty(str))
	{
		string strObj(str);
		
		for (size a = 0; a < num; a++)
		{
			const char* toTrimItem = toTrim[a];
			
			size toTrimItemLength = strlen(toTrimItem);
			
			while ((strObj.length() >= toTrimItemLength) && 
					areEqual((char*)strObj.substr(strObj.length() - toTrimItemLength, toTrimItemLength).c_str(), toTrimItem))
			{
				strObj = strObj.erase(strObj.length() - toTrimItemLength);
			}
		}
		
		return strObj.c_str();
	}
	
	return NULL;
}

const char* StringUtils::concat(size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	string str;
	const char** items = getArgs(num, (const char**)safe_calloc(num, sizeof(const char**)), argsList);
	
	for (size a = 0; a < num; a++)
	{
		str += items[a];
	}
	
	return toString(str);
}

template<class T>
T* getArgs(size num, T* args, va_list& argsList)
{
	T arg;
	
	for (size a = 0; (a < num) && ((arg = va_arg(argsList, T)) != NULL); a++)
	{
		args[a] = arg;
	}
	
	return args;
}
