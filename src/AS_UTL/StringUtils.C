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

static const char* rcsid = "$Id: StringUtils.C,v 1.6 2011-08-04 18:18:56 mkotelbajcvi Exp $";

#include "StringUtils.h"

bool StringUtils::startsWith(const char* str, size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	return startsWith(str, num, VarUtils::getArgs(num, new const char*[num], argsList));
}

bool StringUtils::startsWith(const char* str, size num, const char** toTest)
{
	string strObj(str);
	
	for (size a = 0; a < num; a++)
	{
		string toTestObj(toTest[a]);
		
		if (strObj.substr(0, toTestObj.length()) == toTestObj)
		{
			return TRUE;
		}
	}
	
	return FALSE;
}

bool StringUtils::endsWith(const char* str, size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	return endsWith(str, num, VarUtils::getArgs(num, new const char*[num], argsList));
}

bool StringUtils::endsWith(const char* str, size num, const char** toTest)
{
	string strObj(str);
	
	for (size a = 0; a < num; a++)
	{
		string toTestObj(toTest[a]);
		
		if (strObj.substr(strObj.length() - toTestObj.length()) == toTestObj)
		{
			return TRUE;
		}
	}
	
	return FALSE;
}

const char* StringUtils::concat(size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	string str;
	const char** items = VarUtils::getArgs(num, (const char**)safe_calloc(num, sizeof(const char**)), argsList);
	
	for (size a = 0; a < num; a++)
	{
		str += items[a];
	}
	
	return toString(str);
}

const char* StringUtils::join(const char* delimiter, size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	return join(delimiter, num, VarUtils::getArgs(num, (const char**)safe_calloc(num, sizeof(const char*)), argsList));
}

const char* StringUtils::join(const char* delimiter, size num, const char** toJoin)
{
	string str;
	
	for (size a = 0; a < num; a++)
	{
		if (!str.empty())
		{
			str += delimiter;
		}
		
		str += toJoin[a];
	}
	
	return toString(str);
}

const char* StringUtils::trim(const char* str, size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	return trim(str, num, VarUtils::getArgs(num, (const char**)safe_calloc(num, sizeof(const char*)), argsList));
}

const char* StringUtils::trim(const char* str, size num, const char** toTrim)
{
	return trimEnd(trimStart(str, num, toTrim), num, toTrim);
}

const char* StringUtils::trimStart(const char* str, size num, ...)
{	
	va_list argsList;
	va_start(argsList, num);
	
	return trimStart(str, num, VarUtils::getArgs(num, (const char**)safe_calloc(num, sizeof(const char*)), argsList));
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
		
		return toString(strObj);
	}
	
	return NULL;
}

const char* StringUtils::trimEnd(const char* str, size num, ...)
{
	va_list argsList;
	va_start(argsList, num);
	
	return trimEnd(str, num, VarUtils::getArgs(num, (const char**)safe_calloc(num, sizeof(const char*)), argsList));
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
		
		return toString(strObj);
	}
	
	return NULL;
}
