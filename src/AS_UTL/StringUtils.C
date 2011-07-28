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

static const char *RCSID = "$Id: StringUtils.C,v 1.2 2011-07-28 11:31:00 mkotelbajcvi Exp $";

bool StringUtils::isBlank(const char* str)
{
	if (!isEmpty(str))
	{
		bool notWhitespace;
		
		for (size a = 0; a < strlen(str); a++)
		{
			notWhitespace = true;
			
			for (size b = 0; b < strlen(WHITESPACE_CHARS); b++)
			{
				if (str[a] == WHITESPACE_CHARS[b])
				{
					notWhitespace = false;
					
					break;
				}
			}
			
			if (notWhitespace)
			{
				return false;
			}
		}
	}
	
	return true;
}

bool StringUtils::isEmpty(const char* str)
{
	return str == NULL;
}

bool StringUtils::isEmpty(std::string str)
{
	return str.length() == 0;
}

bool StringUtils::areEqual(const char* str1, const char* str2)
{
	return (str1 == NULL) ? str2 == NULL : (str2 != NULL) && (strcmp(str1, str2) == 0);
}

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
	if (toTrim == NULL)
	{
		throw ArgumentException("toTrim", "Array of strings to trim must not be null.");
	}
	
	if (!isEmpty(str))
	{
		std::string strObj(str);
		
		for (size a = 0; a < num; a++)
		{
			const char* toTrimItem = toTrim[a];
			
			if (toTrimItem == NULL)
			{
				char error[64];
				sprintf(error, "String to trim must not be null: index="F_UL, a);
				
				throw ArgumentException("toTrim", error);
			}
			
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
	if (toTrim == NULL)
	{
		throw  ArgumentException("toTrim", "Array of strings to trim must not be null.");
	}
	
	if (!isEmpty(str))
	{
		std::string strObj(str);
		
		for (size a = 0; a < num; a++)
		{
			const char* toTrimItem = toTrim[a];
			
			if (toTrimItem == NULL)
			{
				char error[64];
				sprintf(error, "String to trim must not be null: index="F_UL, a);

				throw ArgumentException("toTrim", error);
			}
			
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

const char* StringUtils::toString(int value)
{
	return toString<int>(value);
}

const char* StringUtils::toString(long value)
{
	return toString<long>(value);
}

const char* StringUtils::toString(float value)
{
	return toString<float>(value);
}

const char* StringUtils::toString(double value)
{
	return toString<double>(value);
}

const char* StringUtils::toString(char value)
{
	return toString<char>(value);
}

template<class T>
const char* StringUtils::toString(T value)
{
	std::ostringstream stream(std::ostringstream::out);
	
	stream << value;
	
	return stream.str().c_str();
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
