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

#ifndef STRINGUTILS_H
#define STRINGUTILS_H

static const char* rcsid_STRINGUTILS_H = "$Id: StringUtils.h,v 1.12 2011-08-31 06:49:27 mkotelbajcvi Exp $";

#include <stdarg.h>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "VarUtils.h"

#define NEWLINE '\n'
#define NULL_TERMINATOR '\0'

class StringUtils
{
public:
	inline static bool isLowercase(const char* str)
	{
		for (size_t a = 0; a < strlen(str); a++)
		{
			if (!isLowercase(str[a]))
			{
				return false;
			}
		}
		
		return true;
	}
	
	inline static bool isLowercase(const char character)
	{
		return (character >= 'a') && (character <= 'z');
	}
	
	inline static bool isUppercase(const char* str)
	{
		for (size_t a = 0; a < strlen(str); a++)
		{
			if (!isUppercase(str[a]))
			{
				return false;
			}
		}
		
		return true;
	}
	
	inline static bool isUppercase(const char character)
	{
		return (character >= 'A') && (character <= 'Z');
	}
	
	inline static bool isBlank(string str)
	{
		return isBlank(str.c_str());
	}
	
	inline static bool isBlank(const char* str)
	{
		if (!isEmpty(str))
		{
			for (size_t a = 0; a < strlen(str); a++)
			{
				if (!isspace(str[a]))
				{
					return false;
				}
			}
		}
		
		return true;
	}

	inline static bool isEmpty(const char* str)
	{
		return (str == NULL) || strlen(str) == 0;
	}

	inline static bool areEqual(const char* str1, const char* str2)
	{
		return (str1 != NULL) && (str2 != NULL) && (strcmp(str1, str2) == 0);
	}
	
	inline static string concat(string& buffer, const char* format, ...)
	{
		initArgs(format);
		
		vsprintf((char*)buffer.data(), format, argsList);
		
		buffer.assign(buffer.data());
		
		va_end(argsList);
		
		return buffer;
	}
	
	inline static string join(const char* delimiter, string& buffer, size_t num, ...)
	{
		initArgs(num);
		
		return join(delimiter, buffer, num, VarUtils::getArgs<const char*>(num, argsList));
	}
	
	inline static string join(const char* delimiter, string& buffer, size_t num, const char** toJoin)
	{
		for (size_t a = 0; a < num; a++)
		{
			if (!isEmpty(delimiter) && !buffer.empty())
			{
				buffer += delimiter;
			}
			
			buffer += toJoin[a];
		}
		
		return buffer;
	}
	
	inline static string toString(unsigned value, string& buffer)
	{
		ostringstream stream(ostringstream::out);
		stream.str(buffer);
		
		stream << value;
		
		return stream.str();
	}

	inline static string toString(unsigned long value, string& buffer)
	{
		ostringstream stream(ostringstream::out);
		stream.str(buffer);
		
		stream << value;
		
		return stream.str();
	}

	inline static string toString(int value, string& buffer)
	{
		ostringstream stream(ostringstream::out);
		stream.str(buffer);
		
		stream << value;
		
		return stream.str();
	}
	
	inline static string toString(long value, string& buffer)
	{
		ostringstream stream(ostringstream::out);
		stream.str(buffer);
		
		stream << value;
		
		return stream.str();
	}

	inline static string toString(float value, string& buffer)
	{
		ostringstream stream(ostringstream::out);
		stream.str(buffer);
		
		stream << value;
		
		return stream.str();
	}

	inline static string toString(double value, string& buffer)
	{
		ostringstream stream(ostringstream::out);
		stream.str(buffer);
		
		stream << value;
		
		return stream.str();
	}

	inline static string toString(char value, string& buffer)
	{
		buffer += value;
		
		return buffer;
	}
	
	// TODO: remove
	inline static const char* toString(string value)
	{
		char* str = new char[value.length()];
		
		strcpy(str, value.c_str());
		
		return str;
	}
	
	static vector<string> split(string str, vector<string>& buffer, size_t num, ...);
	static vector<string> split(string str, vector<string>& buffer, size_t num, const char** delimiters);
	
	static vector<size_t> findAll(const char* str, vector<size_t>& buffer, size_t num, ...);
	static vector<size_t> findAll(string str, vector<size_t>& buffer, size_t num, ...);
	static vector<size_t> findAll(const char* str, vector<size_t>& buffer, size_t num, const char** toFind);
	static vector<size_t> findAll(string str, vector<size_t>& buffer, size_t num, const char** toFind);
	
	static bool startsWith(const char* str, size_t num, ...);
	static bool startsWith(string str, size_t num, ...);
	static bool startsWith(const char* str, size_t num, const char** toTest);
	static bool startsWith(string str, size_t num, const char** toTest);
	static bool endsWith(const char* str, size_t num, ...);
	static bool endsWith(string str, size_t num, ...);
	static bool endsWith(const char* str, size_t num, const char** toTest);
	static bool endsWith(string str, size_t num, const char** toTest);
	
	static string trim(string& str, size_t num, ...);
	static string trim(string& str, size_t num, const char** toTrim);
	static string trimStart(string& str, size_t num, ...);
	static string trimStart(string& str, size_t num, const char** toTrim);
	static string trimEnd(string& str, size_t num, ...);
	static string trimEnd(string& str, size_t num, const char** toTrim);
	
private:
	StringUtils()
	{
	}
};

#endif
