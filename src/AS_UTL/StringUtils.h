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

static const char* rcsid_STRINGUTILS_H = "$Id: StringUtils.h,v 1.17 2011-09-06 01:11:56 mkotelbajcvi Exp $";

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "VarUtils.h"

namespace Utility
{
	static const char* DEFAULT_TO_STRING_CLOCK_DELIMITER = " ";
	static const size_t MILLESECONDS_IN_SECOND = 1000;
	static const size_t SECONDS_IN_MINUTE = 60;
	static const size_t MILLESECONDS_IN_MINUTE = MILLESECONDS_IN_SECOND * SECONDS_IN_MINUTE;
	static const size_t MINUTES_IN_HOUR = 60;
	static const size_t MILLESECONDS_IN_HOUR = MILLESECONDS_IN_MINUTE * MINUTES_IN_HOUR;
	static const size_t HOURS_IN_DAY = 24;
	static const size_t MILLESECONDS_IN_DAY = MILLESECONDS_IN_HOUR * HOURS_IN_DAY;
	
	class StringUtils
	{
	public:
		inline static bool isNewline(char character)
		{
			return (character == LINE_FEED) || (character == CARRIAGE_RETURN);
		}
		
		inline static bool isLowercase(string str)
		{
			for (size_t a = 0; a < str.length(); a++)
			{
				if (!islower(str[a]))
				{
					return false;
				}
			}
			
			return true;
		}
		
		inline static bool isUppercase(string str)
		{
			for (size_t a = 0; a < str.length(); a++)
			{
				if (!isupper(str[a]))
				{
					return false;
				}
			}
			
			return true;
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
		
		inline static string& delimit(string& buffer, const string delimiter)
		{
			if (!delimiter.empty() && !buffer.empty())
			{
				buffer += delimiter;
			}
			
			return buffer;
		}
		
		inline static string& join(string& buffer, const string delimiter, size_t num, ...)
		{
			initArgs(num);
			
			return join(buffer, delimiter, num, VarUtils::getArgs<const char*>(num, argsList));
		}
		
		inline static string& join(string& buffer, const string delimiter, size_t num, const char** toJoin)
		{
			for (size_t a = 0; a < num; a++)
			{
				delimit(buffer, delimiter);
				
				buffer += toJoin[a];
			}
			
			return buffer;
		}
		
		template<class T>
		inline static string& join(string& buffer, const string delimiter, T begin, T end)
		{
			for (; begin != end; begin++)
			{
				delimit(buffer, delimiter);
				
				buffer += *begin;
			}
			
			return buffer;
		}
		
		inline static string& toString(string& buffer, clock_t value, bool asElapsed = false, string delimiter = DEFAULT_TO_STRING_CLOCK_DELIMITER)
		{
			double valueDbl = ((asElapsed ? (double)(clock() - value) : (double)value) / (double)CLOCKS_PER_SEC) * 
				MILLESECONDS_IN_SECOND;
			size_t milliseconds, seconds, minutes, hours, days;
			
			days = (size_t)(valueDbl / MILLESECONDS_IN_DAY);
			valueDbl -= (days * MILLESECONDS_IN_DAY);
			
			if (days)
			{
				delimit(buffer, delimiter);
				toString(days, buffer);
				buffer += "d";
			}
			
			hours = (size_t)(valueDbl / MILLESECONDS_IN_HOUR);
			valueDbl -= (hours * MILLESECONDS_IN_HOUR);
			
			if (hours)
			{
				delimit(buffer, delimiter);
				toString(hours, buffer);
				buffer += "h";
			}
			
			minutes = (size_t)(valueDbl / MILLESECONDS_IN_MINUTE);
			valueDbl -= (minutes * MILLESECONDS_IN_MINUTE);
			
			if (minutes)
			{
				delimit(buffer, delimiter);
				toString(minutes, buffer);
				buffer += "m";
			}
			
			seconds = (size_t)(valueDbl / MILLESECONDS_IN_SECOND);
			valueDbl -= (seconds * MILLESECONDS_IN_SECOND);
			
			if (seconds)
			{
				delimit(buffer, delimiter);
				toString(seconds, buffer);
				buffer += "s";
			}
			
			milliseconds = (size_t)valueDbl;
			
			if (milliseconds || buffer.empty())
			{
				delimit(buffer, delimiter);
				toString(milliseconds, buffer);
				buffer += "ms";
			}
			
			return buffer;
		}
		
		inline static string& toString(char value, string& buffer)
		{
			buffer += value;
			
			return buffer;
		}
		
		template<class T>
		inline static string& toString(T value, string& buffer)
		{
			ostringstream stream(ostringstream::out);
			
			stream << value;
			
			buffer += stream.str();
			
			return buffer;
		}
		
		static vector<string>& split(string str, vector<string>& buffer, size_t num, ...);
		static vector<string>& split(string str, vector<string>& buffer, const char* delimiter);
		static vector<string>& split(string str, vector<string>& buffer, size_t num, const char** delimiters);
		
		static vector<size_t>& findAll(string str, vector<size_t>& buffer, size_t num, ...);
		static vector<size_t>& findAll(string str, vector<size_t>& buffer, const char* toFind);
		static vector<size_t>& findAll(string str, vector<size_t>& buffer, size_t num, const char** toFind);
		
		static bool startsWith(string str, size_t num, ...);
		static bool startsWith(string str, const char* toTest);
		static bool startsWith(string str, size_t num, const char** toTest);
		static bool endsWith(string str, size_t num, ...);
		static bool endsWith(string str, const char* toTest);
		static bool endsWith(string str, size_t num, const char** toTest);
		
		static string& trim(string& str, size_t num, ...);
		static string& trim(string& str, const char* toTrim);
		static string& trim(string& str, size_t num, const char** toTrim);
		static string& trimStart(string& str, size_t num, ...);
		static string& trimStart(string& str, const char* toTrim);
		static string& trimStart(string& str, size_t num, const char** toTrim);
		static string& trimEnd(string& str, size_t num, ...);
		static string& trimEnd(string& str, const char* toTrim);
		static string& trimEnd(string& str, size_t num, const char** toTrim);
		
	private:
		StringUtils()
		{
		}
	};
}

#endif
