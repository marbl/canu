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

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sstream>
#include <string>

using namespace std;

#include "AS_global.h"
#include "VarUtils.h"

static const char *RCSID_STRINGUTILS_H = "$Id: StringUtils.h,v 1.6 2011-08-04 14:34:41 mkotelbajcvi Exp $";

#define NULL_TERMINATOR '\0'

class StringUtils
{
public:
	inline static bool isBlank(const char* str)
	{
		if (!isEmpty(str))
		{
			for (size a = 0; a < strlen(str); a++)
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
		return str == NULL;
	}

	inline static bool isEmpty(string str)
	{
		return str.length() == 0;
	}

	inline static bool areEqual(const char* str1, const char* str2)
	{
		return (str1 == NULL) ? str2 == NULL : (str2 != NULL) && (strcmp(str1, str2) == 0);
	}
	
	static bool startsWith(const char* str, size num, ...);
	static bool startsWith(const char* str, size num, const char** toTest);
	static bool endsWith(const char* str, size num, ...);
	static bool endsWith(const char* str, size num, const char** toTest);

	static const char* concat(size num, ...);
	
	static const char* join(const char* delimiter, size num, ...);
	static const char* join(const char* delimiter, size num, const char** toJoin);
	
	static const char* trim(const char* str, size num, ...);
	static const char* trim(const char* str, size num, const char** toTrim);
	static const char* trimStart(const char* str, size num, ...);
	static const char* trimStart(const char* str, size num, const char** toTrim);
	static const char* trimEnd(const char* str, size num, ...);
	static const char* trimEnd(const char* str, size num, const char** toTrim);
	
	inline static const char* toString(unsigned value)
	{
		ostringstream stream(ostringstream::out);
		
		stream << value;
		
		return toString(stream.str());
	}

	inline static const char* toString(unsigned long value)
	{
		ostringstream stream(ostringstream::out);
		
		stream << value;
		
		return toString(stream.str());
	}

	inline static const char* toString(int value)
	{
		ostringstream stream(ostringstream::out);
		
		stream << value;
		
		return toString(stream.str());
	}

	inline static const char* toString(long value)
	{
		ostringstream stream(ostringstream::out);
		
		stream << value;
		
		return toString(stream.str());
	}

	inline static const char* toString(float value)
	{
		ostringstream stream(ostringstream::out);
		
		stream << value;
		
		return toString(stream.str());
	}

	inline static const char* toString(double value)
	{
		ostringstream stream(ostringstream::out);
		
		stream << value;
		
		return toString(stream.str());
	}

	inline static const char* toString(char value)
	{
		char* str = new char[2];
		str[0] = value;
		str[1] = NULL_TERMINATOR;
		
		return str;
	}
	
	inline static const char* toString(string value)
	{
		char* str = new char[value.length()];
		
		strcpy(str, value.c_str());
		
		return str;
	}
};

#endif
