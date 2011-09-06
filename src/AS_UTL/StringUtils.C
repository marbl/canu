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

static const char* rcsid = "$Id: StringUtils.C,v 1.12 2011-09-06 09:47:55 mkotelbajcvi Exp $";

#include "StringUtils.h"

using namespace Utility;

vector<string>& StringUtils::split(string str, vector<string>& buffer, const char* delimiter)
{
	vector<size_t> delimiterIndexes;
	findAll(str, delimiterIndexes, delimiter);

	size_t delimiterLength = strlen(delimiter), delimiterIndex, lastDelimiterIndex = 0;
	
	for (size_t a = 0; a < delimiterIndexes.size(); a++)
	{
		delimiterIndex = delimiterIndexes[a];
		
		if ((delimiterIndex - lastDelimiterIndex) == 0)
		{
			buffer.push_back("");
		}
		else
		{
			buffer.push_back(str.substr(lastDelimiterIndex, delimiterIndex - lastDelimiterIndex));
		}
		
		lastDelimiterIndex = delimiterIndex + delimiterLength;
	}
	
	if (lastDelimiterIndex < str.size())
	{
		buffer.push_back(str.substr(lastDelimiterIndex));
	}
	
	return buffer;
}

vector<size_t>& StringUtils::findAll(string str, vector<size_t>& buffer, size_t num, ...)
{
	initArgs(num);
	
	return findAll(str, buffer, num, VarUtils::getArgs<const char*>(num, argsList));
}

vector<size_t>& StringUtils::findAll(string str, vector<size_t>& buffer, const char* toFind)
{
	const char* toFindArr[] = { toFind };
	
	return findAll(str, buffer, 1, toFindArr);
}

vector<size_t>& StringUtils::findAll(string str, vector<size_t>& buffer, size_t num, const char** toFind)
{
	for (size_t a = 0; a < num; a++)
	{
		const char* toFindStr = toFind[a];
		size_t index = string::npos;
		
		do
		{
			index += (index != string::npos) ? strlen(toFindStr) : 0;
			index = str.find(toFindStr, index);
			
			if (index != string::npos)
			{
				buffer.push_back(index);
			}
		}
		while (index != string::npos);
	}
	
	return buffer;
}

bool StringUtils::startsWith(string str, size_t num, ...)
{
	initArgs(num);
	
	return startsWith(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

bool StringUtils::startsWith(string str, const char* toTest)
{
	const char* toTestArr[] = { toTest };
	
	return startsWith(str, 1, toTest);
}

bool StringUtils::startsWith(string str, size_t num, const char** toTest)
{
	for (size_t a = 0; a < num; a++)
	{
		string toTestStr(toTest[a]);
		
		if (str.substr(0, toTestStr.length()) == toTestStr)
		{
			return true;
		}
	}
	
	return false;
}

bool StringUtils::endsWith(string str, size_t num, ...)
{
	initArgs(num);
	
	return endsWith(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

bool StringUtils::endsWith(string str, const char* toTest)
{
	const char* toTestArr[] = { toTest };
	
	return endsWith(str, 1, toTest);
}

bool StringUtils::endsWith(string str, size_t num, const char** toTest)
{
	for (size_t a = 0; a < num; a++)
	{
		string toTestStr(toTest[a]);
		
		if (str.substr(str.length() - toTestStr.length()) == toTestStr)
		{
			return true;
		}
	}
	
	return false;
}

string& StringUtils::trim(string& str, size_t num, ...)
{
	initArgs(num);
	
	return trim(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

string& StringUtils::trim(string& str, const char* toTrim)
{
	const char* toTrimArr[] = { toTrim };
	
	return trim(str, 1, toTrimArr);
}

string& StringUtils::trim(string& str, size_t num, const char** toTrim)
{
	trimStart(str, num, toTrim);
	
	return trimEnd(str, num, toTrim);
}

string& StringUtils::trimStart(string& str, size_t num, ...)
{	
	initArgs(num);
	
	return trimStart(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

string& StringUtils::trimStart(string& str, const char* toTrim)
{
	const char* toTrimArr[] = { toTrim };
	
	return trimStart(str, 1, toTrimArr);
}

string& StringUtils::trimStart(string& str, size_t num, const char** toTrim)
{
	if (!str.empty())
	{
		for (size_t a = 0; a < num; a++)
		{
			string toTrimStr(toTrim[a]);
			
			while ((str.length() >= toTrimStr.length()) && 
					(str.substr(0, toTrimStr.length()) == toTrimStr))
			{
				str.erase(0, toTrimStr.length());
			}
		}
	}
	
	return str;
}

string& StringUtils::trimEnd(string& str, size_t num, ...)
{
	initArgs(num);
	
	return trimEnd(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

string& StringUtils::trimEnd(string& str, const char* toTrim)
{
	const char* toTrimArr[] = { toTrim };
	
	return trimEnd(str, 1, toTrimArr);
}

string& StringUtils::trimEnd(string& str, size_t num, const char** toTrim)
{
	if (!str.empty())
	{
		for (size_t a = 0; a < num; a++)
		{
			string toTrimStr(toTrim[a]);
			
			while ((str.length() >= toTrimStr.length()) && 
					(str.substr(str.length() - toTrimStr.length(), toTrimStr.length()) == toTrimStr))
			{
				str.erase(str.length() - toTrimStr.length());
			}
		}
	}
	
	return str;
}
