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

static const char* rcsid = "$Id: StringUtils.C,v 1.9 2011-08-30 23:09:51 mkotelbajcvi Exp $";

#include "StringUtils.h"

vector<size_t> StringUtils::findAll(const char* str, size_t num, ...)
{
	initArgs(num);
	
	return findAll(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

vector<size_t> StringUtils::findAll(const char* str, size_t num, const char** toFind)
{
	string strObj(str);
	vector<size_t> found;
	
	for (size_t a = 0; a < num; a++)
	{
		const char* toFindItem = toFind[a];
		size_t index = string::npos;
		
		do
		{
			index += (index != string::npos) ? strlen(toFindItem) : 0;
			index = strObj.find(toFindItem, index);
			
			if (index != string::npos)
			{
				found.push_back(index);
			}
		}
		while (index != string::npos);
	}
	
	return found;
}

bool StringUtils::startsWith(string str, size_t num, ...)
{
	initArgs(num);
	
	return startsWith(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

bool StringUtils::startsWith(const char* str, size_t num, ...)
{
	initArgs(num);
	
	return startsWith(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

bool StringUtils::startsWith(string str, size_t num, const char** toTest)
{
	return startsWith(str.c_str(), num, toTest);
}

bool StringUtils::startsWith(const char* str, size_t num, const char** toTest)
{
	string strObj(str);
	
	for (size_t a = 0; a < num; a++)
	{
		string toTestObj(toTest[a]);
		
		if (strObj.substr(0, toTestObj.length()) == toTestObj)
		{
			return TRUE;
		}
	}
	
	return FALSE;
}

bool StringUtils::endsWith(string str, size_t num, ...)
{
	initArgs(num);
	
	return endsWith(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

bool StringUtils::endsWith(const char* str, size_t num, ...)
{
	initArgs(num);
	
	return endsWith(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

bool StringUtils::endsWith(string str, size_t num, const char** toTest)
{
	return endsWith(str.c_str(), num, toTest);
}

bool StringUtils::endsWith(const char* str, size_t num, const char** toTest)
{
	string strObj(str);
	
	for (size_t a = 0; a < num; a++)
	{
		string toTestObj(toTest[a]);
		
		if (strObj.substr(strObj.length() - toTestObj.length()) == toTestObj)
		{
			return TRUE;
		}
	}
	
	return FALSE;
}

string StringUtils::trim(string& str, size_t num, ...)
{
	initArgs(num);
	
	return trim(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

string StringUtils::trim(string& str, size_t num, const char** toTrim)
{
	trimStart(str, num, toTrim);
	
	return trimEnd(str, num, toTrim);
}

string StringUtils::trimStart(string& str, size_t num, ...)
{	
	initArgs(num);
	
	return trimStart(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

string StringUtils::trimStart(string& str, size_t num, const char** toTrim)
{
	if (!str.empty())
	{
		for (size_t a = 0; a < num; a++)
		{
			string toTrimItem(toTrim[a]);
			
			while ((str.length() >= toTrimItem.length()) && 
					(str.substr(0, toTrimItem.length()) == toTrimItem))
			{
				str.erase(0, toTrimItem.length());
			}
		}
	}
	
	return str;
}

string StringUtils::trimEnd(string& str, size_t num, ...)
{
	initArgs(num);
	
	return trimEnd(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

string StringUtils::trimEnd(string& str, size_t num, const char** toTrim)
{
	if (!str.empty())
	{
		for (size_t a = 0; a < num; a++)
		{
			string toTrimItem(toTrim[a]);
			
			while ((str.length() >= toTrimItem.length()) && 
					(str.substr(str.length() - toTrimItem.length(), toTrimItem.length()) == toTrimItem))
			{
				str = str.erase(str.length() - toTrimItem.length());
			}
		}
	}
	
	return str;
}
