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

static const char* rcsid = "$Id: StringUtils.C,v 1.8 2011-08-10 20:25:15 mkotelbajcvi Exp $";

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

bool StringUtils::startsWith(const char* str, size_t num, ...)
{
	initArgs(num);
	
	return startsWith(str, num, VarUtils::getArgs<const char*>(num, argsList));
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

bool StringUtils::endsWith(const char* str, size_t num, ...)
{
	initArgs(num);
	
	return endsWith(str, num, VarUtils::getArgs<const char*>(num, argsList));
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

const char* StringUtils::trim(const char* str, size_t num, ...)
{
	initArgs(num);
	
	return trim(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

const char* StringUtils::trim(const char* str, size_t num, const char** toTrim)
{
	return trimEnd(trimStart(str, num, toTrim), num, toTrim);
}

const char* StringUtils::trimStart(const char* str, size_t num, ...)
{	
	initArgs(num);
	
	return trimStart(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

const char* StringUtils::trimStart(const char* str, size_t num, const char** toTrim)
{
	if (!isEmpty(str))
	{
		string strObj(str);
		
		for (size_t a = 0; a < num; a++)
		{
			string toTrimItem(toTrim[a]);
			
			while ((strObj.length() >= toTrimItem.length()) && 
					(strObj.substr(0, toTrimItem.length()) == toTrimItem))
			{
				strObj.erase(0, toTrimItem.length());
			}
		}
		
		return toString(strObj);
	}
	
	return NULL;
}

const char* StringUtils::trimEnd(const char* str, size_t num, ...)
{
	initArgs(num);
	
	return trimEnd(str, num, VarUtils::getArgs<const char*>(num, argsList));
}

const char* StringUtils::trimEnd(const char* str, size_t num, const char** toTrim)
{
	if (!isEmpty(str))
	{
		string strObj(str);
		
		for (size_t a = 0; a < num; a++)
		{
			string toTrimItem(toTrim[a]);
			
			while ((strObj.length() >= toTrimItem.length()) && 
					(strObj.substr(strObj.length() - toTrimItem.length(), toTrimItem.length()) == toTrimItem))
			{
				strObj = strObj.erase(strObj.length() - toTrimItem.length());
			}
		}
		
		return toString(strObj);
	}
	
	return NULL;
}
