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

static const char* rcsid = "$Id: FileUtils.C,v 1.9 2011-09-06 17:14:15 mkotelbajcvi Exp $";

#include "FileUtils.h"

// TEMP: stopgap commit
#ifdef __APPLE__
ssize_t getline(char** buffer, size_t* bufferSize, FILE* stream)
{
	return -1;
}
#endif

using namespace Utility;

string& FileUtils::getPath(string& buffer, size_t num, ...)
{
	initArgs(num);
	
	return getPath(buffer, num, VarUtils::getArgs<const char*>(num, argsList));
}

string& FileUtils::getPath(string& buffer, size_t num, const char** pathParts)
{
	string pathPartStr;
	
	for (size_t a = 0; a < num; a++)
	{
		pathPartStr = pathParts[a];
		
		if (a != 0)
		{
			pathPartStr = StringUtils::trimStart(pathPartStr, 1, PATH_DELIMITER);
		}
		
		if (a < (num - 1))
		{
			pathPartStr = StringUtils::trimEnd(pathPartStr, 1, PATH_DELIMITER);
		}
		
		pathParts[a] = pathPartStr.c_str();
	}
	
	return StringUtils::join(buffer, PATH_DELIMITER_STR, num, pathParts);
}
