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

static const char* rcsid = "$Id: FileUtils.C,v 1.6 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "FileUtils.h"

using namespace Utility;

/*
char* FileUtils::readLine(FILE* file, char* buffer, size_t bufferSize, bool includeNewline)
{
	if (file == NULL)
	{
		throw ArgumentException("File handle must not be null.");
	}
	
	if (feof(file))
	{
		return NULL;
	}
	
	fgets(buffer, bufferSize, file);
	
	size_t bufferSizeUsed = strlen(buffer);
	
	if (!includeNewline && (bufferSizeUsed > 0) && (buffer[bufferSizeUsed - 1] == NEWLINE))
	{	
		buffer[bufferSizeUsed - 1] = NULL_TERMINATOR;
	}
	
	return buffer;
}
*/

string& FileUtils::readLine(FILE* file, string& buffer, bool includeNewline)
{
	// TODO: implement
	
	return buffer;
}

void FileUtils::writeLine(FILE* file, string& buffer)
{
	// TODO: implement
}

bool FileUtils::canRead(FILE* file)
{
	// TODO: implement
	
	return false;
}

bool FileUtils::canWrite(FILE* file)
{
	// TODO: implement
	
	return false;
}

bool FileUtils::isDirectory(string path)
{
	return isType(path, S_IFDIR);
}

bool FileUtils::isFifo(string path)
{
	return isType(path, S_IFIFO);
}

bool FileUtils::isFile(string path)
{
	return isType(path, S_IFREG);
}

bool FileUtils::isLink(string path)
{
	return isType(path, S_IFLNK);
}

bool FileUtils::isSocket(string path)
{
	return isType(path, S_IFSOCK);
}

bool FileUtils::isType(string path, mode_t type)
{
	Stats stats;
	
	return getStats(path, stats) && ((stats.st_mode & S_IFMT) == type);
}

string& FileUtils::getStats(string path, Stats& stats, string& buffer)
{
	if (!getStats(path, stats))
	{
		buffer += ErrorUtils::getError();
	}
	
	return buffer;
}

bool FileUtils::getStats(string path, Stats& stats)
{
	return exists(path) && (stat(path.c_str(), &stats) != -1);
}

bool FileUtils::exists(string path)
{
	return isAccessible(path, F_OK);
}

string& FileUtils::isAccessible(string path, int accessFlag, string& buffer)
{
	if (!isAccessible(path, accessFlag))
	{
		buffer += ErrorUtils::getError();
	}
	
	return buffer;
}

bool FileUtils::isAccessible(string path, int accessFlag)
{
	return isValidPath(path) && (access(path.c_str(), accessFlag) != -1);
}

bool FileUtils::isValidPath(string path)
{
	return path.size() <= FILENAME_MAX;
}

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
	
	return StringUtils::join(PATH_DELIMITER_STR, buffer, num, pathParts);
}
