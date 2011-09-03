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

#ifndef FILEUTILS_H
#define FILEUTILS_H

static const char* rcsid_FILEUTILS_H = "$Id: FileUtils.h,v 1.8 2011-09-03 01:29:50 mkotelbajcvi Exp $";

#include <fcntl.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "ArgumentException.h"
#include "ErrorUtils.h"
#include "IOException.h"
#include "StringUtils.h"
#include "VarUtils.h"

namespace Utility
{
	static const char PATH_DELIMITER = '/';
	static const char* PATH_DELIMITER_STR = "/";
	
	typedef struct stat Stats;
	
	typedef enum StreamMode
	{
		READ, WRITE, APPEND, READ_WRITE, WRITE_CREATE, APPEND_END
	};
	
	class FileUtils
	{
	public:
		inline static string& readLine(FILE* stream, string& buffer, bool includeNewline = false)
		{
			if (!canRead(stream))
			{
				string errorStr;
				
				throw IOException("Stream cannot be read from: " + ErrorUtils::getError(stream, errorStr));
			}
			
			char* bufferData = (char*)buffer.data();
			size_t bufferSize = buffer.capacity();
			ssize_t lineSize = getline(&bufferData, &bufferSize, stream);
			
			if (lineSize != -1)
			{
				if (!includeNewline)
				{
					while ((bufferData[lineSize - 1] == '\n') || (bufferData[lineSize - 1] == '\r'))
					{
						lineSize--;
					}
				}
				
				buffer.assign(bufferData, lineSize);
			}
			
			return buffer;
		}

		inline static void writeLine(FILE* stream, string& buffer, bool includeNewline = true)
		{
			if (!canWrite(stream))
			{
				string errorStr;
				
				throw IOException("Stream cannot be written to: " + ErrorUtils::getError(stream, errorStr));
			}
			
			fputs(buffer.c_str(), stream);
			
			if (includeNewline && !StringUtils::isNewline(buffer[buffer.size() - 1]))
			{
				fputc(NEWLINE, stream);
			}
		}
		
		inline static void close(FILE* stream, bool ignoreErrors = true)
		{
			if (isValidStream(stream))
			{
				fclose(stream);
				
				if (!ignoreErrors && ErrorUtils::hasError(stream))
				{
					string errorStr;
					
					throw IOException("Error while closing stream: " + ErrorUtils::getError(stream, errorStr));
				}
			}
		}
		
		inline static FILE* openRead(string path, bool binary = false)
		{
			return open(path, READ, binary);
		}
		
		inline static FILE* openWrite(string path, bool binary = false)
		{
			return open(path, WRITE, binary);
		}
		
		inline static FILE* open(string path, StreamMode mode, bool binary = false)
		{
			if (((mode == READ) || (mode == READ_WRITE)) && !isReadable(path))
			{
				string errorStr;
				
				throw IOException("Path is not readable (" + ErrorUtils::getError(errorStr) + "): " + path);
			}
			
			if (((mode == WRITE) || (mode == APPEND) || (mode == READ_WRITE) || (mode == WRITE_CREATE) || (mode == APPEND_END)) && 
				!isWriteable(path))
			{
				string errorStr;
				
				throw IOException("Path is not writeable (" + ErrorUtils::getError(errorStr) + "): " + path);
			}
			
			string modeStr;
			
			FILE* stream = fopen(path.c_str(), getMode(modeStr, mode, binary).c_str());
			
			if (stream == NULL)
			{
				string errorStr;
				
				throw IOException("Unable to open stream (" + ErrorUtils::getError(stream, errorStr) + "): path=" + path + ", mode=" + modeStr);
			}
			
			return stream;
		}
		
		inline static string& getMode(string& buffer, StreamMode mode, bool binary = false)
		{
			switch (mode)
			{
				case READ:
					buffer += 'r';
					break;
				
				case WRITE:
					buffer += 'w';
					break;
					
				case APPEND:
					buffer += 'a';
					break;
					
				case READ_WRITE:
					buffer += "r+";
					break;
					
				case WRITE_CREATE:
					buffer += "w+";
					break;
				
				case APPEND_END:
					buffer += "a+";
					break;
			}
			
			if (binary)
			{
				buffer += 'b';
			}
			
			return buffer;
		}
		
		inline static bool canRead(FILE* stream)
		{
			if (feof(stream))
			{
				return false;
			}
			else
			{
				int accessMode = getAccessMode(stream);
				
				return ((accessMode == O_RDONLY) || (accessMode == O_RDWR));
			}
		}

		inline static bool canWrite(FILE* stream)
		{
			if (feof(stream))
			{
				return false;
			}
			else
			{
				int accessMode = getAccessMode(stream);
				
				return ((accessMode == O_WRONLY) || (accessMode == O_RDWR));
			}
		}
		
		inline static int getAccessMode(FILE* stream)
		{
			int status = getStatus(stream);
			
			return (status != -1) ? (status & O_ACCMODE) : status;
		}
		
		inline static int getStatus(FILE* stream)
		{
			return isValidStream(stream) ? fcntl(fileno(stream), F_GETFL, 0) : 0;
		}

		inline static bool isDirectory(string path)
		{
			return isType(path, S_IFDIR);
		}

		inline static bool isFifo(string path)
		{
			return isType(path, S_IFIFO);
		}

		inline static bool isFile(string path)
		{
			return isType(path, S_IFREG);
		}

		inline static bool isLink(string path)
		{
			return isType(path, S_IFLNK);
		}

		inline static bool isSocket(string path)
		{
			return isType(path, S_IFSOCK);
		}

		inline static bool isType(string path, mode_t type)
		{
			Stats stats;
			
			return getStats(path, stats) && ((stats.st_mode & S_IFMT) == type);
		}
		
		inline static bool getStats(FILE* stream, Stats& stats)
		{
			return isValidStream(stream) && (fstat(fileno(stream), &stats) != -1);
		}

		inline static bool getStats(string path, Stats& stats)
		{
			return exists(path) && (stat(path.c_str(), &stats) != -1);
		}

		inline static bool isExecutable(string path)
		{
			return isAccessible(path, X_OK);
		}
		
		inline static bool isReadable(string path)
		{
			return isAccessible(path, R_OK);
		}
		
		inline static bool isWriteable(string path)
		{
			return isAccessible(path, W_OK);
		}
		
		inline static bool exists(string path)
		{
			return isAccessible(path, F_OK);
		}

		inline static bool isAccessible(string path, int accessFlag)
		{
			return isValidPath(path) && (access(path.c_str(), accessFlag) != -1);
		}

		inline static int isValidStream(FILE* stream)
		{
			return stream != NULL;
		}
		
		inline static bool isValidPath(string path)
		{
			return path.size() <= FILENAME_MAX;
		}
		
		static string& getPath(string& buffer, size_t num, ...);
		static string& getPath(string& buffer, size_t num, const char** pathParts);
		
	private:
		FileUtils()
		{
		}
	};
}

#endif
