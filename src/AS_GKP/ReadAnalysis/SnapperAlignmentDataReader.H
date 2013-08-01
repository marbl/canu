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

#ifndef SNAPPERALIGNMENTDATAREADER_H
#define SNAPPERALIGNMENTDATAREADER_H

static const char* rcsid_SNAPPERALIGNMENTDATAREADER_H = "$Id$";

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AlignmentDataReader.h"
#include "DataException.h"
#include "FileUtils.h"
#include "ReadAlignment.h"
#include "StringUtils.h"

using namespace Utility;

namespace ReadAnalysis
{
	static const char* SNAPPER_READ_ALIGNMENT_START = "sim4begin";
	static const char* SNAPPER_READ_ALIGNMENT_END = "sim4end";
	
	static const char* SNAPPER_READ_INFO_FORMAT = F_U64"["F_U64"-"F_U16I"-"F_U16I"] "F_U16I"["F_U64I"-"F_U64I"] <"F_U64I"-"F_U16I"-"F_U16I"-"F_STRI">";
	static const char* SNAPPER_SANGER_READ_DEFINITION_FORMAT = 
		"edef="F_U32","F_U32I" mate="F_U32","F_U32I" lib="F_STRI","F_U32I" clr="F_STRI","F_U64I","F_U64I" deleted="F_CI;
	static const char* SNAPPER_ILLUMINA_READ_DEFINITION_FORMAT = 
		"edef"F_STRI":"F_U16I":"F_U16I":"F_U32":"F_STRI;
	static const char* SNAPPER_GENOME_DEFINITION_FORMAT = "ddef="F_STRI;
	static const char* SNAPPER_ALIGNMENT_INFO_FORMAT = F_U64"-"F_U64" ("F_U64"-"F_U64") <"F_U64I"-"F_U16I"-"F_U16">";
	
	class SnapperAlignmentDataReader : public AlignmentDataReader
	{
	public:
		SnapperAlignmentDataReader();
		virtual ~SnapperAlignmentDataReader();
		
	protected:
		inline static string& readLine(FILE* stream, string& line, size_t& lineNum)
		{
			line = string();
			
			FileUtils::readLine(stream, line);
			
			lineNum++;
			
			return line;
		}
		
		static void processStart(string line, size_t lineNum);
		static void processReadInfo(ReadAlignment* readAlign, string line, size_t lineNum);
		static void processReadDefinition(ReadAlignment* readAlign, string line, size_t lineNum);
		static void processGenomeDefinition(ReadAlignment* readAlign, string line, size_t lineNum);
		static void processAlignmentInfo(ReadAlignment* readAlign, string line, size_t lineNum);
		static void processReadSequence(ReadAlignment* readAlign, string line, size_t lineNum);
		static void processGenomeSequence(ReadAlignment* readAlign, string line, size_t lineNum);
		static void processSequence(vector<char>& sequence, string line, size_t lineNum);
		static void processEnd(string line, size_t lineNum);
		
		virtual void processData();
	};
}

#endif
