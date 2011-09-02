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

#ifndef ALIGNMENTERROR_H
#define ALIGNMENTERROR_H

static const char* rcsid_ALIGNMENTERROR_H = "$Id: AlignmentError.h,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>

using namespace std;

#include "AS_global.h"
#include "AlignmentErrorType.h"

namespace ReadAnalysis
{
	class AlignmentError
	{
	public:
		AlignmentError(AS_IID iid = 0, AlignmentErrorType type = UNKNOWN);
		
		AS_IID getIid()
		{
			return this->iid;
		}

		AlignmentErrorType getType()
		{
			return this->type;
		}

		void setIid(AS_IID iid)
		{
			this->iid = iid;
		}

		void setType(AlignmentErrorType type)
		{
			this->type = type;
		}
		
	protected:
		AS_IID iid;
		AlignmentErrorType type;
	};
}

#endif
