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

#ifndef STREAMCAPTURE_H
#define STREAMCAPTURE_H

static const char* rcsid_STREAMCAPTURE_H = "$Id: StreamCapture.h,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;

#include "AS_global.h"
#include "IllegalStateException.h"

namespace Utility
{
	class StreamCapture
	{
	public:
		StreamCapture();
		StreamCapture(ios& stream);
		~StreamCapture();
		
		void startCapture();
		string stopCapture();
		string getCaptured();
		
		ios* getStream()
		{
			return this->stream;
		}
		
		void setStream(ios& stream)
		{
			this->stopCapture();
			
			this->stream = &stream;
		}
		
	protected:
		ios* stream;
		streambuf* originalBuffer;
		stringbuf* captureBuffer;
		bool capturing;
	};
}

#endif
