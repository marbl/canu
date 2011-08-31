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

#ifndef ILLEGALSTATEEXCEPTION_H
#define ILLEGALSTATEEXCEPTION_H

static const char* rcsid_ILLEGALSTATEEXCEPTION_H = "$Id: IllegalStateException.h,v 1.3 2011-08-31 06:49:27 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using namespace std;

#include "AS_UTL_alloc.h"
#include "ExceptionUtils.h"
#include "RuntimeException.h"
#include "StringUtils.h"

class IllegalStateException : public RuntimeException
{
public:
	IllegalStateException(const char* message = NULL, RuntimeException* cause = NULL) throw();
	IllegalStateException(string message = string(), RuntimeException* cause = NULL) throw();
	
protected:
	virtual void initialize() throw();
};

#endif
