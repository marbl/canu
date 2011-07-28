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

#ifndef ARGUMENT_EXCEPTION_H_
#define ARGUMENT_EXCEPTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string>

#include "RuntimeException.h"

static const char* RCSID_ARGUMENT_EXCEPTION_H_ = "$Id: ArgumentException.h,v 1.2 2011-07-28 11:31:00 mkotelbajcvi Exp $";

class ArgumentException : public RuntimeException
{
public:
	ArgumentException(char* name = NULL, char* message = NULL) throw();
	
	virtual const char* what() const throw();

protected:
	char* name;
};

#endif
