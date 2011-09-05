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

static const char* rcsid = "$Id: BaseAlignment.C,v 1.3 2011-09-05 21:23:26 mkotelbajcvi Exp $";

#include "BaseAlignment.h"

using namespace ReadAnalysis;

BaseAlignment::BaseAlignment(size_t position, size_t errorTypeReserveSize)
{
	this->position = position;
	
	this->initErrorType(MISMATCH, errorTypeReserveSize);
	this->initErrorType(INSERTION, errorTypeReserveSize);
	this->initErrorType(DELETION, errorTypeReserveSize);
}

void BaseAlignment::addError(AlignmentError* error)
{
	this->errors[error->getType()]->push_back(error);
}

void BaseAlignment::initErrorType(AlignmentErrorType type, size_t reserveSize)
{
	this->errors[type] = new vector<AlignmentError*>();
	this->errors[type]->reserve(reserveSize);
}
