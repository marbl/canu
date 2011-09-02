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

static const char* rcsid = "$Id: BaseAlignment.C,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "BaseAlignment.h"

using namespace ReadAnalysis;

BaseAlignment::BaseAlignment(size_t position, size_t readReserveSize, size_t errorTypeReserveSize)
{
	this->position = position;
	this->errorTypeReserveSize = errorTypeReserveSize;
	
	this->reads.reserve(readReserveSize);
}

void BaseAlignment::addError(AlignmentError error)
{
	this->getErrors(error.getType()).push_back(error);
}

vector<AlignmentError>& BaseAlignment::getErrors(AlignmentErrorType type)
{
	vector<AlignmentError>* typeErrors = (type != UNKNOWN) && this->hasErrors(type) ? &this->errors[type] : NULL;
	
	if (typeErrors == NULL)
	{
		typeErrors = new vector<AlignmentError>();
		
		if (type != UNKNOWN)
		{
			typeErrors->reserve(this->errorTypeReserveSize);
			
			this->errors[type] = *typeErrors;
		}
		else
		{
			typeErrors->reserve(this->getNumErrors(MISMATCH) + this->getNumErrors(INSERTION) + this->getNumErrors(DELETION));
			
			copy(this->errors[MISMATCH].begin(), this->errors[MISMATCH].end(), typeErrors->end());
			copy(this->errors[INSERTION].begin(), this->errors[INSERTION].end(), typeErrors->end());
			copy(this->errors[DELETION].begin(), this->errors[DELETION].end(), typeErrors->end());
		}
	}
	
	return *typeErrors;
}

size_t BaseAlignment::getNumErrors(AlignmentErrorType type)
{
	if (type == UNKNOWN)
	{
		size_t numErrors = 0;
		
		if (this->hasErrors(MISMATCH))
		{
			numErrors += this->errors[MISMATCH].size();
		}
		
		if (this->hasErrors(INSERTION))
		{
			numErrors += this->errors[INSERTION].size();
		}
		
		if (this->hasErrors(DELETION))
		{
			numErrors += this->errors[DELETION].size();
		}
		
		return numErrors;
	}
	else
	{
		return (this->hasErrors(type) ? this->errors[type].size() : 0);
	}
}

bool BaseAlignment::hasErrors(AlignmentErrorType type)
{
	return !this->errors.empty() && ((type == UNKNOWN) || (this->errors.count(type) != 0) && !this->errors[type].empty());
}
