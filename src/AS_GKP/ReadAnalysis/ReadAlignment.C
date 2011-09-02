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

static const char* rcsid = "$Id: ReadAlignment.C,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "ReadAlignment.h"

using namespace ReadAnalysis;

ReadAlignment::ReadAlignment(size_t index)
{
	this->index = index;
	this->iid = 0;
	this->length = 0;
	this->alignment.first = this->alignment.second = 0;
	
	this->sequence.reserve(DEFAULT_READ_ALIGNMENT_SEQUENCE_RESERVE_SIZE);
	this->genomeSequence.reserve(DEFAULT_READ_ALIGNMENT_SEQUENCE_RESERVE_SIZE);
}
