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

static const char* rcsid = "$Id: StreamCapture.C,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "StreamCapture.h"

using namespace Utility;

StreamCapture::StreamCapture()
{
	this->stream = NULL;
	this->originalBuffer = NULL;
	this->captureBuffer = NULL;
	this->capturing = false;
}

StreamCapture::StreamCapture(ios& stream)
{
	this->stream = &stream;
	this->originalBuffer = NULL;
	this->captureBuffer = NULL;
	this->capturing = false;
}

StreamCapture::~StreamCapture()
{
	if (this->capturing)
	{
		this->stopCapture();
	}
}

void StreamCapture::startCapture()
{
	if (this->stream == NULL)
	{
		throw IllegalStateException("Stream to capture must be set.");
	}
	
	if (!this->capturing)
	{
		this->captureBuffer = new stringbuf();
		this->originalBuffer = this->stream->rdbuf(this->captureBuffer);
		
		this->capturing = true;
	}
}

string StreamCapture::stopCapture()
{
	if (this->capturing)
	{
		this->stream->rdbuf(this->originalBuffer);
		
		this->capturing = false;
	}
	
	return this->getCaptured();
}

string StreamCapture::getCaptured()
{
	return (this->captureBuffer != NULL) ? this->captureBuffer->str() : string();
}
