
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
static char CM_ID[] 
= "$Id: AS_CGB_util.c,v 1.5 2007-05-14 09:27:10 brianwalenz Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_util.c
 * 
 * Utility functions for the CGB (Unitigger)
 * subsystem of the Celera WGS genome assembler.
 *
 * Author: Jason Miller
 *********************************************************************/

/*********************************************************************/
/* System include files */
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/*********************************************************************/

/*********************************************************************/
/* Local include files */
#include "AS_CGB_util.h"
/************************************************************************/

/* ---------------------------------------- */
// Safely get the length of a Proto-IO Message delta encoding.
// Later, define DELTA as a class that reports length.
// Assumes DELTA never encodes the NUL character.
// Move this function to the AS_MSG subsystem.
size_t getDeltaLength (signed char * deltaCode) {
    char * uDelta;
    size_t deltaLen;
    uDelta = (char *) deltaCode;
    deltaLen = strlen(uDelta);
    return deltaLen;
}

/* ---------------------------------------- */
// Safely cast these opaque types.
// Note Tnes is specific to the CGB subsystem.
// Later, change assert() to error reporting.
Tnes AS_CGB_SafeCast_cdsInt8_to_Tnes (int8 dataIn) {
    Tnes dataOut;
    dataOut = (Tnes) dataIn;
    assert (dataIn == (int8) dataOut);
    return dataOut;
}
