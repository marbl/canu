
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
#ifndef AS_PER_ENCODESEQUENCEQUALITY_H
#define AS_PER_ENCODESEQUENCEQUALITY_H
/*************************************************************************
 Module:  AS_PER_encodeSequenceQuality

 Description:

 Functions for encoding sequence and quality values one char per seq/quality pair.

 *************************************************************************/

/* RCS Info
 * $Id: AS_PER_encodeSequenceQuality.h,v 1.4 2005-03-22 19:49:20 jason_miller Exp $
 * $Revision: 1.4 $
 *
 */

/* encodeSequenceQuality:
     Compress sequence and quality data to 1 byte per combined value;
*/
int encodeSequenceQuality
( char *encoded, char *sequence, char *quality, uint hasQuality);

/* decodeSequenceQuality:
     Extract sequence and quality data from 1 byte per combined value
*/
int decodeSequenceQuality
( char *encoded, int encodeLength, char *sequence, char *quality, 
  uint hasQuality);


#endif
