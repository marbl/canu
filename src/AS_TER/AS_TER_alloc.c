
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
/**********************************************************************
 Module:      AS_TER
 Description: This file contains functions to allocate
              UIDs from the UID server
 Assumptions:
**********************************************************************/

static char CM_ID[] = "$Id: AS_TER_alloc.c,v 1.1.1.1 2004-04-14 13:53:43 catmandew Exp $";

#include "AS_TER_alloc.h"
#include "AS_MSG_pmesg.h"

uint64 AS_TER_uidStart = 4711;


int32 get_uids(uint64 blockSize, uint64 *interval, int32 real)
{
  interval[0] = 4711;
  interval[1] = blockSize;
  interval[2] = 4711+2*blockSize;
  interval[3] = blockSize;
  return UID_CODE_OK;
}


int32 get_next_uid(uint64 *uid, int32 real){
  *uid = AS_TER_uidStart++;
  return UID_CODE_OK;
}

/*--------------------------------------------------------------------*/
/*  MISC routines */
/*--------------------------------------------------------------------*/

void check_environment(){
  return;
}

void set_start_uid(uint64 s) {
    AS_TER_uidStart = s;
}


