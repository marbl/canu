
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
static char CM_ID[] = "$Id: AS_CVT_uid_interface.c,v 1.4 2005-03-22 19:48:52 jason_miller Exp $";

/*********************************************************************/
// headers
/*********************************************************************/
// standard headers
#include <stdio.h>

// project headers
#include "AS_global.h"
#include "AS_CVT_uid_interface.h"

// function that sets up UID server interaction
UIDInteractorp CreateUIDInteractor( cds_uint64 block_size)
{
  UIDInteractorp ui;

  ui = (UIDInteractorp) malloc(sizeof(UIDInteractor));
  ui->currUID = 1;

  return ui;
}


void DestroyUIDInteractor( UIDInteractorp ui )
{
  if( ui )
    free( ui );
}


static int GetUIDBlock( UIDInteractorp ui )
{
  return 0;
}


// function to get a UID
int GetUID( UIDInteractorp ui, cds_uint64 * uid )
{
  *uid = ui->currUID++;
  return 0;
}
