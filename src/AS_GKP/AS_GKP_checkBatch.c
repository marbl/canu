
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

static char CM_ID[] = "$Id: AS_GKP_checkBatch.c,v 1.13 2007-08-24 15:29:48 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"

int Check_BatchMesg(BatchMesg          *bat_mesg){
  GateKeeperBatchRecord  gkpb = {0};

  clearGateKeeperBatchRecord(&gkpb);

  if (bat_mesg->eaccession == 0) {
    fprintf(errorFP, "# BAT Error: Batch has zero or no UID; can't add it.\n");
    return(GATEKEEPER_FAILURE);
  }
  if (getGatekeeperUIDtoIID(gkpStore, bat_mesg->eaccession, NULL) != 0) {
    fprintf(errorFP, "# BAT Error: Batch "F_UID" exists, can't add it again.\n", bat_mesg->eaccession);
    return(GATEKEEPER_FAILURE);
  }

  gkpb.batchUID       = bat_mesg->eaccession;
  strncpy(gkpb.name,    bat_mesg->name, 255);
  strncpy(gkpb.comment, bat_mesg->comment, 255);

  appendIndexStore(gkpStore->bat, &gkpb);
  setGatekeeperUIDtoIID(gkpStore, bat_mesg->eaccession, getLastElemStore(gkpStore->bat), AS_IID_BAT);

  return GATEKEEPER_SUCCESS;
}


