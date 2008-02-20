
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

static char const *rcsid = "$Id: AS_GKP_checkBatch.c,v 1.17 2008-02-20 10:53:30 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"

int Check_BatchMesg(BatchMesg          *bat_mesg){
  GateKeeperBatchRecord  gkpb;

  clearGateKeeperBatchRecord(&gkpb);

  gkpStore->gkp.batInput++;

  if (AS_UID_isDefined(bat_mesg->eaccession) == FALSE) {
    AS_GKP_reportError(AS_GKP_BAT_ZERO_UID);
    gkpStore->gkp.batErrors++;
    return(1);
  }

  if (getGatekeeperUIDtoIID(gkpStore, bat_mesg->eaccession, NULL) != 0) {
    AS_GKP_reportError(AS_GKP_BAT_EXISTS, AS_UID_toString(bat_mesg->eaccession));
    gkpStore->gkp.batErrors++;
    return(1);
  }

  gkpb.batchUID       = bat_mesg->eaccession;
  strncpy(gkpb.name,    bat_mesg->name, 255);
  strncpy(gkpb.comment, bat_mesg->comment, 255);

  appendIndexStore(gkpStore->bat, &gkpb);
  setGatekeeperUIDtoIID(gkpStore, bat_mesg->eaccession, getLastElemStore(gkpStore->bat), AS_IID_BAT);

  gkpStore->gkp.batLoaded++;

  return(0);
}


