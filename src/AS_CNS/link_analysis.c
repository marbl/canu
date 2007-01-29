
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
/* $Id: link_analysis.c,v 1.7 2007-01-29 20:41:08 brianwalenz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA_MSG.h"
#include "MultiAlignStore_CNS.h"
#include "Statistics_CNS.h"
#include "AS_UTL_ID_store.h"
#include "AS_PER_gkpStore.h"

int main(int argc, char *argv[])
{ GenericMesg *pmesg;
  IntScaffoldMesg *isf;
  IntScreenMatch *mat;
  IntContigPairs *pairs;
  IntConConMesg *contig;
  IntUnitigMesg *unitig;
  IntMultiPos *frag1;
  IntMultiPos *frag2;
  ReadStructp rsp = new_ReadStruct();
  int num_frag1;
  int num_frag2;
  int ma1_len;
  int ma2_len;
  VA_TYPE(int) *placed;
  int scaffold=-1;
  int num_pairs;
  int i;
  int isplaced = 1;
  FragStoreHandle frag_store;
  GateKeeperStore gkp_store;
  GateKeeperLinkStore gkpl_store;
  GateKeeperLinkRecord link;
  GateKeeperRecord gkpf;
      char   string[1000];
  ID_Arrayp  tig_iids;
  ID_Arrayp  tig_iids_found;
  cds_int64  this_id;
  IntFragment_ID frag_id;
  IntFragment_ID linked_id;
  int    num_uids;
  Fragment_ID accf;
  Fragment_ID accl;
  FILE *iumfile;
  int in_unitig;
  iumfile = fopen(argv[4],"r");
  ReadProtoMesg_AS(iumfile,&pmesg);
  unitig = pmesg->m;
  tig_iids = AllocateID_Array( unitig->num_frags );
  tig_iids_found = AllocateID_Array( unitig->num_frags );
  for( i = 0; i < unitig->num_frags - 1; i++ ) {
      AppendToID_Array( tig_iids, (cds_uint64) unitig->f_list[i].ident, 0 );
  }
  AppendToID_Array( tig_iids, (cds_uint64) unitig->f_list[unitig->num_frags-1].ident, 1 );
  frag_store = openFragStore(argv[1], "rb");
  gkp_store = openGateKeeperStore(argv[2],"r");
  gkpl_store = openGateKeeperLinkStore(argv[3],"r");
  GateKeeperLinkRecordIterator iterator;

  for (i = 0; i < unitig->num_frags;i ++ ) {
    frag_id = unitig->f_list[i].ident;
    getGateKeeperStore(gkp_store, frag_id, &gkpf);
    if ( gkpf.linkHead ) {    
      CreateGateKeeperLinkRecordIterator(gkpl_store, gkpf.linkHead, frag_id, &iterator);
      while ( NextGateKeeperLinkRecordIterator(&iterator, &link) ) {
        if ( link.frag1 == frag_id )  { 
          linked_id = link.frag2;
        } else if ( link.frag2 == frag_id )  {
          linked_id = link.frag1;
        } else {
          fprintf(stderr,"Link didn't make sense. id: %d\n",frag_id);
          exit(1);
       }
       getFragStore(frag_store,frag_id,FRAG_S_ALL,rsp);
       getAccID_ReadStruct(rsp, &accf);
       getFragStore(frag_store,linked_id,FRAG_S_ALL,rsp);
       getAccID_ReadStruct(rsp, &accl);
       if ( FindID_ArrayID( tig_iids, linked_id) > -1 ) {
          in_unitig = 1;
       } else {
          in_unitig = 0;
       }
       fprintf(stderr,"%d %lu %d %lu distance: %lu in_unitig: %d\n", frag_id, accf, linked_id, accl,link.distance,in_unitig);
    }
   }
 }
 exit (0);
}
