
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
/* $Id: post_analysis.c,v 1.2 2004-09-23 20:25:21 mcschatz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "AS_UTL_ID_store.h"


int main(int argc, char *argv[])
{ GenericMesg *pmesg;
 IntScaffoldMesg *isf;
 MesgReader   reader;
 IntScreenMatch *mat;
 IntContigPairs *pairs;
 IntConConMesg *contig;
 IntUnitigMesg *unitig;
 MultiAlignT *ma;
 MultiAlignT *ma1;
 MultiAlignT *ma2;
 MultiAlignStoreT *cstore = CreateMultiAlignStoreT(0);
 MultiAlignStoreT *ustore = CreateMultiAlignStoreT(0);
 IntMultiPos *frag1;
 IntMultiPos *frag2;
 int num_frag1;
 int num_frag2;
 int ma1_len;
 int ma2_len;
 int scaffold=-1;
 int num_pairs;
 int i;
 int isplaced = 1;
 FragStoreHandle frag_store;
 FragStoreHandle bactig_store;
 FILE *pcs = NULL;
 FILE *pfs = NULL;
 FILE *dcs = NULL;
 FILE *dfs = NULL;
 FILE *out = NULL;
 FILE *sublist = NULL;
 char buffer[256];
 char *sublist_file;
 char *output_file;
 ID_Arrayp  tig_iids;
 ID_Arrayp  tig_iids_found;
 cds_int64  this_id;
 int do_all = 0;
 out = fopen("post_analysis.out","w");
 reader = InputFileType_AS( stdin );
 frag_store = openFragStore(argv[1], "rb");
 if ( argc == 5 ) {
   bactig_store = openFragStore(argv[2],"rb");
   sublist_file = argv[3];
   output_file = argv[4];
 } else {
   bactig_store = NULL;
   sublist_file = argv[2];
   output_file = argv[3];
 }
 if ( sublist_file[0] == 'A' ) { do_all = 1;}
 sprintf(buffer,"%s.pcs",output_file);
 pcs = fopen(buffer,"w");
 sprintf(buffer,"%s.pfs",output_file);
 pfs = fopen(buffer,"w");
 assert(pfs && pcs );
   
 if ( !do_all ) {
   char   string[1000];
   int    num_uids;
   sublist = fopen(sublist_file,"r");
   if( sublist == NULL )
     {
       fprintf( stderr, "Failed to open list file %s for reading.\n", argv[2] );
       exit(1);
     }
   num_uids = 0;
   while( fgets( string, 1000, sublist ) )
     {
       num_uids++;
     }
   rewind( sublist );
   tig_iids = AllocateID_Array( num_uids );
   tig_iids_found = AllocateID_Array( num_uids );
   if( tig_iids == NULL || tig_iids_found == NULL ) return 1;
   for( this_id = 0; this_id < num_uids - 1; this_id++ )
     {
       fgets( string, 1000, sublist );
       AppendToID_Array( tig_iids, STR_TO_UID(string, NULL, 10), 0 );
     }
   fgets( string, 1000, sublist );
   AppendToID_Array( tig_iids, STR_TO_UID(string, NULL, 10), 1 );
  
   fclose( sublist );
 }

 while (reader(stdin,&pmesg) != EOF){
   if (pmesg->t ==MESG_ISF)  {
     contig = (IntConConMesg *) pmesg->m;
     if( do_all || (this_id = FindID_ArrayID( tig_iids, contig->iaccession)) > -1 ) {
       if ( ! do_all ) AppendToID_Array( tig_iids_found, contig->iaccession, 1 );
       ma = CreateMultiAlignTFromICM(contig, contig->iaccession,  0);
        
       if (contig->placed == AS_PLACED) {
	 CollectStats(ma, frag_store, bactig_store, pcs, pfs,READSTRUCT_LATEST);
       } 
         
       //      PrintMultiAlignT(out,ma,frag_store,0,0);
       fflush(out);
       if ( ! do_all ) {
	 WriteProtoMesg_AS(out,pmesg);
       }
       fflush(out);
       DeleteMultiAlignT(ma);
     }
   }
   if (pmesg->t ==MESG_ISF)  {
     break;
   }
 }
 fclose(pcs);
 fclose(pfs);
 exit (0);
}
