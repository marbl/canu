
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
/* TestBench for Reading Frag files and assigning internal IDs */
/* Also useful for converting .urc files to .ofg files         */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_Hash.h"

#define NOT_PRESENT 0L

uint32 lookupdst(HashTable_AS *htbl_dst, Distance_ID dst){
  uint32 id = (uint32)LookupInHashTable_AS(htbl_dst, 
					   (void *)dst, sizeof(uint64));
  if(id == 0)
    fprintf(stderr,"LookupDST %s!!!  Distance %d Key %lu\n",
	  (id?"Succeeded":"Failed"), id, dst);
  return id;

}

uint32 lookupfrg(HashTable_AS *htbl_frg, Fragment_ID frg){
  uint32 id = (uint32)LookupInHashTable_AS(htbl_frg, (void *)frg, sizeof(uint64));
  if(id == 0)
    fprintf(stderr,"LookupFrag %s!!!  Fragment %d Key %lu\n",
	  (id?"Succeeded":"Failed"), id, frg);

  return id;

}

main(int argc, char **argv){
  int ndst=1,nsfr=1,novl=1, nsfg = 1;
  enum Mesg_types imesgtype;
  GenericMesg   *pmesg;
  InternalDistMesg idt_mesg;
  DistanceMesg  *dst_mesg;
  OverlapMesg   *ovl_mesg;
  OFGMesg  ofg_mesg;
  OFragGuide ofg_guide;
  FragMesg   *frg_mesg;
  ScreenedFragMesg *sfg_mesg;
  HashTable_AS *htbl_frg;
  HashTable_AS *htbl_dst;

  FILE *fin,*fout;
  
  if(argc != 2){
    fprintf(stderr,"Usage: %s <name>\nReads <name>.urc, Outputs <name>.ofg\n\n", argv[0]);
    exit(1);
  }


  htbl_frg = CreateHashTable_uint64_AS(200000);
  htbl_dst = CreateHashTable_uint64_AS(1000);

  {
    char FileName[256];
    sprintf(FileName,"%s.urc", argv[1]);
    fin  = fopen(FileName, "r");
    if(!fin){
      fprintf(stderr,"* Couldn't open %s\n", FileName);
      exit(1);
    }
    sprintf(FileName,"%s.ofg", argv[1]);
    fout = fopen(FileName,"w");
    if(!fout){
      fprintf(stderr,"* Couldn't open %s\n", FileName);
      exit(1);
    }
  }
  while( EOF != read_protomesg(fin, &pmesg)) {

    imesgtype = pmesg->t;
#ifdef DEBUG    
    if(imesgtype != EOF) {
      write_protomesg(fout,pmesg);
      fflush(fout);
    }
#endif
    /* mesg = dup_protomesg(pmesg); */

    if( imesgtype == MESG_DST) {
      int insert;
      /*  Distance record--skip for now */
      dst_mesg = pmesg->m;
      fprintf(stderr,"* Distance Message %d\n", ndst);
      sfg_mesg = pmesg->m;
#ifdef DEBUG
      fprintf(stderr,"Read DST message with accession %ld\n", dst_mesg->accession);
#endif
      insert = InsertInHashTable_AS(htbl_dst, (void *)dst_mesg->accession, 
				    sizeof(uint64), (void *)ndst);
      if(insert == 0)
	fprintf(stderr,"Insert %s!!!  Distance %d Key %lu\n",
		(insert?"Succeeded":"Failed"), ndst, dst_mesg->accession);

      ndst++;
    }
    if( imesgtype == MESG_SFG) {
      /* Put the record where it belongs in the array.
	 This array is indexed by the overlaps. */
      int insert;
      sfg_mesg = pmesg->m;
#ifdef DEBUG
      fprintf(stderr,"Read message with accession %ld\n", sfg_mesg->accession);
#endif
      insert = InsertInHashTable_AS(htbl_frg, (void *)sfg_mesg->accession, 
				    sizeof(uint64), (void *)nsfg);
      if(insert == 0)
	fprintf(stderr,"Insert %s!!!  Fragment %d Key %lu\n",
	      (insert?"Succeeded":"Failed"), nsfg, sfg_mesg->accession);

      nsfg++;
    }

    if( imesgtype == MESG_FRG) {
      /* Put the record where it belongs in the array.
	 This array is indexed by the overlaps. */
      int insert;
      frg_mesg = pmesg->m;
#ifdef DEBUG
      fprintf(stderr,"Read message with accession %ld\n", frg_mesg->accession);
#endif
    }

    if( imesgtype == MESG_OVL) {
      ovl_mesg = pmesg->m;
      novl++;
    }
  } 
  rewind(fin);
  /**************** Phase 2 *******************/


  while( EOF != read_protomesg(fin, &pmesg)) {

    imesgtype = pmesg->t;
#ifdef DEBUG
    if(imesgtype != EOF) {
      write_protomesg(fout,pmesg);
      fflush(fout);
    }
#endif
    /* mesg = dup_protomesg(pmesg); */

    if( imesgtype == MESG_DST) {
      /*  Distance record--skip for now */
      uint64 id;

      dst_mesg = pmesg->m;
      idt_mesg.action = dst_mesg->action;
      idt_mesg.eaccession = dst_mesg->accession;
      idt_mesg.iaccession = lookupdst(htbl_dst, dst_mesg->accession);
      idt_mesg.median = dst_mesg->median;
      idt_mesg.delta = dst_mesg->delta;
      fprintf(stderr,"* Distance Message %d\n", dst_mesg->accession);

      pmesg->m = &idt_mesg;
      pmesg->t = MESG_IDT;
      write_protomesg(fout,pmesg);      
    }
    if( imesgtype == MESG_FRG) {
      /* Put the record where it belongs in the array.
	 This array is indexed by the overlaps. */
      int insert;
      sfg_mesg = pmesg->m;
      assert(0);
      nsfr++;
    }
    if( imesgtype == MESG_SFG) {
      /* Put the record where it belongs in the array.
	 This array is indexed by the overlaps. */
      int insert;
      uint32 id;
      sfg_mesg = pmesg->m;
      ofg_mesg.action = sfg_mesg->action;
      ofg_mesg.eaccession = sfg_mesg->accession;
      ofg_mesg.type = sfg_mesg->type;
      ofg_mesg.source = sfg_mesg->source;
      ofg_mesg.screened = sfg_mesg->screened;
      ofg_mesg.entry_time = sfg_mesg->entry_time;
      ofg_mesg.clear_rng = sfg_mesg->clear_rng;
      switch(ofg_mesg.type){
	case AS_READ:
	case AS_B_READ:
	case AS_EXTR:
	case AS_TRNR:
	    ofg_mesg.u.read.emate  = sfg_mesg->u.read.mate;
	    ofg_mesg.u.read.edistance  = sfg_mesg->u.read.distance;
	  if(sfg_mesg->u.read.mate != NOT_PRESENT){
	    ofg_mesg.u.read.imate  = lookupfrg(htbl_frg,sfg_mesg->u.read.mate);
	    ofg_mesg.u.read.idistance  = lookupdst(htbl_dst, sfg_mesg->u.read.distance);
	  }else{
	    ofg_mesg.u.read.imate  = NOT_PRESENT;
	    ofg_mesg.u.read.idistance  = NOT_PRESENT;
	  }	    
	  break;
      case AS_REREAD:
	ofg_mesg.u.reread.efrag = sfg_mesg->u.reread.frag;
	ofg_mesg.u.reread.ifrag = lookupfrg(htbl_frg, sfg_mesg->u.reread.frag);
	break;
      case AS_GUIDE:
	ofg_guide.eneighbor = sfg_mesg->u.guide->neighbor;
	ofg_guide.edistance = sfg_mesg->u.guide->distance;
	ofg_guide.ineighbor = lookupfrg(htbl_frg, sfg_mesg->u.guide->neighbor);
	ofg_guide.idistance = lookupdst(htbl_dst, sfg_mesg->u.guide->distance);
	ofg_mesg.u.guide = &ofg_guide;
	break;
      default:
	assert(0);
      }
#ifdef DEBUG
      fprintf(stderr,"Read message with accession %ld\n", sfg_mesg->accession);
#endif
      ofg_mesg.iaccession = lookupfrg(htbl_frg,sfg_mesg->accession);

      pmesg->m = &ofg_mesg;
      pmesg->t = MESG_OFG;

#ifdef DEBUG
      fprintf(stderr," mate = %u, %ld  distance = %u, %ld\n",
	      ofg_mesg.u.read.imate, ofg_mesg.u.read.emate,
	      ofg_mesg.u.read.idistance, ofg_mesg.u.read.edistance);
#endif
      write_protomesg(fout, pmesg);
    }

    if( imesgtype == MESG_OVL) {
      ovl_mesg = pmesg->m;
      novl++;
    }
  } 
  fclose(fin);


  fclose(fout);
  DeleteHashTable_AS(htbl_frg);
  DeleteHashTable_AS(htbl_dst);

  return 0;
}
