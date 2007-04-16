
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Instititue.
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
static char CM_ID[]= "$Id: AS_MSG_pmesg2.c,v 1.2 2007-04-16 15:35:41 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>

#include "AS_MSG_pmesg_internal.h"
#include "AS_PER_gkpStore.h"

static
void *
Read_LIB_Mesg(FILE *fin) {
}

static
void
Write_LIB_Mesg(FILE *fout,void *mesg) {
}

static
void
Clear_LIB_Mesg(void *mesg,int typ) {
}



static
void *
Read_Frag_Mesg(FILE *fin,int frag_class) {
  static FragMesg fmesg;
  char   ch;

  assert((frag_class == MESG_FRG) || (frag_class == MESG_IFG) || (frag_class == MESG_OFG));

  fmesg.version        = 2;

  fmesg.library_uid    = 0;
  fmesg.plate_uid      = 0;
  fmesg.plate_location = 0;
  fmesg.is_random      = 1;
  fmesg.status_code    = 'G';
  fmesg.clear_vec.bgn  = 0;
  fmesg.clear_vec.end  = 0;
  fmesg.clear_qlt.bgn  = 0;
  fmesg.clear_qlt.end  = 0;
  
  if (frag_class == MESG_FRG) { 
    GET_TYPE(ch,"act:%c","action");
    fmesg.action = (ActionType) ch;
    GET_FIELD(fmesg.eaccession,"acc:"F_UID,"accession field");
  } else {
    GET_TYPE(ch,"act:%1[AD]","action");
    fmesg.action = (ActionType) ch;
    GET_PAIR(fmesg.eaccession,fmesg.iaccession,"acc:("F_UID","F_IID")","accession field pair");
  } 

  GET_FIELD(fmesg.is_random,"rnd:"F_U32,"is_random field");
  GET_TYPE(fmesg.status_code,"%1[GBUWXVEIR]","status code");

  if (frag_class == MESG_FRG) { 
    GET_FIELD(fmesg.library_uid,"lib:"F_UID,"library accession field");
  } else {
    GET_PAIR(fmesg.library_uid,fmesg.library_iid,"lib:("F_UID","F_IID")","library accession field pair");
  } 

  GET_FIELD(fmesg.plate_uid,"pla:"F_UID,"plate_uid field");
  GET_FIELD(fmesg.plate_location,"loc:"F_U32,"plate_location field");

  fmesg.source   = NULL;
  fmesg.sequence = NULL;
  fmesg.quality  = NULL;
  fmesg.hps      = NULL;

  if (fmesg.action == AS_ADD) { 
    fmesg.source   = (char *) GetText("src:",fin,FALSE);

    GET_FIELD(fmesg.entry_time,"etm:"F_TIME_T,"time field");

    if( frag_class != MESG_OFG ) {
      fmesg.sequence = (char *) GetText("seq:",fin,TRUE);
      fmesg.quality  = (char *) GetText("qlt:",fin,TRUE);
      fmesg.hps      = (char *) GetText("hps:",fin,TRUE);
    }

    GET_PAIR(fmesg.clear_rng.bgn,fmesg.clear_rng.end,"clr:"F_COORD","F_COORD,"clear range field");
    GET_PAIR(fmesg.clear_vec.bgn,fmesg.clear_rng.end,"clv:"F_COORD","F_COORD,"vector clear range field");
    GET_PAIR(fmesg.clear_qlt.bgn,fmesg.clear_rng.end,"clq:"F_COORD","F_COORD,"quality clear range field");

    fmesg.source   = AS_MSG_globals->MemBuffer + ((long) (fmesg.source));
    if( frag_class != MESG_OFG ) {
      // Convert from an index to a pointer.
      fmesg.sequence = AS_MSG_globals->MemBuffer + ((long) (fmesg.sequence));
      fmesg.quality  = AS_MSG_globals->MemBuffer + ((long) (fmesg.quality));
      fmesg.hps      = AS_MSG_globals->MemBuffer + ((long) (fmesg.hps));
    }
  }  //  action is AS_ADD
  GET_EOM;
  return ((void *) (&fmesg));
}

static void *Read_FRG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,MESG_FRG); }

static void *Read_IFG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,MESG_IFG); }

static void *Read_OFG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,MESG_OFG); }


static
void
Write_Frag_Mesg(FILE *fout,void *vmesg,int frag_class) {
  FragMesg *mesg = (FragMesg *) vmesg;

  assert((frag_class == MESG_FRG) || (frag_class == MESG_IFG) || (frag_class == MESG_OFG));

  fprintf(fout,"{%s\n",MessageTypeName[frag_class]);
  fprintf(fout,"act:%c\n",mesg->action);
  if (frag_class == MESG_FRG)
    fprintf(fout,"acc:"F_UID"\n",mesg->eaccession);
  else
    fprintf(fout,"acc:("F_UID","F_IID")\n",mesg->eaccession,mesg->iaccession);

  fprintf(fout,"rnd:%d\n",mesg->is_random);
  fprintf(fout,"sta:%c\n",mesg->status_code);

  if (frag_class == MESG_FRG)
    fprintf(fout,"lib:"F_UID"\n",mesg->library_uid);
  else
    fprintf(fout,"lib:"F_UID","F_IID"\n",mesg->library_uid,mesg->library_iid);

  fprintf(fout,"pla:"F_UID"\n",mesg->plate_uid);
  fprintf(fout,"loc:"F_U32"\n",mesg->plate_location);

  if (mesg->action == AS_ADD) {
    PutText(fout,"src:",mesg->source,FALSE);
    fprintf(fout,"etm:"F_TIME_T"\n",mesg->entry_time);
    if( frag_class != MESG_OFG ) {
      PutText(fout,"seq:",mesg->sequence,TRUE);
      PutText(fout,"qlt:",mesg->quality,TRUE);
      PutText(fout,"hps:",mesg->hps,TRUE);
    }
    fprintf(fout,"clr:"F_COORD","F_COORD"\n",mesg->clear_rng.bgn,mesg->clear_rng.end);
    fprintf(fout,"clv:"F_COORD","F_COORD"\n",mesg->clear_vec.bgn,mesg->clear_vec.end);
    fprintf(fout,"clq:"F_COORD","F_COORD"\n",mesg->clear_qlt.bgn,mesg->clear_qlt.end);
  }

  fprintf(fout,"}\n");
}

static void Write_FRG_Mesg(FILE *fout,void *mesg)
{ Write_Frag_Mesg(fout,mesg,MESG_FRG); }

static void Write_IFG_Mesg(FILE *fout,void *mesg)
{ Write_Frag_Mesg(fout,mesg,MESG_IFG); }

static void Write_OFG_Mesg(FILE *fout,void *mesg)
{ Write_Frag_Mesg(fout,mesg,MESG_OFG); }




static
void *
Read_LKG_Mesg(FILE *fin) {
  static LinkMesg lmesg;

  char ch; 
  GET_TYPE(ch,"act:%c","action");
  lmesg.action = (ActionType) ch;
  lmesg.type = AS_MATE;
  lmesg.entry_time = 0;
  lmesg.link_orient = AS_READ_ORIENT_UNKNOWN;
  GET_FIELD(lmesg.frag1,"fg1:"F_UID,"fragment 1 field");
  GET_FIELD(lmesg.frag2,"fg2:"F_UID,"fragment 2 field");
  lmesg.distance = 0;
  GET_EOM;
  return (&lmesg);
}

static
void
Write_LKG_Mesg(FILE *fout,void *mesg) {
  LinkMesg *lmesg = (LinkMesg *)mesg;
  fprintf(fout,"{LKG\n");
  fprintf(fout,"act:%c\n",lmesg->action);
  fprintf(fout,"fg1:"F_UID"\n",lmesg->frag1);
  fprintf(fout,"fg2:"F_UID"\n",lmesg->frag2);
  fprintf(fout,"}\n");
}




void AS_MSG_setFormatVersion2(void) {
  AS_MSG_callrecord  *ct = AS_MSG_globals->CallTable;

  //  Set us back to format 1
  AS_MSG_setFormatVersion1();

  //  VersionMesg VER doesn't change from format 1.
  //  AuditRecord and AuditLine don't change from format 1.

  //  LinkMesg LKG and DistanceMesg DST do not exist in version 2.
  //  The internal versions do exist.
  //
  ct[MESG_LKG].header  = "";
  ct[MESG_LKG].reader  = NULL;
  ct[MESG_LKG].writer  = NULL;
  ct[MESG_LKG].clearer = NULL;
  ct[MESG_LKG].size    = 0l;

  ct[MESG_DST].header  = "";
  ct[MESG_DST].reader  = NULL;
  ct[MESG_DST].writer  = NULL;
  ct[MESG_DST].clearer = NULL;
  ct[MESG_DST].size    = 0l;

  //  The LibraryMesg LIB is new.
  //
  ct[MESG_LIB].header  = "{LIB";
  ct[MESG_LIB].reader  = Read_LIB_Mesg;
  ct[MESG_LIB].writer  = Write_LIB_Mesg;
  ct[MESG_LIB].clearer = Clear_LIB_Mesg;
  ct[MESG_LIB].size    = sizeof(LibraryMesg);

  //  The FragMesg FRG, IFG and OFG messages are updated.

  ct[MESG_FRG].reader  = Read_FRG_Mesg;
  ct[MESG_FRG].writer  = Write_FRG_Mesg;

  ct[MESG_IFG].reader  = Read_IFG_Mesg;
  ct[MESG_IFG].writer  = Write_IFG_Mesg;

  ct[MESG_OFG].reader  = Read_OFG_Mesg;
  ct[MESG_OFG].writer  = Write_OFG_Mesg;

  //  The LinkMesg LKG is updated.

  ct[MESG_LKG].reader  = Read_LKG_Mesg;
  ct[MESG_LKG].writer  = Write_LKG_Mesg;

}

