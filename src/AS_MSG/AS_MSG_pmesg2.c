
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
static char CM_ID[]= "$Id: AS_MSG_pmesg2.c,v 1.4 2007-04-26 14:07:03 brianwalenz Exp $";

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
  static LibraryMesg  lmesg;  // statics are initialized to zero, important for realloc() below
  char                ch;
  long                featureoffset = 0;

  GET_TYPE(ch,"act:%c","action");
  lmesg.action = (ActionType) ch;

  GET_FIELD(lmesg.eaccession,"acc:"F_UID,"accession field");

  if (lmesg.action == AS_UPDATE) {
    GET_FIELD(lmesg.mean,"mea:%f","mean field");
    GET_FIELD(lmesg.stddev ,"std:%f","stddev field");
  }

  lmesg.num_features = 0;

  if ((lmesg.action == AS_ADD) || (lmesg.action == AS_IGNORE)) {
    GET_TYPE(ch,"ori:%c","orientation");
    lmesg.link_orient = (OrientType) ch;

    GET_FIELD(lmesg.mean,"mea:%f","mean field");
    GET_FIELD(lmesg.stddev ,"std:%f","stddev field");

    lmesg.source   = (char *)GetText("src:",fin,FALSE);

    GET_FIELD(lmesg.num_features ,"nft:%d","number of features");
    featureoffset = GetText("fea:",fin,TRUE);

    //  Adjust source to be a pointer instead of an offset
    lmesg.source = AS_MSG_globals->MemBuffer + (long)lmesg.source;
  }  //  End of AS_ADD & AS_IGNORE

  GET_EOM;


  //  Munge the feature string into a list of features
  //
  if (lmesg.num_features > 0) {
    int      i;
    char    *fb = AS_MSG_globals->MemBuffer + featureoffset;
    char    *fn = fb;

    lmesg.features = (char **)safe_realloc(lmesg.features, sizeof(char *) * lmesg.num_features);
    lmesg.values   = (char **)safe_realloc(lmesg.values,   sizeof(char *) * lmesg.num_features);

    for (i=0; i<lmesg.num_features; i++) {
      //  get rid of spaces in the label
      while (isspace(*fb))
        fb++;

      lmesg.features[i] = fb;

      //  Look for the '='
      while (*fb != '=')
        fb++;

      //  strip whitespace at the end of the label
      fn = fb-1;
      while (isspace(*fn)) {
        *fn = 0;
        fn--;
      }

      //  skip over the =, and any white space
      fb++;
      while (isspace(*fb))
        fb++;

      lmesg.values[i] = fb;
      
      //  Look for the end (a new line)
      while ((*fb != '\n') && (*fb != '\r'))
        fb++;

      //  strip whitespace at the end of the label
      fn = fb-1;
      while (isspace(*fn)) {
        *fn = 0;
        fn--;
      }

      fprintf(stderr, "GOT:   F'%s' V'%s'\n", lmesg.features[i], lmesg.values[i]);
    }
  }  //  num_features > 0


  return(&lmesg);
}

static
void
Write_LIB_Mesg(FILE *fout,void *mesg) {
  LibraryMesg  *lmesg = (LibraryMesg *)mesg;
  int           i;

  fprintf(fout,"{LIB\n");
  fprintf(fout,"act:%c\n", lmesg->action);
  fprintf(fout,"acc:"F_UID"\n", lmesg->eaccession);

  if (lmesg->action == AS_UPDATE) {
    fprintf(fout,"mea:%.3f\n",lmesg->mean);
    fprintf(fout,"std:%.3f\n",lmesg->stddev);
  }

  if ((lmesg->action == AS_ADD) || (lmesg->action == AS_IGNORE)) {
    fprintf(fout,"ori:%c\n",lmesg->link_orient);
    fprintf(fout,"mea:%.3f\n",lmesg->mean);
    fprintf(fout,"std:%.3f\n",lmesg->stddev);

    PutText(fout,"src:",lmesg->source,FALSE);

    fprintf(fout,"nft:%d\n",lmesg->num_features);
    fprintf(fout,"fea:\n");
    for (i=0; i<lmesg->num_features; i++)
      fprintf(fout,"%s=%s\n", lmesg->features[i], lmesg->values[i]);
    fprintf(fout,".\n");
    fprintf(fout,"}\n");
  }
}



static
void *
Read_Frag_Mesg(FILE *fin,int frag_class) {
  static FragMesg fmesg;
  char   ch;

  assert((frag_class == MESG_FRG) || (frag_class == MESG_IFG) || (frag_class == MESG_OFG));

  fmesg.version        = 2;

  fmesg.library_uid    = 0;
  fmesg.library_iid    = 0;
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

  GET_FIELD(fmesg.is_random,"rnd:"F_U32,"is_random");
  GET_TYPE(fmesg.status_code,"sta:%1[GBUWXVEIR]","status code");

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

  if ((fmesg.action == AS_ADD) || (fmesg.action == AS_IGNORE)) { 
    char  *line;
    int    b, e;

    fmesg.source   = (char *) GetText("src:",fin,FALSE);

    if( frag_class != MESG_OFG ) {
      fmesg.sequence = (char *) GetText("seq:",fin,TRUE);
      fmesg.quality  = (char *) GetText("qlt:",fin,TRUE);
      fmesg.hps      = (char *) GetText("hps:",fin,TRUE);
    }

    //  Special handling for clear ranges -- the vector and quality clear are optional.

    fmesg.clear_vec.bgn = 1;
    fmesg.clear_vec.end = 0;

    fmesg.clear_qlt.bgn = 1;
    fmesg.clear_qlt.end = 0;

    fmesg.clear_rng.bgn = 1;
    fmesg.clear_rng.end = 0;

    line = GetLine(fin, TRUE);
    if(sscanf(line,"clv:"F_COORD","F_COORD,&b,&e)==2){
      fmesg.clear_vec.bgn = b;
      fmesg.clear_vec.end = e;
      line = GetLine(fin, TRUE);
    }
    if(sscanf(line,"clq:"F_COORD","F_COORD,&b,&e)==2){
      fmesg.clear_qlt.bgn = b;
      fmesg.clear_qlt.end = e;
      line = GetLine(fin, TRUE);
    }
    if(sscanf(line,"clr:"F_COORD","F_COORD,&b,&e)==2){
      fmesg.clear_rng.bgn = b;
      fmesg.clear_rng.end = e;
    } else {
      MfieldError("final clear range field");	
    }

    // Convert from an index to a pointer.
    //
    fmesg.source   = AS_MSG_globals->MemBuffer + ((long) (fmesg.source));
    if( frag_class != MESG_OFG ) {
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

  if ((mesg->action == AS_ADD) || (mesg->action == AS_IGNORE)) {
    PutText(fout,"src:",mesg->source,FALSE);
    if( frag_class != MESG_OFG ) {
      PutText(fout,"seq:",mesg->sequence,TRUE);
      PutText(fout,"qlt:",mesg->quality,TRUE);
      PutText(fout,"hps:",mesg->hps,TRUE);
    }
    fprintf(fout,"clr:"F_COORD","F_COORD"\n",mesg->clear_rng.bgn,mesg->clear_rng.end);

    if (mesg->clear_vec.bgn <= mesg->clear_vec.end)
      fprintf(fout,"clv:"F_COORD","F_COORD"\n",mesg->clear_vec.bgn,mesg->clear_vec.end);

    if (mesg->clear_qlt.bgn <= mesg->clear_qlt.end)
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
  lmesg.link_orient = AS_READ_ORIENT_UNKNOWN;
  GET_FIELD(lmesg.frag1,"frg:"F_UID,"fragment 1 field");
  GET_FIELD(lmesg.frag2,"frg:"F_UID,"fragment 2 field");
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
  fprintf(fout,"frg:"F_UID"\n",lmesg->frag1);
  fprintf(fout,"frg:"F_UID"\n",lmesg->frag2);
  fprintf(fout,"}\n");
}




void AS_MSG_setFormatVersion2(void) {
  AS_MSG_callrecord  *ct = AS_MSG_globals->CallTable;

  //  Set us back to format 1
  AS_MSG_setFormatVersion1();

  //  VersionMesg VER doesn't change from format 1.
  //  AuditRecord and AuditLine don't change from format 1.

  //  DistanceMesg DST does not exist in version 2.
  //  The internal version (IDT) does exist.
  //
  ct[MESG_DST].header  = "";
  ct[MESG_DST].reader  = NULL;
  ct[MESG_DST].writer  = NULL;
  ct[MESG_DST].size    = 0l;

  //  The LibraryMesg LIB is new.
  //
  ct[MESG_LIB].header  = "{LIB";
  ct[MESG_LIB].reader  = Read_LIB_Mesg;
  ct[MESG_LIB].writer  = Write_LIB_Mesg;
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

