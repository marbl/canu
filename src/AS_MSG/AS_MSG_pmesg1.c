
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
static char CM_ID[]= "$Id: AS_MSG_pmesg1.c,v 1.7 2007-04-28 08:46:22 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>

#include "AS_MSG_pmesg_internal.h"
#include "AS_PER_gkpStore.h"

//  reads old and new AFG message (with and w/o chaff field)
#define AFG_BACKWARDS_COMPATIBLE
#define IAF_BACKWARDS_COMPATIBLE


/******************** INPUT ROUTINES ***************************/

static void *Read_Dist_Mesg(FILE *fin, int external)
{ static InternalDistMesg dmesg;

  if (external)
    { char ch;
      GET_TYPE(ch,"act:%c","action");
      dmesg.action = (ActionType) ch;
      GET_FIELD(dmesg.eaccession,"acc:"F_UID,"accession field")
        }
  else
    { char ch;
      GET_TYPE(ch,"act:%1[ADIUR]","action");  //  R, redefine, for compatibility; same as update
      dmesg.action = (ActionType) ch;
      GET_PAIR(dmesg.eaccession,dmesg.iaccession,"acc:("F_UID","F_IID")","accession field pair");
    }
  if (dmesg.action == 'R')
    dmesg.action = AS_UPDATE;
  if ((dmesg.action == AS_ADD) || (dmesg.action == AS_UPDATE) || (dmesg.action == AS_IGNORE))
    { GET_FIELD(dmesg.mean,"mea:%f","mean field");
      GET_FIELD(dmesg.stddev ,"std:%f","stddev field");
    }
  GET_EOM;
  return ((void *) (&dmesg));
}

static void *Read_DST_Mesg(FILE *fin)
{ return Read_Dist_Mesg(fin,1); }

static void *Read_IDT_Mesg(FILE *fin)
{ return Read_Dist_Mesg(fin,0); }

// This is tricky, since we are trying to update the fields
// of a structure that may be realloced.  Store offsets in
// temps, and stuff them in in one go.

static int Read_ADL_Struct(long last, FILE *fin)
{ AuditLine mesg; // temporary
  long crnt;
  AuditLine *newMesg;

  mesg.name    = (char *)GetString("who:",fin);
  GET_FIELD(mesg.complete,"ctm:"F_TIME_T,"completion field");
  mesg.version = (char *) GetString("vsn:",fin);
  mesg.comment = (char *)GetText("com:",fin, FALSE);
  mesg.next = (AuditLine *)last;
  // this looks SUSPICIOUS
  GET_EOM;

  // Now allocate space for this guy, and copy him to the allocated space
  crnt = MoreSpace(sizeof(AuditLine),8);
  newMesg = (AuditLine *)(AS_MSG_globals->MemBuffer + crnt);
  *newMesg = mesg;

  return crnt;
}

static void *Read_ADT_Mesg(FILE *fin)
{ static AuditMesg amesg;
  AuditLine *cptr, *tail;
  long       last,  crnt;

  /* First build up list (in reverse order) using indices and not pointers */

  last = crnt = -1;
  while (strncmp(GetLine(fin,TRUE),"{ADL",4) == 0)
    { 
      last =  Read_ADL_Struct(last,fin);
    }
  if (AS_MSG_globals->CurLine[0] != '.')
    MgenError("Expecting end of ADL list");
  GET_EOM;

  /* Traverse again, reversing list order and converting indices to ptrs. */

  //  fprintf(stderr,"* crnt = %lx last = %lx \n", crnt, last);

  tail = NULL;

  crnt = last;

  while (crnt >= 0)
    { cptr = (AuditLine *) (AS_MSG_globals->MemBuffer + crnt);
      crnt = (long) (cptr->next);
      cptr->next    = tail;
      cptr->name    = AS_MSG_globals->MemBuffer + ((long) (cptr->name));
      cptr->version = AS_MSG_globals->MemBuffer + ((long) (cptr->version));
      cptr->comment = AS_MSG_globals->MemBuffer + ((long) (cptr->comment));
      
      tail = cptr;
    }
  amesg.list = tail;

  return (&amesg);
}

static void *Read_VER_Mesg(FILE *fin) {
  static VersionMesg vmesg;
  GET_FIELD(vmesg.version,"ver:"F_U32,"version field");
  GET_EOM;
  switch (vmesg.version) {
    case 1:
      AS_MSG_setFormatVersion1();
      break;
    case 2:
      AS_MSG_setFormatVersion2();
      break;
    default:
      fprintf(stderr,"ERROR: Unknown version "F_U32".\n", vmesg.version);
      assert((vmesg.version == 1) ||
             (vmesg.version == 2));
      break;
  }
  return(&vmesg);
}


static void *Read_Frag_Mesg(FILE *fin, int frag_class)
{ static FragMesg fmesg;
  char   ch;
  time_t entry_time;

  assert((frag_class == MESG_FRG) || (frag_class == MESG_IFG) || (frag_class == MESG_OFG));

  fmesg.version        = 1;
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

  fmesg.source   = NULL;
  fmesg.sequence = NULL;
  fmesg.quality  = NULL;
  fmesg.hps      = NULL;

  if ((fmesg.action == AS_ADD) || (fmesg.action == AS_IGNORE)) { 
    if (frag_class == MESG_FRG){
      GET_TYPE(ch,"typ:%c","type");
      fmesg.type = (FragType) ch;
      // We want to succeed on all reads, and let the gatekeeper do its stuff
    }else{
      GET_TYPE(ch,"typ:%1[RXELTFUCBW]","type");
      fmesg.type = (FragType) ch;
    }

    fmesg.source   = (char *) GetText("src:",fin,FALSE);

    //  Unused
    GET_FIELD(entry_time,"etm:"F_TIME_T,"time field");

    if( frag_class != MESG_OFG ) {
      fmesg.sequence = (char *) GetText("seq:",fin,TRUE);
      fmesg.quality  = (char *) GetText("qlt:",fin,TRUE);
    }
    GET_PAIR(fmesg.clear_rng.bgn,fmesg.clear_rng.end,"clr:"F_COORD","F_COORD,"clear range field");
    fmesg.source   = AS_MSG_globals->MemBuffer + ((long) (fmesg.source));
    if( frag_class != MESG_OFG ) {
      // Convert from an index to a pointer.
      fmesg.sequence = AS_MSG_globals->MemBuffer + ((long) (fmesg.sequence));
      fmesg.quality  = AS_MSG_globals->MemBuffer + ((long) (fmesg.quality));
    }
  }  //  action is AS_ADD or AS_IGNORE
  GET_EOM;
  return ((void *) (&fmesg));
}

static void *Read_FRG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,MESG_FRG); } 

static void *Read_IFG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,MESG_IFG); }

static void *Read_OFG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,MESG_OFG); }


static void *Read_OVL_Mesg(FILE *fin)
{ static OverlapMesg omesg;
  int    idx;
  char ch;
  
  GET_FIELD(omesg.aifrag,"afr:"F_IID,"a-fragment field");
  GET_FIELD(omesg.bifrag,"bfr:"F_IID,"b-fragment field");
  GET_TYPE(ch,"ori:%1[NAIO]","orientation");
  omesg.orientation = (OrientType) ch;
  GET_TYPE(ch,"olt:%1[DCSXdc]","overlap");
  omesg.overlap_type = (OverlapType) ch;
  GET_FIELD(omesg.ahg,"ahg:"F_COORD,"a-hang field");
  GET_FIELD(omesg.bhg,"bhg:"F_COORD,"b-hang field");
  GET_FIELD(omesg.quality,"qua:%f","quality field");
  GET_FIELD(omesg.min_offset,"mno:"F_COORD,"min-offset field");
  GET_FIELD(omesg.max_offset,"mxo:"F_COORD,"max-offset field");
  GET_FIELD(omesg.polymorph_ct,"pct:"F_S32,"poly-count field");

  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("delta tag label");
  idx = MoreSpace(2*AS_FRAG_MAX_LEN,1);
  omesg.delta = (signed char *) (AS_MSG_globals->MemBuffer + idx);

  { int i, n;     /* Read a delta item (only one of its kind) */
    char *t, *u;

    i = 0;
    while ((t = GetLine(fin,TRUE))[0] != '.')
      while (1)
        { n = strtol(t,&u,10);
          if (u == t) break;
          t = u;
          if (! isspace((int)*t))
            MgenError("Delta is not a sequence of digits");
          omesg.delta[i++] = n;
        }
    omesg.delta[i] = 0;
  }

  GET_EOM;
  return ((void *) (&omesg));
}

static void *Read_LKG_Mesg(FILE *fin)
{ static LinkMesg lmesg;
  char ch;
  time_t entry_time;
  GET_TYPE(ch,"act:%c","action");
  lmesg.action = (ActionType) ch;
  GET_TYPE(ch,"typ:%c","link");
  lmesg.type = (LinkType) ch;
  GET_FIELD(lmesg.frag1,"fg1:"F_UID,"fragment 1 field");
  GET_FIELD(lmesg.frag2,"fg2:"F_UID,"fragment 2 field");
  if ((lmesg.action == AS_ADD) || (lmesg.action == AS_IGNORE))
    {
      //  Unused
      GET_FIELD(entry_time,"etm:"F_TIME_T,"entry time field");
      GET_FIELD(lmesg.distance,"dst:"F_UID,"distance field");
      GET_TYPE(ch,"ori:%1[NAIOU]","link orientation");
      lmesg.link_orient = (OrientType) ch;
    }
  GET_EOM;
  return ((void *) (&lmesg));
}

static void *Read_UOM_Mesg(FILE *fin)
{ static UnitigOverlapMesg	mesg;
  char ch; 
  GET_FIELD(mesg.chunk1,"ck1:"F_IID,"chunk 1 id field");
  GET_FIELD(mesg.chunk2,"ck2:"F_IID,"chunk 2 id field");
  GET_TYPE(ch,"ori:%1[NAIO]","orientation");
  mesg.orient = (ChunkOrientationType) ch;
  GET_TYPE(ch,"ovt:%1[NOTCIMXdcYZ]","overlap type");
  mesg.overlap_type = (UnitigOverlapType) ch;

#ifdef AS_ENABLE_SOURCE
  {
    long		sindx;
    sindx = GetText("src:",fin,FALSE);
    mesg.source = AS_MSG_globals->MemBuffer + sindx;
  }
#endif
  GET_FIELD(mesg.best_overlap_length,"len:"F_COORD,"best overlap");
  GET_FIELD(mesg.min_overlap_length,"min:"F_COORD,"min overlap");
  GET_FIELD(mesg.max_overlap_length,"max:"F_COORD,"max overlap");
  GET_FIELD(mesg.quality,"qua:%f","quality field");
  GET_EOM;
  return ((void *) (&mesg));
}

static void Read_IMP_Mesg(FILE *fin, long indx)
{
  IntMultiPos	*imp;
  int		 i;
  int32		 n, *delta;
  long		 tindx;
  char		*line, *u;
  char ch;
  
  imp = (IntMultiPos *) (AS_MSG_globals->MemBuffer + indx);
  GET_TYPE(ch,"typ:%1[RXTELUFSUcBCG]","multipos$");
  imp->type = (FragType) ch;
  GET_FIELD(imp->ident,"mid:"F_IID,"multipos id");
  GET_FIELD(imp->contained,"con:"F_IID,"contained id");
#ifdef NEW_UNITIGGER_INTERFACE
  GET_FIELD(imp->ident2,"bid:"F_IID,"multipos id");
#endif
#ifdef AS_ENABLE_SOURCE
  imp->sourceInt = -1;
#endif
  GET_PAIR(imp->position.bgn,imp->position.end,"pos:"F_COORD","F_COORD,"position field");
#ifdef NEW_UNITIGGER_INTERFACE
  GET_FIELD(imp->ahang,"ahg:"F_S32,"ahang");
  GET_FIELD(imp->bhang,"bhg:"F_S32,"bhang");
#endif
  GET_FIELD(imp->delta_length,"dln:"F_S32,"delta length");
  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  if (imp->delta_length > 0) {
    tindx = MoreSpace(sizeof(int32)*imp->delta_length,8);
    imp = (IntMultiPos *) (AS_MSG_globals->MemBuffer + indx);	// in case of realloc
    imp->delta = (int32 *) tindx;
    delta = (int32 *) (AS_MSG_globals->MemBuffer + (long) imp->delta);
    i = 0;
    while (i < imp->delta_length) {
      line = GetLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  else
    imp->delta = NULL;
  GET_EOM;
  return;
}

static void Read_IMV_Mesg(FILE *fin, long indx)
{
  IntMultiVar *imv = (IntMultiVar *) (AS_MSG_globals->MemBuffer + indx);

  GET_PAIR(imv->position.bgn,imv->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(imv->num_reads,"nrd:"F_S32,"number of reads"); 
  GET_FIELD(imv->num_conf_alleles,"nca:"F_S32,"number of confirmed alleles");
  GET_FIELD(imv->anchor_size,"anc:"F_S32,"anchor size");
  GET_FIELD(imv->var_length,"len:"F_S32,"length field");
  GET_FIELD(imv->curr_var_id,"vid:"F_S32,"current VAR record id");
  GET_FIELD(imv->phased_var_id,"pid:"F_S32,"phased VAR record id");
  imv->nr_conf_alleles = (char *)GetText("nra:",fin,FALSE);
  imv->weights         = (char *)GetText("wgt:",fin,FALSE);
  imv->var_seq =         (char *)GetText("seq:",fin,FALSE);
  GET_EOM;
  return;
}

static void Read_VAR_Mesg(FILE *fin, long indx)
{
  IntMultiVar *smv = (IntMultiVar *) (AS_MSG_globals->MemBuffer + indx);

  GET_PAIR(smv->position.bgn,smv->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(smv->num_reads,"nrd:"F_S32,"number of reads");
  GET_FIELD(smv->num_conf_alleles,"nca:"F_S32,"number of confirmed alleles");
  GET_FIELD(smv->anchor_size,"anc:"F_S32,"anchor size");
  GET_FIELD(smv->var_length,"len:"F_S32,"length field");
  GET_FIELD(smv->curr_var_id,"vid:"F_S32,"current VAR record id");
  GET_FIELD(smv->phased_var_id,"pid:"F_S32,"phased VAR record id");
  smv->nr_conf_alleles = (char *)GetText("nra:",fin,FALSE);
  smv->weights         = (char *)GetText("wgt:",fin,FALSE);
  smv->var_seq         = (char *)GetText("seq:",fin,FALSE);
  GET_EOM;
  return;
}

static void Read_IUP_Mesg(FILE *fin, long indx)
{
  IntUnitigPos		*iup;
  int			i;
  int32			n, *delta;
  long			tindx;
  char			*line, *u;
  char ch;
  
  iup = (IntUnitigPos *) (AS_MSG_globals->MemBuffer + indx);
  GET_TYPE(ch,"typ:%1[URSPsX]","unitigpos type");
  iup->type = (UnitigType) ch;
  GET_FIELD(iup->ident,"lid:"F_IID,"unitigpos id");
  GET_PAIR(iup->position.bgn,iup->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(iup->delta_length,"dln:"F_S32,"delta length");
  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  if (iup->delta_length > 0) {
    tindx = MoreSpace(sizeof(int32)*iup->delta_length,8);
    iup = (IntUnitigPos *) (AS_MSG_globals->MemBuffer + indx);	// in case of realloc
    iup->delta = (int32 *) tindx;
    delta = (int32 *) (AS_MSG_globals->MemBuffer + (long) iup->delta);
    i = 0;
    while (i < iup->delta_length) {
      line = GetLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  else
    iup->delta = NULL;
  GET_EOM;
  return;
}

static void *Read_IUM_Mesg(FILE *fin)
{ static IntUnitigMesg		mesg;
  int				i;
  long				cindx, qindx, mpindx, indx;
# ifdef AS_ENABLE_SOURCE
  long				sindx;
# endif
  char ch;

  GET_FIELD(mesg.iaccession,"acc:"F_IID,"accession field");
# ifdef AS_ENABLE_SOURCE
  sindx = GetText("src:",fin,FALSE);
# endif
  GET_FIELD(mesg.coverage_stat,"cov:%f","coverage stat");
  GET_TYPE(ch,"sta:%1[UCNSX]","status");
  mesg.status = (UnitigStatus) ch;
  GET_FIELD(mesg.a_branch_point,"abp:"F_COORD,"a branch point");
  GET_FIELD(mesg.b_branch_point,"bbp:"F_COORD,"b branch point");
  GET_FIELD(mesg.length,"len:"F_COORD,"length field");
  cindx = GetText("cns:",fin,TRUE);
  qindx = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:"F_S32,"forced booleon");
  GET_FIELD(mesg.num_frags,"nfr:"F_S32,"num frags field");
  if (mesg.num_frags > 0) {
    indx = mpindx = MoreSpace(mesg.num_frags*sizeof(IntMultiPos),8);
    /* Why use 8 boundary above? Damned if I know - ela */
    for (i=0; i < mesg.num_frags; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{IMP",4) != 0)
	MgenError("Expecting IMP record");
      Read_IMP_Mesg(fin, indx);
      indx += sizeof(IntMultiPos);
    }
    mesg.f_list = (IntMultiPos *) (AS_MSG_globals->MemBuffer + mpindx);
  }
  else
    mesg.f_list = NULL;
  GET_EOM;
  mesg.consensus = AS_MSG_globals->MemBuffer + cindx;
  mesg.quality  = AS_MSG_globals->MemBuffer + qindx;

  assert(strlen(mesg.consensus) == strlen(mesg.quality));
  assert((strlen(mesg.consensus) == mesg.length) ||
	 (strlen(mesg.consensus) == 0) );

  for (i=0; i < mesg.num_frags; ++i) {
    if (mesg.f_list[i].delta_length > 0)
      mesg.f_list[i].delta = (int32 *) (AS_MSG_globals->MemBuffer+(long) mesg.f_list[i].delta);
  }
# ifdef AS_ENABLE_SOURCE
  mesg.source = AS_MSG_globals->MemBuffer + sindx;
# endif
  return ((void *) (&mesg));
}

static void *Read_IUL_Mesg(FILE *fin)
{ static IntUnitigLinkMesg	mesg;
  int				i,size;
  long				indx;
  IntMate_Pairs			*imp;
  char ch;
  
  GET_FIELD(mesg.unitig1,"ut1:"F_IID,"unitig 1 field");
  GET_FIELD(mesg.unitig2,"ut2:"F_IID,"unitig 2 field");
  GET_TYPE(ch,"ori:%1[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
  mesg.overlap_type = (UnitigOverlapType) ch;
  GET_FIELD(mesg.is_possible_chimera,"ipc:"F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  GET_TYPE(ch,"sta:%1[APBCU]","placement status");
  mesg.status = (PlacementStatusType) ch;
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    indx = MoreSpace(sizeof(IntMate_Pairs)*size,8);
    imp = mesg.jump_list = (IntMate_Pairs *) (AS_MSG_globals->MemBuffer + indx);
    for (i=0; i < size; ++i) {
      GET_TRIPLE(imp->in1,imp->in2,ch,F_IID","F_IID",%1[MSBRYT]","mate pair");
      imp->type = (LinkType)ch;
      ++imp;
    }
  }
  else
    mesg.jump_list = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_ICL_Mesg(FILE *fin)
{ static IntContigLinkMesg	mesg;
  int				i,size;
  long				indx;
  IntMate_Pairs			*imp;
  char ch;
  
  GET_FIELD(mesg.contig1,"co1:"F_IID,"contig 1 field");
  GET_FIELD(mesg.contig2,"co2:"F_IID,"contig 2 field");
  GET_TYPE(ch,"ori:%1[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
  mesg.overlap_type = (UnitigOverlapType) ch;
  GET_FIELD(mesg.is_possible_chimera,"ipc:"F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  GET_TYPE(ch,"sta:%1[APBCU]","placement status");
  mesg.status = (PlacementStatusType) ch;
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    indx = MoreSpace(sizeof(IntMate_Pairs)*size,8);
    imp = mesg.jump_list = (IntMate_Pairs *) (AS_MSG_globals->MemBuffer + indx);
    for (i=0; i < size; ++i) {
      GET_TRIPLE(imp->in1,imp->in2,ch, F_IID","F_IID",%1[MSBRYT]","mate pair");
      imp->type = (LinkType)ch;
      ++imp;
    }
  }
  else
    mesg.jump_list = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_ISL_Mesg(FILE *fin)
{ static InternalScaffoldLinkMesg	mesg;
  int				i,size;
  long				indx;
  IntMate_Pairs			*imp;
  char ch;
  
  GET_FIELD(mesg.iscaffold1,"sc1:"F_IID,"scaffold 1 field");
  GET_FIELD(mesg.iscaffold2,"sc2:"F_IID,"scaffold 2 field");
  GET_TYPE(ch,"ori:%1[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  assert(size > 0);
  indx = MoreSpace(sizeof(IntMate_Pairs)*size,8);
  imp = mesg.jump_list = (IntMate_Pairs *) (AS_MSG_globals->MemBuffer + indx);
  for (i=0; i < size; ++i) {
    GET_TRIPLE(imp->in1,imp->in2,ch, F_IID","F_IID",%1[MSBRYT]","mate pair");
    imp->type = (LinkType) ch;
    ++imp;
  }
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_AFG_Mesg(FILE *fin)
{ static AugFragMesg		mesg;
  char *line;
  int32 i,j;
  char ch;
  
  GET_PAIR(mesg.eaccession,mesg.iaccession,"acc:("F_UID","F_IID")","accession field");
  GET_TYPE(ch,"mst:%1[ZGCLSONHADEURF]","mate status");
  mesg.mate_status = (MateStatType) ch;
  GET_FIELD(mesg.chimeric,"chi:"F_S32,"chimeric flag");
#ifdef AFG_BACKWARDS_COMPATIBLE
  line = GetLine(fin, TRUE);
  if(sscanf(line,"cha:"F_S32,&i)==1){
    mesg.chaff=i;
    line = GetLine(fin, TRUE);
  }
  else
    mesg.chaff=0;

  if(sscanf(line,"clr:"F_COORD","F_COORD,&i,&j)==2){
    mesg.clear_rng.bgn=i;
    mesg.clear_rng.end=j;
  }
  else
    MfieldError("chaff flag/clear range");	
#else
  GET_FIELD(mesg.chaff,"cha:"F_S32,"chaff flag");
  GET_PAIR(mesg.clear_rng.bgn,mesg.clear_rng.end,"clr:"F_COORD","F_COORD,"clear range");
#endif
  GET_EOM;
  return ((void *) (&mesg));
}

static void Read_ICP_Mesg(FILE *fin, IntContigPairs *icp)
{
  char ch;
  GET_FIELD(icp->contig1,"ct1:"F_IID,"contig 1 id");
  GET_FIELD(icp->contig2,"ct2:"F_IID,"contig 2 id");
  GET_FIELD(icp->mean,"mea:%f","mean distance");
  GET_FIELD(icp->stddev,"std:%f","standard deviation");
  GET_TYPE(ch,"ori:%1[NAIOU]","link orientation");
  icp->orient = (ChunkOrientationType) ch;
  GET_EOM;
  return;
}

static void *Read_ISF_Mesg(FILE *fin)
{ static IntScaffoldMesg	mesg;
  int				i, num;
  long				indx;
  IntContigPairs		*icp;

  GET_FIELD(mesg.iaccession,"acc:"F_IID,"ISF accession");
  GET_FIELD(mesg.num_contig_pairs,"noc:"F_S32,"number of contigs");
  num = MAX(1,mesg.num_contig_pairs);
  if (num > 0) {
    indx = MoreSpace(num*sizeof(IntContigPairs),8);
    icp = mesg.contig_pairs = (IntContigPairs *) (AS_MSG_globals->MemBuffer + indx);
    for (i=0; i < num; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{ICP",4) != 0)
	MgenError("Expecting ICP record");
      Read_ICP_Mesg(fin,icp);
      ++icp;
    }
  }
  else
    mesg.contig_pairs = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_IMD_Mesg(FILE *fin)
{ static IntMateDistMesg	mesg;
  long				indx;
  int				i;

  GET_FIELD(mesg.refines,"ref:"F_IID,"distance id");
  GET_FIELD(mesg.mean,"mea:%f","mean distance");
  GET_FIELD(mesg.stddev,"std:%f","standard deviation");
  GET_FIELD(mesg.min,"min:"F_COORD,"min distance");
  GET_FIELD(mesg.max,"max:"F_COORD,"max distance");
  GET_FIELD(mesg.num_buckets,"buc:"F_S32,"number of buckets");
  if (strncmp(GetLine(fin,TRUE),"his:",4) != 0)
    MgenError("Expecting his field");
  if (mesg.num_buckets > 0) {
    indx = MoreSpace(mesg.num_buckets*sizeof(int32),8);
    mesg.histogram = (int32 *) (AS_MSG_globals->MemBuffer + indx);
    for (i=0; i < mesg.num_buckets; ++i)
      GET_FIELD(mesg.histogram[i],F_S32,"histogram entry");
  }
  else
    mesg.histogram = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_ICM_Mesg(FILE *fin)
{ static IntConConMesg		mesg;
  long	 cindx, qindx, mpindx, upindx, indx, uindx, vindx, vpindx;
  int				i;
  char ch;
  
  GET_FIELD(mesg.iaccession,"acc:"F_IID,"accession number");
  GET_TYPE(ch,"pla:%1[PU]","placed flag");
  mesg.placed = (ContigPlacementStatusType) ch;
  GET_FIELD(mesg.length,"len:"F_COORD,"contig length");
  cindx = GetText("cns:",fin,TRUE);
  qindx = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:"F_S32,"forced flag");
  GET_FIELD(mesg.num_pieces,"npc:"F_S32,"number of pieces");
  GET_FIELD(mesg.num_unitigs,"nou:"F_S32,"number of unitigs");
  GET_FIELD(mesg.num_vars,"nvr:"F_S32,"num vars field");
  /* Why use 8 boundary above & below? Damned if I know - ela */
  vindx = vpindx = MoreSpace(mesg.num_vars   *sizeof(IntMultiVar),8); 
  indx  = mpindx = MoreSpace(mesg.num_pieces *sizeof(IntMultiPos),8);
  uindx = upindx = MoreSpace(mesg.num_unitigs*sizeof(IntUnitigPos),8);

  mesg.v_list = NULL;
  mesg.pieces = NULL;
  mesg.unitigs = NULL;

  if (mesg.num_vars > 0) {
    mesg.v_list = (IntMultiVar *) (AS_MSG_globals->MemBuffer + vpindx);
    for (i=0; i < mesg.num_vars; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{IMV",4) != 0)
        MgenError("Expecting IMV record");
      Read_IMV_Mesg(fin, vindx);
      vindx += sizeof(IntMultiVar);
    }
  }

  if (mesg.num_pieces > 0) {
    mesg.pieces = (IntMultiPos *) (AS_MSG_globals->MemBuffer + mpindx);
    for (i=0; i < mesg.num_pieces; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{IMP",4) != 0)
        MgenError("Expecting IMP record");
      Read_IMP_Mesg(fin, indx);
      indx += sizeof(IntMultiPos);
    }
  }  

  if (mesg.num_unitigs > 0) {
    mesg.unitigs = (IntUnitigPos *) (AS_MSG_globals->MemBuffer + upindx);
    for (i=0; i < mesg.num_unitigs; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{IUP",4) != 0)
	MgenError("Expecting IUP record");
      Read_IUP_Mesg(fin, uindx);
      uindx += sizeof(IntUnitigPos);
    }
  }

  GET_EOM;

  mesg.consensus = AS_MSG_globals->MemBuffer + cindx;
  mesg.quality = AS_MSG_globals->MemBuffer + qindx;

  if (mesg.num_vars > 0)
    mesg.v_list = (IntMultiVar *) (AS_MSG_globals->MemBuffer + vpindx);
  for (i=0; i < mesg.num_vars; ++i) {
    mesg.v_list[i].nr_conf_alleles = AS_MSG_globals->MemBuffer + (long) mesg.v_list[i].nr_conf_alleles;
    mesg.v_list[i].weights = AS_MSG_globals->MemBuffer + (long) mesg.v_list[i].weights;
    mesg.v_list[i].var_seq = AS_MSG_globals->MemBuffer + (long) mesg.v_list[i].var_seq;
  }

  if (mesg.num_pieces > 0)
    mesg.pieces = (IntMultiPos *) (AS_MSG_globals->MemBuffer + mpindx);
  for (i=0; i < mesg.num_pieces; ++i) {
    if (mesg.pieces[i].delta_length > 0)
      mesg.pieces[i].delta = (int32 *) (AS_MSG_globals->MemBuffer + (long) mesg.pieces[i].delta);
  }
  if (mesg.num_unitigs > 0)
    mesg.unitigs = (IntUnitigPos *) (AS_MSG_globals->MemBuffer + upindx);
  for (i=0; i < mesg.num_unitigs; ++i) {
    if (mesg.unitigs[i].delta_length > 0)
      mesg.unitigs[i].delta = (int32 *) (AS_MSG_globals->MemBuffer + (long) mesg.unitigs[i].delta);
  }
  return ((void *) (&mesg));
}


static void *Read_IAF_Mesg(FILE *fin)
{ static IntAugFragMesg		mesg;
  char ch;
  
  GET_FIELD(mesg.iaccession,"acc:"F_IID,"accession field");
  GET_TYPE(ch,"typ:%1[RXELTFSUCBWG]","type");
  mesg.type = (FragType) ch;
  GET_FIELD(mesg.chimeric,"chi:"F_S32,"chimeric flag");

#ifdef IAF_BACKWARDS_COMPATIBLE
  {
    int i, j;
    char *line;

    line = GetLine(fin, TRUE);
    if(sscanf(line,"cha:"F_S32,&i)==1)
      {
        mesg.chaff=i;
        line = GetLine(fin, TRUE);
      }
    else
      mesg.chaff=0;
    if(sscanf(line,"clr:"F_COORD","F_COORD,&i,&j)==2)
      {
        mesg.clear_rng.bgn=i;
        mesg.clear_rng.end=j;
      }
    else
      MgenError("IAF: choked on cha: or chr:");
  }
#else
  GET_FIELD(mesg.chaff,"cha:"F_S32,"chaff flag");
  GET_PAIR(mesg.clear_rng.bgn,mesg.clear_rng.end,"clr:"F_COORD","F_COORD,"clear range");
#endif
  GET_TYPE(ch,"mst:%1[ZGCLSONHADEURF]","mate status");
  mesg.mate_status = (MateStatType) ch;
  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_EOF_Mesg(FILE *fin)
{
  static EndOfFileMesg mesg;
  int commentidx;
  time_t entry_time;

  GET_FIELD(mesg.status,"sta:"F_S32,"status field");
  GET_FIELD(entry_time,"crt:"F_TIME_T,"entry time field");  //  Unused
  commentidx = GetText("com:",fin,FALSE);
  mesg.comment = AS_MSG_globals->MemBuffer + commentidx;
  GET_EOM;

  return ((void *) &mesg);
}



/* Genome snapshot input routines */
/**********************************/

static void Read_MPS_Mesg(FILE *fin, long indx)
{
  SnapMultiPos		*imp;
  int			i;
  int32			n, *delta;
  long			tindx;
  char			*line, *u;
  char ch;
  
  imp = (SnapMultiPos *) (AS_MSG_globals->MemBuffer + indx);
  GET_TYPE(ch,"typ:%1[RXTEFUSLuBG]","multipos type");
  imp->type = (FragType) ch;
  GET_FIELD(imp->eident,"mid:"F_UID,"multipos id");
#ifdef AS_ENABLE_SOURCE
  tindx = GetText("src:",fin,FALSE);
  imp = (SnapMultiPos *) (AS_MSG_globals->MemBuffer + indx);	// in case of realloc
  imp->source = (char *) tindx;
#endif
  GET_PAIR(imp->position.bgn,imp->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(imp->delta_length,"dln:"F_S32,"delta length");
  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  if (imp->delta_length > 0) {
    tindx = MoreSpace(sizeof(int32)*imp->delta_length,8);
    imp = (SnapMultiPos *) (AS_MSG_globals->MemBuffer + indx);	// in case of realloc
    imp->delta = (int32 *) tindx;
    delta = (int32 *) (AS_MSG_globals->MemBuffer + (long) imp->delta);
    i = 0;
    while (i < imp->delta_length) {
      line = GetLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  else
    imp->delta = NULL;
  GET_EOM;
  return;
}

static void Read_UPS_Mesg(FILE *fin, long indx)
{
  UnitigPos		*iup;
  int			i;
  int32			n, *delta;
  long			tindx;
  char			*line, *u;
  char ch;
  
  iup = (UnitigPos *) (AS_MSG_globals->MemBuffer + indx);
  GET_TYPE(ch,"typ:%1[URSPs]","unitigpos type");
  iup->type = (UnitigType) ch;
  GET_FIELD(iup->eident,"lid:"F_UID,"unitig id");
  GET_PAIR(iup->position.bgn,iup->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(iup->delta_length,"dln:"F_S32,"delta length");
  if (strncmp(GetLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  if (iup->delta_length > 0) {
    tindx = MoreSpace(sizeof(int32)*iup->delta_length,8);
    iup = (UnitigPos *) (AS_MSG_globals->MemBuffer + indx);	// in case of realloc
    iup->delta = (int32 *) tindx;
    delta = (int32 *) (AS_MSG_globals->MemBuffer + (long) iup->delta);
    i = 0;
    while (i < iup->delta_length) {
      line = GetLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  else
    iup->delta = NULL;
  GET_EOM;
  return;
}


static void Read_CTP_Mesg(FILE *fin, SnapContigPairs *icp)
{
  char ch;
  GET_FIELD(icp->econtig1,"ct1:"F_UID ,"contig 1 id");
  GET_FIELD(icp->econtig2,"ct2:"F_UID ,"contig 2 id");
  GET_FIELD(icp->mean,"mea:%f","mean distance");
  GET_FIELD(icp->stddev,"std:%f","standard deviation");
  GET_TYPE(ch,"ori:%1[NAIOU]","link orientation");
  icp->orient = (ChunkOrientationType) ch;
  GET_EOM;
  return;
}


static void *Read_UTG_Mesg(FILE *fin)
{ static SnapUnitigMesg		mesg;
  int				i;
  long				cindx, qindx, mpindx, indx;
#ifdef AS_ENABLE_SOURCE
  long				sindx;
#endif
  char ch;

  GET_PAIR(mesg.eaccession,mesg.iaccession,"acc:("F_UID","F_IID")","accession field");
#ifdef AS_ENABLE_SOURCE
  sindx = GetText("src:",fin,FALSE);
#endif
  GET_FIELD(mesg.coverage_stat,"cov:%f","coverage stat");
  GET_TYPE(ch,"sta:%1[UCNSX]","status");
  mesg.status = (UnitigStatus) ch;
  GET_FIELD(mesg.a_branch_point,"abp:"F_COORD,"a branch point");
  GET_FIELD(mesg.b_branch_point,"bbp:"F_COORD,"b branch point");
  GET_FIELD(mesg.length,"len:"F_COORD,"length field");
  cindx = GetText("cns:",fin,TRUE);
  qindx = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:"F_S32,"forced booleon");
  GET_FIELD(mesg.num_frags,"nfr:"F_S32,"num frags field");
  if (mesg.num_frags > 0) {
    indx = mpindx = MoreSpace(mesg.num_frags*sizeof(SnapMultiPos),8);
    /* Why use 8 boundary above? Damned if I know - ela */
    for (i=0; i < mesg.num_frags; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{MPS",4) != 0)
	MgenError("Expecting MPS record");
      Read_MPS_Mesg(fin, indx);
      indx += sizeof(SnapMultiPos);
    }
    mesg.f_list = (SnapMultiPos *) (AS_MSG_globals->MemBuffer + mpindx);
  }
  else
    mesg.f_list = NULL;
  GET_EOM;
  mesg.consensus = AS_MSG_globals->MemBuffer + cindx;
  mesg.quality  = AS_MSG_globals->MemBuffer + qindx;
  for (i=0; i < mesg.num_frags; ++i) {
#ifdef AS_ENABLE_SOURCE
    mesg.f_list[i].source = AS_MSG_globals->MemBuffer + (long) mesg.f_list[i].source;
#endif
    if (mesg.f_list[i].delta_length > 0)
      mesg.f_list[i].delta = (int32 *) (AS_MSG_globals->MemBuffer+(long) mesg.f_list[i].delta);
  }
#ifdef AS_ENABLE_SOURCE
  mesg.source = AS_MSG_globals->MemBuffer + sindx;
#endif
  return ((void *) (&mesg));
}


static void *Read_ULK_Mesg(FILE *fin)
{ static SnapUnitigLinkMesg	mesg;
  int				i,size;
  long				indx;
  SnapMate_Pairs	       	*imp;
  char ch;

  GET_FIELD(mesg.eunitig1,"ut1:"F_UID,"unitig 1 field");
  GET_FIELD(mesg.eunitig2,"ut2:"F_UID,"unitig 2 field");
  GET_TYPE(ch,"ori:%1[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
  mesg.overlap_type = (UnitigOverlapType) ch;
  GET_FIELD(mesg.is_possible_chimera,"ipc:"F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  GET_TYPE(ch,"sta:%1[APBCU]","placement status");
  mesg.status = (PlacementStatusType) ch;
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    indx = MoreSpace(sizeof(SnapMate_Pairs)*size,8);
    imp = mesg.jump_list = (SnapMate_Pairs *) (AS_MSG_globals->MemBuffer + indx);
    for (i=0; i < size; ++i) {
      GET_TRIPLE(imp->in1,imp->in2,ch,F_UID","F_UID",%1[MSBRYT]","mate pair");
      imp->type = (LinkType)ch;
      ++imp;
    }
  }
  else
    mesg.jump_list = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_CCO_Mesg(FILE *fin)
{ static SnapConConMesg		mesg;
  long	 cindx, qindx, mpindx, upindx, indx, uindx, vindx, vpindx;
  int  	 i;
  char   ch;
  
  GET_PAIR(mesg.eaccession,mesg.iaccession,"acc:("F_UID","F_IID")","accession number");
  GET_TYPE(ch,"pla:%1[PU]","placed flag");
  mesg.placed = (ContigPlacementStatusType) ch;
  GET_FIELD(mesg.length,"len:"F_COORD,"contig length");
  cindx = GetText("cns:",fin,TRUE);
  qindx = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:"F_S32,"forced flag");
  GET_FIELD(mesg.num_pieces,"npc:"F_S32,"number of pieces");
  GET_FIELD(mesg.num_unitigs,"nou:"F_S32,"number of unitigs");
  GET_FIELD(mesg.num_vars,"nvr:"F_S32,"number of vars");
  /* Why use 8 boundary above & below? Damned if I know - ela */
  vindx = vpindx = MoreSpace(mesg.num_vars  *sizeof(IntMultiVar),8);
  indx  = mpindx = MoreSpace(mesg.num_pieces*sizeof(SnapMultiPos),8);
  uindx =upindx = MoreSpace(mesg.num_unitigs*sizeof(UnitigPos),8);

  mesg.vars = NULL;
  mesg.unitigs = NULL;

  if (mesg.num_vars > 0) {
    mesg.vars = (IntMultiVar *) (AS_MSG_globals->MemBuffer + vpindx);
    for (i=0; i < mesg.num_vars; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{VAR",4) != 0)
        MgenError("Expecting VAR record");
      Read_VAR_Mesg(fin, vindx);
      vindx += sizeof(IntMultiVar);
    }
  }

  for (i=0; i < mesg.num_pieces; ++i) {
    if (strncmp(GetLine(fin,TRUE),"{MPS",4) != 0)
      MgenError("Expecting MPS record");
    Read_MPS_Mesg(fin, indx);
    indx += sizeof(SnapMultiPos);
  }

  if (mesg.num_unitigs > 0) {
    mesg.unitigs  = (UnitigPos *) (AS_MSG_globals->MemBuffer + upindx);
    for (i=0; i < mesg.num_unitigs; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{UPS",4) != 0)
	MgenError("Expecting UPS record");
      Read_UPS_Mesg(fin, uindx);
      uindx += sizeof(UnitigPos);
    }
  }

  GET_EOM;

  mesg.consensus = AS_MSG_globals->MemBuffer + cindx;
  mesg.quality = AS_MSG_globals->MemBuffer + qindx;

  if (mesg.num_vars > 0)
    mesg.vars = (IntMultiVar *) (AS_MSG_globals->MemBuffer + vpindx);
  for (i=0; i < mesg.num_vars; ++i) {
    mesg.vars[i].nr_conf_alleles = AS_MSG_globals->MemBuffer + (long) mesg.vars[i].nr_conf_alleles;
    mesg.vars[i].weights = AS_MSG_globals->MemBuffer + (long) mesg.vars[i].weights;
    mesg.vars[i].var_seq = AS_MSG_globals->MemBuffer + (long) mesg.vars[i].var_seq;
  } 

  if (mesg.num_pieces > 0)
    mesg.pieces = (SnapMultiPos *) (AS_MSG_globals->MemBuffer + mpindx);
  for (i=0; i < mesg.num_pieces; ++i) {
#ifdef AS_ENABLE_SOURCE
    mesg.pieces[i].source = AS_MSG_globals->MemBuffer + (long) mesg.pieces[i].source;
#endif
    if (mesg.pieces[i].delta_length > 0)
      mesg.pieces[i].delta = (int32 *) (AS_MSG_globals->MemBuffer + (long) mesg.pieces[i].delta);
  }

  if (mesg.num_unitigs > 0)
    mesg.unitigs = (UnitigPos *) (AS_MSG_globals->MemBuffer + upindx);
  for (i=0; i < mesg.num_unitigs; ++i) {
    if (mesg.unitigs[i].delta_length > 0)
      mesg.unitigs[i].delta = (int32 *) (AS_MSG_globals->MemBuffer + (long) mesg.unitigs[i].delta);
  }
  
  return ((void *) (&mesg));
}



static void *Read_CLK_Mesg(FILE *fin)
{ static SnapContigLinkMesg	mesg;
  int				i,size;
  long				indx;
  SnapMate_Pairs		*imp;
  char ch;
  // char*                         line;
  
  GET_FIELD(mesg.econtig1,"c%*1[ot]1:"F_UID,"unitig 1 field");
  GET_FIELD(mesg.econtig2,"c%*1[ot]2:"F_UID,"unitig 2 field");

  GET_TYPE(ch,"ori:%1[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
  GET_TYPE(ch,"ovt:%1[NOTCIMXYZ]","overlap type");
  mesg.overlap_type = (UnitigOverlapType) ch;
  GET_FIELD(mesg.is_possible_chimera,"ipc:"F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  GET_TYPE(ch,"sta:%1[APBCU]","placement status");
  mesg.status = (PlacementStatusType) ch;
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    indx = MoreSpace(sizeof(SnapMate_Pairs)*size,8);
    imp = mesg.jump_list = (SnapMate_Pairs *) (AS_MSG_globals->MemBuffer + indx);
    for (i=0; i < size; ++i) {
      GET_TRIPLE(imp->in1,imp->in2,ch,F_UID","F_UID",%1[MSBRYT]","mate pair");
      imp->type = (LinkType)ch;
      ++imp;
    }
  }
  else
    mesg.jump_list = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_SLK_Mesg(FILE *fin)
{ static SnapScaffoldLinkMesg	mesg;
  int				i,size;
  long				indx;
  SnapMate_Pairs		*imp;
  char ch;
  
  GET_FIELD(mesg.escaffold1,"sc1:"F_UID,"scaffold 1 field");
  GET_FIELD(mesg.escaffold2,"sc2:"F_UID,"scaffold 2 field");

  GET_TYPE(ch,"ori:%1[NAOI]","orientation");
  mesg.orientation = (ChunkOrientationType) ch;
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  if (strncmp(GetLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  assert(size > 0) ;
  indx = MoreSpace(sizeof(SnapMate_Pairs)*size,8);
  imp = mesg.jump_list = (SnapMate_Pairs *) (AS_MSG_globals->MemBuffer + indx);
  for (i=0; i < size; ++i) {
    GET_TRIPLE(imp->in1,imp->in2,ch,F_UID","F_UID",%1[MSBRYT]","mate pair");
    imp->type = (LinkType)ch;
    ++imp;
  }
  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_SCF_Mesg(FILE *fin)
{ static SnapScaffoldMesg	mesg;
  int				i, num;
  long				indx;
  SnapContigPairs		*icp;

  GET_PAIR(mesg.eaccession,mesg.iaccession,"acc:("F_UID","F_IID")","accession number");
  GET_FIELD(mesg.num_contig_pairs,"noc:"F_S32,"number of contigs");
  num = MAX(mesg.num_contig_pairs, 1);
  if (num > 0) {
    indx = MoreSpace(num*sizeof(SnapContigPairs),8);
    icp = mesg.contig_pairs = (SnapContigPairs *) (AS_MSG_globals->MemBuffer + indx);
    for (i=0; i < num; ++i) {
      if (strncmp(GetLine(fin,TRUE),"{CTP",4) != 0)
	MgenError("Expecting CTP record");
      Read_CTP_Mesg(fin,icp);
      ++icp;
    }
  }
  else
    mesg.contig_pairs = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_DSC_Mesg(FILE *fin)
{ static SnapDegenerateScaffoldMesg	mesg;

  GET_FIELD(mesg.eaccession,"acc:"F_UID,"scaffold ID");
  GET_FIELD(mesg.econtig,"ctg:"F_UID,"contig ID");
  GET_EOM;
  return ((void *) (&mesg));
}

static void *Read_IDS_Mesg(FILE *fin)
{ static IntDegenerateScaffoldMesg	mesg;

  GET_FIELD(mesg.icontig,"ctg:"F_IID,"contig ID");
  GET_EOM;
  return ((void *) (&mesg));
}



static void *Read_MDI_Mesg(FILE *fin)
{ static SnapMateDistMesg	mesg;
  long				indx;
  int				i;

  GET_PAIR(mesg.erefines,mesg.irefines,"ref:("F_UID","F_IID")","distance id");
  GET_FIELD(mesg.mean,"mea:%f","mean distance");
  GET_FIELD(mesg.stddev,"std:%f","standard deviation");
  GET_FIELD(mesg.min,"min:"F_COORD,"min distance");
  GET_FIELD(mesg.max,"max:"F_COORD,"max distance");
  GET_FIELD(mesg.num_buckets,"buc:"F_S32,"number of buckets");
  if (strncmp(GetLine(fin,TRUE),"his:",4) != 0)
    MgenError("Expecting his field");
  if (mesg.num_buckets > 0) {
    indx = MoreSpace(mesg.num_buckets*sizeof(int32),8);
    mesg.histogram = (int32 *) (AS_MSG_globals->MemBuffer + indx);
    for (i=0; i < mesg.num_buckets; ++i)
      GET_FIELD(mesg.histogram[i],F_S32,"histogram entry");
  }
  else
    mesg.histogram = NULL;
  GET_EOM;
  return ((void *) (&mesg));
}


static void *Read_BAT_Mesg(FILE *fin){
  static BatchMesg mesg;
  int nameidx, commentidx;
  time_t entry_time;
  nameidx = GetString("bna:",fin);
  GET_FIELD(entry_time,"crt:"F_TIME_T,"entry time field");  //  Unused
  GET_FIELD(mesg.eaccession,"acc:"F_UID,"accession number");
  commentidx = GetText("com:",fin, FALSE);

  mesg.comment = AS_MSG_globals->MemBuffer + commentidx;
  mesg.name = AS_MSG_globals->MemBuffer + nameidx;

  GET_EOM;
  return ((void *) (&mesg));
}

/******************** OUTPUT ROUTINES ***************************/

/*  Routine to output each type of proto-IO message. */


static void Write_Dist_Mesg(FILE *fout, void *vmesg, int external)
{ InternalDistMesg *mesg = (InternalDistMesg *) vmesg;

  fprintf(fout,"{%s\n",(external?"DST":"IDT"));
  fprintf(fout,"act:%c\n",mesg->action);
  if (external)
    fprintf(fout,"acc:"F_UID"\n",mesg->eaccession);
  else
    fprintf(fout,"acc:("F_UID","F_IID")\n",mesg->eaccession,mesg->iaccession);
  if (mesg->action != AS_DELETE)
    { fprintf(fout,"mea:%.3f\n",mesg->mean);
      fprintf(fout,"std:%.3f\n",mesg->stddev);
    }
  fprintf(fout,"}\n");
}

static void Write_DST_Mesg(FILE *fout, void *vmesg)
{ Write_Dist_Mesg(fout,vmesg,1); }

static void Write_IDT_Mesg(FILE *fout, void *vmesg)
{ Write_Dist_Mesg(fout,vmesg,0); }

static void Write_ADL_Struct(FILE *fout, AuditLine *mesg)
{ fprintf(fout,"{ADL\n");
  fprintf(fout,"who:%s\n",mesg->name);
  fprintf(fout,"ctm:"F_TIME_T"\n",mesg->complete);
  fprintf(fout,"vsn:%s\n",mesg->version);
  PutText(fout,"com:",mesg->comment,FALSE);
  fprintf(fout,"}\n");
}

static void Write_ADT_Mesg(FILE *fout, void *vmesg)
{ AuditMesg *mesg = (AuditMesg *) vmesg;
  AuditLine *a;

  fprintf(fout,"{ADT\n");
  for (a = mesg->list; a != NULL; a = a->next)
    Write_ADL_Struct(fout,a);
  fprintf(fout,".\n");
  fprintf(fout,"}\n");
}

static void Write_VER_Mesg(FILE *fout, void *vmesg) {
  VersionMesg *mesg = (VersionMesg *) vmesg;

  fprintf(fout,"{VER\n");
  fprintf(fout,"ver:"F_U32"\n", mesg->version);
  fprintf(fout,"}\n");
}
         
static void Write_LKG_Mesg(FILE *fout, void *vmesg)
{ LinkMesg *mesg = (LinkMesg *) vmesg;

  fprintf(fout,"{LKG\n");
  fprintf(fout,"act:%c\n",mesg->action);
  fprintf(fout,"typ:%c\n",(char) mesg->type);
  fprintf(fout,"fg1:"F_UID"\n",mesg->frag1);
  fprintf(fout,"fg2:"F_UID"\n",mesg->frag2);
  if((mesg->action == AS_ADD) || (mesg->action == AS_IGNORE))
    { fprintf(fout,"etm:0\n");
      fprintf(fout,"dst:"F_UID"\n",mesg->distance);
      fprintf(fout,"ori:%c\n",mesg->link_orient);
    }
  fprintf(fout,"}\n");
}

static void Write_Frag_Mesg(FILE *fout, void *vmesg, int frag_class) {
  FragMesg *mesg = (FragMesg *) vmesg;

  assert((frag_class == MESG_FRG) || (frag_class == MESG_IFG) || (frag_class == MESG_OFG));

  fprintf(fout,"{%s\n",MessageTypeName[frag_class]);
  fprintf(fout,"act:%c\n",mesg->action);
  if (frag_class == MESG_FRG)
    fprintf(fout,"acc:"F_UID"\n",mesg->eaccession);
  else
    fprintf(fout,"acc:("F_UID","F_IID")\n",mesg->eaccession,mesg->iaccession);

  if ((mesg->action == AS_ADD) || (mesg->action == AS_IGNORE)) {
    fprintf(fout,"typ:%c\n",(char) mesg->type);
    PutText(fout,"src:",mesg->source,FALSE);
    fprintf(fout,"etm:0\n");
    if( frag_class != MESG_OFG ) {
      PutText(fout,"seq:",mesg->sequence,TRUE);
      PutText(fout,"qlt:",mesg->quality,TRUE);
    }
    fprintf(fout,"clr:"F_COORD","F_COORD"\n", mesg->clear_rng.bgn,mesg->clear_rng.end);
  }

  fprintf(fout,"}\n");
}

static void Write_FRG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,MESG_FRG); }

static void Write_IFG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,MESG_IFG); }

static void Write_OFG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,MESG_OFG); }


static void Write_OVL_Mesg(FILE *fout, void *vmesg)
{ OverlapMesg *omesg = (OverlapMesg *) vmesg;
  int i;

  fprintf(fout,"{OVL\n");
  fprintf(fout,"afr:"F_IID"\n",omesg->aifrag);
  fprintf(fout,"bfr:"F_IID"\n",omesg->bifrag);
  fprintf(fout,"ori:%c\n",omesg->orientation);
  fprintf(fout,"olt:%c\n",omesg->overlap_type);
  fprintf(fout,"ahg:"F_COORD"\n",omesg->ahg);
  fprintf(fout,"bhg:"F_COORD"\n",omesg->bhg);
  fprintf(fout,"qua:%.6f\n",omesg->quality);
  fprintf(fout,"mno:"F_COORD"\n",omesg->min_offset);
  fprintf(fout,"mxo:"F_COORD"\n",omesg->max_offset);
  fprintf(fout,"pct:"F_S32"\n",omesg->polymorph_ct);
  fprintf(fout,"del:\n");
  for (i = 0; omesg->delta[i] != AS_ENDOF_DELTA_CODE; i++)
    fprintf(fout,"%4d%c",omesg->delta[i], (i%15 == 14) ? '\n' : ' ');
  fprintf(fout,"\n");
  fprintf(fout,".\n");
  fprintf(fout,"}\n");
}

static void Write_UOM_Mesg(FILE *fout, void *vmesg)
{ UnitigOverlapMesg *mesg = (UnitigOverlapMesg *) vmesg;

  fprintf(fout,"{UOM\n");
  fprintf(fout,"ck1:"F_IID"\n",mesg->chunk1);
  fprintf(fout,"ck2:"F_IID"\n",mesg->chunk2);
  fprintf(fout,"ori:%c\n",mesg->orient);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
#ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mesg->source,FALSE);
#endif
  fprintf(fout,"len:"F_COORD"\n",mesg->best_overlap_length);
  fprintf(fout,"min:"F_COORD"\n",mesg->min_overlap_length);
  fprintf(fout,"max:"F_COORD"\n",mesg->max_overlap_length);
  fprintf(fout,"qua:%.6f\n",mesg->quality);
  fprintf(fout,"}\n");
  return;
}

static void Write_IMP_Mesg(FILE *fout, IntMultiPos *mlp)
{ int i;

  fprintf(fout,"{IMP\n");
  fprintf(fout,"typ:%c\n",(char) mlp->type);
  fprintf(fout,"mid:"F_IID"\n",mlp->ident);
  fprintf(fout,"con:"F_IID"\n",mlp->contained);
#ifdef NEW_UNITIGGER_INTERFACE
  fprintf(fout,"bid:"F_IID"\n",mlp->ident2);
#endif
  fprintf(fout,"pos:"F_COORD","F_COORD"\n",
          mlp->position.bgn,mlp->position.end);
#ifdef NEW_UNITIGGER_INTERFACE
  fprintf(fout,"ahg:"F_S32"\n",mlp->ahang);
  fprintf(fout,"bhg:"F_S32"\n",mlp->bhang); 
#endif
  fprintf(fout,"dln:"F_S32"\n",mlp->delta_length);
  fprintf(fout,"del:\n");
  if (mlp->delta_length > 0 ) {
    for(i=0; i < mlp->delta_length; i++) {
      fprintf(fout,F_S32"%c", mlp->delta[i], (i%20 == 19) ? '\n' : ' ');
    }
    if (mlp->delta_length%20 != 0) fprintf(fout,"\n");
  }
  fprintf(fout,"}\n");
  return;
}

static void Write_IMV_Mesg(FILE *fout, IntMultiVar *imv)
{
  fprintf(fout,"{IMV\n");
  fprintf(fout,"pos:"F_COORD","F_COORD"\n",imv->position.bgn,imv->position.end);
  fprintf(fout,"nrd:"F_S32"\n",imv->num_reads);
  fprintf(fout,"nca:"F_S32"\n",imv->num_conf_alleles);
  fprintf(fout,"anc:"F_S32"\n",imv->anchor_size);
  fprintf(fout,"len:"F_S32"\n",imv->var_length);
  fprintf(fout,"vid:"F_S32"\n",imv->curr_var_id);
  fprintf(fout,"pid:"F_S32"\n",imv->phased_var_id);
  PutText(fout,"nra:",imv->nr_conf_alleles,FALSE);
  PutText(fout,"wgt:",imv->weights,FALSE);
  PutText(fout,"seq:",imv->var_seq,FALSE);
  fprintf(fout,"}\n");
  return;
}

static void Write_VAR_Mesg(FILE *fout, IntMultiVar *smv)
{
  fprintf(fout,"{VAR\n");
  fprintf(fout,"pos:"F_COORD","F_COORD"\n",smv->position.bgn,smv->position.end);
  fprintf(fout,"nrd:"F_S32"\n",smv->num_reads);
  fprintf(fout,"nca:"F_S32"\n",smv->num_conf_alleles);
  fprintf(fout,"anc:"F_S32"\n",smv->anchor_size);
  fprintf(fout,"len:"F_S32"\n",smv->var_length);
  fprintf(fout,"vid:"F_S32"\n",smv->curr_var_id);
  fprintf(fout,"pid:"F_S32"\n",smv->phased_var_id);
  PutText(fout,"nra:",smv->nr_conf_alleles,FALSE);
  PutText(fout,"wgt:",smv->weights,FALSE);
  PutText(fout,"seq:",smv->var_seq,FALSE);
  fprintf(fout,"}\n");
  return;
}

static void Write_IUP_Mesg(FILE *fout, IntUnitigPos *up)
{ int i;

  fprintf(fout,"{IUP\n");
  fprintf(fout,"typ:%c\n",(char) up->type);
  fprintf(fout,"lid:"F_IID"\n",up->ident);
  fprintf(fout,"pos:"F_COORD","F_COORD"\n",up->position.bgn,up->position.end);
  fprintf(fout,"dln:"F_S32"\n",up->delta_length);
  fprintf(fout,"del:\n");
  if (up->delta_length > 0 ) {
    for(i=0; i < up->delta_length; i++)
      fprintf(fout,F_S32"%c",up->delta[i], (i%20 == 19) ? '\n' : ' ');
    if (up->delta_length%20 != 0) fprintf(fout,"\n");
  }
  fprintf(fout,"}\n");
  return;
}

static void Write_IUM_Mesg(FILE *fout, void *vmesg)
{ IntUnitigMesg *mesg = (IntUnitigMesg *) vmesg;
  int			i;

  assert(strlen(mesg->consensus) == strlen(mesg->quality));
  assert((strlen(mesg->consensus) == mesg->length) || 
	 (strlen(mesg->consensus) == 0) );
  fprintf(fout,"{IUM\n");
  fprintf(fout,"acc:"F_IID"\n",mesg->iaccession);
# ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mesg->source,FALSE);
# endif
  fprintf(fout,"cov:%.3f\n",mesg->coverage_stat);
  fprintf(fout,"sta:%c\n",mesg->status);
  fprintf(fout,"abp:"F_COORD"\n",mesg->a_branch_point);
  fprintf(fout,"bbp:"F_COORD"\n",mesg->b_branch_point);
  fprintf(fout,"len:"F_COORD"\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:"F_S32"\n",mesg->forced);
  fprintf(fout,"nfr:"F_S32"\n",mesg->num_frags);
  for (i=0; i < mesg->num_frags; ++i)
    Write_IMP_Mesg(fout,&(mesg->f_list[i]));
  fprintf(fout,"}\n");
  return;
}

static void Write_IUL_Mesg(FILE *fout, void *vmesg)
{ IntUnitigLinkMesg *mesg = (IntUnitigLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{IUL\n");
  fprintf(fout,"ut1:"F_IID"\n",mesg->unitig1);
  fprintf(fout,"ut2:"F_IID"\n",mesg->unitig2);
  fprintf(fout,"ori:%c\n",mesg->orientation);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  fprintf(fout,"ipc:"F_S32"\n",mesg->is_possible_chimera);
  fprintf(fout,"gui:"F_S32"\n",mesg->includes_guide);
  fprintf(fout,"mea:%.3f\n",mesg->mean_distance);
  fprintf(fout,"std:%.3f\n",mesg->std_deviation);
  fprintf(fout,"num:"F_S32"\n",mesg->num_contributing);
  fprintf(fout,"sta:%c\n",mesg->status);
  npairs = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --npairs;
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout,F_IID","F_IID",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_ICL_Mesg(FILE *fout, void *vmesg)
{ IntContigLinkMesg *mesg = (IntContigLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{ICL\n");
  fprintf(fout,"co1:"F_IID"\n",mesg->contig1);
  fprintf(fout,"co2:"F_IID"\n",mesg->contig2);
  fprintf(fout,"ori:%c\n",mesg->orientation);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  fprintf(fout,"ipc:"F_S32"\n",mesg->is_possible_chimera);
  fprintf(fout,"gui:"F_S32"\n",mesg->includes_guide);
  fprintf(fout,"mea:%.3f\n",mesg->mean_distance);
  fprintf(fout,"std:%.3f\n",mesg->std_deviation);
  fprintf(fout,"num:"F_S32"\n",mesg->num_contributing);
  fprintf(fout,"sta:%c\n",mesg->status);
  npairs = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --npairs;
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout,F_IID","F_IID",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_ISL_Mesg(FILE *fout, void *vmesg)
{ InternalScaffoldLinkMesg *mesg = (InternalScaffoldLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{ISL\n");
  fprintf(fout,"sc1:"F_IID"\n",mesg->iscaffold1);
  fprintf(fout,"sc2:"F_IID"\n",mesg->iscaffold2);
  fprintf(fout,"ori:%c\n",mesg->orientation);
  fprintf(fout,"gui:"F_S32"\n",mesg->includes_guide);
  fprintf(fout,"mea:%.3f\n",mesg->mean_distance);
  fprintf(fout,"std:%.3f\n",mesg->std_deviation);
  fprintf(fout,"num:"F_S32"\n",mesg->num_contributing);
  npairs = mesg->num_contributing;
  assert(npairs > 0);
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout,F_IID","F_IID",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_AFG_Mesg(FILE *fout, void *vmesg)
{ AugFragMesg *mesg = (AugFragMesg *) vmesg;
  
  fprintf(fout,"{AFG\n");
  fprintf(fout,"acc:("F_UID","F_IID")\n",mesg->eaccession,mesg->iaccession);
  fprintf(fout,"mst:%c\n",mesg->mate_status);
  fprintf(fout,"chi:"F_S32"\n",mesg->chimeric);
  fprintf(fout,"cha:"F_S32"\n",mesg->chaff);
  fprintf(fout,"clr:"F_COORD","F_COORD"\n", mesg->clear_rng.bgn,mesg->clear_rng.end);
  fprintf(fout,"}\n");
  return;
}

static void Write_ICP_Mesg(FILE *fout, IntContigPairs *mesg)
{
  fprintf(fout,"{ICP\n");
  fprintf(fout,"ct1:"F_IID"\n",mesg->contig1);
  fprintf(fout,"ct2:"F_IID"\n",mesg->contig2);
  fprintf(fout,"mea:%.3f\n",mesg->mean);
  fprintf(fout,"std:%.3f\n",mesg->stddev);
  fprintf(fout,"ori:%c\n",mesg->orient);
  fprintf(fout,"}\n");
  return;
}

static void Write_ISF_Mesg(FILE *fout, void *vmesg)
{ IntScaffoldMesg *mesg = (IntScaffoldMesg *) vmesg;
  int		i;
  int num = MAX(1, mesg->num_contig_pairs);
  fprintf(fout,"{ISF\n");
  fprintf(fout,"acc:"F_IID"\n", mesg->iaccession);
  fprintf(fout,"noc:"F_S32"\n",mesg->num_contig_pairs);
  for (i=0; i < num; ++i)
    Write_ICP_Mesg(fout,&mesg->contig_pairs[i]);
  fprintf(fout,"}\n");
  return;
}

static void Write_IMD_Mesg(FILE *fout, void *vmesg)
{ IntMateDistMesg *mesg = (IntMateDistMesg *) vmesg;
  int		i;

  fprintf(fout,"{IMD\n");
  fprintf(fout,"ref:"F_IID"\n",mesg->refines);
  fprintf(fout,"mea:%.3f\n",mesg->mean);
  fprintf(fout,"std:%.3f\n",mesg->stddev);
  fprintf(fout,"min:"F_COORD"\n",mesg->min);
  fprintf(fout,"max:"F_COORD"\n",mesg->max);
  fprintf(fout,"buc:"F_S32"\n",mesg->num_buckets);
  fprintf(fout,"his:\n");
  for (i=0; i < mesg->num_buckets; ++i)
    fprintf(fout,F_S32"\n",mesg->histogram[i]);
  fprintf(fout,"}\n");
  return;
}

static void Write_ICM_Mesg(FILE *fout, void *vmesg)
{ IntConConMesg *mesg = (IntConConMesg *) vmesg;
  int		i;

  fprintf(fout,"{ICM\n");
  fprintf(fout,"acc:"F_IID"\n",mesg->iaccession);
  fprintf(fout,"pla:%c\n",mesg->placed);
  fprintf(fout,"len:"F_COORD"\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:"F_S32"\n",mesg->forced);
  fprintf(fout,"npc:"F_S32"\n",mesg->num_pieces);
  fprintf(fout,"nou:"F_S32"\n",mesg->num_unitigs);
  fprintf(fout,"nvr:"F_S32"\n",mesg->num_vars);
  for (i=0; i < mesg->num_vars; ++i)
    Write_IMV_Mesg(fout, &mesg->v_list[i]);
  for (i=0; i < mesg->num_pieces; ++i)
    Write_IMP_Mesg(fout, &mesg->pieces[i]);
  for (i=0; i < mesg->num_unitigs; ++i)
    Write_IUP_Mesg(fout, &(mesg->unitigs[i]));
  fprintf(fout,"}\n");
  return;
}


static void Write_IAF_Mesg(FILE *fout, void *vmesg)
{ IntAugFragMesg *mesg = (IntAugFragMesg *) vmesg;
  
  fprintf(fout,"{IAF\n");
  fprintf(fout,"acc:"F_IID"\n",mesg->iaccession);
  fprintf(fout,"typ:%c\n",(char) mesg->type);
  fprintf(fout,"chi:"F_S32"\n",mesg->chimeric);
  fprintf(fout,"cha:"F_S32"\n",mesg->chaff);
  fprintf(fout,"clr:"F_COORD","F_COORD"\n", mesg->clear_rng.bgn,mesg->clear_rng.end);
  fprintf(fout,"mst:%c\n",mesg->mate_status);
  fprintf(fout,"}\n");
  return;
}

/* Genome Snapshot output routines */
/***********************************/


static void Write_UPS_Mesg(FILE *fout, UnitigPos *up)
{ int i;

  fprintf(fout,"{UPS\n");
  fprintf(fout,"typ:%c\n",(char) up->type);
  fprintf(fout,"lid:"F_UID"\n",up->eident);
  fprintf(fout,"pos:"F_COORD","F_COORD"\n",up->position.bgn,up->position.end);
  fprintf(fout,"dln:"F_S32"\n",up->delta_length);
  fprintf(fout,"del:\n");
  if (up->delta_length > 0 ) {
    for(i=0; i < up->delta_length; i++)
      fprintf(fout,F_S32"%c",up->delta[i], (i%20 == 19) ? '\n' : ' ');
    if (up->delta_length%20 != 0) fprintf(fout,"\n");
  }
  fprintf(fout,"}\n");
  return;
}
static void Write_MPS_Mesg(FILE *fout, SnapMultiPos *mlp)
{ int i;

  fprintf(fout,"{MPS\n");
  fprintf(fout,"typ:%c\n",(char) mlp->type);
  fprintf(fout,"mid:"F_UID"\n",mlp->eident);
#ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mlp->source,FALSE);
#endif
  fprintf(fout,"pos:"F_COORD","F_COORD"\n",
          mlp->position.bgn,mlp->position.end);
  fprintf(fout,"dln:"F_S32"\n",mlp->delta_length);
  fprintf(fout,"del:\n");
  if (mlp->delta_length > 0 ) {
    for(i=0; i < mlp->delta_length; i++)
      fprintf(fout,F_S32"%c",mlp->delta[i], (i%20 == 19) ? '\n' : ' ');
    if (mlp->delta_length%20 != 0) fprintf(fout,"\n");
  }
  fprintf(fout,"}\n");
  return;
}


static void Write_UTG_Mesg(FILE *fout, void *vmesg)
{ SnapUnitigMesg *mesg = (SnapUnitigMesg *) vmesg;
  int			i;

  fprintf(fout,"{UTG\n");
  fprintf(fout,"acc:("F_UID","F_IID")\n",
          mesg->eaccession,mesg->iaccession);
#ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mesg->source,FALSE);
#endif
  fprintf(fout,"cov:%.3f\n",mesg->coverage_stat);
  fprintf(fout,"sta:%c\n",mesg->status);
  fprintf(fout,"abp:"F_COORD"\n",mesg->a_branch_point);
  fprintf(fout,"bbp:"F_COORD"\n",mesg->b_branch_point);
  fprintf(fout,"len:"F_COORD"\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:"F_S32"\n",mesg->forced);
  fprintf(fout,"nfr:"F_S32"\n",mesg->num_frags);
  for (i=0; i < mesg->num_frags; ++i)
    Write_MPS_Mesg(fout,&(mesg->f_list[i]));
  fprintf(fout,"}\n");
  return;
}


static void Write_ULK_Mesg(FILE *fout, void *vmesg)
{ SnapUnitigLinkMesg *mesg = (SnapUnitigLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{ULK\n");
  fprintf(fout,"ut1:"F_UID"\n",mesg->eunitig1);
  fprintf(fout,"ut2:"F_UID"\n",mesg->eunitig2);
  fprintf(fout,"ori:%c\n",mesg->orientation);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  fprintf(fout,"ipc:"F_S32"\n",mesg->is_possible_chimera);
  fprintf(fout,"gui:"F_S32"\n",mesg->includes_guide);
  fprintf(fout,"mea:%.3f\n",mesg->mean_distance);
  fprintf(fout,"std:%.3f\n",mesg->std_deviation);
  fprintf(fout,"num:"F_S32"\n",mesg->num_contributing);
  fprintf(fout,"sta:%c\n",mesg->status);
  npairs = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --npairs;
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout, F_UID","F_UID",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}


static void Write_CCO_Mesg(FILE *fout, void *vmesg)
{ SnapConConMesg *mesg = (SnapConConMesg *) vmesg;
  int		i;

  fprintf(fout,"{CCO\n");
  fprintf(fout,"acc:("F_UID","F_IID")\n",mesg->eaccession,mesg->iaccession);
  fprintf(fout,"pla:%c\n",mesg->placed);
  fprintf(fout,"len:"F_COORD"\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:"F_S32"\n",mesg->forced);
  fprintf(fout,"npc:"F_S32"\n",mesg->num_pieces);
  fprintf(fout,"nou:"F_S32"\n",mesg->num_unitigs);
  fprintf(fout,"nvr:"F_S32"\n",mesg->num_vars);
  for (i=0; i < mesg->num_vars; ++i)
    Write_VAR_Mesg(fout, &(mesg->vars[i]));
  for (i=0; i < mesg->num_pieces; ++i)
    Write_MPS_Mesg(fout, &mesg->pieces[i]);
  for (i=0; i < mesg->num_unitigs; ++i)
    Write_UPS_Mesg(fout, &(mesg->unitigs[i]));
  fprintf(fout,"}\n");
  return;
}


static void Write_CLK_Mesg(FILE *fout, void *vmesg)
{ SnapContigLinkMesg *mesg = (SnapContigLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{CLK\n");
  fprintf(fout,"co1:"F_UID"\n",mesg->econtig1);
  fprintf(fout,"co2:"F_UID"\n",mesg->econtig2);
  fprintf(fout,"ori:%c\n",mesg->orientation);
  fprintf(fout,"ovt:%c\n",mesg->overlap_type);
  fprintf(fout,"ipc:"F_S32"\n",mesg->is_possible_chimera);
  fprintf(fout,"gui:"F_S32"\n",mesg->includes_guide);
  fprintf(fout,"mea:%.3f\n",mesg->mean_distance);
  fprintf(fout,"std:%.3f\n",mesg->std_deviation);
  fprintf(fout,"num:"F_S32"\n",mesg->num_contributing);
  fprintf(fout,"sta:%c\n",mesg->status);
  npairs = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --npairs;
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout, F_UID","F_UID",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_SLK_Mesg(FILE *fout, void *vmesg)
{ SnapScaffoldLinkMesg *mesg = (SnapScaffoldLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{SLK\n");
  fprintf(fout,"sc1:"F_UID"\n",mesg->escaffold1);
  fprintf(fout,"sc2:"F_UID"\n",mesg->escaffold2);
  fprintf(fout,"ori:%c\n",mesg->orientation);
  fprintf(fout,"gui:"F_S32"\n",mesg->includes_guide);
  fprintf(fout,"mea:%.3f\n",mesg->mean_distance);
  fprintf(fout,"std:%.3f\n",mesg->std_deviation);
  fprintf(fout,"num:"F_S32"\n",mesg->num_contributing);
  npairs = mesg->num_contributing;
  assert(npairs > 0);
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout, F_UID","F_UID",%c\n",
            mesg->jump_list[i].in1,
            mesg->jump_list[i].in2,
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}


static void Write_CTP_Mesg(FILE *fout, SnapContigPairs *mesg)
{
  fprintf(fout,"{CTP\n");
  fprintf(fout,"ct1:"F_UID"\n",mesg->econtig1);
  fprintf(fout,"ct2:"F_UID"\n",mesg->econtig2);
  fprintf(fout,"mea:%.3f\n",mesg->mean);
  fprintf(fout,"std:%.3f\n",mesg->stddev);
  fprintf(fout,"ori:%c\n",mesg->orient);
  fprintf(fout,"}\n");
  return;
}

static void Write_SCF_Mesg(FILE *fout, void *vmesg)
{ SnapScaffoldMesg *mesg = (SnapScaffoldMesg *) vmesg;
  int		i;
  int num = MAX(1,mesg->num_contig_pairs);
  fprintf(fout,"{SCF\n");
  fprintf(fout,"acc:("F_UID","F_IID")\n",mesg->eaccession,mesg->iaccession);
  fprintf(fout,"noc:"F_S32"\n",mesg->num_contig_pairs);
  for (i=0; i < num; ++i)
    Write_CTP_Mesg(fout,&mesg->contig_pairs[i]);
  fprintf(fout,"}\n");
  return;
}

static void Write_DSC_Mesg(FILE *fout, void *vmesg)
{ SnapDegenerateScaffoldMesg *mesg = (SnapDegenerateScaffoldMesg *) vmesg;
  fprintf(fout,"{DSC\n");
  fprintf(fout,"acc:"F_UID"\n",mesg->eaccession);
  fprintf(fout,"ctg:"F_UID"\n",mesg->econtig);
  fprintf(fout,"}\n");
  return;
}

static void Write_IDS_Mesg(FILE *fout, void *vmesg)
{ IntDegenerateScaffoldMesg *mesg = (IntDegenerateScaffoldMesg *) vmesg;
  fprintf(fout,"{IDS\n");
  fprintf(fout,"ctg:"F_IID"\n",mesg->icontig);
  fprintf(fout,"}\n");
  return;
}


static void Write_MDI_Mesg(FILE *fout, void *vmesg)
{ SnapMateDistMesg *mesg = (SnapMateDistMesg *) vmesg;
  int		i;

  fprintf(fout,"{MDI\n");
  fprintf(fout,"ref:("F_UID","F_IID")\n",mesg->erefines,mesg->irefines);
  fprintf(fout,"mea:%.3f\n",mesg->mean);
  fprintf(fout,"std:%.3f\n",mesg->stddev);
  fprintf(fout,"min:"F_COORD"\n",mesg->min);
  fprintf(fout,"max:"F_COORD"\n",mesg->max);
  fprintf(fout,"buc:"F_S32"\n",mesg->num_buckets);
  fprintf(fout,"his:\n");
  for (i=0; i < mesg->num_buckets; ++i)
    fprintf(fout,F_S32"\n",mesg->histogram[i]);
  fprintf(fout,"}\n");
  return;
}

static void Write_BAT_Mesg(FILE *fout, void *vmesg){
  BatchMesg *mesg = (BatchMesg *)vmesg;
  fprintf(fout,"{BAT\n");
  fprintf(fout,"bna:%s\n",mesg->name);
  fprintf(fout,"crt:0\n");
  fprintf(fout,"acc:"F_UID"\n",mesg->eaccession);
  PutText(fout,"com:",mesg->comment, FALSE);
  fprintf(fout,"}\n");

}

static void Write_EOF_Mesg(FILE *fout, void *vmesg)
{
  EndOfFileMesg *mesg = (EndOfFileMesg *) vmesg;

  fprintf(fout,"{EOF\n");
  fprintf(fout,"sta:"F_S32"\n", mesg->status);
  fprintf(fout,"crt:0\n");
  PutText(fout,"com:", mesg->comment, FALSE);
  fprintf(fout,"}\n");
}




static AS_MSG_callrecord CallTable1[NUM_OF_REC_TYPES + 1] = {
  {"", NULL, NULL, 0l},
  {"{ADT", Read_ADT_Mesg, Write_ADT_Mesg, sizeof(AuditMesg) },
  {"{VER", Read_VER_Mesg, Write_VER_Mesg, sizeof(VersionMesg)  },
  {"{FRG", Read_FRG_Mesg, Write_FRG_Mesg, sizeof(FragMesg)  },
  {"{IFG", Read_IFG_Mesg, Write_IFG_Mesg, sizeof(InternalFragMesg) },
  {"{OFG", Read_OFG_Mesg, Write_OFG_Mesg, sizeof(OFGMesg) },
  {"{LKG", Read_LKG_Mesg, Write_LKG_Mesg, sizeof(LinkMesg) },
  {"", NULL, NULL, 0l },
  {"{DST", Read_DST_Mesg, Write_DST_Mesg, sizeof(DistanceMesg) },
  {"{IDT", Read_IDT_Mesg, Write_IDT_Mesg, sizeof(InternalDistMesg) },
  {"RLIB", NULL, NULL, 0l },  //  RESERVED for Version 2's LIB message
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"{OVL", Read_OVL_Mesg, Write_OVL_Mesg, sizeof(OverlapMesg) },
  {"", NULL, NULL, 0l },
  {"{UOM", Read_UOM_Mesg, Write_UOM_Mesg, sizeof(UnitigOverlapMesg) },
  {"{IUM", Read_IUM_Mesg, Write_IUM_Mesg, sizeof(IntUnitigMesg) },
  {"{IUL", Read_IUL_Mesg, Write_IUL_Mesg, sizeof(IntUnitigLinkMesg) },
  {"{ICL", Read_ICL_Mesg, Write_ICL_Mesg, sizeof(IntContigLinkMesg) },
  {"{AFG", Read_AFG_Mesg, Write_AFG_Mesg, sizeof(AugFragMesg) },
  {"{ISF", Read_ISF_Mesg, Write_ISF_Mesg, sizeof(IntScaffoldMesg) },
  {"{IMD", Read_IMD_Mesg, Write_IMD_Mesg, sizeof(IntMateDistMesg) },
  {"{IAF", Read_IAF_Mesg, Write_IAF_Mesg, sizeof(IntAugFragMesg) },
  {"{UTG", Read_UTG_Mesg, Write_UTG_Mesg, sizeof(SnapUnitigMesg) },
  {"{ULK", Read_ULK_Mesg, Write_ULK_Mesg, sizeof(SnapUnitigLinkMesg) },
  {"{ICM", Read_ICM_Mesg, Write_ICM_Mesg, sizeof(IntConConMesg) },
  {"{CCO", Read_CCO_Mesg, Write_CCO_Mesg, sizeof(SnapConConMesg) },
  {"{CLK", Read_CLK_Mesg, Write_CLK_Mesg, sizeof(SnapContigLinkMesg) },
  {"{SCF", Read_SCF_Mesg, Write_SCF_Mesg, sizeof(SnapScaffoldMesg) },
  {"{MDI", Read_MDI_Mesg, Write_MDI_Mesg, sizeof(SnapMateDistMesg) },
  {"{BAT", Read_BAT_Mesg, Write_BAT_Mesg, sizeof(BatchMesg) },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"{IDS", Read_IDS_Mesg, Write_IDS_Mesg, sizeof(IntDegenerateScaffoldMesg) },
  {"{DSC", Read_DSC_Mesg, Write_DSC_Mesg, sizeof(SnapDegenerateScaffoldMesg) },
  {"{SLK", Read_SLK_Mesg, Write_SLK_Mesg, sizeof(SnapScaffoldLinkMesg) },
  {"{ISL", Read_ISL_Mesg, Write_ISL_Mesg, sizeof(InternalScaffoldLinkMesg) },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"", NULL, NULL, 0l },
  {"{EOF", Read_EOF_Mesg, Write_EOF_Mesg, sizeof(EndOfFileMesg) }
};


void AS_MSG_setFormatVersion1(void) {
  memcpy(AS_MSG_globals->CallTable, CallTable1, sizeof(AS_MSG_callrecord) * (NUM_OF_REC_TYPES + 1));
}

