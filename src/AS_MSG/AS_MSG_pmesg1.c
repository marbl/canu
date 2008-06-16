
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
static char CM_ID[]= "$Id: AS_MSG_pmesg1.c,v 1.26 2008-06-16 22:53:26 brianwalenz Exp $";

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
char
GetIIDIIDMatePairType(AS_IID *IID1, AS_IID *IID2, FILE *fin) {

  ReadLine(fin,TRUE);

  *IID1 = 0;
  *IID2 = 0;

  char *str = AS_MSG_globals->curLine;

  *IID1 = strtoul(str, NULL, 10);

  while (*str != ',')  str++;
  str++;

  *IID2 = strtoul(str, NULL, 10);

  while (*str != ',')  str++;
  str++;

  return(*str);
}


static
char
GetUIDUIDMatePairType(AS_UID *UID1, AS_UID *UID2, FILE *fin) {

  ReadLine(fin,TRUE);

  (*UID1)  = AS_UID_undefined();
  (*UID2)  = AS_UID_undefined();

  char *str = AS_MSG_globals->curLine; 
  char *currLoc = str;

   // get first UID
   while (*currLoc != ',') { currLoc++; }
   *currLoc = '\0';
   (*UID1) = AS_UID_lookup(str, NULL);
   
   // get second UID
   str = ++currLoc;
   while (*currLoc != ',') { currLoc++; }
   *currLoc = '\0';
   (*UID2) = AS_UID_lookup(str, NULL);

   // return the type value
  return(*(++currLoc));
}





/******************** INPUT ROUTINES ***************************/

static
void *
Read_DST_Mesg(FILE *fin) {
  static DistanceMesg dmesg;

  dmesg.action = (ActionType)GetType("act:%c","action", fin);

  dmesg.eaccession = GetUID("acc:",fin);

  if (dmesg.action == 'R')
    dmesg.action = AS_UPDATE;

  if ((dmesg.action == AS_ADD) ||
      (dmesg.action == AS_UPDATE) ||
      (dmesg.action == AS_IGNORE)) {
    GET_FIELD(dmesg.mean,   "mea:%f", "mean field");
    GET_FIELD(dmesg.stddev, "std:%f", "stddev field");
  }

  GetEOM(fin);

  return(&dmesg);
}


static
void *
Read_VER_Mesg(FILE *fin) {
  static VersionMesg vmesg;

  GET_FIELD(vmesg.version,"ver:"F_U32,"version field");
  GetEOM(fin);

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


static
void *
Read_Frag_Mesg(FILE *fin, int frag_class) {
  static FragMesg fmesg;

  assert(frag_class == MESG_FRG);

  fmesg.version        = 1;
  fmesg.library_uid    = AS_UID_undefined();
  fmesg.library_iid    = 0;
  fmesg.plate_uid      = AS_UID_undefined();
  fmesg.plate_location = 0;
  fmesg.is_random      = 1;
  fmesg.status_code    = 'G';
  fmesg.clear_vec.bgn  = 1;  //  Format 1 cannot have vec or qlt; these disable them
  fmesg.clear_vec.end  = 0;
  fmesg.clear_qlt.bgn  = 1;
  fmesg.clear_qlt.end  = 0;
  
  fmesg.action = (ActionType)GetType("act:%c","action", fin);

  fmesg.eaccession = GetUID("acc:",fin);

  fmesg.source   = NULL;
  fmesg.sequence = NULL;
  fmesg.quality  = NULL;
  fmesg.hps      = NULL;

  if ((fmesg.action == AS_ADD) || (fmesg.action == AS_IGNORE)) { 

    // We want to succeed on all reads, and let the gatekeeper do its stuff
    fmesg.type = (FragType)GetType("typ:%c","type", fin);

    fmesg.source   = (char *) GetText("src:",fin,FALSE);

    ReadLine(fin, TRUE);  //  unused "entry time field" etm:

    fmesg.sequence = (char *) GetText("seq:",fin,TRUE);
    fmesg.quality  = (char *) GetText("qlt:",fin,TRUE);

    GET_PAIR(fmesg.clear_rng.bgn,fmesg.clear_rng.end,"clr:"F_COORD","F_COORD,"clear range field");

  }  //  action is AS_ADD or AS_IGNORE
  GetEOM(fin);
  return(&fmesg);
}

static void *Read_FRG_Mesg(FILE *fin)
{ return Read_Frag_Mesg(fin,MESG_FRG); } 



static void *Read_OVL_Mesg(FILE *fin)
{ static OverlapMesg omesg;
  int    idx;
  
  GET_FIELD(omesg.aifrag,"afr:"F_IID,"a-fragment field");
  GET_FIELD(omesg.bifrag,"bfr:"F_IID,"b-fragment field");

  omesg.orientation = (OrientType)GetType("ori:%1[NAIO]","orientation", fin);
  omesg.overlap_type = (OverlapType)GetType("olt:%1[DCSXdc]","overlap", fin);

  GET_FIELD(omesg.ahg,"ahg:"F_COORD,"a-hang field");
  GET_FIELD(omesg.bhg,"bhg:"F_COORD,"b-hang field");
  GET_FIELD(omesg.quality,"qua:%f","quality field");
  GET_FIELD(omesg.min_offset,"mno:"F_COORD,"min-offset field");
  GET_FIELD(omesg.max_offset,"mxo:"F_COORD,"max-offset field");
  GET_FIELD(omesg.polymorph_ct,"pct:"F_S32,"poly-count field");

  if (strncmp(ReadLine(fin,TRUE),"del:",4) != 0)
    MgenError("delta tag label");

  omesg.delta = (signed char *)GetMemory(2*AS_FRAG_MAX_LEN);

  {
    int i, n;     /* Read a delta item (only one of its kind) */
    char *t, *u;

    i = 0;
    while ((t = ReadLine(fin,TRUE))[0] != '.')
      while (1) {
        n = strtol(t,&u,10);
        if (u == t) break;
        t = u;
        if (! isspace((int)*t))
          MgenError("Delta is not a sequence of digits");
        omesg.delta[i++] = n;
      }
    omesg.delta[i] = 0;
  }

  GetEOM(fin);
  return(&omesg);
}

static
void *
Read_LKG_Mesg(FILE *fin) {
  static LinkMesg lmesg;

  lmesg.action = (ActionType)GetType("act:%c","action", fin);
  lmesg.type   = (LinkType)  GetType("typ:%c","link", fin);

  lmesg.frag1 = GetUID("fg1:",fin);
  lmesg.frag2 = GetUID("fg2:",fin);

  if ((lmesg.action == AS_ADD) || (lmesg.action == AS_IGNORE)) {
    ReadLine(fin, TRUE);  //  unused "entry time field" etm:
    lmesg.distance    = GetUID("dst:",fin);
    lmesg.link_orient = (OrientType)GetType("ori:%1[NAIOU]","link orientation", fin);
  }
  GetEOM(fin);
  return(&lmesg);
}

static void *Read_UOM_Mesg(FILE *fin)
{ static UnitigOverlapMesg	mesg;
  GET_FIELD(mesg.chunk1,"ck1:"F_IID,"chunk 1 id field");
  GET_FIELD(mesg.chunk2,"ck2:"F_IID,"chunk 2 id field");
  
  mesg.orient = (ChunkOrientationType)GetType("ori:%1[NAIO]","orientation", fin);
  mesg.overlap_type = (UnitigOverlapType)GetType("ovt:%1[NOTCIMXdcYZ]","overlap type", fin);

#ifdef AS_ENABLE_SOURCE
  mesg.source = GetText("src:",fin,FALSE);
#endif
  GET_FIELD(mesg.best_overlap_length,"len:"F_COORD,"best overlap");
  GET_FIELD(mesg.min_overlap_length,"min:"F_COORD,"min overlap");
  GET_FIELD(mesg.max_overlap_length,"max:"F_COORD,"max overlap");
  GET_FIELD(mesg.quality,"qua:%f","quality field");
  GetEOM(fin);
  return(&mesg);
}

static
void
Read_IMP_Mesg(FILE *fin, IntMultiPos *imp) {
  int		 i;
  int32		 n, *delta;
  char		*line, *u;
  
  imp->type = (FragType)GetType("typ:%1[RXTELUFSUcBCG]","multipos$", fin);
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
  if (strncmp(ReadLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  imp->delta = NULL;
  if (imp->delta_length > 0) {
    imp->delta = (int32 *)GetMemory(sizeof(int32) * imp->delta_length);
    i = 0;
    while (i < imp->delta_length) {
      line = ReadLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	imp->delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  GetEOM(fin);
}

static
void
Read_IMV_Mesg(FILE *fin, IntMultiVar *imv) {

  GET_PAIR(imv->position.bgn,imv->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(imv->num_reads,"nrd:"F_S32,"number of reads"); 
  GET_FIELD(imv->num_conf_alleles,"nca:"F_S32,"number of confirmed alleles");
  GET_FIELD(imv->min_anchor_size,"anc:"F_S32,"minimal anchor size");
  GET_FIELD(imv->var_length,"len:"F_S32,"length field");
  GET_FIELD(imv->curr_var_id,"vid:"F_S32,"current VAR record id");
  GET_FIELD(imv->phased_var_id,"pid:"F_S32,"phased VAR record id");
  imv->nr_conf_alleles = GetText("nra:",fin,FALSE);
  imv->weights         = GetText("wgt:",fin,FALSE);
  imv->var_seq         = GetText("seq:",fin,FALSE);
  imv->conf_read_iids  = GetText("rid:",fin,FALSE);
  GetEOM(fin);
}

static
void
Read_VAR_Mesg(FILE *fin, IntMultiVar *smv) {

  GET_PAIR(smv->position.bgn,smv->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(smv->num_reads,"nrd:"F_S32,"number of reads");
  GET_FIELD(smv->num_conf_alleles,"nca:"F_S32,"number of confirmed alleles");
  GET_FIELD(smv->min_anchor_size,"anc:"F_S32,"minimal anchor size");
  GET_FIELD(smv->var_length,"len:"F_S32,"length field");
  GET_FIELD(smv->curr_var_id,"vid:"F_S32,"current VAR record id");
  GET_FIELD(smv->phased_var_id,"pid:"F_S32,"phased VAR record id");
  smv->nr_conf_alleles = GetText("nra:",fin,FALSE);
  smv->weights         = GetText("wgt:",fin,FALSE);
  smv->var_seq         = GetText("seq:",fin,FALSE);
  smv->conf_read_iids  = GetText("rid:",fin,FALSE);
  GetEOM(fin);
}

static
void
Read_IUP_Mesg(FILE *fin, IntUnitigPos *iup) {
  int			i;
  int32			n, *delta;
  char			*line, *u;
  
  iup->type = (UnitigType)GetType("typ:%1[URSPsX]","unitigpos type", fin);
  GET_FIELD(iup->ident,"lid:"F_IID,"unitigpos id");
  GET_PAIR(iup->position.bgn,iup->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(iup->delta_length,"dln:"F_S32,"delta length");
  if (strncmp(ReadLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  iup->delta = NULL;
  if (iup->delta_length > 0) {
    iup->delta = (int32 *)GetMemory(sizeof(int32)*iup->delta_length);
    i = 0;
    while (i < iup->delta_length) {
      line = ReadLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	iup->delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  GetEOM(fin);
}

static
void *
Read_IUM_Mesg(FILE *fin) {
  static IntUnitigMesg		mesg;
  int				i;

  GET_FIELD(mesg.iaccession,"acc:"F_IID,"accession field");
# ifdef AS_ENABLE_SOURCE
  mesg.source = GetText("src:",fin,FALSE);
# endif
  GET_FIELD(mesg.coverage_stat,"cov:%f","coverage stat");  
  mesg.status = (UnitigStatus)GetType("sta:%1[UCNSX]","status", fin);

  // flag for handling unitig
  mesg.unique_rept = (UnitigFUR)GetType("fur:%1[XUR]","unique_rept", fin);

  ReadLine(fin, TRUE);  //  unused "a branch point" abp:
  ReadLine(fin, TRUE);  //  unused "b branch point" bbp:

  GET_FIELD(mesg.length,"len:"F_COORD,"length field");
  mesg.consensus = GetText("cns:",fin,TRUE);
  mesg.quality   = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:"F_S32,"forced booleon");
  GET_FIELD(mesg.num_frags,"nfr:"F_S32,"num frags field");

  mesg.f_list = NULL;

  if (mesg.num_frags > 0) {
    mesg.f_list = (IntMultiPos *)GetMemory(mesg.num_frags * sizeof(IntMultiPos));

    for (i=0; i < mesg.num_frags; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{IMP",4) != 0)
	MgenError("Expecting IMP record");
      Read_IMP_Mesg(fin, mesg.f_list + i);
    }
  }

  GetEOM(fin);

  assert(strlen(mesg.consensus) == strlen(mesg.quality));

  return(&mesg);
}

static
void *
Read_IUL_Mesg(FILE *fin) {
  static IntUnitigLinkMesg	mesg;
  int				i,size;
  
  GET_FIELD(mesg.unitig1,"ut1:"F_IID,"unitig 1 field");
  GET_FIELD(mesg.unitig2,"ut2:"F_IID,"unitig 2 field");
  mesg.orientation = (ChunkOrientationType)GetType("ori:%1[NAOI]","orientation", fin);
  mesg.overlap_type = (UnitigOverlapType)GetType("ovt:%1[NOTCIMXYZ]","overlap type", fin);
  GET_FIELD(mesg.is_possible_chimera,"ipc:"F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  mesg.status = (PlacementStatusType)GetType("sta:%1[APBCU]","placement status", fin);
  if (strncmp(ReadLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  mesg.jump_list = NULL;
  if (size > 0) {
    mesg.jump_list = (IntMate_Pairs *)GetMemory(sizeof(IntMate_Pairs)*size);
    for (i=0; i < size; ++i) {
      IntMate_Pairs *imp = mesg.jump_list + i;
      //GET_TRIPLE(imp->in1,imp->in2,ch,F_IID","F_IID",%1[MSBRYT]","mate pair");
      imp->type = (LinkType)GetIIDIIDMatePairType(&imp->in1, &imp->in2, fin);;
    }
  }
  GetEOM(fin);
  return(&mesg);
}

static
void *
Read_ICL_Mesg(FILE *fin) {
  static IntContigLinkMesg	mesg;
  int				i,size;
  
  GET_FIELD(mesg.contig1,"co1:"F_IID,"contig 1 field");
  GET_FIELD(mesg.contig2,"co2:"F_IID,"contig 2 field");
  mesg.orientation = (ChunkOrientationType)GetType("ori:%1[NAOI]","orientation", fin);
  mesg.overlap_type = (UnitigOverlapType)GetType("ovt:%1[NOTCIMXYZ]","overlap type", fin);
  GET_FIELD(mesg.is_possible_chimera,"ipc:"F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  mesg.status = (PlacementStatusType)GetType("sta:%1[APBCU]","placement status", fin);
  if (strncmp(ReadLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  mesg.jump_list = NULL;
  if (size > 0) {
    mesg.jump_list = (IntMate_Pairs *)GetMemory(sizeof(IntMate_Pairs)*size);
    for (i=0; i < size; ++i) {
      IntMate_Pairs *imp = mesg.jump_list + i;
      //GET_TRIPLE(imp->in1,imp->in2,ch, F_IID","F_IID",%1[MSBRYT]","mate pair");
      imp->type = (LinkType)GetIIDIIDMatePairType(&imp->in1, &imp->in2, fin);
    }
  }
  GetEOM(fin);
  return(&mesg);
}

static
void *
Read_ISL_Mesg(FILE *fin) {
  static InternalScaffoldLinkMesg	mesg;
  int				i,size;
  
  GET_FIELD(mesg.iscaffold1,"sc1:"F_IID,"scaffold 1 field");
  GET_FIELD(mesg.iscaffold2,"sc2:"F_IID,"scaffold 2 field");
  mesg.orientation = (ChunkOrientationType)GetType("ori:%1[NAOI]","orientation", fin);
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  if (strncmp(ReadLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  assert(size > 0);
  mesg.jump_list = (IntMate_Pairs *)GetMemory(sizeof(IntMate_Pairs)*size);
  for (i=0; i < size; ++i) {
    IntMate_Pairs *imp = mesg.jump_list + i;
    //GET_TRIPLE(imp->in1,imp->in2,ch, F_IID","F_IID",%1[MSBRYT]","mate pair");
    imp->type = (LinkType)GetIIDIIDMatePairType(&imp->in1, &imp->in2, fin);
  }
  GetEOM(fin);
  return(&mesg);
}

static void *Read_AFG_Mesg(FILE *fin)
{ static AugFragMesg		mesg;
  char *line;
  
  mesg.eaccession = GetUIDIID("acc:",&mesg.iaccession,fin);

  mesg.mate_status = (MateStatType)GetType("mst:%1[ZGCLSONHADEURF]","mate status", fin);

  GET_FIELD(mesg.chimeric_NOTUSED,"chi:"F_S32,"chimeric flag");
  GET_FIELD(mesg.chaff,"cha:"F_S32,"chaff flag");
  GET_PAIR(mesg.clear_rng.bgn,mesg.clear_rng.end,"clr:"F_COORD","F_COORD,"clear range");
  GetEOM(fin);
  return(&mesg);
}

static void *Read_AMP_Mesg(FILE *fin)
{ static AugMatePairMesg	mesg;
  
  mesg.fragment1 = GetUID("frg:",fin);
  mesg.fragment2 = GetUID("frg:",fin);
  mesg.mate_status = (MateStatType)GetType("mst:%1[ZGCLSONHADEURF]","mate status", fin);
  GetEOM(fin);
  return(&mesg);
}

static void Read_ICP_Mesg(FILE *fin, IntContigPairs *icp)
{
  GET_FIELD(icp->contig1,"ct1:"F_IID,"contig 1 id");
  GET_FIELD(icp->contig2,"ct2:"F_IID,"contig 2 id");
  GET_FIELD(icp->mean,"mea:%f","mean distance");
  GET_FIELD(icp->stddev,"std:%f","standard deviation");
  icp->orient = (ChunkOrientationType)GetType("ori:%1[NAIOU]","link orientation", fin);
  GetEOM(fin);
}

static
void *
Read_ISF_Mesg(FILE *fin) {
  static IntScaffoldMesg	mesg;
  int				i, num;
  IntContigPairs		*icp;

  GET_FIELD(mesg.iaccession,"acc:"F_IID,"ISF accession");
  GET_FIELD(mesg.num_contig_pairs,"noc:"F_S32,"number of contigs");
  num = MAX(1,mesg.num_contig_pairs);
  if (num > 0) {
    icp = mesg.contig_pairs = (IntContigPairs *)GetMemory(num*sizeof(IntContigPairs));
    for (i=0; i < num; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{ICP",4) != 0)
	MgenError("Expecting ICP record");
      Read_ICP_Mesg(fin,icp);
      ++icp;
    }
  }
  else
    mesg.contig_pairs = NULL;
  GetEOM(fin);
  return(&mesg);
}

static
void *
Read_IMD_Mesg(FILE *fin) {
  static IntMateDistMesg	mesg;
  int				i;

  GET_FIELD(mesg.refines,"ref:"F_IID,"distance id");
  GET_FIELD(mesg.mean,"mea:%f","mean distance");
  GET_FIELD(mesg.stddev,"std:%f","standard deviation");
  GET_FIELD(mesg.min,"min:"F_COORD,"min distance");
  GET_FIELD(mesg.max,"max:"F_COORD,"max distance");
  GET_FIELD(mesg.num_buckets,"buc:"F_S32,"number of buckets");
  if (strncmp(ReadLine(fin,TRUE),"his:",4) != 0)
    MgenError("Expecting his field");
  if (mesg.num_buckets > 0) {
    mesg.histogram = (int32 *)GetMemory(mesg.num_buckets*sizeof(int32));
    for (i=0; i < mesg.num_buckets; ++i)
      GET_FIELD(mesg.histogram[i],F_S32,"histogram entry");
  }
  else
    mesg.histogram = NULL;
  GetEOM(fin);
  return(&mesg);
}

static
void *
Read_ICM_Mesg(FILE *fin) {
  static IntConConMesg		mesg;
  int				i;
  
  GET_FIELD(mesg.iaccession,"acc:"F_IID,"accession number");
  mesg.placed = (ContigPlacementStatusType)GetType("pla:%1[PU]","placed flag", fin);
  GET_FIELD(mesg.length,"len:"F_COORD,"contig length");
  mesg.consensus = GetText("cns:",fin,TRUE);
  mesg.quality   = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:"F_S32,"forced flag");
  GET_FIELD(mesg.num_pieces,"npc:"F_S32,"number of pieces");
  GET_FIELD(mesg.num_unitigs,"nou:"F_S32,"number of unitigs");
  GET_FIELD(mesg.num_vars,"nvr:"F_S32,"num vars field");

  mesg.v_list = NULL;
  mesg.pieces = NULL;
  mesg.unitigs = NULL;

  if (mesg.num_vars > 0) {
    mesg.v_list = (IntMultiVar *)GetMemory(mesg.num_vars   *sizeof(IntMultiVar));
    for (i=0; i < mesg.num_vars; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{IMV",4) != 0)
        MgenError("Expecting IMV record");
      Read_IMV_Mesg(fin, mesg.v_list + i);
    }
  }

  if (mesg.num_pieces > 0) {
    mesg.pieces = (IntMultiPos *)GetMemory(mesg.num_pieces *sizeof(IntMultiPos));
    for (i=0; i < mesg.num_pieces; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{IMP",4) != 0)
        MgenError("Expecting IMP record");
      Read_IMP_Mesg(fin, mesg.pieces + i);
    }
  }  

  if (mesg.num_unitigs > 0) {
    mesg.unitigs = (IntUnitigPos *)GetMemory(mesg.num_unitigs*sizeof(IntUnitigPos));
    for (i=0; i < mesg.num_unitigs; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{IUP",4) != 0)
	MgenError("Expecting IUP record");
      Read_IUP_Mesg(fin, mesg.unitigs + i);
    }
  }

  GetEOM(fin);

  return(&mesg);
}


static void *Read_IAF_Mesg(FILE *fin)
{ static IntAugFragMesg		mesg;
  
  GET_FIELD(mesg.iaccession,"acc:"F_IID,"accession field");
  mesg.type = (FragType)GetType("typ:%1[RXELTFSUCBWG]","type", fin);
  GET_FIELD(mesg.chimeric_NOTUSED,"chi:"F_S32,"chimeric flag");
  GET_FIELD(mesg.chaff,"cha:"F_S32,"chaff flag");
  GET_PAIR(mesg.clear_rng.bgn,mesg.clear_rng.end,"clr:"F_COORD","F_COORD,"clear range");
  mesg.mate_status = (MateStatType)GetType("mst:%1[ZGCLSONHADEURF]","mate status", fin);
  GetEOM(fin);
  return(&mesg);
}


static void *Read_IAM_Mesg(FILE *fin)
{ static IntAugMatePairMesg	mesg;
  
  GET_FIELD(mesg.fragment1,"frg:"F_IID,"accession field");
  GET_FIELD(mesg.fragment2,"frg:"F_IID,"accession field");
  mesg.mate_status = (MateStatType)GetType("mst:%1[ZGCLSONHADEURF]","mate status", fin);
  GetEOM(fin);
  return(&mesg);
}


static
void *
Read_EOF_Mesg(FILE *fin) {
  static EndOfFileMesg mesg;
  time_t entry_time;

  GET_FIELD(mesg.status,"sta:"F_S32,"status field");
  ReadLine(fin, TRUE);  //  unused "entry time field" crt:
  mesg.comment = GetText("com:",fin,FALSE);
  GetEOM(fin);

  return (&mesg);
}



/* Genome snapshot input routines */
/**********************************/

static
void
Read_MPS_Mesg(FILE *fin, SnapMultiPos *imp) {
  int			i;
  int32			n;
  char			*line, *u;

  imp->type = (FragType)GetType("typ:%1[RXTEFUSLuBG]","multipos type", fin);

  imp->eident = GetUID("mid:",fin);

#ifdef AS_ENABLE_SOURCE
  imp->source = GetText("src:",fin,FALSE);
#endif
  GET_PAIR(imp->position.bgn,imp->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(imp->delta_length,"dln:"F_S32,"delta length");
  if (strncmp(ReadLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  imp->delta = NULL;
  if (imp->delta_length > 0) {
    imp->delta = (int32 *)GetMemory(sizeof(int32) * imp->delta_length);
    i = 0;
    while (i < imp->delta_length) {
      line = ReadLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	imp->delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  GetEOM(fin);
}

static
void
Read_UPS_Mesg(FILE *fin, UnitigPos *iup) {
  int			i;
  int32			n;
  char			*line, *u;
  
  iup->type = (UnitigType)GetType("typ:%1[URSPs]","unitigpos type", fin);

  iup->eident = GetUID("lid:",fin);

  GET_PAIR(iup->position.bgn,iup->position.end,"pos:"F_COORD","F_COORD,"position field");
  GET_FIELD(iup->delta_length,"dln:"F_S32,"delta length");
  if (strncmp(ReadLine(fin,TRUE),"del:",4) != 0)
    MgenError("Missing del: field");
  iup->delta = NULL;
  if (iup->delta_length > 0) {
    iup->delta = (int32 *)GetMemory(sizeof(int32)*iup->delta_length);
    i = 0;
    while (i < iup->delta_length) {
      line = ReadLine(fin,TRUE);
      n = (int32) strtol(line,&u,10);
      while (u != line) {
	line = u;
	iup->delta[i++] = n;
        n = (int32) strtol(line,&u,10);
      }
    }
  } 
  GetEOM(fin);
}


static
void
Read_CTP_Mesg(FILE *fin, SnapContigPairs *icp) {

  icp->econtig1 = GetUID("ct1:",fin);
  icp->econtig2 = GetUID("ct2:",fin);

  GET_FIELD(icp->mean,"mea:%f","mean distance");
  GET_FIELD(icp->stddev,"std:%f","standard deviation");
  icp->orient = (ChunkOrientationType)GetType("ori:%1[NAIOU]","link orientation", fin);
  GetEOM(fin);
}


static
void *
Read_UTG_Mesg(FILE *fin) {
  static SnapUnitigMesg		mesg;
  int				i;

  mesg.eaccession = GetUIDIID("acc:",&mesg.iaccession,fin);

#ifdef AS_ENABLE_SOURCE
  mesg.source = GetText("src:",fin,FALSE);
#endif
  GET_FIELD(mesg.coverage_stat,"cov:%f","coverage stat");
  mesg.status = (UnitigStatus)GetType("sta:%1[UCNSX]","status", fin);

  ReadLine(fin, TRUE);  //  unused "a branch point" abp:
  ReadLine(fin, TRUE);  //  unused "b branch point" bbp:

  GET_FIELD(mesg.length,"len:"F_COORD,"length field");
  mesg.consensus = GetText("cns:",fin,TRUE);
  mesg.quality   = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:"F_S32,"forced booleon");
  GET_FIELD(mesg.num_frags,"nfr:"F_S32,"num frags field");

  mesg.f_list = NULL;
  if (mesg.num_frags > 0) {
    mesg.f_list = (SnapMultiPos *)GetMemory(mesg.num_frags*sizeof(SnapMultiPos));

    for (i=0; i < mesg.num_frags; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{MPS",4) != 0)
	MgenError("Expecting MPS record");
      Read_MPS_Mesg(fin, mesg.f_list + i);
    }
  }

  GetEOM(fin);

  return(&mesg);
}


static void *Read_ULK_Mesg(FILE *fin) {
  static SnapUnitigLinkMesg	mesg;
  int				i,size;

  mesg.eunitig1 = GetUID("ut1:",fin);
  mesg.eunitig2 = GetUID("ut2:",fin);

  mesg.orientation = (ChunkOrientationType)GetType("ori:%1[NAOI]","orientation", fin);
  mesg.overlap_type = (UnitigOverlapType)GetType("ovt:%1[NOTCIMXYZ]","overlap type", fin);
  GET_FIELD(mesg.is_possible_chimera,"ipc:"F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  mesg.status = (PlacementStatusType)GetType("sta:%1[APBCU]","placement status", fin);
  if (strncmp(ReadLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;
  if (size > 0) {
    mesg.jump_list = (SnapMate_Pairs *)GetMemory(sizeof(SnapMate_Pairs)*size);
    for (i=0; i < size; ++i) {
      SnapMate_Pairs *imp = mesg.jump_list + i;
      imp->type = (LinkType)GetUIDUIDMatePairType(&imp->in1, &imp->in2, fin);  //  valid MSBRYT
    }
  }
  else
    mesg.jump_list = NULL;
  GetEOM(fin);
  return(&mesg);
}


static void *Read_CCO_Mesg(FILE *fin)
{ static SnapConConMesg		mesg;
  int  	 i;
  
  mesg.eaccession = GetUIDIID("acc:",&mesg.iaccession,fin);
  mesg.placed = (ContigPlacementStatusType)GetType("pla:%1[PU]","placed flag", fin);
  GET_FIELD(mesg.length,"len:"F_COORD,"contig length");
  mesg.consensus = GetText("cns:",fin,TRUE);
  mesg.quality   = GetText("qlt:",fin,TRUE);
  GET_FIELD(mesg.forced,"for:"F_S32,"forced flag");
  GET_FIELD(mesg.num_pieces,"npc:"F_S32,"number of pieces");
  GET_FIELD(mesg.num_unitigs,"nou:"F_S32,"number of unitigs");
  GET_FIELD(mesg.num_vars,"nvr:"F_S32,"number of vars");

  mesg.vars = NULL;
  mesg.pieces = NULL;
  mesg.unitigs = NULL;

  if (mesg.num_vars > 0) {
    mesg.vars = (IntMultiVar *)GetMemory(mesg.num_vars  *sizeof(IntMultiVar));
    for (i=0; i < mesg.num_vars; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{VAR",4) != 0)
        MgenError("Expecting VAR record");
      Read_VAR_Mesg(fin, mesg.vars + i);
    }
  }

  if (mesg.num_pieces > 0) {
    mesg.pieces = (SnapMultiPos *)GetMemory(mesg.num_pieces * sizeof(SnapMultiPos));
    for (i=0; i < mesg.num_pieces; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{MPS",4) != 0)
        MgenError("Expecting MPS record");
      Read_MPS_Mesg(fin, mesg.pieces + i);
    }
  }

  if (mesg.num_unitigs > 0) {
    mesg.unitigs  = (UnitigPos *)GetMemory(mesg.num_unitigs*sizeof(UnitigPos));
    for (i=0; i < mesg.num_unitigs; ++i) {
      if (strncmp(ReadLine(fin,TRUE),"{UPS",4) != 0)
	MgenError("Expecting UPS record");
      Read_UPS_Mesg(fin, mesg.unitigs + i);
    }
  }

  GetEOM(fin);

  return(&mesg);
}



static void *Read_CLK_Mesg(FILE *fin)
{ static SnapContigLinkMesg	mesg;
  int				i,size;
  
  mesg.econtig1 = GetUID("co1:",fin);
  mesg.econtig2 = GetUID("co2:",fin);

  mesg.orientation = (ChunkOrientationType)GetType("ori:%1[NAOI]","orientation", fin);
  mesg.overlap_type = (UnitigOverlapType)GetType("ovt:%1[NOTCIMXYZ]","overlap type", fin);
  GET_FIELD(mesg.is_possible_chimera,"ipc:"F_S32,"warning");
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  mesg.status = (PlacementStatusType)GetType("sta:%1[APBCU]","placement status", fin);
  if (strncmp(ReadLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");

  size = mesg.num_contributing;
  if (mesg.overlap_type != AS_NO_OVERLAP)
    --size;

  mesg.jump_list = NULL;
  if (size > 0) {
    mesg.jump_list = (SnapMate_Pairs *)GetMemory(sizeof(SnapMate_Pairs)*size);
    for (i=0; i < size; ++i) {
      SnapMate_Pairs *imp = mesg.jump_list + i;
      imp->type = (LinkType)GetUIDUIDMatePairType(&imp->in1, &imp->in2, fin);  //  valid MSBRYT
    }
  }

  GetEOM(fin);

  return(&mesg);
}

static void *Read_SLK_Mesg(FILE *fin)
{ static SnapScaffoldLinkMesg	mesg;
  int				i,size;
  
  mesg.escaffold1 = GetUID("sc1:",fin);
  mesg.escaffold2 = GetUID("sc2:",fin);

  mesg.orientation = (ChunkOrientationType)GetType("ori:%1[NAOI]","orientation", fin);
  GET_FIELD(mesg.includes_guide,"gui:"F_S32,"guide flag");
  GET_FIELD(mesg.mean_distance,"mea:%f","mean distance");
  GET_FIELD(mesg.std_deviation,"std:%f","standard deviation");
  GET_FIELD(mesg.num_contributing,"num:"F_S32,"number of links");
  if (strncmp(ReadLine(fin,TRUE),"jls:",4) != 0)
    MgenError("Expecting jls field");
  size = mesg.num_contributing;
  assert(size > 0) ;
  mesg.jump_list = (SnapMate_Pairs *)GetMemory(sizeof(SnapMate_Pairs)*size);
  for (i=0; i < size; ++i) {
    SnapMate_Pairs *imp = mesg.jump_list + i;
    imp->type = (LinkType)GetUIDUIDMatePairType(&imp->in1, &imp->in2, fin);  //  valid MSBRYT
  }

  GetEOM(fin);

  return(&mesg);
}


static void *Read_SCF_Mesg(FILE *fin)
{ static SnapScaffoldMesg	mesg;
  int				i, num;

  mesg.eaccession = GetUIDIID("acc:",&mesg.iaccession,fin);

  GET_FIELD(mesg.num_contig_pairs,"noc:"F_S32,"number of contigs");
  num = MAX(mesg.num_contig_pairs, 1);
  if (num > 0) {
    mesg.contig_pairs = (SnapContigPairs *)GetMemory(num * sizeof(SnapContigPairs));
    for (i=0; i < num; ++i) {
      SnapContigPairs *icp = mesg.contig_pairs + i;
      if (strncmp(ReadLine(fin,TRUE),"{CTP",4) != 0)
	MgenError("Expecting CTP record");
      Read_CTP_Mesg(fin,icp);
    }
  }
  else
    mesg.contig_pairs = NULL;
  GetEOM(fin);
  return(&mesg);
}


static void *Read_MDI_Mesg(FILE *fin)
{ static SnapMateDistMesg	mesg;
  int				i;

  mesg.erefines = GetUIDIID("ref:",&mesg.irefines,fin);
  GET_FIELD(mesg.mean,"mea:%f","mean distance");
  GET_FIELD(mesg.stddev,"std:%f","standard deviation");
  GET_FIELD(mesg.min,"min:"F_COORD,"min distance");
  GET_FIELD(mesg.max,"max:"F_COORD,"max distance");
  GET_FIELD(mesg.num_buckets,"buc:"F_S32,"number of buckets");
  if (strncmp(ReadLine(fin,TRUE),"his:",4) != 0)
    MgenError("Expecting his field");
  if (mesg.num_buckets > 0) {
    mesg.histogram = (int32 *)GetMemory(mesg.num_buckets*sizeof(int32));

    for (i=0; i < mesg.num_buckets; ++i)
      GET_FIELD(mesg.histogram[i],F_S32,"histogram entry");
  }
  else
    mesg.histogram = NULL;
  GetEOM(fin);
  return(&mesg);
}


static void *Read_BAT_Mesg(FILE *fin){
  static BatchMesg mesg;

  mesg.name       = GetString("bna:",fin);
  ReadLine(fin, TRUE);  //  unused "entry time" crt:
  mesg.eaccession = GetUID("acc:",fin);
  mesg.comment    = GetText("com:",fin, FALSE);
  GetEOM(fin);
  return(&mesg);
}

/******************** OUTPUT ROUTINES ***************************/


static void Write_DST_Mesg(FILE *fout, void *vmesg)
{ DistanceMesg *mesg = (DistanceMesg *) vmesg;

  fprintf(fout,"{DST\n");
  fprintf(fout,"act:%c\n",mesg->action);
  fprintf(fout,"acc:%s\n",AS_UID_toString(mesg->eaccession));
  if (mesg->action != AS_DELETE)
    { fprintf(fout,"mea:%.3f\n",mesg->mean);
      fprintf(fout,"std:%.3f\n",mesg->stddev);
    }
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
  fprintf(fout,"fg1:%s\n",AS_UID_toString(mesg->frag1));
  fprintf(fout,"fg2:%s\n",AS_UID_toString(mesg->frag2));
  if((mesg->action == AS_ADD) || (mesg->action == AS_IGNORE))
    { fprintf(fout,"etm:0\n");
      fprintf(fout,"dst:%s\n",AS_UID_toString(mesg->distance));
      fprintf(fout,"ori:%c\n",mesg->link_orient);
    }
  fprintf(fout,"}\n");
}

static void Write_Frag_Mesg(FILE *fout, void *vmesg, int frag_class) {
  FragMesg *mesg = (FragMesg *) vmesg;

  assert(frag_class == MESG_FRG);

  fprintf(fout,"{%s\n",MessageTypeName[frag_class]);
  fprintf(fout,"act:%c\n",mesg->action);
  if (frag_class == MESG_FRG)
    fprintf(fout,"acc:%s\n",AS_UID_toString(mesg->eaccession));
  else
    fprintf(fout,"acc:(%s,"F_IID")\n",AS_UID_toString(mesg->eaccession),mesg->iaccession);

  if ((mesg->action == AS_ADD) || (mesg->action == AS_IGNORE)) {
    fprintf(fout,"typ:%c\n",(char) mesg->type);
    PutText(fout,"src:",mesg->source,FALSE);
    fprintf(fout,"etm:0\n");
    PutText(fout,"seq:",mesg->sequence,TRUE);
    PutText(fout,"qlt:",mesg->quality,TRUE);
    fprintf(fout,"clr:"F_COORD","F_COORD"\n", mesg->clear_rng.bgn,mesg->clear_rng.end);
  }

  fprintf(fout,"}\n");
}

static void Write_FRG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,MESG_FRG); }


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
  fprintf(fout,"anc:"F_S32"\n",imv->min_anchor_size);
  fprintf(fout,"len:"F_S32"\n",imv->var_length);
  fprintf(fout,"vid:"F_S32"\n",imv->curr_var_id);
  fprintf(fout,"pid:"F_S32"\n",imv->phased_var_id);
  PutText(fout,"nra:",imv->nr_conf_alleles,FALSE);
  PutText(fout,"wgt:",imv->weights,FALSE);
  PutText(fout,"seq:",imv->var_seq,FALSE);
  PutText(fout,"rid:",imv->conf_read_iids,FALSE);
  fprintf(fout,"}\n");
  return;
}

static void Write_VAR_Mesg(FILE *fout, IntMultiVar *smv)
{
  fprintf(fout,"{VAR\n");
  fprintf(fout,"pos:"F_COORD","F_COORD"\n",smv->position.bgn,smv->position.end);
  fprintf(fout,"nrd:"F_S32"\n",smv->num_reads);
  fprintf(fout,"nca:"F_S32"\n",smv->num_conf_alleles);
  fprintf(fout,"anc:"F_S32"\n",smv->min_anchor_size);
  fprintf(fout,"len:"F_S32"\n",smv->var_length);
  fprintf(fout,"vid:"F_S32"\n",smv->curr_var_id);
  fprintf(fout,"pid:"F_S32"\n",smv->phased_var_id);
  PutText(fout,"nra:",smv->nr_conf_alleles,FALSE);
  PutText(fout,"wgt:",smv->weights,FALSE);
  PutText(fout,"seq:",smv->var_seq,FALSE);
  PutText(fout,"rid:",smv->conf_read_iids,FALSE);
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
  fprintf(fout,"fur:%c\n",mesg->unique_rept);
  fprintf(fout,"abp:0\n");  //  Unused!
  fprintf(fout,"bbp:0\n");  //  Unused!
  fprintf(fout,"len:"F_COORD"\n",mesg->length);
  PutText(fout,"cns:",mesg->consensus,TRUE);
  PutText(fout,"qlt:",mesg->quality,TRUE);
  fprintf(fout,"for:"F_S32"\n",mesg->forced);
  fprintf(fout,"nfr:"F_S32"\n",mesg->num_frags);
  assert(mesg->num_frags > 0);
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
  fprintf(fout,"acc:(%s,"F_IID")\n",AS_UID_toString(mesg->eaccession),mesg->iaccession);
  fprintf(fout,"mst:%c\n",mesg->mate_status);
  fprintf(fout,"chi:0\n");  //  chimeric_NOTUSED
  fprintf(fout,"cha:"F_S32"\n",mesg->chaff);
  fprintf(fout,"clr:"F_COORD","F_COORD"\n", mesg->clear_rng.bgn,mesg->clear_rng.end);
  fprintf(fout,"}\n");
  return;
}

static void Write_AMP_Mesg(FILE *fout, void *vmesg)
{ AugMatePairMesg *mesg = (AugMatePairMesg *) vmesg;
  
  fprintf(fout,"{AMP\n");
  fprintf(fout,"frg:%s\n",AS_UID_toString(mesg->fragment1));
  fprintf(fout,"frg:%s\n",AS_UID_toString(mesg->fragment2));
  fprintf(fout,"mst:%c\n",mesg->mate_status);
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
  fprintf(fout,"chi:0\n");  //  chimeric_NOTUSED
  fprintf(fout,"cha:"F_S32"\n",mesg->chaff);
  fprintf(fout,"clr:"F_COORD","F_COORD"\n", mesg->clear_rng.bgn,mesg->clear_rng.end);
  fprintf(fout,"mst:%c\n",mesg->mate_status);
  fprintf(fout,"}\n");
  return;
}

static void Write_IAM_Mesg(FILE *fout, void *vmesg)
{ IntAugMatePairMesg *mesg = (IntAugMatePairMesg *) vmesg;
  
  fprintf(fout,"{IAM\n");
  fprintf(fout,"frg:"F_IID"\n",mesg->fragment1);
  fprintf(fout,"frg:"F_IID"\n",mesg->fragment2);
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
  fprintf(fout,"lid:%s\n",AS_UID_toString(up->eident));
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
  fprintf(fout,"mid:%s\n",AS_UID_toString(mlp->eident));
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
  fprintf(fout,"acc:(%s,"F_IID")\n", AS_UID_toString(mesg->eaccession),mesg->iaccession);
#ifdef AS_ENABLE_SOURCE
  PutText(fout,"src:",mesg->source,FALSE);
#endif
  fprintf(fout,"cov:%.3f\n",mesg->coverage_stat);
  fprintf(fout,"sta:%c\n",mesg->status);
  fprintf(fout,"abp:0\n");  //  Unused!
  fprintf(fout,"bbp:0\n");  //  Unused!
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
  fprintf(fout,"ut1:%s\n",AS_UID_toString(mesg->eunitig1));
  fprintf(fout,"ut2:%s\n",AS_UID_toString(mesg->eunitig2));
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
    fprintf(fout,"%s,%s,%c\n",
            AS_UID_toString(mesg->jump_list[i].in1),
            AS_UID_toString(mesg->jump_list[i].in2),
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}


static void Write_CCO_Mesg(FILE *fout, void *vmesg)
{ SnapConConMesg *mesg = (SnapConConMesg *) vmesg;
  int		i;

  fprintf(fout,"{CCO\n");
  fprintf(fout,"acc:(%s,"F_IID")\n",AS_UID_toString(mesg->eaccession),mesg->iaccession);
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
  fprintf(fout,"co1:%s\n",AS_UID_toString(mesg->econtig1));
  fprintf(fout,"co2:%s\n",AS_UID_toString(mesg->econtig2));
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
    fprintf(fout, "%s,%s,%c\n",
            AS_UID_toString(mesg->jump_list[i].in1),
            AS_UID_toString(mesg->jump_list[i].in2),
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}

static void Write_SLK_Mesg(FILE *fout, void *vmesg)
{ SnapScaffoldLinkMesg *mesg = (SnapScaffoldLinkMesg *) vmesg;
  int		i, npairs;

  fprintf(fout,"{SLK\n");
  fprintf(fout,"sc1:%s\n",AS_UID_toString(mesg->escaffold1));
  fprintf(fout,"sc2:%s\n",AS_UID_toString(mesg->escaffold2));
  fprintf(fout,"ori:%c\n",mesg->orientation);
  fprintf(fout,"gui:"F_S32"\n",mesg->includes_guide);
  fprintf(fout,"mea:%.3f\n",mesg->mean_distance);
  fprintf(fout,"std:%.3f\n",mesg->std_deviation);
  fprintf(fout,"num:"F_S32"\n",mesg->num_contributing);
  npairs = mesg->num_contributing;
  assert(npairs > 0);
  fprintf(fout,"jls:\n");
  for (i=0; i < npairs; ++i)
    fprintf(fout, "%s,%s,%c\n",
            AS_UID_toString(mesg->jump_list[i].in1),
            AS_UID_toString(mesg->jump_list[i].in2),
            (char)(mesg->jump_list[i].type));
  fprintf(fout,"}\n");
  return;
}


static void Write_CTP_Mesg(FILE *fout, SnapContigPairs *mesg)
{
  fprintf(fout,"{CTP\n");
  fprintf(fout,"ct1:%s\n",AS_UID_toString(mesg->econtig1));
  fprintf(fout,"ct2:%s\n",AS_UID_toString(mesg->econtig2));
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
  fprintf(fout,"acc:(%s,"F_IID")\n",AS_UID_toString(mesg->eaccession),mesg->iaccession);
  fprintf(fout,"noc:"F_S32"\n",mesg->num_contig_pairs);
  for (i=0; i < num; ++i)
    Write_CTP_Mesg(fout,&mesg->contig_pairs[i]);
  fprintf(fout,"}\n");
  return;
}


static void Write_MDI_Mesg(FILE *fout, void *vmesg)
{ SnapMateDistMesg *mesg = (SnapMateDistMesg *) vmesg;
  int		i;

  fprintf(fout,"{MDI\n");
  fprintf(fout,"ref:(%s,"F_IID")\n",AS_UID_toString(mesg->erefines),mesg->irefines);
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
  fprintf(fout,"acc:%s\n",AS_UID_toString(mesg->eaccession));
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
  {"{BAT", Read_BAT_Mesg, Write_BAT_Mesg, sizeof(BatchMesg) },
  {"{VER", Read_VER_Mesg, Write_VER_Mesg, sizeof(VersionMesg)  },
  {"{DST", Read_DST_Mesg, Write_DST_Mesg, sizeof(DistanceMesg) },
  {"RLIB", NULL, NULL, 0l },  //  RESERVED for Version 2's LIB message
  {"{FRG", Read_FRG_Mesg, Write_FRG_Mesg, sizeof(FragMesg)  },
  {"{LKG", Read_LKG_Mesg, Write_LKG_Mesg, sizeof(LinkMesg) },

  {"{OVL", Read_OVL_Mesg, Write_OVL_Mesg, sizeof(OverlapMesg) },
  {"{UOM", Read_UOM_Mesg, Write_UOM_Mesg, sizeof(UnitigOverlapMesg) },

  {"{IMD", Read_IMD_Mesg, Write_IMD_Mesg, sizeof(IntMateDistMesg) },
  {"{IAF", Read_IAF_Mesg, Write_IAF_Mesg, sizeof(IntAugFragMesg) },
  {"{IAM", Read_IAM_Mesg, Write_IAM_Mesg, sizeof(IntAugMatePairMesg) },
  {"{IUM", Read_IUM_Mesg, Write_IUM_Mesg, sizeof(IntUnitigMesg) },
  {"{IUL", Read_IUL_Mesg, Write_IUL_Mesg, sizeof(IntUnitigLinkMesg) },
  {"{ICM", Read_ICM_Mesg, Write_ICM_Mesg, sizeof(IntConConMesg) },
  {"{ICL", Read_ICL_Mesg, Write_ICL_Mesg, sizeof(IntContigLinkMesg) },
  {"{ISF", Read_ISF_Mesg, Write_ISF_Mesg, sizeof(IntScaffoldMesg) },
  {"{ISL", Read_ISL_Mesg, Write_ISL_Mesg, sizeof(InternalScaffoldLinkMesg) },

  {"{MDI", Read_MDI_Mesg, Write_MDI_Mesg, sizeof(SnapMateDistMesg) },
  {"{AFG", Read_AFG_Mesg, Write_AFG_Mesg, sizeof(AugFragMesg) },
  {"{AMP", Read_AMP_Mesg, Write_AMP_Mesg, sizeof(AugMatePairMesg) },
  {"{UTG", Read_UTG_Mesg, Write_UTG_Mesg, sizeof(SnapUnitigMesg) },
  {"{ULK", Read_ULK_Mesg, Write_ULK_Mesg, sizeof(SnapUnitigLinkMesg) },
  {"{CCO", Read_CCO_Mesg, Write_CCO_Mesg, sizeof(SnapConConMesg) },
  {"{CLK", Read_CLK_Mesg, Write_CLK_Mesg, sizeof(SnapContigLinkMesg) },
  {"{SCF", Read_SCF_Mesg, Write_SCF_Mesg, sizeof(SnapScaffoldMesg) },
  {"{SLK", Read_SLK_Mesg, Write_SLK_Mesg, sizeof(SnapScaffoldLinkMesg) },

  {"{EOF", Read_EOF_Mesg, Write_EOF_Mesg, sizeof(EndOfFileMesg) }
};


void AS_MSG_setFormatVersion1(void) {
  memcpy(AS_MSG_globals->CallTable, CallTable1, sizeof(AS_MSG_callrecord) * (NUM_OF_REC_TYPES + 1));
}

