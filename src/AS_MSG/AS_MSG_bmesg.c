
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
static char CM_ID[] = "$Id: AS_MSG_bmesg.c,v 1.11 2007-01-28 21:52:24 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include "AS_global.h"

#define BIN_CODE 0x0013fa00

#define ROUNDUP(n,u) ((((n)-1)/(u) + 1)*(u))  /* Round n up to nearest
                                                 multiple of u */

#define FREAD(ptr,siz,num,inp) \
if (fread((ptr),(siz),(num),(inp)) != (num)) ReadError(CDS_FTELL((inp)));

#define FWRITE(ptr,siz,num,inp) \
if (fwrite((ptr),(siz),(num),(inp)) != (num)) WriteError();

static int MessageNum;  /* Number of current message */
static int Mcode;       /* Enum type value of current message routine */


static char *MemBuffer = NULL;   /* Memory allocation buffer for messages */
static int   MemMax = -1, MemTop;   /* Memory ceiling and current top */

/* Make sure there is a block of size bytes left in memory buffer starting
   at an index that is a multiple of boundary.  Return the *index* into the
   array, so that the realloc does not blow structures in the process of
   being built.  All pointers in such structures are saved as integers and
   converted to pointers after all allocation has taken place.             */

void MakeSpaceB(int size){
  size_t   newsize=1;
  char *newbufr;

  if(MemMax> size)
    return;

  if(MemMax < 0)
    newsize = 2 * 2048 * 2048; // This may be excessive, but this code is BRITTLE!
  else
    newsize = 2 * size;

      newsize = ROUNDUP(newsize,8);
      newbufr = (char *) realloc(MemBuffer,newsize);
      /*
      if(MemBuffer)
        fprintf(stderr,
                "* Reallocing MemBuffer to size " F_SIZE_T " old:%p  new:%p\n",
                newsize, MemBuffer, newbufr);
      */
      if (newbufr == NULL)
        { fprintf(stderr,"ERROR: Out of memory \n");
	assert(0);
	exit (1);
        }
      MemBuffer = newbufr;
      MemMax    = newsize;
}


static long MoreSpace(int size, int boundary)
{ 
  MemTop = ROUNDUP(MemTop,boundary);

  MakeSpaceB(MemTop + size);
  { int alloc;

    alloc   = MemTop;
    MemTop += size;
    return (alloc);
  }
}

static void ReadError(int pos)
{ 
  // char *junk = NULL;
  fprintf(stderr,"ERROR: Read error at file position %d ",
	  pos);
  fprintf(stderr," at message %d of type %s\n",
                 MessageNum,MessageTypeName[Mcode]);
  fprintf(stderr," Your input file was probably truncated...check its size.\n");
  exit (1);
}

static void WriteError(void)
{ fprintf(stderr,"ERROR: Write error");
  fprintf(stderr," at message %d of type %s\n",
                 MessageNum,MessageTypeName[Mcode]);
  assert(0);
  //  exit (1);
}

   /* Get a string field item */

static long GetString(FILE *fin)
{ long text;
  int len;
 
  FREAD(&len,sizeof(len),1,fin);
  text = MoreSpace(len+1,1);
  FREAD(MemBuffer+text,sizeof(char),len,fin);
  MemBuffer[text+len] = '\0';
  return (text);
}

   /* Output text field item with 3-code field-name "tag". */

static void PutString(FILE *fout, char *text)
{ int len;
  
 if(text != NULL){ 
   len = strlen(text);
   FWRITE(&len,sizeof(len),1,fout);
   FWRITE(text,sizeof(char),len,fout);
 }else{
   len = 0;
   FWRITE(&len,sizeof(len),1,fout);
 }
}

/******************** INPUT ROUTINES ***************************/

static int Read_ADL_Struct(long last, FILE *fin)
{ 
  AuditLine mesg; // temporary
  long crnt;
  AuditLine *newMesg;

  FREAD(&mesg,sizeof(AuditLine),1,fin);
  mesg.name    = (char *) GetString(fin);
  mesg.version = (char *) GetString(fin);
  mesg.comment = (char *) GetString(fin);
  mesg.next = (AuditLine *)last;

  // Now allocate space for this guy, and copy him to the allocated space
  crnt = MoreSpace(sizeof(AuditLine),8);
  newMesg = (AuditLine *)(MemBuffer + crnt);
  *newMesg = mesg;

  return crnt;
}

static void Read_ADT_Mesg(FILE *fin, void *vmesg)
{ AuditMesg *amesg = (AuditMesg *) vmesg;
  AuditLine *cptr, *tail;
  long       last,  crnt;
  int        idx;

  /* First build up list (in reverse order) using indices and not pointers */

  last = crnt = -1;
  FREAD(&idx,sizeof(idx),1,fin);
  while (idx-- > 0)
    { 
      last =  Read_ADL_Struct(last,fin);
    }

  /* Traverse again, reversing list order and converting indices to ptrs. */

  crnt = last;
  tail = NULL;
  while (crnt >= 0)
    { cptr = (AuditLine *) (MemBuffer + crnt);
      crnt = (long) (cptr->next);
      cptr->next    = tail;
      cptr->name    = MemBuffer + ((long) (cptr->name));
      cptr->version = MemBuffer + ((long) (cptr->version));
      cptr->comment = MemBuffer + ((long) (cptr->comment));
      tail = cptr;
    }
  amesg->list = tail;
}

static void Read_Frag_Mesg(FILE *fin, void *vmesg, int frag_class)
{ ScreenedFragMesg *fmesg = (ScreenedFragMesg *) vmesg;

  if (fmesg->action == AS_ADD)
    { fmesg->source   = (char *) GetString(fin);
      if(frag_class < 2)
        {
          fmesg->sequence = (char *) GetString(fin);
          fmesg->quality  = (char *) GetString(fin);
        }
      fmesg->source   = MemBuffer + ((long) (fmesg->source));
      if(frag_class < 2)
        {
          fmesg->sequence = MemBuffer + ((long) (fmesg->sequence));
          fmesg->quality  = MemBuffer + ((long) (fmesg->quality));
        }
    }
}

static void Read_FRG_Mesg(FILE *fin, void *vmesg)
{ Read_Frag_Mesg(fin,vmesg,0); }

static void Read_SFG_Mesg(FILE *fin, void *vmesg)
{ Read_Frag_Mesg(fin,vmesg,1); }

static void Read_OFG_Mesg(FILE *fin, void *vmesg)
{ Read_Frag_Mesg(fin,vmesg,2); }

static void Read_OVL_Mesg(FILE *fin, void *vmesg)
{ OverlapMesg *omesg = (OverlapMesg *) vmesg;

  omesg->delta = (signed char *) GetString(fin);
  omesg->delta = (signed char *) (MemBuffer + ((long) (omesg->delta)));
}


// THIS ROUTINE IS UNSAFE IF THE MEMBUFFER IS REALLOCATED
// 


#ifdef AS_ENABLE_SOURCE
static void Read_UOM_Mesg(FILE *fin, void *vmesg)
{
  UnitigOverlapMesg *mesg = (UnitigOverlapMesg *) vmesg;
  long		sindx;

  sindx = GetString(fin);
  mesg->source = MemBuffer + sindx;
}
static void Read_FOM_Mesg(FILE *fin, void *vmesg)
{
  FragOverlapMesg *mesg = (FragOverlapMesg *) vmesg;
  long		sindx;

  sindx = GetString(fin);
  mesg->source = MemBuffer + sindx;
}
#endif

static void Read_IUM_Mesg(FILE *fin, void *vmesg)
{ 
  IntUnitigMesg *mesg = (IntUnitigMesg *) vmesg;
  long		cindx, qindx, mpindx, indx;
  // long	   dindx;
# ifdef AS_ENABLE_SOURCE
  long		sindx;
# endif
  IntMultiPos	*mlp;
  // int16		*delta;
  int		i;

  cindx = GetString(fin);
  qindx = GetString(fin);
# ifdef AS_ENABLE_SOURCE
  sindx = GetString(fin);
# endif
  if (mesg->num_frags > 0) {
    mpindx = MoreSpace(sizeof(IntMultiPos)*mesg->num_frags, 8);
    for(i=0; i < mesg->num_frags; i++) {
      mlp = ((IntMultiPos *) (MemBuffer + mpindx)) + i;
      FREAD(mlp,sizeof(IntMultiPos),1,fin);
      if (mlp->delta_length > 0) {
	indx = MoreSpace(sizeof(int32)*mlp->delta_length,8);
	/* The next line guards against reallocs */
	mlp = ((IntMultiPos *) (MemBuffer + mpindx)) + i;
	mlp->delta = (int32 *) indx;
	FREAD((int32 *)(MemBuffer + (long)mlp->delta),sizeof(int32),
	      mlp->delta_length,fin);
      }
      else
	mlp->delta = NULL;
      mesg->f_list = (IntMultiPos *) (MemBuffer + mpindx);
    }
  }
  else
    mesg->f_list = NULL;
  mesg->consensus = MemBuffer + cindx;
  mesg->quality  = MemBuffer + qindx;
  for (i=0; i < mesg->num_frags; ++i) {
    if (mesg->f_list[i].delta_length > 0)
      mesg->f_list[i].delta = (int32 *) (MemBuffer+(long)mesg->f_list[i].delta);
  }
# ifdef AS_ENABLE_SOURCE
  mesg->source = MemBuffer + sindx;
# endif
}

static void Read_IUL_Mesg(FILE *fin, void *vmesg)
{
  IntUnitigLinkMesg *mesg = (IntUnitigLinkMesg *) vmesg;
  long		indx;
  int		total;

  total = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --total;
  if (total > 0) {
    indx = MoreSpace(sizeof(IntMate_Pairs)*total,8);
    mesg->jump_list = (IntMate_Pairs *) (MemBuffer + indx);
    FREAD(mesg->jump_list,sizeof(IntMate_Pairs),total,fin);
  }
  else
    mesg->jump_list = NULL;
}

static void Read_ICL_Mesg(FILE *fin, void *vmesg)
{
  IntContigLinkMesg *mesg = (IntContigLinkMesg *) vmesg;
  long		indx;
  int		total;

  total = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --total;
  if (total > 0) {
    indx = MoreSpace(sizeof(IntMate_Pairs)*total,8);
    mesg->jump_list = (IntMate_Pairs *) (MemBuffer + indx);
    FREAD(mesg->jump_list,sizeof(IntMate_Pairs),total,fin);
  }
  else
    mesg->jump_list = NULL;
}

static void Read_ISL_Mesg(FILE *fin, void *vmesg)
{
  InternalScaffoldLinkMesg *mesg = (InternalScaffoldLinkMesg *) vmesg;
  long		indx;
  int		total;

  total = mesg->num_contributing;
  assert(total > 0) ;
  indx = MoreSpace(sizeof(IntMate_Pairs)*total,8);
  mesg->jump_list = (IntMate_Pairs *) (MemBuffer + indx);
  FREAD(mesg->jump_list,sizeof(IntMate_Pairs),total,fin);
}

static void Read_AFG_Mesg(FILE *fin, void *vmesg)
{
  AugFragMesg *mesg = (AugFragMesg *) vmesg;
}

static void Read_ISF_Mesg(FILE *fin, void *vmesg)
{
  IntScaffoldMesg *mesg = (IntScaffoldMesg *) vmesg;
  long		indx;
  int num = MAX(1, mesg->num_contig_pairs);
  if ( num > 0) {
    indx = MoreSpace(sizeof(IntContigPairs)*(num),8);
    mesg->contig_pairs = (IntContigPairs *) (MemBuffer + indx);
    FREAD(mesg->contig_pairs,sizeof(IntContigPairs),num,fin);
  }
  else
    mesg->contig_pairs = NULL;
}

static void Read_IMD_Mesg(FILE *fin, void *vmesg)
{
  IntMateDistMesg *mesg = (IntMateDistMesg *) vmesg;
  long		indx;

  if (mesg->num_buckets > 0) {
    indx = MoreSpace(sizeof(int32)*mesg->num_buckets,8);
    mesg->histogram = (int32 *) (MemBuffer + indx);
    FREAD(mesg->histogram,sizeof(int32),mesg->num_buckets,fin);
  }
  else
    mesg->histogram = NULL;
}




/* Snap shot input routines */
/****************************/

static void Read_UTG_Mesg(FILE *fin, void *vmesg)
{ 
  SnapUnitigMesg *mesg = (SnapUnitigMesg *) vmesg;
  long		cindx, qindx, mpindx, indx;
  // long	   dindx;
  #ifdef AS_ENABLE_SOURCE
  long		sindx;
  #endif
  SnapMultiPos	*mlp;
  // int16		*delta;
  int		i;

  cindx = GetString(fin);
  qindx = GetString(fin);
  #ifdef AS_ENABLE_SOURCE
  sindx = GetString(fin);
  #endif
  if (mesg->num_frags > 0) {
    mpindx = MoreSpace(sizeof(SnapMultiPos)*mesg->num_frags, 8);
    for(i=0; i < mesg->num_frags; i++) {
      mlp = ((SnapMultiPos *) (MemBuffer + mpindx)) + i;
      FREAD(mlp,sizeof(SnapMultiPos),1,fin);
      #ifdef AS_ENABLE_SOURCE
      indx = GetString(fin);
      mlp = ((SnapMultiPos *) (MemBuffer + mpindx)) + i; // in case realloc
      mlp->source = (char *) indx;
      #endif
      if (mlp->delta_length > 0) {
	indx = MoreSpace(sizeof(int32)*mlp->delta_length,8);
	/* The next line guards against reallocs */
	mlp = ((SnapMultiPos *) (MemBuffer + mpindx)) + i;
	mlp->delta = (int32 *) indx;
	FREAD((int32 *)(MemBuffer + (long)mlp->delta),sizeof(int32),
	      mlp->delta_length,fin);
      }
      else
	mlp->delta = NULL;
      mesg->f_list = (SnapMultiPos *) (MemBuffer + mpindx);
    }
  }
  else
    mesg->f_list = NULL;
  mesg->consensus = MemBuffer + cindx;
  mesg->quality  = MemBuffer + qindx;
  for (i=0; i < mesg->num_frags; ++i) {
    #ifdef AS_ENABLE_SOURCE
    mesg->f_list[i].source = MemBuffer + (long) mesg->f_list[i].source;
    #endif
    if (mesg->f_list[i].delta_length > 0)
      mesg->f_list[i].delta = (int32 *) (MemBuffer+(long)mesg->f_list[i].delta);
  }
  #ifdef AS_ENABLE_SOURCE
  mesg->source = MemBuffer + sindx;
  #endif
}


static void Read_ULK_Mesg(FILE *fin, void *vmesg)
{
  SnapUnitigLinkMesg *mesg = (SnapUnitigLinkMesg *) vmesg;
  long		indx;
  int		total;

  total = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --total;
  if (total > 0) {
    indx = MoreSpace(sizeof(SnapMate_Pairs)*total,8);
    mesg->jump_list = (SnapMate_Pairs *) (MemBuffer + indx);
    FREAD(mesg->jump_list,sizeof(SnapMate_Pairs),total,fin);
  }
  else
    mesg->jump_list = NULL;
}



static void Read_CCO_Mesg(FILE *fin, void *vmesg)
{
  SnapConConMesg *mesg = (SnapConConMesg *) vmesg;
  SnapMultiPos	*mlp;
  UnitigPos *up;
  IntMultiVar *imv;

  int		i;
  int32		*delta;
  long		cindx, qindx, pindx, uindx, indx, vindx;

  cindx = GetString(fin);
  qindx = GetString(fin);
  pindx = MoreSpace(sizeof(SnapMultiPos)*mesg->num_pieces,8);
  vindx = MoreSpace(sizeof(IntMultiVar)*mesg->num_vars,8);
  uindx = MoreSpace(sizeof(UnitigPos)*mesg->num_unitigs,8);
  if (mesg->num_pieces > 0) {
    for (i=0; i < mesg->num_pieces; ++i) {
      mlp = ((SnapMultiPos *) (MemBuffer + pindx)) + i;
      FREAD(mlp,sizeof(SnapMultiPos),1,fin);
      #ifdef AS_ENABLE_SOURCE
      indx = GetString(fin);
      mlp = ((SnapMultiPos *) (MemBuffer + pindx)) + i; // in case of realloc
      mlp->source = (char *) indx;
      #endif
      if (mlp->delta_length > 0) {
	indx = MoreSpace(sizeof(int32)*mlp->delta_length,8);
	mlp = ((SnapMultiPos *) (MemBuffer + pindx)) + i; // in case of realloc
	mlp->delta = (int32 *) indx;
	delta = (int32 *) (MemBuffer + (long) mlp->delta);
	FREAD(delta,sizeof(int32),mlp->delta_length,fin);
      }
      else
	mlp->delta = NULL;
    }
    mesg->pieces = (SnapMultiPos *) (MemBuffer + pindx);
  }
  else
    mesg->pieces = NULL;

  if (mesg->num_vars > 0) {
    for (i=0; i < mesg->num_vars; ++i) {
      imv = ((IntMultiVar *) (MemBuffer + vindx)) + i;
      FREAD(imv,sizeof(IntMultiVar),1,fin);
    }
    mesg->vars   = (IntMultiVar *) (MemBuffer + vindx);
  }
  else
    mesg->vars   = NULL;


  if (mesg->num_unitigs > 0) {
    for (i=0; i < mesg->num_unitigs; ++i) {
      up = ((UnitigPos *) (MemBuffer + uindx)) + i;
      FREAD(up,sizeof(UnitigPos),1,fin);
      if (up->delta_length > 0) {
	indx = MoreSpace(sizeof(int32)*up->delta_length,8);
	up = ((UnitigPos *) (MemBuffer + uindx)) + i; // in case of realloc
	up->delta = (int32 *) indx;
	delta = (int32 *) (MemBuffer + (long) up->delta);
	FREAD(delta,sizeof(int32),up->delta_length,fin);
      }
      else
	up->delta = NULL;
    }
    mesg->unitigs = (UnitigPos *) (MemBuffer + uindx);
  }
  else
    mesg->unitigs = NULL;



  mesg->consensus = MemBuffer + cindx;
  mesg->quality = MemBuffer + qindx;
  
  if(mesg->num_pieces > 0){
    mesg->pieces = (SnapMultiPos *) (MemBuffer + pindx);
    for (i=0; i < mesg->num_pieces; ++i) {
#ifdef AS_ENABLE_SOURCE
      mesg->pieces[i].source = MemBuffer + (long) mesg->pieces[i].source;
#endif
      if (mesg->pieces[i].delta_length > 0)
	mesg->pieces[i].delta = 
	  (int32 *) (MemBuffer + (long) mesg->pieces[i].delta);
    }
  }
  if(mesg->num_unitigs > 0){
    mesg->unitigs = (UnitigPos *) (MemBuffer + uindx);
    for (i=0; i < mesg->num_unitigs; ++i) {
      if (mesg->unitigs[i].delta_length > 0)
	mesg->unitigs[i].delta = 
	  (int32 *) (MemBuffer + (long) mesg->unitigs[i].delta);
    }
  }
}



static void Read_CLK_Mesg(FILE *fin, void *vmesg)
{
  SnapContigLinkMesg *mesg = (SnapContigLinkMesg *) vmesg;
  long		indx;
  int		total;

  total = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --total;
  if (total > 0) {
    indx = MoreSpace(sizeof(SnapMate_Pairs)*total,8);
    mesg->jump_list = (SnapMate_Pairs *) (MemBuffer + indx);
    FREAD(mesg->jump_list,sizeof(SnapMate_Pairs),total,fin);
  }
  else
    mesg->jump_list = NULL;
}


static void Read_SLK_Mesg(FILE *fin, void *vmesg)
{
  SnapScaffoldLinkMesg *mesg = (SnapScaffoldLinkMesg *) vmesg;
  long		indx;
  int		total;

  total = mesg->num_contributing;
  assert(total > 0) ;
  indx = MoreSpace(sizeof(SnapMate_Pairs)*total,8);
  mesg->jump_list = (SnapMate_Pairs *) (MemBuffer + indx);
  FREAD(mesg->jump_list,sizeof(SnapMate_Pairs),total,fin);
}



static void Read_SCF_Mesg(FILE *fin, void *vmesg)
{
  SnapScaffoldMesg *mesg = (SnapScaffoldMesg *) vmesg;
  long		indx;
  int num = MAX(1, mesg->num_contig_pairs);

  if (num > 0) {
    indx = MoreSpace(sizeof(SnapContigPairs)*(mesg->num_contig_pairs),8);
    mesg->contig_pairs = (SnapContigPairs *) (MemBuffer + indx);
    FREAD(mesg->contig_pairs,sizeof(SnapContigPairs),num,
	  fin);
  }
  else
    mesg->contig_pairs = NULL;
}


static void Read_MDI_Mesg(FILE *fin, void *vmesg)
{
  SnapMateDistMesg *mesg = (SnapMateDistMesg *) vmesg;
  long		indx;

  if (mesg->num_buckets > 0) {
    indx = MoreSpace(sizeof(int32)*mesg->num_buckets,8);
    mesg->histogram = (int32 *) (MemBuffer + indx);
    FREAD(mesg->histogram,sizeof(int32),mesg->num_buckets,fin);
  }
  else
    mesg->histogram = NULL;
}

static void Read_BAT_Mesg(FILE *fin, void *vmesg)
{
  BatchMesg *mesg = (BatchMesg *) vmesg;
  mesg->name = (char *)GetString(fin);
  mesg->comment = (char *)GetString(fin);
  mesg->name = MemBuffer + (size_t)mesg->name;
  mesg->comment = MemBuffer + (size_t)mesg->comment;
}

static void Read_IBA_Mesg(FILE *fin, void *vmesg)
{
  InternalBatchMesg *mesg = (InternalBatchMesg *) vmesg;
  mesg->name = (char *)GetString(fin);
  mesg->comment = (char *)GetString(fin);
  mesg->name = MemBuffer + (size_t)mesg->name;
  mesg->comment = MemBuffer + (size_t)mesg->comment;

}
static void Read_BAC_Mesg(FILE *fin, void *vmesg)
{
  BacMesg *mesg = (BacMesg *) vmesg;
  int indx;

  if(mesg->num_bactigs > 0){
    indx = MoreSpace(sizeof(BactigMesg) * mesg->num_bactigs, 8);
    mesg->bactig_list = (BactigMesg *)(MemBuffer + indx);
    FREAD(mesg->bactig_list,sizeof(BactigMesg),mesg->num_bactigs,  fin);
  }
  indx = GetString(fin);
  mesg->source = (char *)(MemBuffer + indx);
}

static void Read_IBC_Mesg(FILE *fin, void *vmesg)
{
  InternalBacMesg *mesg = (InternalBacMesg *) vmesg;
  int indx;

  if(mesg->num_bactigs > 0){
    indx = MoreSpace(sizeof(InternalBactigMesg) * mesg->num_bactigs, 8);
    mesg->bactig_list = (InternalBactigMesg *)(MemBuffer + indx);
    FREAD(mesg->bactig_list,sizeof(BactigMesg),mesg->num_bactigs,  fin);
  }
  indx = GetString(fin);
  mesg->source = (char *)(MemBuffer + indx);
}

static void Read_EOF_Mesg(FILE *fin, void *vmesg)
{
  EndOfFileMesg *mesg = (EndOfFileMesg *) vmesg;

  mesg->comment = (char *) GetString(fin);
  mesg->comment = MemBuffer + (long) (mesg->comment);
}


/******************** OUTPUT ROUTINES ***************************/

static void Write_ADL_Struct(FILE *fout, AuditLine *mesg)
{ FWRITE(mesg,sizeof(AuditLine),1,fout);
  PutString(fout,mesg->name);
  PutString(fout,mesg->version);
  PutString(fout,mesg->comment);
}

static void Write_ADT_Mesg(FILE *fout, void *vmesg)
{ AuditMesg *mesg = (AuditMesg *) vmesg;
  AuditLine *a;
  int        len;

  len = 0;
  for (a = mesg->list; a != NULL; a = a->next)
    len += 1;
  FWRITE(&len,sizeof(len),1,fout);
  for (a = mesg->list; a != NULL; a = a->next)
    Write_ADL_Struct(fout,a);
}

static void Write_Frag_Mesg(FILE *fout, void *vmesg, int frag_class)
{ ScreenedFragMesg *mesg = (ScreenedFragMesg *) vmesg;

  if (mesg->action == AS_ADD)
    { PutString(fout,mesg->source);
      if (frag_class < 2)
        {
          PutString(fout,mesg->sequence);
          PutString(fout,mesg->quality);
        }
    }
}

static void Write_FRG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,0); }

static void Write_SFG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,1); }

static void Write_OFG_Mesg(FILE *fout, void *vmesg)
{ Write_Frag_Mesg(fout,vmesg,2); }

static void Write_OVL_Mesg(FILE *fout, void *vmesg)
{ OverlapMesg *omesg = (OverlapMesg *) vmesg;

  PutString(fout,(char *) (omesg->delta));
}



#ifdef AS_ENABLE_SOURCE
static void Write_UOM_Mesg(FILE *fout, void *vmesg)
{
  UnitigOverlapMesg *mesg = (UnitigOverlapMesg *) vmesg;

  PutString(fout,mesg->source);
}
static void Write_FOM_Mesg(FILE *fout, void *vmesg)
{
  FragOverlapMesg *mesg = (FragOverlapMesg *) vmesg;

  PutString(fout,mesg->source);
}
#endif

static void Write_IUM_Mesg(FILE *fout, void *vmesg)
{ 
  IntUnitigMesg *mesg = (IntUnitigMesg *) vmesg;
  int i;

  PutString(fout,mesg->consensus);
  PutString(fout,mesg->quality);
# ifdef AS_ENABLE_SOURCE
  PutString(fout,mesg->source);
# endif
  for (i=0; i<mesg->num_frags; i++) {
    FWRITE(&mesg->f_list[i],sizeof(IntMultiPos),1,fout);
    if (mesg->f_list[i].delta_length > 0)
      FWRITE(mesg->f_list[i].delta,sizeof(int32),
	     mesg->f_list[i].delta_length,fout);
  }
}

static void Write_IUL_Mesg(FILE *fout, void *vmesg)
{
  IntUnitigLinkMesg *mesg = (IntUnitigLinkMesg *) vmesg;
  int		total;

  total = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --total;
  if (total > 0)
    FWRITE(mesg->jump_list,sizeof(IntMate_Pairs),total,fout);
}

static void Write_ICL_Mesg(FILE *fout, void *vmesg)
{
  IntContigLinkMesg *mesg = (IntContigLinkMesg *) vmesg;
  int		total;

  total = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --total;
  if (total > 0)
    FWRITE(mesg->jump_list,sizeof(IntMate_Pairs),total,fout);
}

static void Write_ISL_Mesg(FILE *fout, void *vmesg)
{
  InternalScaffoldLinkMesg *mesg = (InternalScaffoldLinkMesg *) vmesg;
  int		total;

  total = mesg->num_contributing;
  assert(total > 0);
  FWRITE(mesg->jump_list,sizeof(IntMate_Pairs),total,fout);
}

static void Write_AFG_Mesg(FILE *fout, void *vmesg)
{
  AugFragMesg *mesg = (AugFragMesg *) vmesg;
}

static void Write_ISF_Mesg(FILE *fout, void *vmesg)
{
  IntScaffoldMesg *mesg = (IntScaffoldMesg *) vmesg;
  int num = MAX(1, mesg->num_contig_pairs);

  if (num > 0)
    FWRITE(mesg->contig_pairs,sizeof(IntContigPairs),num,
	   fout);
}

static void Write_IMD_Mesg(FILE *fout, void *vmesg)
{
  IntMateDistMesg *mesg = (IntMateDistMesg *) vmesg;

  if (mesg->num_buckets > 0)
    FWRITE(mesg->histogram,sizeof(int32),mesg->num_buckets,fout);
}



/* Genome Snapshot output routines */
/***********************************/


static void Write_UTG_Mesg(FILE *fout, void *vmesg)
{ 
  SnapUnitigMesg *mesg = (SnapUnitigMesg *) vmesg;
  int i;

  PutString(fout,mesg->consensus);
  PutString(fout,mesg->quality);
  #ifdef AS_ENABLE_SOURCE
  PutString(fout,mesg->source);
  #endif
  for (i=0; i<mesg->num_frags; i++) {
    FWRITE(&mesg->f_list[i],sizeof(SnapMultiPos),1,fout);
    #ifdef AS_ENABLE_SOURCE
    PutString(fout,mesg->f_list[i].source);
    #endif
    if (mesg->f_list[i].delta_length > 0)
      FWRITE(mesg->f_list[i].delta,sizeof(int32),
	     mesg->f_list[i].delta_length,fout);
  }
}


static void Write_ULK_Mesg(FILE *fout, void *vmesg)
{
  SnapUnitigLinkMesg *mesg = (SnapUnitigLinkMesg *) vmesg;
  int		total;

  total = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --total;
  if (total > 0)
    FWRITE(mesg->jump_list,sizeof(SnapMate_Pairs),total,fout);
}


static void Read_ICM_Mesg(FILE *fin, void *vmesg)
{
  IntConConMesg *mesg = (IntConConMesg *) vmesg;
  IntMultiPos	*mlp;
  IntUnitigPos	*iup;
  IntMultiVar   *imv;
  int		 i;
  int32		*delta;
  long		 cindx, qindx, pindx, uindx, indx, vindx;

  cindx = GetString(fin);
  qindx = GetString(fin);
  pindx = MoreSpace(sizeof(IntMultiPos)*mesg->num_pieces,8);
  vindx = MoreSpace(sizeof(IntMultiVar)*mesg->num_vars,8);
  uindx = MoreSpace(sizeof(IntUnitigPos)*mesg->num_unitigs,8);

  if (mesg->num_pieces > 0) {
    for (i=0; i < mesg->num_pieces; ++i) {
      mlp = ((IntMultiPos *) (MemBuffer + pindx)) + i;
      FREAD(mlp,sizeof(IntMultiPos),1,fin);
      if (mlp->delta_length > 0) {
	indx = MoreSpace(sizeof(int32)*mlp->delta_length,8);
	mlp = ((IntMultiPos *) (MemBuffer + pindx)) + i; // in case of realloc
	mlp->delta = (int32 *) indx;
	delta = (int32 *) (MemBuffer + (long) mlp->delta);
	FREAD(delta,sizeof(int32),mlp->delta_length,fin);
      }
      else
	mlp->delta = NULL;
    }
    mesg->pieces = (IntMultiPos *) (MemBuffer + pindx);
  }
  else
    mesg->pieces = NULL;

  if (mesg->num_vars > 0) {
    for (i=0; i < mesg->num_vars; ++i) {
      imv = ((IntMultiVar *) (MemBuffer + vindx)) + i;
      FREAD(imv,sizeof(IntMultiVar),1,fin);
    }
    mesg->v_list = (IntMultiVar *) (MemBuffer + vindx);
  }
  else
    mesg->v_list = NULL;

  if (mesg->num_unitigs > 0) {
    for (i=0; i < mesg->num_unitigs; ++i) {
      iup = ((IntUnitigPos *) (MemBuffer + uindx)) + i;
      FREAD(iup,sizeof(IntUnitigPos),1,fin);
      if (iup->delta_length > 0) {
	indx = MoreSpace(sizeof(int32)*iup->delta_length,8);
	iup = ((IntUnitigPos *) (MemBuffer + uindx)) + i; // in case of realloc
	iup->delta = (int32 *) indx;
	delta = (int32 *) (MemBuffer + (long) iup->delta);
	FREAD(delta,sizeof(int32),iup->delta_length,fin);
      }
      else
	iup->delta = NULL;
    }
    mesg->unitigs = (IntUnitigPos *) (MemBuffer + uindx);
  }
  else
    mesg->unitigs = NULL;



  mesg->consensus = MemBuffer + cindx;
  mesg->quality = MemBuffer + qindx;


  if(mesg->num_pieces > 0){
    mesg->pieces = (IntMultiPos *)(MemBuffer + pindx);
  for (i=0; i < mesg->num_pieces; ++i) {
    if (mesg->pieces[i].delta_length > 0){
      mesg->pieces[i].delta = 
		      (int32 *) (MemBuffer + (long) mesg->pieces[i].delta);
      assert(mesg->pieces[i].delta[0] > -9999);

    }
  }
  }
  if(mesg->num_unitigs > 0){
    mesg->unitigs = (IntUnitigPos *)(MemBuffer + uindx);
  for (i=0; i < mesg->num_unitigs; ++i) {
    if (mesg->unitigs[i].delta_length > 0){

      mesg->unitigs[i].delta = 
		      (int32 *) (MemBuffer + (long) mesg->unitigs[i].delta);
      assert(mesg->unitigs[i].delta[0] > -9999);
    }
  }
  }
}
static void Write_ICM_Mesg(FILE *fout, void *vmesg)
{
  IntConConMesg *mesg = (IntConConMesg *) vmesg;
  int		i;

  PutString(fout,mesg->consensus);
  PutString(fout,mesg->quality);
  for (i=0; i < mesg->num_pieces; ++i) {
    FWRITE(&mesg->pieces[i],sizeof(IntMultiPos),1,fout);
    if (mesg->pieces[i].delta_length > 0)
      FWRITE(mesg->pieces[i].delta,sizeof(int32),
	     mesg->pieces[i].delta_length,fout);
  }
  for (i=0; i < mesg->num_vars; ++i) {
    FWRITE(&mesg->v_list[i],sizeof(IntMultiVar),1,fout);
  }
  for (i=0; i < mesg->num_unitigs; ++i) {
    FWRITE(&mesg->unitigs[i],sizeof(IntUnitigPos),1,fout);
    if (mesg->unitigs[i].delta_length > 0)
      FWRITE(mesg->unitigs[i].delta,sizeof(int32),
	     mesg->unitigs[i].delta_length,fout);
  }
}

static void Write_CCO_Mesg(FILE *fout, void *vmesg)
{
  SnapConConMesg *mesg = (SnapConConMesg *) vmesg;
  int		i;

  PutString(fout,mesg->consensus);
  PutString(fout,mesg->quality);
  for (i=0; i < mesg->num_pieces; ++i) {
    FWRITE(&mesg->pieces[i],sizeof(SnapMultiPos),1,fout);
    #ifdef AS_ENABLE_SOURCE
    PutString(fout,mesg->pieces[i].source);
    #endif
    if (mesg->pieces[i].delta_length > 0)
      FWRITE(mesg->pieces[i].delta,sizeof(int32),
	     mesg->pieces[i].delta_length,fout);
  }
  for (i=0; i < mesg->num_vars; ++i) {
    FWRITE(&mesg->vars[i],sizeof(IntMultiVar),1,fout);
  }
  for (i=0; i < mesg->num_unitigs; ++i) {
    FWRITE(&mesg->unitigs[i],sizeof(UnitigPos),1,fout);
    if (mesg->unitigs[i].delta_length > 0)
      FWRITE(mesg->unitigs[i].delta,sizeof(int32),
	     mesg->unitigs[i].delta_length,fout);
  }

}


static void Write_CLK_Mesg(FILE *fout, void *vmesg)
{
  SnapContigLinkMesg *mesg = (SnapContigLinkMesg *) vmesg;
  int		total;

  total = mesg->num_contributing;
  if (mesg->overlap_type != AS_NO_OVERLAP)
    --total;
  if (total > 0)
    FWRITE(mesg->jump_list,sizeof(SnapMate_Pairs),total,fout);
}


static void Write_SLK_Mesg(FILE *fout, void *vmesg)
{
  SnapScaffoldLinkMesg *mesg = (SnapScaffoldLinkMesg *) vmesg;
  int		total;

  total = mesg->num_contributing;
  assert(total > 0);
  FWRITE(mesg->jump_list,sizeof(SnapMate_Pairs),total,fout);
}


static void Write_SCF_Mesg(FILE *fout, void *vmesg)
{
  SnapScaffoldMesg *mesg = (SnapScaffoldMesg *) vmesg;
  int num = MAX(1, mesg->num_contig_pairs);

  if (num > 1)
    FWRITE(mesg->contig_pairs,sizeof(SnapContigPairs),num,
	   fout);
}


static void Write_MDI_Mesg(FILE *fout, void *vmesg)
{
  SnapMateDistMesg *mesg = (SnapMateDistMesg *) vmesg;

  if (mesg->num_buckets > 0)
    FWRITE(mesg->histogram,sizeof(int32),mesg->num_buckets,fout);
}

static void Write_BAT_Mesg(FILE *fout, void *vmesg)
{
  BatchMesg *mesg = (BatchMesg *) vmesg;
  PutString(fout, mesg->name);
  PutString(fout, mesg->comment);

}
static void Write_IBA_Mesg(FILE *fout, void *vmesg)
{
  InternalBatchMesg *mesg = (InternalBatchMesg *) vmesg;
  PutString(fout, mesg->name);
  PutString(fout, mesg->comment);
}
static void Write_BAC_Mesg(FILE *fout, void *vmesg)
{
  BacMesg *mesg = (BacMesg *) vmesg;

  if(mesg->num_bactigs > 0)
    FWRITE(mesg->bactig_list,sizeof(InternalBactigMesg),mesg->num_bactigs,  fout);
  PutString(fout,mesg->source);

}

static void Write_IBC_Mesg(FILE *fout, void *vmesg)
{
  InternalBacMesg *mesg = (InternalBacMesg *) vmesg;
  if(mesg->num_bactigs > 0)
    FWRITE(mesg->bactig_list,sizeof(InternalBactigMesg),mesg->num_bactigs,  fout);
  PutString(fout,mesg->source);

}

static void Write_EOF_Mesg(FILE *fout, void *vmesg)
{
  EndOfFileMesg *mesg = (EndOfFileMesg *) vmesg;
  PutString(fout, mesg->comment);
}


/******************** EXTERNAL ENTRY POINTS ***************************/

/*  Routines to duplicate the second-level parts of a message.  */


typedef struct {
  void (*reader)(FILE *, void *);
  void (*writer)(FILE *, void *);
  int  size;
} callrecord;

static callrecord CallTable[] = {
  { NULL, NULL, 0},
  { Read_ADT_Mesg, Write_ADT_Mesg, sizeof(AuditMesg) },
  { Read_FRG_Mesg, Write_FRG_Mesg, sizeof(FragMesg)  },
  { Read_FRG_Mesg, Write_FRG_Mesg, sizeof(InternalFragMesg) },
  { Read_SFG_Mesg, Write_SFG_Mesg, sizeof(ScreenedFragMesg) },
  { Read_OFG_Mesg, Write_OFG_Mesg, sizeof(OFGMesg) },
  { NULL,          NULL,           sizeof(LinkMesg) },
  { NULL,          NULL,           sizeof(InternalLinkMesg) },
  { NULL,          NULL,           sizeof(DistanceMesg) },
  { NULL,          NULL,           0l },
  { NULL,          NULL,           0l },
  { NULL,          NULL,           sizeof(InternalDistMesg) },
  { NULL,          NULL,           sizeof(InternalDistMesg) },
  { Read_OVL_Mesg, Write_OVL_Mesg, sizeof(OverlapMesg) },
  { NULL,          NULL,           sizeof(BranchMesg) },
#ifdef AS_ENABLE_SOURCE
  { Read_UOM_Mesg, Write_UOM_Mesg, sizeof(UnitigOverlapMesg) },
#else
  { NULL,          NULL,           sizeof(UnitigOverlapMesg) },
#endif
  { Read_IUM_Mesg, Write_IUM_Mesg, sizeof(IntUnitigMesg) },
  { Read_IUL_Mesg, Write_IUL_Mesg, sizeof(IntUnitigLinkMesg) },
  { Read_ICL_Mesg, Write_ICL_Mesg, sizeof(IntContigLinkMesg) },
  { Read_AFG_Mesg, Write_AFG_Mesg, sizeof(AugFragMesg) },
  { Read_ISF_Mesg, Write_ISF_Mesg, sizeof(IntScaffoldMesg) },
  { Read_IMD_Mesg, Write_IMD_Mesg, sizeof(IntMateDistMesg) },
  { NULL,          NULL,           sizeof(IntAugFragMesg) },
  { Read_UTG_Mesg, Write_UTG_Mesg, sizeof(SnapUnitigMesg) },
  { Read_ULK_Mesg, Write_ULK_Mesg, sizeof(SnapUnitigLinkMesg) },
  { Read_ICM_Mesg, Write_ICM_Mesg, sizeof(IntConConMesg) },
  { Read_CCO_Mesg, Write_CCO_Mesg, sizeof(SnapConConMesg) },
  { Read_CLK_Mesg, Write_CLK_Mesg, sizeof(SnapContigLinkMesg) },
  { Read_SCF_Mesg, Write_SCF_Mesg, sizeof(SnapScaffoldMesg) },
  { Read_MDI_Mesg, Write_MDI_Mesg, sizeof(SnapMateDistMesg) },
  { Read_BAT_Mesg, Write_BAT_Mesg, sizeof(BatchMesg) },
  { Read_IBA_Mesg, Write_IBA_Mesg, sizeof(InternalBatchMesg) },
  { Read_BAC_Mesg, Write_BAC_Mesg, sizeof(BacMesg) },
  { Read_IBC_Mesg, Write_IBC_Mesg, sizeof(InternalBacMesg) },
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { NULL,NULL, sizeof(IntDegenerateScaffoldMesg)},
  { NULL,NULL, sizeof(SnapDegenerateScaffoldMesg)},
  { Read_SLK_Mesg,Write_SLK_Mesg, sizeof(SnapScaffoldLinkMesg)},
  { Read_ISL_Mesg,Write_ISL_Mesg, sizeof(InternalScaffoldLinkMesg)},
#ifdef AS_ENABLE_SOURCE
  { Read_FOM_Mesg, Write_FOM_Mesg, sizeof(FragOverlapMesg) },
#else
  { NULL,          NULL,           sizeof(FragOverlapMesg) },
#endif
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { NULL, NULL, 0l },
  { Read_EOF_Mesg, Write_EOF_Mesg, sizeof(EndOfFileMesg) }
};

MesgReader InputFileType_AS(FILE *fin)
{ int c;

  c = fgetc(fin);
  ungetc(c,fin);
  if (c == 0)
    return (ReadBinaryMesg_AS);

  return (ReadProtoMesg_AS);
}

int ReadBinaryMesg_AS(FILE *fin, GenericMesg **pmesg)
{ void (*reader)(FILE *, void *);

  static GenericMesg ReadMesg;
  static void *BiggestMesg;
  int firstime = (CDS_FTELL(fin) == 0);

  Mcode = 0;
  if (firstime)
    { int msp, code, t;

    //firstime = 0;
      msp = 0;
      for (t = 1; t <= NUM_OF_REC_TYPES; t++)
        if (CallTable[t].size > msp)
          msp = CallTable[t].size;
      BiggestMesg = (void *) malloc(msp);
      MessageNum  = 0;

      FREAD(&code,sizeof(code),1,fin);
      if (code != BIN_CODE)
        { fprintf(stderr,"ERROR: Input file is not a binary message file\n");
          exit (1);
        }
    }

  MessageNum += 1;
  *pmesg = &ReadMesg;

  MemTop = 0;
  if (fread(&Mcode,sizeof(Mcode),1,fin) != 1)
    return (EOF);

  if (!(Mcode >= 0 && Mcode < NUM_OF_REC_TYPES + 1)){
    fprintf(stderr,"* Mcode = %d\n", Mcode);
    assert(0);
  }

  FREAD(BiggestMesg,CallTable[Mcode].size,1,fin);
  reader = CallTable[Mcode].reader;
  if (reader != NULL)
    reader(fin,BiggestMesg);

  ReadMesg.t = (MessageType) Mcode;
  ReadMesg.m = BiggestMesg;
  ReadMesg.s = MemTop;
  return (0);
}

MesgWriter OutputFileType_AS(OutputType type)
{

  if(type == AS_BINARY_OUTPUT)
    return (WriteBinaryMesg_AS);

  return (WriteProtoMesg_AS);
}

int WriteBinaryMesg_AS(FILE *fout, GenericMesg *pmesg)
{ void (*writer)(FILE *, void *);
  int    code;

  int firstime = (CDS_FTELL(fout) == 0);

  Mcode = pmesg->t;
  if (firstime)
    { // firstime = 0;
      code = BIN_CODE;
      FWRITE(&code,sizeof(code),1,fout);
    }

  FWRITE(&Mcode,sizeof(Mcode),1,fout);
  FWRITE(pmesg->m,CallTable[pmesg->t].size,1,fout);
  writer = CallTable[pmesg->t].writer;
  if (writer != NULL)
    writer(fout,pmesg->m);
  return (0);
}

void ResetBinary_AS(void){
   free(MemBuffer);
   MemBuffer = NULL;   /* Memory allocation buffer for messages */
   MemMax = -1;
   MemTop = 0;;   /* Memory ceiling and current top */
}

