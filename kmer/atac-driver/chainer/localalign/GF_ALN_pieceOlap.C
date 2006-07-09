// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Clark Mobarry
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "GF_ALN_local.h"


typedef struct {
  char *aseg;
  char *bseg;
} PAIRALIGN;



/* safely copy a substring of a string into static space which is enlarged
   as needed */

static
int
safe_substr(char **seg,int *segspace,const char *seq,int beg,int end){

  if(*segspace<end-beg+1){
    *segspace=2*(end-beg)+1;
    *seg=(char*)ckrealloc(*seg,sizeof(char)*(*segspace));
  }
  strncpy(*seg,seq+beg,end-beg);
  (*seg)[end-beg]='\0';
  return(strlen(*seg)==end-beg);
}



/* construct a trace (with AS_ALN_OKNAlign) for the first local segment,
   copying the result into its own static location (so that we
   can call AS_ALN_OKNAlign on the second local segment without losing the
   result */

static int *get_trace(const char *aseq, const char *bseq,Local_Overlap *O,int piece,
                      int which){
  static char *aseg=NULL, *bseg=NULL;
  static int asegspace=0,bsegspace=0;
  static int *segtrace[2], tracespace[2]={0,0};
  int alen,blen;
  int spnt, *tmptrace;
  int segdiff;
  int i;
  int iret;

  assert(which==0||which==1);

  if(segtrace[which]==NULL){
    tracespace[which]=100;
    segtrace[which]=(int*)ckalloc(sizeof(int)*tracespace[which]);
  }

  iret = safe_substr(&aseg,&asegspace,aseq,O->chain[piece].piece.abpos,
                     O->chain[piece].piece.aepos);
  if(iret == 0){
    fprintf(stderr,"EXCEPTION get_trace: For aseg: len(aseg)=%d, len(bseg)=%d, alen=%d, blen=%d\n",
            strlen(aseg),strlen(bseg), alen,blen);
    return NULL;
  }
  iret = safe_substr(&bseg,&bsegspace,bseq,O->chain[piece].piece.bbpos,
                     O->chain[piece].piece.bepos);
  if(iret == 0){
    fprintf(stderr,"EXCEPTION get_trace: For bseg: len(aseg)=%d, len(bseg)=%d, alen=%d, blen=%d\n",
            strlen(aseg),strlen(bseg), alen,blen);
    return NULL;
  }

  alen=O->chain[piece].piece.aepos-O->chain[piece].piece.abpos;
  blen=O->chain[piece].piece.bepos-O->chain[piece].piece.bbpos;

  //printf("get_trace: len(aseg)=%d, len(bseg)=%d, alen=%d, blen=%d\n",
  //       strlen(aseg),strlen(bseg), alen,blen);

  /* get trace for the segment from AS_ALN_OKNAlign */
  spnt=0; 
  /* subtract because Gene likes to point to one before string start */
  aseg--; 
  bseg--;
  segdiff=(int)((O->chain[piece].piece.aepos-O->chain[piece].piece.abpos)
		*(1.5*O->chain[piece].piece.error)    +10);
  
  tmptrace=AS_ALN_OKNAlign(aseg,alen,bseg,blen,&spnt,segdiff);

  if(spnt!=0){
    if(spnt>0){
      O->chain[piece].agap+=spnt;
      O->chain[piece].piece.abpos+=spnt;
      i=0; 
      while(tmptrace[i]!=0){
	if(tmptrace[i]<0){
	  tmptrace[i]+=spnt;
	}
	i++;
      }
    } else {
      O->chain[piece].bgap+=-spnt;
      O->chain[piece].piece.bbpos+=-spnt;
      i=0; 
      while(tmptrace[i]!=0){
	if(tmptrace[i]>0){
	  tmptrace[i]+=spnt;
	}
	i++;
      }
    }
  }      
  aseg++; /* restore because need to know where memory block is allocated,
  		 and so that next time around strncpy will work right! */
  bseg++;
  i=0;
  while(tmptrace[i]!=0){
    segtrace[which][i]=tmptrace[i];
    i++;
    if(i==tracespace[which]){
      tracespace[which]*=2;
      segtrace[which]=(int*)ckrealloc(segtrace[which],
                                      sizeof(int)*tracespace[which]);
    }
  }
  segtrace[which][i]=0;
  return(segtrace[which]);

}







static void safe_add_to_seg(char **seg,int pos,char c,int *len){
  if(pos==*len){
    (*len)=(*len)*2;
    *seg=(char*)ckrealloc(*seg,sizeof(char)*((*len)+1));
  }
  (*seg)[pos]=c;
}



static PAIRALIGN *construct_pair_align(const char *aseq,
                                       const char *bseq,
                                       Local_Overlap *O,
                                       int piece,
                                       int *trace,
                                       int which){
  static char *aseg[2]={NULL,NULL},*bseg[2]={NULL,NULL};
  static int alen[2]={0,0},blen[2]={0,0};
  static PAIRALIGN pairalign[2];

  int starta,startb;
  int offseta,offsetb;
  int tpos,apos,bpos;


  if(aseg[which]==NULL){
    alen[which]=blen[which]=1000;
    aseg[which]=(char*)ckalloc((alen[which]+1)*sizeof(char));
    bseg[which]=(char*)ckalloc((blen[which]+1)*sizeof(char));
  }
  starta=offseta=O->chain[piece].piece.abpos;
  startb=offsetb=O->chain[piece].piece.bbpos;
  tpos=0;
  apos=0;
  bpos=0;

  while(trace[tpos]!=0){
    if(trace[tpos]<0){
      for(;offseta<-trace[tpos]+starta-1;apos++,offseta++){
	safe_add_to_seg(&(aseg[which]),apos,aseq[offseta],&(alen[which]));
      }
      safe_add_to_seg(&(aseg[which]),apos,'-',&(alen[which]));
      apos++;
    } else {
      for(;offsetb<trace[tpos]+startb-1;bpos++,offsetb++){
	safe_add_to_seg(&(bseg[which]),bpos,bseq[offsetb],&(blen[which]));
      }
      safe_add_to_seg(&(bseg[which]),bpos,'-',&(blen[which]));
      bpos++;
    }
    tpos++;
  }
  for(;offseta<O->chain[piece].piece.aepos;apos++,offseta++){
    safe_add_to_seg(&(aseg[which]),apos,aseq[offseta],&(alen[which]));
  }
  for(;offsetb<O->chain[piece].piece.bepos;bpos++,offsetb++){
    safe_add_to_seg(&(bseg[which]),bpos,bseq[offsetb],&(blen[which]));
  }

  assert(offseta==O->chain[piece].piece.aepos);
  assert(offsetb==O->chain[piece].piece.bepos);
  assert(offseta-O->chain[piece].piece.abpos+
	 offsetb-O->chain[piece].piece.bbpos+
	 tpos ==
	 apos+bpos);
  safe_add_to_seg(&(aseg[which]),apos,'\0',&(alen[which]));
  safe_add_to_seg(&(bseg[which]),bpos,'\0',&(blen[which]));

  pairalign[which].aseg=aseg[which];
  pairalign[which].bseg=bseg[which];
  return(pairalign+which);
}


static PAIRALIGN *get_align(const char *aseq,const char *bseq,Local_Overlap *O,int piece,
                            int which){
  int *trace=get_trace(aseq,bseq,O,piece,which);
  if(trace == NULL) return NULL;

  PAIRALIGN *pairalign = construct_pair_align(aseq,bseq,O,piece,trace,which);

  return(pairalign);
}












void fix_overlapping_pieces(const char *aseq, const char *bseq,
			    Local_Overlap *O,int piece0, int piece1){

  //int i,j;
  int pos=0;
  PAIRALIGN *pair_align1,*pair_align2;

  int offseta1,offsetb1,offseta2,offsetb2;
  int bestend1a,bestend1b,bestbeg2a,bestbeg2b;
  int into1,into2,bestinto2=0;
  int errs1,errs2,minerrs;


  assert(O->chain[piece0].piece.aepos>=O->chain[piece1].piece.abpos||
	 O->chain[piece0].piece.bepos>=O->chain[piece1].piece.bbpos);

  assert(O->chain[piece0].piece.aepos<=O->chain[piece1].piece.aepos);
  assert(O->chain[piece0].piece.bepos<=O->chain[piece1].piece.bepos);

  /* create alignments for the two segments */  

  pair_align1=get_align(aseq,bseq,O,piece0,0);
  pair_align2=get_align(aseq,bseq,O,piece1,1);

  if(pair_align1 == NULL || pair_align2 == NULL){
    fprintf(stderr,"EXCEPTION pair_align1=%p pair_align2=%p\n", pair_align1, pair_align2);
    fprintf(stderr,"EXCEPTION while fixing gap(%d,%d) (%d,%d)---(%d,%d) vs. gap(%d,%d)  (%d,%d)---(%d,%d)\n",
            O->chain[piece0].agap,	O->chain[piece0].bgap,
            O->chain[piece0].piece.abpos,O->chain[piece0].piece.bbpos,
            O->chain[piece0].piece.aepos,O->chain[piece0].piece.bepos,
            O->chain[piece1].agap,	O->chain[piece1].bgap,
            O->chain[piece1].piece.abpos,O->chain[piece1].piece.bbpos,
            O->chain[piece1].piece.aepos,O->chain[piece1].piece.bepos);
  }
  
  if(pair_align1 == NULL){
    fprintf(stderr,"EXCEPTION Fixing by pseudo-deleting piece0.\n");
    O->chain[piece0].agap=0;
    O->chain[piece0].bgap=0;
    if(piece0>0){
      O->chain[piece0].piece.abpos=O->chain[piece0-1].piece.aepos;
      O->chain[piece0].piece.aepos=O->chain[piece0-1].piece.aepos;
      O->chain[piece0].piece.bbpos=O->chain[piece0-1].piece.bepos;
      O->chain[piece0].piece.bepos=O->chain[piece0-1].piece.bepos;
    } else {
      O->chain[piece0].piece.abpos=0;
      O->chain[piece0].piece.aepos=0;
      O->chain[piece0].piece.bbpos=0;
      O->chain[piece0].piece.bepos=0;
    }
    O->chain[piece1].agap = O->chain[piece1].piece.abpos - O->chain[piece0].piece.aepos;
    O->chain[piece1].bgap = O->chain[piece1].piece.bbpos - O->chain[piece0].piece.bepos;
    return;
  }
  if(pair_align2 == NULL){
    fprintf(stderr,"EXCEPTION Fixing by pseudo-deleting piece1.\n");
    O->chain[piece1].agap=0;
    O->chain[piece1].bgap=0;
    O->chain[piece1].piece.abpos = O->chain[piece0].piece.aepos;
    O->chain[piece1].piece.aepos = O->chain[piece0].piece.aepos;
    O->chain[piece1].piece.bbpos = O->chain[piece0].piece.bepos;
    O->chain[piece1].piece.bepos = O->chain[piece0].piece.bepos;
    
    if(piece1+1<=O->num_pieces){
      O->chain[piece1+1].agap=O->chain[piece1+1].piece.abpos - O->chain[piece0].piece.aepos;
      O->chain[piece1+1].bgap=O->chain[piece1+1].piece.bbpos - O->chain[piece0].piece.bepos;
    }
    return;
  }


  /* if, in finding the alignments, we shift the ends of the
     alignment of the first segment to after the starts of the
     alignment of the second segment, then the overlap has been
     resolved, so we do nothing more */

  if(!(O->chain[piece0].piece.aepos>=O->chain[piece1].piece.abpos||
       O->chain[piece0].piece.bepos>=O->chain[piece1].piece.bbpos)){

    return;
  }

  /* if, in finding the alignments, we shift the end of the
     alignment of the second segment to before the start of the
     alignment of the first segment, then the second is contained
     in the first and we need to do something exceptional;
     the most heuristic, but consistent with the practice elsewhere
     in the local overlapper, is to pseudo-delete the second segment */

  if(!(O->chain[piece0].piece.aepos<=O->chain[piece1].piece.aepos)||
     !(O->chain[piece0].piece.bepos<=O->chain[piece1].piece.bepos)){

    O->chain[piece1].agap=0;
    O->chain[piece1].bgap=0;
    O->chain[piece1].piece.abpos=O->chain[piece0].piece.aepos;
    O->chain[piece1].piece.aepos=O->chain[piece0].piece.aepos;
    O->chain[piece1].piece.bbpos=O->chain[piece0].piece.bepos;
    O->chain[piece1].piece.bepos=O->chain[piece0].piece.bepos;
    
    if(piece1+1<=O->num_pieces){
      O->chain[piece1+1].agap=O->chain[piece1+1].piece.abpos-
        O->chain[piece0].piece.aepos;
      O->chain[piece1+1].bgap=O->chain[piece1+1].piece.bbpos-
        O->chain[piece0].piece.bepos;
    }
    return;
  }


  /* if, in finding the alignments, we shift the start of the
     alignment of the first segment to after the start of the
     alignment of the second segment, then the first is contained
     in the second and we need to do something exceptional;
     the most heuristic, but consistent with the practice elsewhere
     in the local overlapper, is to pseudo-delete the first segment */

  if(O->chain[piece0].piece.abpos>O->chain[piece1].piece.abpos||
     O->chain[piece0].piece.bbpos>O->chain[piece1].piece.bbpos){

    O->chain[piece0].agap=0;
    O->chain[piece0].bgap=0;
    if(piece0>0){
      O->chain[piece0].piece.abpos=O->chain[piece0-1].piece.aepos;
      O->chain[piece0].piece.aepos=O->chain[piece0-1].piece.aepos;
      O->chain[piece0].piece.bbpos=O->chain[piece0-1].piece.bepos;
      O->chain[piece0].piece.bepos=O->chain[piece0-1].piece.bepos;
    } else {
      O->chain[piece0].piece.abpos=0;
      O->chain[piece0].piece.aepos=0;
      O->chain[piece0].piece.bbpos=0;
      O->chain[piece0].piece.bepos=0;
    }
    O->chain[piece1].agap=O->chain[piece1].piece.abpos-
      O->chain[piece0].piece.aepos;
    O->chain[piece1].bgap=O->chain[piece1].piece.bbpos-
      O->chain[piece0].piece.bepos;

    return;
  }

  /* find start of region for evaluation in first alignment */
  /* when done,
     offseta1 and offsetb1 should be the offsets into the sequences
     such that they correspond to a column in the alignment of the
     first segment and that column contains the first possible
     overlap with the second segment */

  offseta1=O->chain[piece0].piece.abpos;
  offsetb1=O->chain[piece0].piece.bbpos;
  into1=0;
  while(offseta1<O->chain[piece1].piece.abpos&&
	offsetb1<O->chain[piece1].piece.bbpos){
    assert(pair_align1->aseg[into1]!='\0');
    assert(pair_align1->bseg[into1]!='\0');
    if(pair_align1->aseg[into1]!='-')offseta1++;
    if(pair_align1->bseg[into1]!='-')offsetb1++;
    into1++;
  }


  //  if(pair_align1->aseg[into1-1]!='-')offseta1--;
  //  if(pair_align1->bseg[into1-1]!='-')offsetb1--;

  /* count mismatches in the second alignment */

  into2=0;
  errs2=0;
  while(pair_align2->aseg[into2]!='\0'){
    assert(pair_align2->bseg[into2]!='\0');
    if(pair_align2->aseg[into2]!=pair_align2->bseg[into2]){
      errs2++;
    }
    into2++;
  }

  /* initialize solution variables and auxiliaries */
  into2=0;
  errs1 = (pair_align1->aseg[into1]!=pair_align1->bseg[into1] ? 1 : 0);
  minerrs=errs2;
  offseta2=O->chain[piece1].piece.abpos;
  offsetb2=O->chain[piece1].piece.bbpos;
  bestend1a=offseta1 - (pair_align1->aseg[into1-1]!='-' ? 1 : 0);
  bestend1b=offsetb1 - (pair_align1->bseg[into1-1]!='-' ? 1 : 0);
  bestbeg2a=offseta2;
  bestbeg2b=offsetb2;

  /* while there is potential overlap still to come ... */

  while(pair_align1->aseg[into1]!='\0'&&pair_align2->aseg[into2]!='\0'){

    // Once, we did the following assert, assuming that the alignment
    // of pair_align2 would not run out before pair_align1, since otherwise
    // there would be a containment or some such that shouldn't happen;
    // But, as luck would have it, alignment trimming quirks etc can
    // make it happen.  So ... no more assert
    //
    // assert(pair_align2->aseg[into2]!='\0');

    /* while a position in the second segment is no greater than
       the position in the first segment, 
       check for mismatch in second segment,
       counting errors,
       incrementing the sequence position counters as appropriate;
       advance the second segment
       position */
  
    while(offseta1>=offseta2||offsetb1>=offsetb2){
      errs2-= (pair_align2->aseg[into2]!=pair_align2->bseg[into2] ? 1 : 0);
      offseta2+= ( pair_align2->aseg[into2]!='-'  ? 1 : 0 );
      offsetb2+= ( pair_align2->bseg[into2]!='-'  ? 1 : 0 );
      into2++;
      if(pair_align2->aseg[into2]=='\0'){
	break;
      }
      //      assert(pair_align2->aseg[into2]!='\0');
      //      assert(pair_align2->bseg[into2]!='\0');
    }

    if(errs1+errs2<=minerrs&&
       pair_align1->aseg[into1]==pair_align1->bseg[into1]){
      minerrs=errs1+errs2;
      bestend1a=offseta1 /* -(pair_align1->aseg[into1-1]!='-' ? 1 : 0 )*/;
      bestend1b=offsetb1 /* -(pair_align1->bseg[into1-1]!='-' ? 1 : 0 )*/;
      bestbeg2a=offseta2;
      bestbeg2b=offsetb2;
      bestinto2=into2;
    }

    /* while the positions in the first segment are no greater than
       the positions in the second segment, 
       check for mismatch in first segment,
       counting errors,
       incrementing the sequence position counters as appropriate;
       advance the first segment
       position */
  
    while(offseta1<offseta2&&offsetb1<offsetb2){
      offseta1+= ( pair_align1->aseg[into1]!='-'  ? 1 : 0 );
      offsetb1+= ( pair_align1->bseg[into1]!='-'  ? 1 : 0 );
      into1++;
      errs1+= (pair_align1->aseg[into1]!=pair_align1->bseg[into1] ? 1 : 0);

      if(pair_align1->aseg[into1]=='\0'){
	break;
      }
    }
  }

  if(bestend1a<O->chain[piece0].piece.aepos)
    bestend1a++;
  if(bestend1b<O->chain[piece0].piece.bepos)
    bestend1b++;
  O->chain[piece0].piece.aepos=bestend1a;
  O->chain[piece0].piece.bepos=bestend1b;
  O->chain[piece1].piece.abpos=bestbeg2a;
  O->chain[piece1].piece.bbpos=bestbeg2b;
  O->chain[piece1].agap=bestbeg2a-bestend1a;
  O->chain[piece1].bgap=bestbeg2b-bestend1b;

  assert(O->chain[piece1].agap>=0);
  assert(O->chain[piece1].bgap>=0);
  assert(O->chain[piece1].agap==0||O->chain[piece1].bgap==0);


  // now, adjust the beginning of the second piece to skip any mismatches
  while(pair_align2->aseg[bestinto2]!=pair_align2->bseg[bestinto2]&&
	pair_align2->aseg[bestinto2]!='\0'){
    bestbeg2a += ( pair_align2->aseg[bestinto2]!='-'  ? 1 : 0 );
    bestbeg2b += ( pair_align2->bseg[bestinto2]!='-'  ? 1 : 0 );
    bestinto2++;
  }
  O->chain[piece1].piece.abpos=bestbeg2a;
  O->chain[piece1].piece.bbpos=bestbeg2b;
  O->chain[piece1].agap=bestbeg2a-bestend1a;
  O->chain[piece1].bgap=bestbeg2b-bestend1b;

  assert(O->chain[piece1].piece.abpos<=O->chain[piece1].piece.aepos);
  assert(O->chain[piece1].piece.bbpos<=O->chain[piece1].piece.bepos);
}
