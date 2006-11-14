
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "AS_global.h"
#include "CA_ALN_local.h"
#include "AS_ALN_aligners.h"
#include "AS_ALN_local.h"

//  Setting this to 4 causes a crash in
//     print_piece (O=0x57179d0, piece=0, aseq=0x5815c5f "", bseq=0x580c9cf "") at AS_ALN_loverlapper.c:186
//
#define DEBUG_LOCALOVL 0
#if DEBUG_LOCALOVL > 3
#define DEBUG_TRACE
#else
#undef DEBUG_TRACE
#endif

// to print calling params and global variables inside Local_Overlap_AS() ...
#undef DEBUG_PARAMS


#define CHECK_FINAL_QUALITY

//maximum number of matching segments that can be pieced together
int MaxGaps= 3;

//maximum allowed mismatch at end of overlap
int MaxBegGap= 200;

//maximum allowed mismatch at end of overlap
int MaxEndGap= 200;

//biggest gap internal to overlap allowed
int MaxInteriorGap=400;

//whether to treat the beginning of the b fragment and 
//   the end of the a fragment as allowed to have more error
int asymmetricEnds=0;

//amount of mismatch at end of the overlap that can cause
//   an overlap to be rejected
int MaxFreeFlap=20;

// set bubblePoppingVersion to 1 to get the original behavior
//    forcing positive ahang or set it to 0 to allow negative ahangs
int bubblePoppingVersion=0;

// set useSizeToOrderBlocks to 1 to get block mismatches resolved such that
//    the smaller block comes first:
//
//       .........AA------.........
//       .........--BBBBBB........
//
// or set it to 0 to get resolution with the A block always first
int useSizeToOrderBlocks = 1;

//global variable holding largest block mismatch of last returned overlap
int max_indel_AS_ALN_LOCOLAP_GLOBAL;


#define ENDHACK
#ifdef ENDHACK
  int ENDGAPHACK=3;
  int MINCORESEG=30;
#endif


#define MINUS_ONE 0   /* whether "zero" is one before the string or the start */

#undef AFFINE_QUALITY   /* overlap diff and length reported in affine terms */


/* safely copy a substring into memory pointed to by a static pointer */
static void safe_suffix(char **Dest,int *DestLen,char *src,int start){
  int len=strlen(src+start);
  if(*DestLen==0||*DestLen<len+1){
    *DestLen=len+1;
    if(*Dest==NULL){
      *Dest=(char*)ckalloc(sizeof(char)*(*DestLen));
    } else {
      *Dest=(char*)ckrealloc(*Dest,sizeof(char)*(*DestLen));
    }
  }
  strcpy(*Dest,src+start);
  assert(strlen(*Dest)==len);
}

/* print alignment of a "piece" -- one local alignment in the overlap chain*/
static void print_piece(Local_Overlap *O,int piece,char *aseq,char *bseq){
  int alen,blen,segdiff,spnt,epnt;
  static char *aseg,*bseg;
  static int aseglen=0,bseglen=0, *segtrace;

  alen=O->chain[piece].piece.aepos-O->chain[piece].piece.abpos;
  blen=O->chain[piece].piece.bepos-O->chain[piece].piece.bbpos;

    /* make sure there is (persistant) space for the strings */
    if(aseglen<alen+1){
      aseglen=2*(alen+1);
      if(aseg==NULL){
	aseg=(char*)ckalloc(sizeof(char)*aseglen);
      } else {
	aseg=(char*)ckrealloc(aseg,sizeof(char)*aseglen);
      }
    }
    if(bseglen<blen+1){
      bseglen=2*(blen+1);
      if(bseg==NULL){
	bseg=(char*)ckalloc(sizeof(char)*bseglen);
      } else {
	bseg=(char*)ckrealloc(bseg,sizeof(char)*bseglen);
      }
    }

    /* copy the segments */
    strncpy(aseg,aseq+O->chain[piece].piece.abpos,alen);
    aseg[alen]='\0';
    strncpy(bseg,bseq+O->chain[piece].piece.bbpos,blen);
    bseg[blen]='\0';

    assert(strlen(bseg)==blen);
    assert(strlen(aseg)==alen);

    /* guesstimate the required number of diagonals/edits to consider to
       get optimal alignment */
    segdiff=1+(int)((O->chain[piece].piece.aepos     
		     -O->chain[piece].piece.abpos)   
		    *1.5*O->chain[piece].piece.error);


    /* get trace for the segment from AS_ALN_OKNAlign */
    spnt=0; 
    /* subtract from aseg, bseg because Gene likes to index from 1, not 0 */
#ifdef OKNAFFINE
    epnt=0;
    segtrace=AS_ALN_OKNAffine(aseg-1,alen,bseg-1,blen,&spnt,&epnt,segdiff);
#else
    segtrace=AS_ALN_OKNAlign(aseg-1,alen,bseg-1,blen,&spnt,segdiff);
#endif

#if DEBUG_LOCALOVL > 4
    { int p=0,pos=0,neg=0;
    while(segtrace[p]!=0){
      printf("trace[%d] = %d\n",p,segtrace[p]);
      if(segtrace[p]>0){
	pos++;
      } else {
	neg++;
      }
      p++;
    }
      printf("alen %d blen %d pos %d neg %d spnt %d epnt %d\n",alen,blen,pos,neg,spnt,epnt);
    }
#endif

    if(epnt!=0){
      if(epnt<0){
	aseg[alen+epnt]='\0';
      } else {
	bseg[blen-epnt]='\0';
      }
    }
    if(spnt<0){
      int i=0;
      printf("(B sequence on top (2nd place in src); spnt=%d!!!\n",spnt);
      while(segtrace[i]!=0){segtrace[i++]*=-1;}
      PrintAlign(stdout,-spnt,0,bseg,aseg,segtrace);
      i=0; while(segtrace[i]!=0){segtrace[i++]*=-1;}
    } else {
      PrintAlign(stdout,spnt,0,aseg,bseg,segtrace);
    }
  

}


/* Create a trace to be interpreted as with DP_Compare_AS, but based
   on a Local_Overlap record.  A Local_Segment within the overlap
   will be aligned using AS_ALN_OKNAlign(), generating a subtrace.  Subtraces,
   with their indices appropriately adjusted, will be spliced together
   by an encoding of the gaps between segments; for now, we'll simply insert
   gaps as follows:

      A "gap" with x bases in A and y bases in B will become a section
      of the alignment x+y positions long, with the A fragment first
      and the B fragment second (with '=' indicating an aligned match):

               ====AAAAAAAAAA--------------======
               ====----------BBBBBBBBBBBBBB======

   Obviously, a more compact treatment is possible, but this makes clear
   the presumptive blocks involved in the mismatch; also, Karin says that it
   will make Consensus happy.  One slight change 

   Assumptions:  (do they matter?)
   By the usual conventions, the ahang should be nonnegative, the 
   bhang negative only if the ahang is positive, and both sequences
   should be in the forward orientation.

*/

int *AS_Local_Trace(Local_Overlap *O, char *aseq, char *bseq){
  static int *TraceBuffer=NULL;
  int i,j,k,segdiff,*segtrace;
  int lastgood=-1;
  static int allocatedspace=0;
  int tracespace=0;
  static char *aseg=NULL,*bseg=NULL;

  static int aseglen=0,bseglen=0;
  int abeg=0,bbeg=0; /* begining of segment; overloaded */
  int tracep=0; /* index into TraceBuffer */
  int spnt=0; /* to pass to AS_ALN_OKNAlign */
#ifdef OKNAFFINE
  int epnt=0; /* to pass to AS_ALN_OKNAffine */
#endif
  int alen,blen;
  int ahang = 0;


#define ADJUSTCHAINEND
#ifdef ADJUSTCHAINEND
  { int n=O->num_pieces;
    assert(O->num_pieces>0);
    O->chain[n].piece.abpos=O->chain[n].agap+O->chain[n-1].piece.aepos;
    O->chain[n].piece.bbpos=O->chain[n].bgap+O->chain[n-1].piece.bepos;
  }
    
#endif

  /*Estimate length required to store trace*/
  tracespace=0;
  tracespace+=abs(O->begpos)+abs(O->endpos);
  for(i=0;i<=O->num_pieces;i++){
    tracespace+=MAX(O->chain[i].agap,1);
    tracespace+=MAX(O->chain[i].bgap,1);
    tracespace+=(int)((O->chain[i].piece.aepos
		       -O->chain[i].piece.abpos)
		      *1.5*O->chain[i].piece.error);
    tracespace+=1000;
  }

  /*(Re)allocate space for the trace as necessary;
    Note that this is persistent storage so ...
       ... it doesn't need to get allocated on every call
       ... it shouldn't get freed
       ... results stored here need to be copied elsewhere if they
           are to be saved
  */
  if(allocatedspace<tracespace){
    allocatedspace=2*tracespace;
    if(TraceBuffer==NULL){
      TraceBuffer=(int*)ckalloc(sizeof(int)*allocatedspace);
    } else {
      TraceBuffer=(int*)ckrealloc(TraceBuffer,sizeof(int)*allocatedspace);
    }
  }
  
  /* for each Local_Overlap chain[i].piece, 
     need to handle the gap at the beginning and
     (for all but the final piece) the associated aligned segment */
  for(i=0;i<=O->num_pieces;i++){

#if DEBUG_LOCALOVL > 3
    if(i<O->num_pieces){
      printf("Top of main loop in AS_Local_Trace, segment %d alignment:\n",i);
      print_piece(O,i,aseq,bseq);
    }
#endif

    /* if conditions indicate the segment was deleted in previous loop,
       skip! */
    if(O->chain[i].agap==0 &&
       O->chain[i].bgap==0 &&
       O->chain[i].piece.abpos==O->chain[i].piece.aepos &&
       O->chain[i].piece.bbpos==O->chain[i].piece.bepos){
      continue;
    }

    /* guesstimate the required number of diagonals/edits to consider to
       get optimal alignment */
    segdiff=1+(int)((O->chain[i].piece.aepos     
		     -O->chain[i].piece.abpos)   
		    *1.5*O->chain[i].piece.error);



    /* Building an alignment/trace under the usual assumptions does not allow
       a given position in one sequence to simultaneously align to two or more
       positions in the other sequence.  However, the Find_Local_Overlap()
       routine can chain together local alignment segments that overlap.

       In order to make the local overlaps compatible with everything else
       we do, we need to trim back the overlaps.  
    
       Since we will "output" this segment at the end of the loop,
       we need to fix its overlap with the following segment in this
       cycle through the loop  
    */

#define WARNING_DELETE_PIECE

    k=i+1;
    while(k<O->num_pieces){

      /* if conditions indicate the segment was deleted previously, skip! */
      if(O->chain[k].agap==0 &&
	 O->chain[k].bgap==0 &&
	 O->chain[k].piece.abpos==O->chain[k].piece.aepos &&
	 O->chain[k].piece.bbpos==O->chain[k].piece.bepos){
	k++;
	continue;
      }

#if DEBUG_LOCALOVL > 4
fprintf(stdout,"Checking for overlaps of %d to %d\n",i,k);
    printf("Currently, segment %d alignment:\n",k);
    print_piece(O,k,aseq,bseq);

#endif

      if(O->chain[k].piece.abpos<O->chain[i].piece.aepos||
	O->chain[k].piece.bbpos<O->chain[i].piece.bepos){

#if DEBUG_LOCALOVL > 3
fprintf(stdout,"Need to fix overlaps of %d to %d\n",i,k);
fprintf(stdout," gap(%d,%d) (%d,%d)-(%d,%d) (piece %d)\n"
	       " gap(%d,%d) (%d,%d)-(%d,%d) (piece %d)\n",
	O->chain[i].agap,
	O->chain[i].bgap,
	O->chain[i].piece.abpos,O->chain[i].piece.bbpos,
	O->chain[i].piece.aepos,O->chain[i].piece.bepos,
	i,
	O->chain[k].agap,
	O->chain[k].bgap,
	O->chain[k].piece.abpos,O->chain[k].piece.bbpos,
	O->chain[k].piece.aepos,O->chain[k].piece.bepos,
	k);
#endif

	/* handle possibility of the first segment being
	 contained within the second;originally simply asserted
         against this; now, try to handle by deleting first segment */

	if(O->chain[i].piece.abpos>O->chain[k].piece.abpos||
	   O->chain[i].piece.bbpos>O->chain[k].piece.bbpos){

#if DEBUG_LOCALOVL > 0
#ifdef WARNING_DELETE_PIECE
	  fprintf(stdout,"Warning: deleting contained piece of local_overlap: a projection of (%d,%d)-(%d,%d) (piece %d) contains the same projection of (%d,%d)-(%d,%d) (piece %d)\n",
		  O->chain[k].piece.abpos,O->chain[k].piece.bbpos,
		  O->chain[k].piece.aepos,O->chain[k].piece.bepos,
		  k,
		  O->chain[i].piece.abpos,O->chain[i].piece.bbpos,
		  O->chain[i].piece.aepos,O->chain[i].piece.bepos,
		  i);
#endif
#endif

#undef FIXDELBYADVANCINGTONEXT
#ifdef FIXDELBYADVANCINGTONEXT
	  if(i>0){
	    O->chain[i].agap = O->chain[k].piece.abpos-
	      O->chain[lastgood].piece.aepos;
	    O->chain[i].bgap = O->chain[k].piece.bbpos-
	      O->chain[lastgood].piece.bepos;
	  } else {
	    O->chain[i].agap = O->chain[k].piece.abpos;
	    O->chain[i].bgap = O->chain[k].piece.bbpos;
	  }
	  O->chain[i].piece.abpos=O->chain[k].piece.abpos;
	  O->chain[i].piece.aepos=O->chain[k].piece.abpos;
	  O->chain[i].piece.bbpos=O->chain[k].piece.bbpos;
	  O->chain[i].piece.bepos=O->chain[k].piece.bbpos;
	  O->chain[k].agap=0;
	  O->chain[k].bgap=0;
#else
	  O->chain[i].agap=0;
	  O->chain[i].bgap=0;
	  if(lastgood>=0){
	    O->chain[i].piece.abpos=O->chain[lastgood].piece.aepos;
	    O->chain[i].piece.aepos=O->chain[lastgood].piece.aepos;
	    O->chain[i].piece.bbpos=O->chain[lastgood].piece.bepos;
	    O->chain[i].piece.bepos=O->chain[lastgood].piece.bepos;
	  } else {
	    O->chain[i].piece.abpos=0;
	    O->chain[i].piece.aepos=0;
	    O->chain[i].piece.bbpos=0;
	    O->chain[i].piece.bepos=0;
          }
	  O->chain[k].agap=O->chain[k].piece.abpos-
	    O->chain[i].piece.aepos;
	  O->chain[k].bgap=O->chain[k].piece.bbpos-
	    O->chain[i].piece.bepos;
	  if(lastgood<0){
//printf("Shrinking gaps for segment %d\n",k);
	    O->chain[k].agap--;
	    O->chain[k].bgap--;
	  }
#endif

        } else 	/* otherwise, check for 2nd piece contained within first */
	  if(O->chain[i].piece.aepos>O->chain[k].piece.aepos||
	   O->chain[i].piece.bepos>O->chain[k].piece.bepos){

#define FIX_OLAP_OFFSET 0 /* hack to handle different position systems*/

	/* if the next piece is completely within current piece, 
	   effectively remove it */


#if DEBUG_LOCALOVL > 0
#ifdef WARNING_DELETE_PIECE
fprintf(stdout,"Warning: deleting contained piece of local_overlap: a projection of (%d,%d)-(%d,%d) (piece %d) contains the same projection of (%d,%d)-(%d,%d) (piece %d)\n",
	O->chain[i].piece.abpos,O->chain[i].piece.bbpos,
	O->chain[i].piece.aepos,O->chain[i].piece.bepos,
	i,
	O->chain[k].piece.abpos,O->chain[k].piece.bbpos,
	O->chain[k].piece.aepos,O->chain[k].piece.bepos,
	k);
#endif
#endif

#ifdef FIXDELBYADVANCINGTONEXT
	  O->chain[k].agap = 0;
	  O->chain[k].piece.abpos=O->chain[i].piece.aepos;
	  O->chain[k].piece.aepos=O->chain[i].piece.aepos;
	  O->chain[k].bgap = 0;
	  O->chain[k].piece.bbpos=O->chain[i].piece.bepos;
	  O->chain[k].piece.bepos=O->chain[i].piece.bepos;
	  if(k+1<=O->num_pieces){
	    O->chain[k+1].agap=O->chain[k+1].piece.abpos-
	      O->chain[k-1].piece.aepos;
	    O->chain[k+1].bgap=O->chain[k+1].piece.bbpos-
	      O->chain[k-1].piece.bepos;
	  }
#else
	  O->chain[k].agap = 0;
	  O->chain[k].bgap = 0;
	  O->chain[k].piece.abpos=O->chain[i].piece.aepos;
	  O->chain[k].piece.aepos=O->chain[i].piece.aepos;
	  O->chain[k].piece.bbpos=O->chain[i].piece.bepos;
	  O->chain[k].piece.bepos=O->chain[i].piece.bepos;
	  if(k+1<=O->num_pieces){
	    int l;
	    l=k-1;
	    while(O->chain[l].agap==0 &&
		  O->chain[l].bgap==0 &&
		  O->chain[l].piece.abpos==O->chain[l].piece.aepos &&
		  O->chain[l].piece.bbpos==O->chain[l].piece.bepos){
	      l--;
	      assert(l>=0);
	    }
	      
#if DEBUG_LOCALOVL > 3
	    fprintf(stdout,"Resetting gaps for segment %d to point to end of segment %d\n",
		    k+1,l);
#endif

	    O->chain[k+1].agap=O->chain[k+1].piece.abpos-
	      O->chain[l].piece.aepos;
	    O->chain[k+1].bgap=O->chain[k+1].piece.bbpos-
	      O->chain[l].piece.bepos;
	  }

#endif

	/* else, fix the overlap */
	} else {
	
#if DEBUG_LOCALOVL > 3
	  fprintf(stdout,"Fixing overlap; was:\n"
		  "gap(%d,%d) (%d,%d)------(%d,%d) (piece %d)\n"
		  "gap(%d,%d) (%d,%d)------(%d,%d) (piece %d)\n",
		  O->chain[i].agap,O->chain[i].bgap,
		  O->chain[i].piece.abpos,O->chain[i].piece.bbpos,
		  O->chain[i].piece.aepos,O->chain[i].piece.bepos,
		  i,
		  O->chain[k].agap,O->chain[k].bgap,
		  O->chain[k].piece.abpos,O->chain[k].piece.bbpos,
		  O->chain[k].piece.aepos,O->chain[k].piece.bepos,
		  k);
#endif
	  fix_overlapping_pieces(aseq+FIX_OLAP_OFFSET,
				 bseq+FIX_OLAP_OFFSET,
				 O,i,k);


	  // if the second piece disappeared
	  if(O->chain[k].piece.abpos==O->chain[k].piece.aepos||
	     O->chain[k].piece.bbpos==O->chain[k].piece.bepos){

#if DEBUG_LOCALOVL > 3
	    fprintf(stdout,"Second segment size now 0!\n");
#endif
	    O->chain[k].agap = 0;
	    O->chain[k].bgap = 0;
	    O->chain[k].piece.abpos=O->chain[i].piece.aepos;
	    O->chain[k].piece.aepos=O->chain[i].piece.aepos;
	    O->chain[k].piece.bbpos=O->chain[i].piece.bepos;
	    O->chain[k].piece.bepos=O->chain[i].piece.bepos;
	    if(k+1<=O->num_pieces){
	      int l;
	      l=k-1;
	      while(O->chain[l].agap==0 &&
		    O->chain[l].bgap==0 &&
		    O->chain[l].piece.abpos==O->chain[l].piece.aepos &&
		    O->chain[l].piece.bbpos==O->chain[l].piece.bepos){
		l--;
		assert(l>=0);
	      }
	      
#if DEBUG_LOCALOVL > 3
	      fprintf(stdout,"Resetting gaps for segment %d to point to end of segment %d\n",
		      k+1,l);
#endif

	      O->chain[k+1].agap=O->chain[k+1].piece.abpos-
		O->chain[l].piece.aepos;
	      O->chain[k+1].bgap=O->chain[k+1].piece.bbpos-
		O->chain[l].piece.bepos;
	    }
	  } else {
	    // if the first piece disappeared
	    if (O->chain[i].piece.abpos==O->chain[i].piece.aepos||
		O->chain[i].piece.bbpos==O->chain[i].piece.bepos){
#ifdef FIXDELBYADVANCINGTONEXT
	      assert(0)/*uncovered situation; cut and paste from above to fix */
#endif

#if DEBUG_LOCALOVL > 3
		fprintf(stdout,"First segment size now 0!\n");
#endif
		O->chain[i].agap=0;
	      O->chain[i].bgap=0;
	      if(lastgood>=0){
		O->chain[i].piece.abpos=O->chain[lastgood].piece.aepos;
		O->chain[i].piece.aepos=O->chain[lastgood].piece.aepos;
		O->chain[i].piece.bbpos=O->chain[lastgood].piece.bepos;
		O->chain[i].piece.bepos=O->chain[lastgood].piece.bepos;
	      } else {
		O->chain[i].piece.abpos=0;
		O->chain[i].piece.aepos=0;
		O->chain[i].piece.bbpos=0;
		O->chain[i].piece.bepos=0;
	      }
	      
	      { // need to adjust gaps for all segments between i and k (inclusive of k)
		int l;
		for(l=i+1;l<=k;l++){
		  O->chain[l].agap=O->chain[l].piece.abpos-
		    O->chain[i].piece.aepos;
		  O->chain[l].bgap=O->chain[l].piece.bbpos-
		    O->chain[i].piece.bepos;
		  if(lastgood<0){
		    //printf("Shrinking gaps for segment %d\n",l);
		    O->chain[l].agap--;
		    O->chain[l].bgap--;
		  }
		}
	      }
	    }
	  }
	}


#if DEBUG_LOCALOVL > 3
	  fprintf(stdout,"now:\n"
		  "gap(%d,%d) (%d,%d)------(%d,%d) (piece %d)\n"
		  "gap(%d,%d) (%d,%d)------(%d,%d) (piece %d)\n",
		  O->chain[i].agap,O->chain[i].bgap,
		  O->chain[i].piece.abpos,O->chain[i].piece.bbpos,
		  O->chain[i].piece.aepos,O->chain[i].piece.bepos,
		  i,
		  O->chain[k].agap,O->chain[k].bgap,
		  O->chain[k].piece.abpos,O->chain[k].piece.bbpos,
		  O->chain[k].piece.aepos,O->chain[k].piece.bepos,
		  k);
#endif

      }

      k++;
    }
  
    /* if conditions indicate the segment was deleted previously, skip! */
    if(O->chain[i].agap==0 &&
       O->chain[i].bgap==0 &&
       O->chain[i].piece.abpos==O->chain[i].piece.aepos &&
       O->chain[i].piece.bbpos==O->chain[i].piece.bepos){
      continue;
    }
  
    /* set up positions before which gaps are inserted to handle
       the gap portion of a chain piece */

    /* put gaps before beginning of aligned piece (but after the portion
       of aseq in the gap); location is relative to the beginning of the
       alignment (i.e., ignores ahang worth of positions)  */

    if(i!=O->num_pieces){
      abeg=O->chain[i].piece.abpos;
    } else {
      assert(lastgood>=0&&lastgood<O->num_pieces);
      abeg=O->chain[lastgood].piece.aepos+
	O->chain[i].agap + MINUS_ONE ;
    }
  
    /*handle boundary case to prevent gaps preceding the b sequence*/
    if((i==0||lastgood<0)&&O->chain[i].bgap>0){
      assert(O->chain[i].agap>=0);

      #undef DEBUG_AHANG_CHANGE
      #ifdef DEBUG_AHANG_CHANGE
            printf("Reset ahang, abeg, agap and bgap\n");
            printf("originally: %d %d %d %d\n",
      	     O->begpos,
      	     abeg,
      	     O->chain[i].agap,
      	     O->chain[i].bgap);
      #endif

      if(bubblePoppingVersion){
	  O->begpos=O->chain[i].piece.abpos-1;
	  assert(O->begpos>=0);
	  O->chain[i].agap=MINUS_ONE;
	  O->chain[i].bgap=O->chain[i].piece.bbpos-1;
      }else{
	if(O->begpos>=0){
	  O->begpos=O->chain[i].piece.abpos-1;
	  assert(O->begpos>=0);
	  O->chain[i].agap=MINUS_ONE;

#define WE_BELIEVE_BGAP_IS_RIGHT
#ifdef WE_BELIEVE_BGAP_IS_RIGHT
	  assert( ( i==0&&	O->chain[i].bgap==O->chain[i].piece.bbpos-1)
	  	  ||( i>0&&lastgood<0&&O->chain[i].bgap==O->chain[i].piece.bbpos-1))
	  	    ;
#else
	  if ( ! 	  ( ( i==0&&	O->chain[i].bgap==O->chain[i].piece.bbpos-1)
			    ||( i>0&&lastgood<0&&O->chain[i].bgap==O->chain[i].piece.bbpos-1))) {
	    fprintf(stderr,"LOVERLAPPER: Problem adjusting first bgap\n");
	    fprintf(stderr,"lastgood %d piece %d bgap %d bbpos %d begpos %d\n",
		    lastgood, 
		    i, 
		    O->chain[i].bgap,
		    O->chain[i].piece.bbpos,
		    O->begpos);
	    fprintf(stderr,">Aseq\n%s\n",aseq+1);
	    fprintf(stderr,">Bseq\n%s\n",bseq+1);
	  }
#endif


	  if(lastgood<0){
	    O->chain[i].bgap=O->chain[i].piece.bbpos-1;
	  }
	} else {
	  if(i==0){
	    O->begpos-=O->chain[i].bgap;
	  } else{
	    O->begpos=-O->chain[i].bgap;
	  }
	  O->chain[i].bgap=0;
	}
      }

      #ifdef DEBUG_AHANG_CHANGE
            printf("afterwards: %d %d %d %d\n",
	     O->begpos,
	     abeg,
	     O->chain[i].agap,
	     O->chain[i].bgap);
      #endif
    }

    /* now prevent gaps at end of A sequence */
//    if(i==O->num_pieces&&O->endpos>=0){
//      O->endpos+=O->chain[i].bgap;
//      O->chain[i].bgap=0;
//    }

    /* now make sure that end mismatches are treated by tucking
       the shorter tail into a gap before the longer tail,
       or, ifdef FORCEPOSITIVEBHANG, by tucking the A tail into
       a gap before the B tail */
    if(i==O->num_pieces){
      if(O->endpos>=0){
	O->endpos+=O->chain[i].bgap;
	O->chain[i].bgap=0;
      }else {
#undef FORCEPOSITIVEBHANG
#ifndef FORCEPOSITIVEBHANG
	O->endpos-=O->chain[i].agap;
	abeg-=O->chain[i].agap;
	O->chain[i].agap=0;
#else
	O->chain[i].agap-=O->endpos;
	O->endpos=O->chain[i].bgap;
	O->chain[i].bgap=0;
#endif
      }	
    }

    if(i==0||lastgood<0){
      ahang=O->begpos;
    }

    /* put gaps before the portion of bseq in the gap; for the first
       piece, this means before position 0 */
    if(i==0||lastgood<0){
//      bbeg=1-min(ahang,0);
      bbeg=1-min(ahang,0);
    } else {
      assert(lastgood<O->num_pieces);
      bbeg=O->chain[lastgood].piece.bepos /*-min(ahang,0)*/;
    }



    /* now insert the right number of gaps! */


    ///////////////////////////////////////


    if(i==O->num_pieces){

      if(O->endpos<0){
	O->chain[i].agap+=-O->endpos;
	O->endpos=0;
      } else {
	O->chain[i].bgap+=O->endpos;
	O->endpos=0;
      }

      if(O->chain[i].agap <= O->chain[i].bgap){
	O->endpos=O->chain[i].bgap;
	O->chain[i].bgap=0;
      }else{
	O->endpos=-O->chain[i].agap;
	O->chain[i].agap=0;
      }

    }

    if(O->chain[i].agap <= O->chain[i].bgap || ! useSizeToOrderBlocks ){

      /* start by putting len(agap) gaps before the chunk of B in the gap */

      for(j=0;j<O->chain[i].agap
#if MINUS_ONE == 1
	    -(i==0)   /* this offset deals with possible bug in Local_Overlap*/
	    +(i==O->num_pieces&&
	      O->chain[i].agap>0) /* and this one too! */
#endif
	    ;j++){
	TraceBuffer[tracep++]=bbeg;
      }


      /* then put len(bgap) gaps before the chunk of A in the gap */

      for(j=0;j<O->chain[i].bgap
#if MINUS_ONE == 1
	    -(i==0)   /* this offset deals with possible bug in Local_Overlap*/
	    +(i==O->num_pieces&&
	      O->chain[i].bgap>0) /* and this one too! */
#endif
	    ;j++){
	TraceBuffer[tracep++]=-abeg;
      }


    } else { // if the bgap is smaller,

      abeg-=O->chain[i].agap;
      bbeg+=O->chain[i].bgap;

      //      fprintf(stderr,"SWAPPING BLOCK MISMATCH GAPS\n");

      /* start by putting len(bgap) gaps before the chunk of A in the gap */

      for(j=0;j<O->chain[i].bgap   ;j++){
	TraceBuffer[tracep++]=-abeg;
      }

      /* then put len(agap) gaps before the chunk of B in the gap */

      for(j=0;j<O->chain[i].agap ;j++){
	TraceBuffer[tracep++]=bbeg;
      }

    }

    ///////////////////////////////////////  

    /* if last piece, there is no aligned segment */

    if(i==O->num_pieces)break;


    /* set bbeg to beginning of aligned segment for piece */

    abeg=O->chain[i].piece.abpos ;
    bbeg=O->chain[i].piece.bbpos /* -min(ahang,0) ??*/;


    /* set lengths of segments */

    alen=O->chain[i].piece.aepos-abeg; /* check +1?? */
    blen=O->chain[i].piece.bepos-bbeg; /* check +1?? */


    /* create strings for just the parts of the sequences in the
       aligned segment */


    /* make sure there is (persistant) space for the strings */
    if(aseglen<alen+1){
      aseglen=2*(alen+1);
      if(aseg==NULL){
	aseg=(char*)ckalloc(sizeof(char)*aseglen);
      } else {
	aseg=(char*)ckrealloc(aseg,sizeof(char)*aseglen);
      }
    }
    if(bseglen<blen+1){
      bseglen=2*(blen+1);
      if(bseg==NULL){
	bseg=(char*)ckalloc(sizeof(char)*bseglen);
      } else {
	bseg=(char*)ckrealloc(bseg,sizeof(char)*bseglen);
      }
    }


    /* copy the segments */
    strncpy(aseg,aseq+abeg,alen);
    aseg[alen]='\0';
    strncpy(bseg,bseq+bbeg,blen);
    bseg[blen]='\0';

    assert(strlen(bseg)==blen);
    assert(strlen(aseg)==alen);

    /* guesstimate the required number of diagonals/edits to consider to
       get optimal alignment */
    segdiff=1+(int)((O->chain[i].piece.aepos     
		     -O->chain[i].piece.abpos)   
		    *1.5*O->chain[i].piece.error);


    /* get trace for the segment from AS_ALN_OKNAlign */
    spnt=0; 
    /* subtract from aseg, bseg because Gene likes to index from 1, not 0 */
#ifdef OKNAFFINE
    epnt=0;
    segtrace=AS_ALN_OKNAffine(aseg-1,alen,bseg-1,blen,&spnt,&epnt,segdiff);

    if(epnt!=0){
      if(epnt>0){ /* throwing away some of B segment */
	O->chain[i+1].bgap+=epnt;
	O->chain[i].piece.bepos-=epnt;
	assert(O->chain[i].piece.bbpos<=O->chain[i].piece.bepos);
      } else {
	O->chain[i+1].agap-=epnt;
	O->chain[i].piece.aepos+=epnt;
	assert(O->chain[i].piece.abpos<=O->chain[i].piece.aepos);
      }
    }
#else
    segtrace=AS_ALN_OKNAlign(aseg-1,alen,bseg-1,blen,&spnt,segdiff);
#endif
    if(spnt>0){
      O->chain[i].agap+=spnt;
      O->chain[i].piece.abpos+=spnt;
    } else {
      O->chain[i].bgap-=spnt;
      O->chain[i].piece.bbpos-=spnt;
    }


#ifdef DEBUG_TRACE
    if(spnt<0){
      int i=0;
      printf("(B sequence on top (2nd place in src); spnt=%d!!!\n",spnt);
      while(segtrace[i]!=0){segtrace[i++]*=-1;}
      PrintAlign(stdout,-spnt,0,bseg,aseg,segtrace);
      i=0; while(segtrace[i]!=0){segtrace[i++]*=-1;}
    } else {
      PrintAlign(stdout,spnt,0,aseg,bseg,segtrace);
    }
#endif

    assert(segtrace!=NULL);


#undef PRINT_SUMMARY
#ifdef PRINT_SUMMARY
    fprintf(stdout,"  Match A[%d,%d] to B[%d,%d]\n",
	    O->chain[i].piece.abpos,
	    O->chain[i].piece.aepos,
	    O->chain[i].piece.bbpos,
	    O->chain[i].piece.bepos);
#endif


    /* Now copy the segment trace into master trace, adjusting positions */
    j=0;

    if(spnt<0){
      int ctr;
      for(ctr=0;ctr<abs(spnt);ctr++){
	TraceBuffer[tracep++]=-abeg;
      }
    } else {
      int ctr;
      for(ctr=0;ctr<spnt;ctr++){
	TraceBuffer[tracep++]=bbeg;
      }
    }

#ifdef PRINT_SUMMARY
    fprintf(stdout,"\t");
#endif
    while(segtrace[j]!=0){
#ifdef PRINT_SUMMARY
      fprintf(stdout," %d",segtrace[j]+(segtrace[j]<0 ? 1 +MAX(0,spnt) : -1-MAX(0,-spnt)));
#endif
      if(segtrace[j]<0){
	TraceBuffer[tracep++]=-abeg+segtrace[j++]+1 /* -MAX(0,spnt) ?? */;
      } else {
	TraceBuffer[tracep++]=bbeg+segtrace[j++]-1 /* +MAX(0,-spnt) ?? */;
      }
    }

#ifdef PRINT_SUMMARY
      fprintf(stdout,"\n");
#undef VERBOSE_SUMMARY
#ifdef VERBOSE_SUMMARY
    if(spnt<0){
      int i=0;
      printf("(B sequence on top (2nd place in src); spnt=%d!!!\n",spnt);
      while(segtrace[i]!=0){segtrace[i++]*=-1;}
      PrintAlign(stdout,-spnt,0,bseg,aseg,segtrace);
      i=0; while(segtrace[i]!=0){segtrace[i++]*=-1;}
    } else {
      PrintAlign(stdout,spnt,0,aseg,bseg,segtrace);
    }
#endif
#endif


    /* set lastgood to this segment */

    lastgood=i;

    /* and back to the top of the loop for another overlap piece */
  }		      

  /* terminate the trace */
  TraceBuffer[tracep]=0;

  /* really non-robust coding on allocated(trace)space, so let's at least test
     for a violation! */
  assert(tracep<allocatedspace);

#if 0
  for(i=0;i<tracep;i++){
    printf("%d .",TraceBuffer[i]);
  }
  printf("\n");
#endif
#undef DEBUG_TRACE_CONSTR
#ifdef DEBUG_TRACE_CONSTR
  if(O->begpos>=0){
    PrintAlign(stdout,O->begpos,O->endpos,
	       aseq+1,bseq+1,TraceBuffer);
  }else{
    for(i=0;i<tracep;i++){
      TraceBuffer[i]*=-1;
    }
    PrintAlign(stdout,-O->begpos,-O->endpos,
	       bseq+1,aseq+1,TraceBuffer);
    for(i=0;i<tracep;i++){
      TraceBuffer[i]*=-1;
    }
  }
  printf("End of Trace construction!\n");
#endif

  return(TraceBuffer);
}


/* Given fragments a and b, find the best overlap between them using
   local alignment code.  The purpose is to permit "bubbles" in the fragment
   graph due to multi-base polymorphisms to be smoothed out, such that a 
   unitig can be constructed.

   The method relies on Myers' local overlap code.  MORE DETAILS????

   The local overlap code returns the best overlap available as a chain of
   matching segments.

   The current function is concerned with three things:
   
   1) Providing a wrapper around the local alignment code so that users
   can carry out operations on the level of the "message" (IFM, OVL, etc).

   2) Determine whether the best possible overlap is good enough to constitute
   an overlap for the purposes of unitigging.  A preliminary set of criteria
   are given below.

   3) Combine bits of existing code to provide a full alignment specification
   for the local overlapper; one can debate whether large gaps should be 
   encoded in an edit trace, but for the purposes of integration with the
   existing code base, the current function will aim to specify an
   alignment in terms of an ahang, a bhang and a trace, which will get
   converted into the core elements of an overlap message just as in
   DP_Compare_AS.  

   PRELIMINARY DECISION CRITERIA:

     - if an overall %-match of less than erate with length MinLen[=40]
       exists, then the overlap should be allowed regardless of the
       following additional criteria
     - at least one segment of at least MinLen bases
     - if we treat an affine event as a single error (and the length involved
       to likewise be a single base) for the purpose of calculating an
       error rate, the error rate should be no more than erate.
     - no more than MaxGaps [=2?] unaligned regions ("gaps") in the overlap;
       this also entails no more than MaxGaps+1 matching segments
       
     Additional criteria that we might want to use:

     - no more than MaxEndGap [=20?] unaligned bases at either end; this
       is to prevent a classic branchpoint out of a repeat into distinct
       regions from being overlapped; on the other hand, the requirement
       that all fragments in a bubble go back together might mean that
       this filter is unnecessary.
     - only one end may end in a gap; having both ends
       be in polymorphic regions should be quite unlikely, and it is
       not currently our concern to handle crappy fragment ends (= spurs,
       handled by other means).


   AS_FIND_LOCAL_OVERLAP:
      Just find a good alignment to the boundary of the d.p. matrix.  From
      this extrapolate a rough overlap relationship without an alignment
      and return the result (if there is one).
   AS_FIND_LOCAL_ALIGN:
      For this option, further go to the trouble of computing the alignment
      and store it in the overlap message.
   AS_FIND_LOCAL_ALIGN_NO_TRACE
      A hack by SAK...don't compute delta encoding.

   As with DP_Compare_AS, the resulting overlap message occupies
   memory which is persistently owned (reused) by the function; it
   must be copied if it is to be retained beyond the given call.

*/


OverlapMesg *Local_Overlap_AS(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where)
{ int   alen,  blen;
  char *aseq, *bseq;
  int NumSegs=0;
  int seglen,sumseglen=0;
  int noCoreSeg=1;
  int i;
  int forcenext=0;

  static char *Ausable=NULL, *Busable=NULL;
  static int AuseLen=0, BuseLen=0;

#ifdef ENDHACK
  int coreseglen=min(MINCORESEG,minlen);
#else
  int coreseglen=minlen;
#endif

  double avgerror=0.;

  static OverlapMesg QVBuffer;

  Local_Segment *local_results=NULL;
  Local_Overlap *O=NULL;
  
  assert(what == AS_FIND_LOCAL_ALIGN || what == AS_FIND_LOCAL_ALIGN_NO_TRACE
	 || what == AS_FIND_LOCAL_OVERLAP );

  aseq = a->sequence;  /* Setup sequence access */
  bseq = b->sequence;
  alen = strlen(aseq);
  blen = strlen(bseq);


#ifdef DEBUG_PARAMS
  fprintf(stderr,
	  "Local_Overlap_AS params:\n"
	  "MaxGaps %d MaxBegGap %d MaxEndGap %d MaxInteriorGap %d asymmetricEnds %d\n" 
	  "MaxFreeFlap %d bubblePoppingVersion %d useSizeToOrderBlocks %d\n"
	  "beg %d end %d opposite %d erate %f minlen %d\n"
	  ">a\n%s\n>b\n%s\n",
	  MaxGaps,MaxBegGap,MaxEndGap,MaxInteriorGap,asymmetricEnds,
	  MaxFreeFlap,bubblePoppingVersion,useSizeToOrderBlocks,
	  beg,end,opposite,erate,minlen,
	  aseq,bseq);
	  
#endif



  aseq -= MINUS_ONE;
  bseq -= MINUS_ONE;

  if (opposite){                 /* Compare in opposite orientation. */
    Complement_Fragment_AS(b);
  }


  /* now generate the substrings of the input sequences which contain the
     portion usable according to the user's begin/end parameters */
#define EXPAND_BAND 10 /* amount to pad ends of  user's begin..end range */
#define STRETCH 1.5  /* amount of "stretching" allowed in bounding 
			overlap end */
#define BIGPAD 500 /* amount of extra slop allowed in bounding overlap end */


  // not really using the "suffix" part of next calls; just a safe string copy
  safe_suffix(&Ausable,&AuseLen,aseq+MINUS_ONE,0);
  safe_suffix(&Busable,&BuseLen,bseq+MINUS_ONE,0);

  for(i=0;i<MAX(0,beg-EXPAND_BAND);i++){
    Ausable[i]='N';
  }
  assert(alen>=i);
  for(i=0;i<MAX(0,-end-EXPAND_BAND);i++){
    Busable[i]='N';
  }
  assert(blen>=i);

  for(i=MAX(0,end)+(blen+min(0,end))*STRETCH+BIGPAD;i<alen;i++){
    Ausable[i]='N';
  }

  i=-min(0,beg)+(alen-MAX(0,beg))*STRETCH+BIGPAD;
  assert(i>=0);

  for(;i<blen;i++){
    Busable[i]='N';
  }


  /* Notes on handling of reverse complement overlap searching:
     1) The implementation here uses Find_Local_Segments and 
        Find_Local_Overlaps only in the forward direction; the interaction
	of these two routines on reverse orientation searches is non-obvious
	[even Gene agrees], so we avoid it.
     2) Thus, it is the responsibility of this routine to reverse the B
        sequence prior to running these routines.
  */



  local_results=Find_Local_Segments(Ausable-MINUS_ONE,strlen(Ausable),
				    Busable-MINUS_ONE,strlen(Busable),
				    LOCAL_FORW, 8, erate*2., &NumSegs);

#if DEBUG_LOCALOVL > 0
  fprintf(stdout,"Find_Local_Segments found %d.\n",NumSegs);
#if DEBUG_LOCALOVL > 2
  { int i;
  for(i=0;i<NumSegs;i++){
    printf("Segment %d: (%d,%d)[----------](%d,%d), error %e\n",i,
	   local_results[i].abpos,
	   local_results[i].bbpos,
	   local_results[i].aepos,
	   local_results[i].bepos,
	   local_results[i].error);
  }
  }
#endif
#endif

  if(NumSegs==0){
    #if DEBUG_LOCALOVL > 0 
    fprintf(stdout,"No segments found by Find_Local_Segments()\n");
    #endif
    goto nooverlap;    
  }

  O=Find_Local_Overlap(alen,blen,
		       0 /*comp==0 -> fwd orientation */,
		       0 /*nextbest==0 -> best overlap*/,
		       local_results,NumSegs,
		       minlen
#ifdef ENDHACK
		       -2*ENDGAPHACK
#endif
		        ,1.);



#if DEBUG_LOCALOVL > 3
  if(O!=NULL){
    for(i=0;i<O->num_pieces;i++){
      printf("overlap segment %d: (%d,%d)[---------](%d,%d)\n",
	     i,
	     O->chain[i].piece.abpos,
	     O->chain[i].piece.bbpos,
	     O->chain[i].piece.aepos,
	     O->chain[i].piece.bepos
	     );
    }
  }
#endif




  if(O!=NULL){
  next_overlap:
#ifndef CHECK_FINAL_QUALITY
    while(O->begpos > end||O->begpos<beg||O->indif/(double)O->length>erate||forcenext){
#else

#undef NONFINAL_BAND_TEST
#ifdef NONFINAL_BAND_TEST
    while(O->begpos > end||O->begpos<beg||forcenext){
#else
    if(forcenext){
#endif

#endif
      free(O);
      forcenext=0;
      //      printf("Looking at lower-scoring overlaps\n");
      O=Find_Local_Overlap(alen,blen,
			 0 /*comp==0 -> fwd orientation */,
			 1 /*nextbest==1 -> next best */,
			 local_results,NumSegs,
			 minlen
#ifdef ENDHACK
			 -2*ENDGAPHACK
#endif
			 ,1.);
      if(O==NULL)goto nooverlap;

#if DEBUG_LOCALOVL > 3
    for(i=0;i<O->num_pieces;i++){
      printf("overlap segment %d: (%d,%d)[---------](%d,%d)\n",
	     i,
	     O->chain[i].piece.abpos,
	     O->chain[i].piece.bbpos,
	     O->chain[i].piece.aepos,
	     O->chain[i].piece.bepos
	     );
    }
#endif

    }
  }

#if DEBUG_LOCALOVL > 1
  if(O!=NULL){
    fprintf(stdout,"Overlap found.\n");
#if DEBUG_LOCALOVL > 3
    Print_Local_Overlap_withAlign(stdout,O,Ausable-MINUS_ONE,Busable-MINUS_ONE);
#else
    Print_Local_Overlap(stdout,O,5);
#endif
  }
#endif


  if(O==NULL){
    #if DEBUG_LOCALOVL > 0 
    fprintf(stdout,"No overlap out of Find_Local_Overlap()\n");
    #endif
    goto nooverlap;    
  }


  if(O->diffs/(double)O->length >= erate){ 
    /* global version fails */

    if(asymmetricEnds){

      /* Test for Mike's clear-range-extended overlaps */

      /* For the beginning of the overlap, the overlap should be rejected if:
	      - the first N bases of the B sequence are not aligned to A
	      - N > MaxBegGap (which should be set to the number of bases
	         the clear range was extended, plus a bit of slop)
	      - there is a "flap" (unaligned bases) at the beginning of
	        the alignment of sufficient size; this distinguishes cases
		which simply have a negative ahang from cases where there is
		a problematic mismatch */

      if(O->chain[0].bgap>MaxFreeFlap&&
	 O->chain[0].bgap+MAX(0,-O->begpos)>MaxBegGap){
#if DEBUG_LOCALOVL > 0 
	fprintf(stdout,"Overlap rejected because: begin flap too big %d\n",
		    O->chain[0].bgap+MAX(0,-O->begpos)>MaxBegGap);
#endif
	forcenext=1;
	goto next_overlap;
      }

      /* For the end of the overlap, as above but reverse orientation */
      
      if(O->chain[O->num_pieces].agap>MaxFreeFlap&&
	 O->chain[O->num_pieces].agap+MAX(0,-O->endpos)>MaxEndGap){

#if DEBUG_LOCALOVL > 0 
	fprintf(stdout,"Overlap rejected because: end flap too big %d\n",
		O->chain[O->num_pieces].agap+MAX(0,-O->endpos)>MaxEndGap);
#endif
	
	forcenext=1;
	goto next_overlap;
	
      }


    } else {  //  asymmetricEnds == 0


      if(min(O->chain[0].agap,
	     O->chain[0].bgap)>MaxBegGap){
#if DEBUG_LOCALOVL > 0 
	fprintf(stdout,"Overlap rejected because: begin gap too big %d\n",
		min(O->chain[0].agap,
		    O->chain[0].bgap));
#endif
	forcenext=1;
	goto next_overlap;
      }
      
      if(min(O->chain[O->num_pieces].agap,
	     O->chain[O->num_pieces].bgap)>MaxEndGap){
#if DEBUG_LOCALOVL > 0 
	fprintf(stdout,"Overlap rejected because: end gap too big %d\n",
		min(O->chain[O->num_pieces].agap,
		    O->chain[O->num_pieces].bgap));
#endif
	
	forcenext=1;
	goto next_overlap;
	
      }


    }  //  if asymmetricEnds


#ifndef CHECK_FINAL_QUALITY
    if(O->num_pieces-
	     (O->chain[0].type==LOCAL_BOUNDARY||
	      O->chain[0].agap+
	      O->chain[0].bgap==2)
	     +(O->chain[O->num_pieces].type!=LOCAL_BOUNDARY) 
       > MaxGaps   /* too many gaps */){
      #if DEBUG_LOCALOVL > 0 
      fprintf(stdout,"Overlap rejected because: too many pieces %d\n",
	      O->num_pieces-
	     (O->chain[0].type==LOCAL_BOUNDARY||
	      O->chain[0].agap+
	      O->chain[0].bgap==2)+
	      (O->chain[O->num_pieces].type!=LOCAL_BOUNDARY) );
      #endif

      forcenext=1;
      goto next_overlap;

    }

#endif

  }



  #undef KEEP_OLD_SUMSEGLEN_BUG
  #ifndef KEEP_OLD_SUMSEGLEN_BUG
  sumseglen=0;
  avgerror=0.;
  #endif



  for(i=0;i<O->num_pieces;i++){

    if(i!=0&&
       MAX(O->chain[i].agap,
	   O->chain[i].bgap)>MaxInteriorGap){
#if DEBUG_LOCALOVL > 0 
      fprintf(stdout,"Overlap rejected because: interior gap %d too big %d\n",
	      i,
	      MAX(O->chain[i].agap,
		  O->chain[i].bgap));
#endif

      forcenext=1;
      goto next_overlap;
      
    }

    /* segments assumed to be plus-sense wrt A sequence */
    assert(O->chain[i].piece.aepos-
	   O->chain[i].piece.abpos>0); 

    seglen=1+((O->chain[i].piece.aepos-
	    O->chain[i].piece.abpos
	    +abs(O->chain[i].piece.bepos-
		 O->chain[i].piece.bbpos))
           /2);

    if(seglen>=coreseglen)noCoreSeg=0;
    sumseglen+=seglen;
    avgerror+=O->chain[i].piece.error*seglen;
  }

  //Test for at least one segment bigger than minlen
  if(noCoreSeg){

#ifdef ENDHACK    
    /* if none found so far, might need to allow for a segment with short
       end gaps due to a mismatch (given +1/-3 scoring, 1 mismatch kills
       2 matches, for example) */
    int last=O->num_pieces-1;
    int begseglen = 1+
      (O->chain[i].piece.aepos-
       O->chain[i].piece.abpos
       +abs(O->chain[i].piece.bepos-
	    O->chain[i].piece.bbpos)
       )/2
      +min(ENDGAPHACK,O->chain[0].agap);
    int endseglen = 1+
      (O->chain[last].piece.aepos-
       O->chain[last].piece.abpos
       +abs(O->chain[last].piece.bepos-
	    O->chain[last].piece.bbpos)
       )/2
      +min(ENDGAPHACK,O->chain[last+1].agap);
    if(last==0)begseglen+=min(ENDGAPHACK,O->chain[1].agap);
    if(last==0)endseglen+=min(ENDGAPHACK,O->chain[0].agap);

    /* if even allowing some slop at ends we can't find a big enough segment, 
       quit! */
    if(MAX(begseglen,endseglen)<minlen){	      
#endif

      #if DEBUG_LOCALOVL > 0 
      fprintf(stdout,"Overlap rejected because no anchoring segment (min = %d)\n",
	      minlen);
      #endif
      
      forcenext=1;
      goto next_overlap;

#ifdef ENDHACK
    }
#endif
  }

#ifndef CHECK_FINAL_QUALITY
  avgerror+=O->num_pieces;
  sumseglen+=O->num_pieces;
  if(O->chain[0].agap<=1||
     O->chain[0].bgap<=1){
    avgerror-=1.;
    sumseglen-=1;
  }
  if(O->chain[O->num_pieces].agap>0||
     O->chain[O->num_pieces].bgap>0){
    avgerror+=1.;
    sumseglen+=1;
  }
  avgerror/=sumseglen;
  //Test for quality
  if(avgerror>erate){
    #if DEBUG_LOCALOVL > 0 
    fprintf(stdout,"Overlap rejected because: "
	    "affine-adjusted error rate %f > erate %f\n",avgerror,erate);
    #endif
    forcenext=1;
    goto next_overlap;
  } else {
    //    printf("avgerror came to %f\n",avgerror);
  }
#ifdef AFFINE_QUALITY
  O->length=(int)(sumseglen+.5);
  O->diffs=(int)(avgerror*sumseglen+.5);
#endif
#endif

  QVBuffer.quality=avgerror;

  /* Now that you have finish point, first compute complete alignment
     by calling AS_ALN_OKNAlign or extrapolate start point as dictated by
     `what'.  Then build an overlap record modelling the overlap.   */

  { int ahang, bhang;
    int *trace;

    //printf("Begpos for OverlapMesg == %d\n",O->begpos);

    
    {   int i;


	#if MINUS_ONE == 0 
	// coordinates from Find_Local routines will be one off from
	// those expected by the trace routines, so adjust them!
	//	O->begpos+= (O->begpos>=0 ? 1 : -1);
	//	O->endpos+= (O->endpos>=0 ? 1 : -1);
	for(i=0;i<=O->num_pieces;i++){
	  if(i<O->num_pieces){
	    O->chain[i].piece.abpos++;
	    O->chain[i].piece.bbpos++;
	    O->chain[i].piece.aepos++;
	    O->chain[i].piece.bepos++;
	  }
	}


	#endif
	//AS_Local_Trace assumes string pointer one before start of string!
	trace = AS_Local_Trace(O,Ausable-1,Busable-1);

#ifdef CHECK_FINAL_QUALITY
	{//TEST NUMBER OF GAPS AT SEGMENT BOUNDARIES
	  int final_num_pieces=0;
	  for(i=0;i<O->num_pieces;i++){
	    if(O->chain[i].piece.abpos < O->chain[i].piece.aepos) final_num_pieces++;
	  }
	  if(O->chain[0].bgap==0){final_num_pieces--;}
	  if(O->chain[O->num_pieces].agap>0){final_num_pieces++;}
	  if(final_num_pieces> MaxGaps){
	    #if DEBUG_LOCALOVL > 0 
	    fprintf(stdout,"Overlap rejected because: too many pieces %d\n",
		    O->num_pieces-
		    (O->chain[0].type==LOCAL_BOUNDARY||
		     O->chain[0].agap+
		     O->chain[0].bgap==2)+
		    (O->chain[O->num_pieces].type!=LOCAL_BOUNDARY) );
            #endif

	    forcenext=1;
	    goto next_overlap;
	  }
	}
#endif


#if DEBUG_LOCALOVL > 1
	fprintf(stdout,"Overlap after trimming:\n");
#if DEBUG_LOCALOVL > 3
	#if MINUS_ONE == 0 
	// coordinates from Find_Local routines will be one off from
	// those expected by the trace routines, so adjust them!
	//	O->begpos+= (O->begpos>=0 ? 1 : -1);
	//	O->endpos+= (O->endpos>=0 ? 1 : -1);
	for(i=0;i<=O->num_pieces;i++){
	  if(i<O->num_pieces){
	    O->chain[i].piece.abpos--;
	    O->chain[i].piece.bbpos--;
	    O->chain[i].piece.aepos--;
	    O->chain[i].piece.bepos--;
	  }
	}
	#endif
	for(i=0;i<O->num_pieces;i++){
	  printf("overlap segment %d: (%d,%d)[---------](%d,%d) gap(%d,%d)\n",
		 i,
		 O->chain[i].piece.abpos,
		 O->chain[i].piece.bbpos,
		 O->chain[i].piece.aepos,
		 O->chain[i].piece.bepos,
		 O->chain[i].agap,
		 O->chain[i].bgap
		 );
	}
	printf("overlap end gaps (%d,%d)\n",
	       O->chain[O->num_pieces].agap,
	       O->chain[O->num_pieces].bgap
	       );
	Print_Local_Overlap_withAlign(stdout,O,Ausable-MINUS_ONE,Busable-MINUS_ONE);
#else
	Print_Local_Overlap(stdout,O,5);
#endif


#endif

    }

    ahang=O->begpos; /* N.B.: can be changed in AS_Local_Trace*/
    //printf("ahang for OverlapMesg set to %d\n",ahang);

    bhang=O->endpos;

#define FIXENDGAPS
#ifdef FIXENDGAPS
    { 

      int i=0;
      int j=0;
      int changeahang=0;
      int changebhang=0;
      int fullLenA = strlen(a->sequence);
      int fullLenB = strlen(b->sequence);
#undef DEBUG_FINAL_TRACE
#ifdef DEBUG_FINAL_TRACE
      char c;
      printf("Trace (lens %d %d):",fullLenA,fullLenB);
#endif
      while(trace[i]!=0){
#ifdef DEBUG_FINAL_TRACE
       	c='*';
#endif

	if(trace[i]<-fullLenA){
	  changebhang++;
	} else if (trace[i]>fullLenB){
	  changebhang--;
	} else if (trace[i]==-1){
	  changeahang--;
	} else if (trace[i]==1){
	  changeahang++;
	} else {
#ifdef DEBUG_FINAL_TRACE
	  c=' ';
#endif
	  trace[j++]=trace[i];
	}
#ifdef DEBUG_FINAL_TRACE
	printf(" %c%d",c,trace[i]);
#endif
	i++;
      }
#ifdef DEBUG_FINAL_TRACE
      printf("\n");
#endif
      trace[j]=0;
      ahang+=changeahang;
      bhang+=changebhang;
    }
#endif  /* FIXENDGAPS */

#ifndef NONFINAL_BAND_TEST
    if( ahang > end|| ahang < beg ){
      forcenext=1;
      goto next_overlap;
    }
#endif

    if (ahang < 0 || (ahang == 0 && bhang > 0))
      { if (bhang >= 0)
          QVBuffer.overlap_type = AS_CONTAINMENT;
        else
          QVBuffer.overlap_type = AS_DOVETAIL;
        QVBuffer.ahg = -ahang;
        QVBuffer.bhg = -bhang;
        if (opposite)
          QVBuffer.orientation = AS_OUTTIE;
        else
          QVBuffer.orientation = AS_NORMAL;
        QVBuffer.aifrag = b->iaccession;
        QVBuffer.bifrag = a->iaccession;
        if (trace != NULL)
          { int i;
            for (i = 0; trace[i] != 0; i++)
              trace[i] = -trace[i];
          }
      }
    else
      { if (bhang <= 0)
          QVBuffer.overlap_type = AS_CONTAINMENT;
        else
          QVBuffer.overlap_type = AS_DOVETAIL;
        QVBuffer.ahg = ahang;
        QVBuffer.bhg = bhang;
        if (opposite)
          QVBuffer.orientation = AS_INNIE;
        else
          QVBuffer.orientation = AS_NORMAL;
        QVBuffer.aifrag = a->iaccession;
        QVBuffer.bifrag = b->iaccession;
      }
    assert(trace!=NULL);
    QVBuffer.delta = Pack_Alignment_AS(trace,QVBuffer.ahg);

    *where = ahang;

#undef FINALTRACEDEBUG
#ifdef FINALTRACEDEBUG
    { int i;
    for (i = 0; trace[i] != 0; i++)
      printf(" %d",trace[i]);
    printf(" trace elements %d\n",i);
    }
#endif



  }



#ifdef CHECK_FINAL_QUALITY
  { 
    int alen,blen,del,sub,ins,affdel,affins,blockdel,blockins;
    double errRate,errRateAffine;
    int AFFINEBLOCKSIZE=4;

    //////////////////
    // TEST FINAL FLAPS
    //////////////////

    if(asymmetricEnds){

      /* Test for Mike's clear-range-extended overlaps */

      /* For the beginning of the overlap, the overlap should be rejected if:
	      - the first N bases of the B sequence are not aligned to A
	      - N > MaxBegGap (which should be set to the number of bases
	         the clear range was extended, plus a bit of slop)
	      - there is a "flap" (unaligned bases) at the beginning of
	        the alignment of sufficient size; this distinguishes cases
		which simply have a negative ahang from cases where there is
		a problematic mismatch */

      // if there is a flap of sufficient size 
      if(MAX(O->chain[0].agap,O->chain[0].bgap)>MaxFreeFlap&&
	 // and the the unaligned amount of B is too long ...
	 O->chain[0].bgap+MAX(0,-O->begpos)>MaxBegGap){
#if DEBUG_LOCALOVL > 0 
	fprintf(stdout,"Overlap rejected because: begin flap too big %d\n",
		    O->chain[0].bgap+MAX(0,-O->begpos)>MaxBegGap);
#endif
	forcenext=1;
	goto next_overlap;
      }

      /* For the end of the overlap, as above but reverse orientation */
  
      // if there is a flap of sufficient size 
      if(MAX(O->chain[O->num_pieces].agap,O->chain[O->num_pieces].bgap)>MaxFreeFlap&&
	 // and the the unaligned amount of A is too long ...
	 O->chain[O->num_pieces].agap+MAX(0,-O->endpos)>MaxEndGap){

#if DEBUG_LOCALOVL > 0 
	fprintf(stdout,"Overlap rejected because: end flap too big %d\n",
		O->chain[O->num_pieces].agap+MAX(0,-O->endpos)>MaxEndGap);
#endif
	
	forcenext=1;
	goto next_overlap;
	
      }


    } else {  //  asymmetricEnds == 0


      if(min(O->chain[0].agap+MAX(O->begpos,0),
	     O->chain[0].bgap+MAX(-O->begpos,0))>MaxBegGap){
#if DEBUG_LOCALOVL > 0 
	fprintf(stdout,"Overlap rejected because: begin gap too big %d\n",
		min(O->chain[0].agap,
		    O->chain[0].bgap));
#endif
	forcenext=1;
	goto next_overlap;
      }
      
      if(min(O->chain[O->num_pieces].agap+MAX(-O->endpos,0),
	     O->chain[O->num_pieces].bgap+MAX(O->endpos,0))>MaxEndGap){
#if DEBUG_LOCALOVL > 0 
	fprintf(stdout,"Overlap rejected because: end gap too big %d\n",
		min(O->chain[O->num_pieces].agap,
		    O->chain[O->num_pieces].bgap));
#endif
	
	forcenext=1;
	goto next_overlap;
	
      }


    }  //  if asymmetricEnds


    //////////

    if (opposite){
      Complement_Fragment_AS(b);
    }

    Analyze_Affine_Overlap_AS(a,b,&QVBuffer,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
			      &affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE,&max_indel_AS_ALN_LOCOLAP_GLOBAL);

    //    {
    //      Print_Overlap_AS(stdout, a, b, &QVBuffer);
    //    }

    if (opposite){
      Complement_Fragment_AS(b);
    }

    errRate       = (sub+   ins+   del) / (double)(alen+      ins);
    errRateAffine = (sub+affins+affdel) / (double)(alen-del+affins+affdel);


    if(errRateAffine>erate){
      #if DEBUG_LOCALOVL > 0 

      fprintf(stdout,"Overlap rejected because: "
	      "affine-adjusted error rate %f > erate %f\n",errRateAffine,erate);
      #endif
      forcenext=1;
      goto next_overlap;
    }
    //    printf("avgerror came to %f\n",avgerror);
    #ifdef AFFINE_QUALITY
    QVBuffer.quality = errRateAffine;
    #endif

    //    printf("\n\nInner: Alen %d, Blen %d, del %d, sub %d, ins %d\n"
    //	   " affdel %d, affins %d, blockdel %d, blockins %d\n",
    //	   alen,blen,del,sub,ins,
    //	   affdel,affins,blockdel,blockins);
    //    printf("Simple mismatch rate %f\n",errRate);
    //    printf("Affine mismatch rate %f\n",errRateAffine);

  }
#else
    Analyze_Affine_Overlap_AS(a,b,&QVBuffer,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
			      &affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE, 0L);
#endif


  if (opposite){
    Complement_Fragment_AS(b);
  }
  if(O!=NULL)free(O);

  return (&QVBuffer);

nooverlap:
  if (opposite)
    Complement_Fragment_AS(b);
  if(O!=NULL)free(O);
  return ((OverlapMesg*)NULL);
  
}
