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

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)


/*  Handle Local_Overlap pieces (local segments) which overlap,
 *  by trimming them back until they abut, 
 *  such that the number of mismatches is minimized
 */
void fix_overlapping_pieces(const char *aseq,
                            const char *bseq,
                            Local_Overlap *O,
                            int piece0,
                            int piece1);




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


// set useSizeToOrderBlocks to 1 to get block mismatches resolved such that
//    the smaller block comes first:
//
//       .........AA------.........
//       .........--BBBBBB........
//
// or set it to 0 to get resolution with the A block always first
int useSizeToOrderBlocks = 1;



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

int *AS_Local_Trace(Local_Overlap *O, const char *aseq, const char *bseq){
  static int *TraceBuffer=NULL;
  int i,j,k,segdiff,*segtrace;
  int lastgood=-1;
  static int allocatedspace=0;
  int tracespace=0;
  static char *aseg=NULL,*bseg=NULL;

  //  Not computing traces generates slight differences.  Why?
  //
  const int computeTraceFlag = 0;

  static int aseglen=0,bseglen=0;
  int abeg=0,bbeg=0; /* begining of segment; overloaded */
  int tracep=0; /* index into TraceBuffer */
  int spnt=0; /* to pass to AS_ALN_OKNAlign */


  {
    int n=O->num_pieces;
    assert(O->num_pieces>0);
    O->chain[n].piece.abpos=O->chain[n].agap+O->chain[n-1].piece.aepos;
    O->chain[n].piece.bbpos=O->chain[n].bgap+O->chain[n-1].piece.bepos;
  }

  if(computeTraceFlag){
    /*Estimate length required to store trace*/
    tracespace=0;
    tracespace+=abs(O->begpos)+abs(O->endpos);
    for(i=0;i<=O->num_pieces;i++){
      tracespace+=max(O->chain[i].agap,1);
      tracespace+=max(O->chain[i].bgap,1);
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
  }  //  computeTraceFlag
  
  /* for each Local_Overlap chain[i].piece, 
     need to handle the gap at the beginning and
     (for all but the final piece) the associated aligned segment */
  for(i=0;i<=O->num_pieces;i++){

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

      if(O->chain[k].piece.abpos<O->chain[i].piece.aepos||
         O->chain[k].piece.bbpos<O->chain[i].piece.bepos){

	/* handle possibility of the first segment being
           contained within the second;originally simply asserted
           against this; now, try to handle by deleting first segment */

	if(O->chain[i].piece.abpos>O->chain[k].piece.abpos||
	   O->chain[i].piece.bbpos>O->chain[k].piece.bbpos){


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

        } else 	/* otherwise, check for 2nd piece contained within first */
	  if(O->chain[i].piece.aepos>O->chain[k].piece.aepos||
             O->chain[i].piece.bepos>O->chain[k].piece.bepos){

            /* if the next piece is completely within current piece, 
               effectively remove it */

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
	      

              O->chain[k+1].agap=O->chain[k+1].piece.abpos-
                O->chain[l].piece.aepos;
              O->chain[k+1].bgap=O->chain[k+1].piece.bbpos-
                O->chain[l].piece.bepos;
            }

            /* else, fix the overlap */
          } else {

            fix_overlapping_pieces(aseq,
                                   bseq,
                                   O,i,k);


            // if the second piece disappeared
            if(O->chain[k].piece.abpos==O->chain[k].piece.aepos||
               O->chain[k].piece.bbpos==O->chain[k].piece.bepos){

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
	      
                O->chain[k+1].agap=O->chain[k+1].piece.abpos-
                  O->chain[l].piece.aepos;
                O->chain[k+1].bgap=O->chain[k+1].piece.bbpos-
                  O->chain[l].piece.bepos;
              }
            } else {
              // if the first piece disappeared
              if (O->chain[i].piece.abpos==O->chain[i].piece.aepos||
                  O->chain[i].piece.bbpos==O->chain[i].piece.bepos){

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
              }
            }
          }
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
	O->chain[i].agap;
    }

    /*handle boundary case to prevent gaps preceding the b sequence*/
    if((i==0||lastgood<0)&&O->chain[i].bgap>0){
      assert(O->chain[i].agap>=0);

      if(O->begpos>=0){
        O->begpos=O->chain[i].piece.abpos-1;
        assert(O->begpos>=0);
        O->chain[i].agap=0;

        //  Instead of asserting, an ifdef previously printed stuff
        //  out and continued happily along.
        //
        assert( ( i==0&&	O->chain[i].bgap==O->chain[i].piece.bbpos-1)
                ||( i>0&&lastgood<0&&O->chain[i].bgap==O->chain[i].piece.bbpos-1)) ;

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
	O->endpos-=O->chain[i].agap;
	abeg-=O->chain[i].agap;
	O->chain[i].agap=0;
      }	
    }

    /* put gaps before the portion of bseq in the gap; for the first
       piece, this means before position 0 */
    if(i==0 || lastgood<0){
      bbeg = 1-min(O->begpos,0);
    } else {
      assert(lastgood<O->num_pieces);
      bbeg = O->chain[lastgood].piece.bepos;
    }

    /* now insert the right number of gaps! */

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


    if (computeTraceFlag) {
      if(O->chain[i].agap <= O->chain[i].bgap || ! useSizeToOrderBlocks ){
        /* start by putting len(agap) gaps before the chunk of B in the gap */
        for(j=0; j<O->chain[i].agap ;j++)
          TraceBuffer[tracep++]=bbeg;

        /* then put len(bgap) gaps before the chunk of A in the gap */
        for(j=0; j<O->chain[i].bgap ;j++)
          TraceBuffer[tracep++]=-abeg;
      } else { // if the bgap is smaller,
        abeg-=O->chain[i].agap;
        bbeg+=O->chain[i].bgap;

        /* start by putting len(bgap) gaps before the chunk of A in the gap */
        for(j=0;j<O->chain[i].bgap   ;j++)
          TraceBuffer[tracep++]=-abeg;

        /* then put len(agap) gaps before the chunk of B in the gap */
        for(j=0;j<O->chain[i].agap ;j++)
          TraceBuffer[tracep++]=bbeg;
      }
    } else {
      //  Not computing traces!
      if(O->chain[i].agap <= O->chain[i].bgap || ! useSizeToOrderBlocks ){
      } else {
        abeg-=O->chain[i].agap;
        bbeg+=O->chain[i].bgap;
      }
    }


    ///////////////////////////////////////  

    /* if last piece, there is no aligned segment */

    if(i==O->num_pieces)break;

    /* set bbeg to beginning of aligned segment for piece */

    abeg=O->chain[i].piece.abpos;
    bbeg=O->chain[i].piece.bbpos;

    /* set lengths of segments */
  
    int alen=O->chain[i].piece.aepos-abeg; /* check +1?? */
    int blen=O->chain[i].piece.bepos-bbeg; /* check +1?? */

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
  
    /* CMM we do not need the trace computed */
    if (computeTraceFlag) {

      /* copy the segments */

      strncpy(aseg, aseq+abeg, alen);
      strncpy(bseg, bseq+bbeg, blen);

      aseg[alen] = 0;
      bseg[blen] = 0;

      if (((int)strlen(bseg) != blen) ||
          ((int)strlen(aseg) != alen)) {
        fprintf(stderr,"EXCEPTION strlen(aseg)=%d alen=%d abeg=%d\n", strlen(aseg), alen, abeg);
        fprintf(stderr,"EXCEPTION strlen(bseg)=%d blen=%d bbeg=%d\n", strlen(bseg), blen, bbeg);

        fprintf(stderr,"EXCEPTION aseg=<%s>\n", aseg);
        fprintf(stderr,"EXCEPTION bseg=<%s>\n", bseg);

        fprintf(stderr,"EXCEPTION aseq=<%s>\n", aseq + 1);
        fprintf(stderr,"EXCEPTION bseq=<%s>\n", bseq + 1);

        return NULL; // Return an exceptional value.
      }

      /* guesstimate the required number of diagonals/edits to consider to
         get optimal alignment */
      segdiff = 1 + (int)((O->chain[i].piece.aepos - O->chain[i].piece.abpos) * 1.5 * O->chain[i].piece.error);


      /* get trace for the segment from AS_ALN_OKNAlign */
      spnt=0; 
      /* subtract from aseg, bseg because Gene likes to index from 1, not 0 */
      segtrace=AS_ALN_OKNAlign(aseg-1,alen,bseg-1,blen,&spnt,segdiff);
      // This adjusts the beginning coordinates so that segment is
      // consistent with the back-trace.

      if(spnt>0){
        O->chain[i].agap+=spnt;
        O->chain[i].piece.abpos+=spnt;
      } else {
        O->chain[i].bgap-=spnt;
        O->chain[i].piece.bbpos-=spnt;
      }

      /* get trace for the segment from AS_ALN_OKNAffine */
      /* Seems like it should be a good idea, but doesn't work as well
         as we might expect! */
      //bpnt=0; 
      //epnt=0; 
      //segtrace=AS_ALN_OKNAffine(aseg,alen,bseg,blen,&bpnt,&epnt,segdiff);

      assert(segtrace!=NULL);

      /* Now copy the segment trace into master trace, adjusting positions */
      j=0;

      if(spnt<0){
        for(int ctr=0;ctr<abs(spnt);ctr++)
          TraceBuffer[tracep++]=-abeg;
      } else {
        for(int ctr=0;ctr<spnt;ctr++)
          TraceBuffer[tracep++]=bbeg;
      }

      while(segtrace[j]!=0){
        if(segtrace[j]<0){
          TraceBuffer[tracep++]=-abeg+segtrace[j++]+1 /* -max(0,spnt) ?? */;
        } else {
          TraceBuffer[tracep++]=bbeg+segtrace[j++]-1 /* +max(0,-spnt) ?? */;
        }
      }
    }  //  computeTraceFlag
    
    /* set lastgood to this segment */

    lastgood=i;

    /* and back to the top of the loop for another overlap piece */
  }  

  /* terminate the trace */
  if(TraceBuffer != NULL) {
    TraceBuffer[tracep]=0;

    if (tracep >= allocatedspace)
      fprintf(stderr,"ERROR memory is already corrupted in %s at %d.\n", __FILE__, __LINE__);
    assert(tracep < allocatedspace);
  }

  return(TraceBuffer);
}
