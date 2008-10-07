
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "CA_ALN_local.h"
#include "AS_ALN_aligners.h"
#include "CA_ALN_scafcomp.h"

#define OVL_ERATE   .99
#define MIN_SEG_LEN  40
#define SEG_ERATE   .10

/* Debug conditional compilation flags */

#undef  DEBUG_COMPARE
#undef  DEBUG_RAYSHOOT
#undef  DEBUG_ALIGN
#undef  DEBUG_ANALYSIS
#undef  DEBUG_LOCAL
#undef  DEBUG_OPTTAIL
#undef  DEBUG_SEGORDER

#define ALLOW_NO_OVERLAP_INTERLEAVINGS



typedef struct CO_tag {
  Segment       *seg;
  int            best;
  struct CO_tag *trace;
  struct CO_tag *Alink;
  struct CO_tag *Blink;
} COvlps;


typedef struct _ScafLap_tag {
  struct _ScafLap_tag *next;    /* NULL-terminated linked-list of said */
  Scaffold            *ascaf;   /* A scaffold */
  Scaffold            *bscaf;   /* B scaffold */
  int                  asnum;   /* A scaffold number in assembly */
  int                  bsnum;   /* B scaffold number in assembly */
  int                  abase;   /* #(in assembly) of 1st contig in A scaffold */
  int                  bbase;   /* #(in assembly) of 1st contig in B scaffold */
  Overlap             *overlap; /* Overlap between scaffolds (sequence and
                                   trace fields are NULL!)              */
  int                  score;   /* Sum of lengths contig overlaps involved
                                   in this scaffold overlap.         */
  float                erate;   /* Error rate of matching parts. */
  int                  firm;    /* Alignment satisfies scaffold constraints */
  int                  D_delta; /* Width of band containing all contig
                                   overlaps of this scaffold overlap */
  int                  D_var;   /* Sum of all gap variation between
                                   contigs involved in the overlap */
  int                  A_delta; /* Sum of all A-scaffold segments that should
                                   have been but were not aligned     */
  int                  B_delta; /* Sum of all B-scaffold segments that should
                                   have been but were not aligned     */
  Segment             *seglist; /* List of contig overlaps constituting the
                                   scaffold overlap                      */
} Scaffold_Overlap;


int     MaxAlign  = -1;
int     MaxBucket = -1;
COvlps  *CtgOvls  = NULL;
COvlps **ABuckets = NULL;
COvlps **BBuckets = NULL;




typedef struct ival_tag {
  double beg;
  double end;
  COvlps *traceback;
} interval;

typedef struct ilist_tag {
  interval ival;
  struct ilist_tag *prev;
  struct ilist_tag *next;
}
interval_list;


static
interval_list *
cleanup_ilist(interval_list *list){
  interval_list *link;
  int fromTail=0;

  if(list==NULL)return NULL;

  // we should be on one end or the other
  assert( list->next==NULL || list->prev==NULL);

  // if next is null, we are at end (else beginning)
  fromTail = ( list->next==NULL ? 1 : 0 );

  // trace through the list, deleting as we go
  while(list!=NULL){
    link = (fromTail ? list->prev : list->next);
    safe_free(list);
    list=link;
  }
  return NULL;
}

static
interval_list *
add_to_ilist(interval_list *tail, interval to_add);

static
interval_list *
add_to_ilist_special(interval_list *tail, interval to_add){
  interval_list *curr,*next=NULL,*to_delete_tail=NULL;

  // we should only be here if certain conditions are met:
  assert(tail!=NULL);
  assert(to_add.beg < tail->ival.beg);

  fprintf(stderr,"Inefficient addition into interval list!\n");

  curr=tail;

  // back up until to_add would be a simple append
  while(curr!=NULL&&curr->ival.beg>to_add.beg){
    next=curr;
    curr=curr->prev;
  }
  assert(next!=NULL);

  // undo backwards link from next, so we can do a simple cleanup later
  next->prev=NULL;
  if (curr)
    curr->next=NULL;

  if (curr)
    curr->next=NULL;

  // for all elements in chain headed by next, append them serially,
  // unless the rest of the chain can simply be tacked on,
  // in which case we need to prepare for a special cleanup
  curr=add_to_ilist(curr,to_add);
  while(next!=NULL){
    if(next->ival.beg >= curr->ival.end){
      curr->next=next;
      to_delete_tail=next->prev;
      if(to_delete_tail!=NULL){
	to_delete_tail->next=NULL;
      }
      next->prev=curr;
      break;
    }
    to_delete_tail=next;
    curr=add_to_ilist(curr,next->ival);
    next=next->next;
  }
  cleanup_ilist(to_delete_tail);

  while(curr->next!=NULL){
    curr=curr->next;
  }

  fprintf(stderr,"Inefficient addition into interval list - finished.\n");

  assert(curr->prev!=curr);
  return(curr);

}

static
interval_list *
add_to_ilist(interval_list *tail, interval to_add){
  interval_list *il;

  assert(to_add.beg<=to_add.end);

  // if this interval is being added to an existing (nonempty) list,
  if (tail!=NULL){

    // mostly, segments will be added in an order, and constrained in a fashion,
    // that allow us to consider only the last segment and the new one to check for overlaps.
    //
    // Ideally, this would never happen and we could:
    //      assert( tail->ival.beg <= to_add.beg);
    //
    // But so far that is not strictly maintained.  Additions violating these criteria can be handled, but it is clearly less efficient,
    // especially if applied to a simple linked list rather than something that can be searched in O(log N) rather than O(N).
    // Currently, we use the linked list, so let us make it a special case:

    if(tail->ival.beg > to_add.beg){
      return(add_to_ilist_special(tail,to_add));
    }

    // if the intervals exactly match ...
    if(tail->ival.beg == to_add.beg && tail->ival.end == to_add.end){

      // nothing will be appended: we either discard new or replace tail

      // if tail is at least as good,
      if( to_add.traceback==NULL || (tail->ival.traceback != NULL && tail->ival.traceback->best >= to_add.traceback->best ) ){

	//do nothing
	assert(tail->ival.beg<=tail->ival.end);
	assert(tail->prev!=tail);
	return tail;

      } else { // new is better

	// replace tail entirely, in place
	tail->ival.traceback=to_add.traceback;
	assert(tail->ival.beg<=tail->ival.end);
	assert(tail->prev!=tail);
	return tail;

      }

    }

    // special case when new segment is contained within old
    if(tail->ival.end>to_add.end){
      if(to_add.traceback==NULL|| (tail->ival.traceback != NULL && tail->ival.traceback->best >= to_add.traceback->best ) ){
	assert(tail->prev!=tail);
	return(tail);
      } else {
	interval bonus=tail->ival;
	bonus.beg=to_add.end;
	tail->ival.end=to_add.beg;
	if(tail->ival.beg==tail->ival.end){
	  tail->ival=to_add;
	} else {
	  tail = add_to_ilist(tail,to_add);
	}
	if(bonus.beg<bonus.end){
	  tail = add_to_ilist(tail,bonus);
	}
      }
      assert(tail->prev!=tail);
      return(tail);
    }

    // if there is either a propoer overlap or tail is a subset of to_add
    if(tail->ival.end > to_add.beg){

      // if tail is a subset of new
      if( tail->ival.beg == to_add.beg && tail->ival.end < to_add.end){

	// if tail is better, scorewise
	if(to_add.traceback == NULL || (tail->ival.traceback != NULL && tail->ival.traceback->best > to_add.traceback->best)){

	  // adjust new to append just part that sticks out past tail
	  to_add.beg=tail->ival.end;
	  assert(to_add.beg<=to_add.end);
	} else {

	  // replace tail in place
	  tail->ival.end=to_add.end;
	  tail->ival.traceback=to_add.traceback;
	  assert(tail->ival.beg<=tail->ival.end);
	  assert(tail->prev!=tail);
	  return tail;

	}

      } else { // tail and new have a proper overlap

	// if tail is better, scorewise
	if(to_add.traceback == NULL || ( tail->ival.traceback != NULL && tail->ival.traceback->best > to_add.traceback->best)){

	  // adjust new
	  to_add.beg=tail->ival.end;
	  assert(to_add.beg<=to_add.end);
	} else {

	  // adjust tail
	  tail->ival.end=to_add.beg;
	  assert(tail->ival.beg<=tail->ival.end);

	}

      }
    }
  }

  il = (interval_list*) safe_malloc(sizeof(interval_list));
  il->ival.beg=to_add.beg;
  il->ival.end=to_add.end;
  il->ival.traceback=to_add.traceback;
  il->prev=tail;
  il->next=NULL;

  if(tail!=NULL){
    tail->next=il;
  }

  assert(il->prev!=il);
  return il;
}


// Project_across_Agap_one_interval():
// Project a single interval accessible on left side of A-gap across the A-gap
//   need input accessibility interval
//        A gap variance
//        B gap info
//        number of sigma allowed
//   modify accessibility interval to provide corresponding interval on the right side of the A-gap
//   return whether we came out the bottom -- i.e. found a solution
//
// Important note: if an interval starts (or ends) in a B-gap, the interpretation is that the amount of stretch
// in gap size above or below the point is only the relevant fraction of the gap's stretchiness.
//

static
int
Project_across_Agap_one_interval(interval *inoutIval,COvlps **bestTerm, double Agap_length, double Agap_var,Scaffold *B, int varwin){
  double top=inoutIval->beg;
  double bot=inoutIval->end;
  int terminal=0;

  {
    // find highest location across the gap
    // This happens when
    //   - Agap is smallest
    //   - Bgaps are biggest
    double remainingToUseUp;
    double BctgBeg,BctgEnd;
    int i,j;
    double Bmin;
    double fraction_of_gap;

    Bmin=top;

    // figure out the minimum distance to be accounted for by the projection


    // what SHOULD we do if minimum length of gap is negative?

#undef ALLOW_NEG_GAP_BACKUP
#ifdef ALLOW_NEG_GAP_BACKUP
    remainingToUseUp = Agap_length -  Agap_var * varwin;
    if(remainingToUseUp < 0){

      { // if Bmin falls in a gap, figure out how to adjust for fact that a
	// position in the gap factors in the stretchiness of the gap itself ...
	j=0;
	double ctgend;
	while(j<B->num_gaps && B->ctgs[j+1].lft_end<=Bmin){
	  j++;
	}
	ctgend=B->ctgs[j].lft_end+B->ctgs[j].length;
	if(ctgend<Bmin){
	  double gapFrac,adjustment;
	  assert(j<B->num_gaps);
	  gapFrac = (Bmin-ctgend)/B->gaps[j].gap_length;
	  assert(gapFrac>=0);
	  adjustment = gapFrac*(B->gaps[j].gap_length+B->gaps[j].gap_var*varwin);
	  if(adjustment <= -remainingToUseUp){
	    remainingToUseUp += adjustment;
	    Bmin = ctgend;
	  } else {
	    gapFrac = -remainingToUseUp/(B->gaps[j].gap_length+B->gaps[j].gap_var*varwin);
	    assert(gapFrac>=0);
	    remainingToUseUp=0;
	    Bmin -= gapFrac * B->gaps[j].gap_length;
	    assert(Bmin>ctgend);
	  }
	}
      }
      assert(remainingToUseUp<=0);

      Bmin+=remainingToUseUp;
      if(Bmin<0){
	Bmin=0;
      }
      top=Bmin;
      remainingToUseUp=0;
    }
#else
    remainingToUseUp = Agap_length -  Agap_var * varwin;
    if (remainingToUseUp < 0) {
      fprintf(stderr, "ALLOW_NEG_GAP_BACKUP would have fixed a negative remainingToUseUp of %f; we just set it to zero now!\n", remainingToUseUp);
      remainingToUseUp = 0;
    }
#endif

    // once past contigs and gaps clearly above the interval,
    // contigs or fractions thereof use up their own length
    // B gaps use up their mean PLUS maximum stretch
    //
    // So, process a contig+gap at a time

    i=0;
    while ( remainingToUseUp > 0 && i < B->num_gaps){
      double ctgend;
      // skip contigs and gaps completely above the relevant interval
      if(B->ctgs[i+1].lft_end<=Bmin){
	i++;
	continue;
      }

      ctgend=B->ctgs[i].lft_end+B->ctgs[i].length;

      // handle contig portion
      if(Bmin<ctgend){
	remainingToUseUp -= ctgend-Bmin;
	Bmin=ctgend;

	// if we end within the contig, we are done
	if(remainingToUseUp<=0){
	  Bmin+=remainingToUseUp;
	  remainingToUseUp=0;
	  break;
	}
      }

      // now use the gap

      // if not full gap, figure out what fraction
      fraction_of_gap = 1.;
      if(Bmin!=ctgend){
	fraction_of_gap -= ((double)( Bmin-ctgend)) /(double) B->gaps[i].gap_length;
	assert(fraction_of_gap>=0);
      }
      remainingToUseUp -=  fraction_of_gap * (B->gaps[i].gap_length + B->gaps[i].gap_var*varwin);
      Bmin = B->ctgs[i+1].lft_end;
      if(remainingToUseUp < 0 ) {
	if(0){
	  Bmin = remainingToUseUp;
	  if(Bmin < ctgend){
	    Bmin = ctgend;
	  }
	  remainingToUseUp=0;
	} else {
	  double frac = -((double) remainingToUseUp)/(double)(B->gaps[i].gap_length + B->gaps[i].gap_var*varwin) ;
	  assert(frac>=0);
	  Bmin = B->ctgs[i+1].lft_end- (frac*B->gaps[i].gap_length);
	  remainingToUseUp=0;
	}
	break;
      }
      i++;
    }
    assert( remainingToUseUp == 0 || i == B->num_gaps );

    if(remainingToUseUp>0){

      //handle the final contig

      assert(Bmin==B->ctgs[i].lft_end||(Bmin==top&&top>B->ctgs[i].lft_end));

      // is this conditional gratuitous?
      if(Bmin<B->ctgs[i].lft_end+B->ctgs[i].length){

	remainingToUseUp -= B->ctgs[i].lft_end+B->ctgs[i].length-Bmin;
	Bmin=B->ctgs[i].lft_end+B->ctgs[i].length;

	// if we end within the contig, we are done
	if(remainingToUseUp<=0){
	  Bmin+=remainingToUseUp;
	  remainingToUseUp=0;
	} else { // we cannot help but come out the bottom

	  if(*bestTerm==NULL || ( inoutIval->traceback!=NULL && (*bestTerm)->best < inoutIval->traceback->best) ){
	    *bestTerm = inoutIval->traceback;
	  }
	  terminal=1;
	}

      }

    }
    top = Bmin;
  }

  {
    // find lowest location across the gap
    // This happens when
    //   - Agap is largest
    //   - Bgaps are smallest
    double remainingToUseUp;
    double BctgBeg,BctgEnd;
    int i,j;
    double Bmax;

    Bmax=bot;

    // figure out the maximum distance to be accounted for by the projection
    remainingToUseUp = MAX(0, Agap_length + Agap_var * varwin);
    // what SHOULD we do if max length of gap is negative?
    assert(remainingToUseUp>=0);

    // once past contigs and gaps clearly above the interval,
    // contigs or fractions thereof use up their own length
    // B gaps use up their mean MINUS maximum stretch
    //
    // So, process a contig+gap at a time

    i=0;
    while ( remainingToUseUp > 0 && i < B->num_gaps){
      double ctgend;
      double fraction_of_gap;

      // skip contigs and gaps completely above the relevant interval
      if(B->ctgs[i+1].lft_end<=Bmax){
	i++;
	continue;
      }

      ctgend=B->ctgs[i].lft_end+B->ctgs[i].length      ;

      // handle contig portion
      if(Bmax<ctgend){
	remainingToUseUp -= ctgend-Bmax;
	Bmax=ctgend;

	// if we end within the contig, we are done
	if(remainingToUseUp<=0){
	  Bmax+=remainingToUseUp;
	  remainingToUseUp=0;
	  break;
	}
      }

      // now use the gap

      // if not full gap, figure out what fraction
      if(Bmax!=ctgend){
	fraction_of_gap = 1.-((double)( Bmax-ctgend)) /(double) B->gaps[i].gap_length;
	assert(fraction_of_gap>=0 && fraction_of_gap <= 1);
      } else {
	fraction_of_gap = 1.;
      }

      remainingToUseUp -=  fraction_of_gap * (B->gaps[i].gap_length - B->gaps[i].gap_var*varwin);

      Bmax = B->ctgs[i+1].lft_end;
      if(remainingToUseUp < 0 ) {
	assert(B->gaps[i].gap_length - B->gaps[i].gap_var*varwin >= 0) ;
	if(0){
	  Bmax -= remainingToUseUp;
	  assert(Bmax >= ctgend);
	  remainingToUseUp=0;
	} else {
	  double frac = -((double) remainingToUseUp)/(double)(B->gaps[i].gap_length - B->gaps[i].gap_var*varwin) ;
	  assert(frac>=0);
	  Bmax = B->ctgs[i+1].lft_end - (frac*B->gaps[i].gap_length);
	  remainingToUseUp=0;
	}
	break;
      }
      i++;
    }
    assert( remainingToUseUp == 0 || i == B->num_gaps );

    if(remainingToUseUp>0){

      //handle the final contig

      assert(Bmax==B->ctgs[i].lft_end||(Bmax==bot&&bot>B->ctgs[i].lft_end));

      // is this conditional gratuitous?
      if(Bmax<B->ctgs[i].lft_end+B->ctgs[i].length){

	remainingToUseUp -= (B->ctgs[i].lft_end + B->ctgs[i].length - Bmax);
	Bmax=B->ctgs[i].lft_end+B->ctgs[i].length;

	// if we end within the contig, we are done
	if(remainingToUseUp<=0){
	  Bmax+=remainingToUseUp;
	  remainingToUseUp=0;
	} else { // we cannot help but come out the bottom
	  terminal=1;
	}
      }

    }

    bot = Bmax;
  }

  inoutIval->beg=top;
  inoutIval->end=bot;
  assert(top<=bot);
  return terminal;
}

// Project_across_Agap():
// Project a set of intervals accessible on left side of A-gap across the A-gap
//   need input list of accessibility intervals
//        A gap variance
//        B gap info
//        number of sigma allowed
//   construct accessibility interval list to provide corresponding interval on the right side of the A-gap
//   return whether we came out the bottom -- i.e. found a solution
//
// To process the left edge of an A-gap (list of accessible intervals)
// For each vertical interval,
//   If terminal "gap"
//     If best so far, update best
//   Else
//     project across A-gap
//     If go off bottom,
//       If best so far, update best

static
int
Project_across_Agap(interval_list *inList,interval_list **outList, COvlps **bestTerm,double Agap_length, double Agap_var,Scaffold *B, int varwin){
  interval_list *curr=inList;
  int terminal=0;

  while(curr!=NULL){
    int thisterm = Project_across_Agap_one_interval(&(curr->ival),bestTerm,Agap_length,Agap_var,B,varwin);
    if(thisterm){
      terminal=1;
      if(*bestTerm==NULL || ( curr->ival.traceback!=NULL && curr->ival.traceback->best > (*bestTerm)->best) ){
	*bestTerm=curr->ival.traceback;
      }
    }
    *outList=add_to_ilist(*outList,curr->ival);
    curr=curr->next;
  }
  while(*outList!=NULL&&(*outList)->prev!=NULL)*outList=(*outList)->prev;
  return terminal;
}


static
int
Process_seg(COvlps *af,COvlps **bestTerm,
            interval_list **nextAgapList,
            interval_list **BgapLists,
            int whichA,int whichB, Scaffold *A, Scaffold*B,int varwin){
  interval to_add;
  int terminal=0;

  assert(af->seg->a_contig==whichA);
  assert(af->seg->b_contig==whichB);

  to_add.traceback=af;

  if(af->seg->overlap->endpos>0){ // we come out in middle of B contig ... makes portion of right  edge of B contig 'accessible'


    // coordinates are global, relative to entire B scaffold

    to_add.beg = B->ctgs[whichB].lft_end + (B->ctgs[whichB].length-af->seg->overlap->endpos);
    assert(to_add.beg>=B->ctgs[whichB].lft_end);
    to_add.end=to_add.beg;

    *nextAgapList=add_to_ilist(*nextAgapList,to_add);

    // if terminal, we should deal with it at a higher level

  }else { // we come out in middle (or at end) of A contig ... makes portion of bottom edge of contig  accessible

    // contig coordinates are local (i.e. relative to A contig)
    // endpos is negative
    to_add.beg=A->ctgs[whichA].length + af->seg->overlap->endpos;
    assert(to_add.beg >= 0);
    to_add.end=to_add.beg;

    BgapLists[whichB]=add_to_ilist(BgapLists[whichB],to_add);

    // if terminal, process accordingly
    if(whichB==B->num_gaps){
      terminal=1;
      if( *bestTerm==NULL || (*bestTerm)->best < af->best ){
	*bestTerm = af;
      }
    }

  }

  return(terminal);
}


static
int
Process_Agap_one_accessible_interval(interval curr,COvlps **bestTerm,
                                     interval_list ** nextAgapList,
                                     interval_list ** BgapLists,
                                     int whichA,Scaffold *A, Scaffold *B,int varwin){
  int i=0;
  double top=curr.beg;
  double bot=curr.end;
  int terminal=0;
  while(i<B->num_gaps && bot>B->ctgs[i].lft_end){
    double ctgend;
    // skip contigs and gaps completely above the relevant interval
    if(B->ctgs[i+1].lft_end <= top){
      i++;
      continue;
    }
    ctgend=B->ctgs[i].lft_end+B->ctgs[i].length;

    // process B contig, if overlapped

    // To process the top edge of a B contig (in the context of a particular A contig)
    // For each accessible (horizontal) interval
    //   For each AxB overlap segment
    //     If the segment starts in the interval,
    //       If the segment ends on the A-gap (comes out in middle of B contig), add to following A-gap's left-edge accessible list
    //       Else, add exit (in middle of A contig) to bottom edge of B contig's accessible list

    if(top<ctgend){
      double from=MAX(top,B->ctgs[i].lft_end);
      double to=MIN(bot,ctgend);
      assert(from >= B->ctgs[i].lft_end);

      from-=B->ctgs[i].lft_end;
      to-=B->ctgs[i].lft_end;

      // {from, to} are positions along the contig, in local (contig) coordinates
      {
	COvlps *afing,*af;
	afing=ABuckets[whichA];
	// while there are segments on far side of horizontal gap, advance until we are in or past the relevant b-contig
	while (afing != NULL && afing->seg->b_contig < i)
	  afing = afing->Alink;
	// over all segments in contigs [whichA,i]
	for (af = afing; af != NULL && af->seg->b_contig == i; af = af->Alink){
	  // if the segment begins in the accessible interval
	  if(af->seg->overlap->begpos <= 0 ) { //ahang is nonpositive, i.e. segment starts in B contig
	    double pnt = - af->seg->overlap->begpos;
	    if (from <= pnt && pnt <= to) {
	      int score = af->seg->overlap->length - MAX(0,-af->seg->overlap->begpos) - MAX(0,-af->seg->overlap->endpos);
	      assert(score>0);
	      if (curr.traceback != NULL)
		score += curr.traceback->best;
	      if (score > af->best){
		af->best  = score;
		af->trace = curr.traceback;
	      }
	      // now process the segment's projections
	      if(Process_seg(af,bestTerm,nextAgapList,BgapLists,whichA,i,A,B,varwin)){
		terminal=1;
	      }
	    }
	  }
	}
      }
    }

    // process following gap, if overlapped
    {
      double from, to;
      from=MAX(top,ctgend);
      to=MIN(bot,ctgend+B->gaps[i].gap_length);

      // if there is an interval...
      if(from<to){

	// determine portion of A contig that can be reached
	double left, right, amountOff_A_contig=0;

	// leftmost point along A contig is:
	//     [ fraction of B gap that is below 'to' value ] * [ maximum height of B gap ]
	//N.B. Things are screwy here if the gap can be negative
	{
	  double frac;
	  frac= ((double)(B->ctgs[i+1].lft_end-to))/(double)B->gaps[i].gap_length;
	  assert(frac>=0);
	  left = MAX(0, frac * ((double) B->gaps[i].gap_length - B->gaps[i].gap_var*varwin) );
	}

	// rightmost point is:
	//      [ fraction of B gap that is below 'from' value ] * [ minimum height of B gap ]
	right = (((double)(B->ctgs[i+1].lft_end-from))/(double) B->gaps[i].gap_length * (B->gaps[i].gap_length + B->gaps[i].gap_var*varwin) );
	assert(right>=left);
	// ... though if this is off the end of the contig, we must correct
	if( right > A->ctgs[whichA].length){
	  amountOff_A_contig = right - A->ctgs[whichA].length;
	  right = A->ctgs[whichA].length;
	}

	// now find and process any segments starting in this interval ...

	if(right>=left){ // this will be false if 'left' is off right end of contig ... we have not adjusted this away
	  COvlps *afing,*af;
	  afing=ABuckets[whichA];
	  // while there are segments on far side of horizontal gap, advance until we are in or past the relevant b-contig
	  while (afing != NULL && afing->seg->b_contig <= i)
            afing = afing->Alink;
	  // over all segments in contigs [whichA,i]
	  for (af = afing; af != NULL && af->seg->b_contig == i+1; af = af->Alink){
	    // if the segment comes begins in the accessible interval
	    if(af->seg->overlap->begpos>0){
	      if (left <= af->seg->overlap->begpos && af->seg->overlap->begpos <= right) {
		int score = af->seg->overlap->length - MAX(0,-af->seg->overlap->begpos) - MAX(0,-af->seg->overlap->endpos);
		assert(score>0);
		if (curr.traceback != NULL)
		  score += curr.traceback->best;
		if (score > af->best){
		  af->best  = score;
		  af->trace = curr.traceback;
		}
		// now process the segment's projections
		if(Process_seg(af,bestTerm,nextAgapList,BgapLists,whichA,i+1,A,B,varwin)){
		  terminal=1;
		}
	      }
	    }
	  }
	}

	// now, if there is any projection onto the far A gap, add interval to the A gap
	if( amountOff_A_contig > 0){
	  double fracTop, fracBot;
	  double fartop,farbot;
	  interval ival;

	  fracTop = ((double)amountOff_A_contig) /(double) (B->gaps[i].gap_length + B->gaps[i].gap_var*varwin);

	  // things get weird when the gap can be negative,
	  // so let us sanity check:
	  // it seems that if left is past end of contig,
	  // (recalling that left is set based on smallest possible
	  // gap length), this should only happen when smallest
	  // possible gap length is positive
	  assert( ( left - A->ctgs[whichA].length <= 0 ) || B->gaps[i].gap_length - B->gaps[i].gap_var*varwin > 0 );


	  fracBot = ((double) MAX(0,left - A->ctgs[whichA].length )) / (B->gaps[i].gap_length - B->gaps[i].gap_var*varwin);

	  assert(fracBot>=0&&fracBot<=fracTop&&fracTop<=1);

	  ival.beg=B->ctgs[i+1].lft_end-(fracTop*B->gaps[i].gap_length);
	  ival.end=B->ctgs[i+1].lft_end-(fracBot*B->gaps[i].gap_length);
	  ival.traceback=curr.traceback;

	  if(ival.end<B->ctgs[i+1].lft_end){
	    assert( B->gaps[i].gap_length - B->gaps[i].gap_var*varwin > 0 );
	  }

	  *nextAgapList = add_to_ilist(*nextAgapList,ival);
	}

      }
    }
    i++;
  }

  // process final B contig, if overlapped, and test for terminal
  if( bot > B->ctgs[i].lft_end ){
    assert(i==B->num_gaps);
    double ctgend=B->ctgs[i].lft_end+B->ctgs[i].length;
    // if there is an overlap
    if(top<ctgend){
      double from=MAX(top,B->ctgs[i].lft_end);
      double to=MIN(bot,ctgend);
      assert(from<=to);
      to-=B->ctgs[i].lft_end;
      from-=B->ctgs[i].lft_end;
      // {from, to} are positions along the contig, in local (contig) coordinates
      {
	COvlps *afing,*af;
	afing=ABuckets[whichA];
	// while there are segments on far side of horizontal gap, advance until we are in or past the relevant b-contig
	while (afing != NULL && afing->seg->b_contig < i)
	  afing = afing->Alink;
	// over all segments in contigs [whichA,i]
	for (af = afing; af != NULL && af->seg->b_contig == i; af = af->Alink){
	  // if the segment begins in the accessible interval
	  if(af->seg->overlap->begpos <= 0 ) { //ahang is nonpositive, i.e. segment starts in B contig
	    double pnt = -af->seg->overlap->begpos;
	    if (from <= pnt && pnt <= to) {
	      int score = af->seg->overlap->length - MAX(0,-af->seg->overlap->begpos) - MAX(0,-af->seg->overlap->endpos);
	      assert(score>0);
	      if (curr.traceback != NULL)
		score += curr.traceback->best;
	      if (score > af->best){
		af->best  = score;
		af->trace = curr.traceback;
	      }
	      // now process the segment's projections
	      if(Process_seg(af,bestTerm,nextAgapList,BgapLists,whichA,i,A,B,varwin)){
		terminal=1;
	      }
	    }
	  }
	}
      }
    }
    // if we go off the bottom ...
    if(bot>ctgend){
      terminal=1;
      if(*bestTerm==NULL || ( curr.traceback!=NULL && (*bestTerm)->best < curr.traceback->best) ){
	*bestTerm = curr.traceback;
      }
    }
  }
  return terminal;
}


static
int
Process_Agap_accessible_intervals(interval_list **accessList,
                                  interval *contigTopEdgeAccess,
                                  interval_list **nextAgapList,
                                  interval_list **BgapLists,
                                  COvlps **bestTerm,
                                  int whichA, Scaffold *A, Scaffold *B, int varwin){


  // To process the right edge of an A-gap ( list of accessible intervals )
  // For each (vertical) interval,
  //   For each B contig it overlaps
  //     If any AxB segments start in middle of the B contig
  //       If exit point is into A-gap
  //         Add exit point to A-gap's left-edge accessible intervals
  //       If exit point is into B-gap, i.e. out of middle of A contig
  //         Add exit point to B contig's bottom edge's accessible intervals
  //   For each B gap it overlaps,
  //     Add interval at left edge of following A-gap that can be reached to that A-gap's accessible intervals
  //     Add interval on A contig that is accessible top edge of following B contig to that
  //   If interval falls off the bottom, declare success

  interval_list *curr=NULL;
  int terminal =0;

  // handle any accessible segments along top edge
  if(contigTopEdgeAccess!=NULL){

    COvlps *afing,*af;
    double from,to;

    from=contigTopEdgeAccess->beg;
    to=contigTopEdgeAccess->end;
    afing=ABuckets[whichA];
    // over all segments in contigs [whichA,0]
    for (af = afing; af != NULL && af->seg->b_contig == 0; af = af->Alink){
      // if the segment begins in the accessible interval
      if(af->seg->overlap->begpos >= 0 ) { //ahang is nonpositive, i.e. segment starts in A contig
	double pnt = af->seg->overlap->begpos;
	if (from <= pnt && pnt <= to) {
	  af->best = af->seg->overlap->length - MAX(0,-af->seg->overlap->begpos) - MAX(0,-af->seg->overlap->endpos);
	  assert(contigTopEdgeAccess->traceback==NULL);
	  af->trace = NULL;
	  // now process the segment's projections
	  if(Process_seg(af,bestTerm,nextAgapList,BgapLists,whichA,0,A,B,varwin)){
	    terminal=1;
	  }
	}
      }
    }

  }


  curr=*accessList;
  while(curr!=NULL){
    if(Process_Agap_one_accessible_interval(curr->ival,bestTerm,nextAgapList,BgapLists,whichA,A,B,varwin)){
      terminal=1;
    }
    curr=curr->next;
  }
  *accessList = cleanup_ilist(*accessList);
  return terminal;
}


static
int
Project_segments_across_Bgaps(COvlps **bestTerm, interval_list **BgapLists, interval_list **nextAgapList, int whichA, Scaffold *A, Scaffold *B, int varwin){
  int terminal=0;


  // To process the bottom edge of a B contig (in the context of a particular A contig)
  // For each accessible interval (actually, point)
  //   Project and add accessible interval on next B contig
  //   Project and add accessible interval on next A-gap's left edge

  int i;
  for(i=0;i<B->num_gaps;i++){
    interval_list *curr = BgapLists[i];
    assert(curr==NULL||curr->next==NULL); /* we are going to try to process backwards, as this adds A intervals in the right order ... */
    while(curr!=NULL){
      assert(curr->ival.traceback != NULL);
      if( curr->ival.traceback->seg->overlap->endpos < 0 ) {  // i.e. we come out the bottom
	// here, we consider either ended up at A gap or hitting another segment...

	double from = (A->ctgs[whichA].length+curr->ival.traceback->seg->overlap->endpos) + (B->gaps[i].gap_length-B->gaps[i].gap_var*varwin);
	double to = (A->ctgs[whichA].length+curr->ival.traceback->seg->overlap->endpos) + (B->gaps[i].gap_length+B->gaps[i].gap_var*varwin);

	// {from, to} are positions along the contig, in local (contig) coordinates
	{
	  COvlps *afing,*af;
	  afing=ABuckets[whichA];
	  // while there are segments on far side of horizontal gap, advance until we are in or past the relevant b-contig
	  while (afing != NULL && afing->seg->b_contig < i+1)
	    afing = afing->Alink;
	  // over all segments in contigs [whichA,i]
	  for (af = afing; af != NULL && af->seg->b_contig == i+1; af = af->Alink){
	    // if the segment begins in the accessible interval

	    if(af->seg->overlap->begpos >= 0 ) { //ahang is nonnegative, i.e. segment starts in A contig
	      double pnt = af->seg->overlap->begpos;
	      if (from <= pnt && pnt <= to) {
		int score = af->seg->overlap->length - MAX(0,-af->seg->overlap->begpos) - MAX(0,-af->seg->overlap->endpos);
		assert(score>0);
		if (curr->ival.traceback != NULL)
		  score += curr->ival.traceback->best;
		if (score > af->best){
		  af->best  = score;
		  af->trace = curr->ival.traceback;
		}
		// now process the segment's projections
		if(Process_seg(af,bestTerm,nextAgapList,BgapLists,whichA,i+1,A,B,varwin)){
		  terminal=1;
		}
	      }
	    }
	  }
	}
	if(from<to && to>A->ctgs[whichA].length){
	  interval to_add;
	  double gapFrac;

	  gapFrac =
	    ((double) (to-A->ctgs[whichA].length))/
	    (double) (B->gaps[i].gap_length + varwin*B->gaps[i].gap_var);
	  to_add.beg = B->ctgs[i+1].lft_end-gapFrac*B->gaps[i].gap_length;

	  to_add.end = B->ctgs[i+1].lft_end;
	  if(from > A->ctgs[whichA].length){
	    assert(B->gaps[i].gap_length - varwin*B->gaps[i].gap_var > 0);

	    gapFrac =
	      ((double) (from - A->ctgs[whichA].length)) /
	      (double) (B->gaps[i].gap_length - varwin*B->gaps[i].gap_var);

	    to_add.end -= gapFrac*B->gaps[i].gap_length;
	  }

	  to_add.traceback=curr->ival.traceback;

	  *nextAgapList = add_to_ilist(*nextAgapList,to_add);
	}

      } else { // we come out the side into A gap, in middle of B contig
	interval to_add;
	to_add.beg = B->ctgs[i].lft_end+B->ctgs[i].length - curr->ival.traceback->seg->overlap->endpos;
	to_add.end = to_add.beg;
	to_add.traceback = curr->ival.traceback;
	*nextAgapList = add_to_ilist(*nextAgapList,to_add);
      }

      curr=curr->prev;
    }

    BgapLists[i]=cleanup_ilist(BgapLists[i]);
  }

  // handle bottom edge (terminal cases) ...
  if(BgapLists[i]!=NULL){
    interval_list *curr=BgapLists[i];
    while(curr!=NULL){
      assert(curr->ival.traceback != NULL);

      if(curr->ival.traceback->seg->overlap->endpos <= 0 ) {  // i.e. we come out the bottom
	if(*bestTerm==NULL || ( curr->ival.traceback!=NULL && (*bestTerm)->best < curr->ival.traceback->best) ){
	  *bestTerm=curr->ival.traceback;
	}
	terminal=1;
      } else { // we come out in A gap
	interval to_add;
	to_add.beg = B->ctgs[i].lft_end+B->ctgs[i].length - curr->ival.traceback->seg->overlap->endpos;
	to_add.end = to_add.beg;
	to_add.traceback = curr->ival.traceback;
	*nextAgapList = add_to_ilist(*nextAgapList,to_add);
      }

      curr=curr->prev;
    }
    BgapLists[i]=cleanup_ilist(BgapLists[i]);
  }

  return terminal;

}


static
int
ProjectFromTopEdge(interval_list **thisAlist, COvlps **bestTerm,int whichA,Scaffold *A, Scaffold *B,int varwin,int bandbeg,int bandend){

  int terminal=0;
  double low,high,top,bot;
  double minX,maxX;
  double gapFrac;

  assert(*thisAlist==NULL);

  // find relevant interval in A-scaffold coordinates:
  low=MAX(bandbeg,A->ctgs[whichA].lft_end+A->ctgs[whichA].length);
  high=MIN(bandend,A->ctgs[whichA+1].lft_end);

  // if none, we are done
  if(low>=high){
    return 0;
  }

  // find the minimum and maximum amount into the B-scaffold we have to travel ...
  gapFrac = ((double)(A->ctgs[whichA+1].lft_end-low))/(double) A->gaps[whichA].gap_length;
  assert(gapFrac>=0&&gapFrac<=1);
  maxX = (gapFrac * A->gaps[whichA].gap_length + A->gaps[whichA].gap_var * varwin);
  gapFrac = ((double)(A->ctgs[whichA+1].lft_end-high))/(double) A->gaps[whichA].gap_length;
  assert(gapFrac>=0&&gapFrac<=1);
  minX = MAX(0,(gapFrac * A->gaps[whichA].gap_length - A->gaps[whichA].gap_var * varwin));

  // for the minimum, turn this into a position within the B scaffold, taking
  // gap stretchiness into account ...
  top=0;
  {

    double intoB=0;

    // what SHOULD we do if minimum length of gap is negative?

    // process a contig+gap at a time

    int i=0;
    while ( intoB < minX && i < B->num_gaps){

      intoB+=B->ctgs[i].length;
      if(intoB>=minX){
	top = B->ctgs[i].lft_end+B->ctgs[i].length-(intoB-minX);
	break;
      }

      intoB+=(B->gaps[i].gap_length+B->gaps[i].gap_var*varwin);

      if(intoB>=minX){
	gapFrac=((double)(intoB-minX))/(double) (B->gaps[i].gap_length+B->gaps[i].gap_var*varwin);
	assert(gapFrac>=0);
	top=B->ctgs[i+1].lft_end-(B->gaps[i].gap_length*gapFrac);
	break;
      }
      i++;
    }
    if(i==B->num_gaps){
      intoB+=B->ctgs[i].length;
      if(intoB>=minX){
	top = B->ctgs[i].lft_end+B->ctgs[i].length-(intoB-minX);
      } else {
	// top edge of interval out bottom of Scaffold B -- terminal, with no interval to add to A-gap access list
	return(1);
      }
    }
  }

  // for the maximum, turn this into a position within the B scaffold, taking
  // gap stretchiness into account ...
  {

    double intoB=0;

    // what SHOULD we do if minimum length of gap is negative?

    // process a contig+gap at a time

    int i=0;
    bot=0;
    while ( intoB < maxX && i < B->num_gaps){

      intoB+=B->ctgs[i].length;
      if(intoB>=maxX){
	bot = B->ctgs[i].lft_end+B->ctgs[i].length-(intoB-maxX);
	break;
      }

      intoB+=(B->gaps[i].gap_length - B->gaps[i].gap_var*varwin);

      if(intoB>=maxX){

	// no trouble with negative gap lengths here?
	assert(B->gaps[i].gap_length - B->gaps[i].gap_var*varwin > 0);

	gapFrac=((double)(intoB-maxX))/(double) (B->gaps[i].gap_length - B->gaps[i].gap_var*varwin);

	assert(gapFrac>=0);
	bot=B->ctgs[i+1].lft_end - (B->gaps[i].gap_length*gapFrac);
	break;
      }
      i++;
    }
    if(i==B->num_gaps){
      intoB+=B->ctgs[i].length;
      if(intoB>=maxX){
	bot = B->ctgs[i].lft_end+B->ctgs[i].length-(intoB-maxX);
      } else {
	// bottom edge of interval out bottom of Scaffold B -- terminal, with no interval to add to A-gap access list
	terminal=1;
	bot = B->ctgs[i].lft_end+B->ctgs[i].length;
      }
    }
  }
  assert(top<=bot);
  {
    interval to_add;
    to_add.beg=top;
    to_add.end=bot;
    to_add.traceback=NULL;
    *thisAlist = add_to_ilist(*thisAlist,to_add);
  }
  return terminal;
}

Segment *Align_Scaffold(Segment *seglist, int numsegs, int varwin,
                        Scaffold *AF, Scaffold *BF, int *best,
                        int bandbeg, int bandend) {

  // Simple initialization stuff at the top; comments on interesting parts further down ...

  interval_list *thisAlist=NULL,*nextAlist=NULL,**Blists=NULL;
  interval topEdgeAccess, *contigTopEdgeAccessPtr=NULL;
  int i;
  COvlps *bestTerm=NULL;
  int term=0;

  *best = -1;  // if not otherwise modified, signal no solution found

  Blists=(interval_list**)safe_malloc(sizeof(interval_list)*(BF->num_gaps+1));

  // setup (copied directly from Align_Scaffold() ):

  assert(bandbeg<=bandend);
  assert(bandend<=AF->length);
  assert(bandbeg>=-(BF->length));

  if (numsegs > MaxAlign)
    { MaxAlign = (int)(1.3*numsegs + 100);
    CtgOvls  = (COvlps *) safe_realloc(CtgOvls,sizeof(COvlps)*MaxAlign);
    }

  if (AF->num_gaps + BF->num_gaps + 2 > MaxBucket)
    { MaxBucket = (int)(1.3*(AF->num_gaps + BF->num_gaps + 2) + 100);
    ABuckets  = (COvlps **) safe_realloc(ABuckets,sizeof(COvlps *)*MaxBucket);
    }
  BBuckets  = ABuckets + (AF->num_gaps+1);

  { int i,c;
  Segment *s;


  for (i = 0; i <= AF->num_gaps; i++)
    ABuckets[i] = NULL;
  for (i = 0; i <= BF->num_gaps; i++)
    BBuckets[i] = NULL;

  c = numsegs;
  for (s = seglist; s != NULL; s = s->next)
    { c -= 1;
    CtgOvls[c].seg = s;
    CtgOvls[c].best = -1;
    CtgOvls[c].trace = NULL;

#ifdef DEBUG_SEGORDER
    fprintf(stderr,"CtgOvls[%d] actg: %d bctg: %d\n",
            c,CtgOvls[c].seg->a_contig,
            CtgOvls[c].seg->b_contig);
#endif
    // push segment onto Alink list; this needs to result in all
    // segments involving s->a_contig being linked together,
    // and the order of the elements should be such that
    // s->b_contig <= s->Alink->b_contig (if s->Alink != NULL)

    CtgOvls[c].Alink = ABuckets[s->a_contig];
    ABuckets[s->a_contig] = CtgOvls+c;
    if(ABuckets[s->a_contig]->Alink!=NULL)
      assert(ABuckets[s->a_contig]->seg->b_contig <= ABuckets[s->a_contig]->Alink->seg->b_contig);

    // original code did something similar for BBuckets and Blink,
    }


  // push segment onto Blink list; this needs to result in all
  // segments involving s->b_contig being linked together,
  // and the order of the elements should be such that
  // s->a_contig <= s->Blink->a_contig (if s->Blink != NULL)

  for(i=AF->num_gaps;i>=0;i--){
    COvlps *co;
    co = ABuckets[i];
    while(co!=NULL){
      co->Blink = BBuckets[co->seg->b_contig];
      BBuckets[co->seg->b_contig] = co;
      if(co->Blink!=NULL)
        assert(co->seg->a_contig <= co->Blink->seg->a_contig);
      co=co->Alink;
    }
  }

  }

#ifdef DEBUG_ALIGN
  { Segment *s;
  COvlps  *c;
  int      i;

  fprintf(stderr,"\nAlign Scaffolds\n\n  Seglist:\n");
  for (s = seglist; s != NULL; s = s->next)
    fprintf(stderr,"    (%d,%d)\n",s->a_contig,s->b_contig);
  fprintf(stderr,"\n  A-Buckets:\n");
  for (i = 0; i <= AF->num_gaps; i++)
    { fprintf(stderr,"    %2d:",i);
    for (c = ABuckets[i]; c != NULL; c = c->Alink)
      fprintf(stderr," %d",c->seg->b_contig);
    fprintf(stderr,"\n");
    }
  fprintf(stderr,"\n  B-Buckets:\n");
  for (i = 0; i <= BF->num_gaps; i++)
    { fprintf(stderr,"    %2d:",i);
    for (c = BBuckets[i]; c != NULL; c = c->Blink)
      fprintf(stderr," %d",c->seg->a_contig);
    fprintf(stderr,"\n");
    }
  fprintf(stderr,"\n");
  }
#endif

  // Interesting stuff:
  // need to be able to:
  // append intervals == add_to_ilist()
  // handle a reached segment == Process_seg()
  // process the left edge of an A-gap (from list of accessible intervals to list) == Project_across_Agap()
  // process the right edge of an A-gap (from list of accessible intervals ) == Process_Agap_accessible_intervals()
  //     this includes - reaching segments on left edge of A contig and processing appropriately
  //                   - reaching segments on top edge of B contig and processing approp.
  //                   - reaching left side of following Agap and processing approp.
  // process the bottom edge of a B-gap: list of accessible intervals along an A contig == Project_segments_across_Bgaps()

  // Given these, top level control is:
  // Initialize 0'th A-gap's right edge with accessible interval based on banding
  // For each A-gap
  //   Initialize B-contigs' accessibility (really, accessibility of A-contig as relevant to each B-contig)
  //     For first B-contig, initialize accessible portion of A-contig based on banding
  //     Null for all others
  //   Process right edge of A-gap (modifying B-contigs' accessibility lists and following A-gap's left-edge accessibility list)
  //     Check for terminal solutions out bottom
  //   For each B contig,
  //     Process accessibility list (modifying following B contigs' lists and also following A-gap's left-edge accessibility list)
  //     If final contig,
  //       Check for terminal solutions
  //   If not final gap,
  //     Process left edge of following A-gap (creating accessibility list for right edge of that gap)
  //   Else
  //     Check for terminal solutions

  // Room for improvement: For a given A-gap, if we could process the A-contig-bottom-exit-points for each B contig
  // at just the right time, then interval additions would be guaranteed monotonic so that add_to_ilist() could always
  // just worry about the tail of the current list.  This would sometimes require processing of the bottom edge inside the processing
  // if an A-gap accessible interval (between processing the portion on a contig and the portion in the following gap) and sometimes
  // require processing between A-gap accessible intervals.


  // negative gap lengths cause serious complications.  Get rid of them!
  // Also, there is some trouble with rounding errors on coordinates
  {
    int i;
    int extra=0;
    int extraB=0;
    int into;

    i=0;
    into=0;
    extra=0;

    while(i<AF->num_gaps){
      int diff=AF->ctgs[i+1].lft_end+extra-(AF->ctgs[i].lft_end+AF->ctgs[i].length);
      into+=AF->ctgs[i].length;
      if(diff<1){
        // if a band value occurs before the gap starts, then adjust
        if(bandbeg>AF->ctgs[i].length+AF->ctgs[i].lft_end){
          bandbeg+=1-diff;
        }
        if(bandend>AF->ctgs[i].length+AF->ctgs[i].lft_end){
          bandend+=1-diff;
        }
        extra+=1-diff;
        AF->gaps[i].gap_length=1;
      }
      AF->ctgs[i+1].lft_end+=extra;
      AF->gaps[i].gap_length=MAX(1,diff);
      assert(AF->ctgs[i].lft_end+AF->ctgs[i].length+AF->gaps[i].gap_length == AF->ctgs[i+1].lft_end);
      into+=AF->gaps[i].gap_length;
      i++;
    }

    i=0;
    into=0;
    extra=0;
    while(i<BF->num_gaps){
      int diff=BF->ctgs[i+1].lft_end+extra-(BF->ctgs[i].lft_end+BF->ctgs[i].length);
      into+=BF->ctgs[i].length;
      if(diff < 1 ){
        // if a band value occurs before the gap starts, then adjust
        if(-bandbeg>BF->ctgs[i].length+BF->ctgs[i].lft_end){
          bandbeg-=1-diff;
        }
        if(-bandend>BF->ctgs[i].length+BF->ctgs[i].lft_end){
          bandend-=1-diff;
        }
        extra+=1-diff;
        BF->gaps[i].gap_length=1;
      }
      BF->ctgs[i+1].lft_end += extra;
      BF->gaps[i].gap_length=MAX(diff,1);
      assert(BF->ctgs[i].lft_end+BF->ctgs[i].length+BF->gaps[i].gap_length == BF->ctgs[i+1].lft_end);
      into+=BF->gaps[i].gap_length;
      i++;
    }

  }

  // Initialize along left edge of first A contig
  if(bandbeg<0){
    interval startupAgap;
    int top=-MIN(0,bandend);
    int bot=-MIN(0,bandbeg);
    int beg;
    int end;
    int intoB=0;
    int i;
    i=0;
    while(i<BF->num_gaps){
      intoB+=BF->ctgs[i].length;
      if(intoB>top){
	beg=top;
	break;
      }
      intoB+=BF->gaps[i].gap_length;
      if(intoB>top){
	double gapFrac = ((double) (top - (BF->ctgs[i].lft_end+BF->ctgs[i].length))) /
          (double) (BF->gaps[i].gap_length+varwin*BF->gaps[i].gap_var);
	assert(gapFrac>=0);
	beg = BF->ctgs[i].lft_end+BF->ctgs[i].length + gapFrac * BF->gaps[i].gap_length;
	break;
      }
    }
    if(i==BF->num_gaps){
      beg=top;
    }
    i=0;
    intoB=0;
    while(i<BF->num_gaps){
      intoB+=BF->ctgs[i].length;
      if(intoB>bot){
	end=bot;
	break;
      }
      intoB+=BF->gaps[i].gap_length;
      if(intoB>bot){
	double gapFrac = ((double)(intoB-bot)) / (double)(BF->gaps[i].gap_length+varwin*BF->gaps[i].gap_var);
	assert(gapFrac>=0);
	end = BF->ctgs[i+1].lft_end - gapFrac * BF->gaps[i].gap_length;
	break;
      }
      i++;
    }
    if(i==BF->num_gaps){
      end=bot;
    }
    assert(end>=beg);
    startupAgap.end=end;
    startupAgap.beg=beg;
    startupAgap.traceback=NULL;
    thisAlist = add_to_ilist(NULL,startupAgap);
  }

  // For each A-gap
  for(i=0;i<=AF->num_gaps;i++){
    int j;
    //   Initialize B-contig bottom-edge accessibility (really, accessibility of A-contig as relevant to each B-contig) to NULL
    for(j=0;j<=BF->num_gaps;j++){
      Blists[j]=NULL; /* was: cleanup_ilist(Blists[j]);*/
    }

    //     For first B-contig, initialize accessible portion of A-contig based on banding
    {
      int left,right;
      left = MAX(bandbeg,AF->ctgs[i].lft_end);
      right = MIN(bandend,AF->ctgs[i].lft_end+AF->ctgs[i].length);
      left-=AF->ctgs[i].lft_end;
      right-=AF->ctgs[i].lft_end;
      if(left<=right){
	topEdgeAccess.beg=left;
	topEdgeAccess.end=right;
	topEdgeAccess.traceback=NULL;
	contigTopEdgeAccessPtr=&topEdgeAccess;
      }else {
	contigTopEdgeAccessPtr=NULL;
      }
    }

    // thisAlist was constructed so that we always have the tail; trace backwards to get the head, as processing is assumed to be from head to tail
    while(thisAlist!=NULL&&thisAlist->prev!=NULL){
      thisAlist=thisAlist->prev;
    }
    //   Process right edge of A-gap (modifying B-contigs' accessibility lists and following A-gap's left-edge accessibility list)
    if(Process_Agap_accessible_intervals(&thisAlist,
					 contigTopEdgeAccessPtr,
					 &nextAlist,
					 Blists,
					 &bestTerm,
					 i,
					 AF,
					 BF,
					 varwin)){
      term=1;
    }

    //   For each B contig,
    //     Process accessibility list (modifying following B contigs' lists and also following A-gap's left-edge accessibility list), and check for terminal solutions out bottom
    for(j=0;j<=BF->num_gaps;j++){
      if(Blists[j]==NULL)continue;
      interval_list *tmp=Blists[j];
      assert(Blists[j]->next==NULL); /* we want to process these guys from tail to head, as this will potentially introduce A gap intervals from top to bottom */
      if(Project_segments_across_Bgaps(&bestTerm,Blists,&nextAlist,i,AF,BF,varwin)){
	term=1;
      }
    }

    // we are done with current interval list; reset to prepare for next gap
    thisAlist=cleanup_ilist(thisAlist);

    if( i < AF->num_gaps ){
      // process top edge of following a-gap (creating accessibility list for right edge of that gap)
      if( ProjectFromTopEdge(&thisAlist,&bestTerm,i,AF,BF,varwin,bandbeg,bandend)){
	term=1;
      }
      //     Process left edge of following A-gap (creating accessibility list for right edge of that gap)
      if(nextAlist!=NULL){
	while(nextAlist->prev!=NULL)nextAlist=nextAlist->prev;
	if( Project_across_Agap(nextAlist,&thisAlist,&bestTerm,AF->gaps[i].gap_length,AF->gaps[i].gap_var,BF,varwin) ){
	  term=1;
	}
	nextAlist=cleanup_ilist(nextAlist);
      }
    }

  }

  // For final A contig, check for terminal solutions
  while(nextAlist!=NULL){
    term=1;
    if(bestTerm==NULL || (nextAlist->ival.traceback!=NULL && bestTerm->best < nextAlist->ival.traceback->best )){
      bestTerm=nextAlist->ival.traceback;
    }
    nextAlist=nextAlist->prev; // we start from tail and work backwards
  }

  // Now, we need to use a solution, if any was found
  if(term){

    Segment *s, *r;
    COvlps  *c;

    // if found, and there was a best terminal segment ...
    if(bestTerm!=NULL){

      // get its score
      *best=bestTerm->best;

      // set things up to protect essential segments from being freed
      c= bestTerm;
      while(c!=NULL){
	c->seg->alow = - (c->seg->alow+1);
	c=c->trace;
      }

    } else {
      *best=0;
    }

    // free inessential segments and restore essential segments
    for (s = seglist; s != NULL; s = r)
      { r = s->next;
      if (s->alow >= 0)
	{ safe_free(s->overlap);
	safe_free(s);
	}
      else
	s->alow = - (s->alow+1);
      }

    // invert the list to set up seglist
    r = NULL;
    for (c = bestTerm; c != NULL; c = c->trace){
      s = c->seg;
      s->next = r;
      r = s;
    }

    seglist = r;

  } else {

    // see various notes towards end of obsolete/align_scaffold_old
    // for conditions on freeing elements

    seglist=NULL;

  }

  return (seglist);
}
