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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>

#include "GF_ALN_local.h"

#define  max(x,y)        ((x<y) ? (y):(x))
#define  min(x,y)        ((x>y) ? (y):(x))


#define MIN_USABLE 3  /* Smallest subpart of a segment usable for chaining */


int MIN_ALIGNED_COLS=30; /* minimum length of a local overlap in the following
			    sense: an overlap is defined by a set of segments;
			    each segment has a length (ignoring minor troubles
			    in determining length in the presence of small
			    indels); let the sum of the overlap for this
			    purpose be the sum of the lengths of the segments
			    that make up the overlap */

#define BIG_INT 0x7FFFFFFF

typedef struct {
  int  start;
  int  base;
  int  segment;
  int  best;
} Candidate;

typedef struct {
  Local_Segment *item;
  int isadd;
} Event;

typedef struct {
  int value;
  int source;
  int start;
  int colsAligned;
} TraceElement;

/*** AVL-TREE LIST ROUTINES ***/

static void OutOfMemory(char *where)
{ fprintf(stderr,"COMPARE_LOCAL: Out of memory (%s)\n",where);
  exit (1);
}

typedef struct _AVLnode {
  int               RC, LN;
  short             H;
  struct _AVLnode  *L, *R;
  Candidate         V;
} AVLnode;

#define CP(v) ((v)->L->LN)

static AVLnode *freept;
static AVLnode *NIL;

#define INC  AVLinc
#define DEC  AVLdec
#define SEL  AVLselect
#define RNK  AVLrank
#define ADD  AVLinsert
#define DEL  AVLdelete

static void AVLinit(void)
{ freept  = NULL;
  NIL     = (AVLnode *) malloc(sizeof(AVLnode));
  if (NIL == NULL)
    OutOfMemory("Candidate list");
  NIL->LN = NIL->RC = 1;
  NIL->H  = 0;
  NIL->V.base = BIG_INT;
  NIL->V.best = BIG_INT;
}

static AVLnode *AVLinc(AVLnode *v)
{ v->RC++;
  return (v);
}

static void AVLdec(AVLnode *v)
{ v->RC--;
  if (v->RC == 0)
    { DEC(v->L);
      DEC(v->R);
      v->L = freept;
      freept = v;
    }
}

static int AVLlength(AVLnode *v)
{ DEC(v);return (v->LN - 1); }

static AVLnode *NEW(AVLnode *l, Candidate *x, AVLnode *r)
{ AVLnode *v;
  int b;

  if (freept == NULL)
    { v = (AVLnode *) malloc(sizeof(AVLnode));
      if (v == NULL)
        OutOfMemory("Candidate list");
    }
  else
    { v = freept;
      freept = v->L;
    }

  v->RC = 1;
  v->V  = *x;
  v->L  = l;
  v->R  = r;
  v->LN = l->LN + r->LN;
  v->H  = (l->H < r->H ? r->H : l->H) + 1;

  b = v->V.base;
  if (v->L->V.best < b)
    b = v->L->V.best;
  if (v->R->V.best < b)
    b = v->R->V.best;
  v->V.best = b;

  return (v);
}

static AVLnode *BAL(AVLnode *l, Candidate *x, AVLnode *r)
{ AVLnode *t;
  if (l->H - r->H >= -1 && l->H - r->H <= 1)
    t = NEW(INC(l),x,INC(r));
  else if (l->H > r->H)
    if (l->L->H >= l->R->H)
      t = NEW(INC(l->L),&(l->V),NEW(INC(l->R),x,INC(r)));
    else
      t = NEW(NEW(INC(l->L),&(l->V),INC(l->R->L)),&(l->R->V),
              NEW(INC(l->R->R),x,INC(r)));
  else
    if (r->R->H >= r->L->H)
      t = NEW(NEW(INC(l),x,INC(r->L)),&(r->V),INC(r->R));
    else
      t = NEW(NEW(INC(l),x,INC(r->L->L)),&(r->L->V),
              NEW(INC(r->L->R),&(r->V),INC(r->R)));
  DEC(l); DEC(r);
  return (t);
}

static Candidate *AVLselect(AVLnode *v, int k)
{ Candidate *x;
  if (k < CP(v))
    x = SEL(INC(v->L),k);
  else if (k > CP(v))
    x = SEL(INC(v->R),k-CP(v));
  else
    x = &(v->V);
  DEC(v);
  return (x);
}

static int AVLrank(AVLnode *v, int pos)
{ int k;
  if (v == NIL)
    k = 0;
  else if (pos < v->V.start)
    k = RNK(INC(v->L),pos);
  else
    k = CP(v) + RNK(INC(v->R),pos);
  DEC(v);
  return (k);
}

static AVLnode *AVLminprf(AVLnode *v, int hgh, int bst)
{ AVLnode *r;
  int      b = 0;

  if (v == NIL)
    r = v;
  else if (hgh < v->V.start)
    r = AVLminprf(INC(v->L),hgh,bst);
  else
    { if (v->L->V.best < bst)
        b = v->L->V.best;
      if (v->V.base < bst)
        b = v->V.base;
      r = AVLminprf(INC(v->R),hgh,b);
      if (r->V.base < bst)
        bst = r->V.base;
      if (v->V.base < bst)
        { r = v; bst = v->V.base; }
      if (v->L->V.best < bst)
        r = AVLminprf(INC(v->L),hgh,bst);
    }
  DEC(v);
  return (r);
}

static AVLnode *AVLminsuf(AVLnode *v, int low, int bst)
{ AVLnode *r;
  int      b = 0;

  if (v == NIL)
    r = v;
  else if (low > v->V.start)
    r = AVLminsuf(INC(v->R),low,bst);
  else
    { if (v->R->V.best < bst)
        b = v->R->V.best;
      if (v->V.base < bst)
        b = v->V.base;
      r = AVLminsuf(INC(v->L),low,b);
      if (r->V.base < bst)
        bst = r->V.base;
      if (v->V.base < bst)
        { r = v; bst = v->V.base; }
      if (v->R->V.best < bst)
        r = AVLminsuf(INC(v->R),low,bst);
    }
  DEC(v);
  return (r);
}

static AVLnode *AVLminrng(AVLnode *v, int low, int hgh)
{ AVLnode *r, *t;

  if (v == NIL)
    r = v;
  else if (hgh < v->V.start)
    r = AVLminrng(INC(v->L),low,hgh);
  else if (low > v->V.start)
    r = AVLminrng(INC(v->R),low,hgh);
  else
    { r = v;
      t = AVLminprf(INC(v->R),hgh,r->V.base);
      if (t->V.base < r->V.base)
        r = t;
      t = AVLminsuf(INC(v->L),low,r->V.base);
      if (t->V.base < r->V.base)
        r = t;
    }
  DEC(v);
  return (r);
}

static AVLnode *AVLinsert(AVLnode *v, int k, Candidate *x)
{ AVLnode *t;
  if (v == NIL)
    t = BAL(INC(NIL),x,INC(NIL));
  else if (k < CP(v))
    t = BAL(ADD(INC(v->L),k,x),&(v->V),INC(v->R));
  else
    t = BAL(INC(v->L),&(v->V),ADD(INC(v->R),k-CP(v),x));
  DEC(v);
  return (t);
}

static AVLnode *AVLdelete(AVLnode *v, int k)
{ AVLnode *t;
  if (v->L == NIL && v->R == NIL)
    t = INC(NIL);
  else if (k <= CP(v) && v->L != NIL)
    if (k == CP(v))
      t = BAL(DEL(INC(v->L),k-1),SEL(INC(v->L),k-1),INC(v->R));
    else
      t = BAL(DEL(INC(v->L),k),&(v->V),INC(v->R));
  else
    if (k == CP(v))
      t = BAL(INC(v->L),SEL(INC(v->R),1),DEL(INC(v->R),1));
    else
      t = BAL(INC(v->L),&(v->V),DEL(INC(v->R),k-CP(v)));
  DEC(v);
  return (t);
}



static int SSORT(const void *l, const void *r)
{ Event *x, *y;
  int ax, ay, bx, by;

  x = (Event *) l;
  y = (Event *) r;
  if (x->isadd)
    { ax = x->item->abpos;
      bx = x->item->bbpos;
    }
  else
    { ax = x->item->aepos;
      bx = x->item->bepos;
    }
  if (y->isadd)
    { ay = y->item->abpos;
      by = y->item->bbpos;
    }
  else
    { ay = y->item->aepos;
      by = y->item->bepos;
    }
  if (ax < ay)
    return (-1);
  else if (ax > ay)
    return (1);
  else if (x->isadd != y->isadd)
    return (x->isadd - y->isadd);
  else
    return (bx - by);
}


static void convert_segs(Local_Segment *Segs,int NumSegs,int comp, int Alen,int Blen)
{ int i;   /* Mark and reverse all complemented local segs */

 if (comp)
   for (i = 0; i < NumSegs; i++)
     { Segs[i].bbpos = Blen - Segs[i].bbpos;
     Segs[i].bepos = Blen - Segs[i].bepos;
     }

 for (i = 0; i < NumSegs; i++)
   if (Segs[i].bbpos > Segs[i].bepos)
     { int x;
     x = Segs[i].bbpos;
     Segs[i].bbpos = Segs[i].bepos;
     Segs[i].bepos = x;
     Segs[i].score = -Segs[i].score-1;
     }
}

static void restore_segs(Local_Segment *Segs,int NumSegs,int comp,int Alen,int Blen)  
{ int i;   /* Unmark and reverse all complemented local segs */

 for (i = 0; i < NumSegs; i++)
   if (Segs[i].score < 0)
     { int x;
     x = Segs[i].bbpos;
     Segs[i].bbpos = Segs[i].bepos;
     Segs[i].bepos = x;
     Segs[i].score = -Segs[i].score-1;
     }

 if (comp)
   { for (i = 0; i < NumSegs; i++)
     { Segs[i].bbpos = Blen - Segs[i].bbpos;
     Segs[i].bepos = Blen - Segs[i].bepos;
     }
   }
}


Local_Overlap *Find_Local_Overlap(int Alen, int Blen, int comp, int nextbest,
                                  Local_Segment *Segs, int NumSegs,
                                  int MinorThresh, double GapThresh)
{ static Candidate Cvals;
  static int MaxTrace = -1;
  static TraceElement *Trace = NULL;
  static Event        *EventList;
  Local_Overlap *Descriptor;
  Local_Chain   *Chain;

  if (NumSegs == 0) return (NULL);

  if (nextbest)
    { if (Trace == NULL) return (NULL);
      convert_segs(Segs,NumSegs,comp,Alen,Blen);
      goto Gen_Overlap;
    }

  if (MaxTrace < 0)
    AVLinit();

  if (NumSegs > MaxTrace)
    { MaxTrace  = (int)(1.3*NumSegs) + 500;
    Trace = (TraceElement *)
      realloc(Trace,(sizeof(Event)+2*sizeof(TraceElement))*MaxTrace);
    if (Trace == NULL)
      OutOfMemory("Overlap Trace Array");
    EventList = (Event *) (Trace + MaxTrace);
    { // We have to make sure that EventList is aligned on an appropriate boundary.
      // It is derived from Trace which has looser alignment constraints.
      long address = (long)EventList; 
      // By convention "long" int is big as the size of a pointer
      long offset = (address % sizeof(void *));  
      int pad = sizeof(void *) - offset; 
      // This is how much we need to add to get things aligned.
      if(offset){
	//	fprintf(stderr,"* Eventlist is %p adding %d up to ", EventList, pad);
	EventList = (Event *)(((char *)EventList) + pad);
	//	fprintf(stderr," %p\n", EventList);
      }
    }
    }

  convert_segs(Segs,NumSegs,comp,Alen,Blen);

  { int i;
    for (i = 0; i < NumSegs; i++)
      { EventList[2*i].item = Segs+i;
        EventList[2*i].isadd = 1;
        EventList[2*i+1].item = Segs+i;
        EventList[2*i+1].isadd = 0;
      }
  }

  qsort(EventList,2*NumSegs,sizeof(Event),SSORT);

  { int e;
    AVLnode *elist, *ilist, *olist;

    elist = AVLinc(NIL);
    ilist = AVLinc(NIL);
    olist = AVLinc(NIL);

    for (e = 0; e < 2*NumSegs; e++)
      { int i, bb, be, ab, ae;
        double err;

        /* Determine least gapped path to i'th segment */

        i = EventList[e].item - Segs;
        bb = Segs[i].bbpos;
        be = Segs[i].bepos;
        ab = Segs[i].abpos;
        ae = Segs[i].aepos;
	err = Segs[i].error;

        if (EventList[e].isadd)  /* Segment begins */
          { int clen, best, srce;

	  // this definition of best differs from the original (below)
	  // it is designed to encourage global alignment
          //
            best = ab+bb;                 /* Best from boundary */

            //best = ab;                 /* Best from boundary */
            //if (best > bb)
            //  best = bb;
            //best *= 2;


            srce = -1;
            clen = AVLlength(AVLinc(elist));

            { int p;    /* Examine bests from elist */

              p = AVLrank(AVLinc(elist),bb);  /* Best @ start of seg */
              if (p > 0)
                { Candidate *cand;
                  int        altr;
  
                  cand = AVLselect(AVLinc(elist),p);
                  altr = cand->base + (ab + bb);
                  if (altr < best)
                    { best = altr;
                      srce = cand->segment;
                    }
                }

              while (++p <= clen)         /* Bests @ midpoints of seg */
                { Candidate *cand;
                  int        altr;
            
                  cand = AVLselect(AVLinc(elist),p);
                  if (cand->start > be - MIN_USABLE) break;
                  altr = cand->base + 2*cand->start + (ab - bb);
                  if (altr < best)
                    { best = altr;
                      srce = cand->segment;
                    }
                }
            }

            /* Examine bests from ilist and olist */

            { AVLnode *m;
              int bdiag, ldiag, altr;

              bdiag = bb - ab;
              ldiag = bdiag + ((ae-ab) - MIN_USABLE);
              m = AVLminprf(AVLinc(ilist),bdiag,BIG_INT);
              if (m != NIL)
                { altr = m->V.base + bdiag;
                  if (altr < best)
		    { srce = m->V.segment;
                      best = altr;
                    }
                }
              m = AVLminrng(AVLinc(olist),-ldiag,-bdiag);
              if (m != NIL)
                { altr = m->V.base - bdiag;
                  if (altr < best)
		    { srce = m->V.segment;
                      best = altr;
                    }
                }
            }

            /* Record best linkage for segment */

            Trace[i].value  = best;
            Trace[i].source = srce;

	    Trace[i].colsAligned = (int)((1.-err)*(double)(min(ae-ab,be-bb)+1));
            if (srce >= 0){
              Trace[i].start = Trace[srce].start;
	      Trace[i].colsAligned += Trace[srce].colsAligned;
	    } else
              Trace[i].start = i;

            /* Add segment to ilist and olist */

            { int p, d;

              d = be - ae;
              Cvals.segment = i;

              Cvals.start = d;
              Cvals.base  = best - d;
              p = AVLrank(AVLinc(ilist),d);
              ilist = AVLinsert(ilist,p,&Cvals);  

              d = -d;
              Cvals.start = d;
              Cvals.base  = best - d;
              p = AVLrank(AVLinc(olist),d);
              olist = AVLinsert(olist,p,&Cvals);
	    }
	  }

        else  /* Segment ends */
          { int best, clen;

            best = Trace[i].value;
            clen = AVLlength(AVLinc(elist));

            /* Add candidate (if any) created by i'th segment */

            { Candidate *cand;
              int p, off;

              off = be + Segs[i].aepos;
              p = AVLrank(AVLinc(elist),be);
              if (p != 0)
                cand = AVLselect(AVLinc(elist),p);
              if (p == 0 || best < cand->base + off)
                { p += 1;
                  while (p <= clen)
                    { cand = AVLselect(AVLinc(elist),p);
                      if (cand->base + off < best) break;
                      elist = AVLdelete(elist,p);
                      clen -= 1;
                    }
                  p -= 1;
                  if (p > 0)
                    { cand = AVLselect(AVLinc(elist),p);
		      if (cand->start == be)
                        elist = AVLdelete(elist,p--);
                    }
                  Cvals.start   = be;
                  Cvals.base    = best - off;
                  Cvals.segment = i;
                  elist = AVLinsert(elist,p,&Cvals);  
                }
            }

            /* Remove candidates from ilist and olist */

            { int p, d;

              d = be-ae;
              p = AVLrank(AVLinc(ilist),d);
              while (AVLselect(AVLinc(ilist),p)->segment != i)
                p -= 1;
              ilist = AVLdelete(ilist,p);

              p = AVLrank(AVLinc(olist),-d);
              while (AVLselect(AVLinc(olist),p)->segment != i)
                p -= 1;
              olist = AVLdelete(olist,p);
            }
          } 

      }

    AVLdec(elist);
    AVLdec(ilist);
    AVLdec(olist);
  }



Gen_Overlap:

  { int i, npiece;
    int best, end, beg;

    best = BIG_INT;   /* Determine best overall overlap */
    end  = -1;
    for (i = 0; i < NumSegs; i++){
      //      if (Trace[i].start >= 0)
      if (Trace[i].start >= 0&&Trace[i].colsAligned >= MIN_ALIGNED_COLS) {
        int sfx;

	//  this definition of sfx differs from the original (below)
	//  it is designed to encourage global alignment
        //
        sfx = Alen - Segs[i].aepos + Blen - Segs[i].bepos;

        //sfx = Alen - Segs[i].aepos;
        //if (Blen - Segs[i].bepos < sfx)
        //  sfx = Blen - Segs[i].bepos;
        //sfx *= 2;

        //  The "- 2 * Trace[i].colsAligned" makes us encourage longer alignments
        //
        if (Trace[i].value + sfx - 2*Trace[i].colsAligned < best) {
          best = Trace[i].value - 2*Trace[i].colsAligned + sfx;
          end  = i;
        }
      }  
    }

    if (end < 0) {
      restore_segs(Segs,NumSegs,comp,Alen,Blen);
      return (NULL);
    }

    beg = Trace[end].start;

    /* How many segments in the best overlap? */

    npiece = 0;
    for (i = end; i >= 0; i = Trace[i].source)
      npiece += 1;



    /* Allocate result data structures in a single memory block */

    Descriptor = (Local_Overlap *) malloc(sizeof(Local_Overlap) +
                                          (npiece+1)*sizeof(Local_Chain));
    if (Descriptor == NULL)
      OutOfMemory("Overlap descriptor");
    Chain = (Local_Chain *) (Descriptor + 1);

    /* Fill out the description of the chain */

    { int n;
      n = npiece;
      for (i = end; i >= 0; i = Trace[i].source)
        Chain[--n].piece = Segs[i];
    }
    

#define ALLOW_DUP_SEGS_IN_NEXT  /* allow all but the first segment to be
				   used in later attempts */
#ifndef ALLOW_DUP_SEGS_IN_NEXT
    for (i = 0; i < NumSegs; i++)
      if (Trace[i].start == beg)
        Trace[i].start = -1; /* this seems to prevent reuse of segments
				in subsequent calls, and/or if we reject
				this segment as too noisy and jump back
				to the top */
#else
    Trace[end].start = -1; /* this seems to prevent reuse of segments
				in subsequent calls, and/or if we reject
				this segment as too noisy and jump back
				to the top */
  #define REUSE_CURRENT_LAST_AS_NONTERMINAL_SEG
  #ifndef REUSE_CURRENT_LAST_AS_NONTERMINAL_SEG
    for (i = 0; i < NumSegs; i++)
      if (Trace[i].source == beg)
        Trace[i].source = -1; 
  #endif

#endif


    {// The last segment doesn't describe an alignment, only a gap. Initialize it to reasonable values
      Local_Segment *lastseg = &Chain[npiece].piece;
      lastseg->abpos = lastseg->bbpos = -1;
      lastseg->aepos = lastseg->bepos = -1;
      lastseg->ldiag = lastseg->hdiag = -1;
      lastseg->score = -1;
      lastseg->error = -1.0;
    }

    { int gl;
    /* there's basically a bug here: abpos = 1 means starts at first char of A;
       so, agap should be 0, but gets set to 1; i.e., every first gap size 
       gets set to one too many; but there's existing code that relies on
       this fact, so leave it alone for now */
      gl = Chain[0].piece.abpos;
      if (gl > Chain[0].piece.bbpos)
        gl = Chain[0].piece.bbpos;
      Chain[0].agap = gl;
      Chain[0].bgap = gl;
    }

    for (i = 1; i < npiece; i++)
      { Chain[i].agap = Chain[i].piece.abpos - Chain[i-1].piece.aepos;
        Chain[i].bgap = Chain[i].piece.bbpos - Chain[i-1].piece.bepos;
      }

    { int gl;
      gl = Alen - Chain[npiece-1].piece.aepos;
      if (gl > Blen - Chain[npiece-1].piece.bepos)
        gl = Blen - Chain[npiece-1].piece.bepos;
      Chain[npiece].agap = gl;
      Chain[npiece].bgap = gl;
    }

    for (i = 0; i <= npiece; i++)
      { if (abs(Chain[i].agap) <= MinorThresh)
          { if (abs(Chain[i].bgap) <= MinorThresh)
              { if (Chain[i].agap != 0 || Chain[i].bgap != 0)
                  Chain[i].type = LOCAL_MINOR;
                else
                  Chain[i].type = LOCAL_BOUNDARY;
              }
            else if (Chain[i].bgap < 0)
              Chain[i].type = LOCAL_REPEAT;
            else if (Chain[i].bgap > 4*Chain[i].agap)
              Chain[i].type = LOCAL_INDEL;
            else
              Chain[i].type = LOCAL_DISAGREE;
          }
        else if (Chain[i].agap < 0)
          { if (Chain[i].bgap < MinorThresh)
              Chain[i].type = LOCAL_REPEAT;
            else
              Chain[i].type = LOCAL_REPnDEL; 
          }
        else
          { if (abs(Chain[i].bgap) < MinorThresh)
              if (Chain[i].agap > 4*Chain[i].bgap)
                Chain[i].type = LOCAL_INDEL;
              else
                Chain[i].type = LOCAL_DISAGREE;
            else if (Chain[i].bgap < 0)
              Chain[i].type = LOCAL_REPnDEL;
            else
              Chain[i].type = LOCAL_DISAGREE;
          }
      }

    /* Fill out overlap descriptor */

    Descriptor->num_pieces = npiece;
    Descriptor->score      = best;
    Descriptor->chain      = Chain;
    Descriptor->comp       = comp;
    
    { Local_Segment *sg;
      int            ln;

      Descriptor->indif = 0;

      for (i = 0; i < npiece; i++)
        { sg = &(Chain[i].piece);
          ln = ((sg->aepos - sg->abpos) + (sg->bepos - sg->bbpos)) / 2;
          if (i > 0 && Chain[i-1].piece.error < sg->error)
            { if (Chain[i].agap < Chain[i].bgap)
                { if (Chain[i].agap < 0)
                    ln += Chain[i].agap;
                }
              else
                { if (Chain[i].bgap < 0)
                    ln += Chain[i].bgap;
                }
            }
          if (i < npiece-1 && Chain[i+1].piece.error <= sg->error) 
            { if (Chain[i+1].agap < Chain[i+1].bgap)
                { if (Chain[i+1].agap < 0)
                    ln += Chain[i+1].agap;
                }
              else
                { if (Chain[i+1].bgap < 0)
                    ln += Chain[i+1].bgap;
                }
            }
          if (ln > 0)
            Descriptor->indif += (int)(ln * sg->error);
        }
    }

    Descriptor->diffs = Descriptor->indif;

    for (i = 0; i <= npiece; i++)
      { int d; 
        if (Chain[i].agap < 0 || Chain[i].bgap < 0)
          d = abs( (Chain[i].piece.bbpos - Chain[i].piece.abpos) -
                   (Chain[i-1].piece.bepos - Chain[i-1].piece.aepos));
        else
          { d = Chain[i].agap;
            if (d < Chain[i].bgap) d = Chain[i].bgap;
          }
        Descriptor->diffs += d;
      }

    { int overa, overb;
      overa = (Chain[npiece-1].piece.aepos + Chain[npiece].agap)
            - (Chain[0].piece.abpos - Chain[0].agap);
      overb = (Chain[npiece-1].piece.bepos + Chain[npiece].bgap)
            - (Chain[0].piece.bbpos - Chain[0].bgap);
      Descriptor->length = (overa + overb) / 2;
    }

    Descriptor->begpos = Chain[0].piece.abpos - Chain[0].piece.bbpos;
    Descriptor->endpos = (Blen - Chain[npiece-1].piece.bepos)
                       - (Alen - Chain[npiece-1].piece.aepos);

    for (i = 0; i < npiece; i++)
      if (Chain[i].piece.score < 0)
        { int x;
          x = Chain[i].piece.bbpos;
          Chain[i].piece.bbpos = Chain[i].piece.bepos;
          Chain[i].piece.bepos = x;
          Chain[i].piece.score = - Chain[i].piece.score-1;
          Chain[i].reversed = 1;
        }
      else
        Chain[i].reversed = 0;
  }

  restore_segs(Segs,NumSegs,comp,Alen,Blen); /* undo comp and rc changes */

  return (Descriptor);
}
