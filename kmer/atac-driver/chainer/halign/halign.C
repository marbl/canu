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
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "halign.h"

// from Liliana's utils.h
#define min(x,y)      ((x)<=(y) ? (x):(y))
#define max(x,y)      ((x)>=(y) ? (x):(y))

typedef 
enum { DEL, INS, SUB, MOV, NOP} OpType;

typedef 
struct EditScript_t {
        OpType op_type; /* SUB, MOV, INS, or DEL */
        int num;        /* Number of operations */
        struct EditScript_t *next;
} EditScript_t;

static void Free_script(EditScript_t *scr) {
   EditScript_t *head;

   while ((head=scr)!=NULL) {
      scr = scr->next;
      free(head);
   }
}

static int  diff(const char *,const char *,int,int,int *,int *,int *,int *,
		 int,int,int,int,int,int,int,int,int,EditScript_t **,EditScript_t **);
static EditScript_t *Condense_script(EditScript_t *,int);
static void convertScriptToAlignment (H_Alignment_t *, const EditScript_t *);
static EditScript_t *new_editscr(OpType const,int const,EditScript_t *);

#define START 1  // The position of the first character.
#define GAP_OPEN        1  // These are costs.
#define GAP_EXTEND      1
#define MISMATCH        1

////////////////////////////////////////////////////////////////

static int diff(const char *s1, const char *s2, int len1, int len2,
                int *CC, int *DD, int *RR, int *SS, int g, int h, int x,
                int start_cgap, int end_cgap, int free_start, int free_end,
                int tb, int te, EditScript_t **head, EditScript_t **tail)
{
   int i, j, s, t, c, e, tmp;
   const char *a, *b;
   int mincost, mintype, midi, midj;
   EditScript_t *tmp_head, *tmp_tail;


   if (len1==0 && len2==0) {

       *head = *tail = NULL;
       return 0;

   } else if (len2 == 0) {
       int tmpb, tmpe;

       *head = *tail = new_editscr(DEL,len1,NULL);
       tmpb = (len1<=start_cgap) ? 0 : tb+h*(len1-start_cgap);
       tmpe = (len1<=end_cgap) ? 0 : te+h*(len1-end_cgap);

       return min(tmpb,tmpe);

   } else if (len1 == 0) {

       *head = *tail = new_editscr(INS,len2,NULL);
       return ((free_start || free_end) ? 0 : g+(len2*h));

   } else if (len1 == 1) {
       int tmpcost;
       char ch;

       /* insert B, delete A; or delete A, insert B */
       mincost = (start_cgap ? 0 : min(tb,te)+h) +
                 ((free_start || free_end) ? 0 : g+len2*h);
       mintype = 2; midj = (free_start ? len2 : 0);

       /* ... or insert some B, substitute A, insert the rest of B */
       for (j=0, ch=*s1; j<len2; j++) {
           tmpcost = (ch==s2[j]) ? 0 : x;
           if (!free_start && j) tmpcost += g+j*h;
           if (!free_end && j+1<len2) tmpcost += g+(len2-j-1)*h;
           if (tmpcost<mincost) { mincost = tmpcost; mintype = 1; midj = j; }
       }
       if (mintype==2) {
           /* delete A */
           if (free_start) {
               *tail = new_editscr(DEL,1,NULL);
               *head = new_editscr(INS,len2,*tail);
           } else {
               *tail = new_editscr(INS,len2,NULL);
               *head = new_editscr(DEL,1,*tail);
           }
       } else {  /* substitute A */
           EditScript_t *aux;
           *tail = (midj<len2-1) ? new_editscr(INS,len2-midj-1,NULL) : NULL;
           aux = new_editscr((ch==s2[midj] ? MOV : SUB),1,*tail);
           if (*tail==NULL) *tail = aux;
           *head = (midj>0) ? new_editscr(INS,midj,aux) : aux;
       }
       return mincost;

   } else {
       int tmph, tmpg;

       midi = (int)(len1/2);

       /* compute CC and DD in the forward phase */
       tmph = free_start ? 0 : h;
       tmpg = free_start ? 0 : g;
       for (CC[0]=0, t=tmpg, j=1; j<=len2; j++) {
           /* if free_start, allow gap-free ends in the genomic sequence */
           CC[j] = DD[j] = t = t+tmph;
           DD[j] += tmpg;
       }

       for (a=s1, i=1; i<=midi; i++, a++) {

           s = CC[0];
           CC[0] = c = t = max(i-start_cgap,0)*h + (i>start_cgap)*tb;

           e = t + g;
           for (b=s2, j=1; j<=len2; j++, b++) {

               e = min(e, c+g) + h;

               DD[j] = (j==len2 && i>=len1-end_cgap+1) ?
                       min(DD[j], CC[j]) : (min(DD[j]+(i==start_cgap+1)*g, CC[j]+g) + h);

               c = min(DD[j], min(e, s+x*(*a!=*b)));
               s = CC[j]; CC[j] = c;
           }
       }
       DD[0] = CC[0];

       /* compute RR and SS in the reverse phase */
       tmph = free_end ? 0 : h;
       tmpg = free_end ? 0 : g;
       for (RR[len2]=0, t=tmpg, j=len2-1; j>=0; --j) {
            /* if free_end, allow gap-free ends in the genomic sequence */
            RR[j] = SS[j] = t = t+tmph;
            SS[j] += tmpg;
       }

       for (a=s1+len1-1, i=len1-1; i>=midi; --i, --a) {

           s = RR[len2];
           RR[len2] = c = t = max((len1-end_cgap)-i,0)*h + (i<len1-end_cgap)*te;

           e = t + g;
           for (b=s2+len2-1, j=len2-1; j>=0; --j, --b) {

               e = min(e, c+g) + h;

               SS[j] = (j==0 && i<start_cgap) ?
                       min(SS[j], RR[j]) : (min(SS[j]+(i==len1-end_cgap-1)*g, RR[j]+g) + h);
               c = min(SS[j], min(e, s+x*(*a!=*b)));
               s = RR[j]; RR[j] = c;


           }
       }
       SS[len2] = RR[len2];
   }

   /* find midj that minimizes the sum */
   /* special cases: columns 0 and len2 */
   midj = 0;
   if (CC[0]+RR[0]<=DD[0]+SS[0]-g*(midi>start_cgap)) {
       mincost = CC[0]+RR[0];
       mintype = 1;
   } else {
       mincost = DD[0]+SS[0]-g*(midi>start_cgap);
       mintype = 2;
   }

   for (j=1; j<len2; j++) {
       tmp = min(CC[j]+RR[j],DD[j]+SS[j]-g);
       if (mincost > tmp) {
           mincost = tmp;
           midj = j;
           mintype = (tmp == CC[j]+RR[j]) ? 1:2;
       }
   }

   tmp = min(CC[len2]+RR[len2],DD[len2]+SS[len2]-g*(midi<len1-end_cgap));
   if (mincost > tmp) {
      mincost = tmp;
      midj = len2;
      mintype = (tmp==CC[len2]+RR[len2]) ? 1:2;
   }
   /* compute recursively in the two subregions */
   if (mintype==1) {
      int cost1, cost2;

      cost1 =
      diff(s1,s2,midi,midj,CC,DD,RR,SS,g,h,x,
           min(start_cgap,midi),max(end_cgap-len1+midi,0),
           free_start,0,tb,g,head,&tmp_tail);
      cost2 =
      diff(s1+midi,s2+midj,len1-midi,len2-midj,CC,DD,RR,SS,g,h,x,
           max(0,start_cgap-midi),min(end_cgap,len1-midi),
           0,free_end,g,te,&tmp_head,tail);
      if (*head)
           tmp_tail->next = tmp_head;
      else
           *head = tmp_head;
      assert(NULL != *tail);

   } else {
      EditScript_t *aux;
      int cost1, cost2;

      cost1 =
      diff(s1,s2,midi-1,midj,CC,DD,RR,SS,g,h,x,
           min(start_cgap,midi-1),max(end_cgap-len1+midi-1,0),
           free_start,0,tb,0,head,&tmp_tail);
      aux = new_editscr(DEL,2,NULL);
      if (*head)
           tmp_tail->next = aux;
      else
           tmp_tail = *head = aux;
      cost2 =
      diff(s1+midi+1,s2+midj,len1-midi-1,len2-midj,CC,DD,RR,SS,g,h,x,
           max(0,start_cgap-midi-1), min(end_cgap,len1-midi-1),
           0,free_end,0,te,&tmp_head,tail);
      aux->next = tmp_head;
      if (*tail==NULL) *tail = aux;
   }

   return mincost;
}


static EditScript_t *new_editscr(OpType const op_type, int const num, EditScript_t * next) {
  EditScript_t *escr = (EditScript_t *) malloc(sizeof(EditScript_t));
  assert(NULL != escr);
  
  escr->op_type = op_type;
  escr->num     = num;
  escr->next    = next;
  
  return escr;
}






/* Condense_script - merge contiguous operations of the same type together, but do not
   free memory; remove the leftmost dummy script (Ascript[0][0]) */
static EditScript_t *Condense_script(EditScript_t *head, int free_flag)
{
   EditScript_t *tp, *tp1;

   tp = head;
   while (tp && tp->next) {
      while (((tp1 = tp->next) != NULL) && (tp->op_type == tp1->op_type)) {
            tp->num = tp->num + tp1->num;
            tp->next = tp1->next;
            if (free_flag) free(tp1);  /* NEW */
      }
      if (tp->next->op_type==NOP) { /* bypass the 'dummy' marker */
          if (free_flag) free(tp->next);
          tp->next = NULL;
      }
      tp = tp->next;
   }

   return head; 
}

static void Flip_script(struct EditScript_t **script)
{
   struct EditScript_t *ep, *ahead, *behind;

   ahead = *script;
   ep = NULL;
   while (ahead!=NULL) {
          behind = ep;
          ep = ahead;
          ahead = ahead->next;
          ep->next = behind;
  }
  *script = ep;
}


// This is the EST-free version provided by Liliana. -- Jason
// Now it gets S from the EditScript_t instead of by argument.
void convertScriptToAlignment(H_Alignment_t* aln, const EditScript_t * tp)  {
  int arySize = 0;
  int * const aryPtr = aln->scriptAsArray;
  int * const basePtr = aryPtr + 1;
  int prev_op = NOP; 
  
  if (tp != NULL) {
    assert (tp->op_type != prev_op); // Jason says 1st one can't be NOP, right?
    while (tp != NULL) {
      if (tp->op_type == prev_op) {
	// Condense any repeat operations.
	// For example, condense this script: "DEL 5; DEL 3"
	// into this array representation: "DEL 8".
	// Assume 1st op is never NOP.
	*(basePtr+arySize-1) += tp->num;
      } else {
	*(basePtr+arySize) = tp->op_type;
	*(basePtr+arySize+1) = tp->num;
	arySize += 2;
      }
      prev_op = tp->op_type;
      tp = tp->next;
    }
  }
  (*aryPtr) = arySize;
}

////////////////////////////////////////////////////////////////

void Free_align(H_Alignment_t* aln_ptr)
{
   free(aln_ptr->scriptAsArray);
   free(aln_ptr);
}

void halignStart
(char const * const s1, 
 char const * const s2, 
 int const offset1, // Sequence coordinates are base-based, starting from 0
 int const offset2, // but start from 1 in Liliana's code.
 int const len1, 
 int const len2, 
 H_Alignment_t ** alignment
 ) {
   int *CC, *DD, *RR, *SS, g, h, x;
   int start_cgap, end_cgap, free_start, free_end, score, script_len;
   EditScript_t *Script_head=NULL, *Script_tail=NULL, *tp;
   H_Alignment_t * new_alignment;
   EditScript_t * scriptAsLinkedList;

   assert(strlen(s1)>0);
   assert(strlen(s2)>0);

   g    = GAP_OPEN;
   h    = GAP_EXTEND;
   x    = MISMATCH;

   free_start = free_end = 0;

   start_cgap = 0; // max(est->uexons[startxi].f2-est->exons[startxi].f2,0);
   end_cgap = 0; // max(est->exons[xi].t2-est->uexons[xi].t2,0);

   CC = (int *)malloc(4*(len2+1)*sizeof(int));
   assert(NULL != CC);
   DD = CC+len2+1;
   RR = DD+len2+1;
   SS = RR+len2+1;


   score = diff(s1,s2,len1,len2,CC,DD,RR,SS,g,h,x,
		start_cgap,end_cgap,free_start,free_end,g,g,&Script_head,&Script_tail);
   //fprintf(stderr,"s1=<%s>\n",s1);
   //fprintf(stderr,"s2=<%s>\n",s2);
   assert(NULL != Script_head);
   assert(NULL != Script_tail);
   Script_tail->next = new_editscr(NOP,1,NULL); /* attach a dummy marker */

   /* generate the alignment(s) */
   new_alignment = (H_Alignment_t *)malloc(sizeof(H_Alignment_t));
   assert(NULL != new_alignment);
   new_alignment->offset1 = offset1+START; // Convert from zero to one start sequence.
   new_alignment->offset2 = offset2+START; // Convert from zero to one start sequence.
   new_alignment->len1 = len1;
   new_alignment->len2 = len2;
   new_alignment->score = score;
   new_alignment->first = 1;
   //   new_alignment->next_script = NULL;

   /* Condense the script in block operations; this will modify the script 
      and will cut off the dummy edit op from the linked list */
   scriptAsLinkedList = Condense_script(Script_head,1);
   /* not necessary to flip the script */

   for (tp=Script_head, script_len=0; tp; script_len++, tp=tp->next);
   new_alignment->scriptAsArray
     = (int *)malloc((2*script_len+1)*sizeof(int));
   assert(NULL != new_alignment->scriptAsArray); 

   //   S2A(new_alignment,est,startxi,xi);
   convertScriptToAlignment (new_alignment, scriptAsLinkedList);

   Free_script(scriptAsLinkedList);  /* NEW */
   scriptAsLinkedList = NULL;

   if(*alignment != NULL){
     Free_align(*alignment);
   }
   *alignment = new_alignment;

   free(CC); // free(DD); free(RR); free(SS);
}

void printUngappedAlign(H_Alignment_t* aln)
{
  int *lastS, *endS;
  int i, j, b1, l1, b2, l2;
  int nmat = 0;

  i = aln->offset1;
  j = aln->offset2;

  lastS = aln->scriptAsArray +1; 
  endS  = aln->scriptAsArray + aln->scriptAsArray[0];

  while (lastS <= endS) {
     switch (*lastS) {
        case SUB:  

        case MOV:  nmat = (*lastS==MOV) ? *(lastS+1) : 0;
                   b1 = i; b2 = j; lastS++; i += *lastS; j += *lastS;
                   l1 = i-b1; l2 = j-b2;
                   lastS++;
               
                   while (lastS<=endS && (*lastS==SUB || *lastS==MOV)) {
                      nmat += (*lastS==MOV) ? *(lastS+1) : 0;
                      lastS++;
                      i += *lastS; j += *lastS;
                      l1 = i-b1; l2 = j-b2;
                      lastS++; 
                   }
                   // nmat is the number of matches.
                   // These lines should become a list of gap filling matches.
                   printf("%d %d %d %d %d\n", b1-START, b2-START, l1, l2, nmat);  // or SAVE

                   break;

        case INS:  j += *(++lastS); ++lastS;

                   break;
 
        case DEL:  i += *(++lastS); ++lastS;

                   break;

        default :  fprintf(stderr, "Unrecognized opcode in alignment.\n");
                   exit(1);

     }
  }

  return;
}

void printUngappedAlignSharpEnds(H_Alignment_t* aln)
{
  int *lastS, *endS;
  int i, j, b1, l1, b2, l2;
  int nmat = 0;

  i = aln->offset1;
  j = aln->offset2;

  lastS = aln->scriptAsArray +1;
  endS  = aln->scriptAsArray + aln->scriptAsArray[0];

  while (lastS <= endS) {
     switch (*lastS) {
        case SUB:
          ++lastS; i += *(lastS); j += *(lastS); ++lastS;
          break;

        case MOV:
          nmat = *(lastS+1);
          b1 = i; b2 = j; lastS++; i += *lastS; j += *lastS;
          l1 = i-b1; l2 = j-b2;
          lastS++;
          
          while (lastS<=endS && (*lastS==SUB || *lastS==MOV)) {
            nmat += (*lastS==MOV) ? *(lastS+1) : 0;
            lastS++;
            i += *lastS; j += *lastS;
            if (*(lastS-1) == MOV) {
              l1 = i-b1; l2 = j-b2;
            }
            lastS++;
          }
          printf("%d %d %d %d %d\n", b1-START, b2-START, l1, l2, nmat);  // or SAVE
          break;
          
        case INS:
          j += *(++lastS); ++lastS;
          break;

        case DEL:
          i += *(++lastS); ++lastS;
          break;

        default :
          fprintf(stderr, "Unrecognized opcode in alignment.\n");
          exit(1);

     }
  }

  return;
}

void printUngappedAlignSharpEndsOnConsole
(H_Alignment_t* aln,
 char const * const seq1, char const * const seq2, int const opt)
{
  int *lastS, *endS;
  int b1, l1, b2, l2;
  int nmat = 0;
  int lastc=0;
  int seg = 1;

  int i = aln->offset1;
  int j = aln->offset2;
  int const len1 = aln->len1;
  int const len2 = aln->len2;

  lastS = aln->scriptAsArray +1;
  endS  = aln->scriptAsArray + aln->scriptAsArray[0];

  if(opt == 1) { printf("%.*s\n", len1, seq1);}
  if(opt == 2) { printf("%.*s\n", len2, seq2);}

  while (lastS <= endS) {
     switch (*lastS) {
        case SUB:
          ++lastS; i += *(lastS); j += *(lastS); ++lastS;
          break;

        case MOV:
          nmat = *(lastS+1);
          b1 = i; b2 = j; lastS++; i += *lastS; j += *lastS;
          l1 = i-b1; l2 = j-b2;
          lastS++;
          
          while (lastS<=endS && (*lastS==SUB || *lastS==MOV)) {
            nmat += (*lastS==MOV) ? *(lastS+1) : 0;
            lastS++;
            i += *lastS; j += *lastS;
            if (*(lastS-1) == MOV) {
              l1 = i-b1; l2 = j-b2;
            }
            lastS++;
          }
          
          if(opt == 0) {
            printf("%d %d %d %d %d\n", b1-START, b2-START, l1, l2, nmat);  // or SAVE
          } else {
            int bgn=0;
            switch(opt){
            case 1:
              { bgn = b1-1;}
              break;
            case 2:
              { bgn = b2-1;}
              break;
            default:
              assert(0);
            }
            { int ic; int len=bgn-lastc; for(ic=0;ic<len;ic++){printf("-");};lastc += len;}
            { int ic; int len=l1; for(ic=0;ic<len;ic++){printf("%d",(seg % 10));};lastc += len;seg +=1;}
          }
          break;
          
     case INS:  j += *(++lastS); ++lastS;
       break;
       
     case DEL:  i += *(++lastS); ++lastS;
       break;
       
     default :  fprintf(stderr, "Unrecognized opcode in alignment.\n");
       exit(1);
       
     }
  }
  if(opt != 0) { printf("\n");}
  
  return;
}

int iterateUngappedAlignSharpEnds
( H_Alignment_t* aln,
  int& bgn1,
  int& bgn2,
  int& len1,
  int& len2,
  int& nmatInSeg)
{
  // Returns zero when exhasted.
  // Returns one when the args are valid.

  static int *lastS, *endS;
  static int i, j;
  
  nmatInSeg = 0;

  if(aln == NULL) return 0; // not valid output

  if(aln->first){
    aln->first = 0;
    i = aln->offset1;
    j = aln->offset2;
    
    lastS = aln->scriptAsArray + 1;
    endS  = aln->scriptAsArray + aln->scriptAsArray[0];
  }

  while (lastS <= endS) {
    int b1, l1, b2, l2;
    int nmat;
    switch (*lastS) {
    case SUB:
      ++lastS; i += *(lastS); j += *(lastS); ++lastS;
      break;
      
    case MOV:
      nmat = *(lastS+1);
      nmatInSeg ++;
      b1 = i; b2 = j; lastS++; i += *lastS; j += *lastS;
      l1 = i-b1; l2 = j-b2;
      lastS++;
      
      while (lastS<=endS && (*lastS==SUB || *lastS==MOV)) {
	nmat += (*lastS==MOV) ? *(lastS+1) : 0;
	nmatInSeg += (*lastS==MOV) ? *(lastS+1) : 0;
	lastS++;
	i += *lastS; j += *lastS;
	if (*(lastS-1) == MOV) {
	  l1 = i-b1; l2 = j-b2;
	}
	lastS++;
      }
      //printf("%d %d %d %d %d\n", b1-START, b2-START, l1, l2, nmatInSeg);  // or SAVE
      bgn1=b1-START; bgn2=b2-START; len1=l1; len2=l2;
      //return (lastS <= endS); // again
      return 1; // valid output
      break;
      
    case INS:
      j += *(++lastS); ++lastS;
      break;
      
    case DEL:
      i += *(++lastS); ++lastS;
      break;
      
    default :
      fprintf(stderr, "Unrecognized opcode in alignment.\n");
      exit(1);
      
    }
  }
  return 0; // not valid output
}

#if 0
#include <string>
#include <time>

int main(int argc, char *argv[])
{
   char *seq1, *seq2;
   int   len1, len2;
   int   offset1, offset2;
   H_Alignment_t* aln_ptr;
   // Sequence coordinates are base-based, starting from 0
   halign(seq1+offset1, // This is the first base in the comparison.
          seq2+offset2,
          offset1, offset2,
          len1, len2,
          &aln_ptr);
   
   printUngappedAlign(aln_ptr);
   printUngappedAlignSharpEnds(aln_ptr);

   printUngappedAlignSharpEndsOnConsole(aln_ptr, seq1, seq2, 0);
   printUngappedAlignSharpEndsOnConsole(aln_ptr, seq1, seq2, 1);
   printUngappedAlignSharpEndsOnConsole(aln_ptr, seq1, seq2, 2);
   Free_align(aln_ptr); // Must call for each halign() but after printing output.

   exit(0);
}

#endif
