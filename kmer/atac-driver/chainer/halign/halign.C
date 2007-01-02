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

#include "halign.H"

#define min(x,y)      ((x)<=(y) ? (x):(y))
#define max(x,y)      ((x)>=(y) ? (x):(y))

typedef enum { DEL, INS, SUB, MOV } OpType;

#define START           1  // The position of the first character.
#define GAP_OPEN        1  // These are costs.
#define GAP_EXTEND      1
#define MISMATCH        1

class EditScript {
public:
  EditScript(OpType op, int nm, EditScript *nx) {
    op_type = op;
    num     = nm;
    next    = nx;
  };

  OpType        op_type;     // SUB, MOV, INS, or DEL
  int           num;         // Number of operations
  EditScript *next;
};


void
convertScriptToAlignment(EditScript    *head,
                         H_Alignment_t *aln) {

  EditScript *tp  = head;
  EditScript *tp1;
  int         scriptLen = 3;


  //  Three sets of comments, from three different sources:
  //
  // Condense_script - merge contiguous operations of the same type
  // together.  Remove the leftmost dummy script (Ascript[0][0])
  //
  // Condense the script in block operations; this will modify the
  // script and will cut off the dummy edit op from the linked list
  //
  // Condense any repeat operations.  For example, condense this
  // script: "DEL 5; DEL 3" into this array representation: "DEL 8".

  while (tp && tp->next) {
    scriptLen += 2;

    while (((tp1 = tp->next) != NULL) &&
           (tp->op_type == tp1->op_type)) {
      tp->num  += tp1->num;
      tp->next  = tp1->next;
      free(tp1);
    }

    tp = tp->next;
  }

  //  Allocate space for the alignment
  //
  if (aln->scriptAsArrayMax <= scriptLen) {
    free(aln->scriptAsArray);

    aln->scriptAsArrayMax = scriptLen;
    aln->scriptAsArray    = (int *)malloc(scriptLen * sizeof(int));

    assert(aln->scriptAsArray != NULL); 
  }

  aln->scriptAsArray[0] = 0;

  //  Convert.
  //
  int          arySize  = 0;
  EditScript  *tpdel    = 0L;

  tp = head;
  while (tp != NULL) {
    aln->scriptAsArray[++arySize] = tp->op_type;
    aln->scriptAsArray[++arySize] = tp->num;

    tpdel   = tp;
    tp      = tp->next;
    delete tpdel;
  }

  aln->scriptAsArray[0] = arySize;
}



static
int
diff(const char *s1,
     const char *s2,
     int len1,
     int len2,
     int *CC,
     int *DD,
     int *RR,
     int *SS,
     int g,
     int h,
     int x,
     int start_cgap,
     int end_cgap,
     int free_start,
     int free_end,
     int tb,
     int te,
     EditScript **head,
     EditScript **tail) {

  int i, j, s, t, c, e, tmp;
  const char *a, *b;
  int mincost, mintype, midi, midj;
  EditScript *tmp_head = 0L;
  EditScript *tmp_tail = 0L;

  if (len1==0 && len2==0) {

    *head = *tail = NULL;
    return 0;

  } else if (len2 == 0) {
    int tmpb, tmpe;

    *head = *tail = new EditScript(DEL,len1,NULL);
    tmpb = (len1 <= start_cgap) ? 0 : tb+h*(len1-start_cgap);
    tmpe = (len1 <= end_cgap)   ? 0 : te+h*(len1-end_cgap);

    return min(tmpb,tmpe);

  } else if (len1 == 0) {

    *head = *tail = new EditScript(INS,len2,NULL);
    return ((free_start || free_end) ? 0 : g+(len2*h));

  } else if (len1 == 1) {
    int tmpcost;
    char ch;

    /* insert B, delete A; or delete A, insert B */
    mincost = (start_cgap ? 0 : min(tb,te)+h) + ((free_start || free_end) ? 0 : g+len2*h);
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
        *tail = new EditScript(DEL,1,NULL);
        *head = new EditScript(INS,len2,*tail);
      } else {
        *tail = new EditScript(INS,len2,NULL);
        *head = new EditScript(DEL,1,*tail);
      }
    } else {  /* substitute A */
      EditScript *aux;
      *tail = (midj<len2-1) ? new EditScript(INS,len2-midj-1,NULL) : NULL;
      aux = new EditScript((ch==s2[midj] ? MOV : SUB),1,*tail);
      if (*tail==NULL) *tail = aux;
      *head = (midj>0) ? new EditScript(INS,midj,aux) : aux;
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

    cost1 = diff(s1,
                 s2,
                 midi,
                 midj,
                 CC,
                 DD,
                 RR,
                 SS,
                 g,
                 h,
                 x,
                 min(start_cgap,midi),
                 max(end_cgap-len1+midi,0),
                 free_start,
                 0,
                 tb,
                 g,
                 head,
                 &tmp_tail);
    cost2 = diff(s1+midi,
                 s2+midj,
                 len1-midi,
                 len2-midj,
                 CC,
                 DD,
                 RR,
                 SS,
                 g,
                 h,
                 x,
                 max(0,start_cgap-midi),
                 min(end_cgap,len1-midi),
                 0,
                 free_end,
                 g,
                 te,
                 &tmp_head,
                 tail);

    if (*head)
      tmp_tail->next = tmp_head;
    else
      *head = tmp_head;

    assert(NULL != *tail);

  } else {
    EditScript *aux;
    int cost1, cost2;

    cost1 = diff(s1,
                 s2,
                 midi-1,
                 midj,
                 CC,
                 DD,
                 RR,
                 SS,
                 g,
                 h,
                 x,
                 min(start_cgap,midi-1),
                 max(end_cgap-len1+midi-1,0),
                 free_start,
                 0,
                 tb,
                 0,
                 head,
                 &tmp_tail);

    aux = new EditScript(DEL,2,NULL);

    if (*head)
      tmp_tail->next = aux;
    else
      tmp_tail = *head = aux;

    cost2 = diff(s1+midi+1,
           s2+midj,
           len1-midi-1,
           len2-midj,
           CC,
           DD,
           RR,
           SS,
           g,
           h,
           x,
           max(0,start_cgap-midi-1),
           min(end_cgap,len1-midi-1),
           0,
           free_end,
           0,
           te,
           &tmp_head,
           tail);
    aux->next = tmp_head;

    if (*tail==NULL)
      *tail = aux;
  }

  return mincost;
}



void
halignStart(char  *s1,
            char  *s2, 
            H_Alignment_t *alignment) {

  int const offset1 = 0; // Sequence coordinates are base-based, starting from 0
  int const offset2 = 0; // but start from 1 in Liliana's code.

  if ((s1[0] == 0) || (s2[0] == 0))
    return;

  int len1 = strlen(s1);
  int len2 = strlen(s2);

  int start_cgap = 0;
  int end_cgap   = 0;
  int free_start = 0;
  int free_end   = 0;
  int score      = 0;

  EditScript *Script_head=NULL;
  EditScript *Script_tail=NULL;

  int *CC = (int *)malloc(4 * (len2+1) * sizeof(int));
  assert(NULL != CC);

  score = diff(s1,
               s2,
               len1,
               len2,
               CC,
               CC+1*(len2+1),
               CC+2*(len2+1),
               CC+3*(len2+1),
               GAP_OPEN,
               GAP_EXTEND,
               MISMATCH,
               start_cgap,
               end_cgap,
               free_start,
               free_end,
               GAP_OPEN,
               GAP_OPEN,
               &Script_head,
               &Script_tail);

  free(CC);

  assert(NULL != Script_head);
  assert(NULL != Script_tail);

  Script_tail->next = NULL;

  convertScriptToAlignment(Script_head, alignment);

  alignment->offset1 = offset1+START; // Convert from zero to one start sequence.
  alignment->offset2 = offset2+START; // Convert from zero to one start sequence.
  alignment->len1    = len1;
  alignment->len2    = len2;
  alignment->score   = score;
  alignment->first   = 1;
}


int
iterateUngappedAlignSharpEnds(H_Alignment_t *aln,
                              int &bgn1,
                              int &bgn2,
                              int &len1,
                              int &len2,
                              int &nmatInSeg) {

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
        ++lastS;
        i += *(lastS);
        j += *(lastS);
        ++lastS;
        break;
      
      case MOV:
        nmat = *(lastS+1);
        nmatInSeg ++;
        b1 = i;
        b2 = j;
        lastS++;
        i += *lastS;
        j += *lastS;
        l1 = i-b1;
        l2 = j-b2;
        lastS++;
      
        while (lastS<=endS && (*lastS==SUB || *lastS==MOV)) {
          nmat      += (*lastS==MOV) ? *(lastS+1) : 0;
          nmatInSeg += (*lastS==MOV) ? *(lastS+1) : 0;
          lastS++;
          i += *lastS;
          j += *lastS;
          if (*(lastS-1) == MOV) {
            l1 = i-b1;
            l2 = j-b2;
          }
          lastS++;
        }

        bgn1=b1-START;
        bgn2=b2-START;
        len1=l1;
        len2=l2;

        return 1; // valid output
        break;
      
      case INS:
        j += *(++lastS);
        ++lastS;
        break;
      
      case DEL:
        i += *(++lastS);
        ++lastS;
        break;

      default :
        fprintf(stderr, "Unrecognized opcode in alignment.\n");
        exit(1);
        break;
    }
  }
  return 0; // not valid output
}
