#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include "bri.h"

//  This is Liliana Florea's sequencing error simulator.  Bri hacked
//  it to use a real RNG, and to make it work from leaff.

typedef struct edit_script {
  int    optype;
  int    num;
  struct edit_script *next;
} EditScript_t;

typedef struct align {
  int offset, len;
  EditScript_t *script;
} Align_t;



//  This guy is provided by leaff
extern mt_s *mtctx;

//  RAND returns x numbers, starting at number y.
//
#define RAND(x,y) (int)((y) + (mtRandom32(mtctx) % (x)))

#define max(x,y) ((x)>=(y) ? (x):(y))
#define min(x,y) ((x)<=(y) ? (x):(y))

#define MOV	3
#define SUB	2
#define INS	1
#define DEL	0



EditScript_t *
new_script(int optype, int num, EditScript_t *next) {
  EditScript_t *newtp = (EditScript_t *)malloc(sizeof(EditScript_t));

  newtp->optype = optype;
  newtp->num    = num;
  newtp->next   = next;

  return newtp;
}



/* DEL(pos), SUB(pos) - modifY position pos; INS - insert right before pos */
void
insert(Align_t *aln, int in_pos, int in_optype) {
  int i, num, optype;
  EditScript_t *t, *tp;

  //fprintf(stderr, "Modify script op=%d pos=%d\n", in_optype, in_pos);

  for (t=aln->script, i=0, tp=NULL; t; tp=t, t=t->next) {
    num    = t->num;
    optype = t->optype;

    switch (optype) {
      case INS:
        if (in_pos==i+1) {
          if (tp) 
            tp->next = new_script(in_optype, 1, tp->next);
          else 
            aln->script = new_script(in_optype, 1, aln->script);
          return;
        }
        break;

      case DEL:
        i += num;
        break;

      case SUB:
            
      case MOV:
        if (i<in_pos && in_pos<=i+num) {
          int l = (in_optype==INS) ? (in_pos-i) : (in_pos-i-1);
          int r = (in_optype==INS) ? (num-l) : (num-l-1);
          if (l && l!=num) {
            t->num = l; tp = t;
            tp->next = new_script(in_optype, 1, tp->next);
            tp = tp->next;
            tp->next = new_script(optype, r, tp->next);
          } else if (!l) {
            if (tp)
              tp->next = new_script(in_optype, 1, t);
            else 
              aln->script = new_script(in_optype, 1, aln->script);
            if (in_optype!=INS) t->num -= 1;
          } else {
            tp = t;
            tp->next = new_script(in_optype, 1, tp->next);
            if (in_optype!=INS) t->num -= 1;
          }
          return;
        }
        i += num;
        break;

      default:
        fprintf(stderr, "Unrecognized optype (%d).\n", in_optype);
        break;
    }
  }

  //fprintf(stderr, "Failed to modify sequence (%d,%d).\n", in_optype, in_pos);
}


void
print_simseq(char *seq, char *hdr, Align_t *aln, double P, int CUT, int COPY) {
  int   k, e;
  char *s;
  char  let_4[4]  = {'A','C','G','T'};
  char  let_3A[3] = {'C','G','T'};
  char  let_3C[3] = {'A','G','T'};
  char  let_3G[3] = {'A','C','T'};
  char  let_3T[3] = {'A','C','G'};
  EditScript_t *t;

  while ((*hdr) && !isspace(*hdr))
    fprintf(stdout, "%c", *hdr++);

  fprintf(stdout, ":seq=%d:copy=%d %d-%d %1.2f\n", CUT+1, COPY+1, aln->offset, aln->offset+aln->len-1, P);

  s = seq + aln->offset-1;

  for (t=aln->script; t; t=t->next) {
    switch (t->optype) {
      case INS:
        for (k=0; k<t->num; k++) {
          e = RAND(4,0);
          fprintf(stdout, "%c", let_4[e]);
        }
        break;

      case DEL:
        s += t->num;
        break;

      case SUB:
        for (k=0; k<t->num; k++) {
          e = RAND(3,0);
          if (*s=='A') fprintf(stdout, "%c", let_3A[e]);
          else if (*s=='C') fprintf(stdout, "%c", let_3C[e]);
          else if (*s=='G') fprintf(stdout, "%c", let_3G[e]);
          else if (*s=='T') fprintf(stdout, "%c", let_3T[e]);
          else fprintf(stdout, "%c", 'A');
          s++;
        }
        break;
                     
      case MOV:
        for (k=0; k<t->num; k++) {
          fprintf(stdout, "%c", *s); s++;
        }
        break;

      default:
        fprintf(stderr, "Unrecognized optype (%d).\n", t->optype);
        break;
    } 
  }
  fprintf(stdout, "\n");
}    




void
simseq(char *seq, char *hdr, int len, int N, int L, int C, double P) {
  Align_t       align;
  int           i, j, k;
  int           start;
  EditScript_t *s;

  for (i=0; i<N; i++) {
    /* generate a new sequence of length min(len,N) */
    start = RAND((len-L+1),1);

    /* now create in_C non-identical copies */ 
    for (j=0; j<C; j++) {
      /* generate a 'trivial' script for the sequence */

      align.offset = start;
      align.len    = L;
      align.script = new_script(MOV,L,NULL);

      for (k=0; k<L*P; k++) {
        int optype = RAND(3,0);
        int pos    = RAND(L,1);

        insert(&align, pos, optype);
      }

      print_simseq(seq, hdr, &align, P, i, j);

      while (align.script) {
        s = align.script;
        align.script = s->next;
        free(s);
      }
    }
  }
}
