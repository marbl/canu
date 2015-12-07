
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    kmer/leaff/simseq.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2004-JUN-24 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2008-SEP-22 to 2009-JUN-13
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "mt19937ar.H"

//  This is Liliana Florea's sequencing error simulator.  Bri hacked
//  it to use a real RNG, and to make it work from leaff.

typedef struct edit_script {
  uint32  optype;
  uint32  num;
  struct edit_script *next;
} EditScript_t;

typedef struct align {
  uint32 offset, len;
  EditScript_t *script;
} Align_t;



//  This guy is provided by leaff
extern mtRandom MT;

//  RAND returns x numbers, starting at number y.
//
#define RAND(x,y) (int)((y) + (MT.mtRandom32() % (x)))

#define max(x,y) ((x)>=(y) ? (x):(y))
#define min(x,y) ((x)<=(y) ? (x):(y))

#define MOV 3
#define SUB 2
#define INS 1
#define DEL 0



EditScript_t *
new_script(int optype, uint32 num, EditScript_t *next) {
  EditScript_t *newtp = new EditScript_t;

  newtp->optype = optype;
  newtp->num    = num;
  newtp->next   = next;

  return(newtp);
}



/* DEL(pos), SUB(pos) - modifY position pos; INS - insert right before pos */
void
insert(Align_t *aln, uint32 in_pos, uint32 in_optype) {
  uint32 i, num, optype;
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
          uint32 l = (in_optype==INS) ? (in_pos-i) : (in_pos-i-1);
          uint32 r = (in_optype==INS) ? (num-l) : (num-l-1);
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
print_simseq(char *seq, char *hdr, Align_t *aln, double P, uint32 CUT, uint32 COPY) {
  uint32   k, e;
  char *s;
  char  let_4[4]  = {'A','C','G','T'};
  char  let_3A[3] = {'C','G','T'};
  char  let_3C[3] = {'A','G','T'};
  char  let_3G[3] = {'A','C','T'};
  char  let_3T[3] = {'A','C','G'};
  EditScript_t *t;

  fprintf(stdout, ">");

  while ((*hdr) && !isspace(*hdr))
    fprintf(stdout, "%c", *hdr++);

  fprintf(stdout, ":seq=%d:copy=%d:loc=%d-%d:err=%1.2f\n", CUT+1, COPY+1, aln->offset, aln->offset+aln->len-1, P);

  s = seq + aln->offset-1;

  for (t=aln->script; t; t=t->next) {
    if (*s == 0)
      break;

    switch (t->optype) {
      case INS:
        for (k=0; k<t->num; k++) {
          e = RAND(4,0);
          fprintf(stdout, "%c", let_4[e]);
        }
        break;

      case DEL:
        while (*s && t->num) {
          s++;
          t->num--;
        }
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
          if (*s == 0) {
            k = t->num;
          } else {
            fprintf(stdout, "%c", *s);
            s++;
          }
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
simseq(char *seq, char *hdr, uint32 len, uint32 N, uint32 L, uint32 C, double P) {
  Align_t       align;
  uint32        i, j, k;
  uint32        start;
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
        uint32 optype = RAND(3,0);
        uint32 pos    = RAND(L,1);

        insert(&align, pos, optype);
      }

      print_simseq(seq, hdr, &align, P, i, j);

      while (align.script) {
        s = align.script;
        align.script = s->next;
        delete s;
      }
    }
  }
}
