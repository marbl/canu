#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bri++.H"
#include "fasta.H"
#include "fasta-cache.H"
#include "sim4polish.h"

//  This code takes basic sim4db format polishes and recomputes the
//  alignments and scores.  Required in the input polishes are the EST
//  id, genomic id, exon coordinates and an orientation.

void
align(const char *string1,
      const char *string2,
      const int len1,
      const int len2,
      char *alnline1,
      char *alnline2);

int
main(int argc, char **argv) {

  //  Load all the sequences.  We really do need them all the ESTs in
  //  core, since they probably aren't in a useful sorted order.  You
  //  can probably figure out a way to get rid of the FastACache for
  //  the GEN.  Doing so will reduce memory usage by about 50%.

  FastACache *EST = 0L;
  FastACache *GEN = 0L;
  int         mergeTolerancePerc = 0;
  int         mergeToleranceBase = 0;
  int         statsOnly          = 0;

  //  Statistics on the exon merge

  int   mergedExons   = 0;
  int   mergedMatches = 0;

  int   numcdnagaps        = 0;
  int   nummatcheswithgaps = 0;

  FILE *mergeLog      = 0L;

  int     arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-merge", 2) == 0) {
      mergeTolerancePerc = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-b", 2) == 0) {
      mergeToleranceBase = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-M", 2) == 0) {
      mergeLog = fopen(argv[++arg], "w");
    } else if (strncmp(argv[arg], "-e", 2) == 0) {
      if (statsOnly)
        EST = new FastACache(argv[++arg], 1000, false, true);  //  debugging only!
      else 
        EST = new FastACache(argv[++arg],    0, true);
    } else if (strncmp(argv[arg], "-g", 2) == 0) {
      GEN = new FastACache(argv[++arg],    1, false, true);
    } else if (strncmp(argv[arg], "-q", 2) == 0) {
      statsOnly = 1;
    }
    arg++;
  }

  if ((statsOnly == 0) && (!EST || !GEN)) {
    fprintf(stderr, "usage: %s [-merge percent-tolerance] [-M merge-log] [-q] -e est.fasta -g genome.fasta < polishes > somewhere\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "       Polishes _MUST_ be sorted by genomic index.\n");
    fprintf(stderr, "       If not, performance will be worse than atrocious.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       percent-tolerance -- merge exons separated by gap if\n");
    fprintf(stderr, "       the cDNA and genomic gaps differ by less than p percent.\n");
    fprintf(stderr, "       A value of 5 means 5%\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -q: Don't actually do the work, just count the statistics\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    exit(1);
  }


  char   *s1 = new char [16 * 1024 * 1024];
  char   *s2 = new char [16 * 1024 * 1024];
  int     l1 = 0;
  int     l2 = 0;

  speedCounter  *C = new speedCounter("%12.0f polishes -- %12.0f polishes/second\r",
                                      1.0, 0xf, true);

  sim4polish *p = 0L;
  while ((p = s4p_readPolish(stdin)) != 0L) {

    //  If we have a mergeTolerance, merge adjacent exons that are
    //  separated my approximately equal sized cDNA and genomic gaps.
    //
    //  Possible a better way to do this is to check if the identity
    //  of the missing region is decent, too.


    //  Remember the id/cv of this guy for the log
    //
    double id = 0.0;
    double cv = 0.0;
    if (mergeLog) {
      id = s4p_percentIdentity(p);
      cv = s4p_percentCoverage(p);
    }

    int   merged = 0;
    int   gapped = 0;

    if ((mergeTolerancePerc > 0) || (mergeToleranceBase > 0)) {
      for (u32bit i=1; i<p->numExons; i++) {
        int cgap = p->exons[i].estFrom - p->exons[i-1].estTo;
        int ggap = p->exons[i].genFrom - p->exons[i-1].genTo;

        bool  mergeGap = false;

        //  New method -- check if the gaps are within 20bp of each other
        //
        int diff = cgap - ggap;
        if (diff < 0)
          diff = -diff;

        if (diff < mergeToleranceBase)
          mergeGap = true;


        //  Original method -- cehck if the gaps are within 10% of each other
        //
        int ctol = cgap * (100 + mergeTolerancePerc);
        int gtol = ggap * (100 + mergeTolerancePerc);

        cgap *= 100;
        ggap *= 100;

        if (((cgap < ggap) && (ctol > ggap)) ||
            ((ggap < cgap) && (gtol > cgap)))
          mergeGap = true;

        if (cgap > 1) {
          numcdnagaps++;
          gapped++;
        }

        if ((cgap > 1) && (mergeGap)) {

          //  Merge i and i-1 if adding in the tolerance makes either
          //  the cgap or the ggap longer than the other gap.  i.e., the
          //  cgap was shorter, but including the tolerance makes it
          //  longer, so they're about the same size.

          if (mergeLog)
            fprintf(mergeLog, "MERGE: %4d-%4d (%6.2f,%6.2f) %4d-%4d and %8d-%8d (%6.2f,%6.2f) %8d-%8d\n",
                    p->exons[i-1].estFrom, p->exons[i-1].estTo,
                    cgap / 100.0, ctol / 100.0,
                    p->exons[i].estFrom, p->exons[i].estTo,
                    p->exons[i-1].genFrom, p->exons[i-1].genTo,
                    ggap / 100.0, gtol / 100.0,
                    p->exons[i].genFrom, p->exons[i].genTo);

          //  merge exons
          p->exons[i-1].estTo = p->exons[i].estTo;
          p->exons[i-1].genTo = p->exons[i].genTo;

          //  delete this exon
          s4p_deleteExon(p, i);

          //  Do it again!
          i--;

          merged++;
          mergedExons++;
        }
      }

      if (merged)
        mergedMatches++;
      if (gapped)
        nummatcheswithgaps++;
    }


    //  For each exon, generate an alignment


    if (statsOnly == 0) {
      p->estLen   = EST->getSequence(p->estID)->sequenceLength();
      p->estPolyA = 0;
      p->estPolyT = 0;

      for (u32bit i=0; i<p->numExons; i++) {
        l1 = p->exons[i].estTo - p->exons[i].estFrom;
        l2 = p->exons[i].genTo - p->exons[i].genFrom;

        strncpy(s1, EST->getSequence(p->estID)->sequence() + p->exons[i].estFrom, l1);
        strncpy(s2, GEN->getSequence(p->genID)->sequence() + p->exons[i].genFrom, l2);

        if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
          strncpy(s1, EST->getSequence(p->estID)->sequence() + p->estLen - p->exons[i].estTo, l1);
          reverseComplementSequence(s1, l1);
        }

        p->exons[i].estAlignment = (char *)malloc(sizeof(char) * (l1+l2+1));
        p->exons[i].genAlignment = (char *)malloc(sizeof(char) * (l1+l2+1));

        align(s1, s2,
              l1, l2,
              p->exons[i].estAlignment,
              p->exons[i].genAlignment);
      }

      //  There isn't an intron after the last exon.  Force it.
      //
      p->exons[p->numExons-1].intronOrientation = SIM4_INTRON_NONE;

      s4p_updateAlignmentScores(p);
      s4p_printPolish(stdout, p, 0);
    }

    if (merged) {
      fprintf(mergeLog, "MERGED\tEST\t%d\tfrom\t%8.3f\t%8.3f\tto\t%8.3f\t%8.3f\n",
              p->estID, id, cv, s4p_percentIdentity(p), s4p_percentCoverage(p));
    }

    s4p_destroyPolish(p);

    C->tick();
  }

  if ((mergeTolerancePerc > 0) || (mergeToleranceBase > 0)) {
    fprintf(stderr, "FOUND:  %d gaps in %d matches.\n", numcdnagaps, nummatcheswithgaps);
    fprintf(stderr, "MERGED: %d gaps in %d matches.\n", mergedExons, mergedMatches);
  }

  delete GEN;
  delete EST;
}


////////////////////////////////////////
//
//  The rest of this file is Liliana Florea's halign.
//
////////////////////////////////////////

#define MININT     (-99999)

#define DEL 0
#define INS 1
#define SUB 2

#define min(x,y)      ((x)<=(y) ? (x):(y))
#define max(x,y)      ((x)>=(y) ? (x):(y))

typedef struct edit_script {
  int  op_type;   /* SUB, INS or DEL */
  int  num;        /* Number of operations */
  struct edit_script *next;
} edit_script;

typedef struct edit_script_list {
  int          offset1, offset2;
  int          len1, len2;
  int          score;
  int          first;
  edit_script *script;
} edit_script_list;

static void sim4_align_path(const char *,const char *,int,int,int,int,int,edit_script**,edit_script**,int,int);
static int  sim4_align_get_dist(const char *,const char *,int, int, int, int, int);

static void Free_script(edit_script *);

static int snake(const char *,const char *,int, int, int, int);
static int rsnake(const char *,const char *,int, int, int, int, int,int);

void
align(const char *seq1,
      const char *seq2,
      const int   len1,
      const int   len2,
      char *alnline1,
      char *alnline2) {
  edit_script      *head, *tail, *tp;
  int               diff, i;
  char             *a, *b, *at, *bt;

  diff = sim4_align_get_dist(seq1, seq2, 0, 0, len1, len2, len1+len2);

  sim4_align_path(seq1, seq2, 0, 0, len1, len2, diff, &head, &tail, len1, len2);
        
  /* generate the alignment(s) */
  alnline1[0] = '\0';
  alnline2[0] = '\0';

  a = (char *)seq1; at = alnline1;
  b = (char *)seq2; bt = alnline2;

  for (tp=head; tp; tp=tp->next) {
    switch (tp->op_type) {
      case SUB:
        for (i=0; i<tp->num; i++) {
          if (*a==*b) {
            *at = tolower(*a);
            *bt = tolower(*b);
          } else {
            *at = toupper(*a);
            *bt = toupper(*b);
          }
          a++; b++; at++; bt++;
        }
        break; 

      case INS:
        for (i=0; i<tp->num; i++) {
          *at = '-'; *bt = toupper(*b);
          b++; at++; bt++;
        }
        break;  

      case DEL:
        for (i=0; i<tp->num; i++) {
          *bt = '-'; *at = toupper(*a);
          a++; at++; bt++;
        }
        break;

      default:
        fprintf(stderr, "Unrecognized op_type in script. %d\n", tp->op_type);
        exit(0);
    }
  }
  *at = '\0';
  *bt = '\0';

  Free_script(head);
}


int
sim4_align_get_dist(const char *seq1, const char *seq2, int i1, int j1, int i2, int j2, int limit) {
  int *last_d, *temp_d;
  int goal_diag, ll, uu;
  int c, k, row;
  int start, lower, upper;

  /* Compute the boundary diagonals */
  start = j1 - i1;
  lower = max(j1-i2, start-limit);
  upper = min(j2-i1, start+limit);
  goal_diag = j2-i2;

  if (goal_diag > upper || goal_diag < lower) {
    (void)fprintf(stderr, "The two sequences are not really similar.\n");
    (void)fprintf(stderr, "Please try an exact aligning method.\n");
    exit(1);
  }

  /* Allocate space for forward vectors */
  last_d = (int *)malloc((upper-lower+1)*sizeof(int)) - lower;
  temp_d = (int *)malloc((upper-lower+1)*sizeof(int)) - lower;

  /* Initialization */
  for (k=lower; k<=upper; ++k) last_d[k] = MININT;
  last_d[start] = snake(seq1, seq2, start, i1, i2, j2);

  if (last_d[goal_diag] >= i2) {
    /* Free working vectors */
    free(last_d+lower);
    free(temp_d+lower);
    return 0;
  }

  for (c=1; c<=limit; ++c) {
    ll = max(lower,start-c); uu = min(upper, start+c);
    for (k=ll; k<=uu; ++k) {
      if (k == ll)
        row = last_d[k+1]+1;    /* DELETE */
      else if (k == uu)
        row = last_d[k-1];      /* INSERT */
      else if ((last_d[k]>=last_d[k+1]) &&
               (last_d[k]+1>=last_d[k-1]))
        row = last_d[k]+1;      /*SUBSTITUTE */
      else if ((last_d[k+1]+1>=last_d[k-1]) &&
               (last_d[k+1]>=last_d[k]))
        row = last_d[k+1]+1;    /* DELETE */
      else
        row = last_d[k-1];      /* INSERT */

      temp_d[k] = snake(seq1,seq2,k,row,i2,j2);
    }

    for (k=ll; k<=uu; ++k) last_d[k] = temp_d[k];

    if (last_d[goal_diag] >= i2) {
#ifdef STATS
      (void)fprintf(stderr, "get_dist = %d\n",c);
#endif

      /* Free working vectors */
      free(last_d+lower);
      free(temp_d+lower);
      return c;
    }
  }

  /* Ran out of distance limit */
  return -1;
}

void
sim4_align_path(const char *seq1, const char *seq2, int i1, int j1, int i2, int j2, int dist, edit_script **head, edit_script **tail, int M, int N) {

  int     *last_d, *temp_d,       /* forward vectors */
    *rlast_d, *rtemp_d;     /* backward vectors */

  edit_script *head1, *tail1, *head2, *tail2;
  int midc, rmidc;
  int start, lower, upper;
  int rstart, rlower, rupper;
  int c, k, row;
  int mi, mj, tmp, ll, uu;
  char flag;

  *head = *tail = NULL;

  /* Boundary cases */
  if (i1 == i2) {
    if (j1 == j2) *head = NULL;
    else {
      head1 = (edit_script *) malloc(sizeof(edit_script));
      head1->op_type = INS;
      head1->num = j2-j1;
      head1->next = NULL;
      *head = *tail = head1;
    }
    return;
  }

  if (j1 == j2) {
    head1 = (edit_script *) malloc(sizeof(edit_script));
    head1->op_type = DEL;
    head1->num = i2-i1;
    head1->next = NULL;
    *head = *tail = head1;
    return;
  }

  if (dist <= 1) {
    start = j1-i1;
    if (j2-i2 == j1-i1) {
      head1 = (edit_script *) malloc(sizeof(edit_script));
      head1->op_type = SUB;
      head1->num = i2-i1;
      head1->next = NULL;
      *head = *tail = head1;
    } else if (j2-j1 == i2-i1+1) {

      tmp = snake(seq1,seq2,start,i1,i2,j2);
      if (tmp>i1) {
        head1 = (edit_script *) malloc(sizeof(edit_script));
        head1->op_type = SUB;
        head1->num = tmp-i1;
        *head = head1;
      }
      head2 = (edit_script *) malloc(sizeof(edit_script));
      head2->op_type = INS;
      head2->num = 1;

      if (*head) head1->next = head2;
      else *head = head2;
      *tail = head2;
      head2->next = NULL;

      if (i2-tmp) {
        head1 = head2;
        *tail = head2 = (edit_script *)malloc(sizeof(edit_script));
        head2->op_type = SUB;
        head2->num = i2-tmp;
        head2->next = NULL;
        head1->next = head2;
      }
    } else if (j2-j1+1 == i2-i1) {

      tmp = snake(seq1,seq2,start,i1,i2,j2);
      if (tmp>i1) {
        head1 = (edit_script *) malloc(sizeof(edit_script));
        head1->op_type = SUB;
        head1->num = tmp-i1;
        *head = head1;
      }
      head2 = (edit_script *) malloc(sizeof(edit_script));
      head2->op_type = DEL;
      head2->num = 1;

      if (*head) head1->next = head2;
      else *head = head2;
      *tail = head2;
      head2->next = NULL;

      if (i2>tmp+1) {
        head1 = head2;
        *tail = head2 = (edit_script *)malloc(sizeof(edit_script));
        head2->op_type = SUB;
        head2->num = i2-tmp-1;
        head2->next = NULL;
        head1->next = head2;
      }
    } else {
      (void)fprintf(stderr,
                    "align.c: warning: something wrong when aligning.");
    }
    return;
  }

  /* Divide the problem at the middle cost */
  midc = dist/2;
  rmidc = dist - midc;

  /* Compute the boundary diagonals */
  start = j1 - i1;
  lower = max(j1-i2, start-midc);
  upper = min(j2-i1, start+midc);
  rstart = j2-i2;
  rlower = max(j1-i2, rstart-rmidc);
  rupper = min(j2-i1, rstart+rmidc);

  /* Allocate space for forward vectors */
  last_d = (int *)malloc((upper-lower+1)*sizeof(int)) - lower;
  temp_d = (int *)malloc((upper-lower+1)*sizeof(int)) - lower;

  for (k=lower; k<=upper; k++) last_d[k] = -1;
  last_d[start] = snake(seq1,seq2,start,i1,i2,j2);

  /* Forward computation */
  for (c=1; c<=midc; ++c) {
    ll = max(lower,start-c);
    uu = min(upper,start+c);
    for (k=ll; k<=uu; ++k) {
      if (k == ll) {
        /* DELETE : down from (k+1,c-1) */
        row = last_d[k+1]+1;
      } else if (k == uu) {
        /* INSERT : right from (k-1,c-1) */
        row = last_d[k-1];
      } else if ((last_d[k]>=last_d[k+1]) &&
                 (last_d[k]+1>=last_d[k-1])) {
        /* SUBSTITUTE */
        row = last_d[k]+1;
      } else if ((last_d[k+1]+1>=last_d[k-1]) &&
                 (last_d[k+1]>=last_d[k])) {
        /* DELETE */
        row = last_d[k+1]+1;
      } else {
        /* INSERT */
        row = last_d[k-1];
      }

      temp_d[k] = snake(seq1,seq2,k,row,i2,j2);
    }
    for (k=ll; k<=uu; ++k)
      last_d[k] = temp_d[k];
  }

  /* Allocate space for backward vectors */
  rlast_d = (int *)malloc((rupper-rlower+1)*sizeof(int)) - rlower;
  rtemp_d = (int *)malloc((rupper-rlower+1)*sizeof(int)) - rlower;

  for (k=rlower; k<=rupper; k++) rlast_d[k] = i2+1;
  rlast_d[rstart] = rsnake(seq1,seq2,rstart,i2,i1,j1,M,N);

  /* Backward computation */
  for (c=1; c<=rmidc; ++c) {
    ll = max(rlower,rstart-c);
    uu = min(rupper,rstart+c);
    for (k=ll; k<=uu; ++k) {
      if (k == ll) {
        /* INSERT : left from (k+1,c-1) */
        row = rlast_d[k+1];
      } else if (k == uu) {
        /* DELETE : up from (k-1,c-1) */
        row = rlast_d[k-1]-1;
      } else if ((rlast_d[k]-1<=rlast_d[k+1]) &&
                 (rlast_d[k]-1<=rlast_d[k-1]-1)) {
        /* SUBSTITUTE */
        row = rlast_d[k]-1;
      } else if ((rlast_d[k-1]-1<=rlast_d[k+1]) &&
                 (rlast_d[k-1]-1<=rlast_d[k]-1)) {
        /* DELETE */
        row = rlast_d[k-1]-1;
      } else {
        /* INSERT */
        row = rlast_d[k+1];
      }

      rtemp_d[k] = rsnake(seq1,seq2,k,row,i1,j1,M,N);
    }
    for (k=ll; k<=uu; ++k)
      rlast_d[k] = rtemp_d[k];
  }

  /* Find (mi, mj) such that the distance from (i1, j1) to (mi, mj) is
     midc and the distance from (mi, mj) to (i2, j2) is rmidc.
  */

  flag = false;
  mi = i1; mj = j1;
  ll = max(lower,rlower); uu = min(upper,rupper);
  for (k=ll; k<=uu; ++k) {
    if (last_d[k]>=rlast_d[k]) {
      if (last_d[k]-i1>=i2-rlast_d[k]) {
        mi = last_d[k]; mj = k+mi;
      } else {
        mi = rlast_d[k]; mj = k+mi;
      }
      flag = true;

      break;
    }
  }
  free(last_d+lower); free(rlast_d+rlower);
  free(temp_d+lower); free(rtemp_d+rlower);

  if (flag) {
    /* Find a path from (i1,j1) to (mi,mj) */
    sim4_align_path(seq1,seq2,i1,j1,mi,mj,midc,&head1,&tail1,M,N);

    /* Find a path from (mi,mj) to (i2,j2) */
    sim4_align_path(seq1,seq2,mi,mj,i2,j2,rmidc,&head2,&tail2,M,N);

    /* Join these two paths together */
    if (head1) tail1->next = head2;
    else head1 = head2;
  } else {
    (void)fprintf(stderr,
                  "align.c: warning: something wrong when dividing\n");
    head1 = NULL;
  }
  *head = head1;
  if (head2) *tail = tail2;
  else *tail = tail1;
}

static
int
snake(const char *seq1, const char *seq2, int k, int x, int endx, int endy) {
  int y;

  if (x<0) return x;
  y = x+k;
  while (x<endx && y<endy && seq1[x]==seq2[y]) {
    ++x; ++y;
  }
  return x;
}


static
int
rsnake(const char *seq1, const char *seq2, int k, int x, int startx, int starty, int M, int N) {
  int y;

  if (x>M) return x;
  if ((startx<0) || (starty<0))
    (void)printf("TROUBLE!!! startx:  %5d,  starty:  %5d\n",startx, starty);
  if ((x>M) || (x+k>N))
    (void)printf("TROUBLE!!! x:  %5d,  y:  %5d\n",x,x+k);

  y = x+k;
  while (x>startx && y>starty && seq1[x-1]==seq2[y-1]) {
    --x; --y;
  }
  return x;
}


void
Free_script(edit_script *head) {
  edit_script *tp, *tp1;

  tp = head;
  while (tp != NULL) {
    tp1 = tp->next;
    free(tp);
    tp = tp1;
  }
}
