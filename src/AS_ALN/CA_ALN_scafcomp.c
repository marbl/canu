
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

//int VarWindows[] = { 3, 4, 5, 6, 7, 8, 9, 10, 0 };
int VarWindows[] = { 3, 0 };
int MaxWindow = 0;

#define OVL_ERATE   .99
#define MIN_SEG_LEN  40
#define SEG_ERATE   .10

/* Debug conditional compilation flags */

#undef  DIAGNOSTICS
#undef  DEBUG_COMPARE
#undef  DEBUG_RAYSHOOT
#undef  DEBUG_ALIGN
#undef  DEBUG_ANALYSIS
#undef  DEBUG_LOCAL
#undef  DEBUG_OPTTAIL
#undef  DEBUG_SEGORDER
#undef  XFIG
#undef  STATS

#define ALLOW_NO_OVERLAP_INTERLEAVINGS


#define MAX_LINE_LEN   1000
#define MAX_SEQ_LEN  500000

/* Read a sequence of fast files for scaffolds and return a list of these. */

char *ScfCmp_Next_Line(FILE *ifile)
{ static char buffer[MAX_LINE_LEN+2];

  buffer[MAX_LINE_LEN+1] = '\n';
  if (fgets(buffer,MAX_SEQ_LEN+2,ifile) == NULL) return (NULL);
  if (buffer[MAX_LINE_LEN+1] != '\n')
    { fprintf(stderr,"Input line is longer than %d chars\n",MAX_LINE_LEN);
      exit (1);
    }
  return (buffer);
}

static Scaffold *Read_Next_Scaffold(FILE *ifile, int first)
{ static char seqbuf[MAX_SEQ_LEN+1];
  static char *line;

  Scaffold_Gap *gaplist;
  Scaffold_Tig *ctglist;
  Scaffold     *scaf;
  char         *packseq;
  int           numgaps;

  if (first)
    line = ScfCmp_Next_Line(ifile);
  if (line == NULL)
    return (NULL);

  if (line[0] != '>')
    { fprintf(stderr,"Header line should begin with '>'\n");
      exit (1);
    }

  { char *locate, *eptr;
    int   i;

    locate = strstr(line,"stddev:");
    if (locate == NULL)
      { fprintf(stderr,"Expecting gap variance in header line\n");
        exit (1);
      }
    locate += 7;

    numgaps = strtol(locate,&eptr,10);
    if (eptr == locate+7)
      { fprintf(stderr,"Expecting # of gaps after 'stddev:'\n");
        exit (1);
      }

    gaplist = (Scaffold_Gap *) malloc(sizeof(Scaffold_Gap)*numgaps);
    ctglist = (Scaffold_Tig *) malloc(sizeof(Scaffold_Tig)*(numgaps+1));

    locate = eptr;
    for (i = 0; i < numgaps; i++)
      { gaplist[i].gap_var = strtod(locate,&eptr);
        if (eptr == locate)
          { fprintf(stderr,"Expecting # of gaps after 'stddev:'\n");
            exit (1);
          }
        locate = eptr;
      }
  }
       
  { char *p;
    int   len;

    p = seqbuf;
    while ((line = ScfCmp_Next_Line(ifile)) != NULL)
      { if (line[0] == '>') break;
        len = strlen(line)-1;
        if ((p-seqbuf) + len > MAX_SEQ_LEN)
          { fprintf(stderr,"Sequence is longer than %d chars\n",
                           MAX_SEQ_LEN);
            exit (1);
          }
        strcpy(p,line);
        p += len;
      }
    *p = '\0';
  }

#ifdef DEBUG_READ
  { int len, i;

    len = strlen(seqbuf);
    for (i = 0; i < len; i += 60)
      fprintf(stderr,"%.60s\n",seqbuf+i);
  }
#endif

  { int gaps, plen;
    int i, j, p;
  
    gaps = 0;
    plen = 0;
    for (i = 0; seqbuf[i] != '\0'; i++)
      if (seqbuf[i] == 'n')
        { for (j = i+1; seqbuf[j] != '\0'; j++)
            if (seqbuf[j] != 'n')
              break;
          gaps += 1;
          i = j-1;
        }
      else
        plen += 1;

    if (gaps != numgaps)
      { fprintf(stderr,"Header gap count (%d) and N-gap count (%d) differ\n",
                       numgaps,gaps);
        exit (1);
      }

    packseq = (char *) malloc(plen+1);

    gaps = 0;
    p = 0;
    ctglist[gaps].insert_pnt = 0;
    ctglist[gaps].lft_end    = 0;
    for (i = 0; seqbuf[i] != '\0'; i++)
      if (seqbuf[i] == 'n')
        { ctglist[gaps].length = p - ctglist[gaps].insert_pnt;
          plen = 1;
          for (j = i+1; seqbuf[j] != '\0'; j++)
            if (seqbuf[j] != 'n')
              break;
            else
              plen += 1;
          gaplist[gaps].gap_length = plen;
          gaps += 1;
          ctglist[gaps].insert_pnt = p;
          ctglist[gaps].lft_end    = ctglist[gaps-1].lft_end
                                   + (ctglist[gaps-1].length + plen);
          i = j-1;
        }
      else
        packseq[p++] = seqbuf[i];
    packseq[p] = '\0';
    ctglist[gaps].length = p - ctglist[gaps].insert_pnt;
  }

  scaf = (Scaffold *) malloc(sizeof(Scaffold));
  scaf->num_gaps = numgaps;
  scaf->gaps = gaplist;
  scaf->ctgs = ctglist;
  scaf->length = ctglist[numgaps].lft_end + ctglist[numgaps].length;
  scaf->packed_seq = packseq;

  return (scaf);
}

Scaffold_List *Read_Multi_Fasta(char *fname)
{ FILE          *ifile;
  Scaffold_List *list = NULL, *fing;
  Scaffold      *S;
  int            first;

  if ((ifile = fopen(fname,"r")) == NULL)
    { fprintf(stderr,"Cannot open file '%s'\n",fname);
      exit (1);
    }

  fing = NULL;
  first = 1;
  while ((S = Read_Next_Scaffold(ifile,first)) != NULL)
    { if (fing == NULL)
        list = fing = (Scaffold_List *) malloc(sizeof(Scaffold_List));
      else
        fing = fing->next = (Scaffold_List *) malloc(sizeof(Scaffold_List));
      fing->scaffold = S;
      first = 0;
    }
  if (fing != NULL)
    fing->next = NULL;
  else
    list = NULL;

  fclose(ifile);

  return (list);
}

void Complement_Scaffold(Scaffold *S)
{ int i, j, plen, slen;
  Scaffold_Gap T;
  Scaffold_Tig U;

  plen = S->ctgs[S->num_gaps].insert_pnt + S->ctgs[S->num_gaps].length;
  slen = S->length;
  if(S->packed_seq!=NULL)  Complement_Seq(S->packed_seq);
  j = S->num_gaps-1;
  for (i = 0; i < j; i++, j--)
    { T = S->gaps[i];
      S->gaps[i] = S->gaps[j];
      S->gaps[j] = T;
    }
  j = S->num_gaps;
  for (i = 0; i < j; i++, j--)
    { U = S->ctgs[i];
      S->ctgs[i] = S->ctgs[j];
      S->ctgs[j] = U;
      S->ctgs[i].insert_pnt = plen - (S->ctgs[i].insert_pnt+S->ctgs[i].length);
      S->ctgs[j].insert_pnt = plen - (S->ctgs[j].insert_pnt+S->ctgs[j].length);
      S->ctgs[i].lft_end    = slen - (S->ctgs[i].lft_end+S->ctgs[i].length);
      S->ctgs[j].lft_end    = slen - (S->ctgs[j].lft_end+S->ctgs[j].length);
    }
  if (i == j)
    { S->ctgs[j].insert_pnt = plen - (S->ctgs[j].insert_pnt+S->ctgs[j].length);
      S->ctgs[j].lft_end    = slen - (S->ctgs[j].lft_end+S->ctgs[j].length);
    }
}

void Free_Scaffold(Scaffold *S)
{ free(S->gaps);
  free(S->ctgs);
  free(S->packed_seq);
  free(S);
}

static AggregateString *Build_Indexable_String(Scaffold_List *S)
{ int  *idx, *sfn;
  char *str;
  int   len, ctg;
  AggregateString *agg;

  { Scaffold_List *a;
    Scaffold      *s;
    int i;

    len = 0;
    ctg = 0;
    for (a = S; a != NULL; a = a->next)
      { s = a->scaffold;
        for (i = 0; i <= s->num_gaps; i++)
          { len += s->ctgs[i].length + 6;
            ctg += 1;
          }
      }
  }

  { int as;

    as = sizeof(int);
    as = ( (sizeof(AggregateString) + as - 1) / as ) * as;
    agg = (AggregateString *)
             malloc(as + 2*ctg*sizeof(int) + len + 1);
    if (agg == NULL)
      { fprintf(stderr,"Out of memory (aggregate string)\n");
        exit (1);
      }
    idx = (int *) (((char *) agg) + as);
    sfn = idx + ctg;
    str = (char *) (sfn + ctg);
  }

  { Scaffold_List *a;
    Scaffold      *s;
    int   i, c, n;
    char *p;

    c = n = 0;
    p = str;
    for (a = S; a != NULL; a = a->next)
      { s = a->scaffold;
        for (i = 0; i <= s->num_gaps; i++)
          { strcpy(p,"nnnnnn");
            p += 6;
            strncpy(p,s->packed_seq + s->ctgs[i].insert_pnt,
                      s->ctgs[i].length);
            p += s->ctgs[i].length;
            sfn[c] = n;
            idx[c] = p-str;
            c += 1;
          }
        n += 1;
      }
  }

  agg->bigstr = str;
  agg->parts  = idx;
  agg->scafs  = sfn;
  agg->slen   = len;
  agg->plen   = ctg;
  return (agg);
}

static int Locate_Local(AggregateString *cmp, int pos)
{ int l, r, m;

  l = 0;
  r = cmp->plen;
  while (l < r)
    { m = (l+r)/2;
      if (pos <= cmp->parts[m])
        r = m;
      else
        l = m+1;
    } 
  return (r);
}

Local_Address *MSORT_addr;
#define USE_ONLY_FORWARD_SEGS
#ifdef USE_ONLY_FORWARD_SEGS
Local_Segment *MSORT_segs;
#endif


int MSORT(const void *l, const void *r)
{ Local_Address *x, *y;

#ifdef USE_ONLY_FORWARD_SEGS
  Local_Segment *X, *Y; 
 int xrev,yrev;
#endif

  x = MSORT_addr + *((int *) l);
  y = MSORT_addr + *((int *) r);
  if (x->ascaf != y->ascaf)
    return (x->ascaf - y->ascaf);
  else if (x->bscaf != y->bscaf)
    return (x->bscaf - y->bscaf);
  else if (x->acntg != y->acntg)
    return (x->acntg - y->acntg);
  else if (x->bcntg != y->bcntg)
      return (x->bcntg - y->bcntg);

#ifdef USE_ONLY_FORWARD_SEGS
  X = MSORT_segs + *((int *) l);
  Y = MSORT_segs + *((int *) r);
  assert(X->abpos<X->aepos);
  assert(Y->abpos<Y->aepos);
  xrev = ( X->bbpos<X->bepos) ? 0 : 1;
  yrev = ( Y->bbpos<Y->bepos) ? 0 : 1;
  return (xrev - yrev);
#endif

  return(0);

}

static Local_Pool *Find_All_Locals(Scaffold_List *AS, Scaffold_List *BS)
{ AggregateString *acmp, *bcmp;
  int             *perm;
  int              nseg;
  Local_Segment   *segs;
  Local_Address   *addr;
  static Local_Pool pool;
 
  acmp = Build_Indexable_String(AS);
  bcmp = Build_Indexable_String(BS);
  segs = Find_Local_Segments(acmp->bigstr,acmp->slen,
                             bcmp->bigstr,bcmp->slen,LOCAL_BOTH,
                             MIN_SEG_LEN,SEG_ERATE,&nseg);

  { int i;

    addr = (Local_Address *) malloc(nseg*sizeof(Local_Address));
    perm = (int *) malloc(nseg*sizeof(int));
    for (i = 0; i < nseg; i++)
      { int ac, bc;

        perm[i] = i;
        addr[i].acntg = ac = Locate_Local(acmp,segs[i].abpos);
        addr[i].bcntg = bc = Locate_Local(bcmp,segs[i].bbpos);
        addr[i].ascaf = acmp->scafs[ac];
        addr[i].bscaf = bcmp->scafs[bc];

        if (ac == 0)
         ac = 6;
        else
          ac = acmp->parts[ac-1] + 6;
        if (bc == 0)
          bc = 6;
        else
          bc = bcmp->parts[bc-1] + 6;
        segs[i].abpos -= ac;
        segs[i].aepos -= ac;
        segs[i].bbpos -= bc;
        segs[i].bepos -= bc;
        segs[i].ldiag -= (ac-bc);
        segs[i].hdiag -= (ac-bc);
      }
  }

  MSORT_addr = addr;
#ifdef USE_ONLY_FORWARD_SEGS
  MSORT_segs = segs;
#endif
  qsort(perm,nseg,sizeof(int),MSORT);
  
#ifdef DEBUG_LOCAL
  { int i;

    fprintf(stderr,"\n  A index\n");
    for (i = 0; i < acmp->plen; i++)
      fprintf(stderr,"   %6d %3d\n",acmp->parts[i],acmp->scafs[i]);
    fprintf(stderr,"\n  B index\n");
    for (i = 0; i < bcmp->plen; i++)
      fprintf(stderr,"   %6d %3d\n",bcmp->parts[i],bcmp->scafs[i]);
    fprintf(stderr,"\nHits\n");
    for (i = 0; i < nseg; i++)
      fprintf(stderr,"%3d: %3d [%6d,%6d] vs [%6d,%6d] to (%d,%d) [%d,%d]\n",
             i,perm[i],segs[i].abpos,segs[i].aepos,
             segs[i].bbpos,segs[i].bepos,addr[i].acntg,addr[i].bcntg,
             addr[i].ascaf,addr[i].bscaf);
  }
#endif

  free(acmp);
  free(bcmp);

  { int a, w;
    Local_Segment t;
    Local_Address u;

    for (a = 0; a < nseg; a++)
      { t = segs[a];
        u = addr[a];
        while ((w = perm[a]) >= 0)
          { if (perm[w] < 0)
              { segs[a] = t;
                addr[a] = u;
              }
            else
              { segs[a] = segs[w];
                addr[a] = addr[w];
              }
            perm[a] = -1;
            a = w;
          }
      }
  }

  free(perm);

#ifdef DEBUG_LOCAL
  { int i;

    fprintf(stderr,"\nSorted Hits\n");
    for (i = 0; i < nseg; i++)
      fprintf(stderr,"%3d: %3d [%6d,%6d] vs [%6d,%6d] to (%d,%d) [%d,%d]\n",
             i,perm[i],segs[i].abpos,segs[i].aepos,
             segs[i].bbpos,segs[i].bepos,addr[i].acntg,addr[i].bcntg,
             addr[i].ascaf,addr[i].bscaf);
  }
#endif

  pool.num_locals = nseg;
  pool.locals     = segs;
  pool.address    = addr;
  return (&pool);
}

Segment *Find_All_Overlaps(Scaffold *AF, Scaffold *BF, Local_Pool *pool,
                                  int as, int bs, int ac, int bc,
                                  int *fing, int *segs, int comp)
{ int      i, alen;
  int      bf;
  Local_Address *addr;
  Local_Segment *locl;
  Segment *list;

  list  = NULL;
  *segs = 0;

  bf   = *fing;
  addr = pool->address;
  locl = pool->locals;
  for (i = 0; i <= AF->num_gaps; i++)
    { int j, blen;
      Segment *last;

      alen = AF->ctgs[i].length;
      last = list;

      for (j = 0; j <= BF->num_gaps; j++)
        { Local_Overlap *ovl;
          int n;

	  // The following looks odd, but remember, if comp==1,
	  // then the scaffold was reversed, but the hits weren't
	  //     --- ALH, 11/7/01
          if (comp)
            blen = BF->ctgs[BF->num_gaps-j].length;
          else
            blen = BF->ctgs[j].length;

#ifdef USE_ONLY_FORWARD_SEGS
	  //if reversed comparison, skip over forward-oriented matches
	  if(comp){
	    while (bf < pool->num_locals &&
		   addr[bf].ascaf == as && addr[bf].bscaf == bs &&
		   addr[bf].acntg == ac+i && addr[bf].bcntg == bc+j&&
		   locl[bf].bbpos < locl[bf].bepos)
	      bf++;
	  }
#endif

          n = bf;
          while (n < pool->num_locals &&
                 addr[n].ascaf == as && addr[n].bscaf == bs &&
                 addr[n].acntg == ac+i && addr[n].bcntg == bc+j){
#ifdef USE_ONLY_FORWARD_SEGS
	    //if forward comparison, stop on reversed match
	    if(!comp&&locl[n].bbpos > locl[bf].bepos)break;
#endif
            n += 1;
	  }

#ifdef DEBUG_LOCAL
          { int x;
            fprintf(stderr,"\nInput to Find_Local_Overlap (%d,%d) [%d,%d] fing = %d\n",
                   ac+i,bc+j,as,bs,bf);
            for (x = bf; x < n; x++)
              fprintf(stderr,"%3d: [%6d,%6d] vs [%6d,%6d] to (%d,%d) [%d,%d]\n",
                     x,locl[x].abpos,locl[x].aepos,
                     locl[x].bbpos,locl[x].bepos,addr[x].acntg,addr[x].bcntg,
                     addr[x].ascaf,addr[x].bscaf);
          }
#endif

          ovl = Find_Local_Overlap(alen,blen,comp,0,locl+bf,n-bf,20,OVL_ERATE);

#ifdef USE_ONLY_FORWARD_SEGS
	  // if forward comparison, skip over reversed matches
	  if(!comp){
	    bf=n;
	    while (bf < pool->num_locals &&
		   addr[bf].ascaf == as && addr[bf].bscaf == bs &&
		   addr[bf].acntg == ac+i && addr[bf].bcntg == bc+j)
	      bf++;
	    n=bf;
	  }
#endif


          while (ovl != NULL)
            { Segment *newseg;

              newseg = (Segment *) malloc(sizeof(Segment));
              newseg->next     = list;
              newseg->a_contig = i;
              newseg->b_contig = j;

              if (ovl->begpos < 0)
                newseg->alow = 0;
              else
                newseg->alow = ovl->begpos;
              if (ovl->begpos > 0)
                newseg->blow = 0;
              else
                newseg->blow = - ovl->begpos;

              if (ovl->endpos > 0)
                newseg->ahgh = alen;
              else
                newseg->ahgh = alen + ovl->endpos;
              if (ovl->endpos < 0)
                newseg->bhgh = blen;
              else
                newseg->bhgh = blen - ovl->endpos;

              newseg->overlap  = ovl;
              list = newseg;
              *segs += 1;

#ifdef DEBUG_LOCAL
              Print_Local_Overlap(stdout,ovl,4);
#endif

              ovl = Find_Local_Overlap(alen,blen,comp,1,
                                       locl+bf,n-bf,20,OVL_ERATE);
            }

          bf  = n;
        }

      // The following looks odd, but remember, if comp==1,
      // then the scaffold was reversed, but the hits weren't
      //     --- ALH, 11/7/01
      if (comp)
        { Segment *r, *s;
          for (r = last; list != last; list = s)
            { s = list->next;
              list->b_contig = BF->num_gaps - list->b_contig;
              list->next = r;
              r = list;
            }
          list = r;
        }
    }

  *fing = bf;
  return (list);
}

void Free_Segments_ScafComp(Segment *seglist)
{ Segment *s, *t;

  for (s = seglist; s != NULL; s = t)
    { t = s->next;
      free(s->overlap);
      free(s);
    }
}

int     MaxAlign  = -1;
int     MaxBucket = -1;
COvlps  *CtgOvls  = NULL;
COvlps **ABuckets = NULL;
COvlps **BBuckets = NULL;
#ifdef XFIG
static double scale_factor;
static FILE *figfile;

#define MAP(x)  ((int) ((x)*scale_factor + 600))
#define UNMAP(y)  (int)(((y)-600.)/scale_factor)

#endif

#ifdef XFIG
static void Draw_Matrix(Segment *seglist, int numsegs,
                 Scaffold *AF, Scaffold *BF,
                 int varwin)
{ Segment *f;
  int i, j;

  if (AF->length > BF->length)
    scale_factor = 6000. / AF->length;
  else
    scale_factor = 6000. / BF->length;
  fprintf(figfile,"#FIG 3.2\n");
  fprintf(figfile,"Landscape\n");
  fprintf(figfile,"Center\n");
  fprintf(figfile,"Inches\n");
  fprintf(figfile,"Letter\n");
  fprintf(figfile,"100.\n");
  fprintf(figfile,"Single\n");
  fprintf(figfile,"-2\n");
  fprintf(figfile,"1200 2\n");
  fprintf(figfile,"2 2 0 3 0 -1 100 0 -1 4.00 0 0 0 0 0 5\n\t");
  fprintf(figfile," %d %d",MAP(0),MAP(0));
  fprintf(figfile," %d %d",MAP(0),MAP(BF->length));
  fprintf(figfile," %d %d",MAP(AF->length),MAP(BF->length));
  fprintf(figfile," %d %d",MAP(AF->length),MAP(0));
  fprintf(figfile," %d %d\n",MAP(0),MAP(0));

  // draw boxes for contigs
  for (i = 0; i <= AF->num_gaps; i++)
    for (j = 0; j <= BF->num_gaps; j++)
      { int xl, yl, xh, yh;
        xl = MAP(AF->ctgs[i].lft_end);
        yl = MAP(BF->ctgs[j].lft_end);
        xh = MAP(AF->ctgs[i].lft_end + AF->ctgs[i].length);
        yh = MAP(BF->ctgs[j].lft_end + BF->ctgs[j].length);
        fprintf(figfile,"2 2 0 2 1 -1 101 0 -1 4.00 0 0 0 0 0 5\n\t");
        fprintf(figfile," %d %d",xl,yl);
        fprintf(figfile," %d %d",xl,yh);
        fprintf(figfile," %d %d",xh,yh);
        fprintf(figfile," %d %d",xh,yl);
        fprintf(figfile," %d %d\n",xl,yl);
      }

  // draw gap symbols?
  for (i = 0; i < AF->num_gaps; i++)
    { int xl, xh, yl;
      xl = MAP(AF->ctgs[i+1].lft_end - varwin * AF->gaps[i].gap_var);
      xh = MAP(AF->ctgs[i+1].lft_end + varwin * AF->gaps[i].gap_var);
      yl = 500;
      fprintf(figfile,"2 1 0 1 1 -1 102 0 -1 4.00 0 0 0 1 1 2\n\t");
      fprintf(figfile,"0 0 1.00 60.00 75.00\n\t");
      fprintf(figfile,"0 0 1.00 60.00 75.00\n\t");
      fprintf(figfile," %d %d",xl,yl);
      fprintf(figfile," %d %d\n",xh,yl);
    }
  // draw gap symbols?
  for (i = 0; i < BF->num_gaps; i++)
    { int yl, yh, xl;
      yl = MAP(BF->ctgs[i+1].lft_end - varwin * BF->gaps[i].gap_var);
      yh = MAP(BF->ctgs[i+1].lft_end + varwin * BF->gaps[i].gap_var);
      xl = 500;
      fprintf(figfile,"2 1 0 1 1 -1 102 0 -1 4.00 0 0 0 1 1 2\n\t");
      fprintf(figfile,"0 0 1.00 60.00 75.00\n\t");
      fprintf(figfile,"0 0 1.00 60.00 75.00\n\t");
      fprintf(figfile," %d %d",xl,yl);
      fprintf(figfile," %d %d\n",xl,yh);
    }

  // draw contig overlaps
  for (f = seglist; f != NULL; f = f->next)
    { int xl, yl, xh, yh;
      xl = MAP(AF->ctgs[f->a_contig].lft_end + f->alow);
      xh = MAP(AF->ctgs[f->a_contig].lft_end + f->ahgh);
      yl = MAP(BF->ctgs[f->b_contig].lft_end + f->blow);
      yh = MAP(BF->ctgs[f->b_contig].lft_end + f->bhgh);
      fprintf(figfile,"2 1 0 2 2 -1 102 0 -1 4.00 0 0 0 0 0 2\n\t");
      fprintf(figfile," %d %d",xl,yl);
      fprintf(figfile," %d %d\n",xh,yh);
    }
}
#endif

#ifdef STATS
static int gapcount, blockcount, linkcount, edgecount;
static int vertcount, boundcount;

void Start_Stats()
{ gapcount = blockcount = linkcount = edgecount = 0;
  vertcount = boundcount = 0;
}

void Print_Stats()
{ fprintf(stderr,"  Block Boundaries Reached  = %3d\n",blockcount);
  fprintf(stderr,"  Gap Portals Reached       = %3d\n",gapcount);
  fprintf(stderr,"  Traverse Calls Made (Bnd) = %3d(%d)\n",linkcount,boundcount);
  fprintf(stderr,"  Edges in Graphs          <= %3d\n",edgecount);
  fprintf(stderr,"  Vertices in Graph         = %3d\n",vertcount);
}
#endif

// BEGINNING OF NEW CODE MARCH 2006

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

interval_list *cleanup_ilist(interval_list *list){
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
    //  munge the memory
    list->prev = 0x87654321;
    list->next = 0x12345678;
    free(list);
    list=link;
  }
  return NULL;
}

interval_list *add_to_ilist(interval_list *tail, interval to_add);

interval_list *add_to_ilist_special(interval_list *tail, interval to_add){
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

interval_list *add_to_ilist(interval_list *tail, interval to_add){
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

  il = (interval_list*) malloc(sizeof(interval_list));
  assert(il!=NULL);
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

int Project_across_Agap_one_interval(interval *inoutIval,COvlps **bestTerm, double Agap_length, double Agap_var,Scaffold *B, int varwin){
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
    remainingToUseUp = max(0, Agap_length + Agap_var * varwin);
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

int Project_across_Agap(interval_list *inList,interval_list **outList, COvlps **bestTerm,double Agap_length, double Agap_var,Scaffold *B, int varwin){

// To process the left edge of an A-gap (list of accessible intervals)
// For each vertical interval,
//   If terminal "gap"
//     If best so far, update best
//   Else 
//     project across A-gap
//     If go off bottom, 
//       If best so far, update best

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

int Process_seg(COvlps *af,COvlps **bestTerm,
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

int Process_Agap_one_accessible_interval(interval curr,COvlps **bestTerm,
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
      double from=max(top,B->ctgs[i].lft_end);
      double to=min(bot,ctgend);
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
	      int score = af->seg->overlap->length - max(0,-af->seg->overlap->begpos) - max(0,-af->seg->overlap->endpos);
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
      from=max(top,ctgend);
      to=min(bot,ctgend+B->gaps[i].gap_length);

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
	  left = max(0, frac * ((double) B->gaps[i].gap_length - B->gaps[i].gap_var*varwin) );
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
		int score = af->seg->overlap->length - max(0,-af->seg->overlap->begpos) - max(0,-af->seg->overlap->endpos);
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


	  fracBot = ((double) max(0,left - A->ctgs[whichA].length )) / (B->gaps[i].gap_length - B->gaps[i].gap_var*varwin); 
	
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
      double from=max(top,B->ctgs[i].lft_end);
      double to=min(bot,ctgend);
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
	      int score = af->seg->overlap->length - max(0,-af->seg->overlap->begpos) - max(0,-af->seg->overlap->endpos);
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

int Process_Agap_accessible_intervals(interval_list **accessList,
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
	  af->best = af->seg->overlap->length - max(0,-af->seg->overlap->begpos) - max(0,-af->seg->overlap->endpos);
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


int Project_segments_across_Bgaps(COvlps **bestTerm, interval_list **BgapLists, interval_list **nextAgapList, int whichA, Scaffold *A, Scaffold *B, int varwin){
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
#if 0
	    if(af->seg->overlap->begpos <= 0 ) { //ahang is nonpositive, i.e. segment starts in B contig
	      double pnt = - af->seg->overlap->begpos;
	      if (from <= pnt && pnt <= to) {
		int score = af->seg->overlap->length - max(0,-af->seg->overlap->begpos) - max(0,-af->seg->overlap->endpos);
		assert(score>0);
		if (curr->ival.traceback != NULL)
		  score += curr->ival.traceback->best;
		if (score > af->best){
		  af->best  = score;
		  af->trace = curr->ival.traceback;
		}
		// now process the segment's projections
		if(Process_seg(af,bestTerm,nextAgapList,BgapLists,whichA,i,A,B,varwin)){
		  terminal=1;
		}
	      }
	    }

#else

	    if(af->seg->overlap->begpos >= 0 ) { //ahang is nonnegative, i.e. segment starts in A contig
	      double pnt = af->seg->overlap->begpos;
	      if (from <= pnt && pnt <= to) {
		int score = af->seg->overlap->length - max(0,-af->seg->overlap->begpos) - max(0,-af->seg->overlap->endpos);
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

#endif

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


int ProjectFromTopEdge(interval_list **thisAlist, COvlps **bestTerm,int whichA,Scaffold *A, Scaffold *B,int varwin,int bandbeg,int bandend){

  int terminal=0;
  double low,high,top,bot;
  double minX,maxX;
  double gapFrac;

  assert(*thisAlist==NULL);

  // find relevant interval in A-scaffold coordinates:
  low=max(bandbeg,A->ctgs[whichA].lft_end+A->ctgs[whichA].length);
  high=min(bandend,A->ctgs[whichA+1].lft_end);

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
  minX = max(0,(gapFrac * A->gaps[whichA].gap_length - A->gaps[whichA].gap_var * varwin)); 
  
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



#ifdef AHANG_BAND_TEST
  assert(bandbeg<=bandend);
  assert(bandend<=AF->length);
  assert(bandbeg>=-(BF->length));
#endif

  if (numsegs > MaxAlign)
    { MaxAlign = (int)(1.3*numsegs + 100);
      CtgOvls  = (COvlps *) realloc(CtgOvls,sizeof(COvlps)*MaxAlign);
      if (CtgOvls == NULL)
        { fprintf(stderr,"Out of memory allocating DP array\n");
          exit (1);
        }
    }

  if (AF->num_gaps + BF->num_gaps + 2 > MaxBucket)
    { MaxBucket = (int)(1.3*(AF->num_gaps + BF->num_gaps + 2) + 100);
      ABuckets  = (COvlps **) realloc(ABuckets,sizeof(COvlps *)*MaxBucket);
      if (ABuckets == NULL)
        { fprintf(stderr,"Out of memory allocating segment sort arrays\n");
          exit (1);
        }
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
     AF->gaps[i].gap_length=max(1,diff);
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
     BF->gaps[i].gap_length=max(diff,1);
     assert(BF->ctgs[i].lft_end+BF->ctgs[i].length+BF->gaps[i].gap_length == BF->ctgs[i+1].lft_end);
     into+=BF->gaps[i].gap_length;
     i++;
   }

 }

  // Initialize along left edge of first A contig
  if(bandbeg<0){
    interval startupAgap;
    int top=-min(0,bandend);
    int bot=-min(0,bandbeg);
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
      left = max(bandbeg,AF->ctgs[i].lft_end);
      right = min(bandend,AF->ctgs[i].lft_end+AF->ctgs[i].length);
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
	{ free(s->overlap);
	free(s);
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



int Link_Horizontal(Scaffold *A, Scaffold *B, int varwin,
                           int i, int j, int low, int hgh, COvlps *source)
{ int k, var;
  int terminal;
  COvlps *afing;
#ifdef XFIG
  int xl, zl, zh, gapped;

  if (low != hgh)
    { xl = A->ctgs[i].lft_end + A->ctgs[i].length + B->ctgs[j].lft_end;
      zl = MAP(xl - hgh > A->ctgs[i].lft_end + A->ctgs[i].length ? 
	       xl - hgh : A->ctgs[i].lft_end + A->ctgs[i].length);
      zh = MAP(xl - low > A->ctgs[i].lft_end + A->ctgs[i].length ? 
	       xl - low : A->ctgs[i].lft_end + A->ctgs[i].length);
      xl = MAP(B->ctgs[j].lft_end);
      gapped = 1;
    }
  else
    { xl  = MAP(A->ctgs[i].lft_end + A->ctgs[i].length);
      zl = MAP(low);
      zh = MAP(hgh);
      gapped = 0;
    }
#endif

#ifdef STATS
  linkcount += 1;
#endif
  if (i == A->num_gaps)
    {
#ifdef DEBUG_ALIGN
      fprintf(stderr,"    Ray reaches A-boundary\n"); 
#endif
      return (1);
    }

  terminal = 0;
  low += A->gaps[i].gap_length;
  hgh += A->gaps[i].gap_length;
  var  = (int)(A->gaps[i].gap_var * varwin);
  afing = ABuckets[i+1];

  for (k = j; k <= B->num_gaps; k++)
    { int x, beg, end;

      x = B->ctgs[k].lft_end;
      if (low-var > x)
        beg = low-var;
      else
        beg = x; 
      x += B->ctgs[k].length;
      if (hgh+var < x)
        end = hgh+var;
      else
        end = x;
      if (beg <= end)
        { COvlps *af;
          int     pnt, score;

#ifdef STATS
          blockcount += 1;
#endif
#ifdef XFIG
          { int xh, yl, yh;

            yl = MAP(beg);
            yh = MAP(end);
            xh  = MAP(A->ctgs[i+1].lft_end);
            fprintf(figfile,"2 3 0 1 4 4  98 0 35 4.01 1 0 0 0 0 5\n\t");
            if (gapped)
              { fprintf(figfile," %d %d",zl,xl);
                fprintf(figfile," %d %d",zh,xl);
                fprintf(figfile," %d %d",xh,yl);
                fprintf(figfile," %d %d",xh,yh);
                fprintf(figfile," %d %d\n",zl,xl);
              }
            else
              { fprintf(figfile," %d %d",xl,zl);
                fprintf(figfile," %d %d",xl,zh);
                fprintf(figfile," %d %d",xh,yh);
                fprintf(figfile," %d %d",xh,yl);
                fprintf(figfile," %d %d\n",xl,zl);
              }
          }
#endif
#ifdef DEBUG_RAYSHOOT
          fprintf(stderr,"    Ray overlaps B%d[%d,%d]\n",k,beg,end);
#endif
          while (afing != NULL && afing->seg->b_contig < k)
            afing = afing->Alink; 
          for (af = afing; af != NULL && af->seg->b_contig == k; af = af->Alink)
            if (af->seg->overlap->begpos <= 0)
              { pnt = B->ctgs[k].lft_end - af->seg->overlap->begpos;
                if (beg <= pnt && pnt <= end)
                  { score = af->seg->overlap->length - max(0,-af->seg->overlap->begpos) - max(0,-af->seg->overlap->endpos);
		  assert(score>0);
                    if (source != NULL)
                      score += source->best;
                    if (score > af->best)
                      { af->best  = score;
                        af->trace = source;
                      }
#ifdef STATS
                    edgecount += 1;
#endif
#ifdef DEBUG_ALIGN
                    fprintf(stderr,"    Finds (%d,%d) hangs(%d,%d)= %d\n",
                           af->seg->a_contig,af->seg->b_contig,af->seg->overlap->begpos,af->seg->overlap->endpos,af->best);
#endif

// If you are reading this, it is possible you are confused as to whether we should 
// recurse on segments reached from Link_Horizontal or Link_Vertical; it appears that
// we should NOT.  Instead, Align_Scaffold proceeds as follows, in three stages:
// (a) determine all the segments that can be reached directly from starting 
// in a gap on the entry border; each such segment has its "best" set to its length
// (b) determine all segments that can be reached by starting from a segment that
// begins on the entry border; again, reachable segments have "best" set to length
// (c) consider all segments with best >= 0 (i.e. those that can be reached from
// the entry border via gaps or segments previously visited); determine where these
// segments can reach--any segment they can reach itself becomes reachable.
// To avoid redundant computation for segments in path tails, once we reach
// a segment, we do not recursively follow it but instead mark it as accessible (i.e.
// set its best > 0).  For this to work, we must be sure to evaluate segments in the
// right order -- it would work to go in increasing order of A contigs, with subsorting
// on increasing B contigs, or vice versa ... but critically we can't go in reverse order
// on either.
//
// As an aside, it appears that repeated evaluation of the same gap intervals is not
// prevented by this scheme ....  Better might be to proceed one A contig or gap at a time.
//
                  }
              }
        }

      if (k < B->num_gaps)
        { int y, l, del;

          if (low-var > x)
            beg = low-var;
          else
            beg = x; 
          y = B->ctgs[k+1].lft_end;
          if (hgh+var < y)
	    if(hgh+var<x){
	      continue;
	    } else {
	      end = hgh+var;
	    }
	  else
            end = y;
          del = (int)(B->gaps[k].gap_var * varwin);
          if (beg <= end)
            { l = A->ctgs[i+1].lft_end;
#ifdef STATS
              gapcount += 1;
#endif
#ifdef XFIG
              { int xh, yl, yh;

                yl = MAP(beg);
                yh = MAP(end);
                xh  = MAP(A->ctgs[i+1].lft_end);
                fprintf(figfile,"2 3 0 1 4 4  98 0 35 4.02 1 0 0 0 0 5\n\t");
                if (gapped)
                  { fprintf(figfile," %d %d",zl,xl);
                    fprintf(figfile," %d %d",zh,xl);
                    fprintf(figfile," %d %d",xh,yl);
                    fprintf(figfile," %d %d",xh,yh);
                    fprintf(figfile," %d %d\n",zl,xl);
                  }
                else
                  { fprintf(figfile," %d %d",xl,zl);
                    fprintf(figfile," %d %d",xl,zh);
                    fprintf(figfile," %d %d",xh,yh);
                    fprintf(figfile," %d %d",xh,yl);
                    fprintf(figfile," %d %d\n",xl,zl);
                  }
              }
#endif
#ifdef DEBUG_RAYSHOOT
              fprintf(stderr,"    Ray overlaps gap B%d:%d to B%d:%d\n",k,beg,k+1,end);
#endif

              if (Link_Vertical(A,B,varwin,
                                i+1,k,l - (end-x),l - (beg-x),source)) 
                terminal = 1;
            }
          else if (beg > y && beg-del < y)
            { l = A->ctgs[i+1].lft_end;
#ifdef STATS
              gapcount += 1;
#endif
#ifdef DEBUG_RAYSHOOT
              fprintf(stderr,"    Ray indirects gap B%d to B%d @ %d\n",k,k+1,beg-del);
#endif
              if (Link_Vertical(A,B,varwin,
                  i+1,k,l - (beg-x),l - (beg-x),source)) 
                terminal = 1;
            }
          var += del;
        }

      else if (end == x)
        { terminal = 1;
#ifdef XFIG
          { int tl, th, yh;

            tl = A->ctgs[i+1].lft_end - ((hgh+var) - end);
            if (tl < A->ctgs[i].lft_end + A->ctgs[i].length)
              tl = A->ctgs[i].lft_end + A->ctgs[i].length;
            yh = MAP(end);
            tl = MAP(tl);
            // th = MAP(A->ctgs[i+1].lft_end);
	    th = A->ctgs[i+1].lft_end - ((low-var) -end);
            if (th > A->ctgs[i+1].lft_end)
              th = A->ctgs[i+1].lft_end ;  /* tl -> th; did this fix anything? */
    	    th = MAP(th);
            fprintf(figfile,"2 3 0 1 4 4  98 0 35 4.03 1 0 0 0 0 5\n\t");
            if (gapped)
              { fprintf(figfile," %d %d",zl,xl);
                fprintf(figfile," %d %d",zh,xl);
                fprintf(figfile," %d %d",th,yh);
                fprintf(figfile," %d %d",tl,yh);
                fprintf(figfile," %d %d\n",zl,xl);
              }
            else
              { fprintf(figfile," %d %d",xl,zl);
                fprintf(figfile," %d %d",xl,zh);
                fprintf(figfile," %d %d",tl,yh);
                fprintf(figfile," %d %d",th,yh);
                fprintf(figfile," %d %d\n",xl,zl);
              }
          }
#endif
#ifdef DEBUG_ALIGN
          fprintf(stderr,"    Ray reaches B-boundary\n"); 
#endif
        }
    }

  return (terminal);
}

int Link_Vertical(Scaffold *A, Scaffold *B, int varwin,
                         int i, int j, int low, int hgh, COvlps *source)
{ int k, var;
  int terminal;
  COvlps *bfing;
#ifdef XFIG
  int yl, zl, zh, gapped;

  if (low != hgh)
    { yl = B->ctgs[j].lft_end + B->ctgs[j].length + A->ctgs[i].lft_end;
      zl = MAP(yl - hgh > B->ctgs[j].lft_end + B->ctgs[j].length ? 
	       yl - hgh : B->ctgs[j].lft_end + B->ctgs[j].length);
      zh = MAP(yl - low > B->ctgs[j].lft_end + B->ctgs[j].length ? 
	       yl - low : B->ctgs[j].lft_end + B->ctgs[j].length);
      yl = MAP(A->ctgs[i].lft_end);
      gapped = 1;
    }
  else
    { yl  = MAP(B->ctgs[j].lft_end + B->ctgs[j].length);
      zl = MAP(low);
      zh = MAP(hgh);
      gapped = 0;
    }
#endif

#ifdef STATS
  linkcount += 1;
#endif
  if (j == B->num_gaps)
    {
#ifdef DEBUG_ALIGN
      fprintf(stderr,"    Ray reaches B-boundary\n"); 
#endif
      return (1);
    }

  terminal = 0;
  low += B->gaps[j].gap_length;
  hgh += B->gaps[j].gap_length;
  var  = (int)(B->gaps[j].gap_var * varwin);
  bfing = BBuckets[j+1];

  for (k = i; k <= A->num_gaps; k++)
    { int x, beg, end;

      x = A->ctgs[k].lft_end;
      if (low-var > x)
        beg = low-var;
      else
        beg = x; 
      x += A->ctgs[k].length;
      if (hgh+var < x)
        end = hgh+var;
      else
        end = x;
      if (beg <= end)
        { COvlps *bf;
          int     pnt, score;

#ifdef STATS
          blockcount += 1;
#endif
#ifdef XFIG
          { int yh, xl, xh;

            xl = MAP(beg);
            xh = MAP(end);
            yh  = MAP(B->ctgs[j+1].lft_end);
            fprintf(figfile,"2 3 0 1 5 5  98 0 35 4.04 1 0 0 0 0 5\n\t");
            if (gapped)
              { fprintf(figfile," %d %d",yl,zl);
                fprintf(figfile," %d %d",yl,zh);
                fprintf(figfile," %d %d",xl,yh);
                fprintf(figfile," %d %d",xh,yh);
                fprintf(figfile," %d %d\n",yl,zl);
              }
            else
              { fprintf(figfile," %d %d",zl,yl);
                fprintf(figfile," %d %d",zh,yl);
                fprintf(figfile," %d %d",xh,yh);
                fprintf(figfile," %d %d",xl,yh);
                fprintf(figfile," %d %d\n",zl,yl);
              }
          }
#endif
#ifdef DEBUG_RAYSHOOT
          fprintf(stderr,"    Ray overlaps A%d[%d,%d]\n",k,beg,end);
#endif
          while (bfing != NULL && bfing->seg->a_contig < k){
#ifdef DEBUG_RAYSHOOT
	    fprintf(stderr,"Skipping over segment involving A%d, ahg:%d\n",
		    bfing->seg->a_contig,
		    bfing->seg->overlap->begpos);
#endif
            bfing = bfing->Blink; 
	  }
          for (bf = bfing; bf != NULL && bf->seg->a_contig == k; bf = bf->Blink){
#ifdef DEBUG_RAYSHOOT
	    fprintf(stderr,"Examining segment involving A%d ahg:%d\n",
		    bf->seg->a_contig,
		    bf->seg->overlap->begpos);
#endif
            if (bf->seg->overlap->begpos >= 0)
              { pnt = A->ctgs[k].lft_end + bf->seg->overlap->begpos;
                if (beg <= pnt && pnt <= end)
                  { score = bf->seg->overlap->length  - max(0,-bf->seg->overlap->begpos) - max(0,-bf->seg->overlap->endpos);
		  assert(score>0);
                    if (source != NULL)
                      score += source->best;
                    if (score > bf->best)
                      { bf->best  = score;
                        bf->trace = source;
                      }
#ifdef STATS
                    edgecount += 1;
#endif
#ifdef DEBUG_ALIGN
                    fprintf(stderr,"    Finds (%d,%d) hangs(%d,%d)= %d\n",
                           bf->seg->a_contig,bf->seg->b_contig,bf->seg->overlap->begpos,bf->seg->overlap->endpos,bf->best);
#endif
                  }
              }
	  }
        }

      if (k < A->num_gaps)
        { int y, l, del;

          if (low-var > x)
            beg = low-var;
          else
            beg = x; 
          y = A->ctgs[k+1].lft_end;
          if (hgh+var < y)
	    if(hgh+var < x){
	      continue;
	    } else {
	      end = hgh+var;
	    }
          else
            end = y;
          del = (int)(A->gaps[k].gap_var * varwin);
          if (beg <= end)
            { l = B->ctgs[j+1].lft_end;
#ifdef STATS
              gapcount += 1;
#endif
#ifdef XFIG
              { int yh, xl, xh;

                xl = MAP(beg);
                xh = MAP(end);
                yh  = MAP(B->ctgs[j+1].lft_end);
                fprintf(figfile,"2 3 0 1 5 5  98 0 35 4.05 1 0 0 0 0 5\n\t");
                if (gapped)
                  { fprintf(figfile," %d %d",yl,zl);
                    fprintf(figfile," %d %d",yl,zh);
                    fprintf(figfile," %d %d",xl,yh);
                    fprintf(figfile," %d %d",xh,yh);
                    fprintf(figfile," %d %d\n",yl,zl);
                  }
                else
                  { fprintf(figfile," %d %d",zl,yl);
                    fprintf(figfile," %d %d",zh,yl);
                    fprintf(figfile," %d %d",xh,yh);
                    fprintf(figfile," %d %d",xl,yh);
                    fprintf(figfile," %d %d\n",zl,yl);
                  }
              }
#endif
#ifdef DEBUG_RAYSHOOT
              fprintf(stderr,"    Ray overlaps gap A%d:%d to A%d:%d\n",k,beg,k+1,end);
#endif
              if (Link_Horizontal(A,B,varwin,
                                  k,j+1,l - (end-x),l - (beg-x),source))
                terminal = 1; 
            }
          else if (beg > y && beg-del < y)
            { l = B->ctgs[j+1].lft_end;
#ifdef STATS
              gapcount += 1;
#endif
#ifdef DEBUG_RAYSHOOT
              fprintf(stderr,"    Ray indirects gap A%d to A%d @ %d\n",k,k+1,beg-del);
#endif
              if (Link_Horizontal(A,B,varwin,
                                  k,j+1,l - (beg-x),l - (beg-x),source))
                terminal = 1;
            }
          var += del;
        }

      else if (end == x)
        { terminal = 1;
#ifdef XFIG
          { int tl, th, xh;

            tl = B->ctgs[j+1].lft_end - ((hgh+var) - end);
            if (tl < B->ctgs[j].lft_end + B->ctgs[j].length)
              tl = B->ctgs[j].lft_end + B->ctgs[j].length;
            xh = MAP(end);
            tl = MAP(tl);
            // th = MAP(B->ctgs[j+1].lft_end);
	    th = B->ctgs[j+1].lft_end - ((low-var) -end);
            if (th > B->ctgs[j+1].lft_end)
              th = B->ctgs[j+1].lft_end ; /* tl -> th : did this fix anything? */
    	    th = MAP(th);
            fprintf(figfile,"2 3 0 1 5 5  98 0 35 4.06 1 0 0 0 0 5\n\t");
            if (gapped)
              { fprintf(figfile," %d %d",yl,zl);
                fprintf(figfile," %d %d",yl,zh);
                fprintf(figfile," %d %d",xh,th);
                fprintf(figfile," %d %d",xh,tl);
                fprintf(figfile," %d %d\n",yl,zl);
              }
            else
              { fprintf(figfile," %d %d",zl,yl);
                fprintf(figfile," %d %d",zh,yl);
                fprintf(figfile," %d %d",xh,tl);
                fprintf(figfile," %d %d",xh,th);
                fprintf(figfile," %d %d\n",zl,yl);
              }
          }
#endif
#ifdef DEBUG_ALIGN
          fprintf(stderr,"    Ray reaches A boundary\n");
#endif
        }
    }

  return (terminal);
}

  

void Free_Scaffold_List(Scaffold_List *SL)
{ Scaffold_List *s, *t;

  for (s = SL; s != NULL; s = t)
    { t = s->next;
      Free_Scaffold(s->scaffold);
      free(s);
    }
}

static int UnAccountedBits(int alow, int ahgh, int blow, int bhgh,
                           int i, int p, int j, int q, 
                           Scaffold *AF, Scaffold *BF)
{ int abg, abe;
  int u, v;
  int delta;

  delta = 0;
  if (alow < ahgh)
    { abg = alow; abe = ahgh; }
  else
    { abg = ahgh; abe = alow; }

#ifdef DEBUG_ANALYSIS
  fprintf(stderr,"  Checking gap [%d,%d] vs [%d,%d]\n",ahgh,alow,bhgh,blow);
  fprintf(stderr,"     Contigs   [%d,%d] vs [%d,%d]\n",p,i,q,j);
#endif

  for (u = p; u <= i; u++)
    { int b, e;
      int x, y;

      b = AF->ctgs[u].lft_end;
      e = b + AF->ctgs[u].length;
      if (abg > b) b = abg;
      if (abe < e) e = abe;
#ifdef DEBUG_ANALYSIS
      fprintf(stderr,"       X [%d,%d] = [%d,%d]\n",
             AF->ctgs[u].lft_end,AF->ctgs[u].lft_end+AF->ctgs[u].length,b,e);
#endif
      if (b >= e && abg != abe) continue;
      if (abg == abe)
        { x = bhgh;
          y = blow;
        }
      else
        { x = (int)(((blow - 1.*bhgh) / (alow - ahgh)) * (b - ahgh) + bhgh);
          y = (int)(((blow - 1.*bhgh) / (alow - ahgh)) * (e - ahgh) + bhgh);
        }
#ifdef DEBUG_ANALYSIS
      fprintf(stderr,"       yields [%d,%d]\n",x,y);
#endif
      if (x > y)
        { int m; m = x; x = y; y = m; }
      for (v = q; v <= j; v++)
        { int s, t;

          s = BF->ctgs[v].lft_end;
          t = s + BF->ctgs[v].length;
          if (s < x) s = x;
          if (t > y) t = y;
#ifdef DEBUG_ANALYSIS
          fprintf(stderr,"       X [%d,%d] = [%d,%d] of size %d\n",
                 BF->ctgs[v].lft_end,
                 BF->ctgs[v].lft_end+BF->ctgs[v].length,
                 s,t,t-s);
#endif
          if (s < t)
            delta += t-s;
        }
    }

  return (delta);
}

Scaffold_Overlap *Analyze_Overlap(Scaffold *AF, Scaffold *BF,
                                         Segment *seglist, int score, int comp)
{ int i, j;
  int diag, dmax = 0, dmin = 0;
  int firstime, lastime;
  int alow, ahgh = 0;
  int blow, bhgh = 0;
  int adelta, bdelta;
  int aone, bone;
  int alst, blst;
  int totdiff, totleng, segdiff;
  Segment *seg;

  Scaffold_Overlap *sovl;

  seg = seglist;
  aone = seg->a_contig;
  bone = seg->b_contig;

  segdiff = 0;
  totdiff = 0;
  totleng = 0;
  adelta = bdelta = 0;
  firstime = 1;
  lastime  = 0;
  i = j = 0;
  while (seg != NULL || ! lastime)
    { int p, q;

      alst = p = i;
      blst = q = j;
      if (seg == NULL)
        { i = AF->num_gaps;
          j = BF->num_gaps;
          if (AF->length - ahgh < BF->length - bhgh)
            { blow = bhgh + (AF->length-ahgh);
              alow = AF->length;
            }
          else
            { alow = ahgh + (BF->length-bhgh);
              blow = BF->length;
            }
        } 
      else
        { i = seg->a_contig;
          j = seg->b_contig;
          if (seg->overlap->begpos > 0)
            { alow = AF->ctgs[i].lft_end + seg->overlap->begpos;
              blow = BF->ctgs[j].lft_end;
            }
          else
            { alow = AF->ctgs[i].lft_end;
              blow = BF->ctgs[j].lft_end - seg->overlap->begpos;
            }
          totdiff += seg->overlap->diffs;
          totleng += seg->overlap->length;

          if (comp)
            { int m;
#ifdef DPC_BASED
              seg->overlap->comp = 1;
#endif
              m = BF->ctgs[j].length - seg->bhgh;
              seg->bhgh = BF->ctgs[j].length - seg->blow;
              seg->blow = m;
            }
        }

      if (firstime)
        { if (alow < blow)
            { ahgh = 0;
              bhgh = blow - alow;
            }
          else
            { ahgh = alow - blow;
              bhgh = 0;
            }
         }

      { int a_bits, b_bits;

        a_bits = UnAccountedBits(blow,bhgh,alow,ahgh,j,q,i,p,BF,AF);
        b_bits = UnAccountedBits(alow,ahgh,blow,bhgh,i,p,j,q,AF,BF);
        adelta += a_bits;
        bdelta += b_bits;
        if (a_bits > b_bits)
          segdiff += a_bits;
        else
          segdiff += b_bits;
      }

      if (seg != NULL)
        { if (seg->overlap->endpos > 0)
            { bhgh = BF->ctgs[j].lft_end
                   + (BF->ctgs[j].length - seg->overlap->endpos);
              ahgh = AF->ctgs[i].lft_end + AF->ctgs[i].length;
            }
          else
            { bhgh = BF->ctgs[j].lft_end + BF->ctgs[j].length;
              ahgh = AF->ctgs[i].lft_end
                   + (AF->ctgs[i].length + seg->overlap->endpos);
            }

          diag = ( (alow - blow) + (ahgh - bhgh) ) / 2;
          if (firstime || diag > dmax)
            dmax = diag;
          if (firstime || diag < dmin)
            dmin = diag;

#ifdef DEBUG_ANALYSIS
          fprintf(stderr,"  diag = %d\n",diag);
          fprintf(stderr,"  [%d,%d]  [%d,%d]\n",alow,ahgh,blow,bhgh);
#endif
        }

      firstime = 0;
      if (seg == NULL)
        lastime = 1;
      else
        seg = seg->next;
    }

  { int      i, diag;
    double   slack;
    Overlap *ovl;

    slack = 0.;
    for (i = aone; i < alst; i++)
      slack += AF->gaps[i].gap_var;
    for (i = bone; i < blst; i++)
      slack += BF->gaps[i].gap_var;

    diag = (dmin + dmax) / 2;

    ovl = (Overlap *) malloc(sizeof(Overlap));
    ovl->begpos = diag;
    ovl->endpos = diag - (AF->length - BF->length);
    ovl->diffs  = totdiff + segdiff;
    ovl->length = (AF->length + BF->length -
                        (abs(diag) + abs(ovl->endpos))) / 2;
    ovl->comp   = comp;
    //    ovl->aseq   = ovl->bseq = NULL;
    //    ovl->trace  = NULL;

    sovl = (Scaffold_Overlap *) malloc(sizeof(Scaffold_Overlap));
    sovl->score   = score;
    sovl->erate   = (1.*ovl->diffs) / ovl->length;
    sovl->ascaf   = AF;
    sovl->bscaf   = BF;
    sovl->D_delta = dmax-dmin;
    sovl->D_var   = (int)slack;
    sovl->A_delta = adelta;
    sovl->B_delta = bdelta;
    sovl->overlap = ovl;
    sovl->seglist = seglist;
  }

  return (sovl);
}
