
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
//#define DEBUG
//#define DEBUGS
/*
   Fragment simulator: generates simulated shotgun data.
 
   Author:            Gene Myers

   Date last Revised: October 12, 1998

   Desciption:
      Got rid of interactive input mode.  Input is now from a file
        parameter and not standard input.
      Got rid of old 4 parameter DNA sequence convention.
      Made it so number of elements can be expressed as a percent of basis.

   * Make it so that # of reads or # of repeats can be expressed as an X
     of genome.
   * Make it so that you can give a series of ramps for error ramp.
   * Make it so that part of repeats can be generated.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_UTL_rand.h"

char FragVersion[] = "$Revision: 1.1.1.1 $";

double drand48();
void   srand48();

void  *malloc();

#define MAX_INT 0x7FFFFFFF
#define FUDGE  0.9999

int   Seed;             /* Seed for random number generation        */

FILE *sfile;            /* Unit to read spec from. */

int   no_elem;          /* Don't print out sub-elements */
int   comments;         /* Don't print out any comments */
int   notation;         /* Emit fragment intervals */
int   printVersion;     /* Print Celsim Version in output */
int   uniform;          /* Use uniform mate pair distance distribution -- default is gaussian */
int   protoIO;          /* Nuke non-actgACTG characters.  If false, replace non-actgnACTGN with n */
//char  commentBuffer[256];


/* >>>> REPETITIVE ELEMENT DEFINITIONS <<<< */

/* An element is either of "type" BASIS or CNCAT, with or without the
   "GLOBAL" bit set.  If of type BASIS then the "len" and "prob" fields
   give the length range and symbol distribution for the basis sequence.
   "rlist" links together all of the element references in the given
   elements definition (in order).  "valstack" is used during sequence
   construction to maintain a stack of current instances for the given
   element.  "instances" points to a list (linked by the "ilist" field
   of the INSTANCE records) of instances of the given element used in
   building the subject strand. 
   Elements loaded from files have the FROMFILE flag set and are GLOBAL
*/

#define BASIS  0
#define CNCAT  1
#define GLOBAL 2
#define CONSTANT 4

typedef struct { struct refr *rlist;
                 int          type;
                 char *seq;            /* If type & CONSTANT, this is valid */
                 int          minlen, maxlen;
                 double       probA, probAC, probACG;
                 struct inst *valstack;
                 struct inst *instances;
               } DEFINITION;

/* Each reference is to an element ("name" in [0,25]).  The remaining
  fields record the orientation probability, mutation rate range,
  repetition factor range, and unique regeneration range, respectively.
  maxins = 0 if regeneration is not to take place.  */

typedef struct refr { struct refr *rlist;
                      int          name;
                      double       orient;
                      double       minerr, maxerr;
                      int          minrep, maxrep;
                      double       xminrp, xmaxrp;
                      int          minins, maxins;
                      double fracturePrefixp, fractureSuffixp, fractureRandomp;
                      double    fractureMinp, fractureMaxp;
                    } REFERENCE;

/* Instances of a given element are linked by their "ilist" field, and
   the "elist" field points to a list of the component instances contained
   in the instance.  The length and unique instance number of the instance
   are given in the fields "length" and "level".  "value" points at the
   string giving the instance sequence (only valid during construction
   save for that of the instance that becomes the subject strand).  "nextval"
   links together the stack of instances used during construction. */

typedef struct inst { struct comp *elist;
                      int          length;
                      int          level;
                      char        *value;
                      struct inst *nextval;
                      struct inst *ilist;
                    } INSTANCE;

/* "refer" and "level" identify the instance constituting this component.
   "length" is the components length (including mutations).  "oprand"
   points at the string value of the instance from which this component
   is derived (valid only during construction).  "forw" gives the orientation
   in which the underlying instance was incorporated, "erate" gives the
   mutation rate, and "errors" the number of mutations introduced.  Finally,
   "begpos" and "endpos" specify the beginning and ending indices within
   the containing instance at which this component was incorporated.  */

typedef enum {invalidFracture = -1, 
	      noFracture , 
	      prefixFracture, 
	      suffixFracture,
              randomFracture
} fractureType;

typedef struct comp { struct comp *elist;
                      int          length;
                      int          level;
                      struct refr *refer;
                      struct inst *source;
                      char        *oprand;
                      int          forw;
                      double       erate;
                      int          errors;
                      int          begpos, endpos;
                      fractureType fType; 
                      int          fbeg,fend;
                    } COMPONENT;
        
DEFINITION Relement[26];

int   Sequence;  /* Element defining string to be fragmented */
FILE *SeqFile;   /* if Sequence < 0 ==> File containing string to be frag'd */


/* >>>> BASIC UTILITIES <<<< */

#define SEGLEN    50    /* Output sequences, SEGLEN chars to a line */
#define DEL_FREQ .333   /* Deletion frequency for mutations */
#define INS_FREQ .333   /* Insertion frequency for mutations */

char SubChar[10][3] = { {'c', 'g', 't'},  /* What chars to substitute */
                       {'a', 'g', 't'},
                       {'a', 'c', 't'},
                       {'a', 'c', 'g'},
                       {'C', 'G', 'T'},
                       {'A', 'G', 'T'},
                       {'A', 'C', 'T'},
                       {'A', 'C', 'G'},
                       {'N', 'N', 'N'},  // n's don't get substituted
                       {'n', 'n', 'n'}};

int Decode[128];     /* Converts chars a,c,g,t to integers 0,1,2,3;
                        and A,C,G,T to integers 4,5,6,7.             */

char FlipChar[128];  /* Gives the Watson-Crick complement */
char TranChar[128];  /* The identity transform */

char *read_seq(FILE *input, int *len);

  /* Allocate space and exit if none */

char *ckalloc(int size)
{ char *m;
  if ((m = (char *) malloc(size)) == NULL)
    { fprintf(stdout,"Out of Memory (frag %d)\n",size);
      exit (2);
    }
  return (m);
}

  /* Generate a character with a coin toss */

char genchar(double pa,double pac,double pacg)
{ double rnd;

  rnd = drand48();
  if (rnd < pac)
    if (rnd < pa)
      return ('a');
    else
      return ('c');
  else
    if (rnd < pacg)
      return ('g');
    else 
      return ('t');
}

  /* Return a substitute for the input character. */

char substitute(int orig)
{ register int i;
  i = (int)((2. + FUDGE)*drand48());
  return (SubChar[Decode[orig]][i]);
}

 
/* >>>> SUBJECT SEQUENCE CONSTRUCTION <<<< */

/*  Place a string at s which is a mutated copy of the string t of
    length len.  The copy is in the direction given by fwd, mutated
    with err changes, ins of which are insertions, and del of which
    are deletions.     */

char *mutate(char *s,char *t,int len,int fwd,int err,int ins,int del)
{ register double x;
  register char *tr;

  len += (del - ins);
  del += ins;
  tr   = TranChar;
  if (!fwd)
    { fwd = -1;
      t  += (len-1);
      tr  = FlipChar;
    }

#ifdef DEBUGS
  fprintf(stderr,"mutate = '");
#endif
  while (len > 0 || err > 0)
    if ((x = err*drand48()) < ins)
      { if ((len+1)*drand48() < err)
          { *s++ = genchar(.25,.50,.75);
            del -= 1;
            ins -= 1;
            err -= 1;
#ifdef DEBUGS
            printf("!%c",(s[-1]-'a')+'A');
#endif
          }
        else
          { *s++ = tr[(int) *t];
            t += fwd;
            len -= 1;
#ifdef DEBUGS
            printf("%c",s[-1]);
#endif
          }
      }
    else if (x < del)
      { if (len*drand48() < err)
          { t += fwd;
            del -= 1;
            err -= 1;
#ifdef DEBUGS
            printf(".%c",(t[-fwd]-'a')+'A');
#endif
          }
        else
          { *s++ = tr[(int) *t];
            t += fwd;
#ifdef DEBUGS
            printf("%c",s[-1]);
#endif
          }
        len -= 1;
      }
    else
      { if (len*drand48() < err)
          { *s++ = substitute(tr[(int) *t]);
            t += fwd;
            err -= 1;
#ifdef DEBUGS
            printf("%c",(s[-1]-'a')+'A');
#endif
          }
        else
          { *s++ = tr[(int) *t];
            t += fwd;
#ifdef DEBUGS
            printf("%c",s[-1]);
#endif
          }
        len -= 1;
      }
#ifdef DEBUGS
  printf("'\n");
#endif
  return (s);
}

/* Build a string for version lev of element e.  This is done by
   recursively building strings for each reference in e's definition
   (if that version has not yet been built), then determining the
   orientation, mutation level, and number of instances for each
   reference (so the size of e's string is known), and then building
   the string for e from these instance descriptions.   */

static int level;

int COMPARE(COMPONENT **a,COMPONENT **b)
{ return ((*a)->forw - (*b)->forw); }

COMPONENT *instantiateReference(REFERENCE *r,INSTANCE *v,INSTANCE *i){
  double rv;
  COMPONENT *x = NULL;
  int actualLength;
  double diceThrow = drand48();
  double indenom;
  
  x = (COMPONENT *) ckalloc(sizeof(COMPONENT));
  x->elist  = NULL;
  x->level  = v->level;
  x->refer  = r;
  x->source = v;
  x->fType = noFracture; /* No Fracture */
  /*  x->oprand = v->value; now a function of fracturing */
  x->erate  = r->minerr + drand48()*(r->maxerr - r->minerr);
  /* We need to compute the actual length of this component, including the 
     effect of fractures.
  */
  
  indenom = 1. - INS_FREQ;
  
  if(r->fracturePrefixp + 
     r->fractureSuffixp + 
     r->fractureRandomp > diceThrow)
  { 
    /* we'll decide what type later */
    
    int offset;
    
    /* Uniform length distribution */
    
    actualLength = (int)
      (v->length  + FUDGE)* 
      (r->fractureMinp + (r->fractureMaxp - r->fractureMinp) * drand48() );
    
#ifdef DEBUGS
    fprintf(stderr,"actualLength = %d v->length = %d fractureMinp = %g fractureMaxp = %g\n",
            actualLength, v->length, r->fractureMinp, r->fractureMaxp);
#endif
    if(diceThrow < r->fracturePrefixp)  /* a prefix fracture */
    {
      x->fType = prefixFracture;
      offset = 0;
      x->fbeg = 0;
#ifdef DEBUGS
      printf("Prefix Fracture of length %d at offset %d\n",
             actualLength,offset);
#endif
    }
    else
    {
      diceThrow -= r->fracturePrefixp;
      if (diceThrow < r->fractureSuffixp)
      {
        /* a suffix fracture */
        x->fType = suffixFracture;
        offset = v->length - actualLength;
        x->fbeg = offset;
#ifdef DEBUGS
        printf("Suffix Fracture of length %d at offset %d\n",
               actualLength,offset);
#endif
      }
      else
      {    /* a random fracture */
        x->fType = randomFracture;
        /* range of offsets for a subfragment of length
           actualLength from a fragment of length v->length */
        offset = (v->length - actualLength + FUDGE) * drand48();
        x->fbeg = offset;
#ifdef DEBUGS
        printf("Random Fracture of length %d at offset %d\n",
               actualLength,offset);
#endif
      }
    }
    x->oprand = v->value + offset; /* get point to offset string */
    
  }
  else
  {
    actualLength = v->length;
    x->oprand = v->value; /* get point to whole string */
  }
  
  /* Sanity */
  assert(actualLength <= v->length &&
         actualLength >= 0);
  
  // This is some kind of trick with float/int roundoff/truncation
  x->errors = rv = actualLength * x->erate;
  if (rv - x->errors > drand48())
  {
    x->errors += 1;
  }
  x->begpos = rv = x->errors * INS_FREQ;
  if (rv - x->begpos > drand48())
  {
    x->begpos += 1;
  }
  if (indenom < 1.e-20)
  {
    x->endpos = 0;
  }
  else
  {
    rv = (x->errors - x->begpos)*DEL_FREQ/indenom;
    x->endpos = rv;
    if (rv - x->endpos > drand48())
      x->endpos += 1;
  }
  x->length = actualLength + x->begpos - x->endpos;
  x->fend   = x->fbeg + actualLength;
  return x;
}

void build_trav(DEFINITION *e, int lev)
{ REFERENCE  *r = NULL;
  DEFINITION *f = NULL;
  INSTANCE   *i = NULL, *u = NULL, *v = NULL;
  COMPONENT  *x = NULL, *y = NULL, **z = NULL;
  char *s;
  int n, l;
  int reps, len, num, ins, avl;
  double xrp;
  INSTANCE uhead;
  int *insters;
  int refLen; /* total length, so far, of a particular reference */

  ins = 0;
  for (r = e->rlist; r != NULL; r = r->rlist)
    if (r->maxins > 0) ins += 1;

  if (ins > 0)
    insters = (int *) malloc(sizeof(int)*ins);
  else
    insters = NULL;

  ins = 0;
  u = &uhead;
  for (r = e->rlist; r != NULL; r = r->rlist)
    { f = Relement + r->name;
      if (r->maxins > 0)
        { len = r->minins + drand48()*(r->maxins - r->minins + FUDGE);
          insters[ins++] = len;
          while (len-- > 0)
            { build_trav(f,l = ++level);
              u->nextval  = f->valstack;
              u = u->nextval;
              f->valstack = u->nextval;
              for (n = 0; n < 26; n++)
                { v = Relement[n].valstack;
                  if (v != NULL && v->level == l)
                    { if (no_elem) free(v->value);
                      Relement[n].valstack = v->nextval;
                    }
                }
            }
        }
      else if (f->valstack == NULL ||
               (f->valstack->level != lev && !(f->type & GLOBAL)))
        { if (f->type & GLOBAL)
            build_trav(f,0);
          else
            build_trav(f,lev);
        }
    }

#ifdef DEBUG
  fprintf(stderr," Build_trav: %c lev = %d\n",'A' + (e-Relement),lev); 
  fprintf(stderr," Build_trav: e->valstack = 0x%p \n", e->valstack);
#endif
  i = (INSTANCE *) ckalloc(sizeof(INSTANCE));
  i->elist = NULL;
  i->level     = lev;
  i->ilist     = e->instances;
  e->instances = i;
  i->nextval   = e->valstack;
  e->valstack  = i;
  if ((e->type & CNCAT) == BASIS)
    i->length = e->minlen + drand48()*(e->maxlen - e->minlen);

  len = 0;
  num = 0;
  ins = 0;
  u = uhead.nextval;
  y = (COMPONENT *) i;
  for (r = e->rlist; r != NULL; r = r->rlist){ 
    f = Relement + r->name;
      if (r->maxins > 0)
        n = insters[ins++];
      else
        n = 1;
      while (n-- > 0)
        { if (r->maxins > 0)
            { v = u;
              u = u->nextval;
            }
          else
            v = f->valstack;


	if (r->minrep >= 0){   /* Fixed number of repeats */
            reps = r->minrep + drand48()*(r->maxrep - r->minrep + FUDGE);
	    refLen = 0;
	    while (reps-- > 0)
            { 
	      x = instantiateReference(r,v,i);
              y = y->elist  = x;
	      refLen += x->length; /* Copies of this reference */
              len += x->length;    /* Total length of this instance */
              num += 1;
#ifdef DEBUG
              fprintf(stderr,"\t%c.%d: o=%d e=%g(%d,%d,%d) length = %d value = 0x%x\n\t%s\n",
                      x->refer->name+'A',x->level,x->forw,x->erate,x->errors,
                      x->begpos,x->endpos,x->length,x->oprand, x->oprand); 
#endif
            }
	} else{ /* Probabilistic number of repeats */
	      int desiredTotalLength;
	      double totalFracturep = ( r->fracturePrefixp +
                                        r->fractureSuffixp +
					r->fractureRandomp);
	      xrp = r->xminrp + drand48()*(r->xmaxrp - r->xminrp);
	      /* Take the fracture probabilities into account in computing the 
		 average length of the instances of this elment */
	      avl =  totalFracturep * (r->fractureMaxp - r->fractureMinp)/2.0 + 
		(1-totalFracturep) * v->length;
	      /* avl now reflects the average length of unmutated instances of
		 the current element */
              avl = avl * (1. +
                     (INS_FREQ-DEL_FREQ) * (r->minerr + r->maxerr) / 2.);
	      desiredTotalLength = (int)(xrp * i->length);
	      refLen = 0;
	      reps = 0;

#ifdef DEBUGS
              fprintf(stderr,"*** Shooting for total length %d, avl %d nominal length = %d\n",
		     desiredTotalLength, avl, v->length);
#endif

	      /* Keep adding until we are one short of, on average, exceeding the desired total length */
	      while(refLen + avl < desiredTotalLength){
#ifdef DEBUGS
		fprintf(stderr,"refLen =%d avl = %d desired = %d\n",
		       refLen, avl, desiredTotalLength);
#endif
		x = instantiateReference(r,v,i);
		y = y->elist  = x;
		refLen += x->length; /* Copies of this reference */
		len += x->length;    /* Total length of this instance */
		num += 1;
		reps++;
#ifdef DEBUGS
              fprintf(stderr,"\t%c.%d: o=%d e=%g(%d,%d,%d) length = %d value = 0x%x %s\n",
                      x->refer->name+'A',x->level,x->forw,x->erate,x->errors,
                      x->begpos,x->endpos,x->length,x->oprand, x->oprand);
#endif

	      }
	      /* Now, probabilistically add another.  Be conservative, since we don't want
	         to go over.  We use the nominal length of the reference to determine the probability,
	         rather than the expected length.  If there are fractures, this may be an
	         overestimate */
#ifdef DEBUGS
	      fprintf(stderr,"*** desired = %d refLen = %d v->length = %d\n",
		     desiredTotalLength, refLen, v->length);
#endif

	      {
		double diceThrow = drand48();
		double pAdd = (double)(desiredTotalLength - refLen)/(double)(v->length);

		if(pAdd > diceThrow){
		x = instantiateReference(r,v,i);
		y = y->elist  = x;
		refLen += x->length; /* Copies of this reference */
		len += x->length;    /* Total length of this instance */
		num += 1;
		reps++;
#ifdef DEBUG
              fprintf(stderr,"\t%c.%d: o=%d e=%g(%d,%d,%d) length = %d value = 0x%x\n",
                      x->refer->name+'A',x->level,x->forw,x->erate,x->errors,
                      x->begpos,x->endpos,x->length,x->oprand); fflush(stdout);
#endif

		}
	      }

#ifdef DEBUG
	      fprintf(stderr,"Added %d reps of total Length %d for target length %d\n",
		     reps, refLen, desiredTotalLength);
#endif
	      
	}
	}
  }

  y->elist = NULL;
  if ((e->type & CNCAT) == CNCAT){
    i->length = len;
  }

  i->value = (char *) ckalloc(i->length+1);


  s = i->value;
  if ((e->type & CNCAT) == CNCAT) {  /* CNCAT */
    for (x = i->elist; x != NULL; x = x->elist)
      { n = s - i->value;
        x->forw = (drand48() <= x->refer->orient);
        s = mutate(s,x->oprand,x->length,x->forw,
                     x->errors,x->begpos,x->endpos);
        x->begpos = n;
        x->endpos = s - i->value;
      }
  }  else { // if(!(e->type & CONSTANT))      /* BASIS */
      if (num > 0){
         z = (COMPONENT **) ckalloc(sizeof(COMPONENT *)*num);
          n = 0;
          for (x = i->elist; x != NULL; x = x->elist)
            { z[n++] = x;
              x->forw = drand48()*MAX_INT;
            }
          qsort(z,
		num,
		sizeof(COMPONENT *),
		(int(*)(const void*, const void*))COMPARE);
          y = (COMPONENT *) i;
          for (n = 0; n < num; n++)
            y = y->elist = z[n];
          y->elist = NULL;
          free(z);
      }
#ifdef DEBUG
      else{
	fprintf(stderr,"* num == 0\n");
      }
#endif
      len = i->length - len;
      x   = i->elist;
      while (len > 0 || num > 0){
        if (drand48() < num/(len+1.)){ 
	    assert(x != NULL);
	    n = s - i->value;
            x->forw = (drand48() <= x->refer->orient);
            s = mutate(s,x->oprand,x->length,x->forw,
                         x->errors,x->begpos,x->endpos);
            x->begpos = n;
            x->endpos = s - i->value;
            x = x->elist;
            num -= 1;
	} else if(e->type & CONSTANT){ /* Simply copy it into the the instance */
	    strcpy(s,e->seq);
	    s += i->length;
	    len -= i->length;
	}else{                       /* Generate a probabilistic character */
	    *s++ = genchar(e->probA,e->probAC,e->probACG);
            len -= 1;
          }
      }
  } /* end BASIS AND NOT CONSTANT */

  *s = '\0';
#ifdef DEBUGS
  fprintf(stderr,"Val = '%s'\n",i->value);
#endif

  if (no_elem)
    for (r = e->rlist; r != NULL; r = r->rlist)
      if (r->maxins > 0)
        { u = uhead.nextval;
          uhead.nextval = u->nextval;
          free(u->value);
        }

  if (insters != NULL) free(insters);
}

void report_trav(INSTANCE *i, int lev, int type, int pos,
                              int rev, int fb, int fe, int ftype)
{ COMPONENT *x = NULL;

  if (rev)
    { COMPONENT *y, *t;

      y = NULL;
      for (x = i->elist; x != NULL; x = t)
        { t = x->elist;
          x->elist = y;
          y = x; 
        }
      i->elist = y;
    }

  for (x = i->elist; x != NULL; x = x->elist)
    { int xb=-1, xe=-1, tb=-1, te=-1;
      int ctype;

      ctype = ftype;
      if (x->fType == randomFracture)
        ctype |= 0x3;
      else if (x->fType == prefixFracture)
        ctype |= 0x2;
      else if (x->fType == suffixFracture)
        ctype |= 0x1;

      xb = x->fbeg;
      xe = x->fend;
      if (x->begpos < fb)
        { tb = fb;
          xb += (fb - x->begpos);
        }
      else
        tb = x->begpos;
      if (x->endpos > fe)
        { te = fe;
          xe -= (x->endpos - fe);
        }
      else
        te = x->endpos;
      if (tb >= te) continue;

      printf("# %*s ",lev,"");
      if (type & 1)
        printf("= ");
      else
        printf("> ");
      printf("%c.%d",x->refer->name+'A',x->level);
      if (rev)
        { if (x->forw)
            printf("r");
          else
            printf("f");
          te = (pos+fe)-te;
          printf(" at %d-%d",te,(pos+fe)-tb);
          tb = te;
        }
      else
        { if (x->forw)
            printf("f");
          else
            printf("r");
          tb = pos+(tb-fb);
          printf(" at %d-%d",tb,pos+(te-fb));
        }
      printf(", mutated %.2f",x->erate);

      if (xe != x->length && (ctype & 0x2) != 0)
        if (xb != 0 && (ctype & 0x1) != 0)
	  printf(", random fracture [%d,%d]", xb, xe);
        else
	  printf(", suffix fracture [%d,%d]", xb, xe);
      else if (xb != 0 && (ctype & 0x1) != 0)
	printf(", prefix fracture [%d,%d]", xb, xe);
      printf(".\n");

      report_trav(x->source,lev+2,Relement[x->refer->name].type,tb,
                  x->forw?rev:1-rev,xb,xe,ctype);
    }

  if (rev)
    { COMPONENT *y, *t;

      y = NULL;
      for (x = i->elist; x != NULL; x = t)
        { t = x->elist;
          x->elist = y;
          y = x; 
        }
      i->elist = y;
    }
}
      

/* Build the subject dna strand. */

void build_seq(void)
{ int n, k, len;
  INSTANCE  *i;

  /* Clear out the valstack and instances for each element */
  for (n = 0; n < 26; n++)
    { Relement[n].valstack  = NULL;
      Relement[n].instances = NULL;
      /* DO NOT clear out the rlist, since this is what defines the
	 sequence grammar 
         *** NO *** Relement[n].rlist = NULL; *** NO ***
       */
    }
  level = 0;
  build_trav(Relement+Sequence,0);

  i = Relement[Sequence].instances;
  if (comments)
    { printf("#  %c.%d (%d)\n",Sequence+'A',i->level,i->length);
      report_trav(i,0,Relement[Sequence].type,0,0,0,i->length,0);
      printf("#\n");
    }

  if (!no_elem && comments)
    { for (n = 0; n < 26; n++)
        for (i = Relement[n].instances; i != NULL; i = i->ilist)
          if (n != Sequence)
            { len = i->length;
              printf("#  %c.%d (%d)\n",n+'A',i->level,len);
              for (k = 0; k < len; k += SEGLEN)
                if (len-k < SEGLEN)
                  printf("#     %.*s\n",len-k,i->value+k);
                else
                  printf("#     %.*s\n",SEGLEN,i->value+k);
            }
      printf("#\n");
    }
}


/* >>>> FRAGMENT GENERATION <<<< */

/*  Fragment generation parameters. */

int     NumFrag;          /* the # of fragments to be sampled         */ 
int     MinLen, MaxLen;   /* the length range for fragments           */
double  Orient;           /* the probability that the given fragment is
		             in the forward direction                 */

double  BegErr, EndErr;   /* the beginning and end of the error int.  */
double  ProbI, ProbID;    /* the probabilities of insertions and
  		             deletions                                */
double  SingleFract;      /* fraction of single (unpaired) reads */
int     MinIns, MaxIns;   /* the length range for for vector inserts  */
double  PairError;        /* fraction of erroneous pairings */

/* Current fragment definition. */

int   StartPos;     /* the start index of the fragment */
int   FragLen;      /* the length of the fragment */
int   OriginalFragLen;      /* the length of the fragment before edits */
int   FragForw;     /* true if the frag is to be in the forward orientation */
char *Frag;         /* the final error-ladden fragment sequence */
char *Qual;         /* the associated quality values */
int   Fragtop;      /* current length of Frag array */

/* Current insert definition. */

int   StartIns;    /* the start index of the insert */
int   InsLen;      /* the length of the insert */
int   InsForw;     /* true if the insert is to be in the forward orientation */

/* Quality value parameters */
int   IntroduceSequenceError = FALSE; /* Used to handle special case of NO sequence error */
float LevelProb = .70;   /* percent of read covered by "level" QV */
float ValleyProb = .7;   /* percent of 'featured' data which start valleys */
int   ValleyDepth = 22;  /* depth of valleys */
int   ValleyBottom;      /* min QV at valley bottom (varies with LevelQV) */
int   HillHeight = 10;   /* height of hills */
int   HillPeak;          /* max QV at peak of hills (varies with LevelQV) */
int   LevelQV;           /* QV of "level" portion ( 70 % of frag ) */
float QVfit1 = -10.0;    /* coef of QV fit to logprob:    */
float QVfit2 =  7.44;    /*    LevelQV = QVfit1*log10(prob)+QVfit2 (+'0') */
char Quality_File_Name[FILENAME_MAX];
FILE *QualityFile;

/* Sample a fragment from the dna strand. */

void newfrag(int dnalen)
{ 
#if 0
  FragLen  = MinLen + (MaxLen - MinLen + FUDGE)*drand48();
  StartPos = (dnalen - FragLen + FUDGE)*drand48();
#endif

  FragLen = GetRand_AS(MinLen, MaxLen, uniform);

  StartPos = GetRand_AS(0, dnalen - FragLen, TRUE);

  if (FragLen >= dnalen)
    { FragLen  = dnalen;
      StartPos = 0;
    }

  if (drand48() < Orient)
    FragForw = 1;
  else
    FragForw = 0;
}

void firstfrag(int dnalen)
{
#if 0
  InsLen  = MinIns + (MaxIns - MinIns + FUDGE)*drand48();
  StartIns = (dnalen - InsLen + FUDGE)*drand48();
#endif
  InsLen = GetRand_AS(MinIns, MaxIns, uniform);
  StartIns = GetRand_AS(0, dnalen - InsLen, TRUE);

  if (InsLen >= dnalen)
    { InsLen  = dnalen;
      StartIns= 0;
    }

  if (drand48() < Orient)
    InsForw = 1;
  else
    InsForw = 0;

  FragLen  = MinLen + (MaxLen - MinLen + FUDGE)*drand48();
  if (FragLen > InsLen) FragLen = InsLen;
  if (InsForw)
    { StartPos = StartIns;
      FragForw = 1;
    }
  else
    { StartPos = StartIns + InsLen - FragLen;
      FragForw = 0;
    }
}

void secondfrag(int dnalen)
{ 
#if 0
  FragLen  = MinLen + (MaxLen - MinLen + FUDGE)*drand48();
#endif
  FragLen = GetRand_AS(MinLen, MaxLen, uniform);

  if (FragLen > InsLen) FragLen = InsLen;
  if (InsForw)
    { StartPos = StartIns + InsLen - FragLen;
      FragForw = 0;
    }
  else
    { StartPos = StartIns;
      FragForw = 1;
    }
}

/* Traverse the input fragment, introducing edits into it. */

void editfrag(char *dna)
{ double rnd, be, de;
  int    ucap;
  register int i, j, t, m;
  register char *xfer;
  static   int call1 = 1;
  int isvalley, valleystat;
  float qvprob;

  static char buf1[SEGLEN+1], buf2[SEGLEN+1];

  if (call1)
    { call1 = 0;
      Fragtop = 2*MaxLen;
      Frag = (char *) ckalloc(Fragtop);
      Qual = (char *) ckalloc(Fragtop);
    }

  if (FragForw)
    { be = BegErr; de = EndErr-be; }
  else
    { be = EndErr; de = BegErr-be;
    }

  dna += StartPos;
  ucap = isupper(*dna);
  i = j = t = 0;
  isvalley = 0;
  valleystat = 0;
  while (1)
    { /* determine the quality value */
      if (IntroduceSequenceError == FALSE || drand48() < LevelProb) { /* remain at level quality */
         Qual[j] = LevelQV;
      } else if (drand48() < .66 ) { /* insert a valley */
         Qual[j] = (int) GetRand_AS(ValleyBottom,LevelQV,TRUE);
      } else {                         /* insert a hill */
         Qual[j] = (int) GetRand_AS(LevelQV,HillPeak,TRUE);
      }
      qvprob = pow(10,(float) (Qual[j]-'0')/(-10));
      //if (drand48() < be + de*i/FragLen)  /* Make an edit. */
      if (IntroduceSequenceError && drand48() < qvprob)  /* Make an edit. */
	{ rnd = drand48();
	  if (rnd < ProbI)                /* Insert a character. */
            { buf1[t] = '-';
              if (ucap)
                buf2[t++] = Frag[j++] = toupper(genchar(.25,.50,.75));
              else
                buf2[t++] = Frag[j++] = genchar(.25,.50,.75);
            }
          else if (i >= FragLen)            /* At end of string => quit */
            break;
	  else if (rnd < ProbID)     /* Delete next character.  */
            { buf1[t] = dna[i++];
              buf2[t++] = '-';
            }
	  else                            /* Substitute next character. */
            { buf1[t] = dna[i];
              buf2[t++] = Frag[j++] = substitute(dna[i++]);
            }
	}
      else if (i >= FragLen)               /* At end of string => quit */
        break;
      else                                 /* Character unedited */
        { buf1[t] = Frag[j++] = dna[i++];
          buf2[t++] = ' ';
        }
      if (j >= Fragtop)
        { xfer = (char *) ckalloc(Fragtop+MaxLen);
          for (j = 0; j < Fragtop; j++)
            xfer[j] = Frag[j];
          free(Frag);
          Frag = xfer;
          Fragtop += MaxLen;
        }
      if (t >= SEGLEN)
        { buf1[t] = buf2[t] = '\0';
          if (comments)
            { printf("#  %s\n",buf1);
              printf("#  %s\n",buf2);
              printf("#\n");
            }
          t = 0;
        }
    }
  if (t > 0)
    { buf1[t] = buf2[t] = '\0';
      if (comments)
        { printf("#  %s\n",buf1);
          printf("#  %s\n",buf2);
          printf("#\n");
        }
    }

  if (!FragForw)
    { i = 0;
      t = j-1;
      while (i <= t)
        { m = Qual[i];
          Qual[i] = Qual[t];
          Qual[t] = m;
          m = Frag[i];
          Frag[i++] = FlipChar[(int) Frag[t]];
          Frag[t--] = FlipChar[m];
        }
    }

  OriginalFragLen = FragLen;	// save position in DNA strand
  FragLen = j;
}

void outputfrag(int index, char c, int isFalseMate)
{ register int i;
/* output sequence */
  printf("> %d%c Interval: [%d,%d] %s %s\n",
	 index,c,
	 StartPos,StartPos+OriginalFragLen,
	  (FragForw?"":" REVERSED "),
	  (isFalseMate?" False Mate ":""));

  for(i = 0; i < FragLen; i++)
    if(!isalpha(Frag[i]) || isspace(Frag[i]))
      fprintf(stderr,"* Non alpha character 0x%d at pos %d\n",
	      Frag[i], i);

  for (i = 0; i < FragLen; i += SEGLEN)
    if (FragLen-i < SEGLEN)
      printf("%.*s\n",FragLen-i,Frag+i);
    else
      printf("%.*s\n",SEGLEN,Frag+i);

  fflush(stdout);
/* output quality */
  fprintf(QualityFile,"> %d%c Interval: [%d,%d] %s %s (Quality Values)\n",
	 index,c,
	 StartPos,StartPos+FragLen,
	  (FragForw?"":" REVERSED "),
	  (isFalseMate?" False Mate ":""));

  for(i = 0; i < FragLen; i++)
    if(isspace(Qual[i]))
      fprintf(stderr,"* Non alpha character 0x%d at pos %d\n",
	      Qual[i], i);

  for (i = 0; i < FragLen; i += SEGLEN)
    if (FragLen-i < SEGLEN)
      fprintf(QualityFile,"%.*s\n",FragLen-i,Qual+i);
    else
      fprintf(QualityFile,"%.*s\n",SEGLEN,Qual+i);
  fflush(QualityFile);
}


/* >>>> INPUTING AND PARSING A SIMULATION SPECIFICATION <<<< */

void swallowComment(void)
{  int c;
/* start of comment */
    while(1){
      c = fgetc(sfile);
	if(c == EOF){
	  fprintf(stderr,"Unexpected EOF in comment\n");
	  exit(1);
	}
      if(c == '\n')
	return;
    }
}


/* Nextsymbol gets the next non blank character.  Any error terminates
   execution.  To make error messages informative in that context, we
   adopt the policy of writing out object descriptions as soon as they
   are correctly input.   */

int nextsymbol(void)
{ register int c;
/* first, suck up spaces */
  do c = fgetc(sfile); while (isspace(c));
  if(c == '#'){
    swallowComment();
    return nextsymbol(); /* return the nextsymbol, possibly skipping whitespace and comments */
  }
  return c;
}

/* Get a postive natural number from the sfile.   The input parameter
   query is true if getnatural is to handle the entire input transaction,
   including re-entry prompts.  Query is always false when in batch mode. */

int getnatural(char *label)
{ int datum;
 
  if (fscanf(sfile," %d",&datum) != 1 || datum <= 0)
    {      fprintf(stderr,"\nFrag Error:\n\t*** %s must be a positive int\n",label);
     fflush(stderr);
      exit(1);
    }
  return(datum);
}

char *getstring(char *label, char *buffer, int length)
{ char tbuffer[2048];
 
  if (fscanf(sfile," %s", tbuffer) != 1 || strlen(tbuffer) >= length)
    { fprintf(stdout,"\nFrag Error:\n\t*** %s must be a string of length less than %d\n",label,length);
      exit(1);
    }
  strcpy(buffer,tbuffer);
  return(buffer);
}

/* Get a probability from the sfile. */

double getprobability(char *label)
{
  double datum;
  if (fscanf(sfile," %lf",&datum) != 1 || datum < 0. || datum > 1.)
    { fprintf(stderr,"\nFrag Error:\n\t*** %s must be a real number in [0,1]\n",label);
      exit(1);
    }
  return(datum);
}

/* Get a number from the sfile. */

int getfactor(char *label, double *real)
{ int value=-1;
  double xvalue;
  register int c;
  int cnt1, cnt2, expo;
  unsigned int frac;

  c = nextsymbol();
  ungetc(c,sfile);
  value = 0;
  if (fscanf(sfile,"%d",&value) != 1 && c != '.')
    { fprintf(stderr,
              "\nFrag Error:\n\t*** Expecting an integer or real number\n");
      exit(1);
    }
  c = getc(sfile);
  if (c == '.' || c == 'e' || c == 'E')
    { if (c == '.')
        { c = getc(sfile);
          ungetc(c,sfile);
          if (isdigit(c))
            fscanf(sfile,"%n%u%n",&cnt1,&frac,&cnt2);
          else
            cnt1 = cnt2 = frac = 0;
          if (cnt1 == cnt2 || frac == 0) 
            xvalue = 0.;
          else
            { while (frac % 10 == 0)
                { frac /= 10;
                  cnt1 += 1;
                }
              xvalue = frac;
              while (cnt1 < cnt2)
                { xvalue /= 10.;
                  cnt1 += 1;
                }
	    }
          xvalue += value;
        }
      else
        { ungetc(c,sfile);
          xvalue = value;
        }
      c = getc(sfile);
      if (c == 'e' || c == 'E')
        { if (fscanf(sfile,"%d",&expo) != 1)
            ungetc(c,sfile);
          else
            xvalue *= pow(10.,1.*expo);
        }
      else
        ungetc(c,sfile);
      if (xvalue < 0. || xvalue > 1.)
        { fprintf(stdout,"\n\t*** %s %% must be a real ",label);
          fprintf(stdout,"number in [0,1]\n");
          exit(1);
        }
      *real = xvalue;
      return (-1);
    }
  else
    { ungetc(c,sfile);
      if (value < 0)
        { fprintf(stdout,
              "\n\t*** %s must be a non-negative int\n",label);
          exit(1);
        }
      return (value);
    }
}

/* Read the string for an element definition from a file */

void getElementFromFile(DEFINITION *e)
{ 
  char tfilename[1024];
  char *tstring;
  int len;
  FILE *efile;

  tstring = getstring("Filename", tfilename, 1024);
  if((efile = fopen(tstring,"r")) == NULL){
    fprintf(stderr,"** Couldn't open file %s\n", tstring);
    exit(1);
  }
  e->seq = read_seq(efile,&len);
#ifdef DEBUG
  fprintf(stderr,"*** read_seq returning(" F_SIZE_T "):\n%s\n\n",
	  strlen(e->seq),e->seq);
#endif
  fclose(efile);
  e->minlen  = len;   /* Establish record fields */
  e->maxlen  = len;
  e->probA   = 0.25;
  e->probAC  = 0.50;
  e->probACG = 0.75;
}

/* Read the string for an element definition from a inline, quoted string */

void getElementFromInlineString(DEFINITION *e)
{ 
  char *buffer;
  char *tstring;
  int len;
  char c;
  int charsDeleted = 0;
  int charsSubstituted = 0;
  int size = 2048;
  buffer = (char*) malloc(size);
  tstring = buffer;
  for(c = nextsymbol(); 
      c != '"' && c != ';' && c != EOF;
      c = nextsymbol()){
    // fprintf(stderr,"* Read char %c\n", c);
    if(!isspace(c))
       switch(tolower(c)){
	  case 'a':
	  case 'c':
	  case 't':
	  case 'g':
	  case '\n':   // these get deleted below
	    *tstring++ = c;
	    break;
	  case 'm':
	  case 'r':
	  case 'w':
	  case 's':
	  case 'y':
	  case 'k':
	  case 'v':
	  case 'h':
	  case 'd':
	  case 'b':
	    if(protoIO == 0)
	      charsSubstituted++;
	    /* Intentional fallthrough */
	  case 'n':  // These are tolerated
	    if(protoIO == 0){
	      *tstring++ = 'n';  // We don't know how to complement these yet
	    }else{
	      charsDeleted++;
	    }
	    break;
	  default:
	    charsDeleted++;
	    break;
	  }

    if(tstring - buffer >= size){
      int oldsize = tstring - buffer;
      size *= 2;
      buffer = (char*)realloc(buffer,size);
      assert(NULL != buffer);
      tstring = buffer + oldsize;
    }
  }
    if(c != '"'){

    fprintf(stderr,"* Error reading inline string from input %s\n",
	    buffer);
    }

    *tstring++ = '\0';

    if(charsDeleted > 0 ||
       charsSubstituted > 0)
      fprintf(stderr,"* Read inline string substituted %d and deleted %d chars\n",
	      charsDeleted, charsSubstituted);

  e->seq = buffer;
  len = strlen(e->seq);
#ifdef DEBUG
  fprintf(stderr,"*** Read (" F_SIZE_T "):\n%s\n\n",
	  strlen(e->seq),e->seq);
#endif
  e->minlen  = len;   /* Establish record fields */
  e->maxlen  = len;
  e->probA   = 0.25;
  e->probAC  = 0.50;
  e->probACG = 0.75;
}

void getbasis(DEFINITION *e)
{ register int c;
  int minl, maxl;
  double pa, pac, pacg;

/* Get size range */

  minl = getnatural("(Min) Element length");
  c = nextsymbol();
  if (c == '-')
    { maxl = getnatural("Max Element length");
      c = nextsymbol();
    }
  else
    maxl = minl;
  if (minl > maxl)
    { fprintf(stdout,"\n\t*** Min length (%d) > Max length (%d)\n",minl,maxl);
      exit (1);
    }

/* Get symbol bias probabilities (if specified) */

  if (c == 'p')
    { if ((c = nextsymbol()) != '(') ungetc(c,sfile);
      pa = getprobability("Odds of A");
      if ((c = nextsymbol()) != ',') ungetc(c,sfile);
      pac = getprobability("Odds of C");
      if ((c = nextsymbol()) != ',') ungetc(c,sfile);
      pacg = getprobability("Odds of G");
      if ((c = nextsymbol()) != ')') ungetc(c,sfile);
      if (pacg > 1.)
        { fprintf(stderr,"\nFragError:\n\t*** Sum (%.2f) of A|C|G odds > 1\n",pacg);
          exit (1);
        }
    }
  else
    { pa = pac = pacg = .25;
      ungetc(c,sfile);
    }

  e->minlen  = minl;   /* Establish record fields */
  e->maxlen  = maxl;
  e->probA   = pa;
  e->probAC  = pac += pa;
  e->probACG = pacg += pac;
}

/* Build and return another reference record for element e's definition. */

REFERENCE *getreference(DEFINITION *e)
{ REFERENCE *r;
  int c, m, minr, maxr = 0, mini, maxi;
  double o, mine, maxe;
  double xminr, xmaxr = 0.;

  double fracturePrefixp, fractureSuffixp, fractureRandomp;
  double fractureMinp, fractureMaxp;

  c = nextsymbol();    /* Get referenced element */
  if (! isupper(c))
    { fprintf(stderr,"\nFrag Error:\n\t*** %c is not an element name\n",c);
      exit (1);
    }
  m = c - 'A';
  if (Relement[m].type < 0)
    { fprintf(stderr,"\nFrag Error:\n\t*** %c referenced before being defined\n",c);
      exit (1);
    }
  c = nextsymbol();

  if (c == 'r')     /* Get orientation (if specified) */
    { o = 0.;
      c = nextsymbol();
    }
  else if (c == 'o')
    { if ((c = nextsymbol()) != '(') ungetc(c,sfile);
      o = getprobability("Orientation odds");
      if ((c = nextsymbol()) == ')') c = nextsymbol();
    }
  else
    o = 1.;
                    /* Get mutation rate/range (if specified) */
  if (c == 'm')
    { if ((c = nextsymbol()) != '(') ungetc(c,sfile);
      mine = maxe = getprobability("(Min) Variability");
      c = nextsymbol();
      if (c == '.' || isdigit(c) || c == ',')
        { if (c != ',') ungetc(c,sfile);
          maxe = getprobability("Max Variability");
          c = nextsymbol();
        }
      if (c == ')') c = nextsymbol();
    }
  else
    mine = maxe = 0.;
  if (mine > maxe)
    {
      fprintf(stderr,"\nFrag Error:\n\t*** Min variability (%lf) > Max variability (%lf)\n",
                     mine,maxe);
      exit (1);
    }

  if (c == 'n')
    { if ((c = nextsymbol()) != '(') ungetc(c,sfile);
      minr = getfactor("(Min) occurence",&xminr);
      if (minr < 0)
        { if ((e->type & CNCAT) != BASIS)
            { fprintf(stderr,"\nFrag Error:\n\t*** %% spec. only allowed in basis def's.\n");
              exit (1);
            }
          xmaxr = xminr;
          c = nextsymbol();
          if (isdigit(c) || c == ',' || c == '.')
            { if (c != ',') ungetc(c,sfile);
              xmaxr = getprobability("Max Occurence %");
              c = nextsymbol();
            }
          if (c == ')') c = nextsymbol();
          if (xminr > xmaxr)
            { fprintf(stderr,"\nFrag Error:\n\t*** Min occurrence %% (%g) > ",xminr);
              fprintf(stderr,"Max occurrence %% (%g)\n",xmaxr);
              exit (1);
            }
        }
      else
        { maxr = minr;
          c = nextsymbol();
          if (isdigit(c) || c == ',')
            { if (c != ',') ungetc(c,sfile);
              maxr = getnatural("Max Occurences");
              c = nextsymbol();
            }
          if (c == ')') c = nextsymbol();
          if (minr > maxr)
            { fprintf(stderr,"\nFrag Error:\n\t*** Min occurrences (%d) > ",minr);
              fprintf(stderr,"Max occurrences (%d)\n",maxr);
              exit (1);
            }
        }
    }
  else
    minr = maxr = 1;
  
  if (c == '!')
    { c = nextsymbol();
      if (c == '(' || isdigit(c))
        { if (c != '(') ungetc(c,sfile);
          if (fscanf(sfile," %d",&mini) != 1 || mini < 0)
            { fprintf(stderr,
                      "\nFrag Error:\n\t*** (Min) Instances must be a non-negative int\n");
              exit(1);
            }
          maxi = mini;
          c = nextsymbol();
          if (isdigit(c) || c == ',')
            { if (c != ',') ungetc(c,sfile);
              maxi = getnatural("Max Instances");
              c = nextsymbol();
            }
          if (c == ')') c = nextsymbol();
        }
      else
        maxi = mini = 1;
    }
  else
    mini = maxi = 0;
  if (mini > maxi)
    { fprintf(stderr,"\nFrag Error:\n\t*** Min instances (%d) > Max instances (%d)\n",
                     mini,maxi);
      exit (1);
    }

  /**  Fracture Stuff:
       f(Fract.prefix, Fract.suffix, Fract.random, Fract.minLength, Fract.maxLength)
  **/


  /* Set defaults */
  fracturePrefixp = fractureSuffixp = fractureRandomp = 0.0;
  fractureMinp = 0.0;
  fractureMaxp = 1.0;
  if (c == 'f' && (c = nextsymbol()) == '('){  /* found 'f(' */
    fracturePrefixp = getprobability("Fracture Prefix  %");
    c = nextsymbol();
    if (isdigit(c) || c == ','){
      if (c != ',') 
	ungetc(c,sfile);
      fractureSuffixp = getprobability("Fracture Suffix  %");
      c = nextsymbol();
      if (isdigit(c) || c == ','){
	if (c != ',') 
	  ungetc(c,sfile);
	fractureRandomp = getprobability("Fracture Random  %");
	c = nextsymbol();
	if (isdigit(c) || c == ','){
	  if (c != ',') 
	    ungetc(c,sfile);
	  fractureMinp = getprobability("Fracture Min  %");
	  c = nextsymbol();
	  if (isdigit(c) || c == ','){
	    if (c != ',') 
	      ungetc(c,sfile);
	    fractureMaxp = getprobability("Fracture Max  %");
	    c = nextsymbol();
	  }
	}
      }
    }
    if (c == ')') c = nextsymbol(); /* swallow the closing paren */

  }
  /* Check sanity of read results */
  if (  fracturePrefixp + fractureSuffixp + fractureRandomp > 1.0){
   fprintf(stderr,"\nFrag Error:\n\t***Total Fracture Probability out of range (%g,%g,%g) => %g\n",
	   fracturePrefixp, fractureSuffixp, fractureRandomp,
	   fracturePrefixp + fractureSuffixp + fractureRandomp   );
   exit (1);
 }





  ungetc(c,sfile);
  
/* Build reference record */

  r = (REFERENCE *) ckalloc(sizeof(REFERENCE));
  r->name   = m;
  r->orient = o;
  r->minerr = mine;
  r->maxerr = maxe;
  r->minrep = minr;
  r->maxrep = maxr;
  r->xminrp = xminr;
  r->xmaxrp = xmaxr;
  r->minins = mini;
  r->maxins = maxi;
  /*  New stuff for fracture support */
  r->fracturePrefixp = fracturePrefixp;
  r->fractureSuffixp = fractureSuffixp;
  r->fractureRandomp = fractureRandomp;
  r->fractureMinp = fractureMinp;
  r->fractureMaxp = fractureMaxp;
  /*  End New stuff for fracture support */
/* Print reference information */

  if (comments)
    { if ((e->type & CNCAT) == BASIS && e->rlist == NULL)
        printf("#    Containing:\n");
      printf("#\t %c:",r->name + 'A');
      printf(" O'odds = %.2f,",r->orient);
      printf(" Mut's = %.2f-%.2f,",r->minerr,r->maxerr);
      if (r->minrep >= 0)
        printf(" %d-%d Rep's,",r->minrep,r->maxrep);
      else
        printf(" %g-%g %% of Basis,",100*r->xminrp,100*r->xmaxrp);
      printf(" %d-%d Gen's,",r->minins,r->maxins);
      
      if(  r->fracturePrefixp + r->fractureSuffixp + r->fractureRandomp >= 1.e-09){
	printf(" Fract's = (%.2f %%p, %.2f %%s, %.2f %%r, length = %g-%g",
	       r->fracturePrefixp, r->fractureSuffixp, r->fractureRandomp,
	       r->fractureMinp, r->fractureMaxp);
      }else{
	printf(" 0%% Fractures");
      }
      printf("\n");
      fflush(stdout);
    }

  return (r);
}


/* Read and parse all element definitions from the input */

int EditMax(int size, double error)
{ int n, q;

  n = size*error;
  if (n < size*error) n += 1;
  q = n*INS_FREQ;
  if (q < n*INS_FREQ) q += 1;
  size += q - ((int) (n*DEL_FREQ));
  return (size);
}

int EditMin(int size, double error)
{ int n, q;

  n = size*error;
  q = n*DEL_FREQ;
  if (q < n*DEL_FREQ) q += 1;
  size += ((int) (n*INS_FREQ)) - q;
  return (size);
}

void getgrammar(void)
{ register int c;
  register DEFINITION *e;
  register REFERENCE  *r, *rlink;
  int name, type, mlen, llen;

  while (1)
                          /* Is there another definition? */
    { if (Sequence < 0)
        { if (comments) printf("# Seed = %d\n#\n",Seed); }
      else
        { c = nextsymbol();
          if (feof(sfile)) break;
          ungetc(c,sfile);
          if (isdigit(c)) break;
        }
      if (comments)
        { printf("# Element: ");
          fflush(stdout);
        }
                          /* Get name and basis (if any) */
      name = c = nextsymbol();
      if (! isupper(name))
        { if (c == '\n')
            fprintf(stderr,"\nFrag Error:\n\t*** Expecting a name\n");
          else
            fprintf(stdout,"\nFrag Error:\n\t*** %c is not an uppercase letter\n",name);
          exit (1);
        }
      else if (Relement[name-'A'].type >= 0)
        { fprintf(stderr,"\nFrag Error:\n\t*** %c is already defined\n",name);
          exit (1);
        }
      e = Relement + (name - 'A');
      c = nextsymbol();
      switch(c){
      case '=':
	type = BASIS;
	break;
      case '~':
	type = GLOBAL;
	break;
      case '<':
	{
	  char c;
	  type = GLOBAL | CONSTANT;
	  c = nextsymbol();
	  if(c == '"'){
	    getElementFromInlineString(e);
	  }else{
	    ungetc(c,sfile);
	    getElementFromFile(e);
	  }
#ifdef DEBUG
	  fprintf(stderr,"* Read e->seq %s\n", e->seq);
#endif
	}
	  break;
      default:
	fprintf(stderr,"\nFrag Error:\n\t*** Expecting an =, ~ or < after %c\n",name);
	exit(1);
      }
      if(!(type & CONSTANT)){
	c = nextsymbol();
	if (isdigit(c))
	  { type = BASIS | type;
          ungetc(c,sfile);
          getbasis(e);
          c = nextsymbol();
	  }
	else
	  { type = CNCAT | type;
          if (Sequence < 0)
            { fprintf(stderr,"\nFrag Error:\n\t*** 1st declaration must have a basis\n");
	    exit(1);
            }
	  }
	ungetc(c,sfile);
      }

      e->type  = type;
      e->rlist = NULL;
      Sequence = name - 'A';

      if (comments)
        { printf("%c",name);
          if (e->type & GLOBAL)
            printf(" (static)");
          printf("\n");
          if ((e->type & CNCAT) == BASIS)
            { 
	      printf("#    Length = [%d,%d]",
                     e->minlen,e->maxlen);
	      if(!(e->type & CONSTANT)){
		printf(", ACGT odds = %.2f/%.2f/%.2f/%.2f\n",
		     e->probA,e->probAC - e->probA,
                     e->probACG - e->probAC, 1.0 - e->probACG);
	      }else{
		printf(", (constant)\n");
	      }
	      
            }
          else
            printf("#    Concatenation of:\n");
          fflush(stdout);
        }

      rlink = (REFERENCE *) e;    /* Get references */
      mlen  = llen = 0;
      while (1)
        { c = nextsymbol();
          if (c == ';')
            { if ((e->type & 1) == CNCAT && e->rlist == NULL)
                { fprintf(stderr,"\nFrag Error:\n\t *** Must define at least one element\n");
                  exit (1);
                }
              break;
            }
          else
            { 
	      int maxlen,minlen;
	      ungetc(c,sfile);
              r = getreference(e);
              maxlen = Relement[r->name].maxlen;
              minlen = Relement[r->name].minlen;
              if (r->minrep >= 0)
                { if (r->maxins)
                    { mlen += EditMax(maxlen,r->maxerr)*r->maxrep*r->maxins;
                      llen += EditMin(minlen,r->minerr)*r->minrep*r->minins;
                    }
                  else
                    { mlen += EditMax(maxlen,r->maxerr)*r->maxrep;
                      llen += EditMin(minlen,r->minerr)*r->minrep;
                    }
                }    /* Probabilistic number of repetitions */
              else
                { int x, m, n, s, l;
                  double pad;

                  x = e->minlen*r->xmaxrp;
                  pad = (1.+(INS_FREQ-DEL_FREQ)*(r->minerr+r->maxerr)/2.);
                  m = maxlen * pad;
                  n = x / m;
                  if (n*m < x) n += 1;
                  s = EditMax(maxlen,r->maxerr)*n;
                  if (n*m == x)
                    l = maxlen-1;
                  else
                    { l = x / (n*pad);
                      m = l*pad;
                      if (m * n == x) l -= 1;
                    }
                  m = minlen*pad;
                  if (l >= minlen)
                    l = EditMax(l,r->maxerr)*(n+1);
                  else if (m * (n+1) <= e->maxlen*r->xmaxrp)
                    l = (EditMax(minlen,r->maxerr)*x)/m;
                  else
                    l = s;
                  if (s < l) s = l;
                  if (r->maxins)
                    mlen += s*r->maxins;
                  else
                    mlen += s;
                }
              rlink = rlink->rlist = r;
              r->rlist = NULL;
            }
        }

      if ((e->type & 1) == BASIS)
        { if (e->minlen < mlen)
            { fprintf(stderr,"\nFrag Error:\n\t*** Min basis length(%d) < ",e->minlen);
              fprintf(stderr,"Possible sum of contained elements");
              fprintf(stderr," (%d)\n",mlen);
              exit (1);
            }
        }
      else
        { e->maxlen = mlen;
          e->minlen = llen;
          mlen = 0;
          for (r = e->rlist; r != NULL; r = r->rlist)
            mlen += r->minrep;
          if (mlen == 0)
            { fprintf(stderr,"\nFrag Error:\n\t*** Potentially empty concatenation\n");
              exit (1);
            }
        }
    }
}

void getfile(void)
{ int i, n, top;
  char *buffer;

  top     = 2048;
  buffer  = (char *) ckalloc(sizeof(char)*(top+1));
  SeqFile = NULL;
  if (comments)
    { printf("# Seed = %d\n#\n",Seed);
      fflush(stdout);
    }
  n = nextsymbol();
  n = nextsymbol();
  i = 0;
  buffer[i++] = n;
  n = fgetc(sfile);
  while ( ! isspace(n))
    { if (i >= top)
        { buffer = (char *) realloc(buffer,sizeof(char)*(2*top + 1));
          if (buffer == NULL)
            { fprintf(stderr,"Frag Error: Out of Memory (frag %d)\n",2*top+1);
              exit (2);
            }
          top *= 2;
        }
      buffer[i++] = n;
      n = fgetc(sfile);
    }
  buffer[i] = '\0';
  SeqFile = fopen(buffer,"r");
  if (SeqFile == NULL)
    { fprintf(stderr,"\nFrag Error:\n\t*** Cannot open file '%s'\n",buffer);
      exit(1);
    }
  if (comments)
    { printf("# Sequence from file: %s\n",buffer);
      fflush(stdout);
    }
}

/*  Read in the input specification, honoring the old format
    if it is present.  */

int getinput(void)
{ register int n;

  for (n = 0; n < 26; n++)
    Relement[n].type = -1;
  Sequence = -1;
  
  if (comments)
    printf("#\n");

  n = nextsymbol();
  if (feof(sfile)) exit (0);
  ungetc(n,sfile);
  if (n == '<')
    getfile();
  else
    getgrammar();

  n = nextsymbol();
  if (feof(sfile))
    { if (comments) fprintf(stdout,"#\n");
      return (1);
    }
  ungetc(n,sfile);

  if (comments)
    { printf("# Fragments:\n");   /* Get fragment parameters */
      printf("#    Number = ");
      fflush(stdout);
    }
  NumFrag = getnatural("Number of fragments");
  if (comments)
    { printf("%d, Length Range = ",NumFrag);
      fflush(stdout);
    }
  MinLen = getnatural("Min fragment length");
  MaxLen = getnatural("Max fragment length");
  if (MinLen > MaxLen)
    { fprintf(stderr,"\n\t*** Min len(%d) > Max len(%d)\n",MinLen,MaxLen);
      exit(1);
    }

  if (comments)
    { printf("[%d,%d], F/R odds = ",MinLen,MaxLen);
      fflush(stdout);
    }
  Orient = getprobability("Orientation Odds");
  if (comments)
    { printf("%.2f/%.2f\n",Orient,1.0-Orient);
      printf("# Edit Characteristics:\n");
      printf("#    Error Ramp = ");
      fflush(stdout);
    }
  BegErr = getprobability("Start error rate");
  if (comments)
    { printf("%.2f->",BegErr);
      fflush(stdout);
    }
  EndErr = getprobability("Final error rate");
  if (comments)
    { printf("%.2f, Ins/Del/Sub Odds = ",EndErr);
      fflush(stdout);
    }
  if(BegErr < 1.e-9 &&
     EndErr < 1.e-9){ // Special case of NO sequencing error
    IntroduceSequenceError = FALSE;
    LevelQV = '0' + 60;
    ValleyBottom = HillPeak = 0;
  }else{
    IntroduceSequenceError = TRUE;
    LevelQV = (int) (QVfit1*log10((double) (BegErr+EndErr)/2 ) + QVfit2) + '0';
    ValleyBottom = LevelQV - ValleyDepth;
    ValleyBottom = (ValleyBottom > '0')?ValleyBottom:'0';  // min allowed QV is '0'
    HillPeak = LevelQV + HillHeight;
    HillPeak = (HillPeak < 108)?HillPeak:108;  // max allowed QV is 108 = 'l' (ell)
  }

  ProbI   = getprobability("Odds of inserting");
  ProbID  = ProbI + getprobability("Odds of deleting");
  if (ProbID > 1.)
    { fprintf(stdout,"\n\t*** Sum (%.2f) of Indel odds > 1\n",ProbID);
      exit(1);
    }
  if (comments)
    { printf("%.2f/%.2f/%.2f\n",ProbI,ProbID-ProbI,1.-ProbID);
      fflush(stdout);
    }

  n = nextsymbol();
  if (n == '\n')
    n = fgetc(sfile);
  ungetc(n,sfile);

  SingleFract = 1.0;
  if (feof(sfile))
    { if (comments) fprintf(stdout,"#\n");
      return (0);
    }

  if (comments)
    { printf("# Dual-End Inserts:\n");   /* Get fragment parameters */
      printf("#    Single Odds = "); /* Get dual-end parameters */
      fflush(stdout);
    }
  SingleFract = getprobability("Odds of single read");
  if (comments)
    { printf("%.2f\n#    Insert Range = ",SingleFract);
      fflush(stdout);
    }
  MinIns  = getnatural("Min insert length");
  MaxIns  = getnatural("Max insert length");
  if (MinIns > MaxIns)
    { fprintf(stderr,"\n\t*** Min insert len(%d) > Max insert len(%d)\n",
                     MinIns,MaxIns);
      exit(1);
    }
  else if (MinIns < MaxLen)
    { fprintf(stderr,"\n\t*** Min insert len(%d) < Max fragment len(%d)\n",
                     MinIns,MaxLen);
      exit(1);
    }
  if (comments)
    { printf("[%d,%d]\n#    Pairing Error Rate = ",MinIns,MaxIns);
      fflush(stdout);
    }
  PairError  = getprobability("Pairing Error Rate");
  if (comments)
    { printf("%.2f\n",PairError);
      printf("#\n");
    }

  return (0);
}

/* Read_seq gets the first FASTA formated sequence from the given input file.
   It is designed to handle lines and strings of arbitrary length.  It
   allocates memory for the sequence and returns a pointer to it.  

   If we encounter characters other than actgACTG we nuke 'em!!! SAK */

#define LBUFLEN 512

char *read_seq(FILE *input, int *len)
{ char  linebuf[LBUFLEN];
  register char *seqbuf;
  register int l, e, bol = 0, top;
  int charsDeleted = 0;
  int charsSubstituted = 0;

  top = 2048;
  seqbuf = (char *) ckalloc(sizeof(char)*top);
  if (fgets(linebuf,LBUFLEN,input) == NULL) return (NULL);
  if (*linebuf == '#')
  {
    do
    {
      if (bol && *linebuf == '>')
        break;
      else if (bol)
      {
        if (*linebuf != '#')
        {
          fprintf(stderr,"Frag read_seq: No FASTA header after a comment\n** %s\n",linebuf);
          exit (1);
        }
        else
          /* Put your comment processing code here */;
      }
      l = strlen(linebuf);
      bol = (linebuf[l-1] == '\n');
    }
    while (fgets(linebuf,LBUFLEN,input) != NULL);
  }

  if (*linebuf != '>')
    { fprintf(stderr,"First line must start with an >-sign\n");
      exit (1);
    }

  do
    { l = strlen(linebuf);
      if (linebuf[l-1] == '\n') break;
    }
  while (fgets(linebuf,LBUFLEN,input) != NULL);

  bol = 1;
  e   = 0;
  while (fgets(linebuf,LBUFLEN,input) != NULL)
    { if (bol && *linebuf == '>')
        break;
      l = strlen(linebuf);

      if (e + l >= top)
        { top = 1.5*(e+l) + 200;
#ifdef USE_CKALLOC
          newbuf = (char *) ckalloc(sizeof(char)*top);
          seqbuf[e] = '\0';
          strcpy(newbuf,seqbuf);
          free(seqbuf);
          seqbuf = newbuf;
#else
	  /* This should be marginally more efficient, since strcpy is relatively costly */
	  seqbuf = (char *) realloc(seqbuf, sizeof(char)*top);
	  assert(NULL != seqbuf);
#endif
        }
      
      // Strip out bogus characters that will confuse frag  SAK
      {
	int s,t, d,sub;
	for(t = 0,s = 0,d = 0,sub=0; s < l; s++){
	  switch(tolower(linebuf[s])){
	  case 'a':
	  case 'c':
	  case 't':
	  case 'g':
	  case '\n':   // these get deleted below
	    seqbuf[e + t++] = linebuf[s];
	    break;
	  case 'm':
	  case 'r':
	  case 'w':
	  case 's':
	  case 'y':
	  case 'k':
	  case 'v':
	  case 'h':
	  case 'd':
	  case 'b':
	    if(protoIO == 0)
	      sub++;
	    /* Intentional fallthrough */
	  case 'n':  // These are tolerated
	    if(protoIO == 0){
	      seqbuf[e + t++] = 'n';  // We don't know how to complement these yet
	    }else{
	      d++;
	    }
	    break;
	  default:
	    d++;
	    break;
	  }
	}
	charsDeleted += d;
	charsSubstituted += sub;
      bol = (linebuf[l-1] == '\n');
      // l - e is the number of non-bogus characters read
      e = (e+l-d) - bol;
      }
    }
  seqbuf[e] = '\0';
  if(charsDeleted > 0 || charsSubstituted > 0)
    fprintf(stderr,"* Read_seq deleted %d bad characters and substituted %d\n",
	    charsDeleted, charsSubstituted);

  {
    int i;
    for(i = 0; i < e; i++)
      {
      if(isspace(seqbuf[i]))
	 fprintf(stderr,"* Non Alpha character 0x%x found at pos %d of seq\n",
		 seqbuf[i],i);
      }
  }
  if (e == 0)
    { fprintf(stderr,"Input sequence is empty\n");
      exit (1);
    }

  *len = e;

  return (seqbuf);
}


/* >>>> MAIN <<<< */

int main(int argc, char *argv[])
{ char *dna, *pname;
  int   idx, len, i;
  int   no_seed, dna_only, illegal;

/* Initialize random number generator with process id. */

  pname = argv[0];

  no_seed  = 1;
  no_elem  = 1;
  comments = 1;
  printVersion = 0;
  notation = 0;
  uniform = 0;
  illegal  = 0;
  protoIO = 0;
  while (argc > 1 && argv[1][0] == '-')
    { if (argv[1][1] == 's')
        { if (argc > 2)
            { Seed    = atoi(argv[2]);
              no_seed = 0;
              argc   -= 1;
              argv   += 1;
            }
        }
      else
        { for (i = 1; argv[1][i] != '\0'; i++)
            if (argv[1][i] == 'p')
              no_elem = 0;
	    else if (argv[1][i] == 'P')
	      protoIO = 1;
            else if (argv[1][i] == 'F')
              comments = 0;
            else if (argv[1][i] == 'N')
              notation = 1;
            else if (argv[1][i] == 'V')
              printVersion = 1;
            else if (argv[1][i] == 'u')
              uniform = 1;
            else if (argv[1][i] == 'n')
              uniform = 0;
            else
              illegal = 1;
          if (i == 1) illegal = 1; 
        }
      argv += 1;
      argc -= 1;
    }

  if (argc != 2 || illegal)
    { fprintf(stderr,"Error: Usage: %s [-s {seed:%%d}] [-npuFNV] specfile\n",pname);
      exit (1);
    }

    fprintf(stderr,"* Using %s distribution seed = %d\n",
	    (uniform?"UNIFORM":"NORMAL"), Seed);

  if ((sfile = fopen(argv[1],"r")) == NULL)
    { fprintf(stderr,"Error: Cannot open spec file %s\n",argv[1]);
      exit (1);
    }

  fprintf(stderr,"* spec file is %s\n", argv[1]);
  { /* open file for quality output */
    char *suffix;
    char *psuffix;
    suffix = strrchr(argv[1],(int)'.');
    psuffix = strchr(argv[1],(int)'.');
    fprintf(stderr,"Input file is %s suffix is %s psuffix is %s\n",argv[1], suffix, psuffix);
    *psuffix = '\0';
    sprintf(Quality_File_Name,"%s.qlt%s",argv[1],suffix);
    fprintf(stderr,"Quality file is %s \n",Quality_File_Name);
    if ((QualityFile = fopen(Quality_File_Name, "w")) == NULL ) { 
      { fprintf(stderr,"Error: Cannot open quality file %s\n",Quality_File_Name);
        exit (1);
      }
    }
  }


  if (comments|| printVersion){
    printf("# Frag %s:(%s)\n", FragVersion, (uniform?"uniform":"gaussian"));
    printf("#\n");
  }
  if (no_seed)
    Seed = getpid();
  srand48(Seed);

/* Setup mutation tables */

  TranChar['a'] = 'a';   FlipChar['a'] = 't';   Decode['a'] = 0;
  TranChar['c'] = 'c';   FlipChar['c'] = 'g';   Decode['c'] = 1;
  TranChar['g'] = 'g';   FlipChar['g'] = 'c';   Decode['g'] = 2;
  TranChar['t'] = 't';   FlipChar['t'] = 'a';   Decode['t'] = 3;
  TranChar['A'] = 'A';   FlipChar['A'] = 'T';   Decode['A'] = 4;
  TranChar['C'] = 'C';   FlipChar['C'] = 'G';   Decode['C'] = 5;
  TranChar['G'] = 'G';   FlipChar['G'] = 'C';   Decode['G'] = 6;
  TranChar['T'] = 'T';   FlipChar['T'] = 'A';   Decode['T'] = 7;
  TranChar['N'] = 'N';   FlipChar['N'] = 'N';   Decode['N'] = 8;
  TranChar['n'] = 'n';   FlipChar['n'] = 'n';   Decode['n'] = 9;
 
/* Read in the specification from sfile. */

  dna_only = getinput();

/* Generate a DNA sequence.  */

  if (Sequence < 0)
    dna = read_seq(SeqFile,&len);
  else
    { build_seq();
      dna = Relement[Sequence].valstack->value;
      len = Relement[Sequence].valstack->length;
    }

/* Sample NumFrag samples from this sequence,
   editing each fragment after sampling,
   and generating appropriate output .
*/

#ifdef DEBUGS
  fprintf(stderr,"dna = '%s'\n",dna);
#endif

  if (dna_only)
    { register int i;
      if (comments)
        { printf("# DNA Sequence:\n");
          printf("#\n");
        }
      printf("> 0\n");
      for (i = 0; i < len; i += SEGLEN)
        if (len-i < SEGLEN)
          printf("%.*s\n",len-i,dna+i);
        else
          printf("%.*s\n",SEGLEN,dna+i);
    }

  else
    { if (comments) printf("# Generated Fragments:\n");
    fprintf(stderr,"* Generating %d fragments\n", NumFrag);
      if (SingleFract >= 1.)
        for (idx = 0; idx < NumFrag; idx++)
          { newfrag(len);
            editfrag(dna);
	    //            printf("> %d\n",idx+1);
            outputfrag(idx+1,' ',FALSE);
          }
      else
        { SingleFract = 2*SingleFract/(1 + SingleFract);
          idx = 1;
          while (idx <= NumFrag)
            if (idx == NumFrag || drand48() < SingleFract)
              { newfrag(len);
                editfrag(dna);
		//                printf("> %d\n",idx);
                outputfrag(idx,' ',FALSE);
                idx += 1;
              }
            else if (drand48() < PairError)
              { newfrag(len);
                editfrag(dna);
		//                if (comments || notation){
		//printf("# Not Related to %dr!\n#\n",idx);
                  // fprintf(stderr,"# Not Related to %dr! PairError = %g\n#\n",idx,PairError);
		//		}
		//               printf("> %df\n",idx);
                outputfrag(idx,'f',TRUE);
                newfrag(len);
                editfrag(dna);
		//                if (comments || notation){
		//                  printf("# Not Related to %df!\n#\n",idx);
                  // fprintf(stderr,"# Not Related to %df\n",idx);
		//		}
		//                printf("> %dr\n",idx);
                outputfrag(idx,'r',TRUE);
                idx += 2;
              }
            else
              { firstfrag(len);
                editfrag(dna);
		//                printf("> %df\n",idx);
                outputfrag(idx,'f',FALSE);
                secondfrag(len);
                editfrag(dna);
		//                printf("> %dr\n",idx);
                outputfrag(idx,'r',FALSE);
                idx += 2;
              }
        }
    }
  /* This is so insure will not detect this as a leak */
  if(Sequence < 0)
    free(dna);
  exit(0);
}
