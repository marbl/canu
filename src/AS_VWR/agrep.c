
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
/* $Id: agrep.c,v 1.4 2005-03-22 19:49:30 jason_miller Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#undef DEBUG

typedef void *regexp;

#define MALLOC(x) malloc(x)
#define FREE(x)   free(x)

#define POUND_AT_FEATURES
#define NUMBER_FEATURE

/* REGULAR EXPRESSION DATA TYPE:
     An e-nfa and work storage to do approximate matching.  Occupies
     at most 60N bytes for a regular expression of length N.
*/

#define BEGIN    0
#define FINISH   1
#define EPSILON  2
#define SYMBOL   3
#define CLASS    4

struct state { int    kind;
               int    data;
               int    mark;
               int    prd1, prd2;
               int    suc1, suc2;
             };

struct a_mach { struct state *start;
                int          *loop;
                char         *cvecs;
                int           size, back;
                int          *work1, *work2;
              };

static struct node *or(void);


/* RE_PARSE:
     Normally returns regexp object used to search for instances of a
     regular expression.  The parameter i_opt is a boolean input
     parameter and should be set if the search is to be case-insensitive.
     Returns NULL and sets most recent error if the regular expression
     is not syntactically valid or if there is insufficient memory.
*/

   /* Parse tree structure */

#define OR      0
#define CAT     1
#define PLUS    2
#define STAR    3
#define OPT     4
#define ATOM_S  5
#define ATOM_C  6

struct node  { int         label;
               struct node *left, *right;
             };

struct leaf  { int    label;
               int    vecoff; };

   /* Space buffer for parse tree: 44*parse_top bytes */

static int parse_top = -1;

static int itop, etop, ctop;

static struct node *interior;
static struct leaf *exterior;
static char        *classes;

   /* Global communication between re_parse & parser routines */

static int  size;     /* Compute size of e-nfa in this var */
static int  i_flag;   /* Case insensitive? */
static char *scan;    /* Current scan pointer */
static char *patbeg;  /* Pointer to first char of pattern */

static int  only_digits;   /* Only digits used in pattern? */
static int  digit_check;   /* Check that only digits are used in pattern */

   /* Synerr sets most recent error and location */

#define RP_ERROR "Missing )"
#define UC_ERROR "Missing ]"
#define NC_ERROR "\\n in class"
#define PT_ERROR "Unexpected end"
#define BE_ERROR "\\ is last char"
#define OM_ERROR out_of_memory
#define NR_ERROR "Range syntax error"
#define NV_ERROR "1st # not <= 2nd"
#define ES_ERROR "Illegal char, try $"

static char **message;
static int   *position;

static char *out_of_memory = "Out of memory";

static void synerr(char *msg)
{ static char bufmess[150];

  sprintf(bufmess,"%s",msg);
  *message  = bufmess;
  *position = scan - patbeg;
}

   /* Parse tree node allocation routines */

static struct node *salloc(int t, int v)
{ register struct leaf *l;

  l = exterior + etop++;
  l->label  = t;
  l->vecoff = v;
  return ((struct node *) l);
}

static struct node *nalloc(int b, struct node *l, struct node *r)
{ register struct node *p;

  if (b <= CAT && r == NULL) return (NULL);
  p = interior + itop++;
  p->label = b;
  p->left  = l;
  p->right = r;
  return (p);
}

   /* Character class construction and allocation routines */

static char vec[16];

static void cclstart(void)
{ register int c;
  for (c = 0; c <= 15; c++)
    vec[c] = 0;
}

static void cclin(register int c)
{ vec[c>>3] |= (1 << c%8);
  if (i_flag)
    { if ( isupper(c))
        { c += 040;
          vec[c>>3] |= (1 << c%8);
        }
      else if (islower(c))
        { c -= 040;
          vec[c>>3] |= (1 << c%8);
        }
    }
}

static int cclcomp(int comp)
{ register int   c;
  register char *res;

  res = classes + ctop;
  ctop += 16;
  if (comp)
    for (c = 0; c <= 15; c++)
      res[c] = (vec[c] ^ 0377);
  else
    for (c = 0; c <= 15; c++)
      res[c] = vec[c];
  res[0] &= 0xFC;

  if (digit_check)
    { for (c = 0; c <= 5; c++)
        if (res[c] != 0) only_digits = 0;
      for (c = 8; c <= 15; c++)
        if (res[c] != 0) only_digits = 0;
      for (c = 2; c <= 7; c++)
        if (res[7] & (1<<c%8)) only_digits = 0;
    }
 
  return (ctop-16);
}

   /* Routines to add special subtrees */

#define NUM_SPECIAL 3
#define SIZE_DIGIT  7
#define SIZE_IDENT  4

static struct node *dig_tree;
static struct node *id_tree;

static struct node *digit(void)
{ register struct node *p;
  register int x;

  if (dig_tree == NULL)
    { cclstart();
      for (x = '0'; x <= '9'; x++)
        cclin(x);
      p = nalloc(STAR,salloc(ATOM_C,cclcomp(0)),NULL);
 
      cclstart();
      for (x = '1'; x <= '9'; x++)
        cclin(x);
      p = nalloc(CAT,salloc(ATOM_C,cclcomp(0)),p);

      dig_tree = nalloc(OR,salloc(ATOM_S,(int) '0'),p);
    }
  return (dig_tree);
}

static struct node *identifier(void)
{ register struct node *p;
  register int x;

  if (id_tree == NULL)
    { cclstart();
      for (x = 'a'; x <= 'z'; x++)
        cclin(x);
      for (x = 'A'; x <= 'Z'; x++)
        cclin(x);
      for (x = '0'; x <= '9'; x++)
        cclin(x);
      cclin('_');
      p = nalloc(STAR,salloc(ATOM_C,cclcomp(0)),NULL);
 
      cclstart();
      for (x = 'a'; x <= 'z'; x++)
        cclin(x);
      for (x = 'A'; x <= 'Z'; x++)
        cclin(x);
      cclin('_');
      id_tree = nalloc(CAT,salloc(ATOM_C,cclcomp(0)),p);
    }
  return (id_tree);
}

  /* Routines to convert number ranges to regular expressions */

static char *g_num1, *g_num2;
static int   g_off,   g_max;

struct node *zero_nine;

static struct node *make_09(void)
{ register int x;

  cclstart();
  for (x = '0'; x <= '9'; x++)
    cclin(x);
  return (salloc(ATOM_C,cclcomp(0)));
}

static struct node *genfil(int lev)
{ register int i;
  register struct node *p;

  size += (g_max-lev);
  if (zero_nine == NULL)
    zero_nine = make_09();
  p = zero_nine;
  for (i = lev+1; i < g_max; i++)
    p = nalloc(CAT,p,zero_nine);
  return (p);
}

static struct node *genabv(int lev)
{ register int c, o, x;
  register struct node *p, *q;

  o = lev - g_off;
  if (o < 0)
    c = '0';
  else
    c = g_num1[o];

  if (lev == g_max-1)
    { size += 1;
      if (c < '9')
        { cclstart();
          for (x = c; x <= '9'; x++)
            cclin(x);
          return (salloc(ATOM_C,cclcomp(0)));
        }
      else
        return (salloc(ATOM_S,c));
    }
  else
    { p = genabv(lev+1);
      if (o >= 0)
        { size += 1;
          p = nalloc(CAT,salloc(ATOM_S,c),p);
        }
      if (c < '9')
        { size += 3;
          if (c < '8')
            { cclstart();
              for (x = c+1; x <= '9'; x++)
                cclin(x);
              q = salloc(ATOM_C,cclcomp(0));
            }
          else
            q = salloc(ATOM_S,'9');
          p = nalloc(OR,p,nalloc(CAT,q,genfil(lev+1)));
        }
      return (p);
    }
}

static struct node *genbel(int lev)
{ register int c, x;
  register struct node *p, *q;

  c = g_num2[lev];

  if (lev == g_max-1)
    { size += 1;
      if (c > '0')
        { cclstart();
          for (x = '0'; x <= c; x++)
            cclin(x);
          return (salloc(ATOM_C,cclcomp(0)));
        }
      else
        return (salloc(ATOM_S,c));
    }
  else
    { size += 1;
      p = nalloc(CAT,salloc(ATOM_S,c),genbel(lev+1));
      if (c > '0')
        { size += 3;
          if (c > '1')
            { cclstart();
              for (x = '0'; x <= c-1; x++)
                cclin(x);
              q = salloc(ATOM_C,cclcomp(0));
            }
          else
            q = salloc(ATOM_S,'0');
          p = nalloc(OR,p,nalloc(CAT,q,genfil(lev+1)));
        }
      return (p);
    }
}
  
static struct node *genrng(int lev)
{ register int c1, c2, o, x;
  register struct node *p, *q;

  o = lev - g_off;
  if (o < 0)
    c1 = '0';
  else
    c1 = g_num1[o];

  c2 = g_num2[lev];
  
  if (c1 == c2)
    { size += 1;
      if (lev+1 >= g_max)
        return(salloc(ATOM_S,c1));
      else
        return (nalloc(CAT,salloc(ATOM_S,c1),genrng(lev+1)));
    }
  else if (lev == g_max-1)
    { size += 1;
      cclstart();
      for (x = c1; x <= c2; x++)
        cclin(x);
      return (salloc(ATOM_C,cclcomp(0)));
    }
  else
    { p = genabv(lev+1);
      if (o >= 0)
        { size += 1;
          p = nalloc(CAT,salloc(ATOM_S,c1),p);
        }
      if (c1 < c2-1)
        { size += 3;
          if (c1 < c2-2)
            { cclstart();
              for (x = c1+1; x <= c2-1; x++)
                cclin(x);
              q = salloc(ATOM_C,cclcomp(0));
            }
          else
            q = salloc(ATOM_S,c1+1);
          p = nalloc(OR,p,nalloc(CAT,q,genfil(lev+1)));
        }
      size += 3;
      p = nalloc(OR,p,nalloc(CAT,salloc(ATOM_S,c2),genbel(lev+1)));
      return (p);
    }
}
          
static struct node *gen_number_pat(char *num1, char *num2)
{ g_num1 = num1;
  g_num2 = num2;
  g_off  = strlen(num1);
  g_max  = strlen(num2);
  g_off  = g_max - g_off;

  return (genrng(0));
}

   /* Recursive decent parse routines for regular expressions */

static struct node *fact(void)
{ struct node *p;
  register char *num1, *num2;
  register int comp, c;

  switch (*scan)

  { case '(':
      scan += 1;
      p = or();
      if (p == NULL) return (NULL);
      if (*scan != ')') 
        { synerr(RP_ERROR);
          return (NULL);
        }
      scan += 1;
      return (p);

    case '.':
      scan += 1;
      size += 1;
      cclstart();
      return (salloc(ATOM_C,cclcomp(1)));

    case '[':
      cclstart();
      comp = 0;
      if (*++scan == '^')
        { scan += 1;
          comp = 1;
        }
      while (*scan != ']')
        { if (*scan == '\\') scan += 1;
          if (*scan == '\0')
            { synerr(UC_ERROR);
              return (NULL);
            }
          if (*scan == '\n')
            { synerr(NC_ERROR);
              return (NULL);
            }
          if (scan[1] == '-' && ((isdigit(*scan) && isdigit(scan[2])) ||
                                 (isupper(*scan) && isupper(scan[2])) ||
                                 (islower(*scan) && islower(scan[2]))))
            { for (c = *scan; c <= scan[2]; c++)
                cclin(c);
              scan += 3;
            }
          else
            cclin(*scan++);
        }
      scan += 1;
      size += 1;
      return (salloc(ATOM_C,cclcomp(comp)));

#ifdef NUMBER_FEATURE

    case '{':
      num1 = ++scan;
      while (isdigit(*scan))
        scan += 1;
      if (scan <= num1)
        { synerr(NR_ERROR);
          return (NULL);
        }
      num2 = ++scan;
      while (isdigit(*scan))
        scan += 1;
      if (*scan != '}' || scan <= num2)
        { synerr(NR_ERROR);
          return (NULL);
        }
      c = num2[-1];
      num2[-1] = '\0';
      *scan    = '\0';
      if (atoi(num1) > atoi(num2))
        { num2[-1] = c;
          *scan    = '}';
          synerr(NV_ERROR);
          return (NULL);
        }
      p = gen_number_pat(num1,num2);
      num2[-1] = c;
      *scan++  = '}';
      return (p);
      
#endif

#ifdef POUND_AT_FEATURES

    case '#':
      size += SIZE_DIGIT;
      scan += 1;
      return (digit());

    case '@':
      size += SIZE_IDENT;
      scan += 1;
      only_digits = 0;
      return (identifier());

#endif

    case '\0':
      synerr(PT_ERROR);
      return (NULL);

    case '\\':
      if (*++scan == '\0')
        { synerr(BE_ERROR);
          return (NULL);
        }
      if (*scan == 't')
        { scan += 1;
          size += 1;
          only_digits = 0;
          return (salloc(ATOM_S,(int) '\t'));
        }
     if (*scan == 'n' || *scan == 'r')
       { synerr(ES_ERROR);
         return (NULL);
       }

    default:
      if (*scan == '\n')
        { synerr(ES_ERROR);
          return (NULL);
        }
      size += 1;
      if (i_flag && isalpha(*scan))
        { cclstart();
          cclin(*scan++);
          return (salloc(ATOM_C,cclcomp(0)));
        }
      else
        { if (! isdigit(*scan)) only_digits = 0;
          return (salloc(ATOM_S,(int) (*scan++)));
        }
  }
}

static struct node *rep(void)
{ struct node *p;
  p = fact();
  while (p != NULL)
    { if (*scan == '?')
        p = nalloc(OPT,p,NULL);
      else if (*scan == '*')
        p = nalloc(STAR,p,NULL);
      else if (*scan == '+')
        p = nalloc(PLUS,p,NULL);
      else
        break;
      scan += 1;
      size += 2;
    }
  return (p);
}

static struct node *con(void)
{ struct node *p;
  p = rep();
  while (p != NULL && *scan != ')'  && *scan != '|' &&
                      *scan != '\n' && *scan != '\0')
    p = nalloc(CAT,p,rep());
  return (p);
}

static struct node *or(void)
{ struct node *p;
  p = con();
  while (p != NULL && (*scan == '|' || *scan == '\n'))
    { scan += 1;
      size += 2;
      p = nalloc(OR,p,con());
    }
  return (p);
}

   /* Build_enfa traverses parse tree to construct an e-nfa for
      recognizing the regular expression.  */

static struct state *NFA;    /* Array of machine states */
static int           zip;    /* Integer code for undefined */

static int          *loop;   /* Array of machines outermost loops */
static int           outer;  /* Am I outside a loop? */
static int           back;   /* Number of outermost loops */

static int build_enfa(register struct node *p,
                      register int cnt, int hat, int dol)
{ register int fs1, fs2;

  switch (p->label)
  { case OR:
      fs1 = build_enfa(p->left ,cnt+1,hat,dol);
      fs2 = build_enfa(p->right,fs1+1,hat,dol);
      NFA[cnt  ].kind = EPSILON;
      NFA[cnt  ].mark = 0;
      NFA[cnt  ].suc1 = cnt+1;
      NFA[cnt  ].suc2 = fs1+1;
      NFA[cnt+1].prd1 = cnt;
      NFA[cnt+1].prd2 = zip;
      NFA[fs1+1].prd1 = cnt;
      NFA[fs1+1].prd2 = zip;
      NFA[fs2+1].kind = EPSILON;
      NFA[fs2+1].mark = 0;
      NFA[fs1  ].suc1 = fs2+1;
      NFA[fs1  ].suc2 = zip;
      NFA[fs2  ].suc1 = fs2+1;
      NFA[fs2  ].suc2 = zip;
      NFA[fs2+1].prd1 = fs1;
      NFA[fs2+1].prd2 = fs2;

      return (fs2+1);
    case CAT:
      fs1 = build_enfa(p->left ,cnt,hat,0);
      fs2 = build_enfa(p->right,fs1+1,0,dol);
      NFA[fs1  ].suc1 = fs1+1;
      NFA[fs1  ].suc2 = zip;
      NFA[fs1+1].prd1 = fs1;
      NFA[fs1+1].prd2 = zip;

      return (fs2);
    case PLUS:
    case STAR:
    case OPT:
      if (outer && p->label != OPT)
        { outer = 0;
          fs1 = build_enfa(p->left,cnt+1,0,0);
          outer = 1;
          loop[++back] = cnt+1;
        }
      else
        fs1 = build_enfa(p->left,cnt+1,hat,dol);
      NFA[cnt  ].kind = EPSILON;
      NFA[cnt  ].mark = 0;
      NFA[cnt  ].suc1 = cnt+1;
      NFA[cnt  ].suc2 = (p->label == PLUS ? zip : fs1+1);
      NFA[cnt+1].prd1 = cnt;
      NFA[cnt+1].prd2 = (p->label == OPT ? zip : fs1);
      NFA[fs1+1].kind = EPSILON;
      NFA[fs1+1].mark = 0;
      NFA[fs1  ].suc1 = fs1+1;
      NFA[fs1  ].suc2 = (p->label == OPT ? zip : cnt+1);
      NFA[fs1+1].prd1 = fs1;
      NFA[fs1+1].prd2 = (p->label == PLUS ? zip : cnt);

      return (fs1+1);
    case ATOM_S:
      if (hat && ((struct leaf *) p)->vecoff == '^')
        { NFA[cnt].kind = BEGIN;
          NFA[cnt].mark = 0;
          NFA[cnt].data = 1;
        }
      else if (dol && ((struct leaf *) p)->vecoff == '$')
        { NFA[cnt].kind = FINISH;
          NFA[cnt].mark = 0;
          NFA[cnt].data = 0;
        }
      else
        { NFA[cnt].kind = SYMBOL;
          NFA[cnt].mark = 0;
          NFA[cnt].data = ((struct leaf *) p)->vecoff;
        }

      return (cnt);
    default /* case ATOM_C */:
      NFA[cnt].kind = CLASS;
      NFA[cnt].mark = 0;
      NFA[cnt].data = ((struct leaf *) p)->vecoff;
      return (cnt);
    }
}

  /* Re_parse top level. */

regexp re_parse(int i_opt, char *pat, char **mesg, int *pos)
{ struct node *p;
  register char *s;
  register int   x, k;
  int patlen, even_size;
  int patplus;
  struct a_mach *mach;
#ifdef DEBUG
  static void show_tree(), show_nfa();
#endif

  patlen  = strlen(pat);

  message  = mesg;
  position = pos;

  patplus = patlen;

#ifdef NUMBER_FEATURE
  x = k = 0;
  for (s = pat; *s != '\0'; s++)
    switch (x)
    { case 0:
        if (*s == '\\')
          x = 1;
        else if (*s == '{')
          x = 2;
        break;
      case 1:
        x = 0;
        break;
      case 2:
        k = 0;
        if (!isdigit(*s)) x = 3;
        break;
      case 3:
        if (*s == '}')
          { if (k > 0) patplus += (k+2)*k-2;
            x = 0;
          }
        if (!isdigit(*s)) x = 0;
        k += 1;
        break;
    }
  zero_nine = NULL;
#endif

#ifdef POUND_AT_FEATURES
  patplus += NUM_SPECIAL;
  dig_tree = id_tree = NULL;
#endif

  if (patplus > parse_top)  /* Prepare parse tree memory buffers */
    { if (parse_top >= 0)
        FREE(interior);
      parse_top = patplus*1.2 + 30;
      s = (char *) MALLOC(parse_top*(sizeof(struct leaf) +
                                     sizeof(struct node) + 
                                     16*sizeof(char)));
      if (s == NULL)
        { synerr(OM_ERROR);
          return (NULL);
        }
      interior = (struct node *) s;
      exterior = (struct leaf *) (s + parse_top*sizeof(struct node));
      classes  = (char *) (s + parse_top*(sizeof(struct node) +
                                          sizeof(struct leaf)));
    }
  etop = 0;
  itop = 0;
  ctop = 0;

  i_flag = i_opt;   /* Parse regular expression */
  size   = 0;
  patbeg = pat;
  scan   = pat;

  p = or();

  if (p != NULL && *scan != '\0')
    { synerr(PT_ERROR);
      p = NULL;
    }
  if (p == NULL) return (NULL);

#ifdef DEBUG
  show_tree(p);
#endif

    /* Allocate e-nfa */

  even_size = size+1;
  if (even_size % 2) even_size += 1;

  s = (char *) MALLOC(sizeof(struct a_mach) + ctop*sizeof(char) +
                      (3*sizeof(int)+sizeof(struct state))*even_size);
  if (s == NULL)
    { synerr(OM_ERROR);
      return (NULL);
    }

  mach = (struct a_mach *) (s + ctop*sizeof(char) +
                            (3*sizeof(int)+sizeof(struct state))*even_size);
  mach->cvecs = (char *) (s + (3*sizeof(int)+sizeof(struct state))*even_size);
  mach->start = (struct state *) (s + 3*sizeof(int)*even_size);
  mach->loop  = (int *) (s + 2*sizeof(int)*even_size);
  mach->work1 = (int *) (s + sizeof(int)*even_size);
  mach->work2 = (int *) (s);

  mach->size = size;     /* Build machine */

  for (x = 0; x < ctop; x++)
    mach->cvecs[x] = classes[x];

  NFA = mach->start;
  zip = size+1;

  loop  = mach->loop;
  back  = 0;
  outer = 1;

  NFA[0].kind = EPSILON;
  NFA[0].mark = 0;
  NFA[0].suc1 = 1;
  NFA[0].suc2 = zip;
  NFA[1].prd1 = 0;
  NFA[1].prd2 = zip;

  x = build_enfa(p,1,1,1);
  mach->back = back;

  NFA[0].prd1 = NFA[0].prd2 = zip;
  NFA[x].suc1 = NFA[x].suc2 = zip;

#ifdef DEBUG
  show_nfa(mach);
#endif

  return ((regexp) mach);
}

/* RE_MATCH:
     Return 0 or 1 according to whether string is recognized by machine
     regexp with at most thresh differences.  For algorithm details
     see Myers & Miller (1989) Bulletin of Mathematical Biology.
*/

static int approx_match(regexp rexpr, char *string, int thresh)
{ int a, j;
  register int x, t, *C, *B;
  register struct state *s;
  register int p;
  register char *cv;
  struct state  *NFA;
  struct a_mach *mach;
  int   *loop, size, back;
  int    b, k;
  int   *E;

  mach  = (struct a_mach *) rexpr;

  NFA   = mach->start;
  loop  = mach->loop;
  size  = mach->size;
  back  = mach->back;
  cv    = mach->cvecs;

  C = mach->work1;
  B = mach->work2;

  C[0] = 0.;
  for (j = 1; j <= size; j++)
    { s = NFA+j;
      t = (s->kind > EPSILON);
      C[j] = C[s->prd1] + t;
      if (s->prd2 < j)
        { if ((t = C[s->prd2] + t) < C[j]) C[j] = t; }
      if (s->kind == FINISH && *string != '\0')
        C[j] = thresh+1;
    }
#ifdef DEBUG
  printf("^:");
  for (j = 0; j <= size; j++)
    printf(" %1d",C[j]);
  if (C[size] <= thresh)
    printf(" Match\n");
  printf("\n");
#else
  if (C[size] <= thresh)
    return (1);
#endif

  while ((a = *string++) != '\0')
    { E = C; C = B; B = E;
      C[0] = 0.;
      for (j = 1; j <= size; j++)
        { s = NFA+j;
          p = s->kind;
          x = C[s->prd1];
          if (s->prd2 < j)
            { if ((t = C[s->prd2]) < x) x = t; }
          if (p > EPSILON)
            { x += 1;
              if ((t = B[j] + 1) < x) x = t;
              if (p == SYMBOL)
                p = (s->data != a);
              else
                p = ((cv[s->data + (a>>3)] & (1 << a%8)) == 0);
              if ((t = B[s->prd1] + p) < x) x = t;
              if (s->prd2 < size)
                { if ((t = B[s->prd2] + p) < x) x = t; }
            }
          else if (p < EPSILON)
            { if (p == BEGIN)
                x = B[j]+1;
              else if (*string != '\0')
                x = thresh+1;
            } 
          C[j] = x;
        }
      for (b = 1; b <= back; b++)
        { j = loop[b];
          k = NFA[j].prd2;
          if ((t = C[k] + (NFA[j].kind != 0)) < C[j])
            { C[j] = t;
              for (j = j+1; j < k; j++)
                { s = NFA+j;
                  x = C[s->prd1];
                  if (s->prd2 < size)
                    { if ((t = C[s->prd2]) < x) x = t; }
                  x += (s->kind != EPSILON);
                  if (x < C[j]) C[j] = x;
                }
            }
        }
#ifdef DEBUG
      printf("%c:",a);
      for (j = 0; j <= size; j++)
        printf(" %1d",C[j]);
      if (C[size] <= thresh)
        printf(" Match\n");
      printf("\n");
#else
      if (C[size] <= thresh)
        return (1);
#endif
    }
#ifdef DEBUG
  printf("Score = %d\n",C[size]);
#endif
  return (0);
}



/* DFA-based exact r.e. matching */

static struct state *aNFA;
static int           azip;
static char         *aclass;

#define SIGMA 128
#define AVESET 20
#define HSIZE 179

struct chain { struct chain *hlink;
               int          *set;
             };

static struct chain *bucket[HSIZE];

static struct chain  hash[SIGMA];

static char dfa[SIGMA][SIGMA];
static int  final[SIGMA];
static int  dalloc;

static int  setvec[AVESET*SIGMA], *settop;
static int *setlimit = setvec + AVESET*SIGMA;

static int *sbot, *dbot;

static int *goset(int s, int c, int *hvalue)
{ register int *src, *dst;
  register int  stp,  gop;
  register unsigned h;

  src = sbot-1;
  for (dst = hash[s].set; *dst != azip; dst++)
    *++src = *dst;

#define XFER  aNFA[*++dst = gop].mark = 1

#define GOTO							\
if (!aNFA[gop].mark)						\
  if (aNFA[gop].kind == EPSILON)				\
    { *++src = gop; XFER; }					\
  else if (aNFA[gop].kind != CLASS)				\
    { if (aNFA[gop].data == c) XFER; }				\
  else								\
    { if (aclass[aNFA[gop].data + (c>>3)] & 1<<c%8) XFER; }

  final[s] = 0;

  dst = dbot-1;
  while (src >= sbot)
    { stp = *src--;
      gop = aNFA[stp].suc1;
      if (gop == azip)
        { final[s] = 1;
          continue;
         }
      GOTO
      gop = aNFA[stp].suc2;
      if (gop == azip) continue;
      GOTO
    }

  h = 0;
  while (dst >= dbot)
    { gop = *dst--;
      if (aNFA[gop].kind != EPSILON)
        h ^= (unsigned) (*++src = gop);
      else
        aNFA[gop].mark = 0;
    }
  h ^= (unsigned) (*++src = 0);
  aNFA[0].mark = 1;
  *hvalue = h%HSIZE;
  return (src);
}

static int hsearch(int hvalue, int size)
{ register struct chain *alter;
  register int *base, *btop;

  for (alter = bucket[hvalue]; alter != NULL; alter = alter->hlink)
    { base = alter->set;
      for (btop = base; *btop != azip; btop++)
        if (!aNFA[*btop].mark) break;
      if (*btop == azip && size == (btop-base)) break;
    }
  if (alter == NULL)
    return (-1);
  else
    return (alter-hash);
}

static void clear(void)
{ register int i;
  for (i = 0; i < HSIZE; i++)
    bucket[i] = NULL;
  for (i = 0; i < SIGMA; i++)
    dfa[0][i] = 0;
  final[0] = -1;
  dalloc = 1;
  settop = setvec+2;
#ifdef DEBUG
  printf("***CLEAR***\n     ");
#endif
}

static int install(int *states, int size, int hvalue)
{ register int n, i;

  if (dalloc == SIGMA || settop+size >= setlimit)
    clear();

  n = dalloc++;
  hash[n].hlink  = bucket[hvalue];
  hash[n].set    = settop;
  bucket[hvalue] = hash + n;

  while (states >= sbot)
    *settop++ = *states--;
  *settop++ = azip;

  for (i = 0; i < SIGMA; i++)
    dfa[n][i] = 0;

  final[n] = -1;

  return (n);
}

static void unmark(register int *states)
{ while (states >= sbot)
    aNFA[*states--].mark = 0;
}

static void dfa_init(struct a_mach *mach)
{ clear();
  sbot   = mach->work1;
  dbot   = mach->work2;
  aNFA   = mach->start;
  azip   = mach->size + 1;
  aclass = mach->cvecs;
  hash[0].set  = setvec;
  setvec[0] = 0;
  setvec[1] = azip;
}


/* RE_FREE:
     Free regexp data type.
*/

static struct a_mach *amach = NULL;

void re_free(regexp rexpr)
{ struct a_mach *mach;

  mach = (struct a_mach *) rexpr;
  if (mach == amach) amach = NULL;
  FREE(mach->work2);
}

/* RE_MATCH:
     Search for a match to regular expression.
*/

int re_match(regexp rexpr, char *string, int thresh)
{ register int *states;
  register int c, s, n;
  int hvalue, size;
  struct a_mach *mach;
#ifdef DEBUG
  register int *ind;
#endif

  if (thresh > 0) return (approx_match(rexpr,string,thresh));

  mach  = (struct a_mach *) rexpr;

  if (mach != amach)
    { dfa_init(mach);
      amach = mach;
    }

#ifdef DEBUG
  printf("\nSTR:     0\n");
#endif

  s = 0;
  c = 1;
  while (1)
    { n = dfa[s][c];
      if (n == 0)
        { states = goset(s,c,&hvalue);
#ifdef DEBUG
          if (c == 1)
            printf("\\01: ");
          else if (c == 0)
            printf(" \\0: ");
          else if (c == '\n')
            printf(" \\n: ");
          else
            printf("  %c: ",c);
#endif
          size = states-sbot+1;
          n = hsearch(hvalue,size);
          if (n >= 0)
#ifdef DEBUG
          { printf("->%3d {%3d}\n",n,hvalue);
#endif
            dfa[s][c] = n;
#ifdef DEBUG
          }
#endif
          else
            { n = install(states,size,hvalue);
              if (s < n) dfa[s][c] = n;
#ifdef DEBUG
              printf("->%3d {%3d} ",n,hvalue);
              s = '[';
              for (ind = sbot; ind <= states; ind++)
                { printf("%c%3d",s,*ind);
                  s = ' ';
                }
              printf("]\n");
#endif
            }
          unmark(states);
        }
#ifdef DEBUG
      else
        { if (c == 1)
            printf("\\01: ");
          else if (c == 0)
            printf(" \\0: ");
          else if (c == '\n')
            printf(" \\n: ");
          else
            printf("  %c: ",c);
          printf("  %3d\n",n);
        }
#endif
      if (final[s]) return (1);
      s = n;
      if (c == 0) break;
      c = *string++;
    }
  if (final[s] < 0)
    { states = goset(s,0,&hvalue);
      unmark(states);
    }
  return (final[s]);
}


/* DEBUG routines for displaying parse trees and e-nfas */

#ifdef DEBUG

static char graph[80];

static char *sym[] = {"__OR","_CAT","PLUS","STAR","_OPT","ATMS","ATMC"};

static void _display(int ind, int sec, struct node *p)
{ printf("%s_%s ",graph,sym[p->label]);
  if (p->label == ATOM_S)
    printf("%c",((struct leaf *) p)->vecoff);
  else if (p->label == ATOM_C)
    { int x;
      for (x = 15; x >= 0; x--)
        printf("%02x",classes[((struct leaf *) p)->vecoff + x] & 0xff);
    }
  printf("\n");
  if (sec) graph[ind] = ' ';
  strcat(graph," |");
  if (p->label <= CAT)
    _display(ind+2,0,p->right);
  if (p->label <= OPT)
    _display(ind+2,1,p->left);
  graph[ind+1] = '\0';
}

static void show_tree(struct node *p)
{ printf("PARSE TREE:\n");
  printf("\n");
  strcpy(graph,"|");
  _display(0,1,p);
  printf("\n");
  fflush(stdout);
}

static void show_nfa(struct a_mach *mach)
{ register int i, x;

  printf("MACHINE:\n");
  printf("\n");
  printf("State table\n");
  for (i = 0; i <= size; i++)
    { printf("%2d: p=%2d %2d s=%2d %2d (%d): ",
            i,mach->start[i].prd1,mach->start[i].prd2,
              mach->start[i].suc1,mach->start[i].suc2,mach->start[i].kind);
      if (mach->start[i].kind == SYMBOL)
        printf("%c",mach->start[i].data);
      else if (mach->start[i].kind == CLASS)
        for (x = 15; x >= 0; x--)
          printf("%02x",mach->cvecs[mach->start[i].data + x] & 0xff);
      printf("\n");
      fflush(stdout);
    }
  printf("\n");
  printf("Loop table\n");
  for (i = 1; i <= back; i++)
    printf("%2d: %2d\n",i,mach->loop[i]);
  printf("\n");
  fflush(stdout);
}

static void show_dfa()
{ register int i, c, n;
  register int *s;
  register struct chain *b;

  printf("\nHASH TABLE\n");
  for (i = 0; i < HSIZE; i++)
    if (bucket[i] != 0)
      { printf("%3d:",i);
        for (b = bucket[i]; b != 0; b = b->hlink)
          printf(" %3d",(int) (b-hash));
        printf("\n");
      }
  fflush(stdout);
  printf("\nDFA ARRAY\n");
  for (i = 0; i < dalloc; i++)
    { printf("%3d:",i);
      for (c = 0; c < SIGMA; c++)
        if ((n=dfa[i][c]) != 0)
          if (isprint(c))
            printf(" (%c,%d)",c,n);
          else
            printf(" (\\%o,%d)",c,n);
      printf("\n");
      fflush(stdout);
      printf("    ");
      for (s = hash[i].set; *s != azip; s++)
        printf(" %3d",*s);
      printf("\n");
      fflush(stdout);
    }
  fflush(stdout);
}

#define BUFLEN 2048

static char buf[BUFLEN];

void main(int argc, char *argv[])
{ regexp mach;
  FILE *ifile;
  int thresh;
  int level;
  int match;
  char message[100];

  mach = re_parse(0,argv[1],message,&level);
  if (mach == NULL)
    { fprintf(stderr,"%d: %s\n",level,message);
      exit (1);
    }

  ifile  = fopen(argv[2],"r");
  thresh = atoi(argv[3]);

  while (fgets(buf,BUFLEN,ifile) != NULL)
    { buf[strlen(buf)-1] = '\0';
      match = re_match(mach,buf,thresh);
      if (match)
        printf("***MATCH***\n");
      else
        printf("***NO MATCH***\n");
    }

  if (thresh == 0) show_dfa();

  re_free(mach);

  exit (0);
}

#endif
