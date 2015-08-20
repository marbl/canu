
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2008-JUN-09 to 2008-JUL-24
 *      are Copyright 2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

//static char sccsid[] = "@(#)qsort.c	8.1 (Berkeley) 6/4/93";
//__FBSDID("$FreeBSD: src/lib/libc/stdlib/qsort.c,v 1.12 2002/09/10 02:04:49 wollman Exp $");

//#include <sys/cdefs.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include <sys/types.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#ifdef __FreeBSD__
#include <pmc.h>
#endif

typedef int		 cmp_t(const void *, const void *);

static inline char	*med3(char *, char *, char *, cmp_t *);
static inline void	 swapfunc(char *, char *, int, int);

#define min(a, b)	(a) < (b) ? a : b

/*
 * Qsort routine from Bentley & McIlroy's "Engineering a Sort Function".
 */
#define swapcode(TYPE, parmi, parmj, n) {       \
    long i = (n) / sizeof (TYPE);               \
    TYPE *pi = (TYPE *) (parmi); 		\
    TYPE *pj = (TYPE *) (parmj); 		\
    do {                                        \
      TYPE	t = *pi;                        \
      *pi++ = *pj;				\
      *pj++ = t;				\
    } while (--i > 0);				\
  }


static inline void
swapfunc(a, b, n, swaptype)
     char *a, *b;
     int n, swaptype;
{
  if(swaptype <= 1)
    swapcode(long, a, b, n)
    else
      swapcode(char, a, b, n)
        }

#define swap(a, b)                              \
  if (swaptype == 0) {				\
    long t = *(long *)(a);			\
    *(long *)(a) = *(long *)(b);		\
    *(long *)(b) = t;                           \
  } else                                        \
    swapfunc(a, b, es, swaptype)

#define vecswap(a, b, n) 	if ((n) > 0) swapfunc(a, b, n, swaptype)

#define	CMP(x, y) (cmp((x), (y)))

static inline char *
med3(char *a, char *b, char *c, cmp_t *cmp)
{
  return CMP(a, b) < 0 ?
    (CMP(b, c) < 0 ? b : (CMP(a, c) < 0 ? c : a ))
    :(CMP(b, c) > 0 ? b : (CMP(a, c) < 0 ? a : c ));
}

/*
 * We use some elaborate condition variables and signalling
 * to ensure a bound of the number of active threads at
 * 2 * maxthreads and the size of the thread data structure
 * to maxthreads.
 */

/* Condition of starting a new thread. */
enum thread_state {
  ts_idle,		/* Idle, waiting for instructions. */
  ts_work,		/* Has work to do. */
  ts_term			/* Asked to terminate. */
};

/* Variant part passed to qsort invocations. */
struct qsort {
  enum thread_state st;	/* For coordinating work. */
  struct common *common;	/* Common shared elements. */
  void *a;		/* Array base. */
  size_t n;		/* Number of elements. */
  pthread_t id;		/* Thread id. */
  pthread_mutex_t mtx_st;	/* For signalling state change. */
  pthread_cond_t cond_st;	/* For signalling state change. */
};

/* Invariant common part, shared across invocations. */
struct common {
  int swaptype;		/* Code to use for swapping */
  size_t es;		/* Element size. */
  cmp_t *cmp;		/* Comparison function */
  int nthreads;		/* Total number of pool threads. */
  int idlethreads;	/* Number of idle threads in pool. */
  int forkelem;		/* Minimum number of elements for a new thread. */
  struct qsort *pool;	/* Fixed pool of threads. */
  pthread_mutex_t mtx_al;	/* For allocating threads in the pool. */
};

static void *qsort_thread(void *p);

/* The multithreaded qsort public interface */

void
qsort_mt(void *a, size_t n, size_t es, cmp_t *cmp, int maxthreads, int forkelem)
{
  struct qsort *qs;
  struct common c;
  int i, islot;
  int    bailout = 1;

  if (n < forkelem)
    goto f1;
  errno = 0;

  if (maxthreads <= 1)
    goto f1;

  /* Try to initialize the resources we need. */
  if (pthread_mutex_init(&c.mtx_al, NULL) != 0)
    goto f1;
  if ((c.pool = (struct qsort *)calloc(maxthreads, sizeof(struct qsort))) ==NULL)
    goto f2;
  for (islot = 0; islot < maxthreads; islot++) {
    qs = &c.pool[islot];
    if (pthread_mutex_init(&qs->mtx_st, NULL) != 0)
      goto f3;
    if (pthread_cond_init(&qs->cond_st, NULL) != 0) {
      pthread_mutex_destroy(&qs->mtx_st);
      goto f3;
    }
    qs->st = ts_idle;
    qs->common = &c;
    if (pthread_create(&qs->id, NULL, qsort_thread, qs) != 0) {
      pthread_mutex_destroy(&qs->mtx_st);
      pthread_cond_destroy(&qs->cond_st);
      goto f3;
    }
  }

  /* All systems go. */
  bailout = 0;

  /* Initialize common elements. */
  c.swaptype = ((char *)a - (char *)0) % sizeof(long) || \
    es % sizeof(long) ? 2 : es == sizeof(long)? 0 : 1;
  c.es = es;
  c.cmp = cmp;
  c.forkelem = forkelem;
  c.idlethreads = c.nthreads = maxthreads;

  /* Hand out the first work batch. */
  qs = &c.pool[0];
  pthread_mutex_lock(&qs->mtx_st);
  qs->a = a;
  qs->n = n;
  qs->st = ts_work;
  c.idlethreads--;
  pthread_cond_signal(&qs->cond_st);
  pthread_mutex_unlock(&qs->mtx_st);

  /*
   * Wait for all threads to finish, and
   * free acquired resources.
   */
 f3:	for (i = 0; i < islot; i++) {
    qs = &c.pool[i];
    if (bailout) {
      pthread_mutex_lock(&qs->mtx_st);
      qs->st = ts_term;
      pthread_cond_signal(&qs->cond_st);
      pthread_mutex_unlock(&qs->mtx_st);
    }
    pthread_join(qs->id, NULL);
    pthread_mutex_destroy(&qs->mtx_st);
    pthread_cond_destroy(&qs->cond_st);
  }
  free(c.pool);
 f2:	pthread_mutex_destroy(&c.mtx_al);
  if (bailout) {
    /* XXX should include a syslog call here */
    fprintf(stderr, "Resource initialization failed; bailing out.\n");
  f1:		qsort(a, n, es, cmp);
  }
}


/*
 * Allocate an idle thread from the pool, lock its
 * mutex, change its state to work, decrease the number
 * of idle threads, and return a
 * pointer to its data area.
 * Return NULL, if no thread is available.
 */
static struct qsort *
allocate_thread(struct common *c)
{
  int i;

  pthread_mutex_lock(&c->mtx_al);
  for (i = 0; i < c->nthreads; i++)
    if (c->pool[i].st == ts_idle) {
      c->idlethreads--;
      c->pool[i].st = ts_work;
      pthread_mutex_lock(&c->pool[i].mtx_st);
      pthread_mutex_unlock(&c->mtx_al);
      return (&c->pool[i]);
    }
  pthread_mutex_unlock(&c->mtx_al);
  return (NULL);
}

/* Thread-callable quicksort. */
static void
qsort_algo(struct qsort *qs)
{
  char *pa, *pb, *pc, *pd, *pl, *pm, *pn;
  long d, r, swaptype, swap_cnt;
  void *a;			/* Array of elements. */
  size_t n, es;			/* Number of elements; size. */
  cmp_t *cmp;
  long nl, nr;
  struct common *c;
  struct qsort *qs2;
  pthread_t id;

  /* Initialize qsort arguments. */
  id = qs->id;
  c = qs->common;
  es = c->es;
  cmp = c->cmp;
  swaptype = c->swaptype;
  a = qs->a;
  n = qs->n;
 top:

  /* From here on qsort(3) business as usual. */
  swap_cnt = 0;
  if (n < 7) {
    for (pm = (char *)a + es; pm < (char *)a + n * es; pm += es)
      for (pl = pm;
           pl > (char *)a && CMP(pl - es, pl) > 0;
           pl -= es)
        swap(pl, pl - es);
    return;
  }
  pm = (char *)a + (n / 2) * es;
  if (n > 7) {
    pl = a;
    pn = (char *)a + (n - 1) * es;
    if (n > 40) {
      d = (n / 8) * es;
      pl = med3(pl, pl + d, pl + 2 * d, cmp);
      pm = med3(pm - d, pm, pm + d, cmp);
      pn = med3(pn - 2 * d, pn - d, pn, cmp);
    }
    pm = med3(pl, pm, pn, cmp);
  }
  swap(a, pm);
  pa = pb = (char *)a + es;

  pc = pd = (char *)a + (n - 1) * es;
  for (;;) {
    while (pb <= pc && (r = CMP(pb, a)) <= 0) {
      if (r == 0) {
        swap_cnt = 1;
        swap(pa, pb);
        pa += es;
      }
      pb += es;
    }
    while (pb <= pc && (r = CMP(pc, a)) >= 0) {
      if (r == 0) {
        swap_cnt = 1;
        swap(pc, pd);
        pd -= es;
      }
      pc -= es;
    }
    if (pb > pc)
      break;
    swap(pb, pc);
    swap_cnt = 1;
    pb += es;
    pc -= es;
  }
  if (swap_cnt == 0) {  /* Switch to insertion sort */
    for (pm = (char *)a + es; pm < (char *)a + n * es; pm += es)
      for (pl = pm;
           pl > (char *)a && CMP(pl - es, pl) > 0;
           pl -= es)
        swap(pl, pl - es);
    return;
  }

  pn = (char *)a + n * es;
  r = min(pa - (char *)a, pb - pa);
  vecswap(a, pb - r, r);
  r = min(pd - pc, pn - pd - es);
  vecswap(pb, pn - r, r);

  nl = (pb - pa) / es;
  nr = (pd - pc) / es;

  /* Now try to launch subthreads. */
  if (nl > c->forkelem && nr > c->forkelem &&
      (qs2 = allocate_thread(c)) != NULL) {
    qs2->a = a;
    qs2->n = nl;
    pthread_cond_signal(&qs2->cond_st);
    pthread_mutex_unlock(&qs2->mtx_st);
  } else if (nl > 0) {
    qs->a = a;
    qs->n = nl;
    qsort_algo(qs);
  }
  if (nr > 0) {
    a = pn - nr * es;
    n = nr;
    goto top;
  }
}

/* Thread-callable quicksort. */
static void *
qsort_thread(void *p)
{
  struct qsort *qs, *qs2;
  int i;
  struct common *c;
  pthread_t id;

  qs = p;
  id = qs->id;
  c = qs->common;
 again:
  /* Wait for work to be allocated. */
  pthread_mutex_lock(&qs->mtx_st);
  while (qs->st == ts_idle)
    pthread_cond_wait(&qs->cond_st, &qs->mtx_st);
  pthread_mutex_unlock(&qs->mtx_st);
  if (qs->st == ts_term) {
    return(NULL);
  }
  assert(qs->st == ts_work);

  qsort_algo(qs);

  pthread_mutex_lock(&c->mtx_al);
  qs->st = ts_idle;
  c->idlethreads++;
  if (c->idlethreads == c->nthreads) {
    for (i = 0; i < c->nthreads; i++) {
      qs2 = &c->pool[i];
      if (qs2 == qs)
        continue;
      pthread_mutex_lock(&qs2->mtx_st);
      qs2->st = ts_term;
      pthread_cond_signal(&qs2->cond_st);
      pthread_mutex_unlock(&qs2->mtx_st);
    }
    pthread_mutex_unlock(&c->mtx_al);
    return(NULL);
  }
  pthread_mutex_unlock(&c->mtx_al);
  goto again;
}
