#include <stdio.h>
#include <stdlib.h>

#include "bri.h"



//  Define PALLOC_DEBUG   to get message whenever the state of palloc's memory map changes
//  Define PALLOC_VERBOSE to get messages about space returned and remaining (one for each palloc call)
//  Define PALLOC_TEST    to include a short main() in palloc.c
//
//#define PALLOC_DEBUG
//#define PALLOC_VERBOSE
//#define PALLOC_TEST




typedef struct pallocroot pallocroot;
typedef struct pallocnode pallocnode;

struct pallocroot {
  size_t       _bs;  //  number of blocks per data element
  pallocnode  *_nl;  //  nodeList
  pallocnode  *_cn;  //  currentNode
};

struct pallocnode {
  size_t        _cp;  //  cuurentPosition
  char         *_dt;  //  data
  pallocnode   *_nx;  //  next pallocnode
};

extern pallocroot _palloc_stuff;


pallocroot  _palloc_stuff = { 128 * 1024 * 1024, NULL, NULL };


static
void *
really_allocate(size_t size) {
  return(malloc(size));
}

void
psetblocksize(size_t size) {
  if (_palloc_stuff._nl == 0L)
    _palloc_stuff._bs = size;
}

size_t
pgetblocksize(void) {
  return(_palloc_stuff._bs);
}


void
pfree(void) {
  pallocnode *n;

  while ((n = _palloc_stuff._nl) != 0L) {
    _palloc_stuff._nl = n->_nx;
    free(n->_dt);
    free(n);
  }

  _palloc_stuff._nl = 0L;
  _palloc_stuff._cn = 0L;
}


void *
palloc(size_t size) {

  //  Make size a multiple of 8
  //
  if (size & 0x7) {
    size >>= 3;
    size++;
    size <<= 3;
  }

  //  If the requested space is larger than our block size, allocate a
  //  new node with the required amount of space.  The new node is
  //  placed on the start of the alloc'd list.
  //
  if (size > _palloc_stuff._bs) {
    pallocnode *n;

    n = (pallocnode *)really_allocate(sizeof(pallocnode));
    n->_cp = size;
    n->_dt = (char *)really_allocate(size);
    n->_nx = _palloc_stuff._nl;

#ifdef PALLOC_DEBUG
    fprintf(stderr, "palloc()-- Custom block for %d bytes at %p.\n", size, n);
#endif

    _palloc_stuff._nl = n;
    if (_palloc_stuff._cn == 0L)
      _palloc_stuff._cn = n;

    return(n->_dt);
  }

  //  Allocate the initial stuff
  //
  if (_palloc_stuff._nl == NULL) {
    _palloc_stuff._nl = (pallocnode *)really_allocate(sizeof(pallocnode));
    _palloc_stuff._cn = _palloc_stuff._nl;

#ifdef PALLOC_DEBUG
    fprintf(stderr, "palloc()-- Inital block of %d bytes at %p.\n", _palloc_stuff._bs, _palloc_stuff._cn);
#endif

    _palloc_stuff._cn->_cp = 0;
    _palloc_stuff._cn->_dt = (char *)really_allocate(_palloc_stuff._bs);
    _palloc_stuff._cn->_nx = NULL;
  }

  //  Need more space?
  //
  if (size + _palloc_stuff._cn->_cp > _palloc_stuff._bs) {
    _palloc_stuff._cn->_nx = (pallocnode *)really_allocate(sizeof(pallocnode));
    _palloc_stuff._cn      = _palloc_stuff._cn->_nx;

#ifdef PALLOC_DEBUG
    fprintf(stderr, "palloc()-- New block or %d bytes at %p.\n", _palloc_stuff._bs, _palloc_stuff._cn);
#endif

    _palloc_stuff._cn->_cp = 0;
    _palloc_stuff._cn->_dt = (char *)really_allocate(_palloc_stuff._bs);
    _palloc_stuff._cn->_nx = NULL;
  }

  //  OK, grab the space, and return it.
  //
  _palloc_stuff._cn->_cp += size;

#ifdef PALLOC_VERBOSE
  fprintf(stderr, "palloc()--   returning %d bytes; %d left.\n",
          size, 
          _palloc_stuff._bs - _palloc_stuff._cn->_cp);
#endif

  return(_palloc_stuff._cn->_dt + _palloc_stuff._cn->_cp - size);
}


#ifdef PALLOC_TEST
static
void
dumppalloc(void) {
  pallocnode *n = _palloc_stuff._nl;

  fprintf(stderr, "0x%016lx bs=%d\n", _palloc_stuff._cn, _palloc_stuff._bs);
  while (n != 0L) {
    fprintf(stderr, "0x%016x: %d\n", n, n->_cp);
    n = n->_nx;
  }
}

int
main(int argc, char **argv) {

  psetblocksize(1024);

  (void)palloc(2048);
  (void)palloc(128);
  (void)palloc(999);
  (void)palloc(8);
  (void)palloc(8);
  (void)palloc(8);
  (void)palloc(8);
  (void)palloc(2056);
  (void)palloc(8);
  (void)palloc(2064);
  (void)palloc(8);
  (void)palloc(2072);
  (void)palloc(8);

  (void)dumppalloc();

  pfree();

  fprintf(stderr, "----------------------------------------\n");

  psetblocksize(10240);

  (void)palloc(2048);
  (void)palloc(128);
  (void)palloc(999);
  (void)palloc(8);
  (void)palloc(8);
  (void)palloc(8);
  (void)palloc(8);
  (void)palloc(2056);
  (void)palloc(8);
  (void)palloc(2064);
  (void)palloc(8);
  (void)palloc(2072);
  (void)palloc(8);

  (void)dumppalloc();

  pfree();

  return(0);
}
#endif
