#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "bri.h"

typedef struct pallocroot pallocroot;
typedef struct pallocnode pallocnode;

//  _dbg: 0 -- print nothing
//        1 -- print block allocation
//        2 -- print all allocations

struct pallocroot {
  size_t       _bs;  //  size of block
  pallocnode  *_nl;  //  nodeList
  pallocnode  *_cn;  //  currentNode
  int          _dbg; //  if set, debug information is printed
};

struct pallocnode {
  size_t        _cp;  //  cuurentPosition
  char         *_dt;  //  data
  pallocnode   *_nx;  //  next pallocnode
};

extern pallocroot _palloc_stuff;

pallocroot  _palloc_stuff = { 128 * 1024 * 1024, NULL, NULL, 0 };

static
void *
really_allocate(size_t size) {
  errno = 0;
  void *ret = malloc(size);
  if (ret == 0L) {
    fprintf(stderr, "palloc()-- can't allocate %lu bytes: %s.\n", size, strerror(errno));
    exit(1);
  }
  return(ret);
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
psetdebug(int on) {
  _palloc_stuff._dbg = on;
}

void
pfree(void) {
  pallocnode *n;
  size_t      r = 0;
  size_t      b = 0;

  while ((n = _palloc_stuff._nl) != 0L) {
    r += n->_cp;
    b++;
    _palloc_stuff._nl = n->_nx;
    free(n->_dt);
    free(n);
  }

  if (_palloc_stuff._dbg > 0)
    fprintf(stderr, "palloc()-- %lu bytes in %lu blocks returned to free store.\n", r, b);

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


  //  Allocate the initial block if it doesn't exist.
  //
  if (_palloc_stuff._nl == NULL) {
    _palloc_stuff._nl = (pallocnode *)really_allocate(sizeof(pallocnode));
    _palloc_stuff._cn = _palloc_stuff._nl;

    if (_palloc_stuff._dbg > 0)
      fprintf(stderr, "palloc()-- Inital block of %lu bytes at %p.\n", _palloc_stuff._bs, _palloc_stuff._cn);

    _palloc_stuff._cn->_cp = 0;
    _palloc_stuff._cn->_dt = (char *)really_allocate(_palloc_stuff._bs);
    _palloc_stuff._cn->_nx = NULL;
  }


  //  If the requested space is larger than our block size, allocate a
  //  new node with the required amount of space.  The new node is
  //  placed on the start of the alloc'd list.
  //
  //  We also place blocks that are bigger than the amount free in the
  //  current block, AND bigger than the amount used in the current
  //  block here.  Since the new block is larger than the free space,
  //  it won't fit in the current block.  Since the new block is
  //  larger than the current block, it is wasteful to throw out the
  //  current block and replace it with a new block.
  //
  //  The tests read:
  //    new block is bigger than our block size
  //    new block won't fit in current block
  //    new block is larger than current block
  //
  if ((size > _palloc_stuff._bs) ||
      ((size > _palloc_stuff._bs - _palloc_stuff._cn->_cp) &&
       (size > _palloc_stuff._cn->_cp))) {
    pallocnode *n;

    n = (pallocnode *)really_allocate(sizeof(pallocnode));
    n->_cp = size;
    n->_dt = (char *)really_allocate(size);
    n->_nx = _palloc_stuff._nl;

    if (_palloc_stuff._dbg > 0)
      fprintf(stderr, "palloc()-- New needs %lu bytes: custom new block at %p.\n",
              size,
              n);

    _palloc_stuff._nl = n;
    if (_palloc_stuff._cn == 0L)
      _palloc_stuff._cn = n;

    return(n->_dt);
  }


  //  Need more space?
  //
  if (size + _palloc_stuff._cn->_cp > _palloc_stuff._bs) {
    _palloc_stuff._cn->_nx = (pallocnode *)really_allocate(sizeof(pallocnode));

    if (_palloc_stuff._dbg > 0)
      fprintf(stderr, "palloc()-- Old block %.3f%% used (%lu bytes remaining), new needs %lu bytes: new block of %lu bytes at %p.\n",
              100.0 * _palloc_stuff._cn->_cp / _palloc_stuff._bs,
              _palloc_stuff._bs - _palloc_stuff._cn->_cp,
              size,
              _palloc_stuff._bs,
              _palloc_stuff._cn->_nx);

    _palloc_stuff._cn      = _palloc_stuff._cn->_nx;
    _palloc_stuff._cn->_cp = 0;
    _palloc_stuff._cn->_dt = (char *)really_allocate(_palloc_stuff._bs);
    _palloc_stuff._cn->_nx = NULL;
  }


  //  OK, grab the space, and return it.
  //
  _palloc_stuff._cn->_cp += size;

  if (_palloc_stuff._dbg > 1)
    fprintf(stderr, "palloc()-- Old block %.3f%% used (%lu bytes remaining): returning %lu bytes.\n",
              100.0 * _palloc_stuff._cn->_cp / _palloc_stuff._bs,
            _palloc_stuff._bs - _palloc_stuff._cn->_cp,
            size);

  return(_palloc_stuff._cn->_dt + _palloc_stuff._cn->_cp - size);
}


void
pdumppalloc(void) {
  pallocnode *n = _palloc_stuff._nl;
  fprintf(stderr, "palloc dump\n");
  fprintf(stderr, "%lu bytes per block\n", _palloc_stuff._bs);
  while (n != 0L) {
    fprintf(stderr, "%p: currentPosition: %8ld bytes used%s\n",
            n, n->_cp, (n == _palloc_stuff._cn) ? ", current block" : "");
    n = n->_nx;
  }
}

