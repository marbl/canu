
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
 *    Brian P. Walenz from 2003-JAN-02 to 2004-MAY-06
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-JAN-13 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-05 to 2011-JAN-10
 *      are Copyright 2005-2008,2011 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "util.h"

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
  void *ret = malloc(size);
  if (ret == 0L) {
    fprintf(stderr, "palloc()-- can't allocate "sizetFMT" bytes: %s.\n", size, strerror(errno));
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

void*
pallochandle(size_t size) {
  pallocroot *root = (pallocroot *)malloc(sizeof(pallocroot));
  if (root == NULL)
    fprintf(stderr, "pallochandle()-- can't allocate a handle!\n"), exit(1);
  if (size == 0)
    size = 128 * 1024 * 1024;
  root->_bs  = size;
  root->_nl  = NULL;
  root->_cn  = NULL;
  root->_dbg = 0;
  return(root);
}


//  Release a palloc handle, does not release the memory in the handle!
void
pfreehandle(void *handle) {
  free((pallocroot *)handle);
}

//  Clear out memory inside the handle.  The handle remains valid after this.
void
pfree2(void *handle) {
  pallocroot  *root = (pallocroot *)handle;
  pallocnode *n;
  size_t      r = 0;
  size_t      b = 0;

  if (root == NULL)
    root = &_palloc_stuff;

  while ((n = root->_nl) != 0L) {
    r += n->_cp;
    b++;
    root->_nl = n->_nx;
    free(n->_dt);
    free(n);
  }

  if (root->_dbg > 0)
    fprintf(stderr, "palloc()-- "sizetFMT" bytes in "sizetFMT" blocks returned to free store.\n", r, b);

  root->_nl = 0L;
  root->_cn = 0L;
}

void
pfree(void) {
  pfree2(&_palloc_stuff);
}


void *
palloc2(size_t size, void *handle) {
  pallocroot  *root = (pallocroot *)handle;

  if (root == NULL)
    root = &_palloc_stuff;

  //  Make size a multiple of 8
  //
  if (size & 0x7) {
    size >>= 3;
    size++;
    size <<= 3;
  }
  if (size == 0)
    return(0L);

  //  Allocate the initial block if it doesn't exist.
  //
  if (root->_nl == NULL) {
    root->_nl = (pallocnode *)really_allocate(sizeof(pallocnode));
    root->_cn = root->_nl;

    if (root->_dbg > 0)
      fprintf(stderr, "palloc()-- Inital block of "sizetFMT" bytes at %p.\n", root->_bs, root->_cn);

    root->_cn->_cp = 0;
    root->_cn->_dt = (char *)really_allocate(root->_bs);
    root->_cn->_nx = NULL;
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
  if ((size > root->_bs) ||
      ((size > root->_bs - root->_cn->_cp) &&
       (size > root->_cn->_cp))) {
    pallocnode *n;

    n = (pallocnode *)really_allocate(sizeof(pallocnode));
    n->_cp = size;
    n->_dt = (char *)really_allocate(size);
    n->_nx = root->_nl;

    if (root->_dbg > 0)
      fprintf(stderr, "palloc()-- New needs "sizetFMT" bytes: custom new block at %p.\n",
              size,
              n);

    root->_nl = n;
    if (root->_cn == 0L)
      root->_cn = n;

    return(n->_dt);
  }


  //  Need more space?
  //
  if (size + root->_cn->_cp > root->_bs) {
    root->_cn->_nx = (pallocnode *)really_allocate(sizeof(pallocnode));

    if (root->_dbg > 0)
      fprintf(stderr, "palloc()-- Old block %.3f%% used ("sizetFMT" bytes remaining), new needs "sizetFMT" bytes: new block of "sizetFMT" bytes at %p.\n",
              100.0 * root->_cn->_cp / root->_bs,
              root->_bs - root->_cn->_cp,
              size,
              root->_bs,
              root->_cn->_nx);

    root->_cn      = root->_cn->_nx;
    root->_cn->_cp = 0;
    root->_cn->_dt = (char *)really_allocate(root->_bs);
    root->_cn->_nx = NULL;
  }


  //  OK, grab the space, and return it.
  //
  root->_cn->_cp += size;

  if (root->_dbg > 1)
    fprintf(stderr, "palloc()-- Old block %.3f%% used ("sizetFMT" bytes remaining): returning "sizetFMT" bytes at %p.\n",
              100.0 * root->_cn->_cp / root->_bs,
            root->_bs - root->_cn->_cp,
            size, root->_cn->_dt + root->_cn->_cp - size);

  return(root->_cn->_dt + root->_cn->_cp - size);
}



void *
palloc(size_t size) {
  return(palloc2(size, &_palloc_stuff));
}


void
pdumppalloc(void *handle) {
  pallocroot  *root = (pallocroot *)handle;
  pallocnode *n = root->_nl;
  fprintf(stderr, "palloc dump\n");
  fprintf(stderr, ""sizetFMT" bytes per block\n", root->_bs);
  while (n != 0L) {
    fprintf(stderr, "%p: currentPosition: "sizetFMT" bytes used%s\n",
            n, n->_cp, (n == root->_cn) ? ", current block" : "");
    n = n->_nx;
  }
}

