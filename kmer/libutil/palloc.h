#include <stdio.h>
#include <stdlib.h>

//  Pac-Man's memory allocator.
//
//  Grabs big chunks of memory, then gives out little pieces.  You can
//  only free ALL memory, not single blocks.
//
//  This is useful when one needs to malloc() tens of millions of
//  things, at which point the overhead of finding a free block is
//  large.
//
void   *palloc(size_t size);
void    pfree(void);

//  The block size can only be changed before the first call to
//  palloc().  Calling psetblocksize() after that has no effect.
//
void    psetblocksize(size_t size);
size_t  pgetblocksize(void);

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

