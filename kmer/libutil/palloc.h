#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

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


//  Define PALLOC_DEBUG   to get message whenever the state of palloc's memory map changes
//  Define PALLOC_VERBOSE to get messages about space returned and remaining (one for each palloc call)
//  Define PALLOC_TEST    to include a short main() in palloc.c
//
//#define PALLOC_DEBUG
//#define PALLOC_VERBOSE
//#define PALLOC_TEST



#ifdef __cplusplus
}
#endif
