#ifndef LIBBRI_H
#define LIBBRI_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#include "alphabet.h"

//  Create the O_LARGEFILE type for open(), if it doesn't already
//  exist (FreeBSD, Tru64)
//
#ifndef O_LARGEFILE
#define O_LARGEFILE    0
#endif


#ifdef __apple
inline
void
nanosleep(struct timespec ts, struct timespec tr) {
  usleep(ts.tv_nsec >> 10);
}
#endif


#include "file.h"
#include "palloc.h"
#include "time.H"

#include "readBuffer.H"
#include "splittowords.H"
#include "bit-operations.H"
#include "bit-packing.H"
#include "bitPackedFile.H"

#include "sim4reader.h"
#include "fastastream.H"
#include "merstream.H"



#endif  //  LIBBRI_H
