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

#ifdef __apple
inline
void
nanosleep(struct timespec ts, struct timespec tr) {
  usleep(ts.tv_nsec >> 10);
}
#endif

#include "file.h"
#include "stat.h"
#include "memory.h"
#include "palloc.h"
#include "time.H"
#include "md5.h"

#include "readBuffer.H"
#include "splittowords.H"
#include "bit-operations.H"
#include "bit-packing.H"
#include "bitPackedFile.H"

#include "fastastream.H"
#include "merstream.H"

#endif  //  LIBBRI_H
