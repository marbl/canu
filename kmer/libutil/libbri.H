#ifndef LIBBRI_H
#define LIBBRI_H

#ifdef __apple
inline
void
nanosleep(struct timespec ts, struct timespec tr) {
  usleep(ts.tv_nsec >> 10);
}
#endif

#include "alphabet.h"
#include "britime.H"
#include "file.h"
#include "stat.h"
#include "memory.h"
#include "palloc.h"
#include "md5.h"

#include "readBuffer.H"
#include "splittowords.H"
#include "bit-operations.H"
#include "bit-packing.H"
#include "bitPackedFile.H"

#include "fastastream.H"
#include "merstream.H"

#endif  //  LIBBRI_H
