#include <pthread.h>

//
//  Sun5.6 and Linux do not work if we define this.
//  BSD and Tru64 don't work if this isn't defined.
//

#if defined(__bsd) || defined(__alpha) || defined(_AIX)

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199506L
#endif

#endif
