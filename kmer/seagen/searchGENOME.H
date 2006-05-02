#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

//  At one time, this was needed for pthread.h or semaphore.h
//typedef unsigned short ushort;

#include <pthread.h>
#include <semaphore.h>

#include "util++.H"
#include "bio++.H"
#include "positionDB.H"
#include "existDB.H"

#include "hitMatrix.H"
#include "searcherState.H"

#include "configuration.H"

