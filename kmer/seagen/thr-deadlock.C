#include "searchGENOME.H"

//  OSF/1 on Compaq Alpha has, in the past, gotten stuck in a deadlock
//  situation allocating memory.  There's lots of debugging stuff at
//  the end if this file.

#ifdef __alpha

//  Define this to kill the process with a vengance instead of
//  gracefully exiting.  exit() tries to free memory, and is thus gets
//  caught in the deadlock -- but is useful for debugging.
//
#define KILL_INSTEAD_OF_EXIT

#ifdef KILL_INSTEAD_OF_EXIT
#include <signal.h>
#endif

u32bit deadlockTested = 0;
u32bit deadlockPassed = 0;

void*
deadlockDetector(void *) {

  fprintf(stderr, "Hello!  I'm a deadlockDetector!\n");

 detectAgain:

  //  Wait for the deadlock checker to reset things
  //
  while ((deadlockTested == 1) || (deadlockPassed == 1))
    sleep(4);

  deadlockTested = 1;
  char *x = new char [16];
  delete [] x;
  deadlockPassed = 1;

  goto detectAgain;

  return(0L);  //  Ignore the warning!
}

void*
deadlockChecker(void *) {

  fprintf(stderr, "Hello!  I'm a deadlockChecker!\n");

 checkAgain:

  //  Wait for the tester to test
  //
  while (deadlockTested == 0)
    sleep(5);

  //  Give it another ten seconds to return
  //
  sleep(5);

  if (deadlockPassed == 0) {
    fprintf(stderr, "\n\n\nESTmapper/search-- Deadlock detected!  Aborting the process!\n\n");
    fflush(stderr);
#ifdef KILL_INSTEAD_OF_EXIT
    kill(getpid(), SIGKILL);
#endif
    exit(1);
  }

  //fprintf(stderr, "Deadlock OK\n");

  //  Reset the testing/checking flags
  //
  deadlockPassed = 0;
  deadlockTested = 0;

  goto checkAgain;

  return(0L);  //  Ignore the warning!
}

#endif  //  _alpha










#ifdef DONT_EVER_ENABLE_THIS

//
//  Here are some notes on what was tried, and the stack trace from a lock.
//  This test failed to find the cause.
//

#define SIZE   (16 * 1024 * 1024)

void*
mallocStressor(void *) {
  struct timespec sleepAmt = { 0, 10000 };
  unsigned long   v = 0;

  fprintf(stderr, "Hello!  I'm a mallocStressor!\n");

 mallocAgain:
  //nanosleep(&sleepAmt, 0L);
  char *x = new char [SIZE];
  for (unsigned int i=SIZE; i--; )
    x[i] = i >> 5;
  for (unsigned int i=SIZE; i--; )
    x[i] |= x[SIZE-i];
  for (unsigned int i=SIZE; i--; )
    v    += x[i];
  delete [] x;

  goto mallocAgain;

  return((void*)v);  //  Ignore the warning!
}

void
main(int argc, char **argv) {
  pthread_attr_t   threadAttr;
  pthread_t        threadID;

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);

  pthread_create(&threadID, &threadAttr, deadlockDetector, 0L);
  pthread_create(&threadID, &threadAttr, deadlockChecker, 0L);

  for (unsigned int i=0; i<16; i++)
    pthread_create(&threadID, &threadAttr, mallocStressor, (void *)i);

  sleep(100);
}

//
//
//  Stack trace #1
//
//

(ladebug) show thread
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
*      1 default thread            blocked         kern usleep SCHED_OTHER  19
      -1 manager thread            blk SCS                     SCHED_RR     19
      -2 null thread for VP 2      running VP 2                null thread  -1
>      2 <anonymous>               blocked         kern usleep SCHED_OTHER  19
      -3 null thread for VP 3      running VP 3                null thread  -1
       3 <anonymous>               blocked         mut 9       SCHED_OTHER  19
      -4 null thread for VP 4      running VP 4                null thread  -1
       4 <anonymous>               blocked         mut 9       SCHED_OTHER  19
      -5 null thread for VP 5      running VP 5                null thread  -1
       5 <anonymous>               blocked         mut 9       SCHED_OTHER  19
       6 <anonymous>               blocked         mut 9       SCHED_OTHER  19
(ladebug) thread 6
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      6 <anonymous>               blocked         mut 9       SCHED_OTHER  19

(ladebug) where
>0  0x3ff805caf3c in __hstTransferRegisters(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x20001e2fc40) in /usr/shlib/libpthread.so
#1  0x3ff805af74c in __osTransferContext(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x20001e2fc40) in /usr/shlib/libpthread.so
#2  0x3ff805a3c50 in __dspTransferContext(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x20001e2fc40) in /usr/shlib/libpthread.so
#3  0x3ff805a12f4 in __dspDispatch(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x20001e2fc40) in /usr/shlib/libpthread.so
#4  0x3ff805ab90c in UnknownProcedure32FromFile8(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x20001e2fc40) in /usr/shlib/libpthread.so
#5  0x3ff805ab2f4 in UnknownProcedure31FromFile8(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x20001e2fc40) in /usr/shlib/libpthread.so
#6  0x3ff805abe3c in UnknownProcedure34FromFile8(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x20001e2fc40) in /usr/shlib/libpthread.so
#7  0x3ff801bf6a0 in UnknownProcedure2FromFile22(0x3ffc01b8400, 0x3ff805ae6d0, 0xfffffffffffffffc, 0x40100000, 0x222384c, 0x20001e2fc40) in /usr/shlib/lib
c.so
#8  0x3ff800cdad4 in malloc(0x3ffc01b8400, 0x3ff805ae6d0, 0xfffffffffffffffc, 0x40100000, 0x222384c, 0x20001e2fc40) in /usr/shlib/libc.so
#9  0x3ff81f300e8 in operator new(0x3ffc01b8400, 0x3ff805ae6d0, 0xfffffffffffffffc, 0x40100000, 0x222384c, 0x20001e2fc40) in /usr/lib/cmplrs/cxx/libcxx.so
#10 0x12000a53c in filter(0x3ffc01b8400, 0x3ff805ae6d0, 0xfffffffffffffffc, 0x40100000, 0x222384c, 0x20001e2fc40) in searchGENOME
#11 0x120008a7c in doSearch(0x3ffc01b8400, 0x3ff805ae6d0, 0xfffffffffffffffc, 0x40100000, 0x222384c, 0x20001e2fc40) in searchGENOME
#12 0x120008e2c in searchThread(0x3ffc01b8400, 0x3ff805ae6d0, 0xfffffffffffffffc, 0x40100000, 0x222384c, 0x20001e2fc40) in searchGENOME
#13 0x3ff805bd2c8 in __thdBase(0x3ffc01b8400, 0x3ff805ae6d0, 0xfffffffffffffffc, 0x40100000, 0x222384c, 0x20001e2fc40) in /usr/shlib/libpthread.so
(ladebug) thread 5  
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      5 <anonymous>               blocked         mut 9       SCHED_OTHER  19

(ladebug) where
>0  0x3ff805caf3c in __hstTransferRegisters(0x3ffc01b8400, 0x3ff805a1ce4, 0x3ff805a3608, 0x20000a0f600, 0x140028000, 0x20000a0f600) in /usr/shlib/libpthre
ad.so
#1  0x3ff805af74c in __osTransferContext(0x3ffc01b8400, 0x3ff805a1ce4, 0x3ff805a3608, 0x20000a0f600, 0x140028000, 0x20000a0f600) in /usr/shlib/libpthread.
so
#2  0x3ff805a3c50 in __dspTransferContext(0x3ffc01b8400, 0x3ff805a1ce4, 0x3ff805a3608, 0x20000a0f600, 0x140028000, 0x20000a0f600) in /usr/shlib/libpthread
.so
#3  0x3ff805a12f4 in __dspDispatch(0x3ffc01b8400, 0x3ff805a1ce4, 0x3ff805a3608, 0x20000a0f600, 0x140028000, 0x20000a0f600) in /usr/shlib/libpthread.so
#4  0x3ff805ab90c in UnknownProcedure32FromFile8(0x3ffc01b8400, 0x3ff805a1ce4, 0x3ff805a3608, 0x20000a0f600, 0x140028000, 0x20000a0f600) in /usr/shlib/lib
pthread.so
#5  0x3ff805ab2f4 in UnknownProcedure31FromFile8(0x3ffc01b8400, 0x3ff805a1ce4, 0x3ff805a3608, 0x20000a0f600, 0x140028000, 0x20000a0f600) in /usr/shlib/lib
pthread.so
#6  0x3ff805abe3c in UnknownProcedure34FromFile8(0x3ffc01b8400, 0x3ff805a1ce4, 0x3ff805a3608, 0x20000a0f600, 0x140028000, 0x20000a0f600) in /usr/shlib/lib
pthread.so
#7  0x3ff801be4f0 in UnknownProcedure12FromFile22(0x3ffc01b8400, 0x3ff805a1ce4, 0x140a275c0, 0x100000, 0x3700498c55, 0x20000a0f600) in /usr/shlib/libc.so
#8  0x3ff800cf2b0 in free(0x3ffc01b8400, 0x3ff805a1ce4, 0x140a275c0, 0x100000, 0x3700498c55, 0x20000a0f600) in /usr/shlib/libc.so
#9  0x3ff81f15a7c in operator delete(0x3ffc01b8400, 0x3ff805a1ce4, 0x140a275c0, 0x100000, 0x3700498c55, 0x20000a0f600) in /usr/lib/cmplrs/cxx/libcxx.so
#10 0x12000b090 in filter(0x3ffc01b8400, 0x3ff805a1ce4, 0x140a275c0, 0x100000, 0x3700498c55, 0x20000a0f600) in searchGENOME
#11 0x120008a7c in doSearch(0x3ffc01b8400, 0x3ff805a1ce4, 0x140a275c0, 0x100000, 0x3700498c55, 0x20000a0f600) in searchGENOME
#12 0x120008dec in searchThread(0x3ffc01b8400, 0x3ff805a1ce4, 0x140a275c0, 0x100000, 0x3700498c55, 0x20000a0f600) in searchGENOME
#13 0x3ff805bd2c8 in __thdBase(0x3ffc01b8400, 0x3ff805a1ce4, 0x140a275c0, 0x100000, 0x3700498c55, 0x20000a0f600) in /usr/shlib/libpthread.so
(ladebug) thread 4
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      4 <anonymous>               blocked         mut 9       SCHED_OTHER  19

(ladebug) where
>0  0x3ff805caf3c in __hstTransferRegisters(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x2000141fc40) in /usr/shlib/libpthread.so
#1  0x3ff805af74c in __osTransferContext(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x2000141fc40) in /usr/shlib/libpthread.so
#2  0x3ff805a3c50 in __dspTransferContext(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x2000141fc40) in /usr/shlib/libpthread.so
#3  0x3ff805a12f4 in __dspDispatch(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x2000141fc40) in /usr/shlib/libpthread.so
#4  0x3ff805ab90c in UnknownProcedure32FromFile8(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x2000141fc40) in /usr/shlib/libpthread.so
#5  0x3ff805ab2f4 in UnknownProcedure31FromFile8(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x2000141fc40) in /usr/shlib/libpthread.so
#6  0x3ff805abe3c in UnknownProcedure34FromFile8(0x3ffc01b8400, 0x3ff805ae6d0, 0x3ff805a33a4, 0x0, 0x0, 0x2000141fc40) in /usr/shlib/libpthread.so
#7  0x3ff801bf6a0 in UnknownProcedure2FromFile22(0x3ffc01b8400, 0x3ff805ae6d0, 0x1, 0x40100000, 0x31c, 0x2000141fc40) in /usr/shlib/libc.so
#8  0x3ff800cdad4 in malloc(0x3ffc01b8400, 0x3ff805ae6d0, 0x1, 0x40100000, 0x31c, 0x2000141fc40) in /usr/shlib/libc.so
#9  0x3ff81f300e8 in operator new(0x3ffc01b8400, 0x3ff805ae6d0, 0x1, 0x40100000, 0x31c, 0x2000141fc40) in /usr/lib/cmplrs/cxx/libcxx.so
#10 0x12000a53c in filter(0x3ffc01b8400, 0x3ff805ae6d0, 0x1, 0x40100000, 0x31c, 0x2000141fc40) in searchGENOME
#11 0x120008a7c in doSearch(0x3ffc01b8400, 0x3ff805ae6d0, 0x1, 0x40100000, 0x31c, 0x2000141fc40) in searchGENOME
#12 0x120008dec in searchThread(0x3ffc01b8400, 0x3ff805ae6d0, 0x1, 0x40100000, 0x31c, 0x2000141fc40) in searchGENOME
#13 0x3ff805bd2c8 in __thdBase(0x3ffc01b8400, 0x3ff805ae6d0, 0x1, 0x40100000, 0x31c, 0x2000141fc40) in /usr/shlib/libpthread.so
(ladebug) thread 3  
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      3 <anonymous>               blocked         mut 9       SCHED_OTHER  19

(ladebug) where
>0  0x3ff805caf3c in __hstTransferRegisters(0x3ffc01b8400, 0x13, 0x3ff805a33a4, 0x0, 0x0, 0x2000283fc40) in /usr/shlib/libpthread.so
#1  0x3ff805af74c in __osTransferContext(0x3ffc01b8400, 0x13, 0x3ff805a33a4, 0x0, 0x0, 0x2000283fc40) in /usr/shlib/libpthread.so
#2  0x3ff805a3c50 in __dspTransferContext(0x3ffc01b8400, 0x13, 0x3ff805a33a4, 0x0, 0x0, 0x2000283fc40) in /usr/shlib/libpthread.so
#3  0x3ff805a12f4 in __dspDispatch(0x3ffc01b8400, 0x13, 0x3ff805a33a4, 0x0, 0x0, 0x2000283fc40) in /usr/shlib/libpthread.so
#4  0x3ff805ab90c in UnknownProcedure32FromFile8(0x3ffc01b8400, 0x13, 0x3ff805a33a4, 0x0, 0x0, 0x2000283fc40) in /usr/shlib/libpthread.so
#5  0x3ff805ab2f4 in UnknownProcedure31FromFile8(0x3ffc01b8400, 0x13, 0x3ff805a33a4, 0x0, 0x0, 0x2000283fc40) in /usr/shlib/libpthread.so
#6  0x3ff805abe3c in UnknownProcedure34FromFile8(0x3ffc01b8400, 0x13, 0x3ff805a33a4, 0x0, 0x0, 0x2000283fc40) in /usr/shlib/libpthread.so
#7  0x3ff801bf6a0 in UnknownProcedure2FromFile22(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/shlib/libc.so
#8  0x3ff800cdad4 in malloc(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/shlib/libc.so
#9  0x3ff81f32050 in UnknownProcedure3FromFile46(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/lib/cmplrs/cxx/libc
xx.so
#10 0x3ff81f34190 in __cxx_v60_dispatch__X4need3new8libcxxso(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/lib/cmp
lrs/cxx/libcxx.so
#11 0x3ff807f29d4 in UnknownProcedure11FromFile0(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/shlib/libexc.so
#12 0x3ff807f2cd8 in exc_dispatch_exception(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/shlib/libexc.so
#13 0x3ff807f39e0 in exc_raise_signal_exception(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/shlib/libexc.so
#14 0x3ff805b9470 in UnknownProcedure8FromFile16(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/shlib/libpthread.so
#15 0x3ff800d0b9c in __sigtramp(0x3ffc01b8400, 0x13, 0x20001927600, 0x40100000, 0x3ff81f34150, 0x2000283fc40) in /usr/shlib/libc.so
#16 0x3ff801be2c0 in UnknownProcedure12FromFile22(0x3ffc0086f90, 0x0, 0x100000, 0x0, 0x3f003d2689, 0xa8003d2689) in /usr/shlib/libc.so
#17 0x3ff800cf2b0 in free(0x3ffc0086f90, 0x0, 0x100000, 0x0, 0x3f003d2689, 0xa8003d2689) in /usr/shlib/libc.so
#18 0x3ff81f15a7c in operator delete(0x3ffc0086f90, 0x0, 0x100000, 0x0, 0x3f003d2689, 0xa8003d2689) in /usr/lib/cmplrs/cxx/libcxx.so
#19 0x12000b090 in filter(0x3ffc0086f90, 0x0, 0x100000, 0x0, 0x3f003d2689, 0xa8003d2689) in searchGENOME
#20 0x120008a7c in doSearch(0x3ffc0086f90, 0x0, 0x100000, 0x0, 0x3f003d2689, 0xa8003d2689) in searchGENOME
#21 0x120008dec in searchThread(0x3ffc0086f90, 0x0, 0x100000, 0x0, 0x3f003d2689, 0xa8003d2689) in searchGENOME
#22 0x3ff805bd2c8 in __thdBase(0x3ffc0086f90, 0x0, 0x100000, 0x0, 0x3f003d2689, 0xa8003d2689) in /usr/shlib/libpthread.so
(ladebug) thread 2
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      2 <anonymous>               blocked         kern usleep SCHED_OTHER  19

(ladebug) where
>0  0x3ff800e5c38 in __usleep_thread(0x20000f15a40, 0x0, 0x0, 0x0, 0x0, 0x356) in /usr/shlib/libc.so
#1  0x3ff801b3314 in __usleep(0x20000f15a40, 0x0, 0x0, 0x0, 0x0, 0x356) in /usr/shlib/libc.so
#2  0x1200091ac in loaderThread(0x20000f15a40, 0x0, 0x0, 0x0, 0x0, 0x356) in searchGENOME
#3  0x3ff805bd2c8 in __thdBase(0x20000f15a40, 0x0, 0x0, 0x0, 0x0, 0x356) in /usr/shlib/libpthread.so
(ladebug) thread 1
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>*     1 default thread            blocked         kern usleep SCHED_OTHER  19

(ladebug) where
>0  0x3ff800e5c38 in __usleep_thread(0x11fffbe88, 0x0, 0x5254, 0x0, 0x0, 0x140002408) in /usr/shlib/libc.so
#1  0x3ff801b3314 in __usleep(0x11fffbe88, 0x0, 0x5254, 0x0, 0x0, 0x140002408) in /usr/shlib/libc.so
#2  0x12000661c in main(0x11fffbe88, 0x0, 0x5254, 0x0, 0x0, 0x140002408) in searchGENOME
#3  0x1200055c8 in __start(0x11fffbe88, 0x0, 0x5254, 0x0, 0x0, 0x140002408) in searchGENOME

//
//
//  Stack trace #2
//
//

(ladebug) show thread
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>*     4 <anonymous>               blocked         kern usleep SCHED_OTHER  19
       1 default thread            blocked         mut 15      SCHED_OTHER  19
      -1 manager thread            blk SCS                     SCHED_RR     19
      -2 null thread for VP 2      running VP 2                null thread  -1
       2 <anonymous>               blocked         mut 15      SCHED_OTHER  19
      -3 null thread for VP 3      running VP 3                null thread  -1
       3 <anonymous>               blocked         mut 15      SCHED_OTHER  19
      -4 null thread for VP 4      running VP 4                null thread  -1
       5 <anonymous>               blocked         mut 15      SCHED_OTHER  19
      -5 null thread for VP 5      running VP 5                null thread  -1
       6 <anonymous>               blocked         mut 15      SCHED_OTHER  19
       7 <anonymous>               blocked         mut 15      SCHED_OTHER  19
       8 <anonymous>               blocked         mut 15      SCHED_OTHER  19



(ladebug) show mutex
Mutex  Name                      State Owner  Pri Type     Waiters (+Count)
------ ------------------------- ----- ------ --- -------- --------------------
     1 Once                                       Normal
     2 debugger client registry                   Normal
     3 VM stats                                   Normal
     4 key creation                               Normal
     5 malloc heap                                Normal
     6 malloc hash                                Normal
     7 malloc cache[0]                            Normal
     8 malloc cache[1]                            Normal
     9 malloc cache[2]                            Normal
    10 malloc cache[3]                            Normal
    11 malloc cache[4]                            Normal
    12 malloc cache[5]                            Normal
    13 malloc cache[6]                            Normal
    14 malloc cache[7]                            Normal
    15 malloc cache[8]           Lock             Normal   6, 7, 8, 1, 5, 2, 3
    16 malloc cache[9]                            Normal
    17 malloc cache[10]                           Normal
    18 malloc cache[11]                           Normal
    19 malloc cache[12]                           Normal
    20 malloc cache[13]                           Normal
    21 malloc cache[14]                           Normal
    22 malloc cache[15]                           Normal
    23 malloc cache[16]                           Normal
    24 malloc cache[17]                           Normal
    25 malloc cache[18]                           Normal
    26 malloc cache[19]                           Normal
    27 malloc cache[20]                           Normal
    28 malloc cache[21]                           Normal
    29 malloc cache[22]                           Normal
    30 malloc cache[23]                           Normal
    31 malloc cache[24]                           Normal
    32 malloc cache[25]                           Normal
    33 malloc cache[26]                           Normal
    34 malloc cache[27]                           Normal
    35 malloc cache[28]                           Normal
    36 brk                                        Normal
    37 exc cr                                     Recurs
    38 exc read rwl                               Normal
    39 VM 0 lookaside                             Normal
    40 VM 1 lookaside                             Normal
    41 VM 2 lookaside                             Normal
    42 VM 3 lookaside                             Normal
    43 VM 4 lookaside                             Normal
    44 VM 5 lookaside                             Normal
    45 VM 6 lookaside                             Normal
    46 VM 0 cache                                 Normal
    47 VM 1 cache                                 Normal
    48 VM 2 cache                                 Normal
    49 Global lock                                Recurs
    50 ldr                                        Recurs
    51 <anonymous>                                Recurs
    52 stderr                                     Recurs
    53 stdout                                     Recurs
    54 <anonymous>                                Recurs
    55 <anonymous>                                Recurs
    56 inputTailMutex(0x14000105                  Normal
    57 queryMatchMutex(0x1400010                  Normal




(ladebug) thread 1
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      1 default thread            blocked         mut 15      SCHED_OTHER  19

(ladebug) where
>0  0x3ff805ba8ac in __hstTransferRegisters(0x20002d47600, 0x0, 0x0, 0x100000000, 0x20002d47c40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#1  0x3ff805acf74 in __osTransferContext(0x20002d47600, 0x0, 0x0, 0x100000000, 0x20002d47c40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#2  0x3ff805a004c in __dspDispatch(0x20002d47600, 0x0, 0x0, 0x100000000, 0x20002d47c40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#3  0x3ff805a94e4 in UnknownProcedure146FromFile0(0x20002d47600, 0x0, 0x0, 0x100000000, 0x20002d47c40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#4  0x3ff805a9bf8 in UnknownProcedure148FromFile0(0x20002d47600, 0x0, 0x0, 0x100000000, 0x20002d47c40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#5  0x3ff801bed30 in UnknownProcedure12FromFile22(0x20002d47600, 0x0, 0x0, 0x0, 0x0, 0x3ff805ac400) in /usr/shlib/libc.so
#6  0x3ff800cf2c0 in free(0x20002d47600, 0x0, 0x0, 0x0, 0x0, 0x3ff805ac400) in /usr/shlib/libc.so
#7  0x3ff81f15a7c in operator delete(0x20002d47600, 0x0, 0x0, 0x0, 0x0, 0x3ff805ac400) in /usr/lib/cmplrs/cxx/libcxx.so
#8  0x3ff81f2f53c in operator delete[](0x20002d47600, 0x0, 0x0, 0x0, 0x0, 0x3ff805ac400) in /usr/lib/cmplrs/cxx/libcxx.so
#9  0x1200073fc in main(0x20002d47600, 0x0, 0x0, 0x0, 0x0, 0x3ff805ac400) in searchGENOME
#10 0x120006088 in __start(0x20002d47600, 0x0, 0x0, 0x0, 0x0, 0x3ff805ac400) in searchGENOME

(ladebug) thread 2
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      2 <anonymous>               blocked         mut 15      SCHED_OTHER  19

(ladebug) where
>0  0x3ff805ba8ac in __hstTransferRegisters(0x2000141f600, 0x0, 0x0, 0x100000000, 0x3ff805a0194, 0x3ff805ac400) in /usr/shlib/libpthread.so
#1  0x3ff805acf74 in __osTransferContext(0x2000141f600, 0x0, 0x0, 0x100000000, 0x3ff805a0194, 0x3ff805ac400) in /usr/shlib/libpthread.so
#2  0x3ff805a004c in __dspDispatch(0x2000141f600, 0x0, 0x0, 0x100000000, 0x3ff805a0194, 0x3ff805ac400) in /usr/shlib/libpthread.so
#3  0x3ff805a94e4 in UnknownProcedure146FromFile0(0x2000141f600, 0x0, 0x0, 0x100000000, 0x3ff805a0194, 0x3ff805ac400) in /usr/shlib/libpthread.so
#4  0x3ff805a9bf8 in UnknownProcedure148FromFile0(0x2000141f600, 0x0, 0x0, 0x100000000, 0x3ff805a0194, 0x3ff805ac400) in /usr/shlib/libpthread.so
#5  0x3ff801bfee0 in UnknownProcedure2FromFile22(0x2000141f600, 0x0, 0x4, 0x0, 0x0, 0x3ff805ac400) in /usr/shlib/libc.so
#6  0x3ff800cdae4 in malloc(0x2000141f600, 0x0, 0x4, 0x0, 0x0, 0x3ff805ac400) in /usr/shlib/libc.so
#7  0x3ff81f300e8 in operator new(0x2000141f600, 0x0, 0x4, 0x0, 0x0, 0x3ff805ac400) in /usr/lib/cmplrs/cxx/libcxx.so
#8  0x3ff81f2f5dc in operator new[](0x2000141f600, 0x0, 0x4, 0x0, 0x0, 0x3ff805ac400) in /usr/lib/cmplrs/cxx/libcxx.so
#9  0x12000a9bc in deadlockDetector(0x2000141f600, 0x0, 0x4, 0x0, 0x0, 0x3ff805ac400) in searchGENOME
#10 0x3ff805c67e0 in __thdBase(0x2000141f600, 0x0, 0x4, 0x0, 0x0, 0x3ff805ac400) in /usr/shlib/libpthread.so

(ladebug) thread 3  
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      3 <anonymous>               blocked         mut 15      SCHED_OTHER  19

(ladebug) where
>0  0x3ff805ba8ac in __hstTransferRegisters(0x20001e2f600, 0x0, 0x0, 0x100000000, 0x0, 0x3ff805ac400) in /usr/shlib/libpthread.so
#1  0x3ff805acf74 in __osTransferContext(0x20001e2f600, 0x0, 0x0, 0x100000000, 0x0, 0x3ff805ac400) in /usr/shlib/libpthread.so
#2  0x3ff805a004c in __dspDispatch(0x20001e2f600, 0x0, 0x0, 0x100000000, 0x0, 0x3ff805ac400) in /usr/shlib/libpthread.so
#3  0x3ff805a94e4 in UnknownProcedure146FromFile0(0x20001e2f600, 0x0, 0x0, 0x100000000, 0x0, 0x3ff805ac400) in /usr/shlib/libpthread.so
#4  0x3ff805a9bf8 in UnknownProcedure148FromFile0(0x20001e2f600, 0x0, 0x0, 0x100000000, 0x0, 0x3ff805ac400) in /usr/shlib/libpthread.so
#5  0x3ff801bfee0 in UnknownProcedure2FromFile22(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libc.so
#6  0x3ff800cdae4 in malloc(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libc.so
#7  0x3ff805be20c in UnknownProcedure0FromFile99(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libpthread.so
#8  0x3ff805be508 in UnknownProcedure1FromFile99(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libpthread.so
#9  0x3ff805be5d0 in UnknownProcedure3FromFile99(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libpthread.so
#10 0x3ff807f369c in UnknownProcedure15FromFile0(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libexc.so
#11 0x3ff807f3a08 in exc_raise_signal_exception(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libexc.so
#12 0x3ff805b5a9c in UnknownProcedure283FromFile0(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libpthread.so
#13 0x3ff800d0bbc in __sigtramp(0x20001e2f600, 0x0, 0x0, 0x0, 0x3ff80119094, 0x3ff805ac400) in /usr/shlib/libc.so
#14 0x3ff800e2158 in __kill(0x27aaf6, 0x6, 0x1000000, 0x0, 0x0, 0x1) in /usr/shlib/libc.so
#15 0x12000aad0 in deadlockChecker(0x27aaf6, 0x6, 0x1000000, 0x0, 0x0, 0x1) in searchGENOME
#16 0x3ff805c67e0 in __thdBase(0x27aaf6, 0x6, 0x1000000, 0x0, 0x0, 0x1) in /usr/shlib/libpthread.so

(ladebug) thread 4
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>*     4 <anonymous>               blocked         kern usleep SCHED_OTHER  19

(ladebug) where
>0  0x3ff800e5e68 in __usleep_thread(0x20002335a30, 0x20002335a28, 0x0, 0x0, 0x12000a690, 0x12000a6c0) in /usr/shlib/libc.so
#1  0x3ff80ba527c in nanosleep(0x20002335a30, 0x20002335a28, 0x0, 0x0, 0x12000a690, 0x12000a6c0) in /usr/shlib/librt.so
#2  0x12000a6ec in loaderThread(0x20002335a30, 0x20002335a28, 0x0, 0x0, 0x12000a690, 0x12000a6c0) in searchGENOME
#3  0x3ff805c67e0 in __thdBase(0x20002335a30, 0x20002335a28, 0x0, 0x0, 0x12000a690, 0x12000a6c0) in /usr/shlib/libpthread.so

(ladebug) thread 5
  Thread Name                      State           Substate    Policy       Pri
  ------ ------------------------- --------------- ----------- ------------ ---
>      5 <anonymous>               blocked         mut 15      SCHED_OTHER  19

(ladebug) where
>0  0x3ff805ba8ac in __hstTransferRegisters(0x20000a0f600, 0x0, 0x0, 0x100000000, 0x20000a0fc40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#1  0x3ff805acf74 in __osTransferContext(0x20000a0f600, 0x0, 0x0, 0x100000000, 0x20000a0fc40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#2  0x3ff805a004c in __dspDispatch(0x20000a0f600, 0x0, 0x0, 0x100000000, 0x20000a0fc40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#3  0x3ff805a94e4 in UnknownProcedure146FromFile0(0x20000a0f600, 0x0, 0x0, 0x100000000, 0x20000a0fc40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#4  0x3ff805a9bf8 in UnknownProcedure148FromFile0(0x20000a0f600, 0x0, 0x0, 0x100000000, 0x20000a0fc40, 0x3ff805ac400) in /usr/shlib/libpthread.so
#5  0x3ff801bfee0 in UnknownProcedure2FromFile22(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/shlib/libc.so
#6  0x3ff800cdae4 in malloc(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/shlib/libc.so
#7  0x3ff81f32050 in UnknownProcedure3FromFile46(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/lib/cmplrs/cxx/libcxx.so
#8  0x3ff81f34190 in __cxx_v60_dispatch__X4need3new8libcxxso(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/lib/cmplrs/cxx/libcxx.so
#9  0x3ff807f29d4 in UnknownProcedure11FromFile0(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/shlib/libexc.so
#10 0x3ff807f2cd8 in exc_dispatch_exception(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/shlib/libexc.so
#11 0x3ff807f39e0 in exc_raise_signal_exception(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/shlib/libexc.so
#12 0x3ff805b5a9c in UnknownProcedure283FromFile0(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/shlib/libpthread.so
#13 0x3ff800d0bbc in __sigtramp(0x20000a0f600, 0x0, 0x2000283f600, 0x2000283d288, 0x3ff81f34150, 0x3ff805ac400) in /usr/shlib/libc.so
#14 0x3ff801beb00 in UnknownProcedure12FromFile22(0x3ffc0087290, 0x0, 0x100000, 0x0, 0x3ff801bead0, 0x3ff801beba4) in /usr/shlib/libc.so
#15 0x3ff800cf2c0 in free(0x3ffc0087290, 0x0, 0x100000, 0x0, 0x3ff801bead0, 0x3ff801beba4) in /usr/shlib/libc.so
#16 0x3ff81f15a7c in operator delete(0x3ffc0087290, 0x0, 0x100000, 0x0, 0x3ff801bead0, 0x3ff801beba4) in /usr/lib/cmplrs/cxx/libcxx.so
#17 0x120009680 in ~encodedQuery(0x3ffc0087290, 0x0, 0x100000, 0x0, 0x3ff801bead0, 0x3ff801beba4) in searchGENOME
#18 0x120009f9c in doSearch(0x3ffc0087290, 0x0, 0x100000, 0x0, 0x3ff801bead0, 0x3ff801beba4) in searchGENOME
#19 0x12000a39c in searchThread(0x3ffc0087290, 0x0, 0x100000, 0x0, 0x3ff801bead0, 0x3ff801beba4) in searchGENOME
#20 0x3ff805c67e0 in __thdBase(0x3ffc0087290, 0x0, 0x100000, 0x0, 0x3ff801bead0, 0x3ff801beba4) in /usr/shlib/libpthread.so

#endif  //  DONT_EVER_ENABLE_THIS
