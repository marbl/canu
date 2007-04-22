#include "snapper2.H"

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
 checkAgain:

  //  Wait for the tester to test
  //
  while (deadlockTested == 0)
    sleep(5);

  //  Give it another chance
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

  //  Reset the testing/checking flags
  //
  deadlockPassed = 0;
  deadlockTested = 0;

  goto checkAgain;

  return(0L);  //  Ignore the warning!
}

#endif
