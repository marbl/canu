  //  If we want to use palloc() in here, we need to have a private
  //  pallochandle, which means that all memory allocation must be
  //  done within the class.
  //

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>

#include "util++.H"


#ifdef WITH_PALLOC

  void  *_pallochandle;

  void  *ckalloc(size_t amount) {
    return(palloc2(amount, _pallochandle));
  };
  void  *ckcalloc(size_t amount) {
    void *x = palloc2(amount, _pallochandle);
    memset(x, 0x00, amount);
    return(x);
  };
  char *strsave(const char *s) {
    char *p = (char *)ckalloc(strlen(s)+1);       /* +1 to hold '\0' */
    return strcpy(p, s);
  }
  void   ckfree(void *) {
  };
  void   mem_init(void) {
    _pallochandle = pallochandle(2 * 1024 * 1024);
  }
  void   mem_end(void) {
    pfree2(_pallochandle);
    pfreehandle(_pallochandle);
  }

#else

  void  *ckalloc(size_t amount) {
#ifdef __APPLE__
    if (amount == 0)
      amount = 16;
#else
    if (amount == 0)
      amount = 8;
#endif

    void *p = malloc(amount);
    if (p == NULL) {
      fprintf(stderr, "Can't allocate "u64bitFMT" bytes.\n", (u64bit)amount);
      kill(getpid(), SIGKILL);
    }

    return(p);
  };
  // strsave - save string s somewhere; return address
  char *strsave(const char *s) {
    char *p = (char *)ckalloc(strlen(s)+1);       /* +1 to hold '\0' */
    return strcpy(p, s);
  }
  void   ckfree(void *P) {
    free(P);
  };
  void   mem_init(void) { ; }
  void   mem_end(void) { ; }

#endif
