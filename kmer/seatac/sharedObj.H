#ifndef SHAREDOBJ_H
#define SHAREDOBJ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>
#include <errno.h>

//  A wrapper around a shared object.
//
//  Responsible for opening, accessing and closing a shared object.
//

class sharedObj {
public:

  sharedObj(char *p) {
    path = new char [strlen(p) + 1];
    strcpy(path, p);

    handle = dlopen(path, RTLD_NOW);
    if (handle == 0L) {
      fprintf(stderr, "ERROR: Failed to open shared object '%s'\n%s\n", path, dlerror());
      exit(1);
    }
  };

  ~sharedObj() {
    if (dlclose(handle)) {
      fprintf(stderr, "WARNING: Failed to close shared object '%s'\n%s\n", path, dlerror());
    }
    delete [] path;
  };

  void   *get(char *name) {
    void *ptr = 0L;

    errno = 0;
    ptr = dlsym(handle, name);
    if (errno) {
      fprintf(stderr, "ERROR: Failed to find symbol '%s' in shared object '%s'\n%s\n", name, path, dlerror());
      exit(1);
    }
    if (ptr == 0L) {
      fprintf(stderr, "ERROR: Symbol '%s' not present in shared object '%s'\n", name, path);
      exit(1);
    }

    return(ptr);
  };

private:
  char *path;
  void *handle;
};


#endif  //  SHAREDOBJ_H
