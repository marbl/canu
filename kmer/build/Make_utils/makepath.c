#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

/** assumes both inputs, if they exist, are directories
    and returns a relative path to the first from the second **/

int main(int argc,char **argv) {
  char path[PATH_MAX+1],subpath[PATH_MAX+1];
  char *p=path,*sp=subpath;
  *path = *subpath = 0;

  if(argc<2) {
    fprintf(stderr,"usage: %s path [subpath]\n",argv[0]);
    exit(1);
  }
  if ( !realpath(argv[1],path) ) {
    fprintf(stderr,"%s: the path '%s' does not exist.\n",argv[0],argv[1]);
    fprintf(stderr,"%s: perhaps you are building from an incomplete CVS?\n",argv[0]);
    printf("\n");
    exit(0);
  }
  /**fprintf(stderr,"   path:%s\n",   path);**/

  /** if subpath is found then it is non-empty otherwise empty **/
  if ( argc==3 && realpath(argv[2],subpath) ) {
    /**fprintf(stderr,"subpath:%s\n",subpath);**/

    for(; *p && *sp && *p==*sp; ++p, ++sp);

    if (*sp==0 && *p=='/') ++p;
    if (*sp=='/' && *p==0) ++sp;

    if(*sp) {
      printf("../") ;
      for( ; *sp ; ++sp)  if(*sp=='/') printf("../") ;
    }
  }

  if (*p)
    printf("%s/\n",p);
  else
    printf("\n");

  exit(0);
}
