#include <stdio.h>
#include <stdlib.h>
#include "../sim4reader.h"


int
main(int argc, char **argv) {
  sim4polish *p;

  while ((p = readPolish(stdin)) != 0L) {
    if (p->comment) {
      fprintf(stdout, "================================================================================\n");
      fprintf(stdout, "p->comment='%s'\n", p->comment);
    }
#if 0
    if ((p->numExons == 11) && (p->percentIdentity < 97)) {
      fprintf(stdout, "================================================================================\n");
      //printPolish(stdout, p);

#if 0
      while (p->numExons > 1) {
        s4p_deleteExon(p, 0);
        fprintf(stdout, "========== %d\n", p->numExons);
        printPolish(stdout, p);
      }
#endif

#if 0
      while (p->numExons > 1) {
        s4p_deleteExon(p, p->numExons-1);
        fprintf(stdout, "========== %d\n", p->numExons);
        printPolish(stdout, p);
      }
#endif

#if 1
      while (p->numExons > 3) {
        s4p_deleteExon(p, 2);
        fprintf(stdout, "========== %d\n", p->numExons);
        printPolish(stdout, p);
      }
#endif

    }
#endif

    destroyPolish(p);
  }
}

