#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

#define   TEST_INTERVAL_LIST
#define   TEST_SIZE  2000
#define   TEST_ITERS 1000

#include "libbri.H"
#include "intervalList.H"

void
fixedTest(void) {
  intervalList   G(10);

  G.addInterval(110);  fprintf(stderr, "Adding %3d -> %3d:\t", 110, 110+10);  G.dump();
  G.addInterval(130);  fprintf(stderr, "Adding %3d -> %3d:\t", 130, 130+10);  G.dump();
  G.addInterval(105);  fprintf(stderr, "Adding %3d -> %3d:\t", 105, 105+10);  G.dump();
  G.addInterval(115);  fprintf(stderr, "Adding %3d -> %3d:\t", 115, 115+10);  G.dump();
  G.addInterval(124);  fprintf(stderr, "Adding %3d -> %3d:\t", 124, 124+10);  G.dump();
  G.addInterval( 50);  fprintf(stderr, "Adding %3d -> %3d:\t",  50,  50+10);  G.dump();
  G.addInterval(200);  fprintf(stderr, "Adding %3d -> %3d:\t", 200, 200+10);  G.dump();
  G.addInterval(150);  fprintf(stderr, "Adding %3d -> %3d:\t", 150, 150+10);  G.dump();
  G.addInterval(205);  fprintf(stderr, "Adding %3d -> %3d:\t", 205, 205+10);  G.dump();
  G.addInterval(195);  fprintf(stderr, "Adding %3d -> %3d:\t", 195, 195+10);  G.dump();
  G.addInterval( 61);  fprintf(stderr, "Adding %3d -> %3d:\t",  61,  61+10);  G.dump();
  G.addInterval( 72);  fprintf(stderr, "Adding %3d -> %3d:\t",  72,  72+10);  G.dump();
  G.addInterval( 65);  fprintf(stderr, "Adding %3d -> %3d:\t",  65,  65+10);  G.dump();
  G.addInterval( 83);  fprintf(stderr, "Adding %3d -> %3d:\t",  83,  83+10);  G.dump();
  G.addInterval( 94);  fprintf(stderr, "Adding %3d -> %3d:\t",  94,  94+10);  G.dump();
  G.addInterval( 75);  fprintf(stderr, "Adding %3d -> %3d:\t",  75,  75+10);  G.dump();
  G.addInterval( 84);  fprintf(stderr, "Adding %3d -> %3d:\t",  84,  84+10);  G.dump();
  G.addInterval(104);  fprintf(stderr, "Adding %3d -> %3d:\t", 104, 104+10);  G.dump();
  G.addInterval(114);  fprintf(stderr, "Adding %3d -> %3d:\t", 114, 114+10);  G.dump();
  G.addInterval(124);  fprintf(stderr, "Adding %3d -> %3d:\t", 124, 124+10);  G.dump();
  G.addInterval(134);  fprintf(stderr, "Adding %3d -> %3d:\t", 134, 134+10);  G.dump();
  G.addInterval(144);  fprintf(stderr, "Adding %3d -> %3d:\t", 144, 144+10);  G.dump();
  G.addInterval( 51);  fprintf(stderr, "Adding %3d -> %3d:\t",  51,  51+10);  G.dump();
  G.addInterval(161);  fprintf(stderr, "Adding %3d -> %3d:\t", 161, 161+10);  G.dump();
  G.addInterval(172);  fprintf(stderr, "Adding %3d -> %3d:\t", 172, 172+10);  G.dump();
  G.addInterval(183);  fprintf(stderr, "Adding %3d -> %3d:\t", 183, 183+10);  G.dump();
  G.addInterval(156);  fprintf(stderr, "Adding %3d -> %3d:\t", 156, 156+10);  G.dump();
  G.addInterval(166);  fprintf(stderr, "Adding %3d -> %3d:\t", 166, 166+10);  G.dump();
  G.addInterval(176);  fprintf(stderr, "Adding %3d -> %3d:\t", 176, 176+10);  G.dump();
  G.addInterval(186);  fprintf(stderr, "Adding %3d -> %3d:\t", 186, 186+10);  G.dump();
  G.addInterval(  0);  fprintf(stderr, "Adding %3d -> %3d:\t",   0,   0+10);  G.dump();
  G.addInterval(  0);  fprintf(stderr, "Adding %3d -> %3d:\t",   0,   0+10);  G.dump();
  G.addInterval(  1);  fprintf(stderr, "Adding %3d -> %3d:\t",   1,   1+10);  G.dump();
  G.addInterval(  2);  fprintf(stderr, "Adding %3d -> %3d:\t",   2,   2+10);  G.dump();
  G.addInterval(300);  fprintf(stderr, "Adding %3d -> %3d:\t", 300, 300+10);  G.dump();
  G.addInterval(320);  fprintf(stderr, "Adding %3d -> %3d:\t", 320, 320+10);  G.dump();
  G.addInterval(280);  fprintf(stderr, "Adding %3d -> %3d:\t", 280, 280+10);  G.dump();
  G.addInterval( 20);  fprintf(stderr, "Adding %3d -> %3d:\t",  20,  20+10);  G.dump();
}


void
main(int argc, char **argv) {

  fixedTest();

  srand48(237831);

loop:


#if 0
  intervalList  *G = new intervalList(10);
  for (u32bit i=0; i<TEST_ITERS; i++) {
    G->addInterval(floor(drand48() * (TEST_SIZE - 10)));
    G->test();
  }
  G->dump();
  delete G;
#endif

  intervalList *A = new intervalList(10);
  intervalList *B = new intervalList(10);
  intervalList *C = new intervalList(10);

  for (u32bit i=0; i<TEST_ITERS; i++) {
    u32bit j = floor(drand48() * (TEST_SIZE - 10));

    C->addInterval(j);
    if (drand48() < 0.5)
      A->addInterval(j);
    else
      B->addInterval(j);
  }

  fprintf(stderr, "A & B ----------------------------------------\n");
  A->dump();
  B->dump();

  A->merge(B);

  fprintf(stderr, "A & C ----------------------------------------\n");
  A->dump();
  C->dump();

  A->compare(C);

  delete A;
  delete B;
  delete C;

  goto loop;
}
