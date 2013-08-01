#include "AS_global.H"
#include "AS_PER_gkpStore.H"

int
main(void) {
  fprintf(stderr, "gkPackedFragment  %d\n", sizeof(gkPackedFragment));
  fprintf(stderr, "gkNormalFragment  %d\n", sizeof(gkNormalFragment));
  fprintf(stderr, "gkStrobeFragment  %d\n", sizeof(gkStrobeFragment));
  fprintf(stderr, "gkFragment        %d\n", sizeof(gkFragment));
}
