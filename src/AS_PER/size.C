#include "AS_global.h"
#include "AS_PER_gkpStore.h"

int
main(void) {
  fprintf(stderr, "gkShortFragment   %d\n", sizeof(gkShortFragment));
  fprintf(stderr, "gkMediumFragment  %d\n", sizeof(gkMediumFragment));
  fprintf(stderr, "gkLongFragment    %d\n", sizeof(gkLongFragment));
  fprintf(stderr, "gkFragment        %d\n", sizeof(gkFragment));
}
