

const char *mainid = "$Id: tgStoreDump.C 6665 2015-01-13 16:48:09Z bri $";

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"


int
main(int argc, char **argv) {
  tgTig  tig;

  gkStore  *gkpStore = new gkStore(argv[1]);

  FILE *F = fopen(argv[2], "r");

  tig.loadFromStream(F);

  fclose(F);

  uint32  displayWidth    = 250;
  uint32  displaySpacing  = 10;
  bool    withQV          = false;
  bool    withDots        = true;

  tig.display(stdout, gkpStore, displayWidth, displaySpacing, withQV, withDots);

  exit(0);
}

