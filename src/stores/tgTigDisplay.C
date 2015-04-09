

const char *mainid = "$Id: tgStoreDump.C 6665 2015-01-13 16:48:09Z bri $";

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"


int
main(int argc, char **argv) {
  tgTig  tig;
  char  *gkpName;
  char  *tigFileName;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigFileName = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if (gkpName == NULL)
    err++;
  if (tigFileName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -t tigFile\n", argv[0]);
    exit(1);
  }

  gkStore  *gkpStore = new gkStore(gkpName);

  FILE *F = fopen(tigFileName, "r");

  tig.loadFromStream(F);

  fclose(F);

  uint32  displayWidth    = 250;
  uint32  displaySpacing  = 10;
  bool    withQV          = false;
  bool    withDots        = true;

  tig.display(stdout, gkpStore, displayWidth, displaySpacing, withQV, withDots);

  exit(0);
}
