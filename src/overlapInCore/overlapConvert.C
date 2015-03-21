


#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"

#include <vector>

using namespace std;


int
main(int argc, char **argv) {
  bool            asCoords = true;

  char           *gkpStoreName = NULL;
  gkStore        *gkpStore = NULL;

  vector<char *>  files;


  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-coords") == 0) {
      asCoords = true;

    } else if (strcmp(argv[arg], "-hangs") == 0) {
      asCoords = false;

    } else if (AS_UTL_fileExists(argv[arg])) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((gkpStoreName == NULL) && (asCoords == true))
    err++;

  if ((err) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] file.ovb[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G             gkpStore (needed for -coords, the default)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -coords        output coordiantes on reads\n");
    fprintf(stderr, "  -hangs         output hangs on reads\n");

    if ((gkpStoreName == NULL) && (asCoords == true))
      fprintf(stderr, "ERROR:  -coords mode requires a gkpStore (-G)\n");

    if (files.size() == 0)
      fprintf(stderr, "ERROR:  no overlap files supplied\n");

    exit(1);
  }

  if (gkpStoreName)
    gkpStore = new gkStore(gkpStoreName);

  char  *ovStr = new char [1024];

  for (uint32 ff=0; ff<files.size(); ff++) {
    ovFile      *of = new ovFile(files[ff], ovFileFull);
    ovsOverlap   ov;

    while (of->readOverlap(&ov))
      fputs(ov.toString(ovStr, gkpStore, asCoords), stdout);

    delete of;

    arg++;
  }

  delete [] ovStr;

  exit(0);
}
