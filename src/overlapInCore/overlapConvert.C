


#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"

int
main(int argc, char **argv) {
  bool      asCoords = true;
  bool      asHangs  = true;

  char     *gkpStoreName = NULL;
  gkStore  *gkpStore = NULL;

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-coords") == 0) {
      asCoords = true;
      asHangs  = false;

    } else if (strcmp(argv[arg], "-hangs") == 0) {
      asCoords = false;
      asHangs  = true;

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (AS_UTL_fileExists(argv[arg])) {
      break;

    } else {
      err++;
    }

    arg++;
  }

  if ((gkpStoreName == NULL) && (asCoords == true))
    err++;

  if ((err) || (arg == argc)) {
    fprintf(stderr, "usage: %s [options] file.ovb[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G             gkpStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -coords        output coordiantes on reads\n");
    fprintf(stderr, "  -hangs         output hangs on reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    if ((gkpStoreName == NULL) && (asCoords == true))
      fprintf(stderr, "ERROR:  -coords mode requires a gkpStore (-G)\n");

    exit(1);
  }

  if (gkpStoreName)
    gkpStore = new gkStore(gkpStoreName);

  while (arg < argc) {
    ovFile      *of = new ovFile(argv[arg], ovFileFull);
    ovsOverlap   ov;

    while (of->readOverlap(&ov)) {
      if (asCoords) {
        fprintf(stdout, "%10"F_U32P" %10"F_U32P"  %c  %6"F_U32P" %6"F_U32P"  %6"F_U32P" %6"F_U32P"  %6.3f\n",
                ov.a_iid, ov.b_iid,
                ov.flipped() ? 'I' : 'N',
                ov.a_bgn(), ov.a_end(gkpStore),
                ov.b_bgn(), ov.b_end(gkpStore),
                ov.erate());
      }
      
      if (asHangs) {
        fprintf(stdout, "%10"F_U32P" %10"F_U32P"  %c  %6"F_S32P" %6"F_S32P"  %6.3f%s\n",
                ov.a_iid, ov.b_iid,
                ov.flipped() ? 'I' : 'N',
                ov.a_hang(), ov.b_hang(),
                ov.erate(),
                (ov.overlapIsDovetail()) ? "" : "  PARTIAL");
      }
    }

    delete of;

    arg++;
  }


  exit(0);
}
