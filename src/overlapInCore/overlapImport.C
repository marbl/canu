


#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"

#include "splitToWords.H"

#include <vector>

using namespace std;

#define  TYPE_NONE    'N'
#define  TYPE_LEGACY  'L'
#define  TYPE_COORDS  'C'
#define  TYPE_HANGS   'H'
#define  TYPE_RAW     'R'



int
main(int argc, char **argv) {
  char                  *gkpStoreName = NULL;
  gkStore               *gkpStore = NULL;

  char                  *ovlFileName = NULL;
  char                  *ovlStoreName = NULL;

  char                   inType = TYPE_NONE;

  vector<char *>         files;


  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      ovlFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-legacy") == 0) {
      inType = TYPE_LEGACY;

    } else if (strcmp(argv[arg], "-coords") == 0) {
      fprintf(stderr, "-coords not implemented.\n"), exit(1);
      inType = TYPE_COORDS;

    } else if (strcmp(argv[arg], "-hangs") == 0) {
      fprintf(stderr, "-hangs not implemented.\n"), exit(1);
      inType = TYPE_HANGS;

    } else if (strcmp(argv[arg], "-raw") == 0) {
      fprintf(stderr, "-raw not implemented.\n"), exit(1);
      inType = TYPE_RAW;

    } else if ((strcmp(argv[arg], "-") == 0) ||
               (AS_UTL_fileExists(argv[arg]))) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (gkpStoreName == NULL)
    err++;
  if (inType == TYPE_NONE)
    err++;

  if ((err) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] ascii-ovl-file-input.[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -G name.gkpStore   path to valid gatekeeper store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "  -o file.ovb        output file name\n");
    fprintf(stderr, "  -O name.ovlStore   output overlap store");
    fprintf(stderr, "\n");
    fprintf(stderr, "Format options:\n");
    fprintf(stderr, "  -legacy            'CA8 overlapStore -d' format\n");
    fprintf(stderr, "  -coords            'overlapConvert -coords' format\n");
    fprintf(stderr, "  -hangs             'overlapConvert -hangs' format\n");
    fprintf(stderr, "  -raw               'overlapConvert -raw' format\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input file can be stdin ('-') or a gz/bz2/xz compressed file.\n");
    fprintf(stderr, "\n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: need to supply a gkpStore (-G).\n");
    if (inType == TYPE_NONE)
      fprintf(stderr, "ERROR: need to supply a format type (-legacy, -coords, -hangs, -raw).\n");

    exit(1);
  }

  if (gkpStoreName)
    gkpStore = new gkStore(gkpStoreName);

  char         *S     = new char [1024];
  splitToWords  W;
  ovsOverlap    ov;

  ovFile       *of    = (ovlFileName  == NULL) ? NULL : new ovFile(ovlFileName, ovFileFullWrite);
  ovStore      *os    = (ovlStoreName == NULL) ? NULL : new ovStore(ovlStoreName, ovStoreWrite);

  for (uint32 ff=0; ff<files.size(); ff++) {
    compressedFileReader   *in = new compressedFileReader(files[ff]);

    fgets(S, 1024, in->file());

    while (!feof(in->file())) {
      W.split(S);

      switch (inType) {
      case TYPE_LEGACY:
        //  Aiid Biid 'I/N' ahang bhang erate erate
        ov.a_iid = W(0);
        ov.b_iid = W(1);

        ov.flipped(W[2][0] == 'I');

        ov.a_hang(W(3));
        ov.b_hang(W(4));

        //  Overlap store reports %error, but we expect fraction error.
        //ov.erate(atof(W[5]);  //  Don't use the original uncorrected error rate
        ov.erate(atof(W[6]) / 100.0);
        break;

      case TYPE_COORDS:
        break;
      
      case TYPE_HANGS:
        break;

      case TYPE_RAW:
        break;

      default:
        break;
      }

      if (of)
        of->writeOverlap(&ov);

      if (os)
        os->writeOverlap(&ov);

      fgets(S, 1024, in->file());
    }

    delete in;
  }

  delete    os;
  delete    of;

  delete [] S;

  delete    gkpStore;

  exit(0);
}
