#include "tapperTag.H"
#include "tapperResult.H"
#include "tapperAlignment.H"
#include "tapperHit.H"
#include "tapperGlobalData.H"
#include "tapperThreadData.H"
#include "tapperComputation.H"

int
main(int argc, char **argv) {
  char     *outName = 0L;
  u32bit    numInputs = 0;
  char      fileName[FILENAME_MAX];

  //  Parse and check the inputs.

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-output", 2) == 0) {
      outName = argv[++arg];

    } else {
      sprintf(fileName, "%s.tapperMappedIndex", argv[arg]);
      numInputs++;

      if (fileExists(fileName) == false)
        err++;
    }
    arg++;
  }
  if ((err) || (numInputs == 0)) {
    fprintf(stderr, "usage: %s -prefix outputprefix inprefix [inprefix ...]\n", argv[0]);
    exit(1);
  }

  //  Open the output file

  tapperResultFile     *out = new tapperResultFile(outName, 'w');

  //  Loop over the inputs, copying to the output.  We could be much
  //  looser here, just blindly copying all records in each file, but
  //  we'll be a little more careful, and copy frag by frag.

  arg=1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-output", 2) == 0) {
      //  Skip the output.
      arg++;
    } else {
      tapperResultFile *inp = new tapperResultFile(argv[arg], 'r');
      tapperResult     *res = new tapperResult;

      while (inp->read(res))
        out->write(res);

      delete inp;
      delete res;
    }

    arg++;
  }

  delete out;

  exit(0);
}
