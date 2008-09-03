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
  u32bit    inputsLen   = 0;
  char     *inputs[8192];

  //  Parse and check the inputs.

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-output", 2) == 0) {
      outName = argv[++arg];

    } else {
      if (tapperResultFile::validResultFile(argv[arg]) == false) {
        fprintf(stderr, "Didn't find tapperResultFile '%s'\n", argv[arg]);
        err++;
      } else {
        inputs[inputsLen++] = argv[arg];
      }
    }
    arg++;
  }
  if ((err) || (inputsLen == 0)) {
    fprintf(stderr, "usage: %s -output out-directory in-directory [in-directory ...]\n", argv[0]);
    exit(1);
  }

  //  Open the output file

  tapperResultFile     *out = new tapperResultFile(outName, 'w');

  //  Loop over the inputs, copying to the output.  We could be much
  //  looser here, just blindly copying all records in each file, but
  //  we'll be a little more careful, and copy frag by frag.

  for (u32bit inputsIdx=0; inputsIdx<inputsLen; inputsIdx++) {
    tapperResultFile *inp = new tapperResultFile(inputs[inputsIdx], 'r');
    tapperResult     *res = new tapperResult;

    while (inp->read(res))
      out->write(res);

    delete inp;
    delete res;
  }

  delete out;

  exit(0);
}
