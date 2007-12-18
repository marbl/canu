#include "util++.H"
#include "bio++.H"

#include "markFile.H"

int
main(int argc, char **argv) {
  char           *outName       = 0L;
  markFileMerger  mm;
  
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];
    } else if (fileExists(argv[arg])) {
      mm.addFile(argv[arg]);
    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      exit(1);
    }

    arg++;
  }

  if (outName == 0L)
    fprintf(stderr, "ERROR: no output file supplied!\n"), exit(1);

  mm.doMerge(outName);

  return(0);
}
