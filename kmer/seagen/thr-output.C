#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "posix.H"
#include "searchGENOME.H"

#if 0
void
statusThread(void *) {
  double finish = 0.0;

  if (config._outputPos > 0)
    finish = (config._numberOfQueries - config._outputPos) / (config._outputPos / (getTime() - config._zeroTime));

  fprintf(stderr, "O:"u32bitFMTW(7)" S:"u32bitFMTW(7)" I:"u32bitFMTW(7)" T:"u32bitFMTW(7)" (%5.1f%%; %8.3f/sec) Finish in %5.2f seconds.\r",
          outputPos,
          inputTail,
          inputHead,
          numberOfQueries,
          100.0 * outputPos / numberOfQueries,
          outputPos / (getTime() - zeroTime),
          finish);
  fflush(stderr);

  double perSec    = outputPos / (getTime() - zeroTime + 0.0000001);

  if      (perSec < 32.0)
    outputMask = 0xf;
  else if (perSec < 256.0)
    outputMask = 0x7f;
  else if (perSec < 1024.0)
    outputMask = 0x1ff;
  else
    outputMask = 0x3ff;
}
#endif



void*
writerThread(void *U, void *Q) {
  encodedQuery  *query = (encodedQuery *)Q;

  //  Write the hits
  //
  if (query->theOutputLength() > 0) {
    errno = 0;
    write(config._outputFile, query->theOutput(), query->theOutputLength());
    if (errno)
      fprintf(stderr, "Couldn't write to the output file '%s'.\n%s\n",
              config._outputFileName, strerror(errno)), exit(1);
  }

  //  Write the query match counts, too!
  //
  if (config._matchCountsFile) {
    char   str[256];

    sprintf(str, u32bitFMT"\n", query->numberOfResults());

    errno = 0;
    write(config._matchCountsFile, str, strlen(str));
    if (errno)
      fprintf(stderr, "Couldn't write to the match counts file '%s'.\n%s\n",
              config._queryMatchFileName, strerror(errno)), exit(1);
  }

  delete query;

  return(0L);
}
