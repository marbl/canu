#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"

//  Writes the location of differences in sim4db alignments.

int
main(int argc, char **argv) {

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-1", 2) == 0) {
    } else if (strncmp(argv[arg], "-2", 2) == 0) {
    } else if (strncmp(argv[arg], "-f", 2) == 0) {
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
    }
    arg++;
  }

  sim4polish *p = new sim4polish(stdin);
  while (p->_numExons > 0) {
    bool    fwd  = (p->_matchOrientation == SIM4_MATCH_FORWARD);

    for (u32bit exon=0; exon<p->_numExons; exon++) {
      sim4polishExon *e = p->_exons + exon;

      //  Parse the alignment to find ungapped blocks

      u32bit  aPos = 0;

      u32bit  qBeg = e->_estFrom - 1;
      u32bit  gBeg = e->_genFrom - 1;

      if (fwd == false)
        qBeg = p->_estLen - e->_estFrom + 1;


      bool  notDone = true;  //  There should be a way to get rid of this stupid variable....
      while (notDone) {
        notDone = ((e->_estAlignment[aPos] != 0) &&
                   (e->_genAlignment[aPos] != 0));

        //  If we find the end of a gapless block, emit a match

        if (e->_estAlignment[aPos] != e->_genAlignment[aPos]) {

          if (fwd) {
            fprintf(stdout, "%s "u32bitFMT" %c ->_ %s "u32bitFMT" %c\n",
                    p->_estDefLine, qBeg, e->_estAlignment[aPos],
                    p->_genDefLine, gBeg, e->_genAlignment[aPos]);
          } else {
            fprintf(stdout, "%s "u32bitFMT" %c ->_ %s "u32bitFMT" %c\n",
                    p->_estDefLine, qBeg, e->_estAlignment[aPos],
                    p->_genDefLine, gBeg, e->_genAlignment[aPos]);
          }
        }

        if (e->_estAlignment[aPos] != '-')
          if (fwd) qBeg++;
          else     qBeg--;
        if (e->_genAlignment[aPos] != '-')
          gBeg++;

        aPos++;

      }
    }

    delete p;
  }

  return(0);
}

