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

  while (!feof(stdin)) {
    sim4polish *p = s4p_readPolish(stdin);

    if (p != 0L) {
      s4p_normalize(p);

      bool    fwd  = (p->matchOrientation == SIM4_MATCH_FORWARD);

      for (u32bit exon=0; exon<p->numExons; exon++) {
        sim4polishExon *e = p->exons + exon;

        //  Parse the alignment to find ungapped blocks

        u32bit  aPos = 0;

        u32bit  qBeg = e->estFrom - 1;
        u32bit  gBeg = e->genFrom - 1;

        if (fwd == false)
          qBeg = p->estLen - e->estTo + 1;


        bool  notDone = true;  //  There should be a way to get rid of this stupid variable....
        while (notDone) {
          notDone = ((e->estAlignment[aPos] != 0) &&
                     (e->genAlignment[aPos] != 0));

          //  If we find the end of a gapless block, emit a match

          if (e->estAlignment[aPos] != e->genAlignment[aPos]) {

            if (fwd) {
              fprintf(stdout, "%s "u32bitFMT" %c -> %s "u32bitFMT" %c\n",
                      p->estDefLine, qBeg, e->estAlignment[aPos],
                      p->genDefLine, gBeg, e->genAlignment[aPos]);
            } else {
              fprintf(stdout, "%s "u32bitFMT" %c -> %s "u32bitFMT" %c\n",
                      p->estDefLine, qBeg, e->estAlignment[aPos],
                      p->genDefLine, gBeg, e->genAlignment[aPos]);
            }
          }

          if (e->estAlignment[aPos] != '-')
            if (fwd) qBeg++;
            else     qBeg--;
          if (e->genAlignment[aPos] != '-')
            gBeg++;

          aPos++;

        }
      }
      s4p_destroyPolish(p);
    }
  }

  return(0);
}

