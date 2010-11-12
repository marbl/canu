#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"
#include "util++.H"

int
main(int argc, char **argv) {
  char     *outPrefix = 0L;
  char      datName[FILENAME_MAX];
  char      gnuName[FILENAME_MAX];
  char      pngName[FILENAME_MAX];
  char      gnuCmd[FILENAME_MAX];
  char     *inName = 0L;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-o", 2) == 0) {
      outPrefix = argv[++arg];

    } else if (strncmp(argv[arg], "-i", 2) == 0) {
      inName = argv[++arg];

    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((inName == 0L) || (outPrefix == 0L) || (err != 0)) {
    fprintf(stderr, "usage: %s -i sim4db -o outputPrefix\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  Creates outputPrefix.dat containing the number of errors at each\n");
    fprintf(stderr, "  base position, and outputPrefix.png the visual representation.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Suggested usage:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  snapper2\n");
    fprintf(stderr, "    -queries Q.fasta\n");
    fprintf(stderr, "    -genomic G.fasta\n");
    fprintf(stderr, "    -positions G.posDB\n");
    fprintf(stderr, "    -aligns\n");
    fprintf(stderr, "    -minmatchidentity 94\n");
    fprintf(stderr, "    -minmatchcoverage 90\n");
    fprintf(stderr, "    -mersize 18\n");
    fprintf(stderr, "    -ignore 500\n");
    fprintf(stderr, "    -numthreads 16\n");
    fprintf(stderr, "    -verbose\n");
    fprintf(stderr, "    -output Q.sim4db\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  pickBestPolish < Q.sim4db > Q.best.sim4db\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  reportAlignmentDifferences\n");
    fprintf(stderr, "    -i Q.best.sim4db\n");
    fprintf(stderr, "    -o Q\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  fprintf(stderr, "Reading input from '%s'\n", inName);
  fprintf(stderr, "Writing output to '%s'\n", outPrefix);

  //  Open output files early, in case they fail.

  errno = 0;

  sprintf(datName, "%s.dat", outPrefix);
  sprintf(gnuName, "%s.gnuplot", outPrefix);
  sprintf(pngName, "%s.png", outPrefix);

  FILE *DAT = fopen(datName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing data: %s\n", datName, strerror(errno)), exit(1);

  FILE *GNU = fopen(gnuName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing gnuplot command: %s\n", gnuName, strerror(errno)), exit(1);

  //  Read matches.

  u32bit            lMax = 10240;
  u32bit            lLen = 0;

  u32bit           *nTot = new u32bit [lMax];
  u32bit           *nIde = new u32bit [lMax];
  u32bit           *nMis = new u32bit [lMax];
  u32bit           *nIns = new u32bit [lMax];
  u32bit           *nDel = new u32bit [lMax];

  memset(nTot, 0, sizeof(u32bit) * lMax);
  memset(nIde, 0, sizeof(u32bit) * lMax);
  memset(nMis, 0, sizeof(u32bit) * lMax);
  memset(nIns, 0, sizeof(u32bit) * lMax);
  memset(nDel, 0, sizeof(u32bit) * lMax);

  sim4polishReader *R = new sim4polishReader(inName);
  sim4polish       *p = 0L;

  while (R->nextAlignment(p)) {
    bool    fwd  = (p->_matchOrientation == SIM4_MATCH_FORWARD);

    for (u32bit exon=0; exon<p->_numExons; exon++) {
      sim4polishExon *e = p->_exons + exon;

      //  Parse the alignment to find ungapped blocks

      u32bit  aPos = 0;                //  Position in the alignment
      u32bit  qPos = e->_estFrom - 1;  //  Actual position in the query sequence
      u32bit  gPos = e->_genFrom - 1;  //  Actual position in the genome sequence

      if (fwd == false)
        qPos = p->_estLen - e->_estFrom + 1;


      bool  notDone = true;  //  There should be a way to get rid of this stupid variable....
      while (notDone) {
        notDone = ((e->_estAlignment[aPos] != 0) &&
                   (e->_genAlignment[aPos] != 0));

        //  If we find the end of a gapless block, emit a match

        if        (e->_estAlignment[aPos] == e->_genAlignment[aPos])
          nIde[qPos]++;

        else if (e->_estAlignment[aPos] == '-')
          nDel[qPos]++;

        else if (e->_genAlignment[aPos] == '-')
          nIns[qPos]++;

        else
          nMis[qPos]++;

        nTot[qPos]++;

        assert(qPos < lMax);

        if (lLen < qPos)
          lLen = qPos;

        //fprintf(stdout, "%s "u32bitFMT" %c ->_ %s "u32bitFMT" %c\n",
        //        p->_estDefLine, qPos, e->_estAlignment[aPos],
        //        p->_genDefLine, gPos, e->_genAlignment[aPos]);

        if (e->_estAlignment[aPos] != '-')
          if (fwd) qPos++;
          else     qPos--;
        if (e->_genAlignment[aPos] != '-')
          gPos++;

        aPos++;
      }
    }
  }


  //  Index
  //  nTot
  //  nIde, percent
  //  nDel,   percent
  //  nIns,   percent
  //  nMis,   percent

  fprintf(DAT, "#idx\tnTot\tnIde\tfrac\tnDel\tfrac\tnIns\tfrac\tnMis\tfrac\tnErr\tfrac\n");
  for (u32bit i=0; i<=lLen; i++)
    fprintf(DAT, "%u\t%u\t%u\t%6.4f\t%u\t%6.4f\t%u\t%6.4f\t%u\t%6.4f\t%u\t%6.4f\n",
            i,
            nTot[i],
            nIde[i], (double)nIde[i] / nTot[i],
            nDel[i], (double)nDel[i] / nTot[i],
            nIns[i], (double)nIns[i] / nTot[i],
            nMis[i], (double)nMis[i] / nTot[i],
            nTot[i] - nIde[i], (double)(nTot[i] - nIde[i]) / nTot[i]);

  fprintf(GNU, "set terminal png\n");
  fprintf(GNU, "set output \"%s\"\n", pngName);
  fprintf(GNU, "set title \"Fraction error per base for '%s'\"\n", inName);
  fprintf(GNU, "set xlabel \"Base position\"\n");
  fprintf(GNU, "set ylabel \"Fraction error\"\n");
  fprintf(GNU, "plot [][0:0.1] \\\n");
  fprintf(GNU, "  \"%s\" using 1:4  with lines title \"nTot\", \\\n", datName);
  fprintf(GNU, "  \"%s\" using 1:6  with lines title \"nDel\", \\\n", datName);
  fprintf(GNU, "  \"%s\" using 1:8  with lines title \"nIns\", \\\n", datName);
  fprintf(GNU, "  \"%s\" using 1:10 with lines title \"nMis\", \\\n", datName);
  fprintf(GNU, "  \"%s\" using 1:12 with lines title \"nErr\"\n", datName);

  fclose(DAT);
  fclose(GNU);

  sprintf(gnuCmd, "gnuplot < %s", gnuName);
  system(gnuCmd);

  return(0);
}

