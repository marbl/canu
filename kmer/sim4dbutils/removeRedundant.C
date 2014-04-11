#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "bio++.H"
#include "sim4.H"
#include "s4p_overlap.H"

//  Remove redundant polishes from an input set.
//
//  Redundancy is defined as two polishes that overlap on the genome.
//  Any amount of overlap is redundant.
//
//  The longest of the overlapping matches is saved.

//#define DEBUGOUT

int
main(int argc, char **argv) {

  if (argc < 2) {
    fprintf(stderr, "usage: %s [-gff3] <polishes-file>\n", argv[0]);
    fprintf(stderr, "(yes, you _must_ give it a file.  stdin is not possible.)\n");
    fprintf(stderr, "WARNING THIS IS PROTOTYPE BROKEN CODE!\n");
    exit(1);
  }

  sim4polishStyle  wstyle = sim4polishStyleDefault;
  sim4polishStyle  rstyle = sim4polishStyleDefault;

  int arg = 1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-gff3") == 0) {
       wstyle = sim4polishGFF3;
    }

    arg++;
  }


  uint32  matchesWithNoOverlap       = 0;
  uint32  matchesWithOverlap         = 0;
  uint32  notPerfectClique           = 0;

  //  Open a polishFile and force the index to build
  //  First find the input file type, with a hack

  sim4polishReader *reader = new sim4polishReader(argv[argc-1]);
  rstyle = reader->getsim4polishStyle();
  delete reader;

  sim4polishFile   *Afile  = new sim4polishFile(argv[argc-1], rstyle);
  Afile->setPosition(0);

  sim4polishWriter *writer = new sim4polishWriter("-", wstyle);

  if (rstyle != wstyle)
    fprintf(stderr, "warning: input format and output format differ.\n");


  //  Ask both for the largest EST iid seen, then iterate over those.
  //
  uint32  largestIID = Afile->maxIID();

  for (uint32 iid=0; iid<largestIID; iid++) {
    sim4polishList *A = Afile->getEST(iid);

    if (A->length() > 0) {

      //  fill out the overlap matrix

      olap_t  **overlap = new olap_t* [A->length()];
      overlap[0] = new olap_t [A->length() * A->length()];
      for (uint32 i=1; i<A->length(); i++)
        overlap[i] = overlap[i-1] + A->length();

      for (uint32 a=0; a<A->length(); a++)
        for (uint32 b=0; b<A->length(); b++)
          if (a == b)
            overlap[a][b] = 0;
          else
            overlap[a][b] = findOverlap((*A)[a], (*A)[b]);

      //  look for guys with no overlaps, print and remove them

      sim4polishList *W = new sim4polishList;

      for (uint32 a=0; a<A->length(); a++) {
        bool nooverlaps = true;

        for (uint32 b=0; b<A->length(); b++)
          if (overlap[a][b])
            nooverlaps = false;

        if (nooverlaps) {
          matchesWithNoOverlap++;

          writer->writeAlignment((*A)[a]);
        } else {
          matchesWithOverlap++;
          W->push(new sim4polish((*A)[a]));
        }
      }


#if 1
      fprintf(stderr, "IID="uint32FMTW(8)" -- overlap:"uint32FMT" noOverlap:"uint32FMT"\r",
              iid, matchesWithOverlap, matchesWithNoOverlap);
      fflush(stderr);
#endif


      //  A is junk, W contains the matches that overlap.

      delete A;
      A = 0L;


      //  Report all the overlaps

#ifdef DEBUGOUT
      for (uint32 a=0; a<W->length(); a++) {
        sim4polish *p = (*W)[a];
        fprintf(stderr, uint32FMTW(3)": "uint32FMTW(3)"--"uint32FMTW(3)"\n",
                iid, p->exons[0].genFrom, p->exons[p->numExons-1].genTo);
      }
#endif



      //  while we have matches in the set of overlapping matches,
      //  find a connected component, check that it is/is not a
      //  clique, and decide which match to keep.

      uint32  *clique      = new uint32 [W->length()];
      uint32   cliqueSize  = 0;
      bool     inserted    = false;
      uint32   *length     = new uint32 [W->length()];

      while (W->length() > 0) {

#ifdef DEBUGOUT
        fprintf(stderr, "IID="uint32FMTW(8)" -- examine "uint32FMT" matches\n",
                iid, W->length());
#endif

        //  Find the length of all the matches in this set

        for (uint32 a=0; a<W->length(); a++) {
          length[a] = 0;
          for (uint32 i=0; i<(*W)[a]->_numExons; i++)
            length[a] += (*W)[a]->_exons[i]._genTo - (*W)[a]->_exons[i]._genFrom + 1;
        }

        //  reconstruct the overlap matrix -- hey, if you want to be
        //  efficient and recover this from the existing one, nobody is
        //  stopping you.
        
        for (uint32 a=0; a<W->length(); a++)
          for (uint32 b=0; b<W->length(); b++)
            if (a == b)
              overlap[a][b] = 0;
            else
              overlap[a][b] = findOverlap((*W)[a], (*W)[b]);

        //  OK, now find the clique/connected component

        for (uint32 i=0; i<W->length(); i++)
          clique[i] = 0;

        clique[0]   = 1;
        cliqueSize  = 1;
        inserted    = true;

        while (inserted) {
          inserted = false;

          //  If a is in the clique, add all it's overlaps

          for (uint32 a=0; a<W->length(); a++) {
            if (clique[a]) {
              for (uint32 b=0; b<W->length(); b++) {
                if ((overlap[a][b]) && (!clique[b])) {
                  clique[b] = 1;
                  cliqueSize++;
                  inserted  = true;
                }
              }
            }
          }
        }

#ifdef DEBUGOUT
        fprintf(stderr, "IID="uint32FMTW(8)" -- examine "uint32FMT" matches, found "uint32FMT" overlapping\n",
                iid, W->length(), cliqueSize);
#endif

        //  Check that it is a clique

        if (cliqueSize > 2) {

          uint32 num = 0;

          for (uint32 a=0; a<W->length(); a++)
            for (uint32 b=0; b<W->length(); b++)
              if (clique[a] && clique[b] && overlap[a][b])
                num++;

          if (num != cliqueSize * (cliqueSize-1)) {
            notPerfectClique++;

            fprintf(stderr, "\nNOT A PERFECT CLIQUE!  Found "uint32FMT" overlaps, wanted "uint32FMT" in the clique.\n",
                    num, cliqueSize * (cliqueSize-1));

            //for (uint32 a=0; a<W->length(); a++)
            //  if (clique[a])
            //    writer->writeAlignment((*W)[a]);
          }
          
        }

        //  Find the longest member, output it

        uint32 longest = 0;
        while (clique[longest] == 0)
          longest++;

        for (uint32 i=0; i<W->length(); i++)
          if ((clique[i]) && (length[longest] < length[i]))
            longest = i;

        writer->writeAlignment((*W)[longest]);

        //  Remove the clique from the set of overlaps

        A = new sim4polishList;
        for (uint32 i=0; i<W->length(); i++) {
          if (clique[i] == 0)
            A->push(new sim4polish((*W)[i]));
        }

        delete W;
        W = A;
        A = 0L;
      }

      delete [] clique;
      delete    W;

      delete [] overlap[0];
      delete [] overlap;
    }

    delete A;
  }

  delete writer;
  delete Afile;

  fprintf(stderr, "\nmatches withOvl:"uint32FMT" withoutOvl:"uint32FMT"\n",
          matchesWithOverlap, matchesWithNoOverlap);
  fprintf(stderr, "not perfect clique:"uint32FMT"\n", notPerfectClique);
}


