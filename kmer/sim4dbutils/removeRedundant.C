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
    fprintf(stderr, "usage: %s [] <polishes-file>\n", argv[0]);
    fprintf(stderr, "(yes, you _must_ give it a file.  stdin is not possible.)\n");
    fprintf(stderr, "WARNING THIS IS PROTOTYPE BROKEN CODE!\n");
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "", 2) == 0) {
    } else if (strncmp(argv[arg], "", 2) == 0) {
    }

    arg++;
  }


  u32bit  matchesWithNoOverlap       = 0;
  u32bit  matchesWithOverlap         = 0;
  u32bit  notPerfectClique           = 0;

  //  Open a polishFile and force the index to build
  //
  sim4polishFile *Afile = new sim4polishFile(argv[argc-1]);
  Afile->setPosition(0);

  sim4polishWriter *writer = new sim4polishWriter("-", sim4polishS4DB);

  //  Ask both for the largest EST iid seen, then iterate over those.
  //
  u32bit  largestIID = Afile->maxIID();

  for (u32bit iid=0; iid<largestIID; iid++) {
    sim4polishList *A = Afile->getEST(iid);

    if (A->length() > 0) {

      //  fill out the overlap matrix

      olap_t  **overlap = new olap_t* [A->length()];
      overlap[0] = new olap_t [A->length() * A->length()];
      for (u32bit i=1; i<A->length(); i++)
        overlap[i] = overlap[i-1] + A->length();

      for (u32bit a=0; a<A->length(); a++)
        for (u32bit b=0; b<A->length(); b++)
          if (a == b)
            overlap[a][b] = 0;
          else
            overlap[a][b] = findOverlap((*A)[a], (*A)[b]);

      //  look for guys with no overlaps, print and remove them

      sim4polishList *W = new sim4polishList;

      for (u32bit a=0; a<A->length(); a++) {
        bool nooverlaps = true;

        for (u32bit b=0; b<A->length(); b++)
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
      fprintf(stderr, "IID="u32bitFMTW(8)" -- overlap:"u32bitFMT" noOverlap:"u32bitFMT"\r",
              iid, matchesWithOverlap, matchesWithNoOverlap);
      fflush(stderr);
#endif


      //  A is junk, W contains the matches that overlap.

      delete A;
      A = 0L;


      //  Report all the overlaps

#ifdef DEBUGOUT
      for (u32bit a=0; a<W->length(); a++) {
        sim4polish *p = (*W)[a];
        fprintf(stderr, u32bitFMTW(3)": "u32bitFMTW(3)"--"u32bitFMTW(3)"\n",
                iid, p->exons[0].genFrom, p->exons[p->numExons-1].genTo);
      }
#endif



      //  while we have matches in the set of overlapping matches,
      //  find a connected component, check that it is/is not a
      //  clique, and decide which match to keep.

      u32bit  *clique      = new u32bit [W->length()];
      u32bit   cliqueSize  = 0;
      bool     inserted    = false;
      u32bit   *length     = new u32bit [W->length()];

      while (W->length() > 0) {

#ifdef DEBUGOUT
        fprintf(stderr, "IID="u32bitFMTW(8)" -- examine "u32bitFMT" matches\n",
                iid, W->length());
#endif

        //  Find the length of all the matches in this set

        for (u32bit a=0; a<W->length(); a++) {
          length[a] = 0;
          for (u32bit i=0; i<(*W)[a]->_numExons; i++)
            length[a] += (*W)[a]->_exons[i]._genTo - (*W)[a]->_exons[i]._genFrom + 1;
        }

        //  reconstruct the overlap matrix -- hey, if you want to be
        //  efficient and recover this from the existing one, nobody is
        //  stopping you.
        
        for (u32bit a=0; a<W->length(); a++)
          for (u32bit b=0; b<W->length(); b++)
            if (a == b)
              overlap[a][b] = 0;
            else
              overlap[a][b] = findOverlap((*W)[a], (*W)[b]);

        //  OK, now find the clique/connected component

        for (u32bit i=0; i<W->length(); i++)
          clique[i] = 0;

        clique[0]   = 1;
        cliqueSize  = 1;
        inserted    = true;

        while (inserted) {
          inserted = false;

          //  If a is in the clique, add all it's overlaps

          for (u32bit a=0; a<W->length(); a++) {
            if (clique[a]) {
              for (u32bit b=0; b<W->length(); b++) {
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
        fprintf(stderr, "IID="u32bitFMTW(8)" -- examine "u32bitFMT" matches, found "u32bitFMT" overlapping\n",
                iid, W->length(), cliqueSize);
#endif

        //  Check that it is a clique

        if (cliqueSize > 2) {

          u32bit num = 0;

          for (u32bit a=0; a<W->length(); a++)
            for (u32bit b=0; b<W->length(); b++)
              if (clique[a] && clique[b] && overlap[a][b])
                num++;

          if (num != cliqueSize * (cliqueSize-1)) {
            notPerfectClique++;

            fprintf(stderr, "\nNOT A PERFECT CLIQUE!  Found "u32bitFMT" overlaps, wanted "u32bitFMT" in the clique.\n",
                    num, cliqueSize * (cliqueSize-1));

            //for (u32bit a=0; a<W->length(); a++)
            //  if (clique[a])
            //    writer->writeAlignment((*W)[a]);
          }
          
        }

        //  Find the longest member, output it

        u32bit longest = 0;
        while (clique[longest] == 0)
          longest++;

        for (u32bit i=0; i<W->length(); i++)
          if ((clique[i]) && (length[longest] < length[i]))
            longest = i;

        writer->writeAlignment((*W)[longest]);

        //  Remove the clique from the set of overlaps

        A = new sim4polishList;
        for (u32bit i=0; i<W->length(); i++) {
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

  fprintf(stderr, "\nmatches withOvl:"u32bitFMT" withoutOvl:"u32bitFMT"\n",
          matchesWithOverlap, matchesWithNoOverlap);
  fprintf(stderr, "not perfect clique:"u32bitFMT"\n", notPerfectClique);
}


