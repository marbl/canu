#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "bio.h"
#include "sim4.H"

//  Kaz Kylheku <kaz@ashi.footprints.net> library.
#include "kazlib/dict.h"
#include "kazlib/except.h"
#include "kazlib/hash.h"
#include "kazlib/list.h"
#include "kazlib/sfx.h"

//  Updates the IID's in a set of polishes.  If a file of deflines (or
//  fasta file) is supplied, the IIDs will match those, otherwise,
//  IIDs are guaranteed to be unique, and in the order of the
//  polishes.

const char *usage =
"usage: %s [-v] [-c x] [-g x] < polishes > polishes\n"
"     -v         Entertain the user\n"
"     -c x       Read cDNA deflines from x\n"
"     -g x       Read genomic deflines from x\n"
"  x is a fasta file, or a list of deflines (fasta file with no sequence)\n";


void
addToDict(dict_t *d, char *n) {
  dnode_t  *node = 0L;
  char     *dcpy = 0L;

  if (n == 0L)
    return;

  seqCache  *F = new seqCache(n);
  seqInCore *S = F->getSequenceInCore();
  while (S) {
    node = (dnode_t *)palloc(sizeof(dnode_t));
    dcpy = (char    *)palloc(sizeof(char) * S->headerLength() + 1);

    strcpy(dcpy, S->header());

    dnode_init(node, (void *)(unsigned long)S->getIID());
    dict_insert(d, node, dcpy);

    delete S;
    S = F->getSequenceInCore();
  }
  delete F;
}


int
headerCompare(const void *a, const void *b) {
  char  *A = *((char **)a);
  char  *B = *((char **)b);

  //fprintf(stderr, "%s -- %s\n", A, B);
  return(strcmp(A, B));
}


int
main(int argc, char **argv) {
  bool   beVerbose = false;
  char  *cDeflines = 0L;
  char  *gDeflines = 0L;

  int arg=1;
  while (arg < argc) {
    
    if        (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;
    } else if (strcmp(argv[arg], "-c") == 0) {
      cDeflines = argv[++arg];
    } else if (strcmp(argv[arg], "-g") == 0) {
      gDeflines = argv[++arg];
    } else {
      fprintf(stderr, "Unknown arg: %s\n", argv[arg]);
    }
    arg++;
  }

  if (isatty(fileno(stdin))) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  //  We parse args, then build the dictionaries, so we can do
  //  any quick error detection first.

  if (beVerbose)
    fprintf(stderr, "%s: Reading genomic deflines from '%s'\n", argv[0], gDeflines);

  dict_t  *g = dict_create(DICTCOUNT_T_MAX, headerCompare);
  addToDict(g, gDeflines);

  if (beVerbose)
    fprintf(stderr, "%s: Reading genomic deflines from '%s'\n", argv[0], cDeflines);

  dict_t  *c = dict_create(DICTCOUNT_T_MAX, headerCompare);
  addToDict(c, cDeflines);

  //  Read all the matches, changing IIDs.  If we find a defline
  //  with no IID, holler and die.

  sim4polish *p;

  //speedCounter  *C = new speedCounter("%12.0f polishes -- %12.0f polishes/second\r",
  //                                    1.0, 0xfff, beVerbose);

  while ((p = s4p_readPolish(stdin)) != 0L) {

    if ((p->estDefLine == 0L) || (p->genDefLine == 0L)) {
      s4p_printPolish(stdout, p, S4P_PRINTPOLISH_NOTVALUABLE);
      fprintf(stderr, "ERROR:  Polish has no deflines!\n");
      exit(1);
    }

    dnode_t *cid = dict_lookup(c, p->estDefLine);
    dnode_t *gid = dict_lookup(g, p->genDefLine);

    if ((cid == 0L) || (gid == 0L)) {
      const char *msg = "both deflines";
      if (cid)  msg = "genomic defline";
      if (gid)  msg = "est defline";

      s4p_printPolish(stdout, p, S4P_PRINTPOLISH_NOTVALUABLE);
      fprintf(stderr, "ERROR:  Couldn't find %s (%p %p) in the dictionary!\n", msg, cid, gid);
      exit(1);
    }

    p->estID = (u32bit)(unsigned long)dnode_get(cid);
    p->genID = (u32bit)(unsigned long)dnode_get(gid);

    s4p_printPolish(stdout, p, S4P_PRINTPOLISH_NOTVALUABLE);
    s4p_destroyPolish(p);

    //C->tick();
  }

  //delete C;

  //  We should clean up our dictionaries, but we just let the OS do it.
  exit(0);
}


