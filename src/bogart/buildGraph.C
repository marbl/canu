
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-FEB-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "tgStore.H"

#include "splitToWords.H"

#include <vector>
using namespace std;



class bestRead {
public:
  bestRead() {
  };
  ~bestRead() {
  };

  uint32  tigID;
  char    tigType;
  uint32  readID;
  uint32  readBgn;
  uint32  readEnd;
};


class bestEdge {
public:
  bestEdge() {
  };
  ~bestEdge() {
  };

  bestRead fr;
  int32    ahang;
  int32    bhang;
  bestRead to;

  bool     flipped;
};


void
loadEdges(char *edgeName, vector<bestEdge> &edges, set<uint32> &vertices) {

  errno = 0;
  FILE  *edgeFile = fopen(edgeName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", edgeName, strerror(errno));

  char   edgeLine[1024];

  fgets(edgeLine, 1024, edgeFile);

  while (!feof(edgeFile)) {
    splitToWords  W(edgeLine);
    bestEdge      E;
    uint32        w = 0;

    assert(W[w++][0] == 't');           //  'tig'
    E.fr.tigID   = W(w++);              //  tigID
    E.fr.tigType = W[w++][0];           //  tig type 'R', 'N', ...
    assert(W[w++][0] == 'r');           //  'read'
    E.fr.readID  = W(w++);              //  readID
    assert(W[w++][0] == 'a');           //  'at'
    E.fr.readBgn = W(w++);              //  bgn-position
    E.fr.readEnd = W(w++);              //  end-position

    E.ahang      = W(w++);              //  a-hang
    E.flipped    = (W[w++][0] == '<');  //  '<' if flipped, '>' if normal
    E.bhang      = W(w++);              //  b-hang

    assert(W[w++][0] == 't');           //  'tig'
    E.to.tigID   = W(w++);              //  tigID
    E.fr.tigType = W[w++][0];           //  tig type
    assert(W[w++][0] == 'r');           //  'read'
    E.to.readID  = W(w++);              //  readID
    assert(W[w++][0] == 'a');           //  'at'
    E.to.readBgn = W(w++);              //  bgn-position
    E.to.readEnd = W(w++);              //  end-position

    edges.push_back(E);

    vertices.insert(E.fr.tigID);
    vertices.insert(E.to.tigID);

    fgets(edgeLine, 1024, edgeFile);
  }

  fprintf(stderr, "Loaded "F_SIZE_T" edges from '%s'.\n", edges.size(), edgeName);
}




tgPosition *
findRead(tgTig *tig, uint32 id) {
  uint32      rr = 0;
  tgPosition *rd = NULL;

  do {
    rd = tig->getChild(rr++);
  } while ((rr < tig->numberOfChildren()) && (rd->ident() != id));

  if (rd->ident() != id) {
    fprintf(stderr, "WARNING:  failed to find read %u in tig %u - ejected?\n", id, tig->tigID());
    rd = NULL;
  }

  return(rd);
}




int
main(int argc, char **argv) {
  char              *gkpName = NULL;
  char              *tigName = NULL;
  int32              tigVers = -1;
  char              *edgesName = NULL;
  char              *graphName = NULL;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      edgesName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      graphName = argv[++arg];

    } else {
      char *s = new char [1024];
      sprintf(s, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (gkpName == NULL)
    err.push_back("No gatekeeper store (-G option) supplied.\n");
  if ((tigName == NULL) || (tigVers == -1))
    err.push_back("No tigStore store (-T option) supplied.\n");
  if (edgesName == NULL)
    err.push_back("No edges file (-E option) supplied.\n");
  if (graphName == NULL)
    err.push_back("No output graph file (-o option) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -G gkpStore -T tigStore tigVersion -E edgesFile ...\n", argv[0]);
    fprintf(stderr, "  -G gkpStore           path to gkpStore\n");
    fprintf(stderr, "  -T tigStore version   path to tigStore\n");
    fprintf(stderr, "  -E edgeFile           path to bogart unused-edges file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o graph.gfa          write to 'graph.gfa'\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //  Open output.

  errno = 0;
  FILE               *graph = fopen(graphName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output graph '%s': %s\n", graphName, strerror(errno)), exit(1);

  //  Open inputs, load the graph.

  gkStore            *gkpStore = gkStore::gkStore_open(gkpName);
  tgStore            *tigStore = new tgStore(tigName, tigVers);

  vector<bestEdge>    edges;
  set<uint32>         vertices;

  loadEdges(edgesName, edges, vertices);

  //  Dump vertcies.

  fprintf(graph, "H\tVN:Z:canu\n");

  for (uint32 tt=0; tt<tigStore->numTigs(); tt++) {
    tgTig  *tig = tigStore->loadTig(tt);

    if (vertices.count(tig->tigID()) > 0)
      fprintf(graph, "S\ttig%08u\t*\tLN:i:%u\tRC:i:%u\n", tig->tigID(), tig->length(), tig->numberOfChildren());

    tigStore->unloadTig(tt, true);
  }

  //  Dump graph.

  char       *cigar = new char [1024 * 1024];

  uint32      nEdgesUnassembled = 0;

  for (uint32 ee=0; ee<edges.size(); ee++) {

    //  Get the tigs for this edge, ignore if either is unassembled.

    tgTig  *frTig   = tigStore->loadTig(edges[ee].fr.tigID);
    uint32  frTigID = frTig->tigID();

    tgTig  *toTig   = tigStore->loadTig(edges[ee].to.tigID);
    uint32  toTigID = toTig->tigID();

    assert(frTig->_class != tgTig_noclass);
    assert(toTig->_class != tgTig_noclass);

    if ((frTig->_class == tgTig_unassembled) ||
        (toTig->_class == tgTig_unassembled)) {
      nEdgesUnassembled++;
      continue;
    }

    //  Find the reads we're using to anchor the tigs together.

    tgPosition  *frRead = findRead(frTig, edges[ee].fr.readID);
    tgPosition  *toRead = findRead(toTig, edges[ee].to.readID);

    if ((frRead == NULL) ||
        (toRead == NULL))
      continue;

    //  Map coordinates from gapped to ungapped.

    uint32  frReadMin = frTig->mapGappedToUngapped(frRead->min());
    uint32  frReadMax = frTig->mapGappedToUngapped(frRead->max());
    uint32  frLen     = frTig->length(false);

    uint32  toReadMin = toTig->mapGappedToUngapped(toRead->min());
    uint32  toReadMax = toTig->mapGappedToUngapped(toRead->max());
    uint32  toLen     = toTig->length(false);

    //
    //  Convert from a read-read overlap to a tig-tig overlap.
    //

    //  Orient tigs based oin the read orientation.
    //    For 'fr', we require that the read always be forward.
    //    For 'to', the overlap dictates the orientation.

    bool  frTigFwd = frRead->isForward();  //  Tig forward if read forward.
    bool  toTigFwd = toRead->isForward();  //  Same, unless...

    if (edges[ee].flipped == true)         //  ...edge is flipped, so flip
      toTigFwd = !toTigFwd;                //  'to' tig.

    //  Cleanup.  Makes skipping an edge much easier.

    tigStore->unloadTig(frTigID, true);  frRead = NULL;
    tigStore->unloadTig(toTigID, true);  toRead = NULL;

    //  Based on tig orientation, find the bgn and end lengths from each read.

    int32   frBgn = (frTigFwd) ? (        frReadMin) : (frLen - frReadMax);
    int32   frEnd = (frTigFwd) ? (frLen - frReadMax) : (        frReadMin);

    int32   toBgn = (toTigFwd) ? (        toReadMin) : (toLen - toReadMax);
    int32   toEnd = (toTigFwd) ? (toLen - toReadMax) : (        toReadMin);

    //fprintf(graph, "hangs0- fr %d-%d to %d-%d ahang %d bhang %d\n",
    //        frBgn, frEnd, toBgn, toEnd, edges[ee].ahang, edges[ee].bhang);

    //  Apply the overlap hangs to find the overlapping regions on the tigs.

    if (edges[ee].ahang < 0)
      toBgn += -edges[ee].ahang;
    else
      frBgn +=  edges[ee].ahang;

    if (edges[ee].bhang < 0)
      frEnd += -edges[ee].bhang;
    else
      toEnd +=  edges[ee].bhang;

    //fprintf(graph, "hangs1- fr %d-%d to %d-%d\n",
    //        frBgn, frEnd, toBgn, toEnd);

    //  The overlap is now between regions toBgn-toEnd and frBgn-frEnd.  Extend this to cover the ends of each tig.
    //
    //     ------------------------------------------
    //                           +++ ||| olap ||| +++
    //                           -------------------------------

    if (toBgn < frBgn) {
      toBgn -= toBgn;
      frBgn -= toBgn;
    } else {
      toBgn -= frBgn;
      frBgn -= frBgn;
    }

    if (toEnd < frEnd) {
      toEnd -= toEnd;
      frEnd -= toEnd;
    } else {
      toEnd -= frEnd;
      frEnd -= frEnd;
    }

    //fprintf(graph, "hangs2- fr %d-%d to %d-%d\n",
    //        frBgn, frEnd, toBgn, toEnd);

    //  Compute the alignment between the two regions, and convert to a cigar string.

    frLen -= (frBgn + frEnd);
    toLen -= (toBgn + toEnd);

    sprintf(cigar, "%dM", (frLen + toLen) / 2);    //  Used to be 'm', Bandage complained about it not being 'M'.

    //  The overlap should now have one of:
    //     frBgn == toEnd == 0 -- to has an overlap to fr
    //     frEnd == toBgn == 0 -- fr has an overlap to to
    //
    //  If not, the overlap is inconsistent with the tigs; it implies the two tigs overlap in their
    //  entirety.

    //  GFA requires that the overlap be between the end of the first read and the start of the second read.
    //  Flip the order if needed.

    if (frTigID == toTigID) {
      fprintf(stderr, "L\ttig%08u\t%c\ttig%08d\t%c\t%s circular\n",
              frTigID, (frTigFwd) ? '+' : '-',
              toTigID, (toTigFwd) ? '+' : '-',
              cigar);
      continue;
    }

    if ((toBgn == 0) && (frEnd == 0)) {
      fprintf(graph, "L\ttig%08u\t%c\ttig%08d\t%c\t%s\n",
              frTigID, (frTigFwd) ? '+' : '-',
              toTigID, (toTigFwd) ? '+' : '-',
              cigar);
      continue;
    }

    if ((frBgn == 0) && (toEnd == 0)) {
      fprintf(graph, "L\ttig%08u\t%c\ttig%08d\t%c\t%s\n",
              toTigID, (toTigFwd) ? '+' : '-',
              frTigID, (frTigFwd) ? '+' : '-',
              cigar);
      continue;
    }

    //  Inconsistent edge.

    fprintf(stderr, "L\ttig%08u\t%c\ttig%08d\t%c\t%s inconsistent\n",
            frTigID, (frTigFwd) ? '+' : '-',
            toTigID, (toTigFwd) ? '+' : '-',
            cigar);
  }

  edges.clear();       //  Make valgrind slightly happier.
  vertices.clear();

  delete [] cigar;
  delete    tigStore;

  gkpStore->gkStore_close();

  exit(0);
}
