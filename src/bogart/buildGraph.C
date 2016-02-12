


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
  uint32  end;
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
};


void
loadEdges(char *edgeName, vector<bestEdge> &edges) {

  errno = 0;
  FILE  *edgeFile = fopen(edgeName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", edgeName, strerror(errno));

  char   edgeLine[1024];

  fgets(edgeLine, 1024, edgeFile);

  while (!feof(edgeFile)) {
    splitToWords  W(edgeLine);
    bestEdge      E;

    E.fr.tigID   = W(1);
    E.fr.tigType = W[2][0];
    E.fr.readID  = W(4);
    E.fr.readBgn = W(6);
    E.fr.readEnd = W(7);
    E.fr.end     = W[8][0];

    E.ahang      = W(10);
    E.ahang      = W(11);

    E.to.tigID   = W(14);
    E.fr.tigType = W[15][0];
    E.to.readID  = W(17);
    E.to.readBgn = W(19);
    E.to.readEnd = W(20);
    E.to.end     = W[21][0];

    edges.push_back(E);

    fgets(edgeLine, 1024, edgeFile);
  }

  fprintf(stderr, "Loaded "F_SIZE_T" edges from '%s'.\n", edges.size(), edgeName);
}




int
main(int argc, char **argv) {
  char              *gkpName = NULL;
  char              *tigName = NULL;
  int32              tigVers = -1;
  char              *edgesName = NULL;

  vector<bestEdge>   edges;

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

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -G gkpStore -T tigStore tigVersion -E edgesFile ...\n", argv[0]);
    fprintf(stderr, "  -G gkpStore           path to gkpStore\n");
    fprintf(stderr, "  -T tigStore version   path to tigStore\n");
    fprintf(stderr, "  -E edgeFile           path to bogart used-best-edges file\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);
  tgStore  *tigStore = new tgStore(tigName, tigVers);


  fprintf(stdout, "H\tVN:Z:bogart\n");

  for (uint32 tt=0; tt<tigStore->numTigs(); tt++) {
    tgTig  *tig = tigStore->loadTig(tt);

    fprintf(stdout, "S\ttig%08u\t*\tLN:i:%u\tRC:i:%u\n", tig->tigID(), tig->length(), tig->numberOfChildren());

    tigStore->unloadTig(tt, true);
  }



  loadEdges(edgesName, edges);

  uint32      edgeTypes[4][4] = {0};
  char       *cigar = new char [1024 * 1024];

  for (uint32 ee=0; ee<edges.size(); ee++) {

    //  Get the tigs for this edge, ignore if either is unassembled.

    tgTig  *frTig = tigStore->loadTig(edges[ee].fr.tigID);
    tgTig  *toTig = tigStore->loadTig(edges[ee].to.tigID);

    assert(frTig->_class != tgTig_noclass);
    assert(toTig->_class != tgTig_noclass);

    if ((frTig->_class == tgTig_unassembled) ||
        (toTig->_class == tgTig_unassembled)) {
#if 0
      fprintf(stderr, "SKIPPING edge %d -- fr tig %d read %d -- to tig %d read %d\n",
              ee,
              edges[ee].fr.tigID, edges[ee].fr.readID,
              edges[ee].to.tigID, edges[ee].to.readID);
#endif
      continue;
    }

    //  Find the reads we're using to anchor the tigs together.

#if 0
    fprintf(stderr, "SEARCHING edge %d -- fr tig %d read %d -- to tig %d read %d\n",
            ee,
            frTig->tigID(), edges[ee].fr.readID,
            toTig->tigID(), edges[ee].to.readID);
#endif

    tgPosition  *frRead = NULL;
    tgPosition  *toRead = NULL;
    uint32       rr;

    rr = 0;
    do {
      frRead = frTig->getChild(rr++);
    } while ((rr < frTig->numberOfChildren()) && (frRead->ident() != edges[ee].fr.readID));

    rr = 0;
    do {
      toRead = toTig->getChild(rr++);
    } while ((rr < toTig->numberOfChildren()) && (toRead->ident() != edges[ee].to.readID));


    if (frRead->ident() != edges[ee].fr.readID)
      fprintf(stderr, "WARNING:  failed to find read %u in tig %u - ejected?\n",
              edges[ee].fr.readID, edges[ee].fr.tigID);
    if (toRead->ident() != edges[ee].to.readID)
      fprintf(stderr, "WARNING:  failed to find read %u in tig %u - ejected?\n",
              edges[ee].to.readID, edges[ee].to.tigID);

    if (frRead->ident() != edges[ee].fr.readID)
      continue;
    if (toRead->ident() != edges[ee].to.readID)
      continue;

    fprintf(stderr, "WORKING -- fr tig %d read %d -- to tig %d read %d\n",
            frTig->tigID(), frRead->ident(),
            toTig->tigID(), toRead->ident());

    //  Decide orientations

    bool  frForward = frRead->isForward();
    bool  toForward = toRead->isForward();

    int32  frBgn = frTig->mapGappedToUngapped(frRead->bgn());
    int32  frEnd = frTig->mapGappedToUngapped(frRead->end());

    int32  toBgn = toTig->mapGappedToUngapped(toRead->bgn());
    int32  toEnd = toTig->mapGappedToUngapped(toRead->end());

    //  Normalize to forward/forward.

    if (frForward) {
      frBgn = frTig->length(false) - frBgn;
      frEnd = frTig->length(false) - frEnd;
    }

    if (toForward) {
      toBgn = toTig->length(false) - toBgn;
      toEnd = toTig->length(false) - toEnd;
    }

    //  Compute an offset based on the read-to-read overlap.

    int32  alignOffset = 0;

    //  Figure out the hangs of each tig.

    int32    fr5 = frBgn - 0;
    int32    fr3 = frTig->length(false) - frEnd;

    int32    to5 = toBgn - 0;
    int32    to3 = toTig->length(false) - toEnd;

    //  Find the min of each of the 5' and 3' hangs.  These tell what portions of each tig should be overlapping.

    int32    min5 = min(fr5, to5);
    int32    min3 = min(fr3, to3);

    //  Extending the read positions by these minimum extensions will then give up (more or less) the region that aligns.

    frBgn -= min5;
    frEnd += min3;

    toBgn -= min5;
    toEnd += min3;

    //  And undo the normalization.

    if (frForward) {
      frBgn = frTig->length(false) - frBgn;
      frEnd = frTig->length(false) - frEnd;
    }

    if (toForward) {
      toBgn = toTig->length(false) - toBgn;
      toEnd = toTig->length(false) - toEnd;
    }

    //  Find an overlap, convert it to a cigar string.

    //  For now, we just take the average of the two lengths.
    int32  frLen = (frBgn < frEnd) ? (frEnd - frBgn) : (frBgn - frEnd);
    int32  toLen = (toBgn < toEnd) ? (toEnd - toBgn) : (toBgn - toEnd);

    sprintf(cigar, "%dm", (frLen + toLen) / 2);

    //  Report the edge.

    fprintf(stdout, "L\ttig%08u\t%c\ttig%08d\t%c\t%s\n",
            frTig->tigID(), frForward ? '+' : '-',
            toTig->tigID(), toForward ? '+' : '-',
            cigar);
    
    tigStore->unloadTig(frTig->tigID(), true);
    tigStore->unloadTig(toTig->tigID(), true);
  }

  delete [] cigar;

#if 0
  for (uint32 ii=0; ii<4; ii++) {
    fprintf(stderr, "%8u%8u%8u%8u\n",
            edgeTypes[ii][0],
            edgeTypes[ii][1],
            edgeTypes[ii][2],
            edgeTypes[ii][3]);
  }
#endif

  delete tigStore;

  gkpStore->gkStore_close();

  exit(0);
}
