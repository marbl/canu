
#include "abAbacus.H"


void
abColumn::CheckBaseCounts(abAbacus *abacus) {
  uint32 counts[256] = {0};
  uint32 nBeads = 0;

  //  Why skip the last column??
  //if (next == -1)
  //  return;

  abBead *bead = abacus->getBead(callID());

  while (bead->downID().isValid()) {
    bead = abacus->getBead(bead->downID());

    counts[ abacus->getBase(bead->ident()) ]++;

    nBeads++;
  }

  if (counts['A'] != GetColumnBaseCount('A'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u A computed %u != stored %u (%u beads)\n", ident().get(), counts['A'], GetColumnBaseCount('A'), nBeads);

  if (counts['C'] != GetColumnBaseCount('C'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u C computed %u != stored %u (%u beads)\n", ident().get(), counts['C'], GetColumnBaseCount('C'), nBeads);

  if (counts['G'] != GetColumnBaseCount('G'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u G computed %u != stored %u (%u beads)\n", ident().get(), counts['G'], GetColumnBaseCount('G'), nBeads);

  if (counts['T'] != GetColumnBaseCount('T'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u T computed %u != stored %u (%u beads)\n", ident().get(), counts['T'], GetColumnBaseCount('T'), nBeads);

  if (counts['-'] != GetColumnBaseCount('-'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%u - computed %u != stored %u (%u beads)\n", ident().get(), counts['-'], GetColumnBaseCount('-'), nBeads);
};
