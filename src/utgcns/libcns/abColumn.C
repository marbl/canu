
#include "abAbacus.H"


void
abColumn::CheckBaseCounts(abAbacus *abacus) {
  uint32 counts[CNS_NALPHABET] = {0};

  //  Why skip the last column??
  //if (next == -1)
  //  return;

  abBead *bead = abacus->getBead(callID());

  while (bead->downID().isValid()) {
    bead = abacus->getBead(bead->downID());

    counts[ abacus->getBase(bead->ident()) ]++;
  }

  if (counts['A'] != GetColumnBaseCount('A'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d A %d != %d\n", ident().get(), counts['A'], GetColumnBaseCount('A'));

  if (counts['C'] != GetColumnBaseCount('C'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d C %d != %d\n", ident().get(), counts['C'], GetColumnBaseCount('C'));

  if (counts['G'] != GetColumnBaseCount('G'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d G %d != %d\n", ident().get(), counts['G'], GetColumnBaseCount('G'));

  if (counts['T'] != GetColumnBaseCount('T'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d T %d != %d\n", ident().get(), counts['T'], GetColumnBaseCount('T'));

  if (counts['-'] != GetColumnBaseCount('-'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d - %d != %d\n", ident().get(), counts['-'], GetColumnBaseCount('-'));
};
