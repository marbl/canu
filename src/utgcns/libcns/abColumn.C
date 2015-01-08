
void
abColumn::CheckBaseCounts(void) {
  uint32 counts[CNS_NALPHABET] = {0};

  //  Why skip the last column??
  //if (next == -1)
  //  return;

  abBead *bead = ab->getBead(call);

  while (bead->down.isValid()) {
    bead = ab->getBead(bead->down);

    counts[ ab->getBase(bead->soffset) ]++;
  }

  if (counts['A'] != GetColumnBaseCount('A'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d A %d != %d\n", lid, counts['A'], GetColumnBaseCount('A'));

  if (counts['C'] != GetColumnBaseCount('C'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d C %d != %d\n", lid, counts['C'], GetColumnBaseCount('C'));

  if (counts['G'] != GetColumnBaseCount('G'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d G %d != %d\n", lid, counts['G'], GetColumnBaseCount('G'));

  if (counts['T'] != GetColumnBaseCount('T'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d T %d != %d\n", lid, counts['T'], GetColumnBaseCount('T'));

  if (counts['-'] != GetColumnBaseCount('-'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d - %d != %d\n", lid, counts['-'], GetColumnBaseCount('-'));
};
