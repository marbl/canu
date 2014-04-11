#ifndef MERTABLE_H
#define MERTABLE_H

//  The obvious simple small mer table, appropriate for large sequences

#error merTable should be unused

class merTable {
public:
  merTable() {
  };
  ~merTable() {
    delete [] merToPositions;
    delete [] positions;
  };

  void         build(seqStream *CS, uint32 merSize, uint32 merSkip=0) {

    //  Allocate the mer table
    //
    uint32  tableSize = uint32ONE << (2*merSize);
    fprintf(stderr, "allocate "uint32FMT" entries for a merTable.\n", tableSize);

    merToPositions = new uint32 [tableSize+1];

    //  First pass, count the number of times we see each mer
    //
    for (uint32 i=0; i<=tableSize; i++)
      merToPositions[i] = 0;

    merStream  MS(merSize, CS);

    while (MS.nextMer(merSkip)) {
      uint64  m = (uint64)MS.theFMer();
      //fprintf(stderr, "add mer "uint64FMT"\n", m);
      merToPositions[m]++;
    }

    //  Convert those to indexes into positions - m[i] is the start of
    //  the locations in positions[] where positions are stored.
    //
    for (uint32 pos=0, val=0, i=0; i<=tableSize; i++) {
      val               = merToPositions[i];
      merToPositions[i] = pos;
      pos              += val;
    }

    //  Allocate space
    //
    fprintf(stderr, "allocate "uint32FMT" entries for positions\n", merToPositions[tableSize]);
    positions = new uint32 [merToPositions[tableSize]];

    //  Second pass, fill in positions
    //
    MS.rewind();

    while (MS.nextMer(merSkip))
      positions[ merToPositions[(uint64)MS.theFMer()]++ ] = MS.thePositionInStream();
  };

  uint32       numberOfPositions(uint64  mer) {
    return(merToPositions[mer+1] - merToPositions[mer]);
  };

  uint32       getPosition(uint64 mer, uint32 index) {
    if (index >= merToPositions[mer+1] - merToPositions[mer])
      return(~uint32ZERO);
    return(merToPositions[mer] + index);
  };

private:
  uint32    *merToPositions;  //  index into positions[]; merToPositions[mer] is the first base in the mer
  uint32    *positions;       //  list of positions for mers, sorted by mer
};

#endif  //  MERTABLE_H
