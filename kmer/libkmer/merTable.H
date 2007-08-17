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

  void         build(seqStream *CS, u32bit merSize, u32bit merSkip=0) {

    //  Allocate the mer table
    //
    u32bit  tableSize = u32bitONE << (2*merSize);
    fprintf(stderr, "allocate "u32bitFMT" entries for a merTable.\n", tableSize);

    merToPositions = new u32bit [tableSize+1];

    //  First pass, count the number of times we see each mer
    //
    for (u32bit i=0; i<=tableSize; i++)
      merToPositions[i] = 0;

    merStream  MS(merSize, CS);

    while (MS.nextMer(merSkip)) {
      u64bit  m = (u64bit)MS.theFMer();
      //fprintf(stderr, "add mer "u64bitFMT"\n", m);
      merToPositions[m]++;
    }

    //  Convert those to indexes into positions - m[i] is the start of
    //  the locations in positions[] where positions are stored.
    //
    for (u32bit pos=0, val=0, i=0; i<=tableSize; i++) {
      val               = merToPositions[i];
      merToPositions[i] = pos;
      pos              += val;
    }

    //  Allocate space
    //
    fprintf(stderr, "allocate "u32bitFMT" entries for positions\n", merToPositions[tableSize]);
    positions = new u32bit [merToPositions[tableSize]];

    //  Second pass, fill in positions
    //
    MS.rewind();

    while (MS.nextMer(merSkip))
      positions[ merToPositions[(u64bit)MS.theFMer()]++ ] = MS.thePositionInStream();
  };

  u32bit       numberOfPositions(u64bit  mer) {
    return(merToPositions[mer+1] - merToPositions[mer]);
  };

  u32bit       getPosition(u64bit mer, u32bit index) {
    if (index >= merToPositions[mer+1] - merToPositions[mer])
      return(~u32bitZERO);
    return(merToPositions[mer] + index);
  };

private:
  u32bit    *merToPositions;  //  index into positions[]; merToPositions[mer] is the first base in the mer
  u32bit    *positions;       //  list of positions for mers, sorted by mer
};

#endif  //  MERTABLE_H
