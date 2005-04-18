
//  Loads a set of matches from a file

class matchList {
public:
  matchList(char *filename);
  ~matchList() {
    delete _seq1;
    delete _seq2;
    delete [] _matches;
  };

  char          _file1[1024];
  char          _file2[1024];

  FastAWrapper *_seq1;
  FastAWrapper *_seq2;

  u32bit        _matchesLen;
  u32bit        _matchesMax;
  match_t      *_matches;
};





matchList::matchList(char *filename) {

  errno = 0;
  FILE   *inFile = fopen(filename, "r");
  if (errno)
    fprintf(stderr, "matchList::matchList()-- failed to load %s: %s\n", filename, strerror(errno)), exit(1);


  //  While loading matches, we compute the mapped length and covered
  //  length.

  //fprintf(stderr, "Loading matches from %s\n", filename);
  
  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in file1 and file2.
  //
  char    inLine[1024];
  readHeader(inLine, inFile, _file1, _file2, 0L);

  //fprintf(stderr, "Opening '%s' for sequence one.\n", _file1);
  //fprintf(stderr, "Opening '%s' for sequence two.\n", _file2);

  //  Open some FastAWrappers for each of the files -- we use these
  //  only to get the length of the sequence.
  //
  _seq1 = new FastAWrapper(_file1);
  _seq2 = new FastAWrapper(_file2);

  _seq1->openIndex();
  _seq2->openIndex();

  _matchesLen = 0;
  _matchesMax = 2 * 1048576;
  _matches    = new match_t [_matchesMax];

  //  For the coverage to work correctly, we need to either have one
  //  intervalList per input sequence, or build a table of the chained
  //  sequence positions.
  //
  u64bit  *offset1 = new u64bit [_seq1->getNumberOfSequences()];
  u64bit  *offset2 = new u64bit [_seq2->getNumberOfSequences()];

  offset1[0] = 1000000;
  for (u32bit i=1; i<_seq1->getNumberOfSequences(); i++)
    offset1[i] = offset1[i-1] + _seq1->sequenceLength(i-1) + 1;

  offset2[0] = 1000000;
  for (u32bit i=1; i<_seq2->getNumberOfSequences(); i++)
    offset2[i] = offset2[i-1] + _seq2->sequenceLength(i-1) + 1;

  u32bit        skippedCount = 0;
  u64bit        skippedLengthA = 0;
  u64bit        skippedLengthB = 0;

  intervalList  intervalA;
  intervalList  intervalB;

  while (!feof(inFile)) {
    if (inLine[0] == 'M') {
      splitToWords  S(inLine);

      if ((S[1][0] == 'u') || (S[1][0] == 'x')) {
        //if ((S[1][0] == 'r')) {
        u32bit  iid1=0, pos1=0, len1=0, ori1=0;
        u32bit  iid2=0, pos2=0, len2=0, ori2=0;
        decodeMatch(S, iid1, pos1, len1, ori1, iid2, pos2, len2, ori2);

        if ((pos1 + len1) > _seq1->sequenceLength(iid1)) {
          chomp(inLine);
          fprintf(stderr, "Match longer than sequence in 1: "u32bitFMT" %s\n", _seq1->sequenceLength(iid1), inLine);
        }

        if ((pos2 + len2) > _seq2->sequenceLength(iid2)) {
          chomp(inLine);
          fprintf(stderr, "Match longer than sequence in 2: "u32bitFMT" %s\n", _seq2->sequenceLength(iid2), inLine);
        }

        if ((iid1 >= _seq1->getNumberOfSequences()) || (iid2 >= _seq2->getNumberOfSequences())) {
          //  Hmmm.  Skip it.
          skippedCount++;
          skippedLengthA += len1;
          skippedLengthB += len2;
        } else {
          intervalA.add(offset1[iid1] + (u64bit)pos1, (u64bit)len1);
          intervalB.add(offset2[iid2] + (u64bit)pos2, (u64bit)len2);

          //  Add it to our list of matches
          //
          if (_matchesLen > _matchesMax) {
            fprintf(stderr, "SORRY!  I don't feel like reallocating matches.  Increase\n");
            fprintf(stderr, "the preallocated size in %s\n", __FILE__);
            exit(1);
          }

          _matches[_matchesLen].matchiid = _matchesLen;
          _matches[_matchesLen].matchuid = strtou32bit(S[2], 0L);

          _matches[_matchesLen].iid1 = iid1;
          _matches[_matchesLen].pos1 = pos1;
          _matches[_matchesLen].len1 = len1;
          _matches[_matchesLen].ori1 = ori1;

          _matches[_matchesLen].iid2 = iid2;
          _matches[_matchesLen].pos2 = pos2;
          _matches[_matchesLen].len2 = len2;
          _matches[_matchesLen].ori2 = ori2;

          _matchesLen++;
        }
      }
    }

    fgets(inLine, 1024, inFile);
  }

  fprintf(stderr, "skipped "u32bitFMT" matches with length "u64bitFMT" and "u64bitFMT"\n",
          skippedCount, skippedLengthA, skippedLengthB);

  fprintf(stderr, "intervalLength A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)intervalA.sumOfLengths(),
          (u64bit)intervalB.sumOfLengths());

  intervalA.merge();
  intervalB.merge();

  fprintf(stderr, "coveredLength  A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)intervalA.sumOfLengths(),
          (u64bit)intervalB.sumOfLengths());
}

