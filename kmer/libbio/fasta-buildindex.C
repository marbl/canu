#include <ctype.h>

#include "bio++.H"

char*
FastAWrapper::indexTypeNames(u32bit indextype) {

  switch (indextype) {
    case FASTA_INDEX_ANY:
      return("any-index");
      break;
    case FASTA_INDEX_ONLY:
      return("index-only");
      break;
    case FASTA_INDEX_PLUS_IDS:
      return("index-plus-names");
      break;
    case FASTA_INDEX_PLUS_DEFLINES:
      return("index-plus-deflines");
      break;
    case FASTA_INDEX_ANY | FASTA_INDEX_MD5:
      return("any-index-md5");
      break;
    case FASTA_INDEX_ONLY | FASTA_INDEX_MD5:
      return("index-only-md5");
      break;
    case FASTA_INDEX_PLUS_IDS | FASTA_INDEX_MD5:
      return("index-plus-names-md5");
      break;
    case FASTA_INDEX_PLUS_DEFLINES | FASTA_INDEX_MD5:
      return("index-plus-deflines-md5");
      break;
  }

  return("Unknown Index Type");
}






void
FastAWrapper::createIndex(u32bit indextype) {


  //  If the user didn't specify which type of index to build, build
  //  just the basic one.
  //
  if ((indextype & FASTA_INDEX_MASK) == FASTA_INDEX_ANY)
    indextype |= FASTA_INDEX_ONLY;

  //  Clear the description structure
  //
  _idxfa_global  theGlobalDesc;

  theGlobalDesc._magic                 = FASTA_MAGICNUMBER;
  theGlobalDesc._version               = FASTA_VERSIONNUMBER;
  theGlobalDesc._indexType             = indextype;
  theGlobalDesc._fastaType             = 0;                  //  XXX:  Where to get this??
  theGlobalDesc._numberOfSequences     = 0;
  theGlobalDesc._fastaFileSize         = 0;
  theGlobalDesc._fastaModificationTime = 0;
  theGlobalDesc._fastaCreationTime     = 0;
  theGlobalDesc._seqlineLength         = 0;
  theGlobalDesc._seqendlLength         = 0;
  theGlobalDesc._fixedWidth            = true;
  theGlobalDesc._squeezedSequences     = true;
  for (u32bit z=0; z<256; z++)
    theGlobalDesc._alphabet[z] = 0;



  //  Remove any existing index -- this will reset permissions if they
  //  are hosed, and generally clean the slate for us.
  //
  errno = 0;
  unlink(_indexname);
  if (errno) {
    if (errno != ENOENT)
      fprintf(stderr, "WARNING:  Can't remove old index '%s' -- index being built anyway\n%s\n", _indexname, strerror(errno));
  }



  //  Before we compute the index, try to open the output file.  This
  //  saves us the compute if we cannot actually open the file
  //  (read-only directory, usually).
  //
  errno = 0;
  int indexfile = open(_indexname,
                       O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                       S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
  if (errno) {
    fprintf(stderr, "FastA::buildIndex() can't open index file '%s' for writing.\n%s\n",
            _indexname, strerror(errno));
    exit(1);
  }



  //  Stat the file to get size and times.  Also checks if the file exists.
  //
  struct stat  fastastat;
  errno = 0;
  if (stat(_filename, &fastastat)) {
    fprintf(stderr, "FastA::buildIndex()-- stat() of file '%s' failed.\n%s\n",
            _filename, strerror(errno));
    exit(1);
  }
  theGlobalDesc._fastaFileSize         = fastastat.st_size;
  theGlobalDesc._fastaModificationTime = fastastat.st_mtime;
  theGlobalDesc._fastaCreationTime     = fastastat.st_ctime;




  //  Describes the sequences themselves
  //
  u32bit                     numSeqs = 0;
  u32bit                     maxSeqs = 1024;
  _idxfa_desc               *theSeqs = new _idxfa_desc [maxSeqs];

  //  Copies of the deflines or names
  //
  u32bit                     theNamesLen = 0;
  u32bit                     theNamesMax = 16 * 1024;
  char                      *theNames    = new char [theNamesMax];

  //  Checksums
  //
  md5_increment_s           *curmd5  = 0L;
  md5_s                     *theMD5s = new md5_s [maxSeqs];

  //  Info about the sequence we just read
  //
  _idxfa_len                 defLen    = 0;
  _idxfa_pos                 defStart  = 0;
  _idxfa_len                 seqLen    = 0;
  _idxfa_pos                 seqStart  = 0;

  //  Copy of the defline for the sequence we just read
  //
  u32bit                     theDefLineLen    = 0;
  u32bit                     theDefLineMax    = 2 * 1024;
  char                      *theDefLine       = new char [theDefLineMax];


  //  Open a buffered input sequence
  //
  readBuffer                 B(_filename);


  //  We need to remember the length of the first sequence line
  //  so we can test if all sequence lines are the same length.
  //  This could be extended to be per-sequence, but that makes
  //  the index 25% larger, and doesn't (IMHO) gain much.
  //
  theGlobalDesc._seqlineLength = ~((_idxfa_len)0);


  while (!B.eof()) {

    //  Skip any whitespace before the defline
    //
    while ((!B.eof()) && isspace(B.get()))
      B.next();


    //  We should be at a '>' character now.  Fail if not.
    //
    if (B.get() != '>') {
      fprintf(stderr, "FastA::buildIndex()-- In file %s, expected '>' at beginning of defline, got '%c' instead.\n",
              _filename, B.get());
      fprintf(stderr, "FastA::buildIndex()-- File position is %ld, sequence number %u\n", B.tell(), numSeqs);

      exit(1);
    }


    //  Save the start of the defline
    //
    defStart = B.tell();
    defLen   = 0;

    //  Count the length of the defline -- counts until the first \r or \n
    //
    //  Also saves a copy of the characters if 'saveDefLine' is enabled.  This
    //  duals as a method of saving just the first word.
    //
    //  theDefLine is NOT zero terminated!
    //
    bool saveDefLine = true;
    if ((theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_ONLY)
      saveDefLine = false;

    theDefLineLen = 0;

    while ((!B.eof()) && (B.get() != '\r') && (B.get() != '\n')) {

      if (((theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_IDS) && (isspace(B.get())))
        saveDefLine = false;

      if (saveDefLine) {
        if (theDefLineLen > theDefLineMax) {
          theDefLineMax *= 2;
          char *newdef = new char [theDefLineMax];
          memcpy(newdef, theDefLine, theDefLineLen);
          delete [] theDefLine;
          theDefLine = newdef;
        }

        theDefLine[theDefLineLen++] = B.get();
      }

      defLen++;

      B.next();
    }



    //  Skip any whitespace between the defline and the start of the sequence
    //
    while ((!B.eof()) && isspace(B.get()))
      B.next();


    //  Save the start of the sequence
    //
    seqStart = B.tell();
    seqLen   = 0;


    //  We'll compare thisLineLen to firstLineLen.  If they differ,
    //  then the global _fixedWidth flag is cleared.
    //
    _idxfa_len                 thisLineLen      = 0;      //  length of the line we've read in, including whitespace
    _idxfa_len                 thisLineSpaces   = 0;      //  How many spaces before hitting end of line?

    bool                       lastLineDisagree = false;
    bool                       multipleLines    = false;

    while ((!B.eof()) &&
           (B.get() != '>')) {

      if (!isspace(B.get())) {
        seqLen++;
        thisLineLen++;
        theGlobalDesc._alphabet[(int)B.get()] = 1;

        //  Add this character to the MD5 hash
        //
        if (theGlobalDesc._indexType & FASTA_INDEX_MD5)
          curmd5 = md5_increment_char(curmd5, B.get());

        //  If we've seen space already, then we have embedded space,
        //  and we're not fixed width or squeezed.
        //
        if (thisLineSpaces) {
          theGlobalDesc._fixedWidth        = false;
          theGlobalDesc._squeezedSequences = false;
        }

        if (multipleLines)
          theGlobalDesc._squeezedSequences = false;

        B.next();

      } else {

        //  Count the number of spaces we see.  If we get more
        //  sequence before the end of line, that's bad.  See above.
        //
        thisLineSpaces++;

        //  If our space is a newline, then set the multipleLines
        //  flag.  If we see more sequence in this sequence, then we
        //  are not squeezed.
        //
        if ((B.get() == '\r') ||
            (B.get() == '\n')) {
          multipleLines = true;

          //  First, count any remaining whitespace on this "line"
          //  (extra linebreak characters, etc).  This is used much
          //  the same as seqlineLength.
          //
          B.next();
          while ((!B.eof()) && (isspace(B.get()))) {
            thisLineSpaces++;
            B.next();
          }

          //  Check the line lengths.
          //
          //  0) if the length of the first line in the file isn't
          //  set, set it.
          //
          //  1) if the *last* line we read in was of a different
          //  length than the first, then set the non-fixed width flag
          //
          //  2) if this line length differs with the first line
          //  length, set the lastLineDisagree flag.  We'll update the
          //  global flag if this turns out to be not the last line.
          //
          if (theGlobalDesc._seqlineLength == ~((_idxfa_len)0)) {
            theGlobalDesc._seqlineLength = thisLineLen;
            theGlobalDesc._seqendlLength = thisLineSpaces;
          }

          if (lastLineDisagree) {
            theGlobalDesc._fixedWidth = false;
          }

          if ((thisLineLen != theGlobalDesc._seqlineLength) ||
              (thisLineSpaces != theGlobalDesc._seqendlLength)) {
            //fprintf(stderr, "thisLineLen=%d / %d   thisLineSpaces=%d / %d\n",
            //        thisLineLen, theGlobalDesc._seqlineLength,
            //        thisLineSpaces, theGlobalDesc._seqendlLength);
            lastLineDisagree = true;
          }


          //  Reset this line length to zero; we're done with it.
          //
          thisLineLen      = 0;
          thisLineSpaces   = 0;
        } else {
          B.next();
        }
      }
    }

    //  Make the index bigger, if needed
    //
    if (numSeqs >= maxSeqs) {
      if (maxSeqs == 0)
        maxSeqs = 16;
      maxSeqs *= 2;

      _idxfa_desc *sa = new _idxfa_desc [ maxSeqs ];
      memcpy(sa, theSeqs, sizeof(_idxfa_desc) * numSeqs);
      delete [] theSeqs;
      theSeqs = sa;

      md5_s *sm = new md5_s [ maxSeqs ];
      memcpy(sm, theMD5s, sizeof(md5_s) * numSeqs);
      delete [] theMD5s;
      theMD5s = sm;
    }

    //  Add the new sequence description to the list.
    //
    theSeqs[numSeqs]._headerStart  = defStart;
    theSeqs[numSeqs]._headerLen    = defLen;
    theSeqs[numSeqs]._seqStart     = seqStart;
    theSeqs[numSeqs]._seqLen       = seqLen;

    //  Add the md5 checksum to the list
    //
    if (theGlobalDesc._indexType & FASTA_INDEX_MD5) {
      md5_increment_finalize(curmd5);
      theMD5s[numSeqs].a = curmd5->a;
      theMD5s[numSeqs].b = curmd5->b;
      md5_increment_destroy(curmd5);
      curmd5 = 0L;
    }

    //  Add the description of the sequence to the list of descriptions
    //
    if (((theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_IDS) ||
        ((theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_DEFLINES)) {
      if (theNamesLen + defLen >= theNamesMax) {
        theNamesMax *= 2;
        char *nd = new char [theNamesMax];
        memcpy(nd, theNames, sizeof(char) * theNamesLen);
        delete [] theNames;
        theNames = nd;
      }
      memcpy(theNames + theNamesLen, theDefLine, theDefLineLen);
      theNamesLen += theDefLineLen;
      theNames[theNamesLen++] = 0;
    }

    numSeqs++;
  }


  //  All done
  //
  theGlobalDesc._numberOfSequences = numSeqs;


  //  XXX:
  //
  //  If the sequences are all on one line (e.g., squeezed) then the
  //  method used to check for lines of the same length fails.
  //  Ideally, we should now check all the _seqLen's in the index to
  //  decide if they are really fixed length or not......but why?
  //
  if (theGlobalDesc._squeezedSequences)
    theGlobalDesc._fixedWidth = false;






  ///////////////////////////////////////
  //
  //  Write the data file; version, number of sequences and sequence
  //  descriptions.
  //
  errno = 0;
  write(indexfile, &theGlobalDesc, sizeof(_idxfa_global));
  if (errno) {
    fprintf(stderr, "FastA::buildIndex() can't write header to index file '%s'.\n%s\n", _filename, strerror(errno));
    exit(1);
  }

  errno = 0;
  write(indexfile, theSeqs, sizeof(_idxfa_desc) * numSeqs);
  if (errno) {
    fprintf(stderr, "FastA::buildIndex() can't write index to index file '%s'.\n%s\n", _filename, strerror(errno));
    exit(1);
  }

  if (((theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_IDS) ||
      ((theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_DEFLINES)) {
    errno = 0;
    write(indexfile, &theNamesLen,  sizeof(u32bit));
    if (errno) {
      fprintf(stderr, "FastA::buildIndex() can't write nameslen to index file '%s'.\n%s\n", _filename, strerror(errno));
      exit(1);
    }

    errno = 0;
    write(indexfile, theNames,  sizeof(char) * theNamesLen);
    if (errno) {
      fprintf(stderr, "FastA::buildIndex() can't write names to index file '%s'.\n%s\n", _filename, strerror(errno));
      exit(1);
    }
  }

  if (theGlobalDesc._indexType & FASTA_INDEX_MD5) {
    errno = 0;
    write(indexfile, theMD5s, sizeof(md5_s) * numSeqs);
    if (errno) {
      fprintf(stderr, "FastA::buildIndex() can't write checksums to index file '%s'.\n%s\n", _filename, strerror(errno));
      exit(1);
    }
  }


  //  Close the index file.
  //
  errno = 0;
  close(indexfile);
  if (errno) {
    fprintf(stderr, "FastA::buildIndex() can't close index file '%s'.\n%s\n", _filename, strerror(errno));
    exit(1);
  }


  delete [] theSeqs;
  delete [] theNames;
  delete [] theDefLine;
}
