#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

#include "bio++.H"

//  If the index exists, check the global information in the index
//  against the version of this code, and the fasta file it is
//  indexing.  If everything checks out OK, just return.  Otherwise,
//  continue and build a new index.
//
//  Returns false if
//    The index doesn't exist
//    The index exists, but is empty
//    The index doesn't match the fasta file and would need to be rebuilt
//
bool
FastAWrapper::isIndexOnDiskCompatible(u32bit indextype, const char *indexname, bool beVerbose) {
  errno = 0;

  // Construct an index filename if one wasn't given.
  //
  if (indexname) {
    _indexname = new char [strlen(indexname) + 1];
    strcpy(_indexname, indexname);
  } else {
    _indexname = new char [strlen(_filename) + 4];
    strcpy(_indexname, _filename);
    strcat(_indexname, "idx");
  }


  //  Stat the index to see if it exists
  struct stat    st;
  stat(_indexname, &st);
  if (errno)
    return(false);

  if (st.st_size == 0)
    //  Index is empty
    return(false);

  //  Stat the -> fasta <- file to get the modification time
  // 
  stat(_filename, &st);
  if (errno) {
    fprintf(stderr, "ERROR: Can't stat '%s'.\n%s\n", _filename, strerror(errno));
    exit(1);
  }

  if (st.st_size == 0) {
    //  file is empty
    fprintf(stderr, "ERROR: '%s' is empty?\n", _filename);
    exit(1);
  }

  int indexfile = open(_indexname, O_RDONLY | O_LARGEFILE);
  if (errno) {
    //  Index doesn't exist??
    fprintf(stderr, "ERROR: Index file exists, but can't be opened for reading?\n%s\n", strerror(errno));
    exit(1);
  }

  //  Opened an index, so read the description.
  //
  _idxfa_global  theGlobalDesc;
  read(indexfile, &theGlobalDesc, sizeof(_idxfa_global));
  if (errno) {
    //  can't read
    fprintf(stderr, "ERROR: Index file exists, but can't read description.\n%s\n", strerror(errno));
    exit(1);
  }

  //  Close the index
  //
  errno = 0;
  close(indexfile);
  if (errno == 0) {
    //  can't close
  }

  //  Check that the magic number and version are correct and that
  //  the modification time of the fasta file is the same as stored
  //  in the index.
  //
  if ((theGlobalDesc._magic                          == FASTA_MAGICNUMBER) &&
      (theGlobalDesc._version                        <= FASTA_VERSIONNUMBER) &&
      (theGlobalDesc._fastaModificationTime          == st.st_mtime) &&
      ((theGlobalDesc._indexType & FASTA_INDEX_MASK) >= (indextype & FASTA_INDEX_MASK)) &&
      ((theGlobalDesc._indexType & FASTA_INDEX_MD5)  >= (indextype & FASTA_INDEX_MD5))) {
    return(true);
  } else {
    if (beVerbose) {
      fprintf(stderr, "WARNING: Index found, but stale or wrong type!\n");

      if (theGlobalDesc._magic != FASTA_MAGICNUMBER)
        fprintf(stderr, "           Magic number incorrect; perhaps not an index file?  got=0x%016lx expected=0x%016lx\n", theGlobalDesc._magic, FASTA_MAGICNUMBER);

      if (theGlobalDesc._version != FASTA_VERSIONNUMBER)
        fprintf(stderr, "           Version number incorrect; got %d, expected %d.\n", theGlobalDesc._version, FASTA_VERSIONNUMBER);

      if (theGlobalDesc._fastaModificationTime != st.st_mtime)
        fprintf(stderr, "           File age not correct; got %ld, expected %ld.\n", theGlobalDesc._fastaModificationTime, st.st_mtime);

      if ((theGlobalDesc._indexType & FASTA_INDEX_MASK) < (indextype & FASTA_INDEX_MASK))
        fprintf(stderr, "           Type of index insufficient; got %s, need %s.\n",
                indexTypeNames(theGlobalDesc._indexType),
                indexTypeNames(indextype));

      if ((theGlobalDesc._indexType & FASTA_INDEX_MD5) < (indextype & FASTA_INDEX_MD5))
        fprintf(stderr, "           MD5 checksums not present.\n");
    }
  }

  return(false);
}


bool
FastAWrapper::isIndexValid(u32bit indextype, const char *indexname) {
  return (isIndexOnDiskCompatible(indextype, indexname, false));
}





void
FastAWrapper::openIndex(u32bit indextype, const char *indexname) {


  //  If we've been told to open an index, but we're not random
  //  access, complain.
  //
  if (_isStreamInput) {
    fprintf(stderr, "FastAWrapper()--  '%s' is a stream and not valid for indexing!\n", _filename);
    return;
  }

  //  If the index is already open, and it's the same type as desired, we
  //  can just return.  Downgrading an index is not allowed -- if you open
  //  a FASTA_INDEX_PLUS_DEFLINES, you're stuck with it.
  //
  if ((_isRandomAccess) && (indextype & FASTA_INDEX_MASK <= _theGlobalDesc._indexType & FASTA_INDEX_MASK))
    return;


  //  If no index, or it's not current, build a new one.
  //
  if (! isIndexOnDiskCompatible(indextype, indexname, true))
    createIndex(indextype);

  //  Index exists, or we wouldn't have made it here.
  //
  _isRandomAccess = true;

  //  Open the index, read it in.
  //
  errno = 0;
  int indexfile = open(_indexname, O_RDONLY | O_LARGEFILE);
  if (errno) {
    fprintf(stderr, "FastAWrapper()-- couldn't open the index '%s'.\n%s\n", _indexname, strerror(errno));
    exit(1);
  }


  //
  //  Read the index
  //
  read(indexfile, &_theGlobalDesc, sizeof(_idxfa_global));
  if (errno) {
    fprintf(stderr, "FastAWrapper()-- couldn't read description from the index '%s'.\n%s\n", _indexname, strerror(errno));
    exit(1);
  }

 _theSeqs = new _idxfa_desc [ _theGlobalDesc._numberOfSequences ];

  read(indexfile, _theSeqs, sizeof(_idxfa_desc) * _theGlobalDesc._numberOfSequences);
  if (errno) {
    fprintf(stderr, "FastAWrapper()-- couldn't read sequence descriptions from the index '%s'.\n%s\n", _indexname, strerror(errno));
    exit(1);
  }

  //  If indextype is FASTA_INDEX_ANY, open the index as reported by
  //  the file.  Otherwise, open the index as specified.
  //
  if ((indextype & FASTA_INDEX_MASK) == FASTA_INDEX_ANY)
    indextype = _theGlobalDesc._indexType;
  else
    _theGlobalDesc._indexType = indextype;



  if (((_theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_IDS) ||
      ((_theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_DEFLINES)) {
    read(indexfile, &_theNamesLen, sizeof(u32bit));
    if (errno) {
      fprintf(stderr, "FastAWrapper()-- couldn't read lengths of names from the index '%s'.\n%s\n", _indexname, strerror(errno));
      exit(1);
    }

    _theNames = new char [sizeof(char) * _theNamesLen];

    read(indexfile, _theNames, sizeof(char) * _theNamesLen);
    if (errno) {
      fprintf(stderr, "FastAWrapper()-- couldn't read names from the index '%s'.\n%s\n", _indexname, strerror(errno));
      exit(1);
    }
  }


  if (_theGlobalDesc._indexType & FASTA_INDEX_MD5) {
    _theMD5s = new md5_s [_theGlobalDesc._numberOfSequences];

    errno = 0;
    read(indexfile, _theMD5s, sizeof(md5_s) * _theGlobalDesc._numberOfSequences);
    if (errno) {
      fprintf(stderr, "FastA::buildIndex() couldn't read checksums from the index '%s'.\n%s\n", _filename, strerror(errno));
      exit(1);
    }
  }


  errno = 0;
  close(indexfile);
  if (errno) {
    fprintf(stderr, "FastAWrapper()-- couldn't close the index '%s'.\n%s\n", _indexname, strerror(errno));
    //exit(1);
  }
}



void
FastAWrapper::printATADescription(FILE *out, char *name) {

    fprintf(out, "! format ata 1.0\n");
    fprintf(out, "S %s %s FASTA DNA ",
            name,
            _filename);

    for (u32bit i=0; i<256; i++)
      if (_theGlobalDesc._alphabet[i])
        fprintf(out, "%c", (char)i);
    
    if      (_theGlobalDesc._squeezedSequences) {
      fprintf(out, " SQUEEZED ");
    } else if (_theGlobalDesc._fixedWidth) {
      fprintf(out, " FIXED ");
    } else {
      fprintf(out, " VARIABLE ");
    }

    fprintf(out, u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            _theGlobalDesc._seqlineLength,
            _theGlobalDesc._seqendlLength,
            _theGlobalDesc._numberOfSequences);

    //
    //  XXX:  assumes these are chromosomes, which is probably wrong.
    //
    //  ch == chromosomes
    //  s  == scaffold
    //  c  == contig
    //

    const char *dumpIndex1 = "G ch %s:%u . 0 1 %s 0 "u32bitFMT" "u64bitFMT" "u32bitFMT" "u32bitFMT" "u64bitFMT" . .\n";
    const char *dumpIndex2 = "G ch %s:%u . 0 1 %s 0 "u32bitFMT" "u64bitFMT" "u32bitFMT" "u32bitFMT" "u64bitFMT" %s . %s\n";

    char *names = _theNames;
    char  seqid[16 * 1024];

    for (u32bit iid=0; iid<_theGlobalDesc._numberOfSequences; iid++) {
      switch (_theGlobalDesc._indexType & FASTA_INDEX_MASK) {
        case FASTA_INDEX_ONLY:
          fprintf(out, dumpIndex1,
                  name,
                  iid,
                  name,
                  _theSeqs[iid]._seqLen,
                  _theSeqs[iid]._seqStart,
                  iid,
                  _theSeqs[iid]._headerLen,
                  _theSeqs[iid]._headerStart);
          break;
        case FASTA_INDEX_PLUS_IDS:
        case FASTA_INDEX_PLUS_DEFLINES:
          strcpy(seqid, names+1);
          for (u32bit x=0; seqid[x]; x++)
            if (whitespaceSymbol[seqid[x]]) {
              seqid[x] = 0;
              break;
            }

          fprintf(out, dumpIndex2,
                  name,
                  iid,
                  name,
                  _theSeqs[iid]._seqLen,
                  _theSeqs[iid]._seqStart,
                  iid,
                  _theSeqs[iid]._headerLen,
                  _theSeqs[iid]._headerStart,
                  seqid,
                  names);

          names = moveToNextName(names, iid);
          break;
      }
    }

}



void
FastAWrapper::printTextDescription(FILE *out) {

  fprintf(out, "/FastAIndex magic=0x%016lx version=0x%08x\n",
          _theGlobalDesc._magic,
          _theGlobalDesc._version);

  fprintf(out, "/indexType             = %s\n", indexTypeNames(_theGlobalDesc._indexType));

  switch (_theGlobalDesc._fastaType) {
    case FASTA_UNDECLARED:
      fprintf(out, "/fastaType             = undeclared\n");
      break;
    default:
      fprintf(out, "/fastaType             = unknown\n");
      break;
  }

  fprintf(out, "/filename              = %s\n",          _filename);
  fprintf(out, "/numberOfSequences     = "u32bitFMT"\n", _theGlobalDesc._numberOfSequences);
  fprintf(out, "/fastaFileSize         = "u64bitFMT"\n", _theGlobalDesc._fastaFileSize);
  fprintf(out, "/fastaModificationTime = "u32bitFMT"\n", _theGlobalDesc._fastaModificationTime);
  fprintf(out, "/fastaCreationTime     = "u32bitFMT"\n", _theGlobalDesc._fastaCreationTime);
  fprintf(out, "/seqlineLength         = "u32bitFMT"\n", _theGlobalDesc._seqlineLength);
  fprintf(out, "/seqendlLength         = "u32bitFMT"\n", _theGlobalDesc._seqendlLength);
  fprintf(out, "/fixedWidth            = %s\n",          _theGlobalDesc._fixedWidth ? "yes" : "no");
  fprintf(out, "/squeezedSequences     = %s\n",          _theGlobalDesc._squeezedSequences ? "yes" : "no");

  fprintf(out, "/alphabet              = ");
  for (u32bit i=0; i<256; i++)
    if (_theGlobalDesc._alphabet[i])
      fprintf(out, "%c", (char)i);
  fprintf(out, "\n");

  char *names = _theNames;

  for (u32bit iid=0; iid<_theGlobalDesc._numberOfSequences; iid++) {
    switch (_theGlobalDesc._indexType & FASTA_INDEX_MASK) {
      case FASTA_INDEX_ONLY:
        fprintf(out, "I "u32bitFMT"\t"u64bitFMT"\t"u32bitFMT"\t"u64bitFMT"\n",
                _theSeqs[iid]._headerLen,
                _theSeqs[iid]._headerStart,
                _theSeqs[iid]._seqLen,
                _theSeqs[iid]._seqStart);
        break;
      case FASTA_INDEX_PLUS_IDS:
      case FASTA_INDEX_PLUS_DEFLINES:
        fprintf(out, "I "u32bitFMT"\t"u64bitFMT"\t"u32bitFMT"\t"u64bitFMT"\t%s\n",
                _theSeqs[iid]._headerLen,
                _theSeqs[iid]._headerStart,
                _theSeqs[iid]._seqLen,
                _theSeqs[iid]._seqStart,
                names);
        names = moveToNextName(names, iid);
        break;
    }
  }
}
