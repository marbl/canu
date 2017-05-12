
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
 *    Brian P. Walenz beginning on 2017-APR-04
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_UTL_fileIO.H"

#include "gfa.H"



template<typename TT>
static
bool
findGFAtokenI(char *features, char *token, TT &value) {
  char *p = NULL;

  p = strstr(features, token);

  if (p == NULL)
    return(false);

  p += strlen(token);  //  Skip over the token...

  if (*p == ':')       //  ...and any :, if the user forgot to include it.
    p++;

  value = (TT)strtoll(p, NULL, 10);

  //fprintf(stderr, "FOUND feature '%s' in '%s' -> '%s' %u\n", token, features, p, value);

  return(true);
}


//  Search for canu-specific names, and convert to tigID's.
//    Allow either 'tig', 'utg' or 'ctg'.
static
uint32
nameToCanuID(char *name) {
  uint32   id = UINT32_MAX;

  if (((name[0] == 't') && (name[1] == 'i') && (name[2] == 'g')) ||
      ((name[0] == 'u') && (name[1] == 't') && (name[2] == 'g')) ||
      ((name[0] == 'c') && (name[1] == 't') && (name[2] == 'g')))
    id = strtoll(name + 3, NULL, 10);

  return(id);
}



gfaSequence::gfaSequence() {
  _name     = NULL;
  _id       = UINT32_MAX;
  _sequence = NULL;
  _features = NULL;
  _length   = 0;
}


gfaSequence::gfaSequence(char *inLine) {
  load(inLine);
}


gfaSequence::gfaSequence(char *name, uint32 id, uint32 len) {
  _name     = new char [strlen(name) + 1];
  _id       = id;
  _sequence = NULL;
  _features = NULL;
  _length   = len;

  strcpy(_name, name);
}


gfaSequence::~gfaSequence() {
  delete [] _name;
  delete [] _sequence;
  delete [] _features;
}


void
gfaSequence::load(char *inLine) {
  splitToWords W(inLine);

  _name     = new char [strlen(W[1]) + 1];
  _id       = UINT32_MAX;
  _sequence = new char [strlen(W[2]) + 1];
  _features = new char [strlen(W[3]) + 1];

  _length   = 0;

  strcpy(_name,     W[1]);
  strcpy(_sequence, W[2]);
  strcpy(_features, W[3]);

  //  Scan the _features for a length.

  findGFAtokenI(_features, "LN:i:", _length);

  //  And any canu ID

  _id = nameToCanuID(_name);
}


void
gfaSequence::save(FILE *outFile) {
  fprintf(outFile, "S\t%s\t%s\tLN:i:%u\n",
          _name,
          _sequence ? _sequence : "*",
          _length);
}


gfaLink::gfaLink() {
  _Aname    = NULL;
  _Aid      = UINT32_MAX;
  _Afwd     = false;

  _Bname    = NULL;
  _Bid      = UINT32_MAX;
  _Bfwd     = false;

  _cigar    = NULL;
  _features = NULL;
}


gfaLink::gfaLink(char *inLine) {
  load(inLine);
}


gfaLink::gfaLink(char *Aname, uint32 Aid, bool Afwd,
                 char *Bname, uint32 Bid, bool Bfwd, char *cigar) {
  _Aname    = new char [strlen(Aname) + 1];
  _Aid      = Aid;
  _Afwd     = Afwd;

  _Bname    = new char [strlen(Bname) + 1];
  _Bid      = Bid;
  _Bfwd     = Bfwd;

  _cigar    = new char [strlen(cigar) + 1];
  _features = NULL;

  strcpy(_Aname,    Aname);
  strcpy(_Bname,    Bname);
  strcpy(_cigar,    cigar);

  _Aid = nameToCanuID(_Aname);    //  Search for canu-specific names, and convert to tigID's.
  _Bid = nameToCanuID(_Bname);
}


gfaLink::~gfaLink() {
  delete [] _Aname;
  delete [] _Bname;
  delete [] _cigar;
  delete [] _features;
}


void
gfaLink::load(char *inLine) {
  splitToWords W(inLine);

  _Aname    = new char [strlen(W[1]) + 1];
  _Aid      = UINT32_MAX;
  _Afwd     = W[2][0] == '+';

  _Bname    = new char [strlen(W[3]) + 1];
  _Bid      = UINT32_MAX;
  _Bfwd     = W[4][0] == '+';

  _cigar    = new char [strlen(W[5]) + 1];

  _features = new char [(W[6]) ? strlen(W[6]) + 1 : 1];

  strcpy(_Aname,    W[1]);
  strcpy(_Bname,    W[3]);
  strcpy(_cigar,    W[5]);
  strcpy(_features, (W[6]) ? W[6] : "");

  _Aid = nameToCanuID(_Aname);    //  Search for canu-specific names, and convert to tigID's.
  _Bid = nameToCanuID(_Bname);
}


void
gfaLink::save(FILE *outFile) {
  fprintf(outFile, "L\t%s\t%c\t%s\t%c\t%s\n",
          _Aname, (_Afwd == true) ? '+' : '-',
          _Bname, (_Bfwd == true) ? '+' : '-',
          (_cigar == NULL) ? "*" : _cigar);
}


void
gfaLink::alignmentLength(int32 &queryLen, int32 &refceLen, int32 &alignLen) {
  char  *cp = _cigar;

  refceLen = 0;  //  Bases on the reference involved in the alignment
  queryLen = 0;  //  Bases on the query     involved in the alignment
  alignLen = 0;  //  Length of the alignment

  if (cp == NULL)
    return;

  if (*cp == '*')
    return;

  do {
    int64  val  = strtoll(cp, &cp, 10);
    char   code = *cp++;

    switch (code) {
      case 'M':  //  Alignment, either match or mismatch
        refceLen += val;
        queryLen += val;
        alignLen += val;
        break;
      case 'I':  //  Insertion to the reference - gap in query
        queryLen += val;
        alignLen += val;
        break;
      case 'D':  //  Deletion from the reference - gap in reference
        refceLen += val;
        alignLen += val;
       break;
      case 'N':  //  Skipped in the reference (e.g., intron)
        refceLen += val;
        fprintf(stderr, "warning - unsupported CIGAR code '%c' in '%s'\n", code, _cigar);
        break;
      case 'S':  //  Soft-clipped from the query - not part of the alignment
        fprintf(stderr, "warning - unsupported CIGAR code '%c' in '%s'\n", code, _cigar);
        break;
      case 'H':  //  Hard-clipped from the query - not part of the alignment, and removed from the read as input
        fprintf(stderr, "warning - unsupported CIGAR code '%c' in '%s'\n", code, _cigar);
        break;
      case 'P':  //  Padding - "silent deletion from padded reference" - ???
        fprintf(stderr, "warning - unsupported CIGAR code '%c' in '%s'\n", code, _cigar);
        break;
      case '=':  //  Alignment, match
        refceLen += val;
        queryLen += val;
        alignLen += val;
        break;
      case 'X':  //  Alignment, mismatch
        refceLen += val;
        queryLen += val;
        alignLen += val;
        break;
      default:
        fprintf(stderr, "unknown CIGAR code '%c' in '%s'\n", code, _cigar);
        break;
    }

  } while (*cp != 0);
}




gfaFile::gfaFile() {
  _header = NULL;
}


gfaFile::gfaFile(char *inFile) {
  _header = NULL;

  if ((inFile[0] == 'H') && (inFile[1] == '\t')) {
    _header = new char [strlen(inFile) + 1];
    strcpy(_header, inFile);
  }

  else {
    loadFile(inFile);
  }
}


gfaFile::~gfaFile() {
  delete [] _header;

  for (uint32 ii=0; ii<_sequences.size(); ii++)
    delete _sequences[ii];

  for (uint32 ii=0; ii<_links.size(); ii++)
    delete _links[ii];
}


bool
gfaFile::loadFile(char *inFile) {
  FILE  *F    = NULL;
  char  *L    = NULL;
  uint32 Llen = 0;
  uint32 Lmax = 0;

  errno = 0;
  F = fopen(inFile, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", inFile, strerror(errno)), exit(1);


  while (AS_UTL_readLine(L, Llen, Lmax, F)) {
    char  type = L[0];

    if (L[1] != '\t')
      fprintf(stderr, "gfaFile::loadFile()-- misformed file; second letter must be tab in line '%s'\n", L), exit(1);

    if      (type == 'H') {
      delete [] _header;
      _header = new char [Llen];
      strcpy(_header, L+2);
    }

    else if (type == 'S') {
      _sequences.push_back(new gfaSequence(L));
    }

    else if (type == 'L') {
      _links.push_back(new gfaLink(L));
    }

    else {
      fprintf(stderr, "gfaFile::loadFile()-- unrecognized line '%s'\n", L), exit(1);
    }
  }

  fclose(F);

  delete [] L;

  fprintf(stderr, "gfa:  Loaded " F_S64 " sequences and " F_S64 " links.\n", _sequences.size(), _links.size());

  return(true);
}




bool
gfaFile::saveFile(char *outFile) {
  FILE  *F = NULL;

  errno = 0;
  F = fopen(outFile, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", outFile, strerror(errno)), exit(1);

  fprintf(F, "H\t%s\n", _header);

  for (uint32 ii=0; ii<_sequences.size(); ii++)
    if (_sequences[ii])
      _sequences[ii]->save(F);

  for (uint32 ii=0; ii<_links.size(); ii++)
    if (_links[ii])
      _links[ii]->save(F);

  fclose(F);

  return(true);
}

