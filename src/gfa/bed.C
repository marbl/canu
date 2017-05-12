
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

#include "bed.H"



//  Search for canu-specific names, and convert to tigID's.
static
uint32
nameToCanuID(char *name) {
  uint32   id = UINT32_MAX;

  if ((name[0] == 't') &&
      (name[1] == 'i') &&
      (name[2] == 'g'))
    id = strtoll(name + 3, NULL, 10);

  if ((name[0] == 'c') &&
      (name[1] == 't') &&
      (name[2] == 'g'))
    id = strtoll(name + 3, NULL, 10);

  if ((name[0] == 'u') &&
      (name[1] == 't') &&
      (name[2] == 'g'))
    id = strtoll(name + 3, NULL, 10);

  return(id);
}



bedRecord::bedRecord() {
  _Aname    = NULL;
  _Aid      = UINT32_MAX;

  _bgn      = UINT32_MAX;
  _end      = 0;

  _Bname    = NULL;
  _Bid      = UINT32_MAX;

  _score    = 0;
  _Bfwd     = false;
}


bedRecord::bedRecord(char *inLine) {
  load(inLine);
}


bedRecord::~bedRecord() {
  delete [] _Aname;
  delete [] _Bname;
}


void
bedRecord::load(char *inLine) {
  splitToWords W(inLine);

  _Aname    = new char [strlen(W[0]) + 1];
  _Aid      = UINT32_MAX;

  _bgn      = W(1);
  _end      = W(2);

  _Bname    = new char [strlen(W[3]) + 1];
  _Bid      = UINT32_MAX;

  _score    = W(4);
  _Bfwd     = W[5][0] == '+';

  strcpy(_Aname,    W[0]);
  strcpy(_Bname,    W[3]);

  _Aid = nameToCanuID(_Aname);    //  Search for canu-specific names, and convert to tigID's.
  _Bid = nameToCanuID(_Bname);
}


void
bedRecord::save(FILE *outFile) {
  fprintf(outFile, "%s\t%d\t%d\t%s\t%u\t%c\n",
          _Aname, _bgn, _end, _Bname, _score, (_Bfwd == true) ? '+' : '-');
}



bedFile::bedFile(char *inFile) {
  loadFile(inFile);
}


bedFile::~bedFile() {
  for (uint32 ii=0; ii<_records.size(); ii++)
    delete _records[ii];
}


bool
bedFile::loadFile(char *inFile) {
  FILE  *F    = NULL;
  char  *L    = NULL;
  uint32 Llen = 0;
  uint32 Lmax = 0;

  errno = 0;
  F = fopen(inFile, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", inFile, strerror(errno)), exit(1);

  while (AS_UTL_readLine(L, Llen, Lmax, F)) {
    _records.push_back(new bedRecord(L));
  }

  fclose(F);

  delete [] L;

  fprintf(stderr, "bed:  Loaded " F_S64 " records.\n", _records.size());

  return(true);
}




bool
bedFile::saveFile(char *outFile) {
  FILE  *F = NULL;

  errno = 0;
  F = fopen(outFile, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", outFile, strerror(errno)), exit(1);

  for (uint32 ii=0; ii<_records.size(); ii++)
    if (_records[ii])
      _records[ii]->save(F);

  fclose(F);

  return(true);
}

