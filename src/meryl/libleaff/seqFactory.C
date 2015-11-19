
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
 *    Brian P. Walenz on 2015-FEB-04
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "seqFactory.H"

#include "fastaFile.H"
#include "fastaStdin.H"
#include "fastqFile.H"
#include "fastqStdin.H"
#include "seqStore.H"
#include "gkStoreFile.H"

seqFactory *seqFactory::me = 0L;


seqFactory::seqFactory() {
  _filesNum = 0;
  _filesMax = 16;
  _files = new seqFile * [_filesMax];

  registerFile(new fastaFile);
  registerFile(new fastaStdin);
  registerFile(new fastqFile);
  registerFile(new fastqStdin);
  registerFile(new seqStore);
  registerFile(new gkStoreFile);
  //registerFile(new sffFile);
}


seqFactory::~seqFactory() {
  for (uint32 i=0; i<_filesNum; i++)
    delete _files[i];
  delete [] _files;
}


void
seqFactory::registerFile(seqFile *f) {
  if (_filesNum >= _filesMax) {
    fprintf(stderr, "seqFactory::registerFile()--  Wow!  You registered lots of files!  Now fix %s at line %d.\n", __FILE__, __LINE__);
    exit(1);
  }
  _files[_filesNum++] = f;
}


seqFile *
seqFactory::openFile(const char *name) {
  seqFile  *n = 0L;

  for (uint32 i=0; i<_filesNum; i++) {
    n = _files[i]->openFile(name);
    if (n)
      return(n);
  }

  fprintf(stderr, "seqFactory::registerFile()--  Cannot determine type of file '%s'.  Tried:\n", name);

  for (uint32 i=0; i<_filesNum; i++)
    fprintf(stderr, "seqFactory::registerFile()--         '%s'\n", _files[i]->getFileTypeName());

  exit(1);
  return(n);
}
