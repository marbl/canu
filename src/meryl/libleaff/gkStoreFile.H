
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
 *    Brian P. Walenz from 2015-FEB-04 to 2015-MAR-17
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef GKSTOREFILE_H
#define GKSTOREFILE_H

#include "seqFile.H"
#include "gkStore.H"

class gkStoreFile : public seqFile {
protected:
  gkStoreFile(const char *filename);
  gkStoreFile();

public:
  ~gkStoreFile();

protected:
  seqFile      *openFile(const char *name);

public:
  uint32        find(const char *sequencename) {
    fprintf(stderr, "gkStoreFile::find()--  Lookup of sequencename '%s' not supported.\n", sequencename);
    assert(0);
    return(0);
  };

  uint32        getSequenceLength(uint32 iid) {
    return(gkp->gkStore_getRead(iid + 1)->gkRead_sequenceLength());
  };

  bool          getSequence(uint32 iid,
                            char *&h, uint32 &hLen, uint32 &hMax,
                            char *&s, uint32 &sLen, uint32 &sMax);
  bool          getSequence(uint32 iid,
                            uint32 bgn, uint32 end, char *s);

private:
  void          clear(void);

  gkStore      *gkp;
  gkReadData    readData;

  friend class seqFactory;
};


#endif  //  GKSTOREFILE_H
