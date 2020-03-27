
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef __GCON_ALIGNMENT_HPP__
#define __GCON_ALIGNMENT_HPP__

///
/// Super-simple alignment representation.  Represents an alignment between two
/// PacBio reads, one of which we're trying to correct.  The read to correct
/// may be either the target or the query, depending on how the alignment was
/// done.
///

#include <string>
#include <stdint.h>

class dagAlignment {
public:
  dagAlignment() {
    start   = 0;
    end     = 0;

    length  = 0;

    qstr    = NULL;
    tstr    = NULL;
  };
  ~dagAlignment() {
    delete [] qstr;
    delete [] tstr;
  };

  void   clear(void) {

    delete [] qstr;
    delete [] tstr;

    start  = 0;
    end    = 0;

    length = 0;

    qstr   = NULL;
    tstr   = NULL;
  };

  void     normalizeGaps(void);

  uint32_t start;   //   1-based!
  uint32_t end;

  uint32_t length;

  char *qstr;
  char *tstr;
};


#endif // __GCON_ALIGNMENT_HPP__
