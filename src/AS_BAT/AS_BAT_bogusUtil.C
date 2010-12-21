
/**************************************************************************
 * Copyright (C) 2010, J Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id: AS_BAT_bogusUtil.C,v 1.8 2010-12-21 19:45:18 brianwalenz Exp $";

#include "AS_BAT_bogusUtil.H"

bool
byFragmentID(const genomeAlignment &A, const genomeAlignment &B) {
  return(A.frgIID < B.frgIID);
}

bool
byGenomePosition(const genomeAlignment &A, const genomeAlignment &B) {
  if (A.genBgn < B.genBgn)
    //  A clearly before B.
    return(true);

  if (A.genBgn > B.genBgn)
    //  A clearly after B.
    return(false);

  if ((A.isRepeat == false) && (B.isRepeat == true))
    //  Start at the same spot, put unique stuff first
    return(true);

  //  Note the second condition above.  If both A and B are repeat=false, byGenomePosition(A,B) and
  //  byGenomePosition(B,A) bth return true, logically impossible.

  return(false);
}




void
addAlignment(vector<genomeAlignment>   &genome,
             int32  frgIID, int32  frgLen,
             int32  frgBgn, int32  frgEnd, bool  isReverse,
             int32  genBgn, int32  genEnd) {
  genomeAlignment  A;

  A.frgIID    = frgIID;
  A.frgLen    = frgLen;
  A.frgBgn    = frgBgn;
  A.frgEnd    = frgEnd;
  A.genBgn    = genBgn;
  A.genEnd    = genEnd;
  A.isReverse = isReverse;
  A.isSpanned = false;
  A.isRepeat  = true;

  assert(A.frgBgn < A.frgEnd);
  assert(A.genBgn < A.genEnd);

  genome.push_back(A);
}


void
loadNucmer(char                      *nucmerName,
           vector<genomeAlignment>   &genome,
           map<string, int32>        &IIDmap,
           vector<string>            &IIDname,
           vector<uint32>            &IIDcount) {
  FILE  *inFile = 0L;
  char   inLine[1024];

  errno = 0;
  inFile = fopen(nucmerName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", nucmerName, strerror(errno)), exit(1);

  //  First FIVE lines is a header
  fgets(inLine, 1024, inFile);
  fgets(inLine, 1024, inFile);
  fgets(inLine, 1024, inFile);
  fgets(inLine, 1024, inFile);
  fgets(inLine, 1024, inFile);

  //  Read the first line.
  fgets(inLine, 1024, inFile);
  chomp(inLine);

  while (!feof(inFile)) {
    splitToWords     W(inLine);
    genomeAlignment  A;
    string           ID = W[12];

    if (IIDmap.find(ID) == IIDmap.end()) {
      IIDname.push_back(ID);
      IIDcount.push_back(0);
      IIDmap[ID] = IIDname.size() - 1;
    }

    //  Unlike snapper, these are already in base-based coords.

    A.frgIID    = IIDmap[ID];
    A.frgLen    = 0;
    A.frgBgn    = W(3);
    A.frgEnd    = W(4);
    A.genBgn    = W(0);
    A.genEnd    = W(1);
    A.isReverse = false;
    A.isSpanned = false;
    A.isRepeat  = true;

    if (A.frgBgn > A.frgEnd) {
      A.frgBgn    = W(4);
      A.frgEnd    = W(3);
      A.isReverse = true;
    }

    IIDcount[A.frgIID]++;

    assert(A.frgBgn < A.frgEnd);
    assert(A.genBgn < A.genEnd);

    genome.push_back(A);

    fgets(inLine, 1024, inFile);
    chomp(inLine);
  }

  fclose(inFile);
}



void
loadSnapper(char                      *snapperName,
            vector<genomeAlignment>   &genome,
            map<string, int32>        &IIDmap,
            vector<string>            &IIDname,
            vector<uint32>            &IIDcount) {
  FILE  *inFile = 0L;
  char   inLine[1024];

  errno = 0;
  inFile = fopen(snapperName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", snapperName, strerror(errno)), exit(1);

  //  Read the first line
  fgets(inLine, 1024, inFile);
  chomp(inLine);

  while (!feof(inFile)) {
    splitToWords     W(inLine);
    genomeAlignment  A;
    string           ID = W[0];

    if (strncmp(inLine, "cDNAid", 6) == 0)
      //  Skip header lines.
      goto nextLine;

    if (IIDmap.find(ID) == IIDmap.end()) {
      IIDname.push_back(ID);
      IIDcount.push_back(0);
      IIDmap[ID] = IIDname.size() - 1;
    }

    //  "+1" -- Convert from space-based coords to base-based coords.

    A.frgIID    = IIDmap[ID];
    A.frgLen    = 0;
    A.frgBgn    = W(3) + 1;
    A.frgEnd    = W(4);
    A.genBgn    = W(6) + 1;
    A.genEnd    = W(7);
    A.isReverse = false;
    A.isSpanned = false;
    A.isRepeat  = true;

    if (A.frgBgn > A.frgEnd) {
      A.frgBgn    = W(4) + 1;
      A.frgEnd    = W(3);
      A.isReverse = true;
    }

    IIDcount[A.frgIID]++;

    assert(A.frgBgn < A.frgEnd);
    assert(A.genBgn < A.genEnd);

    genome.push_back(A);

  nextLine:
    fgets(inLine, 1024, inFile);
    chomp(inLine);
  }

  fclose(inFile);
}



void
loadReferenceSequence(char     *refName,
                      char    *&refhdr,
                      char    *&refseq,
                      int32    &reflen) {

  FILE    *F      = fopen(refName, "r");

  refhdr = new char [1024];
  refseq = new char [16 * 1024 * 1024];

  fgets(refhdr,             1024, F);
  fgets(refseq, 16 * 1026 * 1024, F);

  fclose(F);

  chomp(refhdr);
  chomp(refseq);

  for (uint32 i=0; refhdr[i]; i++) {
    refhdr[i] = refhdr[i+1];
    if (isspace(refhdr[i]))
      refhdr[i] = 0;
  }

  reflen = strlen(refseq);
}

