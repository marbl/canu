
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

#ifndef CLASSIFYMATES_H
#define CLASSIFYMATES_H

static const char *rcsid_CLASSIFYMATES_H = "$Id: classifyMates.H,v 1.2 2011-08-29 20:58:32 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"


//  At 16 bytes, this is bigger than I want.
//
//  We could get down to 8 bytes by losing 'iteration' and 'mateIID', both of which are convenient
//  for development, but not needed for production.
//
//  We could get down to 4 bytes by storing the starting read IID in the file header and counting.
//
//  With 1 billion reads, the full output will be 16GB, which is entirely reasonable for a few
//  releases.
//
class classifyMatesResult {
public:
  classifyMatesResult() {
    fragSpur     = 0;  //  No overlaps on one side of either the read or the mate
    mateSpur     = 0;  //
    fragChimer   = 0;  //  Overlaps on both sides, but one fragment looked chimeric
    mateChimer   = 0;  //
    fragJunction = 0;  //  Overlaps on both sides, but all we found were junctions
    mateJunction = 0;  //
    classified   = 0;  //  Found the mate at the expected distance, correctly oriented.
    limited      = 0;  //  Ran out of time.
    exhausted    = 0;  //  Ran out of places to look.

    innie        = 0;
    distance     = 0;

    isBFS        = 0;
    isDFS        = 0;
    isRFS        = 0;
    iteration    = 0;

    readIID      = 0;
    mateIID      = 0;
  };

public:
  void    print(FILE *out) {
    char   reason[256] = {0};

    if (fragSpur == true)
      strcat(reason, (reason[0]) ? ":a-spur" : "a-spur");
    if (mateSpur == true)
      strcat(reason, (reason[0]) ? ":b-spur" : "b-spur");

    if (fragChimer == true)
      strcat(reason, (reason[0]) ? ":a-chimer" : "a-chimer");
    if (mateChimer == true)
      strcat(reason, (reason[0]) ? ":b-chimer" : "b-chimer");

    if (fragJunction == true)
      strcat(reason, (reason[0]) ? ":a-junction" : "a-junction");
    if (mateJunction == true)
      strcat(reason, (reason[0]) ? ":b-junction" : "b-junction");

    if (limited == true)
      strcat(reason, (reason[0]) ? ":limited" : "limited");

    if (exhausted == true)
      strcat(reason, (reason[0]) ? ":exhausted" : "exhausted");

    if (reason[0] != 0)
      fprintf(out, "Path from "F_U32"/%s to "F_U32"/%s NOT FOUND (%s).\n",
              readIID, (innie == true) ? "5'3'" : "3'5'", 
              mateIID, (innie == true) ? "3'5'" : "5'3'",
              reason);
    else
      fprintf(out, "Path from "F_U32"/%s to "F_U32"/%s found at iter "F_U64" of length "F_U64".\n",
              readIID, (innie == true) ? "5'3'" : "3'5'", 
              mateIID, (innie == true) ? "3'5'" : "5'3'",
              iteration,
              distance);
  };

public:
  uint64  fragSpur     : 1;   //  Result
  uint64  mateSpur     : 1;
  uint64  fragChimer   : 1;
  uint64  mateChimer   : 1;
  uint64  fragJunction : 1;
  uint64  mateJunction : 1;
  uint64  classified   : 1;
  uint64  limited      : 1;
  uint64  exhausted    : 1;

  uint64  innie        : 1;   //  Orientation of solution.
  uint64  distance     : 22;  //  Distance of solution.

  uint64  isBFS        : 1;   //  Algorithm.
  uint64  isDFS        : 1;   //  Algorithm.
  uint64  isRFS        : 1;   //  Algorithm.
  uint64  iteration    : 29;  //  Up to 1m iterations.

  AS_IID  readIID;
  AS_IID  mateIID;
};



class classifyMatesResultFile {
  struct classifyMatesResultFileHeader {
    uint32  magicNumber;

    uint32  isBFS;
    uint32  isRFS;
    uint32  isDFS;

    uint32  iterationLimit;

    uint32  isInnie;
    uint32  minDistance;
    uint32  maxDistance;
  };

public:
  classifyMatesResultFile(char   *filename,
                          char    algorithm,
                          uint32  iterationLimit_,
                          bool    innie_,
                          uint32  minDistance_,
                          uint32  maxDistance_) {
    header.magicNumber    = 0x1a2b0001;
    header.isBFS          = (algorithm == 'b');
    header.isDFS          = (algorithm == 'd');
    header.isRFS          = (algorithm == 'r');
    header.iterationLimit = iterationLimit_;
    header.isInnie        = innie_;
    header.minDistance    = innie_;
    header.maxDistance    = innie_;

    errno = 0;
    F = fopen(filename, "w");
    if (errno)
      fprintf(stderr, "ERROR:  failed to open '%s': %s\n", filename, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, &header, "classifyMatesResultFileHeader", sizeof(classifyMatesResultFileHeader), 1);
  };

  classifyMatesResultFile(char *filename) {
    errno = 0;
    F = fopen(filename, "r");
    if (errno)
      fprintf(stderr, "ERROR:  failed to open '%s': %s\n", filename, strerror(errno)), exit(1);

    AS_UTL_safeRead(F, &header, "classifyMatesResultFileHeader", sizeof(classifyMatesResultFileHeader), 1);

    if (header.magicNumber != 0x1a2b0001)
      fprintf(stderr, "ERROR:  invalid magic number 0x%08x: probably not a classifyMates result file.\n", header.magicNumber), exit(1);
  };

  ~classifyMatesResultFile() {
    fclose(F);
  };


  void  write(classifyMatesResult &result) {
    AS_UTL_safeWrite(F, &result, "classifyMatesResult", sizeof(classifyMatesResult), 1);
  };

  bool  read(classifyMatesResult &result) {
    if (feof(F))
      return(false);
    AS_UTL_safeRead(F, &result, "classifyMatesResult", sizeof(classifyMatesResult), 1);
    return(!feof(F));
  };

private:
  classifyMatesResultFileHeader   header;
  FILE                           *F;
};


#endif  //  CLASSIFYMATES_H
