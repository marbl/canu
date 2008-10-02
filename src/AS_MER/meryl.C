
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2003-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2004,      Brian Walenz
 * Copyright (C) 2005-2007, J. Craig Venter Institute
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

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bio++.H"
#include "meryl.H"

#include "AS_MER_gkpStore_to_FastABase.H"

//
//  This is a SHADOW COPY!  The main really exists in kmer/meryl/meryl.C!
//

int
main(int argc, char **argv) {
  merylArgs   *args = new merylArgs(argc, argv);

  //seqFactory::instance()->registerFile(new gkpStoreFile());
  gkpStoreFile::registerFile();

  switch (args->personality) {
    case 'P':
      estimate(args);
      break;

    case 'B':
      build(args);
      break;

    case 't':
      dumpThreshold(args);
      break;
    case 'c':
      countUnique(args);
      break;
#ifdef PLOT_DISTANCE
      //  Can't do this with big mers
    case 'p':
      plotDistanceBetweenMers(args);
      break;
#endif
    case 'h':
      plotHistogram(args);
      break;

    case PERSONALITY_MIN:
    case PERSONALITY_MINEXIST:
    case PERSONALITY_MAX:
    case PERSONALITY_ADD:
    case PERSONALITY_AND:
    case PERSONALITY_NAND:
    case PERSONALITY_OR:
    case PERSONALITY_XOR:
      multipleOperations(args);
      break;

    case PERSONALITY_SUB:
    case PERSONALITY_ABS:
      binaryOperations(args);
      break;

    case PERSONALITY_LEQ:
    case PERSONALITY_GEQ:
    case PERSONALITY_EQ:
      unaryOperations(args);
      break;

    default:
      args->usage();
      fprintf(stderr, "%s: unknown personality.  Specify -P, -B, -S or -M!\n", args->execName);
      exit(1);
      break;
  }

  delete args;

#ifdef MEMORY_DEBUG
  _dump_allocated_delta(fileno(stdout));
#endif

  return(0);
}
