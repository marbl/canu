
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
 *    Brian P. Walenz from 2008-JUL-24 to 2014-APR-11
 *      are Copyright 2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <unistd.h>
#include <time.h>
#include <math.h>

#include "util++.H"

struct header_s {
  uint64  t1;
  char    s1[570];
  uint64  t2;
};

struct record_s {
  uint64   t1;
  char     s1[123];
};

int
main(int argc, char **argv) {
  header_s      h;
  record_s      r;

  h.t1 = 0x0123456789abcdefllu;
  memset(h.s1, 0x66, 570);
  strcpy(h.s1, "this is the header");
  h.t2 = 0xdeadbeefdeadbeefllu;

  recordFile   *RF = new recordFile("test", sizeof(header_s), sizeof(record_s), 'w');

  memcpy(RF->header(), &h, sizeof(header_s));

  r.t1 = 1;   memset(r.s1, 0x66, 123);  strcpy(r.s1, "record1");
  RF->putRecord(&r);

  r.t1 = 2;   memset(r.s1, 0x66, 123);  strcpy(r.s1, "record2");
  RF->putRecord(&r);

  r.t1 = 3;   memset(r.s1, 0x66, 123);  strcpy(r.s1, "record3");
  RF->putRecord(&r);

  r.t1 = 4;   memset(r.s1, 0x66, 123);  strcpy(r.s1, "record4");
  RF->putRecord(&r);

  r.t1 = 5;   memset(r.s1, 0x66, 123);  strcpy(r.s1, "record5");
  RF->putRecord(&r);

  delete RF;

  RF = new recordFile("test", sizeof(header_s), sizeof(record_s), 'r');

  header_s *hh = (header_s *)RF->header();

  fprintf(stderr, "header t1 "uint64HEX" '%s' t2 "uint64HEX"\n", hh->t1, hh->s1, hh->t2);
  RF->getRecord(&r);  fprintf(stderr, "record "uint64FMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "uint64FMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "uint64FMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "uint64FMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "uint64FMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "uint64FMT" '%s'\n", r.t1, r.s1);

  delete RF;

  return(0);
}
