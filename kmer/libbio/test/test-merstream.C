#include <stdio.h>
#include <stdlib.h>

#include "bri++.H"

//  This is the input file.  Files are tested in test-merstream-speed.C
//
const char *sequence =
">(one)\n"
"TTTTTTTTTT\n"
"AAAAAAAAAACGNTTTTTTTTTTNGGGGGGGGGGNAAAAAAAAANTTTTTTTTTT\n"
">(two)\n"
"TTTTTTTTAA\n"
">(thr)\n"
"TTTTTTTTAC\n"
">(fou)\n"
"T       T       T       T       T\n"
"T       T       T       A       G\n"
">(fiv)\n"
"T T T T T T T T A T C\n"
"\n"
">(six)\n"
"\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n";

//  This is the correct output
//
struct answers {
  u64bit  mer;
  u64bit  pos;
};

//  These numbers seem to be correct, but I haven't rigorously
//  verified each one.
//
answers correct[22] = {
  { 0x00000000000fffff, 0 },
  { 0x00000000000ffffc, 1 },
  { 0x00000000000ffff0, 2 },
  { 0x00000000000fffc0, 3 },
  { 0x00000000000fff00, 4 },
  { 0x00000000000ffc00, 5 },
  { 0x00000000000ff000, 6 },
  { 0x00000000000fc000, 7 },
  { 0x00000000000f0000, 8 },
  { 0x00000000000c0000, 9 },
  { 0x0000000000000000, 10 },
  { 0x0000000000000001, 11 },
  { 0x0000000000000006, 12 },
  { 0x00000000000fffff, 23 },
  { 0x00000000000aaaaa, 34 },
  { 0x00000000000fffff, 55 },
  { 0x00000000000ffff0, 0 },
  { 0x00000000000ffff1, 0 },
  { 0x00000000000ffff2, 0 },
  { 0x00000000000ffff3, 0 },
  { 0x00000000000fffcd, 1 },
  { 0x00000000000fffff, 0 }
};


int
main(int argc, char **argv) {
  merStream  *S = 0L;
  FILE       *F = 0L;
  u32bit      c = 0;
  u32bit      e = 0;

  S = new merStream(10, sequence, strlen(sequence));
  c = 0;
  while (S->nextMer()) {
    if (c > 21)
      fprintf(stderr, "merStream(char): Too many mers in stream.\n"), e++;
    else if ((S->thePosition() != correct[c].pos) || (S->theFMer() != correct[c].mer))
      fprintf(stderr, "merStream(char): "u64bitHEX"/"u64bitFMT" != correct: "u64bitHEX"/"u64bitFMT".\n",
              S->theFMer(), S->thePosition(),
              correct[c].mer, correct[c].pos), e++;
    c++;
  }
  delete S;


  F = fopen("test-merstream.junk", "w");
  fwrite(sequence, sizeof(char), strlen(sequence) + 1, F);
  fclose(F);

  S = new merStream(10, "test-merstream.junk");
  c = 0;
  while (S->nextMer()) {
    if (c > 21)
      fprintf(stderr, "merStream(char): Too many mers in stream.\n"), e++;
    else if ((S->thePosition() != correct[c].pos) || (S->theFMer() != correct[c].mer))
      fprintf(stderr, "merStream(file): "u64bitHEX"/"u64bitFMT" != correct: "u64bitHEX"/"u64bitFMT".\n",
              S->theFMer(), S->thePosition(),
              correct[c].mer, correct[c].pos), e++;
    c++;
  }
  delete S;


  unlink("test-merstream.junk");
}

