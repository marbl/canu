#include <unistd.h>
#include <time.h>
#include <math.h>

#include "util++.H"

struct header_s {
  u64bit  t1;
  char    s1[570];
  u64bit  t2;
};

struct record_s {
  u64bit   t1;
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

  fprintf(stderr, "header t1 "u64bitHEX" '%s' t2 "u64bitHEX"\n", hh->t1, hh->s1, hh->t2);
  RF->getRecord(&r);  fprintf(stderr, "record "u64bitFMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "u64bitFMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "u64bitFMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "u64bitFMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "u64bitFMT" '%s'\n", r.t1, r.s1);
  RF->getRecord(&r);  fprintf(stderr, "record "u64bitFMT" '%s'\n", r.t1, r.s1);

  delete RF;

  return(0);
}
