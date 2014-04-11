#ifndef SIM4B1_S_H
#define SIM4B1_S_H

#define CONTINUOUS_SEED   10
#define SPACED_SEED       11

struct position_t {
  int    begin;
  int    end;
  int    width;
  int    result_shifts;
};

class mss_t {
public:
  mss_t() {};
  mss_t(char seed[32]);
  ~mss_t() {};

  int     type;
  uint64  mask;
  int     masknum;
  int     seedLength;
  int     matchedLength;
  uint64  masks[64];      //  Fails assert in sim4b1_s.C if exceeded
  int     shifts[64];

  uint64  mask_shift(uint64 ecode);
};

#endif

