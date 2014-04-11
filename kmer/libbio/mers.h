#ifndef BIO_MERS_H
#define BIO_MERS_H



inline
uint64
reverseComplementMer(uint32 ms, uint64 fmer) {

  //  The interested reader shall consult bri-bits.h

  //  Reverse the mer
  //
  uint64 rmer = fmer;
  rmer = ((rmer >>  2) & 0x3333333333333333llu) | ((rmer <<  2) & 0xccccccccccccccccllu);
  rmer = ((rmer >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((rmer <<  4) & 0xf0f0f0f0f0f0f0f0llu);
  rmer = ((rmer >>  8) & 0x00ff00ff00ff00ffllu) | ((rmer <<  8) & 0xff00ff00ff00ff00llu);
  rmer = ((rmer >> 16) & 0x0000ffff0000ffffllu) | ((rmer << 16) & 0xffff0000ffff0000llu);
  rmer = ((rmer >> 32) & 0x00000000ffffffffllu) | ((rmer << 32) & 0xffffffff00000000llu);

  //  Complement the bases
  //
  rmer ^= 0xffffffffffffffffllu;

  //  Shift and mask out the bases not in the mer
  //
  rmer >>= 64 - ms * 2;
  rmer  &= uint64MASK(ms * 2);
  return(rmer);
}


//  Used for in seagen/encodedQuery.C (diagnostics) and
//  libbio/kmerhuge.H (in its merToString method).
inline
char *
uint64ToMerString(uint32 ms, uint64 mer, char *str) {
  for (uint32 i=0; i<ms; i++)
    str[ms-i-1] = bitsToLetter[(mer >> (2*i)) & 0x03];
  str[ms] = 0;
  return(str);
}


#if 0
#error this is not used anywhere
inline
uint64
stringToMer(uint32 ms, char *str) {
  uint64  mer = 0L;

  for (uint32 i=0; i<ms; i++) {
    mer <<= 2;
    mer  |= compressSymbol[str[i]];
  }

  return(mer);
}
#endif



#endif  //  BIO_MERS_H
