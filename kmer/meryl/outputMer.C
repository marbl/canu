#include "outputMer.H"

//
//  XXX:  This has a known problem with large counts -- if the highest bit is
//  set, it will write an infinite count.
//
//  XXX:  Why isn't 'prefix' being used?
//

void
outputMer(bitPackedFileWriter *DAT,
          mcDescription&       mcd,
          u64bit               prefix,
          u64bit               check,
          u32bit               count) {
  u64bit v = check & mcd._chckMask;
  u64bit q = count;

  //fprintf(stderr, "0x%016lx 0x%016lx -- %u\n", prefix, check, count);

  if (count == 0) {
    fprintf(stderr, "WARNING: outputMer() told to write a zero count mer!\n");
    fprintf(stderr, "         prefix = 0x%016lx\n", prefix);
    fprintf(stderr, "         check  = 0x%016lx\n", check);
  }

  if (count == 1) {
    //  Count is exactly one.  Special case -- set a flag in the check
    //  bits to denote a count of one.
    //
    v |= u64bitONE << mcd._chckBits;
    DAT->putBits(v, mcd._chckBits + 1);
  } else {
    if (count == 0) {
      //  Zero count; write a zero block, with the last block bit set.
      //
      DAT->putBits(v, mcd._chckBits + 1);
      DAT->putBits(u64bitONE << mcd._chckBits, mcd._chckBits + 1);
    } else {
      //  Count is > 1.
      //
      DAT->putBits(v, mcd._chckBits + 1);

      while (q != u64bitZERO) {

        //  Get the bits to write.  The bits are written in blocks
        //  of size mcd._chckBits, starting with the least significant
        //  block first.
        //
        v = q & mcd._chckMask;

        //  Strip off the bits we will soon write.
        //
        q >>= mcd._chckBits;

        //  If the count is now zero, this will be the last block.
        //  Mark it.
        //
        if (q == u32bitZERO)
          v |= u64bitONE << mcd._chckBits;

        //  Write the bits
        //
        DAT->putBits(v, mcd._chckBits + 1);
      }
    }
  }
}
