#ifndef OUTPUTMER_H
#define OUTPUTMER_H

#include "libbri.H"
#include "mcDescription.H"

void
outputMer(bitPackedFileWriter *DAT,
          mcDescription       *mcd,
          u64bit               prefix,
          u64bit               check,
          u32bit               count);


#endif  //  OUTPUTMER_H
