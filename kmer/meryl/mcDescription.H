#ifndef MCDESCRIPTION_H
#define MCDESCRIPTION_H

#include "libbri.H"

class mcDescription {
public:
  mcDescription() {};
  mcDescription(mcDescription const &m) {
    _merSizeInBases     = m._merSizeInBases;
    _merSizeInBits      = m._merSizeInBits;
    _tableSizeInBits    = m._tableSizeInBits;
    _tableSizeInEntries = m._tableSizeInEntries;
    _chckBits           = m._chckBits;
    _chckMask           = m._chckMask;
    _hashWidth          = m._hashWidth;
    _hashMask           = m._hashMask;
    _actualNumberOfMers = m._actualNumberOfMers;
  };
  ~mcDescription() {};

  u32bit      _merSizeInBases;
  u32bit      _merSizeInBits;

  u32bit      _tableSizeInBits;
  u64bit      _tableSizeInEntries;

  u32bit      _chckBits;
  u64bit      _chckMask;

  u32bit      _hashWidth;
  u64bit      _hashMask;

  u64bit      _actualNumberOfMers;

  u64bit HASH(u64bit a) {
    return((a >> _chckBits) & _hashMask);
  }

  void         write(bitPackedFileWriter *DAT) {
    DAT->putBits(_merSizeInBases,        32);
    DAT->putBits(_tableSizeInBits,       32);
    DAT->putBits(_chckBits,              32);
    DAT->putBits(_hashWidth,             32);
    DAT->putBits(_actualNumberOfMers,    64);
  };

  void         read(bitPackedFileReader *DAT) {
    _merSizeInBases        = (u32bit)DAT->getBits(32);
    _merSizeInBits         = 2 * _merSizeInBases;
    _tableSizeInBits       = (u32bit)DAT->getBits(32);
    _tableSizeInEntries    = u64bitONE << _tableSizeInBits;
    _chckBits              = (u32bit)DAT->getBits(32);
    _chckMask              = u64bitMASK(_chckBits);
    _hashWidth             = (u32bit)DAT->getBits(32);
    _hashMask              = u64bitMASK(_tableSizeInBits);
    _actualNumberOfMers    = DAT->getBits(64);
  };

  void         print(FILE *F) {
    fprintf(F, "merSizeInBases     = "u32bitFMT"\n", _merSizeInBases);
    fprintf(F, "merSizeInBits      = "u32bitFMT"\n", _merSizeInBits);
    fprintf(F, "tableSizeInBits    = "u32bitFMT"\n", _tableSizeInBits);
    fprintf(F, "tableSizeInEntries = "u64bitFMT"\n", _tableSizeInEntries);
    fprintf(F, "chckBits           = "u32bitFMT"\n", _chckBits);
    fprintf(F, "chckMask           = "u64bitHEX"\n", _chckMask);
    fprintf(F, "hashWidth          = "u32bitFMT"\n", _hashWidth);
    fprintf(F, "hashMask           = "u64bitHEX"\n", _hashMask);
    fprintf(F, "actualNumberOfMers = "u64bitFMT"\n", _actualNumberOfMers);
  }
};

#endif  //  MCDESCRIPTION_H
