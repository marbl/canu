#ifndef MARK_FILE_H
#define MARK_FILE_H

#include "util++.H"


class markFile {
public:
  markFile(char *filename);
  ~markFile();

  void    checkTag(void);
  void    writeTag(bool bogus=false);

  u64bit  get(u64bit merNumber) {
    return(_bpf->getBits(merNumber + 256, 1));
  };

  void    set(u64bit merNumber) {
    if (_bpf->getBits(merNumber + 256, 1) == 0) {
      _numberMarked++;
      _haveWritten = true;
      _bpf->putBits(merNumber + 256, 1, 1);

      if (DEBUGlargest < merNumber)
        DEBUGlargest = merNumber;
    }
  };

  void    rewind(void) {
    _bpf->seek(32);
  };

  void    setLastMer(u64bit mernum) {
    _numberOfMers = mernum;
  };
  u64bit  numberOfMers(void) {
    return(_numberOfMers);
  };
  u64bit  numberMarked(void) {
    return(_numberMarked);
  };

private:
  bitPackedFile   *_bpf;

  u64bit DEBUGlargest;

  char             _filename[1024];
  char             _correctMagic[16];

  u64bit           _numberOfMers;
  u64bit           _numberMarked;

  bool             _haveWritten;

  friend class markFileMerger;
};




//  This is a class just so it can be a friend of markFile.
//
//  It needs lowlevel access to the bitPackedFile for efficiency.
//
class markFileMerger {
public:
  markFileMerger() {
    inNameLen     = 0;
    inNameMax     = 10 * 1024;
    inName        = new char     * [inNameMax];
    inFile        = new markFile * [inNameMax];
  };
  ~markFileMerger() {
    for (u32bit f=0; f<inNameLen; f++)
      delete inFile[f];
  };

  void   addFile(char *filetomerge) {
    if (inNameLen >= inNameMax) {
      //  Lazy programmer should figure out what the maximum number of
      //  open files is, and use that instead of 10 * 1024.
      //
      fprintf(stderr, "TOO MANY FILES!  Lazy programmer didn't implement resizing!\n");
      exit(1);
    }

    inName[inNameLen] = filetomerge;
    inFile[inNameLen] = new markFile(inName[inNameLen]);
    inNameLen++;
  };

  void   doMerge(char *outname) {
    markFile *outFile = new markFile(outname);

    //  All files should be the same length, but we should probably be
    //  storing some meta information somewhere.  Well, "probably" is a
    //  bit weak - we have no way of knowing when the stream of bits
    //  ends, so we definitely need to know how long (how many mers)
    //  each is.

    for (u32bit f=0; f<inNameLen; f++)
      fprintf(stderr, "Merging ["u32bitFMTW(2)"] with "u64bitFMT" mers "u64bitFMT" marked.\n",
              f,
              inFile[f]->numberOfMers(),
              inFile[f]->numberMarked());

    u64bit  numMers = inFile[0]->numberOfMers() / 64 + 1;

    //  want to read a word at a time, but bitpackedfile doesn't
    //  really support that nicely (it still does it's bit stuff) and
    //  markFile certainly doesn't support that.

    fprintf(stderr, u64bitFMT" words to merge.\n", numMers);
    speedCounter *SC = new speedCounter(" %8f Mwords (%8.5f Mwords/sec)\n", 1000000, 0xff, true);

    for (u64bit position = 0; position < numMers; position++) {
      u64bit   mergedWord = u64bitZERO;

      for (u32bit f=0; f<inNameLen; f++)
        mergedWord |= inFile[f]->_bpf->getBits(64);

      outFile->_bpf->putBits(mergedWord, 64);

      SC->tick();
    }

    delete SC;

    outFile->writeTag();

    delete outFile;
  };

private:
  u32bit          inNameLen;
  u32bit          inNameMax;
  char          **inName;
  markFile      **inFile;
};



#endif  //  MARK_FILE_H
