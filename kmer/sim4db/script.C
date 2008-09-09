#include "sim4db.H"


//  Parse the command line to create a sim4command object
//
//  [-f|-r] -e ESTid -D GENid GENlo GENhi
//
//    -f  Forward only
//    -r  Reverse only
//    -D  genSeqIID genLo genHi
//    -e  estSeqIID
//
//
char*
getNextScript(u32bit     &ESTiid,
              u32bit     &GENiid, u32bit &GENlo, u32bit &GENhi,
              bool       &doForward,
              bool       &doReverse,
              readBuffer *scriptFile) {

  char x = scriptFile->read();

  //  Skip any white space in the file
  //
  while ((scriptFile->eof() == false) && (whitespaceSymbol[x]))
    x = scriptFile->read();

  //  Exit if we're all done.
  //
  if (scriptFile->eof())
    return(0L);

  u32bit  linePos = 0;
  u32bit  lineMax = 128;
  char   *line    = new char [lineMax];

  //  Copy the line from the readBuffer into our storage
  //
  while ((scriptFile->eof() == false) && (x != '\n')) {
    line[linePos++] = x;
    x = scriptFile->read();
  }
  line[linePos] = 0;

  //  Decode the line
  //
  u32bit         argWords = 0;
  splitToWords   words(line);

  while (words.getWord(argWords)) {
    switch (words.getWord(argWords)[1]) {
      case 'f':
        doForward = true;
        doReverse = false;
        break;
      case 'r':
        doForward = false;
        doReverse = true;
        break;
      case 'D':
        GENiid = strtou32bit(words.getWord(++argWords), 0L);
        GENlo  = strtou32bit(words.getWord(++argWords), 0L);
        GENhi  = strtou32bit(words.getWord(++argWords), 0L);
        break;
      case 'e':
        ESTiid = strtou32bit(words.getWord(++argWords), 0L);
        break;
      default:
        //fprintf(stderr, "Unknown option '%s'\n", words.getWord(argWords));
        break;
    }

    argWords++;
  }

  return(line);
}
