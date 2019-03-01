
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
 *  This file is derived from:
 *
 *    src/utility/decodeBooleanString.C
 *    src/utility/findKeyAndValue.H
 *    src/utility/splitToWords.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-AUG-13
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "strings.H"



uint64
scaledNumber(uint64 n, uint32 div) {

  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;

  return(n);
}


char
scaledUnit(uint64 n, uint32 div) {
  char u = ' ';

  if (n > 9999)  {  n /= div; u = 'k';  }
  if (n > 9999)  {  n /= div; u = 'M';  }
  if (n > 9999)  {  n /= div; u = 'G';  }
  if (n > 9999)  {  n /= div; u = 'T';  }
  if (n > 9999)  {  n /= div; u = 'P';  }
  if (n > 9999)  {  n /= div; u = 'E';  }
  if (n > 9999)  {  n /= div; u = 'Z';  }
  if (n > 9999)  {  n /= div; u = 'Y';  }

  return(u);
}



const char *
scaledName(uint64 n, uint32 div) {
  const char *s = "";

  if (n > 9999)  {  n /= div; s = " thousand";  }
  if (n > 9999)  {  n /= div; s = " million";  }
  if (n > 9999)  {  n /= div; s = " billion";  }
  if (n > 9999)  {  n /= div; s = " trillion";  }
  if (n > 9999)  {  n /= div; s = " quadrillion";  }
  if (n > 9999)  {  n /= div; s = " quintillion";  }
  if (n > 9999)  {  n /= div; s = " sextillion";  }
  if (n > 9999)  {  n /= div; s = " septillion";  }

  return(s);
}



bool
decodeBoolean(char *value) {
  bool ret = false;

  switch (value[0]) {
    case '0':
    case 'f':
    case 'F':
    case 'n':
    case 'N':
      ret = false;
      break;
    case '1':
    case 't':
    case 'T':
    case 'y':
    case 'Y':
      ret = true;
      break;
    default:
      fprintf(stderr, "decodeBoolean()-- unrecognized value '%s'\n", value);
      break;
  }

  return(ret);
}



//  Returns true if a key and value are found.  line is modified.
//  Returns true, with value == NULL, if no delimiter is found.
//  Returns false if the line is blank, or is a comment.
//
bool
KeyAndValue::find(char *line) {

  key_ = NULL;
  val_ = NULL;

  if (line == NULL)
    return(false);

  key_ = line;

  while  (isspace(*key_) == true)        //  Spaces before the key
    key_++;

  if ((iscomment(*key_) == true) ||      //  If we're at a comment right now, there is no key
      (*key_ == 0)) {                    //  and we return failure.
    key_ = NULL;
    val_ = NULL;
    return(false);
  }

  val_ = key_;                           //  We're at the key now

  while ((*val_ != 0) &&
         (isdelimiter(*val_) == false))  //  The key cannot contain a delimiter.
    val_++;

  if (*val_ == 0) {                      //   If at the end of the string, there isn't a
    val_ = NULL;                         //   value, but we'll return true and a key anyway.
    return(true);
  }

  *val_++ = 0;

  while (isdelimiter(*val_) == true) {   //  Spaces or delimiter after the key
    *val_ = 0;
    val_++;
  }

  if (*val_ == 0)                        //  And there is no value, must be a filename.
    return(true);

  char *eol = val_;                      //  We're at the value now

  //  If quoted, all we need to do is find the other quote and stop.
  if ((*val_ == '"') ||
      (*val_ == '\'')) {
    val_++;
    eol++;

    while (*eol != '"')                  //  The value itself.
      eol++;                             //  The value CAN contain delimiters and comment markers.

    *eol = 0;
  }

  //  Otherwise, not quoted.  Find the first comment marker (or eol) then backup to the first non-space.
  else {
    while (iscomment(*eol) == false)     //  The value MUST NOT contain delimiters or comment markers.
      eol++;                             //  But it can contains spaces and other nasties.

    eol--;                               //  Back up off the comment or eol.

    while (isspace(*eol) == true)        //  And keep backing up as long as we're a space.
      eol--;

    eol++;                               //  Move past the last non-space, non-comment

    *eol = 0;                            //  And terminate the value
  }

  return(true);
}



splitToWords::splitToWords(const char *string, splitType type) {
  _wordsLen  = 0;
  _wordsMax  = 0;
  _words     = NULL;

  _charsLen = 0;
  _charsMax = 0;
  _chars    = NULL;

  if (string)
    split(string, type);
}



splitToWords::~splitToWords() {
  delete [] _chars;
  delete [] _words;
}



void
splitToWords::split(const char *line, splitType type) {

  _wordsLen = 0;        //  Initialize to no words
  _charsLen = 0;        //  and no characters.

  if (line == NULL)     //  Bail if there isn't a line to process.
    return;

  //  Count the number of words and chars in the input line, then make
  //  sure there is space for us to store them.

  while (line[_charsLen] != 0)
    if (isSeparator(line[_charsLen++], type))
      _wordsLen++;

  resizeArray(_words, 0, _wordsMax, _wordsLen + 1, resizeArray_doNothing);
  resizeArray(_chars, 0, _charsMax, _charsLen + 1, resizeArray_doNothing);

  //  Clear all the words pointers, and copy the input line to our storage.
  //  This greatly simplifies the loop, as we don't need to worry about
  //  terminating the final word.

  memset(_words, 0,    sizeof(char *) * (_wordsLen + 1));
  memcpy(_chars, line, sizeof(char)   * (_charsLen + 1));

  //  Scan the line copy, converting word separators to NUL bytes.
  //  counting and saving the start of each word in _words.

  _wordsLen = 0;

  for (uint32 st=1, ii=0; ii < _charsLen; ii++) {
    if (isSeparator(line[ii], type)) {      //  If the character is a word
      _chars[ii] = 0;                       //  separator, convert to NUL,
      st         = true;                    //  and flag the next character
    }                                       //  as the start of a new word.

    else if (st) {                          //  Otherwise, if this is the
      _words[_wordsLen++] = _chars + ii;    //  start of a word, make
      st                  = false;          //  a new word.
    }
  }
}



