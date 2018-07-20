
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-21 to 2015-JAN-30
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-FEB-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

class KeyAndValue {
public:

  KeyAndValue(char *line = NULL) {
    find(line);
  };

  ~KeyAndValue() {
  }

  bool     find(char *line);

  char    *key(void)           { return(key_); };
  char    *value(void)         { return(val_); };

  bool     value_bool(void)    { return((val_[0] == 't') || (val_[0] == 'T') || (val_[0] == '1')); };

  int32    value_int32(void)   { return(strtol (val_, NULL, 10)); };
  int64    value_int64(void)   { return(strtoll(val_, NULL, 10)); };

  uint32   value_uint32(void)  { return(strtoul (val_, NULL, 10)); };
  uint64   value_uint64(void)  { return(strtoull(val_, NULL, 10)); };

  float    value_float(void)   { return(strtof(val_, NULL)); };
  double   value_double(void)  { return(strtod(val_, NULL)); };

public:
  bool    iscomment(char c)    { return((c == '!') || (c == '#') || (c == 0)); };
  bool    isdelimiter(char c)  { return((c == ':') || (c == '=') || isspace(c)); };

  char   *key_;
  char   *val_;
};



//  Returns true if a key and value are found.  line is modified.
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
