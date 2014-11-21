


class KeyAndValue {
public:

  KeyAndValue(char *line = NULL) {
    //initialize();
    if (line)
      find(line);
  };

  ~KeyAndValue() {
  }

  bool     find(char *line);

  char    *key(void)    { return(key_); };
  char    *value(void)  { return(val_); };

public:
#if 0
  void    initialize(void) {
    for (int ii=0; ii<256; ii++) {
      iscomment[ii]   = false;
      isdelimiter[ii] = false;
    }

    iscomment['!'] = true;
    iscomment['#'] = true;

    isdelimiter[':'] = true;
    isdelimiter['='] = true;
  };
#endif

  //  Code would be simpler if isdelimiter() included spaces.

  bool    iscomment(char c)    { return((c == '!') || (c == '#')); };
  bool    isdelimiter(char c)  { return((c == ':') || (c == '=') || isspace(c)); };
  
#if 0
  char    iscomment[256];
  char    isdelimiter[256];
#endif

  char   *key_;
  char   *val_;
};



//  Returns true if a key and value are found.
//
bool
KeyAndValue::find(char *line) {

  key_ = NULL;
  val_ = NULL;

  if (line == NULL)
    return(false);

  key_ = line;

  while  (isspace(*key_) == true)    //  Spaces before the key
    key_++;

  if ((iscomment(*key_) == true) ||
      (*key_ == 0)) {
    key_ = NULL;
    val_ = NULL;
    return(false);
  }

  val_ = key_;                       //  We're at the key now

  while (isdelimiter(*val_) == false) //  The key cannot contain a delimiter.
    val_++;

  *val_++ = 0;

  while (isdelimiter(*val_) == true) {  //  Spaces or delimiter after the key
    //*val_ = 0;
    val_++;
  }

  char *eol = val_;                       //  We're at the value now

  bool isQuoted = (*val_ == '"');

  if (isQuoted) {
    //  Quoted, so all we need to do is find the other quote and stop.

    val_++;
    eol++;

    while (*eol != '"')   //  The value itself.
      eol++;              //  The value CAN contain delimiters and comment markers.

    *eol = 0;

  } else {
    //  Not quoted.  Find the first space or comment.
      
    while ((isspace(*eol) == false) &&   //  The value itself.
           (iscomment(*eol) == false))   //  The value MUST NOT contain delimiters or comment markers.
      eol++;                         

    *eol = 0;
  }

  return(true);
}
