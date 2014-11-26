


static
bool
decodeBoolean(char *feature, char *value) {
  bool ret = false;

  //  Decodes a string with 0/1, false/true, no/yes into an integer flag.

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
      fprintf(stderr, "decodeBoolean()-- feature '%s' has unknown boolean value '%s'\n",
              feature, value);
      break;
  }

  return(ret);
}
