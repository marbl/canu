#include "bri++.H"

#include <stdio.h>
#include <string.h>


char*
FastAWrapper::moveToNextName(char *ptr, u32bit iid) {

  //  XXX:
  //
  //  Sadly, if our index has saved only the id, we don't know the
  //  length of this field.  We must blindly search for the start of
  //  the next id.
  //
  //  If we have the whole defline, then we know the length.
  //
  if (_theGlobalDesc._indexType & FASTA_INDEX_MASK == FASTA_INDEX_PLUS_IDS) {
    while (*ptr)
      ptr++;
    ptr++;
  } else {
    ptr += _theSeqs[iid]._headerLen + 1;
  }

  return(ptr);
}



bool
FastAWrapper::find(char *id) {
  u32bit         iid = 0;
  char          *ptr = _theNames;
  char          *res = 0L;

  //
  //  need to trim off the whitespace from the front/end of the id
  //

  //  This is just too ugly.
  //
#if 0
  u32bit         len = 0;

  char *tmpid = new char [ strlen(id) + 1];
  while (*id && isspace(*id))
    id++;

  for (len=0; *id; len++, id++)
    tmpid[len] = *id;

  if (len == 0) {
    delete [] tmpid;
    return(false);
  }

  while (isspace(tmpid[len]))
    len--;

  len++;
  tmpid[len] = 0;

  id = tmpid;
#endif

  while (iid < _theGlobalDesc._numberOfSequences) {
    res = strstr(ptr, id);
    if (res) {
      _filebuffer->seek(_theSeqs[iid]._headerStart);
      _currentSequenceNumber = iid;
      return(true);
    }

    ptr = moveToNextName(ptr, iid);

    iid++;
  }
    
  return(false);
}



