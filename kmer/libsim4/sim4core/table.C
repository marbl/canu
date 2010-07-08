#include "sim4.H"

//  The position of a mer (word) is the position of the last base
//  (base-based).  Note that the sequence starts at position 1.
//
//           11111111112
//  12345678901234567890
//  acgggctactcgaggcta
//
//  First mer is at position 12.
//


void
Sim4::add_word(int ecode, int pos) {
  struct hash_node *h;
  int hval;

  hval = ecode & HASH_SIZE;

  //  Find the word in the hash table
  //
  for (h = hashtable->table[hval]; h; h = h->link)
    if (h->ecode == ecode)
      break;

  //  Didn't find the word?  Add a new one!
  //
  if (h == NULL) {
    h       = hashtable->nodes + hashtable->nodesused++;
    h->link = hashtable->table[hval];
    hashtable->table[hval] = h;

    h->ecode = ecode;
    h->pos   = -1;
  }

  //  Set the position -- this keeps a list of words from high
  //  position to low position.
  //
  hashtable->nextPos[pos] = h->pos;
  h->pos = pos;
}


void
Sim4::bld_table(char *s, int len, mss_t MSS, int type) {
  u64bit ecode;
  int i, j, masked_ecode;
  char *t;

  if (type == PERM) {
    mask = (1 << (2*MSS.seedLength-2)) - 1;   /* LLL 6/16/10 we are setting this for continuous seeds, where 2*seedLength=matchedLength; no effect if seed is spaced */
    hashtable = &phashtable;
    return; 
  }

  /* perform initializations */
  if (type == INIT) {
    mask = (1 << (2*MSS.seedLength-2)) - 1;  /* LLL 6/16/10 we are setting this for continuous seeds, where 2*seedLength=matchedLength; no effect if seed is spaced */

    hashtable = &phashtable;

    if (phashtable.nextPos) {
      delete [] phashtable.nextPos;
      delete [] phashtable.nodes;
    }

    phashtable.nextPos   = new int [len+1];
    phashtable.nodes     = new struct hash_node [len+1];
    phashtable.nodesused = 0;

    for (i=0; i<HASH_SIZE+1; ++i)
      phashtable.table[i] = u64bitZERO;
  } else if (type == TEMP) {
    mask = (1 << (2*MSS.seedLength-2)) - 1;  /* LLL 6/16/10 we are setting this for continuous seeds, where 2*seedLength=matchedLength; no effect if seed is spaced */

    hashtable = &thashtable;

    if (thashtable.nextPos) {
      delete [] thashtable.nextPos;
      delete [] thashtable.nodes;
    }

    thashtable.nextPos   = new int [len+1];
    thashtable.nodes     = new struct hash_node [len+1];
    thashtable.nodesused = 0;

    for (i=0; i<HASH_SIZE+1; ++i)
      thashtable.table[i] = 0L;
  } else {
    fprintf(stderr, "unknown type in bld_table: %d\n", type);
  }

  // skip any word containing an N/X

  int emer;

  //  This is because seq-1 is passed in
  //
  t = s+1;

  if (MSS.type == CONTINUOUS_SEED) {
    for (i=1; (i<=len) && *t; ) {
    restart_c:
      ecode = u64bitZERO;

      for (j=1; (j<MSS.seedLength) && (i<=len) && *t; ++j) {
        emer = encoding[(int)(*t++)];
        i++;

        if (emer < 0)
          goto restart_c;

          ecode <<= 2;
          ecode  |= emer;
        }
   
      for (; (i<=len) && *t; ) {
        emer = encoding[(int)(*t++)];
        i++;

        if (emer < 0)
          goto restart_c;

        ecode  &= mask;
        ecode <<= 2;
        ecode  |= emer;

        add_word(ecode, (int)(t-s-1));
      }
    }
  } else {
    /* SPACED_SEED */
    for (i=1; (i<=len) && *t; ) {
    restart_s: 
      ecode = u64bitZERO;

      for (j=1; (j<MSS.seedLength) && (i<=len) && *t; ++j) {
        emer = encoding[(int)(*t++)];
        i++;

        if (emer < 0)
          goto restart_s;

        ecode <<= 2;
        ecode  |= emer;
      }
    
      for (; (i<=len) && *t; ) {
        emer = encoding[(int)(*t++)];
        i++;

        if (emer < 0)
          goto restart_s;

        ecode  &= MSS.mask;
        ecode <<= 2;
        ecode  |= emer;

#if 1
        // much cheaper...

        for (j=masked_ecode=0; j<MSS.masknum; j++)
          masked_ecode += (ecode & MSS.masks[j]) >> MSS.shifts[j];

        add_word(masked_ecode, (int)(t-s-1));
#else
        add_word(masked_shift(ecode), (int)(t-s-1)); 
#endif
      }
    }
  }
}
