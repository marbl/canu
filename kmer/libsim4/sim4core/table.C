#include "sim4.H"

//
// add_word - add a word to the table of critical words
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

  //  Set the position
  //
  hashtable->nextPos[pos] = h->pos;
  h->pos = pos;
}


void
Sim4::bld_table(char *s, int len, int in_W, int type) {
  int ecode;
  int i, j;
  char *t;

  //fprintf(stdout, "Building table: len=%d type=%d\n", len, type);

  if (type == PERM) {
    mask = (1 << (in_W+in_W-2)) - 1;
    hashtable = &phashtable;
    return; 
  }

  /* perform initializations */
  if (type == INIT) {
    mask = (1 << (in_W+in_W-2)) - 1;

    hashtable = &phashtable;

    if (phashtable.nextPos) {
      delete [] phashtable.nextPos;
      delete [] phashtable.nodes;
    }

    phashtable.nextPos   = new int [len+1];               //(int *)ckalloc((len+1)*sizeof(int));
    phashtable.nodes     = new struct hash_node [len+1];  //(struct hash_node *)ckalloc((len+1) * sizeof(struct hash_node));
    phashtable.nodesused = 0;

    for (i=0; i<HASH_SIZE+1; ++i)
      phashtable.table[i] = 0L;
  } else if (type == TEMP) {
    mask = (1 << (in_W+in_W-2)) - 1;

    hashtable = &thashtable;

    if (thashtable.nextPos) {
      delete [] thashtable.nextPos;
      delete [] thashtable.nodes;
    }

    thashtable.nextPos   = new int [len+1];               //(int *)ckalloc((len+1)*sizeof(int));
    thashtable.nodes     = new struct hash_node [len+1];  //(struct hash_node *)ckalloc((len+1) * sizeof(struct hash_node));
    thashtable.nodesused = 0;

    for (i=0; i<HASH_SIZE+1; ++i)
      thashtable.table[i] = 0L;
  } else {
    fprintf(stderr, "unknown type in bld_table: %d\n", type);
  }

  /* skip any word containing an N/X  */

  int emer;

  t = s+1;

  for (i=1; (i<=len) && *t; ) {
  restart: 
    ecode = 0L;

    for (j=1; (j<in_W) && (i<=len) && *t; ++j) {
      emer = encoding[(int)(*t++)];
      i++;

      if (emer < 0)
        goto restart;

      ecode <<= 2;
      ecode  |= emer;
    }
    
    for (; (i<=len) && *t; ) {
      emer = encoding[(int)(*t++)];
      i++;

      if (emer < 0)
        goto restart;

      ecode  &= mask;
      ecode <<= 2;
      ecode  |= emer;

      add_word(ecode, (int)(t-s-1)); 
    }
  }
}
