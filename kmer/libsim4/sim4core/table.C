#include "sim4.H"

//
// add_word - add a word to the table of critical words
//
void
Sim4::add_word(int ecode, int pos) {
  struct hash_node *h;
  int hval;

  hval = ecode & HASH_SIZE;
  for (h = hashtable->table[hval]; h; h = h->link)
    if (h->ecode == ecode)
      break;
  if (h == NULL) {
    h = hashtable->nodes + hashtable->nodesused++;
    //(struct hash_node *) ckalloc (sizeof(struct hash_node));
    h->link = hashtable->table[hval];
    hashtable->table[hval] = h;
    h->ecode = ecode;
    h->pos   = -1;
  }
  hashtable->nextPos[pos] = h->pos;
  h->pos = pos;
}


void
Sim4::bld_table(char *s, int len, int in_W, int type)
{
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
  t = s+1;
  for (i=1; (i<=len) && *t; ) {
  restart: 
    ecode = 0L;
    for (j=1; (j<in_W) && (i<=len) && *t; ++j) {
      int tmp = encoding[*t++];
      ++i;
      if (tmp<0)
        goto restart;
      ecode = (ecode << 2) + tmp;
    }
    
    for (; (i<=len) && *t; ) {
      int tmp = encoding[*t++]; i++;
      if (tmp<0) goto restart;   
      ecode = ((ecode & mask) << 2) + tmp;
      add_word(ecode, (int)(t-s-1)); 
    }
  }
}
