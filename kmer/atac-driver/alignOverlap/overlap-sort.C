

int
sortMatches1(const void *a, const void *b) {
  const match_t *A = *((const match_t * const *)a);
  const match_t *B = *((const match_t * const *)b);

  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);
  if (A->ori1 > B->ori1)  return(-1);
  if (A->ori1 < B->ori1)  return(1);
  return(0);
}

int
sortMatches2(const void *a, const void *b) {
  const match_t *A = *((const match_t * const *)a);
  const match_t *B = *((const match_t * const *)b);

  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->pos2 < B->pos2)  return(-1);
  if (A->pos2 > B->pos2)  return(1);
  if (A->len2 > B->len2)  return(-1);
  if (A->len2 < B->len2)  return(1);
  if (A->ori2 > B->ori2)  return(-1);
  if (A->ori2 < B->ori2)  return(1);
  return(0);
}



int
spanCompare(const void *a, const void *b) {
  const span_t *A = *((const span_t * const *)a);
  const span_t *B = *((const span_t * const *)b);

  if (A->_iid < B->_iid)  return(-1);
  if (A->_iid > B->_iid)  return(1);
  if (A->_beg < B->_beg)  return(-1);
  if (A->_beg > B->_beg)  return(1);
  return(0);
}
