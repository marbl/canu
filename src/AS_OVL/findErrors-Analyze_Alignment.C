

#include "findErrors.H"






//  Add vote val to G.reads[sub] at sequence position  p
static
void
Cast_Vote(feParameters *G,
          Vote_Value_t val,
          int32        p,
          int32        sub) {

  switch (val) {
    case DELETE:    if (G->reads[sub].vote[p].deletes  < MAX_VOTE)  G->reads[sub].vote[p].deletes++;   break;
    case A_SUBST:   if (G->reads[sub].vote[p].a_subst  < MAX_VOTE)  G->reads[sub].vote[p].a_subst++;   break;
    case C_SUBST:   if (G->reads[sub].vote[p].c_subst  < MAX_VOTE)  G->reads[sub].vote[p].c_subst++;   break;
    case G_SUBST:   if (G->reads[sub].vote[p].g_subst  < MAX_VOTE)  G->reads[sub].vote[p].g_subst++;   break;
    case T_SUBST:   if (G->reads[sub].vote[p].t_subst  < MAX_VOTE)  G->reads[sub].vote[p].t_subst++;   break;
    case A_INSERT:  if (G->reads[sub].vote[p].a_insert < MAX_VOTE)  G->reads[sub].vote[p].a_insert++;  break;
    case C_INSERT:  if (G->reads[sub].vote[p].c_insert < MAX_VOTE)  G->reads[sub].vote[p].c_insert++;  break;
    case G_INSERT:  if (G->reads[sub].vote[p].g_insert < MAX_VOTE)  G->reads[sub].vote[p].g_insert++;  break;
    case T_INSERT:  if (G->reads[sub].vote[p].t_insert < MAX_VOTE)  G->reads[sub].vote[p].t_insert++;  break;
    case NO_VOTE:
      break;
    default :
      fprintf(stderr, "ERROR:  Illegal vote type\n");
      break;
  }
}



//  Return the substitution vote corresponding to  Ch .
static
Vote_Value_t
Matching_Vote(char ch) {

  switch  (ch) {
    case 'a':  return(A_SUBST);  break;
    case 'c':  return(C_SUBST);  break;
    case 'g':  return(G_SUBST);  break;
    case 't':  return(T_SUBST);  break;
  }

  fprintf(stderr, "Matching_Vote()-- invalid letter '%c'\n", ch);

  return(NO_VOTE);
}




//  Analyze the delta-encoded alignment in  delta[0 .. (deltaLen - 1)]
//  between  a_part  and  b_part  and store the resulting votes
//  about the a sequence in  G->reads[sub]. The alignment starts
//   a_offset  bytes in from the start of the a sequence in  G->reads[sub] .
//   a_len  and  b_len  are the lengths of the prefixes of  a_part  and
//   b_part , resp., that align.


void
Analyze_Alignment(Thread_Work_Area_t *wa,
                  char   *a_part, int32 a_len, int32 a_offset,
                  char   *b_part, int32 b_len,
                  int32   sub) {

  assert(a_len >= 0);
  assert(b_len >= 0);

  int32  ct = 0;

  wa->vote[ct].frag_sub   = -1;
  wa->vote[ct].align_sub  = -1;
  wa->vote[ct].vote_val   = A_SUBST;   // Dummy value
  ct++;

  int32  i = 0;
  int32  j = 0;
  int32  p = 0;

  for (int32 k=0;  k<wa->ped.deltaLen;  k++) {
    for (int32 m=1;  m<abs (wa->ped.delta[k]);  m++) {
      if (a_part[i] != b_part[j]) {
        wa->vote[ct].frag_sub = i;
        wa->vote[ct].align_sub = p;

        switch (b_part[j]) {
          case 'a':  wa->vote[ct].vote_val = A_SUBST;  break;
          case 'c':  wa->vote[ct].vote_val = C_SUBST;  break;
          case 'g':  wa->vote[ct].vote_val = G_SUBST;  break;
          case 't':  wa->vote[ct].vote_val = T_SUBST;  break;
          default :
            fprintf(stderr, "ERROR:[1] Bad sequence '%c' 0x%02x)\n", b_part[j], b_part[j]);
            assert(0);
        }

        ct++;
      }

      i++;
      j++;
      p++;
    }

    if (wa->ped.delta[k] < 0) {
      wa->vote[ct].frag_sub = i - 1;
      wa->vote[ct].align_sub = p;

      switch (b_part[j]) {
        case 'a':  wa->vote[ct].vote_val = A_INSERT;  break;
        case 'c':  wa->vote[ct].vote_val = C_INSERT;  break;
        case 'g':  wa->vote[ct].vote_val = G_INSERT;  break;
        case 't':  wa->vote[ct].vote_val = T_INSERT;  break;
        default :
          fprintf(stderr, "ERROR:[2] Bad sequence '%c' 0x%02x)\n", b_part[j], b_part[j]);
          assert(0);
      }

      ct++;

      j++;
      p++;

    } else {
      wa->vote[ct].frag_sub = i;
      wa->vote[ct].align_sub = p;
      wa->vote[ct].vote_val = DELETE;

      ct++;

      i++;
      p++;
    }
  }

  while (i < a_len) {
    if (a_part[i] != b_part[j]) {
      wa->vote[ct].frag_sub = i;
      wa->vote[ct].align_sub = p;

      switch (b_part[j]) {
        case 'a':  wa->vote[ct].vote_val = A_SUBST;  break;
        case 'c':  wa->vote[ct].vote_val = C_SUBST;  break;
        case 'g':  wa->vote[ct].vote_val = G_SUBST;  break;
        case 't':  wa->vote[ct].vote_val = T_SUBST;  break;
        default :
          fprintf(stderr, "ERROR:[3] Bad sequence '%c' 0x%02x)\n", b_part[j], b_part[j]);
          assert(0);
      }

      ct++;
    }

    i++;
    j++;
    p++;
  }

  wa->vote[ct].frag_sub  = i;
  wa->vote[ct].align_sub = p;

  for (int32 i=1; i<=ct; i++) {
    int32  prev_match = wa->vote[i].align_sub - wa->vote[i - 1].align_sub - 1;
    int32  p_lo = (i == 1 ? 0 : wa->G->End_Exclude_Len);
    int32  p_hi = (i == ct ? prev_match : prev_match - wa->G->End_Exclude_Len);

    if (prev_match >= wa->G->Kmer_Len) {
      for (int32 p=0;  p<p_lo;  p++)
        Cast_Vote(wa->G,
                  Matching_Vote(a_part[wa->vote[i-1].frag_sub + p + 1]),
                  a_offset + wa->vote[i-1].frag_sub + p + 1,
                  sub);

      for (int32 p=p_lo;  p<p_hi;  p++) {
        int32 k = a_offset + wa->vote[i-1].frag_sub + p + 1;

        if (wa->G->reads[sub].vote[k].confirmed < MAX_VOTE)
          wa->G->reads[sub].vote[k].confirmed++;
        if ((p < p_hi - 1)  && (wa->G->reads[sub].vote[k].no_insert < MAX_VOTE))
          wa->G->reads[sub].vote[k].no_insert++;
      }

      for (int32 p=p_hi;  p<prev_match;  p++)
        Cast_Vote(wa->G,
                  Matching_Vote(a_part[wa->vote[i - 1].frag_sub + p + 1]),
                  a_offset + wa->vote[i - 1].frag_sub + p + 1,
                  sub);
    }

    // Don't allow consecutive inserts

    if ((i < ct) &&
        ((prev_match > 0) || (wa->vote[i - 1].vote_val <= T_SUBST) || (wa->vote[i].vote_val <= T_SUBST))) {
      int32 next_match = wa->vote[i + 1].align_sub - wa->vote[i].align_sub - 1;

      if (prev_match + next_match >= wa->G->Vote_Qualify_Len)
        Cast_Vote(wa->G,
                  wa->vote[i].vote_val,
                  a_offset + wa->vote[i].frag_sub,
                  sub);
    }
  }

#if 0
  printf(">a_part\n");
  for (j = 0;  a_part[j] != '\0';  j++) {
    putchar (wa->G->reads[sub].wa->vote[a_offset + j].confirmed ? '*' : ' ');
  }
  putchar ('\n');
#endif
}



