


#include "findErrors.H"



#if 0
void
Output_Details(feParameters &G) {

  for  (uint32 i = 0;  i < Num_Frags;  i ++) {
    printf (">%d\n", Lo_Frag_IID + i);

    for  (uint32 j = 0;  Frag[i].sequence[j] != '\0';  j ++)
      printf ("%3d: %c  %3d  %3d | %3d %3d %3d %3d | %3d %3d %3d %3d %3d\n",
              j,
              j >= Frag[i].clear_len ? toupper (Frag[i].sequence[j]) : Frag[i].sequence[j],
              Frag[i].vote[j].confirmed,
              Frag[i].vote[j].deletes,
              Frag[i].vote[j].a_subst,
              Frag[i].vote[j].c_subst,
              Frag[i].vote[j].g_subst,
              Frag[i].vote[j].t_subst,
              Frag[i].vote[j].no_insert,
              Frag[i].vote[j].a_insert,
              Frag[i].vote[j].c_insert,
              Frag[i].vote[j].g_insert,
              Frag[i].vote[j].t_insert);
    }
  }
}
#endif


void 
Output_Corrections(feParameters &G) {
  Correction_Output_t  out;

  errno = 0;
  FILE *fp = fopen(G.outputFileName, "wb");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", G.outputFileName, strerror(errno)), exit(1);


  for (uint32 i=0; i<G.readsLen; i++) {

    out.keep_left   = (G.reads[i].left_degree  < G.Degree_Threshold);
    out.keep_right  = (G.reads[i].right_degree < G.Degree_Threshold);
    out.type        = IDENT;
    out.pos         = 0;
    out.readID      = G.bgnID + i;

    AS_UTL_safeWrite(fp, &out, "correction1", sizeof(Correction_Output_t), 1);

    if (G.reads[i].sequence == NULL)
      // Deleted fragment
      continue;

    for (uint32 j=0; j<G.reads[i].clear_len; j++) {
      if  (G.reads[i].vote[j].confirmed < 2) {
        Vote_Value_t  vote      = DELETE;
        int32         max       = G.reads[i].vote[j].deletes;
        bool          is_change = true;

        if  (G.reads[i].vote[j].a_subst > max) {
          vote      = A_SUBST;
          max       = G.reads[i].vote[j].a_subst;
          is_change = (G.reads[i].sequence[j] != 'a');
        }

        if  (G.reads[i].vote[j].c_subst > max) {
          vote      = C_SUBST;
          max       = G.reads[i].vote[j].c_subst;
          is_change = (G.reads[i].sequence[j] != 'c');
        }

        if  (G.reads[i].vote[j].g_subst > max) {
          vote      = G_SUBST;
          max       = G.reads[i].vote[j].g_subst;
          is_change = (G.reads[i].sequence[j] != 'g');
        }

        if  (G.reads[i].vote[j].t_subst > max) {
          vote      = T_SUBST;
          max       = G.reads[i].vote[j].t_subst;
          is_change = (G.reads[i].sequence[j] != 't');
        }

        int32 haplo_ct  =  ((G.reads[i].vote[j].deletes >= MIN_HAPLO_OCCURS) +
                            (G.reads[i].vote[j].a_subst >= MIN_HAPLO_OCCURS) +
                            (G.reads[i].vote[j].c_subst >= MIN_HAPLO_OCCURS) +
                            (G.reads[i].vote[j].g_subst >= MIN_HAPLO_OCCURS) +
                            (G.reads[i].vote[j].t_subst >= MIN_HAPLO_OCCURS));

        int32 total  = (G.reads[i].vote[j].deletes +
                        G.reads[i].vote[j].a_subst +
                        G.reads[i].vote[j].c_subst +
                        G.reads[i].vote[j].g_subst +
                        G.reads[i].vote[j].t_subst);

        //  The original had a gargantuajn if test (five clauses, all had to be true) to decide if a record should be output.
        //  It was negated into many small tests if we should skip the output.
        //  A side effect is that we can abort a little earlier in two cases (and we don't even bother).

        //  (2 * max > total)
        if (2 * max <= total)
          continue;

        //  (total > 1)
        if (total <= 1)
          continue;

        //  (is_change == true)
        if (is_change == false)
          continue;

        //  ((haplo_ct < 2) || (G.Use_Haplo_Ct == false))
        if ((haplo_ct >= 2) && (G.Use_Haplo_Ct == true))
          continue;

        //  ((G.reads[i].vote[j].confirmed == 0) ||
        //   ((G.reads[i].vote[j].confirmed == 1) && (max > 6)))
        if ((G.reads[i].vote[j].confirmed > 0) &&
            ((G.reads[i].vote[j].confirmed != 1) || (max <= 6)))
          continue;

        //  Otherwise, output.

        out.type       = vote;
        out.pos        = j;

        AS_UTL_safeWrite(fp, &out, "correction2", sizeof(Correction_Output_t), 1);
      }


      if  (G.reads[i].vote[j].no_insert < 2) {
        Vote_Value_t  ins_vote = A_INSERT;
        int32         ins_max  = G.reads[i].vote[j].a_insert;

        if  (ins_max < G.reads[i].vote[j].c_insert) {
          ins_vote = C_INSERT;
          ins_max  = G.reads[i].vote[j].c_insert;
        }

        if  (ins_max < G.reads[i].vote[j].g_insert) {
          ins_vote = G_INSERT;
          ins_max  = G.reads[i].vote[j].g_insert;
        }

        if  (ins_max < G.reads[i].vote[j].t_insert) {
          ins_vote = T_INSERT;
          ins_max  = G.reads[i].vote[j].t_insert;
        }

        int32 ins_haplo_ct = ((G.reads[i].vote[j].a_insert >= MIN_HAPLO_OCCURS) +
                              (G.reads[i].vote[j].c_insert >= MIN_HAPLO_OCCURS) +
                              (G.reads[i].vote[j].g_insert >= MIN_HAPLO_OCCURS) +
                              (G.reads[i].vote[j].t_insert >= MIN_HAPLO_OCCURS));

        int32 ins_total = (G.reads[i].vote[j].a_insert +
                           G.reads[i].vote[j].c_insert +
                           G.reads[i].vote[j].g_insert +
                           G.reads[i].vote[j].t_insert);


        //if  (2 * ins_max > ins_total
        //     && ins_total > 1
        //     && (ins_haplo_ct < 2 || ! Use_Haplo_Ct)
        //     && (G.reads[i].vote[j].no_insert == 0
        //         || (G.reads[i].vote[j].no_insert == 1
        //             && ins_max > 6))) {

        if (2 * ins_max >= ins_total)
          continue;

        if (ins_total <= 1)
          continue;

        if ((ins_haplo_ct >= 2) && (G.Use_Haplo_Ct == true))
          continue;

        if ((G.reads[i].vote[j].no_insert > 0) &&
            ((G.reads[i].vote[j].no_insert != 1) || (ins_max <= 6)))
          continue;

        //  Otherwise, output.

        out.type  = ins_vote;
        out.pos   = j;

        AS_UTL_safeWrite(fp, &out, "correction3", sizeof(Correction_Output_t), 1);
      }
    }
  }

  fclose (fp);
}
