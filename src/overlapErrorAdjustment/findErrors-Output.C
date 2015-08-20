
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "findErrors.H"



void
Output_Details(feParameters *G, uint32 i) {

  fprintf(stderr, ">%d\n", G->bgnID + i);

  for  (uint32 j=0;  G->reads[i].sequence[j] != '\0';  j++)
    fprintf(stderr, "%3d: %c  conf %3d  deletes %3d | subst %3d %3d %3d %3d | no_insert %3d insert %3d %3d %3d %3d\n",
            j,
            j >= G->reads[i].clear_len ? toupper (G->reads[i].sequence[j]) : G->reads[i].sequence[j],
            G->reads[i].vote[j].confirmed,
            G->reads[i].vote[j].deletes,
            G->reads[i].vote[j].a_subst,
            G->reads[i].vote[j].c_subst,
            G->reads[i].vote[j].g_subst,
            G->reads[i].vote[j].t_subst,
            G->reads[i].vote[j].no_insert,
            G->reads[i].vote[j].a_insert,
            G->reads[i].vote[j].c_insert,
            G->reads[i].vote[j].g_insert,
            G->reads[i].vote[j].t_insert);
}


void
Output_Corrections(feParameters *G) {
  Correction_Output_t  out;

  errno = 0;
  FILE *fp = fopen(G->outputFileName, "wb");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", G->outputFileName, strerror(errno)), exit(1);


  for (uint32 i=0; i<G->readsLen; i++) {
    //if (i == 0)
    //  Output_Details(G, i);

    out.keep_left   = (G->reads[i].left_degree  < G->Degree_Threshold);
    out.keep_right  = (G->reads[i].right_degree < G->Degree_Threshold);
    out.type        = IDENT;
    out.pos         = 0;
    out.readID      = G->bgnID + i;

    //fprintf(stderr, "read %d clear_len %d\n", i, G->reads[i].clear_len);
    AS_UTL_safeWrite(fp, &out, "correction1", sizeof(Correction_Output_t), 1);

    if (G->reads[i].sequence == NULL)
      // Deleted fragment
      continue;

    for (uint32 j=0; j<G->reads[i].clear_len; j++) {

      if  (G->reads[i].vote[j].confirmed < 2) {
        Vote_Value_t  vote      = DELETE;
        int32         max       = G->reads[i].vote[j].deletes;
        bool          is_change = true;

        if  (G->reads[i].vote[j].a_subst > max) {
          vote      = A_SUBST;
          max       = G->reads[i].vote[j].a_subst;
          is_change = (G->reads[i].sequence[j] != 'a');
        }

        if  (G->reads[i].vote[j].c_subst > max) {
          vote      = C_SUBST;
          max       = G->reads[i].vote[j].c_subst;
          is_change = (G->reads[i].sequence[j] != 'c');
        }

        if  (G->reads[i].vote[j].g_subst > max) {
          vote      = G_SUBST;
          max       = G->reads[i].vote[j].g_subst;
          is_change = (G->reads[i].sequence[j] != 'g');
        }

        if  (G->reads[i].vote[j].t_subst > max) {
          vote      = T_SUBST;
          max       = G->reads[i].vote[j].t_subst;
          is_change = (G->reads[i].sequence[j] != 't');
        }

        int32 haplo_ct  =  ((G->reads[i].vote[j].deletes >= MIN_HAPLO_OCCURS) +
                            (G->reads[i].vote[j].a_subst >= MIN_HAPLO_OCCURS) +
                            (G->reads[i].vote[j].c_subst >= MIN_HAPLO_OCCURS) +
                            (G->reads[i].vote[j].g_subst >= MIN_HAPLO_OCCURS) +
                            (G->reads[i].vote[j].t_subst >= MIN_HAPLO_OCCURS));

        int32 total  = (G->reads[i].vote[j].deletes +
                        G->reads[i].vote[j].a_subst +
                        G->reads[i].vote[j].c_subst +
                        G->reads[i].vote[j].g_subst +
                        G->reads[i].vote[j].t_subst);

        //  The original had a gargantuajn if test (five clauses, all had to be true) to decide if a record should be output.
        //  It was negated into many small tests if we should skip the output.
        //  A side effect is that we can abort a little earlier in two cases (and we don't even bother).


        //fprintf(stderr, "TEST   read %d position %d type %d -- ", i, j, vote);

        //  (total > 1)
        if (total <= 1) {
          //fprintf(stderr, "FEW   total = %d <= 1\n", total);
          continue;
        }

        //  (2 * max > total)
        if (2 * max <= total) {
          //fprintf(stderr, "WEAK  2*max = %d <= total = %d\n", 2*max, total);
          continue;
        }

        //  (is_change == true)
        if (is_change == false) {
          //fprintf(stderr, "SAME  is_change = %d\n", is_change);
          continue;
        }

        //  ((haplo_ct < 2) || (G->Use_Haplo_Ct == false))
        if ((haplo_ct >= 2) && (G->Use_Haplo_Ct == true)) {
          //fprintf(stderr, "HAPLO haplo_ct=%d >= 2 AND Use_Haplo_Ct = %d\n", haplo_ct, G->Use_Haplo_Ct);
          continue;
        }

        //  ((G->reads[i].vote[j].confirmed == 0) ||
        //   ((G->reads[i].vote[j].confirmed == 1) && (max > 6)))
        if ((G->reads[i].vote[j].confirmed > 0) &&
            ((G->reads[i].vote[j].confirmed != 1) || (max <= 6))) {
          //fprintf(stderr, "INDET confirmed = %d max = %d\n", G->reads[i].vote[j].confirmed, max);
          continue;
        }

        //  Otherwise, output.

        out.type       = vote;
        out.pos        = j;

        //fprintf(stderr, "CORRECT!\n");

        AS_UTL_safeWrite(fp, &out, "correction2", sizeof(Correction_Output_t), 1);
      }  //  confirmed < 2


      if  (G->reads[i].vote[j].no_insert < 2) {
        Vote_Value_t  ins_vote = A_INSERT;
        int32         ins_max  = G->reads[i].vote[j].a_insert;

        if  (ins_max < G->reads[i].vote[j].c_insert) {
          ins_vote = C_INSERT;
          ins_max  = G->reads[i].vote[j].c_insert;
        }

        if  (ins_max < G->reads[i].vote[j].g_insert) {
          ins_vote = G_INSERT;
          ins_max  = G->reads[i].vote[j].g_insert;
        }

        if  (ins_max < G->reads[i].vote[j].t_insert) {
          ins_vote = T_INSERT;
          ins_max  = G->reads[i].vote[j].t_insert;
        }

        int32 ins_haplo_ct = ((G->reads[i].vote[j].a_insert >= MIN_HAPLO_OCCURS) +
                              (G->reads[i].vote[j].c_insert >= MIN_HAPLO_OCCURS) +
                              (G->reads[i].vote[j].g_insert >= MIN_HAPLO_OCCURS) +
                              (G->reads[i].vote[j].t_insert >= MIN_HAPLO_OCCURS));

        int32 ins_total = (G->reads[i].vote[j].a_insert +
                           G->reads[i].vote[j].c_insert +
                           G->reads[i].vote[j].g_insert +
                           G->reads[i].vote[j].t_insert);

        //fprintf(stderr, "TEST   read %d position %d type %d (insert) -- ", i, j, ins_vote);

        if (ins_total <= 1) {
          //fprintf(stderr, "FEW   ins_total = %d <= 1\n", ins_total);
          continue;
        }

        if (2 * ins_max >= ins_total) {
          //fprintf(stderr, "WEAK  2*ins_max = %d <= ins_total = %d\n", 2*ins_max, ins_total);
          continue;
        }

        if ((ins_haplo_ct >= 2) && (G->Use_Haplo_Ct == true)) {
          //fprintf(stderr, "HAPLO ins_haplo_ct=%d >= 2 AND Use_Haplo_Ct = %d\n", ins_haplo_ct, G->Use_Haplo_Ct);
          continue;
        }

        if ((G->reads[i].vote[j].no_insert > 0) &&
            ((G->reads[i].vote[j].no_insert != 1) || (ins_max <= 6))) {
          //fprintf(stderr, "INDET no_insert = %d ins_max = %d\n", G->reads[i].vote[j].no_insert, ins_max);
          continue;
        }

        //  Otherwise, output.

        out.type  = ins_vote;
        out.pos   = j;

        //fprintf(stderr, "INSERT!\n");

        AS_UTL_safeWrite(fp, &out, "correction3", sizeof(Correction_Output_t), 1);
      }  //  insert < 2
    }
  }

  fclose (fp);
}
