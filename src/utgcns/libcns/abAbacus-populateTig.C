

//  Save an abVarRegion into a tgTig tgVariantPosition.

static
int
GetTheMostDistantRead(int32 curr_read_id,
                      int32 nr,
                      int32 **dist_matrix) {
  int32 i, dist_read_id = -1;
  int32 max_dist = -1;

  for (i=0; i<nr; i++) {
    if (i == curr_read_id)
      continue;

    if (max_dist < dist_matrix[curr_read_id][i]) {
      max_dist = dist_matrix[curr_read_id][i];
      dist_read_id = i;
    }
  }

  return dist_read_id;
}



static
void
PopulateVARRecord(bool             is_phased,
                  uint32          *cids,
                  uint32          &vn,
                  uint32          &vm,
                  abVarRead      *&v,
                  abVarRegion      vreg,
                  CNS_Options     *opp,
                  uint32           get_scores,
                  uint32          *conf_read_iids) {

  if (v == NULL)
    v = new abVarRead [vm];

  if (vn == vm)
    resizeArray(v, vn, vm, vm + 10);

  vreg_id++;

  //  v is the tgVariantPosition, see IntMultiVar in AS_MSG_pmesg.H

  v[vn].var_id                 = vreg_id;
  v[vn].phased_id              = is_phased ? vreg_id - 1 : -1;
  v[vn].position.bgn           = vreg.beg;
  v[vn].position.end           = vreg.end;
  v[vn].num_reads              = vreg.nr;
  v[vn].num_alleles            = (vreg.nca < 2) ? 2 : vreg.nca;
  v[vn].num_alleles_confirmed  = vreg.nca;
  v[vn].min_anchor_size        = opp->smooth_win;
  v[vn].var_length             = vreg.end - vreg.beg;

  v[vn].alleles        = new IntVarAllele [v[vn].num_alleles];
  v[vn].var_seq_memory = new char         [v[vn].num_alleles];
  v[vn].read_id_memory = new int32        [v[vn].num_reads];

  v[vn].enc_num_reads = NULL;
  v[vn].enc_weights   = NULL;
  v[vn].enc_var_seq   = NULL;
  v[vn].enc_read_ids  = NULL;

  int32 distant_read_id   = -42;
  int32 distant_allele_id = -42;

  if (vreg.nca < 2) {
    distant_read_id   = GetTheMostDistantRead(vreg.alleles[0].read_ids[0], vreg.nr, vreg.dist_matrix);
    distant_allele_id = vreg.reads[distant_read_id].allele_id;
  }

  char *base  = new char [v[vn].num_alleles];
  int32 shift = vreg.end - vreg.beg + 1;

  for (int32 m=0; m<vreg.end - vreg.beg; m++) {
    memset(base, 0, sizeof(char) * v[vn].num_alleles);

    for (int32 al=v[vn].num_alleles-1; al >=0; al--) {
      if ((al == 0) || (al < vreg.nca)) {
        int32 read_id = vreg.alleles[al].read_ids[0];

        base[al] = vreg.reads[read_id].bases[m];

        if (al == 0) {
          int32   cid      = cids[vreg.beg+m];
          Column *column   = GetColumn(columnStore,cid);
          Bead   *call     = GetBead(beadStore, column->call);
          double  fict_var = 0;
          char    cbase    = 0;

          // Set the consensus quality and base
          BaseCall(cid, 1, fict_var, &vreg, -1, cbase, 0, opp);
          Setchar(sequenceStore, call->soffset, &base[al]);
        }
      } else {
        // vreg.nca < 2 and al == 1
        base[al] = vreg.reads[distant_read_id].bases[m];
      }

      v[vn].var_seq_memory[m + al * shift] = base[al];
    }

    if (get_scores > 0)
      for (uint32 i=0; i<v[vn].num_alleles; i++)
        for (uint32 j=i+1; j<v[vn].num_alleles; j++)
          if (base[i] != base[j])
            NumAAMismatches++;
  }

  delete [] base;

  int32  vso = 0;
  int32  rio = 0;

  for (int32 al=0; al < v[vn].num_alleles; al++) {
    v[vn].alleles[al].num_reads      = vreg.alleles[al].num_reads;
    v[vn].alleles[al].weight         = vreg.alleles[al].weight;
    v[vn].alleles[al].var_seq_offset = vso;
    v[vn].alleles[al].read_id_offset = rio;

    vso += v[vn].var_length + 1;
    rio += v[vn].alleles[al].num_reads;
  }

  int32   tot_num_conf_reads = 0;

  for (int32 m=0; m<vreg.nca; m++)
    tot_num_conf_reads += vreg.alleles[m].num_reads;

  for (int32 rd=0; rd < tot_num_conf_reads; rd++)
    v[vn].read_id_memory[rd] = conf_read_iids[rd];

  for (int32 al=0; al < v[vn].num_alleles; al++)
    if (v[vn].var_seq_memory[al * shift]             == '-' &&
        v[vn].var_seq_memory[al * shift + shift - 2] == '-')
      NumVARStringsWithFlankingGaps++;

  vn++;
}

