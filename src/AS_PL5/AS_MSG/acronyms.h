/* 3 letters alias of struct member
 */
#define PACKAGE(mod) package AS::MSG:: ## mod; \
%MEMBERS = (); \
sub EXISTS { my ($self, $key) = @_; return exists $MEMBERS{$key}; } \
sub FIRSTKEY { my ($self) = @_; my $a = keys %MEMBERS; each %MEMBERS; } \
sub NEXTKEY { my ($self, $key) = @_; each %MEMBERS; } \
sub SCALAR { return 1; } \
$TYPE = AS::MSG::GetMessageTypeFromStruct( #mod)

#define GETTERALIAS(mod,new,old) *swig_ ## new ## _get = *swig_ ## old ## _get; \
*swig_ ## new ## _set = *swig_ ## old ## _set; \
$MEMBERS{ #new } = #old

%perlcode %{
  package AS::MSG;
  use Scalar::Util qw(blessed);

  sub members {
    my $self = shift;
    my $pack = blessed($self);

    no strict 'refs';
    return \%{"${pack}::MEMBERS"};
  }
%}

%perlcode {
  PACKAGE(LibraryMesg);
  GETTERALIAS(LibraryMesg, act, action);
  GETTERALIAS(LibraryMesg, acc, eaccession);
  GETTERALIAS(LibraryMesg, mea, mean);
  GETTERALIAS(LibraryMesg, std, stddev);
  GETTERALIAS(LibraryMesg, src, source);
  GETTERALIAS(LibraryMesg, ori, link_orient);
  GETTERALIAS(LibraryMesg, nft, num_features);
  GETTERALIAS(LibraryMesg, fea, features);
  GETTERALIAS(LibraryMesg, val, values);

  PACKAGE(FragMesg);
  GETTERALIAS(FragMesg, act, action);
  GETTERALIAS(FragMesg, acc, eaccession);
  GETTERALIAS(FragMesg, acc1, eaccession);
  GETTERALIAS(FragMesg, acc2, iaccession);
  GETTERALIAS(FragMesg, typ, type);
  GETTERALIAS(FragMesg, rnd, is_random);
  GETTERALIAS(FragMesg, sta, status_code);
  GETTERALIAS(FragMesg, lib, library_uid);
  GETTERALIAS(FragMesg, pla, plate_uid);
  GETTERALIAS(FragMesg, loc, plate_location);
  GETTERALIAS(FragMesg, src, source);
  GETTERALIAS(FragMesg, seq, sequence);
  GETTERALIAS(FragMesg, qlt, quality);
  GETTERALIAS(FragMesg, hps, hps);
  GETTERALIAS(FragMesg, clr, clear_rng);
  GETTERALIAS(FragMesg, clv, clear_vec);
  GETTERALIAS(FragMesg, clq, clear_qlt);

  PACKAGE(LinkMesg);
  GETTERALIAS(LinkMesg, act, action);
  GETTERALIAS(LinkMesg, typ, type);
  GETTERALIAS(LinkMesg, fg1, frag1);
  GETTERALIAS(LinkMesg, fg2, frag2);
  GETTERALIAS(LinkMesg, dst, distance);
  GETTERALIAS(LinkMesg, ori, link_orient);

  PACKAGE(DistanceMesg);
  GETTERALIAS(DistanceMesg, act, action);
  GETTERALIAS(DistanceMesg, acc, eaccession);
  GETTERALIAS(DistanceMesg, mea, mean);
  GETTERALIAS(DistanceMesg, std, stddev);

  PACKAGE(AuditLine);
  GETTERALIAS(AuditLine, who, name);
  GETTERALIAS(AuditLine, ctm, complete);
  GETTERALIAS(AuditLine, vsn, version);
  GETTERALIAS(AuditLine, com, comment);

  PACKAGE(VersionMesg);
  GETTERALIAS(VersionMesg, ver, version);

  PACKAGE(OverlapMesg);
  GETTERALIAS(OverlapMesg, afr, aifrag);
  GETTERALIAS(OverlapMesg, bfr, bifrag);
  GETTERALIAS(OverlapMesg, ori, orientation);
  GETTERALIAS(OverlapMesg, olt, overlap_type);
  GETTERALIAS(OverlapMesg, ahg, ahg);
  GETTERALIAS(OverlapMesg, bhg, bhg);
  GETTERALIAS(OverlapMesg, qua, quality);
  GETTERALIAS(OverlapMesg, mno, min_offset);
  GETTERALIAS(OverlapMesg, mxo, max_offset);
  GETTERALIAS(OverlapMesg, pct, polymorph_ct);
#ifdef AS_MSG_USE_OVL_DELTA
  GETTERALIAS(OverlapMesg, del, alignment_delta);
#endif

  PACKAGE(UnitigOverlapMesg);
  GETTERALIAS(UnitigOverlapMesg, ck1, chunk1);
  GETTERALIAS(UnitigOverlapMesg, ck2, chunk2);
  GETTERALIAS(UnitigOverlapMesg, ori, orient);
  GETTERALIAS(UnitigOverlapMesg, ovt, overlap_type);
  //GETTERALIAS(UnitigOverlapMesg, src, source);
  GETTERALIAS(UnitigOverlapMesg, len, best_overlap_length);
  GETTERALIAS(UnitigOverlapMesg, min, min_overlap_length);
  GETTERALIAS(UnitigOverlapMesg, max, max_overlap_length);
  GETTERALIAS(UnitigOverlapMesg, qua, quality);

  PACKAGE(IntMultiPos);
  GETTERALIAS(IntMultiPos, typ, type);
  GETTERALIAS(IntMultiPos, mid, ident);
  GETTERALIAS(IntMultiPos, con, contained);
  GETTERALIAS(IntMultiPos, bid, ident2);
  GETTERALIAS(IntMultiPos, pos, position);
  GETTERALIAS(IntMultiPos, ahg, ahang);
  GETTERALIAS(IntMultiPos, bhg, bhang);
  GETTERALIAS(IntMultiPos, dln, delta_length);
  GETTERALIAS(IntMultiPos, del, delta);

  PACKAGE(IntMultiVar);
  GETTERALIAS(IntMultiVar, pos, position);
  GETTERALIAS(IntMultiVar, nrd, num_reads);
  GETTERALIAS(IntMultiVar, nca, num_conf_alleles);
  GETTERALIAS(IntMultiVar, anc, anchor_size);
  GETTERALIAS(IntMultiVar, len, var_length);
  GETTERALIAS(IntMultiVar, vid, curr_var_id);
  GETTERALIAS(IntMultiVar, pid, phased_var_id);
  GETTERALIAS(IntMultiVar, nra, nr_conf_alleles);
  GETTERALIAS(IntMultiVar, wgt, weights);
  GETTERALIAS(IntMultiVar, seq, var_seq);
  GETTERALIAS(IntMultiVar, rid, conf_read_iids);

  PACKAGE(IntUnitgPos);
  GETTERALIAS(IntUnitgPos, typ, type);
  GETTERALIAS(IntUnitgPos, lid, ident);
  GETTERALIAS(IntUnitgPos, pos, position);
  GETTERALIAS(IntUnitgPos, dln, delta_length);
  GETTERALIAS(IntUnitgPos, del, delta);

  PACKAGE(IntUnitigMesg);
  GETTERALIAS(IntUnitigMesg, acc, iaccession);
  //GETTERALIAS(IntUnitigMesg, src, source);
  GETTERALIAS(IntUnitigMesg, cov, coverage_stat);
  GETTERALIAS(IntUnitigMesg, sta, status);
  GETTERALIAS(IntUnitigMesg, len, length);
  GETTERALIAS(IntUnitigMesg, cns, consensus);
  GETTERALIAS(IntUnitigMesg, qlt, quality);
  GETTERALIAS(IntUnitigMesg, for, forced);
  GETTERALIAS(IntUnitigMesg, nfr, num_frags);
  GETTERALIAS(IntUnitigMesg, IMPs, f_list);


  PACKAGE(IntUnitigLinkMesg);
  GETTERALIAS(IntUnitigLinkMesg, ut1, unitig1);
  GETTERALIAS(IntUnitigLinkMesg, ut2, unitig2);
  GETTERALIAS(IntUnitigLinkMesg, ori, orientation);
  GETTERALIAS(IntUnitigLinkMesg, ovt, overlap_type);
  GETTERALIAS(IntUnitigLinkMesg, ipc, is_possible_chimera);
  GETTERALIAS(IntUnitigLinkMesg, gui, includes_guide);
  GETTERALIAS(IntUnitigLinkMesg, mea, mean_distance);
  GETTERALIAS(IntUnitigLinkMesg, std, std_deviation);
  GETTERALIAS(IntUnitigLinkMesg, num, num_contributing);
  GETTERALIAS(IntUnitigLinkMesg, jls, jump_list);

  PACKAGE(IntContigLinkMesg);
  GETTERALIAS(IntContigLinkMesg, co1, contig1);
  GETTERALIAS(IntContigLinkMesg, co2, contig2);
  GETTERALIAS(IntContigLinkMesg, ori, orientation);
  GETTERALIAS(IntContigLinkMesg, ovt, overlap_type);
  GETTERALIAS(IntContigLinkMesg, ipc, is_possible_chimera);
  GETTERALIAS(IntContigLinkMesg, gui, includes_guide);
  GETTERALIAS(IntContigLinkMesg, mea, mean_distance);
  GETTERALIAS(IntContigLinkMesg, std, std_deviation);
  GETTERALIAS(IntContigLinkMesg, num, num_contributing);
  GETTERALIAS(IntContigLinkMesg, sta, status);
  GETTERALIAS(IntContigLinkMesg, jls, jump_list);

  PACKAGE(InternalScaffoldLinkMesg);
  GETTERALIAS(InternalScaffoldLinkMesg, sc1, iscaffold1);
  GETTERALIAS(InternalScaffoldLinkMesg, sc2, iscaffold2);
  GETTERALIAS(InternalScaffoldLinkMesg, ori, orientation);
  GETTERALIAS(InternalScaffoldLinkMesg, gui, includes_guide);
  GETTERALIAS(InternalScaffoldLinkMesg, mea, mean_distance);
  GETTERALIAS(InternalScaffoldLinkMesg, std, std_deviation);
  GETTERALIAS(InternalScaffoldLinkMesg, num, num_contributing);
  GETTERALIAS(InternalScaffoldLinkMesg, jls, jump_list);

  PACKAGE(AugFragMesg);
  GETTERALIAS(AugFragMesg, acc1, eaccession);
  GETTERALIAS(AugFragMesg, acc2, iaccession);
  GETTERALIAS(AugFragMesg, mst, mate_status);
  GETTERALIAS(AugFragMesg, chi, chimeric);
  GETTERALIAS(AugFragMesg, cha, chaff);
  GETTERALIAS(AugFragMesg, clr, clear_rng);

  PACKAGE(IntContigPairs);
  GETTERALIAS(IntContigPairs, ct1, contig1);
  GETTERALIAS(IntContigPairs, ct2, contig2);
  GETTERALIAS(IntContigPairs, mea, mean);
  GETTERALIAS(IntContigPairs, std, stddev);
  GETTERALIAS(IntContigPairs, ori, orient);

  PACKAGE(IntScaffoldMesg);
  GETTERALIAS(IntScaffoldMesg, acc, iaccession);
  GETTERALIAS(IntScaffoldMesg, noc, num_contig_pairs);
  GETTERALIAS(IntScaffoldMesg, ICPs, contig_pairs);

  PACKAGE(IntMateDistMesg);
  GETTERALIAS(IntMateDistMesg, ref, refines);
  GETTERALIAS(IntMateDistMesg, mea, mean);
  GETTERALIAS(IntMateDistMesg, std, stddev);
  GETTERALIAS(IntMateDistMesg, min, min);
  GETTERALIAS(IntMateDistMesg, max, max);
  GETTERALIAS(IntMateDistMesg, buc, num_buckets);
  GETTERALIAS(IntMateDistMesg, his, histogram);

  PACKAGE(IntConConMesg);
  GETTERALIAS(IntConConMesg, acc, iaccession);
  GETTERALIAS(IntConConMesg, pla, placed);
  GETTERALIAS(IntConConMesg, len, length);
  GETTERALIAS(IntConConMesg, cns, consensus);
  GETTERALIAS(IntConConMesg, qlt, quality);
  GETTERALIAS(IntConConMesg, for, forced);
  GETTERALIAS(IntConConMesg, npc, numpieces);
  GETTERALIAS(IntConConMesg, nou, num_unitigs);
  GETTERALIAS(IntConConMesg, nvr, num_vars);
  GETTERALIAS(IntConConMesg, IMVs, v_list);
  GETTERALIAS(IntConConMesg, IMPs, pieces);
  GETTERALIAS(IntConConMesg, IUPs, unitigs);


  PACKAGE(IntAugFragMesg);
  GETTERALIAS(IntAugFragMesg, acc, iaccession);
  GETTERALIAS(IntAugFragMesg, typ, type);
  GETTERALIAS(IntAugFragMesg, chi, chimeric);
  GETTERALIAS(IntAugFragMesg, cha, chaff);
  GETTERALIAS(IntAugFragMesg, clr, clear_rng);
  GETTERALIAS(IntAugFragMesg, mst, mate_status);

  PACKAGE(UnitigPos);
  GETTERALIAS(UnitigPos, typ, type);
  GETTERALIAS(UnitigPos, lid, eident);
  GETTERALIAS(UnitigPos, pos, position);
  GETTERALIAS(UnitigPos, dln, delta_length);
  GETTERALIAS(UnitigPos, del, delta);

  PACKAGE(SnapMultiPos);
  GETTERALIAS(SnapMultiPos, typ, type);
  GETTERALIAS(SnapMultiPos, mid, eident);
  GETTERALIAS(SnapMultiPos, src, source);
  GETTERALIAS(SnapMultiPos, pos, position);
  GETTERALIAS(SnapMultiPos, dln, delta_length);
  GETTERALIAS(SnapMultiPos, del, delta);

  PACKAGE(SnapUnitigMesg);
  GETTERALIAS(SnapUnitigMesg, acc1, eaccession);
  GETTERALIAS(SnapUnitigMesg, acc2, iaccession);
  GETTERALIAS(SnapUnitigMesg, src, source);
  GETTERALIAS(SnapUnitigMesg, cov, coverage_stat);
  GETTERALIAS(SnapUnitigMesg, sta, status);
  GETTERALIAS(SnapUnitigMesg, len, length);
  GETTERALIAS(SnapUnitigMesg, cns, consensus);
  GETTERALIAS(SnapUnitigMesg, qlt, quality);
  GETTERALIAS(SnapUnitigMesg, for, forced);
  GETTERALIAS(SnapUnitigMesg, nfr, num_frags);
  GETTERALIAS(SnapUnitigMesg, MPSs, f_list);

  PACKAGE(SnapUnitigLinkMesg);
  GETTERALIAS(SnapUnitigLinkMesg, ut1, eunitig1);
  GETTERALIAS(SnapUnitigLinkMesg, ut2, eunitig2);
  GETTERALIAS(SnapUnitigLinkMesg, ori, orientation);
  GETTERALIAS(SnapUnitigLinkMesg, ovt, overlap_type);
  GETTERALIAS(SnapUnitigLinkMesg, ipc, is_possible_chimera);
  GETTERALIAS(SnapUnitigLinkMesg, gui, includes_guide);
  GETTERALIAS(SnapUnitigLinkMesg, mea, mean_distance);
  GETTERALIAS(SnapUnitigLinkMesg, std, std_deviation);
  GETTERALIAS(SnapUnitigLinkMesg, num, num_contributing);
  GETTERALIAS(SnapUnitigLinkMesg, sta, status);
  GETTERALIAS(SnapUnitigLinkMesg, jls, jump_list);

  PACKAGE(SnapConConMesg);
  GETTERALIAS(SnapConConMesg, acc1, eaccession);
  GETTERALIAS(SnapConConMesg, acc2, iaccession);
  GETTERALIAS(SnapConConMesg, pla, placed);
  GETTERALIAS(SnapConConMesg, len, length);
  GETTERALIAS(SnapConConMesg, cns, consensus);
  GETTERALIAS(SnapConConMesg, qlt, quality);
  GETTERALIAS(SnapConConMesg, for, forced);
  GETTERALIAS(SnapConConMesg, npc, num_pieces);
  GETTERALIAS(SnapConConMesg, nou, num_unitigs);
  GETTERALIAS(SnapConConMesg, nvr, num_vars);
  GETTERALIAS(SnapConConMesg, VARs, vars);
  GETTERALIAS(SnapConConMesg, MPSs, pieces);
  GETTERALIAS(SnapConConMesg, UPSs, unitigs);

  PACKAGE(SnapContigLinkMesg);
  GETTERALIAS(SnapContigLinkMesg, co1, econtig1);
  GETTERALIAS(SnapContigLinkMesg, co2, econtig2);
  GETTERALIAS(SnapContigLinkMesg, ori, orientation);
  GETTERALIAS(SnapContigLinkMesg, ovt, overlap_type);
  GETTERALIAS(SnapContigLinkMesg, ipc, is_possible_chimera);
  GETTERALIAS(SnapContigLinkMesg, gui, includes_guide);
  GETTERALIAS(SnapContigLinkMesg, mea, mean_distance);
  GETTERALIAS(SnapContigLinkMesg, std, std_deviation);
  GETTERALIAS(SnapContigLinkMesg, num, num_contributing);
  GETTERALIAS(SnapContigLinkMesg, sta, status);
  GETTERALIAS(SnapContigLinkMesg, jls, jump_list);

  PACKAGE(SnapScaffoldLinkMesg);
  GETTERALIAS(SnapScaffoldLinkMesg, sc1, escaffold1);
  GETTERALIAS(SnapScaffoldLinkMesg, sc2, escaffold2);
  GETTERALIAS(SnapScaffoldLinkMesg, ori, orientation);
  GETTERALIAS(SnapScaffoldLinkMesg, gui, includes_guide);
  GETTERALIAS(SnapScaffoldLinkMesg, mea, mean_distance);
  GETTERALIAS(SnapScaffoldLinkMesg, std, std_deviation);
  GETTERALIAS(SnapScaffoldLinkMesg, num, num_contributing);
  GETTERALIAS(SnapScaffoldLinkMesg, jls, jump_list);

  PACKAGE(SnapContigPairs);
  GETTERALIAS(SnapContigPairs, ct1, econtig1);
  GETTERALIAS(SnapContigPairs, ct2, econtig2);
  GETTERALIAS(SnapContigPairs, mea, mean);
  GETTERALIAS(SnapContigPairs, std, stddev);
  GETTERALIAS(SnapContigPairs, ori, orient);

  PACKAGE(SnapScaffoldMesg);
  GETTERALIAS(SnapScaffoldMesg, acc1, eaccession);
  GETTERALIAS(SnapScaffoldMesg, acc2, iaccession);
  GETTERALIAS(SnapScaffoldMesg, noc, num_contig_pairs);
  GETTERALIAS(SnapScaffoldMesg, CTPs, contig_pairs);

  PACKAGE(SnapDegenerateScaffoldMesg);
  GETTERALIAS(SnapDegenerateScaffoldMesg, acc, eaccession);
  GETTERALIAS(SnapDegenerateScaffoldMesg, ctg, econtig);

  PACKAGE(IntDegenerateScaffoldMesg);
  GETTERALIAS(IntDegenerateScaffoldMesg, ctg, icontig);

  PACKAGE(SnapMateDistMesg);
  GETTERALIAS(SnapMateDistMesg, ref1, erefines);
  GETTERALIAS(SnapMateDistMesg, ref2, irefines);
  GETTERALIAS(SnapMateDistMesg, mea, mean);
  GETTERALIAS(SnapMateDistMesg, std, stddev);
  GETTERALIAS(SnapMateDistMesg, min, min);
  GETTERALIAS(SnapMateDistMesg, max, max);
  GETTERALIAS(SnapMateDistMesg, buc, num_buckets);
  GETTERALIAS(SnapMateDistMesg, his, histogram);

  PACKAGE(BatchMesg);
  GETTERALIAS(BatchMesg, bna, name);
  GETTERALIAS(BatchMesg, acc, eaccession);
  GETTERALIAS(BatchMesg, com, comment);

  PACKAGE(EndOfFileMesg);
  GETTERALIAS(EndOfFileMesg, sta, status);
  GETTERALIAS(EndOfFileMesg, com, comment);
}
