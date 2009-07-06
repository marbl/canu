#ifdef SWIGPERL5
%typemap(in) FILE * {
  $1 = PerlIO_findFILE(IoIFP(sv_2io($input)));
}
#endif

/* Typemap for array of struct types and associated getter which extract
 * the relevant information about the array.
 */
#ifdef SWIGPERL5
%typemap(argout) (void **ary_out, int *elt_size, int *len, swig_type_info **ty) {
  AV *ary;
  SV **svs;
  int i;

  if (argvi >= items) {
    EXTEND(sp,1);              /* Extend the stack by 1 object */
  }
  svs = (SV **) malloc(*$3 * sizeof(SV *));
  for(i = 0; i < *$3; i++) {
    svs[i] = sv_newmortal();
    SWIG_MakePtr(svs[i], (void *)((char *)*$1 + i * *$2), *$4, SWIG_SHADOW);
  }
  ary = av_make(*$3, svs);
  free(svs);
  $result = newRV_noinc((SV *)ary);
  sv_2mortal($result);
  argvi++;
}
%typemap(in,numinputs=0) (void **ary_out, int *elt_size, int *len, swig_type_info **ty) (void *j1, int j2, int j3, swig_type_info *j4) {
  $1 = &j1;
  $2 = &j2;
  $3 = &j3;
  $4 = &j4;
}

/* Typemap for the setters for struct types.
 */
%typemap(in,numinputs=1) (void *objout, void **ary_in, int *elt_size,
                          int len, swig_type_info **ty) (void *j1, int j2,
                                                         swig_type_info *j3) {
  AV *ary;

  $2 = &j1;
  $3 = &j2;
  $5 = &j3;

  $1 = ary = (AV *)SvRV($input);
  $4 = av_len((AV *)ary) + 1;
}
%typemap(argout) (void *objout, void **ary_in, int *elt_size,
                  int len, swig_type_info **ty) {
  AV *ary;
  SV *res, **item;
  int i, ret;
  char *raw;
  void *item_ptr;

  if (argvi >= items) {
    EXTEND(sp,1);              /* Extend the stack by 1 object */
  }

  ary = (AV *)$1;
  res = newSV(*$3 * $4);
  raw = SvPVX(res);
  for(i = 0; i < $4; i++) {
    item = av_fetch(ary, i, 0);
    ret = SWIG_ConvertPtr(*item, &item_ptr, *$5, 0);
    if(!SWIG_IsOK(ret)) {
      sv_2mortal(res);
      SWIG_croak("Can't convert object to appropriate type");
    }
    memcpy(raw + i * *$3, item_ptr, *$3);
  }
  *$2 = raw;
  $result = newRV_noinc((SV *) res);
  sv_2mortal($result);
  argvi++;
}
#endif
%define ARRAYOUT(parent,field,length,type)
     %ignore field;
     %extend parent {
  void swig_ ## field ## _get(void **ary_out, int *elt_size, int *len,
                              swig_type_info **ty) {
    *ary_out = $self-> ## field;
    *elt_size = sizeof(type);
    *len = $self-> ## length;
    *ty = SWIGTYPE_p_ ## type;
  }
  void swig_ ## field ## _set(void *objout, void **ary_in, int *elt_size,
                              int len, swig_type_info **ty) {
    *ary_in = (void *)&($self-> ## field);
    *elt_size = sizeof(type);
    *ty = SWIGTYPE_p_ ## type;
    $self-> ## length = len;
  }
     }
%enddef

/* Special case for SCF messages where the number of contig pairs (here noc
 * field, may be equal to 0, in which case there is still 1 element in array.
 */
%define ARRAYOUT1(parent,field,length,type)
     %ignore field;
     %extend parent {
  void swig_ ## field ## _get(void **ary_out, int *elt_size, int *len,
                              swig_type_info **ty) {
    *ary_out = $self-> ## field;
    *elt_size = sizeof(type);
    *len = ($self-> ## length == 0) ? 1 : $self-> ## length;
    *ty = SWIGTYPE_p_ ## type;
  }
  void swig_ ## field ## _set(void *objout, void **ary_in, int *elt_size,
                              int len, swig_type_info **ty) {
    *ary_in = (void *)&($self-> ## field);
    *elt_size = sizeof(type);
    *ty = SWIGTYPE_p_ ## type;
    $self-> ## length = len;
  }
     }
%enddef

#ifdef SWIGPERL5
/* Typemap for getter of strings (char *)
 */
%typemap(in,numinputs=0) (char **str_out) (char *j1) {
  $1 = &j1;
}
%typemap(argout) (char **str_out) {
  $result = SWIG_FromCharPtr((const char *)*$1);
  argvi++;
}
/* Typemap for setter of strings (char *)
 */
%typemap(in,numinputs=1) (void *objout, char ***str_in) (char **j1) {
  $2 = &j1;
  $1 = $input;
}
%typemap(argout) (void *objout, char ***str_in) {
  SV *res;

  res = newSVsv($1);
  **$2 = SvPV_nolen(res);
  $result = newRV_noinc(res);
  sv_2mortal($result);
  argvi++;
}
#endif
%define STRINGOUT(parent,field)
     %ignore field;
     %extend parent {
  void swig_ ## field ## _get(char **str_out) {
    *str_out = $self-> ## field;
  }
  void swig_ ## field ## _set(void *objout, char ***str_in) {
    *str_in = &($self-> ## field);
  }
     }
%enddef

#ifdef SWIGPERL5
/* Typemap for getter of arrays of strings (char **)
 */
%typemap(argout) (char ***ary_out, int *len) {
  AV *ary;
  SV **svs;
  int i;

  if(argvi >= items) {
    EXTEND(sp,1);
  }
  svs = (SV **) malloc(*$2 * sizeof(SV *));
  for(i = 0; i < *$2; i++) {
    svs[i] = sv_newmortal();
    sv_setpv((SV*)svs[i], (*$1)[i]);
  }
  ary = av_make(*$2, svs);
  free(svs);
  $result = newRV((SV*)ary);
  sv_2mortal($result);
  argvi++;
}
%typemap(in,numinputs=0) (char ***ary_out, int *len) (char **j1, int j2) {
  $1 = &j1;
  $2 = &j2;
}

/* Typemap for setter of arrays of strings (char **)
 */
%typemap(in,numinputs=1) (void *objout, char ****ary_in, int len) (char ***j1) {
  $2 = &j1;

  $1 = (AV *)SvRV($input);
  $3 = av_len((AV *)$1) + 1;
}
%typemap(argout) (void *objout, char ***ary_in, int len) {
  AV *ary;
  SV *res;
  int i, size = 0, len;
  char *raw, *str;

  if(argvi >= items) {
    EXTEND(sp,1);
  }

  ary = (AV *)$1;
  for(i = 0; i < $3; i++) {
    size += 1 + SvCUR(*av_fetch($1, i, 0));
  }
  res = newSV(size + $3 * sizeof(char *));
  raw = SvPVX(res);

  size =  $3 * sizeof(char *);
  for(i = 0; i < $3; i++, size += len) {
    str = SvPV(*av_fetch($1, i, 0), len);
    strncpy(raw + size, str, len)
  }

  *$2 = raw;
  $result = newRV_noinc((SV *)res);
  sv_2mortal($result);
  argvi++
}
#endif
%define STRINGARRAYOUT(parent,field,length)
     %ignore field;
     %extend parent {
  void swig_ ## field ## _get(char ***ary_out, int *len) {
    *ary_out = $self-> ## field;
    *len = $self-> ## length;
  }
  void swig_ ## field ## _set(void *objout, char ****ary_in, int len) {
    *ary_in = &$self-> ## field;
    $self-> ## length = len;
  }
     }
%enddef

/* Macro to defines typemap for arrays of integer type.
 */
%define INTTYPEMAP(type)
     // Typemaps for the getter
%typemap(argout) (type **ary_out, int *len) {
  AV *ary;
  SV **svs;
  int i;

  if(argvi >= items) {
    EXTEND(sp,1);
  }
  svs = (SV **) malloc(*$2 * sizeof(SV *));
  for(i = 0; i < *$2; i++) {
    svs[i] = sv_newmortal();
    sv_setiv((SV*)svs[i], (IV)((*$1)[i]));
  }
  ary = av_make(*$2, svs);
  free(svs);
  $result = newRV((SV*)ary);
  sv_2mortal($result);
  argvi++;
}
%typemap(in,numinputs=0) (type **ary_out, int *len) (type *j1, int j2) {
  $1 = &j1;
  $2 = &j2;
}

// Typemaps for the setter
/* The content of the input Perl array is packed into a Perl string
 * and returned. A reference to the returned string must be kept for
 * as long as the message is kept.
 */
%typemap(argout) (void *obj_out, type *ary_in, int len) {
  if(argvi >= items) {
    EXTEND(sp,1);
  }
  $result = newRV_noinc((SV *)$1);
  sv_2mortal($result);
  argvi++;
}
%typemap(in,numinputs=1) (void *obj_out, type *ary_in, int len) {
  int i, nb_items;
  SV **integer;
  SV *res;
  AV *ary;
  type *raw;

  ary = (AV *)SvRV($input);

  nb_items = av_len((AV *)ary) + 1;
  res = newSV(sizeof(type) * nb_items);
  // raw = (type *)malloc(sizeof(type) * nb_items);
  raw = (type *)SvPVX(res);
  for(i = 0; i < nb_items; i++) {
    integer = av_fetch((AV *)ary, i, 0);
    raw[i] = (type)SvIV(*integer);
  }
  $1 = res;
  $2 = raw;
  $3 = nb_items;
}
%enddef

%define INTARRAYOUT(type,parent,field,length)
     %ignore field;
     %extend parent {
  void swig_ ## field ## _get(type **ary_out, int *len) {
    *ary_out = $self-> ## field;
    *len = $self-> ## length;
  }

  void swig_ ## field ## _set(void *obj_out, type *ary_in, int len) {
    $self-> ## field = ary_in;
    $self-> ## length = len;
  }
     }
%enddef

INTTYPEMAP(int32);
INTTYPEMAP(signed char);

/* Those enum types are really characters
 */
%apply char { LinkType, OrientType, FragType, UnitigType, OverlapType,
                ChunkOrientType, LabelType, UnitigOverlapType, DirectionType,
                UnitigStatus, PlacementStatusType, MateStatType }



/* Load basic int types
 */
#define __extension__
%import stdint.h

/* Define getter for array field of structs
 */
ARRAYOUT(IntConConMesg, pieces, num_pieces, IntMultiPos);
ARRAYOUT(IntConConMesg, unitigs, num_unitigs, IntMultiPos);
ARRAYOUT(IntConConMesg, v_list, num_vars, IntMultiVar);
ARRAYOUT(IntUnitigLinkMesg, jump_list, num_contributing, IntMate_Pairs);
ARRAYOUT(IntContigLinkMesg, jump_list, num_contributing, IntMate_Pairs);
ARRAYOUT(InternalScaffoldLinkMesg, jump_list, num_contributing, IntMate_Pairs);
ARRAYOUT1(IntScaffoldMesg, contig_pairs, num_contig_pairs, IntContigPairs);
ARRAYOUT(SnapUnitigMesg, f_list, num_frags, SnapMultiPos);
ARRAYOUT(SnapUnitigMesg, v_list, num_vars, IntMultiVar);
ARRAYOUT(SnapUnitigLinkMesg, jump_list, num_contributing, SnapMate_Pairs);
ARRAYOUT(SnapConConMesg, pieces, num_pieces, SnapMultiPos);
ARRAYOUT(SnapConConMesg, unitigs, num_unitigs, UnitigPos);
ARRAYOUT(SnapConConMesg, vars, num_vars, IntMultiVar);
ARRAYOUT(SnapContigLinkMesg, jump_list, num_contributing, SnapMate_Pairs);
ARRAYOUT(SnapScaffoldLinkMesg, jump_list, num_contributing, SnapMate_Pairs);
ARRAYOUT1(SnapScaffoldMesg, contig_pairs, num_contig_pairs, SnapContigPairs);
ARRAYOUT(IntUnitigMesg, f_list, num_frags, IntMultiPos);
INTARRAYOUT(int32, SnapMateDistMesg, histogram, num_buckets);
INTARRAYOUT(int32, SnapMultiPos, delta, delta_length);
INTARRAYOUT(int32, IntMateDistMesg, histogram, num_buckets);
INTARRAYOUT(int32, UnitigPos, delta, delta_length);
INTARRAYOUT(int32, IntUnitigPos, delta, delta_length);
INTARRAYOUT(int32, IntMultiPos, delta, delta_length);
INTARRAYOUT(int32, MultiPos, delta, delta_length);
#ifdef AS_MSG_USE_OVL_DELTA
INTARRAYOUT(signed char, OverlapMesg, delta, polymorph_ct);
#endif
STRINGARRAYOUT(LibraryMesg, features, num_features);
STRINGARRAYOUT(LibraryMesg, values, num_features);
/* STRINGOUT(BatchMesg, name); */
/* STRINGOUT(BatchMesg, comment); */
/* STRINGOUT(AuditLine, name); */
/* STRINGOUT(AuditLine, version); */
/* STRINGOUT(AuditLine, comment); */
STRINGOUT(LibraryMesg, source);
STRINGOUT(FragMesg, source);
STRINGOUT(FragMesg, sequence);
STRINGOUT(FragMesg, quality);
STRINGOUT(FragMesg, hps);
//STRINGOUT(UnitigOverlapMesg, source);
//STRINGOUT(ChunkMesg, source);
STRINGOUT(IntMultiVar, nr_conf_alleles);
STRINGOUT(IntMultiVar, weights);
STRINGOUT(IntMultiVar, var_seq);
STRINGOUT(IntMultiVar, conf_read_iids);
//STRINGOUT(IntUnitigMesg, source);
STRINGOUT(IntUnitigMesg, consensus);
STRINGOUT(IntUnitigMesg, quality);
STRINGOUT(UnitigMesg, consensus);
STRINGOUT(UnitigMesg, quality);
STRINGOUT(IntConConMesg, consensus);
STRINGOUT(IntConConMesg, quality);
STRINGOUT(SnapUnitigMesg, consensus);
STRINGOUT(SnapUnitigMesg, quality);
STRINGOUT(SnapConConMesg, consensus);
STRINGOUT(SnapConConMesg, quality);
STRINGOUT(EndOfFileMesgTag, comment);

/* Define __cpluplus, otherwise AS_global.h will define the type bool
 * as an int, on which swig choaks.
 */
#define __cplusplus
%import AS_global.h
%include AS_MSG/AS_MSG_pmesg.h
