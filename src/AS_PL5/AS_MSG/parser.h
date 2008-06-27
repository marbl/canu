/* Defines extension to parser.
 */
#ifdef SWIGPERL5
/* Filter arrays. Translate back & forth between strings ("IUM",
 * "FRG", etc.) and numerical constante.  This is crazy long just for
 * a simple array. Got to be a simpler solution!
 */
%typemap(memberin) MessageType[ANY] {
  int i;
  int len = ($1_dim0 < FILTER_SIZE) ? $1_dim0 : FILTER_SIZE;

  for(i = 0; i < len && $input[i]; i++) {
    $1[i] = $input[i];
  }
  if(i < FILTER_SIZE) {
    $1[i] = 0;
  }
}
%typemap(in) MessageType[ANY] (MessageType j1[$1_dim0]) {
  int i, j, len, type;
  STRLEN str_len;
  AV *ary;
  SV **tv;
  char *name;
  char cname[4];

  if(!SvROK($input))
    croak("Argument $argnum is not a reference.");
  if(SvTYPE(SvRV($input)) != SVt_PVAV)
    croak("Argument $argnum is not an array.");
  ary = (AV *)SvRV($input);
  len = (int)av_len(ary) + 1;
  len = (len < FILTER_SIZE) ? len : FILTER_SIZE;

  cname[3] = 0;
  $1 = j1;
  for(i = 0, j = 0; i < len; i++) {
    tv = av_fetch(ary, i, 0);
    name = SvPV(*tv, str_len);
    strncpy(cname, name, (str_len < 3) ? str_len : 3);
    type = GetMessageType(cname);
    if(type > 0 && type < NUM_OF_REC_TYPES) {
      $1[j] = type;
      j++;
    }
  }
  if(j < FILTER_SIZE) {
    $1[j] = 0;
  }
}
%typemap(out) MessageType[ANY] {
  AV *ary;
  int i, len;

  if(argvi >= items) {
    EXTEND(sp,1);
  }
  len = ($1_dim0 < FILTER_SIZE) ? $1_dim0 : FILTER_SIZE;
  ary = newAV();
  av_extend(ary, len);
  for(i = 0; i < len && $1[i]; i++) {
    av_store(ary, i, newSVpv(GetMessageName($1[i]), 0));
  }
  $result = newRV_noinc((SV*)ary);
  sv_2mortal($result);
  argvi++;
}

/* Typemaps for next_message() && message()
 */
%typemap(argout) (OutMesg *nmesg, void **ty, int *own) {
  SV *sv;
  int flags = SWIG_SHADOW | (*$3 ? SWIG_OWNER : 0);

  if(argvi >= items) {
    EXTEND(sp,1);
  }
  if(*$1 && *$2) {
    sv = sv_newmortal();
    SWIG_MakePtr(sv, *$1, (swig_type_info *)*$2, flags);
  } else {
    sv = &PL_sv_undef;
  }
  $result = sv;
  argvi++;
}
%typemap(in,numinputs=0) (OutMesg *nmesg, void **ty, int *own) \
                         (OutMesg j1, void *j2, int j3) {
  $1 = &j1;
  $2 = &j2;
  $3 = &j3;
}
#endif
%extend Parser {
  void initialize();
  void next_message(OutMesg *nmesg, void **ty, int *own);
  void write_message();
  const char *message_type();
}

%include MSG.h
%perlcode %{
  require AS::MSG::Parser;
    %}
