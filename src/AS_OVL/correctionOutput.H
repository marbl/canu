

//  Definitions for our exportable data.

typedef enum {
  IDENT,    //  Just an iid in this record.
  DELETE,
  A_SUBST,
  C_SUBST,
  G_SUBST,
  T_SUBST,
  A_INSERT,
  C_INSERT,
  G_INSERT,
  T_INSERT,
  NO_VOTE,
  EXTENSION
} Vote_Value_t;


typedef  union {
  uint32  keep_left   :  1;    // set true if left overlap degree is low
  uint32  keep_right  :  1;    // set true if right overlap degree is low
  uint32  type        :  4;    //  Vote_Value_t
  uint32  pos         :  26;   //
  uint32  readID;
} Correction_Output_t;


