#include "prefixEditDistance.H"

#include "Binomial_Bound.H"

extern const int32 Edit_Match_Limit_0100[];
extern const int32 Edit_Match_Limit_0200[];
extern const int32 Edit_Match_Limit_0300[];
extern const int32 Edit_Match_Limit_0400[];
extern const int32 Edit_Match_Limit_0500[];
extern const int32 Edit_Match_Limit_0600[];
extern const int32 Edit_Match_Limit_0700[];
extern const int32 Edit_Match_Limit_0800[];
extern const int32 Edit_Match_Limit_0900[];
extern const int32 Edit_Match_Limit_1000[];
extern const int32 Edit_Match_Limit_1100[];
extern const int32 Edit_Match_Limit_1200[];
extern const int32 Edit_Match_Limit_1300[];
extern const int32 Edit_Match_Limit_1400[];
extern const int32 Edit_Match_Limit_1500[];
extern const int32 Edit_Match_Limit_1600[];
extern const int32 Edit_Match_Limit_1700[];
extern const int32 Edit_Match_Limit_1800[];
extern const int32 Edit_Match_Limit_1900[];
extern const int32 Edit_Match_Limit_2000[];
extern const int32 Edit_Match_Limit_2100[];
extern const int32 Edit_Match_Limit_2200[];
extern const int32 Edit_Match_Limit_2300[];
extern const int32 Edit_Match_Limit_2400[];
extern const int32 Edit_Match_Limit_2500[];
extern const int32 Edit_Match_Limit_2600[];
extern const int32 Edit_Match_Limit_2700[];
extern const int32 Edit_Match_Limit_2800[];
extern const int32 Edit_Match_Limit_2900[];
extern const int32 Edit_Match_Limit_3000[];
extern const int32 Edit_Match_Limit_3100[];
extern const int32 Edit_Match_Limit_3200[];
extern const int32 Edit_Match_Limit_3300[];
extern const int32 Edit_Match_Limit_3400[];
extern const int32 Edit_Match_Limit_3500[];
extern const int32 Edit_Match_Limit_3600[];
extern const int32 Edit_Match_Limit_3700[];
extern const int32 Edit_Match_Limit_3800[];
extern const int32 Edit_Match_Limit_3900[];
extern const int32 Edit_Match_Limit_4000[];
extern const int32 Edit_Match_Limit_4100[];
extern const int32 Edit_Match_Limit_4200[];
extern const int32 Edit_Match_Limit_4300[];
extern const int32 Edit_Match_Limit_4400[];
extern const int32 Edit_Match_Limit_4500[];
extern const int32 Edit_Match_Limit_4600[];
extern const int32 Edit_Match_Limit_4700[];
extern const int32 Edit_Match_Limit_4800[];
extern const int32 Edit_Match_Limit_4900[];
extern const int32 Edit_Match_Limit_5000[];


static
const
int32 *
Edit_Match_Limit_Data[50] = {
  Edit_Match_Limit_0100,
  Edit_Match_Limit_0200,
  Edit_Match_Limit_0300,
  Edit_Match_Limit_0400,
  Edit_Match_Limit_0500,
  Edit_Match_Limit_0600,
  Edit_Match_Limit_0700,
  Edit_Match_Limit_0800,
  Edit_Match_Limit_0900,
  Edit_Match_Limit_1000,
  Edit_Match_Limit_1100,
  Edit_Match_Limit_1200,
  Edit_Match_Limit_1300,
  Edit_Match_Limit_1400,
  Edit_Match_Limit_1500,
  Edit_Match_Limit_1600,
  Edit_Match_Limit_1700,
  Edit_Match_Limit_1800,
  Edit_Match_Limit_1900,
  Edit_Match_Limit_2000,
  Edit_Match_Limit_2100,
  Edit_Match_Limit_2200,
  Edit_Match_Limit_2300,
  Edit_Match_Limit_2400,
  Edit_Match_Limit_2500,
  Edit_Match_Limit_2600,
  Edit_Match_Limit_2700,
  Edit_Match_Limit_2800,
  Edit_Match_Limit_2900,
  Edit_Match_Limit_3000,
  Edit_Match_Limit_3100,
  Edit_Match_Limit_3200,
  Edit_Match_Limit_3300,
  Edit_Match_Limit_3400,
  Edit_Match_Limit_3500,
  Edit_Match_Limit_3600,
  Edit_Match_Limit_3700,
  Edit_Match_Limit_3800,
  Edit_Match_Limit_3900,
  Edit_Match_Limit_4000,
  Edit_Match_Limit_4100,
  Edit_Match_Limit_4200,
  Edit_Match_Limit_4300,
  Edit_Match_Limit_4400,
  Edit_Match_Limit_4500,
  Edit_Match_Limit_4600,
  Edit_Match_Limit_4700,
  Edit_Match_Limit_4800,
  Edit_Match_Limit_4900,
  Edit_Match_Limit_5000
};


const char *
toString(pedAlignType at) {
  const char *ret = NULL;

  switch (at) {
    case pedOverlap:  ret = "full-overlap";      break;
    case pedLocal:    ret = "local-overlap";     break;
    case pedGlobal:   ret = "global-alignment";  break;
    default:
      assert(0);
      break;
  }

  return(ret);
}


const char *
toString(pedOverlapType ot) {
  const char *ret = NULL;

  switch (ot) {
    case pedBothBranch:   ret = "both-branch";       break;
    case pedLeftBranch:   ret = "left-branch";       break;
    case pedRightBranch:  ret = "right-branch";      break;
    case pedDovetail:     ret = "dovetail-overlap";  break;
    default:
      assert(0);
      break;
  }

  return(ret);
}


prefixEditDistance::prefixEditDistance(pedAlignType alignType_, double maxErate_) {
  alignType            = alignType_;
  maxErate             = maxErate_;

  MAX_ERRORS             = (1 + (int32)ceil(maxErate * AS_MAX_READLEN));
  ERRORS_FOR_FREE        = 1;
  MIN_BRANCH_END_DIST    = 20;
  MIN_BRANCH_TAIL_SLOPE  = ((maxErate > 0.06) ? 1.0 : 0.20);

  Left_Delta  = new int32  [MAX_ERRORS];
  Right_Delta = new int32  [MAX_ERRORS];

  allocated += 3 * MAX_ERRORS * sizeof(int32);

  Delta_Stack = new int32  [MAX_ERRORS];

  Edit_Space_Lazy = new pedEdit *  [MAX_ERRORS];
  Edit_Array_Lazy = new pedEdit *  [MAX_ERRORS];

  memset(Edit_Space_Lazy, 0, sizeof(pedEdit *) * MAX_ERRORS);
  memset(Edit_Array_Lazy, 0, sizeof(pedEdit *) * MAX_ERRORS);

  allocated += MAX_ERRORS * sizeof (pedEdit);
  allocated += MAX_ERRORS * sizeof (pedEdit);

  int32   dataIndex = (int)ceil(maxErate * 100) - 1;

  if ((dataIndex < 0) || (50 <= dataIndex))
    fprintf(stderr, "prefixEditDistance()-- Invalid maxErate=%f -> dataIndex=%d\n",
            maxErate, dataIndex);
  assert(0 <= dataIndex);
  assert(dataIndex < 50);

#if 0

  //  Use the precomputed values.

  Edit_Match_Limit_Allocation = NULL;
  Edit_Match_Limit            = Edit_Match_Limit_Data[dataIndex];

  fprintf(stderr, "prefixEditDistance()-- Set Edit_Match_Limit to %p; dataIndex=%d 6 = %p\n",
          Edit_Match_Limit, dataIndex, Edit_Match_Limit_0600);

#else

  //  Compute values on the fly.

  Edit_Match_Limit_Allocation = new int32 [MAX_ERRORS + 1];

  for (int32 e=0;  e<= ERRORS_FOR_FREE; e++)
    Edit_Match_Limit_Allocation[e] = 0;

  int Start = 1;

  for (int32 e=ERRORS_FOR_FREE + 1; e<MAX_ERRORS; e++) {
    Start = Binomial_Bound(e - ERRORS_FOR_FREE,
                           maxErate,
                           Start);
    Edit_Match_Limit_Allocation[e] = Start - 1;

    assert(Edit_Match_Limit_Allocation[e] >= Edit_Match_Limit_Allocation[e-1]);
  }

  Edit_Match_Limit = Edit_Match_Limit_Allocation;

#endif

  for (int32 i=0; i <= AS_MAX_READLEN; i++) {
    //Error_Bound[i] = (int32) (i * maxErate + 0.0000000000001);
    Error_Bound[i] = (int32)ceil(i * maxErate);
  }


  //  Value to add for a match in finding branch points.
  //
  //  ALH: Note that maxErate also affects what overlaps get found
  //
  //  ALH: Scoring seems to be unusual: given an alignment of length l with k mismatches, the score
  //  seems to be computed as l + k * error value and NOT (l-k)*match+k*error
  //
  //  I.e. letting x := DEFAULT_BRANCH_MATCH_VAL, the max mismatch fraction p to give a non-negative
  //  score would be p = x/(1-x); conversely, to compute x for a goal p, we have x = p/(1+p).  E.g.
  //
  //  for p=0.06, x = .06 / (1.06) = .0566038
  //  for p=0.35, x = .35 / (1.35) = .259259
  //  for p=0.2,  x = .20 / (1.20) = .166667
  //  for p=0.15, x = .15 / (1.15) = .130435
  //
  //  Value was for 6% vs 35% error discrimination.
  //
  //  Converting to integers didn't make it faster.
  //
  //  Corresponding error value is this value minus 1.0

  Branch_Match_Value = maxErate / (1 + maxErate);

  for (uint32 ii=0; ii<256; ii++)
    tolower[ii] = ::tolower(ii);
};



prefixEditDistance::~prefixEditDistance() {
  delete [] Left_Delta;
  delete [] Right_Delta;

  delete [] Delta_Stack;

  for (uint32 i=0; i<MAX_ERRORS; i++)
    if (Edit_Space_Lazy[i])
      delete [] Edit_Space_Lazy[i];

  delete [] Edit_Space_Lazy;
  delete [] Edit_Array_Lazy;

  delete [] Edit_Match_Limit_Allocation;
};

