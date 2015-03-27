#include "AS_global.H"
#include "overlapInCore.H"

#include "Binomial_Bound.H"  //  liboverlap


Overlap_t
Extend_Alignment(Match_Node_t  *Match,
                 char          *S,     int32   S_Len,
                 char          *T,     int32   T_Len,
                 int32         &S_Lo,  int32   &S_Hi,
                 int32         &T_Lo,  int32   &T_Hi,
                 int32         &Errors,
                 Work_Area_t   *WA);



double AS_OVL_ERROR_RATE = 0.06;

oicParameters  G;

double  Branch_Match_Value = 0.0;
double  Branch_Error_Value = 0.0;


//  Mirrored from overlapInCore.C
void
Initialize_Work_Area(Work_Area_t * WA, int id) {
  uint64  allocated = 0;

  WA->String_Olap_Size  = INIT_STRING_OLAP_SIZE;
  WA->String_Olap_Space = new String_Olap_t [WA->String_Olap_Size];

  WA->Match_Node_Size  = INIT_MATCH_NODE_SIZE;
  WA->Match_Node_Space = new Match_Node_t [WA->Match_Node_Size];

  allocated += WA->String_Olap_Size * sizeof (String_Olap_t);
  allocated += WA->Match_Node_Size  * sizeof (Match_Node_t);

  WA->status     = 0;
  WA->thread_id  = id;

  WA->overlapsLen = 0;
  WA->overlapsMax = 1024 * 1024 / sizeof(ovsOverlap);
  WA->overlaps    = new ovsOverlap [WA->overlapsMax];

  allocated += sizeof(ovsOverlap) * WA->overlapsMax;

  WA->editDist = new prefixEditDistance(G.Doing_Partial_Overlaps, G.maxErate);

  //fprintf(stderr, "Initialize_Work_Area:  MAX_ERRORS=%d ERROR_LIMIT=%d  allocated "F_U64"MB\n", MAX_ERRORS, ERROR_LIMIT, allocated >> 20);
}




#define MAX_LEN 10240

int
main(int argc, char **argv) {

  char  *aStr = new char [MAX_LEN];
  int32  aLen = 0;
  int32  aEnd = 0;

  char  *bStr = new char [MAX_LEN];
  int32  bLen = 0;
  int32  bEnd = 0;

  int32  toEnd = 0;

  Work_Area_t  *WA = new Work_Area_t;

  Initialize_Work_Area(WA, 0);

#if 0
  {
    WA->Edit_Match_Limit = new int32 [AS_MAX_READLEN];

    for  (int i = 0;  i <= ERRORS_FOR_FREE;  i ++)
      WA->Edit_Match_Limit[i] = 0;

    int Start = 1;

    for  (int e = ERRORS_FOR_FREE + 1;  e < MAX_ERRORS;  e ++) {
      Start = Binomial_Bound (e - ERRORS_FOR_FREE, AS_OVL_ERROR_RATE, Start, EDIT_DIST_PROB_BOUND);
      WA->Edit_Match_Limit[e] = Start - 1;
      assert (WA->Edit_Match_Limit[e] >= WA->Edit_Match_Limit[e - 1]);
    }
  }

  {
    WA->Error_Bound = new int32 [AS_MAX_READLEN];

    for  (int i = 0;  i <= AS_MAX_READLEN;  i ++)
      WA->Error_Bound[i] = (int) (i * AS_OVL_ERROR_RATE + 0.0000000000001);
  }
#endif

  FILE *F = fopen(argv[1], "r");


  fgets(aStr, MAX_LEN, F);  //  Header
  fgets(aStr, MAX_LEN, F);  //  Read

  aLen = strlen(aStr);
  aStr[aLen-1] = 0;
  aLen--;


  fgets(bStr, MAX_LEN, F);  //  Header
  fgets(bStr, MAX_LEN, F);  //  Read

  bLen = strlen(bStr);
  aStr[aLen-1] = 0;
  bLen--;


  while (!feof(F)) {
    fprintf(stderr, "Starting!\n");

    Match_Node_t *match  = new Match_Node_t;
    int32         errors = 0;

    int32         aLo = 0, aHi = 0;
    int32         bLo = 0, bHi = 0;

    match->Offset = 100;   //  Begin position in b
    match->Len    = 20;    //  Length of exact match between a and b
    match->Start  = 90;    //  Begin position in a
    match->Next   = 0;     //  Not used

    Overlap_t  ovl = Extend_Alignment(match,
                                      aStr, aLen,
                                      bStr, bLen,
                                      aLo,  aHi,
                                      bLo,  bHi,
                                      errors,
                                      WA);

    fprintf(stderr, "overlapType = %d  overlap  %d-%d  %d-%d  errors = %d  llen = %d rlen = %d\n",
            ovl, aLo, aHi, bLo, bHi, errors, WA->editDist->Left_Delta_Len, WA->editDist->Right_Delta_Len);

#if 0
    int res = Prefix_Edit_Dist(aStr, aLen,
                               bStr, bLen,
                               100,      //  Error Limit
                               &aEnd, &bEnd,
                               &toEnd,
                               WA);

    fprintf(stderr, "res = %d  aEnd =%d  bEnd = %d  toEnd = %d\n", res, aEnd, bEnd, toEnd);
#endif


    fgets(bStr, MAX_LEN, F);  //  Header
    fgets(bStr, MAX_LEN, F);  //  Read

    bLen = strlen(bStr);
    aStr[aLen-1] = 0;
    bLen--;
  }

  return(0);
}


