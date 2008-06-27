    ////////////////////////////////////////
    //
    //  1) Repeat ends - for repeat-ends, examine error rates and
    //  lengths.  (e.g., error rate vs length)
    //

    ////////////////////////////////////////
    //
    //  2) Polymorphism - for non-repeat ends, examine error rates and
    //  lengths.
    //

    ////////////////////////////////////////
    //
    //  3) Short insert check - count number of times a read overlaps
    //  with its mate.  Is it an innie overlap?  Is it a repeat end?
    //

    ////////////////////////////////////////
    //
    //  4) Genome length
    //

    ////////////////////////////////////////
    //
    //  5) Library randomness
    //



void
loadClearLengths(GateKeeperStore *gkp) {

  if (fragClearLength != NULL)
    return;

  FragStream      *frgStream = openFragStream(gkp, FRAG_S_INF);
  fragRecord       fr        = {0};
  uint64           maxIID    = getLastElemFragStore(gkp) + 1;

  fragLibrary     = (AS_IID *)safe_calloc(maxIID, sizeof(AS_IID));
  fragMateIID     = (AS_IID *)safe_calloc(maxIID, sizeof(AS_IID));
  fragClearLength = (uint16 *)safe_calloc(maxIID, sizeof(uint16));

  int  typ = AS_READ_CLEAR_OBT;

#if 0
  //  In general, we need to load a different clear length based on
  //  the type of overlap.  Since we currently only compute stats on
  //  OVL overlaps, this is disabled.

  switch (ovl->dat.ovl.type) {
    case AS_OVS_TYPE_OVL:
      typ = AS_READ_CLEAR_OBT;
      break;
    case AS_OVS_TYPE_OBT:
      typ = AS_READ_CLEAR_OBTINI;
      break;
    case AS_OVS_TYPE_MER:
      typ = AS_READ_CLEAR_UNTRIM;
      break;
    default:
      fprintf(stderr, "Unknown type %d in overlap.\n", ovl->dat.ovl.type);
      exit(1);
      break;
  }
#endif

  while (nextFragStream(frgStream, &fr)) {
    AS_IID  iid = getFragRecordIID(&fr);

    fragLibrary[iid]     = getFragRecordLibraryIID(&fr);
    fragMateIID[iid]     = getFragRecordMateIID(&fr);
    fragClearLength[iid] = (getFragRecordClearRegionEnd  (&fr, typ) -
                            getFragRecordClearRegionBegin(&fr, typ));
  }

  closeFragStream(frgStream);
}




//  a_hang   b_hang     type    label
//  (compared to zero)          (describes A)
//
//  any      any        0x00    all overlaps
//                      %0000
//
//  =        =          0x05    degenerate     ----------
//                      %0101                  ----------
//
//  =        <          0x06    5' containee   ----------
//                      %0110                  ---
//
//  =        >          0x07    5' contained   ---
//                      %0111                  ----------
//
//  <        =          0x09    3' contained          ---
//                      %1001                  ----------
//
//  <        <          0x0a    5' dovetal        ----------
//                      %1010                  ----------
//
//  <        >          0x0b    contained         ----
//                      %1011                  ----------
//
//  >        =          0x0d    3' containee   ----------
//                      %1101                         ---
//
//  >        <          0x0e    containee      ----------
//                      %1110                     ----
//
//  >        >          0x0f    3' dovetail    ----------
//                      %1111                     ----------
//
//  not an overlap -- 0x00
//  degenerate     -- 0x05
//  5' types       -- 0x06, 0x07, 0x0a
//  3' types       -- 0x09, 0x0d, 0x0f
//  C  types       -- 0x0b
//  C  types       -- 0x0e
//  unused types   -- 0x01, 0x02, 0x03, 0x04, 0x08, 0x0c
//
//
//
uint32
computeTypeOfOverlap(OVSoverlap ovl) {
  int32   ah = ovl.dat.ovl.a_hang;
  int32   bh = ovl.dat.ovl.b_hang;
  uint32  tp = 0;

  if (ah == 0)
    tp |= 0x00000004;
  else if (ah < 0)
    tp |= 0x00000008;
  else
    tp |= 0x0000000c;

  if (bh == 0)
    tp |= 0x00000001;
  else if (bh < 0)
    tp |= 0x00000002;
  else
    tp |= 0x00000003;

  return(tp);
}

uint32
overlapTypeIs3prime(uint32 tp) {
  return((tp == 0x06) || (tp == 0x07) || (tp == 0x0a));
}

uint32
overlapTypeIs5prime(uint32 tp) {
  return((tp == 0x09) || (tp == 0x0d) || (tp == 0x0f));
}


//  Swiped from AS_BOG/AS_BOG_BestOverlapGraph.cc::olapLength()
uint32
computeLengthOfOverlap(OVSoverlap ovl) {
  int32   ah = ovl.dat.ovl.a_hang;
  int32   bh = ovl.dat.ovl.b_hang;
  uint32  le = 0;

  if (ah < 0) {
    if (bh < 0)
      le = fragClearLength[ovl.a_iid] + bh;
    else
      le = fragClearLength[ovl.b_iid] + ah - bh;
  } else {
    if (bh < 0)
      le = fragClearLength[ovl.a_iid] + bh - ah;
    else
      le = fragClearLength[ovl.a_iid] - ah;
  }

  return(le);
}



RepeatModel *
computeRepeatModels(OverlapStore *ovs, GateKeeperStore *gkp) {
  int          i;
  uint32       ovl5   = 0;
  uint32       ovl3   = 0;
  AS_IID       lastID = 0;
  OVSoverlap   ovl;

  //  The [0] repeat model is a global model, all the others are
  //  specific to a single library.

  //  Allocate N models, one for each library.
  RepeatModel *rm = (RepeatModel *)safe_calloc(getNumGateKeeperLibraries(gkp) + 1, sizeof(RepeatModel));

  //  Then populate the histograms.
  for (i=0; i <= getNumGateKeeperLibraries(gkp); i++) {
    AS_UTL_histogramAllocate(&rm[i].hist5);
    AS_UTL_histogramAllocate(&rm[i].hist3);
  }

  AS_OVS_resetRangeOverlapStore(ovs);
  while (AS_OVS_readOverlapFromStore(ovs, &ovl, AS_OVS_TYPE_OVL) == TRUE) {

    if (ovl.a_iid != lastID) {

      //  Update the stats for this fragment.
      if (lastID > 0) {
        AS_IID  libIID = fragLibrary[lastID];

        AS_UTL_histogramAdd(&rm[0].hist5, ovl5);
        AS_UTL_histogramAdd(&rm[0].hist3, ovl3);

        AS_UTL_histogramAdd(&rm[libIID].hist5, ovl5);
        AS_UTL_histogramAdd(&rm[libIID].hist3, ovl3);
      }

      ovl5   = 0;
      ovl3   = 0;
      lastID = ovl.a_iid;
    }

    uint32 tp = computeTypeOfOverlap(ovl);

    if (overlapTypeIs5prime(tp))  ovl5++;
    if (overlapTypeIs3prime(tp))  ovl3++;
  }

  //  Update the stats for the last fragment.
  if (lastID > 0) {
    AS_IID  libIID = fragLibrary[lastID];

    AS_UTL_histogramAdd(&rm[0].hist5, ovl5);
    AS_UTL_histogramAdd(&rm[0].hist3, ovl3);

    AS_UTL_histogramAdd(&rm[libIID].hist5, ovl5);
    AS_UTL_histogramAdd(&rm[libIID].hist3, ovl3);
  }

  //  Examine the histograms, build a model.

  for (i=0; i <= getNumGateKeeperLibraries(gkp); i++) {
    char  label[256] = {0};
    char  name[FILENAME_MAX] = {0};
    FILE *file = NULL;

    AS_UTL_histogramCompute(&rm[i].hist5);
    AS_UTL_histogramCompute(&rm[i].hist3);

#warning bogus compute of repeatThreshold
    rm[i].repeatThreshold = ((rm[i].hist5.mode + 2 * rm[i].hist5.mad) +
                             (rm[i].hist3.mode + 2 * rm[i].hist3.mad)) / 2;

    sprintf(name, "%s.repeatmodel.lib.%03d.stats", outputPrefix, i);

    errno = 0;
    file = fopen(name, "w");
    if (errno) {
      fprintf(stderr, "Couldn't open '%s' for write: %s\n", name, strerror(errno));
      exit(1);
    }

    fprintf(file, "repeatThreshold = %d\n", rm[i].repeatThreshold);

    sprintf(label, "Lib "F_IID" 5'", i);
    AS_UTL_histogramShow(&rm[i].hist5, file, label);

    sprintf(label, "Lib "F_IID" 3'", i);
    AS_UTL_histogramShow(&rm[i].hist3, file, label);

    fclose(file);

    sprintf(name,  "%s.repeatmodel.lib.%03d.5prime.dat", outputPrefix, i);
    sprintf(label, "%s.repeatmodel.lib.%03d.5prime.dat", outputPrefix, i);
    AS_UTL_histogramDump(&rm[i].hist5, name, label);

    sprintf(name,  "%s.repeatmodel.lib.%03d.3prime.dat", outputPrefix, i);
    sprintf(label, "%s.repeatmodel.lib.%03d.3prime.dat", outputPrefix, i);
    AS_UTL_histogramDump(&rm[i].hist3, name, label);
  }

#if 1
  //  Dump the repeat thresholds
  fprintf(stderr, "== Repeat Model ==\n");
  fprintf(stderr, "\n");
  for (i=0; i <= getNumGateKeeperLibraries(gkp); i++)
    fprintf(stderr, "repeatThreshold[%2d] = "F_U64"\n", i, rm[i].repeatThreshold);
  AS_UTL_histogramShow(&rm[0].hist5, stderr, "Global 5'");
  AS_UTL_histogramShow(&rm[0].hist3, stderr, "Global 3'");
#endif

  return(rm);
}


//  Returns 0 if neither end is a repeat;
//          1 if both ends are repeats;
//          5 if the 5' end is;
//          3 if the 3' end is.
int
isRepeatEnd(OVSoverlap *ovls, uint64 ovlsLen, RepeatModel *rm) {
  uint64  i, n5=0, n3=0;

  //  If no chance of being a repeat, get out of here.
  if (ovlsLen < rm[0].repeatThreshold)
    return(0);

  for (i=0; i<ovlsLen; i++) {
    uint32  tp = computeTypeOfOverlap(ovls[i]);

    if (overlapTypeIs5prime(tp))
      n5++;

    if (overlapTypeIs3prime(tp))
      n3++;
  }

  if ((n5 >= rm[0].repeatThreshold) &&
      (n3 >= rm[0].repeatThreshold))
    return(1);

  if (n5 >= rm[0].repeatThreshold)
    return(5);

  if (n3 >= rm[0].repeatThreshold)
    return(3);

  return(0);
}

