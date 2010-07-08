#include "sim4.H"

mss_t
Sim4::masks_shifts(char seed[32]) {

    struct positions *MP;
    int    seedLength = strlen(seed);
    int    masknum=0,matchedLength=0;
    int    i, k;
    int    total=0,maskSeedLength=0;
    char   seed_mask[2*seedLength+1];

    for(i=0;i<seedLength;i++){
      if(seed[i] == '0') {
        seed_mask[2*i]   = '0';
        seed_mask[2*i+1] = '0';
      }
      else if(seed[i] == 'x') {
        seed_mask[2*i]   = '0';
        seed_mask[2*i+1] = '1';
      }
      else if(seed[i] == '1') {
        seed_mask[2*i]   = '1';
        seed_mask[2*i+1] = '1';
      }
      else {
        printf("The seed can only contain 0, 1, or x, or any other characters\n");
        exit(1);
      }
    }
    seed_mask[2*seedLength] = '\0';
    maskSeedLength = strlen(seed_mask);

    for(i=0;i<maskSeedLength;i++) {
       if(seed_mask[i] == '1') matchedLength++;
    }
    if(seed_mask[0] == '1')  masknum = 1;

    for(i=1;i<maskSeedLength;i++) {
       if(seed_mask[i] == '1' && seed_mask[i-1]!='1') {
           masknum++;
           continue;
       }
    }

    mss_t MSS; 
    MSS.seedLength = seedLength;
    MSS.masknum    = masknum;
    MSS.matchedLength = matchedLength;
    MSS.mask   = (u64bitONE << (seedLength+seedLength-2)) - 1;
    MSS.masks  = (u64bit*)ckalloc(masknum*sizeof(u64bit));
    MSS.shifts = (int *)ckalloc(masknum*sizeof(int));

#ifdef debug
    printf(u64bitHEX, MSS.mask);
    printf("\n");
#endif

    MP = (struct positions *)ckalloc(masknum*sizeof(struct positions));
    k=0;
    if(seed_mask[0] == '1') MP[masknum-1].end = maskSeedLength - 1;
    for(i=0;i<maskSeedLength-1;i++){
      if(seed_mask[i]!= '1' && seed_mask[i+1] == '1') MP[masknum - k-1].end = maskSeedLength - (i+1) -1 ;
      if(seed_mask[i] == '1' && seed_mask[i+1] != '1') {
          MP[masknum-k-1].begin = maskSeedLength -i-1 ;
          k++;
      }
    }
    if(seed_mask[maskSeedLength-1] == '1') MP[0].begin = 0;
    if(seed_mask[maskSeedLength-1] == '1' && seed_mask[maskSeedLength-2]!= '1' ) MP[0].begin = 0;

    for(i=0;i<masknum;i++){
        MP[i].width = MP[i].end - MP[i].begin + 1;
        total = 0;
        for(k=0;k<i;k++){
          total = total + MP[k].width;
        }
        MP[i].result_shifts = MP[i].begin - total;
    }


    for(i=0;i<MSS.masknum;i++){
        MSS.masks[i]  = ( (u64bitONE << MP[i].width) - u64bitONE) << MP[i].begin;
        MSS.shifts[i] = MP[i].result_shifts;
    }
    ckfree(MP);

    MSS.type = ((2*MSS.seedLength == MSS.matchedLength) ? CONTINUOUS_SEED : SPACED_SEED);

    return MSS;
}

int
Sim4::mask_shift(u64bit ecode, mss_t MSS)
{
     int i=0;
     u64bit masked_ecode=0L;

     for(i=0;i<MSS.masknum; i++){
        masked_ecode += (ecode & MSS.masks[i]) >> MSS.shifts[i];
     }

     return (int)masked_ecode;
}

