#include "sim4.H"
#include "sim4b1_s.H"

mss_t::mss_t(char seed[32]) {
  position_t    MP[64];

  type          = 0;
  mask          = 0;
  masknum       = 0;
  seedLength    = strlen(seed);
  matchedLength = 0;

  int    total=0;
  int    maskSeedLength=0;
  char   seed_mask[2*seedLength+1];

  for (int i=0;i<seedLength;i++){
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

  for (int i=0;i<maskSeedLength;i++) {
    if(seed_mask[i] == '1') matchedLength++;
  }
  if(seed_mask[0] == '1')  masknum = 1;

  for (int i=1;i<maskSeedLength;i++) {
    if(seed_mask[i] == '1' && seed_mask[i-1]!='1') {
      masknum++;
      continue;
    }
  }

  assert(masknum <= 64);

  mask          = (u64bitONE << (seedLength+seedLength-2)) - 1;

#ifdef DEBUG
  printf(u64bitHEX, mask);
  printf("\n");
#endif

  int k=0;
  if(seed_mask[0] == '1') MP[masknum-1].end = maskSeedLength - 1;

  for (int i=0;i<maskSeedLength-1;i++){
    if(seed_mask[i]!= '1' && seed_mask[i+1] == '1') MP[masknum - k-1].end = maskSeedLength - (i+1) -1 ;
    if(seed_mask[i] == '1' && seed_mask[i+1] != '1') {
      MP[masknum-k-1].begin = maskSeedLength -i-1 ;
      k++;
    }
  }
  if(seed_mask[maskSeedLength-1] == '1') MP[0].begin = 0;
  if(seed_mask[maskSeedLength-1] == '1' && seed_mask[maskSeedLength-2]!= '1' ) MP[0].begin = 0;

  for (int i=0;i<masknum;i++){
    MP[i].width = MP[i].end - MP[i].begin + 1;
    total = 0;
    for(k=0;k<i;k++){
      total = total + MP[k].width;
    }
    MP[i].result_shifts = MP[i].begin - total;
  }


  for (int i=0;i<masknum;i++){
    masks[i]  = ( (u64bitONE << MP[i].width) - u64bitONE) << MP[i].begin;
    shifts[i] = MP[i].result_shifts;
  }

  type = ((2*seedLength == matchedLength) ? CONTINUOUS_SEED : SPACED_SEED);
}


u64bit
mss_t::mask_shift(u64bit ecode) {
  u64bit masked_ecode = 0;

  for (int i=0; i<masknum; i++)
    masked_ecode += (ecode & masks[i]) >> shifts[i];

  return(masked_ecode);
}

