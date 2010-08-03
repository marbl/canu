#ifndef SIM4B1_S_H
#define SIM4B1_S_H

#define CONTINUOUS_SEED   10
#define SPACED_SEED       11

typedef struct mss {
    int     type;
    u64bit  mask;
    int     masknum;
    int     seedLength;
    int     matchedLength;
    u64bit *masks;
    int    *shifts;
} mss_t;

struct mss masks_shifts(char *);
int        mask_shift(u64bit, struct mss);

#endif

