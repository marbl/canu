#ifndef ATAC_COMMON_H
#define ATAC_COMMON_H


void   readHeader(char *inLine, FILE *in, char *file1, char *file2, FILE *out);
void   decodeMatch(splitToWords &W,
                   u32bit &iid1, u32bit &pos1, u32bit &len1, u32bit &ori1,
                   u32bit &iid2, u32bit &pos2, u32bit &len2, u32bit &ori2);

#endif  //  ATAC_COMMON_H
