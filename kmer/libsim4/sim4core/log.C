#include "sim4.H"

/* compute expected length of longest exact match .. */
/* K = (int) (log(.75*(double)len1)+log((double)len2))/log(4.0); */
/* .. and throw in a fudge factor */
/* K *= 1.4; */

/* The log4 arrays were computed to mimick the behaviour of the log formula
 * for computing the msp threshold in exon_cores(). For genomic_log4s,
 * entry i stores the value for the length of a genomic sequence
 * for which the contribution to the msp threshold is i/2, i.e.:
 *    1.4*log_4(3/4*len1) = i/2;  
 *
 * Similarly, cDNA_log4s entries store lengths of the cDNA sequence for which
 * the contribution to the msp threshold is i/2, i.e.:
 *    1.4*log_4(len2) = i/2;
 *
 * Both arrays are sorted in increasing order, and can be searched with 
 * binary search.
 */
#define GEN_LOG4_ENTRIES  45
#define CDNA_LOG4_ENTRIES 25

long const
genomic_log4s[]= {1, 2, 3, 5, 9, 15, 26, 42, 70, 114,
                  188, 309, 507, 832, 1365, 1365, 2240, 2240, 3675, 6029,
                  9892, 16231, 26629, 43690, 71681,
                  117606, 192953, 316573, 519392, 852152,
                  1398101, 2293823, 3763409, 6174516, 10130347,
                  16620564, 27268873, 44739242, 73402365, 120429110,
                  197584514, 324171126, 531858072, 872603963, 1431655765 };

long const
cDNA_log4s[]= {1, 1, 2, 4, 7, 11, 19, 32, 52, 86,
               141, 231, 380, 624, 1024, 1680, 2756, 4522, 7419, 12173,
               19972, 32768, 53761, 88204, 144715 };



#if 1
//  The original used to use a binary search; however, with so
//  few entries, it is twice as fast to just brute force search
//  the arrays.
//
//  This makes little changes in the output.
//
int
get_msp_threshold(int len1, int len2) {
  int   i, j;

  //  Find the index of the largest value smaller than our lengths.
  //
  i = 0;
  while (i<GEN_LOG4_ENTRIES) {
    if (genomic_log4s[i] > len1)
      break;
    i++;
  }
  i--;


  j = 0;
  while (j<CDNA_LOG4_ENTRIES) {
    if (cDNA_log4s[j] > len2)
      break;
    j++;
  }
  j--;

  //
  //  XXX:  This looks suspicious!
  //

  if ((i % 2) == 0)
    return(i/2+j/2);

  if ((j % 2) == 0)
    return(i/2+j/2);

  return(i/2+j/2+1);
}

#else

int
find_log_entry(long const *log4s, int n, int len, int offset)
{
  int a;

  a = n/2;  
  if ((len<log4s[a]) && (!a || (len>=log4s[a-1]))) 
    return max(0,(a-1))+offset;
  else if ((len>=log4s[a]) && ((a==n-1) || (len<log4s[a+1]))) 
    return min(n-1,(a+1))+offset;
  else if (len<log4s[a]) 
    return find_log_entry(log4s,a-1,len, offset);   
  else if (len>log4s[a])
    return find_log_entry(log4s+a+1,n-a-1,len, offset+a+1);
  return -1;
}

int
get_msp_threshold(int len1, int len2)
{
  int i, j;

  i = find_log_entry(genomic_log4s, GEN_LOG4_ENTRIES, len1, 0);
  j = find_log_entry(cDNA_log4s, CDNA_LOG4_ENTRIES, len2, 0);

  if (!(i % 2)) return (int)(i/2+j/2);
  else if (!(j % 2)) return (int)(i/2+j/2);
  else return (int)(i/2+j/2+1);
}
#endif








int
Sim4::getMSPthreshold(int defaultK, int l1, int l2) {
  int K;

  if (defaultK <= 0) {
    K = get_msp_threshold(l1, l2);

    /* compensate for the rounding in the log formula */
    if (K >= 0)
      K--;
  } else {
    K = defaultK;
  }

  return(K);
}
