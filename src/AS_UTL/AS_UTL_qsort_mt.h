


void
qsort_mt(void *a,
         size_t n,
         size_t es,
         int (*cmp)(const void *, const void *),
         int maxthreads,
         int forkelem);

#define qsort(A, N, ES, CMP)  qsort_mt((A), (N), (ES), (CMP), 4, 64 * 1024)
