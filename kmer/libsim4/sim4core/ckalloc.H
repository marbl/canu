#ifndef CKALLOC_H
#define CKALLOC_H

extern  void  *ckalloc(size_t amount);
extern  char  *strsave(const char *s);
extern  void   ckfree(void *P);
extern  void   mem_init(void);
extern  void   mem_end(void);

#endif /* CKALLOC_H */
