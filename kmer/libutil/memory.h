#ifndef MEMORY_H
#define MEMORY_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void    *memdup(const void *orig, size_t size);

#ifdef __cplusplus
}
#endif

#endif  //  MEMORY_H
