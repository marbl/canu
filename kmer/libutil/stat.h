#ifndef STAT_H
#define STAT_H

#include <stdio.h>
#include <sys/stat.h>

//
//  Simple wrapper around stat(), lstat() and fstat().
//

#ifdef __cplusplus
extern "C" {
#endif

typedef struct stat stat_s;

stat_s   *stat_onPath(const char *path, stat_s *sb);
stat_s   *stat_onLink(const char *path, stat_s *sb);
stat_s   *stat_onDescriptor(int file, stat_s *sb);
stat_s   *stat_onFile(FILE *F, stat_s *sb);

void      stat_free(stat_s *sb);

int       stat_fileIsPipe(stat_s *sb);
int       stat_fileIsCharacterSpecial(stat_s *sb);
int       stat_fileIsDirectory(stat_s *sb);
int       stat_fileIsBlockSpecial(stat_s *sb);
int       stat_fileIsRegular(stat_s *sb);
int       stat_fileIsSymbolicLink(stat_s *sb);
int       stat_fileIsSocket(stat_s *sb);
int       stat_fileIsWhiteout(stat_s *sb);

uid_t     stat_getUID(stat_s *sb);
gid_t     stat_getGID(stat_s *sb);

double    stat_getAccessTime(stat_s *sb);
double    stat_getModificationTime(stat_s *sb);
double    stat_getStatusTime(stat_s *sb);

off_t     stat_getSize(stat_s *sb);


//  Convenience functions
//
int   fileExists(const char *path);
off_t sizeOfFile(const char *path);

#if 0
int   fileIsDirectory(const char *path);
int   fileIsRegular(const char *path);
#endif

#ifdef __cplusplus
}
#endif

#endif  //  STAT_H
