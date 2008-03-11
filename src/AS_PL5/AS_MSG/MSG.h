#ifndef __MSG_H__
#define __MSG_H__

#include "AS_MSG/AS_MSG_pmesg.h"

#define FILTER_SIZE     16
typedef struct {
  FILE          *in;
  FILE          *out;
  GenericMesg   *pmesg;
  void          *ty;
  int           filter_mode;
  MessageType   include[FILTER_SIZE];
  MessageType   exclude[FILTER_SIZE];
} Parser;
/* Generic type to return messages
 */
typedef void * OutMesg;
typedef void * InMesg;
int GetMessageTypeFromStruct(char *s);
void *GetSwigTypeInfo(int i);
#ifndef SWIG
void Parser_initialize(Parser *parser);
void Parser_next_message(Parser *parser, OutMesg *nmesg, void **ty, int *own);
void Parser_write_message(Parser *parser);
const char *Parser_message_type(Parser *parser);
#endif

#endif /* __MSG_H__ */
