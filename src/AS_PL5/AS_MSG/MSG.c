#include "MSG.h"
#include <string.h>
#include "swigperlrun.h"

/* Map message type to struct name of message.
 * What to do of the SPx (x = o, g, c, 1-9, a, d, ...)?
 * What about IDT? IBI?
 * Mapped to NULL for now.
 */
static const char *MessageTypeStruct[NUM_OF_REC_TYPES + 1] = {
  NULL,
  "BatchMesg",
  "VersionMesg",
  "DistanceMesg",
  "LibraryMesg",
  "FragMesg",
  "LinkMesg",

  "OverlapMesg",

  "UnitigOverlapMesg",

  "IntMateDistMesg",
  "IntAugFragMesg",
  "IntAugMatePairMesg",
  "IntUnitigMesg",
  "IntUnitigLinkMesg",
  "IntConConMesg",
  "IntContigLinkMesg",
  "IntScaffoldMesg",
  "InternalScaffoldLinkMesg",

  "SnapMateDistMesg",
  "AugFragMesg",
  "AugMatePairMesg",
  "SnapUnitigMesg",
  "SnapUnitigLinkMesg",
  "SnapConConMesg",
  "SnapContigLinkMesg",
  "SnapScaffoldMesg",
  "SnapScaffoldLinkMesg",
  "EndOfFileMesg"
};

static swig_type_info *MessageTypeSwigTypeInfo[NUM_OF_REC_TYPES + 1];
static int message_type_info_init = 0;

static void init_message_type_info()
{
  int i;
#define STR_LEN 100
  char str[STR_LEN+1];
  str[STR_LEN] = '\0';

 if(!message_type_info_init) {
    for(i = 0; i < NUM_OF_REC_TYPES; i++) {
      if(MessageTypeStruct[i]) {
        snprintf(str, STR_LEN, "_p_%s", MessageTypeStruct[i]);
        MessageTypeSwigTypeInfo[i] = SWIG_TypeQuery(str);
      } else {
        MessageTypeSwigTypeInfo[i] = NULL;
      }
    }
    message_type_info_init = 1;
  }
}

int GetMessageTypeFromStruct(char *s) {
  int i;

  for(i = 1; i < NUM_OF_REC_TYPES + 1; i++) {
    if(MessageTypeStruct[i])
      if(!strcmp(MessageTypeStruct[i], s))
        return i;
  }
  return -1;
}

void *GetSwigTypeInfo(int i) {
  if(i < 0 || i > NUM_OF_REC_TYPES)
    return NULL;

  return MessageTypeSwigTypeInfo[i];
}

void Parser_initialize(Parser *parser)
{
  init_message_type_info();

  parser->in = stdin;
  parser->out = stdout;
  parser->pmesg = NULL;
  parser->filter_mode = 0;
  parser->include[0] = 0;
  parser->exclude[0] = 0;
}

void Parser_next_message(Parser *parser, OutMesg *nmesg, void **ty, int *own)
{
  int ret, i, found, exclude;

  found = 0;
  exclude = 0;
  while(!found) {
    ret = ReadProtoMesg_AS(parser->in, &parser->pmesg);
    if(ret == EOF) {
      parser->pmesg = NULL;
      *nmesg = *ty = NULL;
      return;
    }
    if(!parser->include[0]) {
      found = 1;
    } else {
      for(i = 0; i < FILTER_SIZE && parser->include[i]; i++) {
        if(parser->pmesg->t == parser->include[i]) {
          found = 1;
          break;
        }
      }
    }
    exclude = 0;
    if(parser->exclude[0] && (found || parser->filter_mode)) {
      for(i = 0; i < FILTER_SIZE && parser->exclude[i]; i++) {
        if(parser->pmesg->t == parser->exclude[i]) {
          found = 0;
          exclude = 1;
          break;
        }
      }
    }
    if(parser->filter_mode && !found && !exclude) {
      WriteProtoMesg_AS(parser->out, parser->pmesg);
    }
  }

  parser->ty = MessageTypeSwigTypeInfo[parser->pmesg->t];
  if(!parser->ty) {
    *nmesg = *ty = NULL;
    SWIG_croak("Unsupported message type");
  }

  *nmesg = parser->pmesg->m;
  *ty = parser->ty;
  *own = 0;
  return;

 fail:
  SWIG_croak_null();
}

const char *Parser_message_type(Parser *parser)
{
  if(parser->pmesg) {
    return GetMessageName(parser->pmesg->t);
  } else {
    return "";
  }
}


void Parser_write_message(Parser *parser)
{
  if(parser->pmesg) {
    WriteProtoMesg_AS(parser->out, parser->pmesg);
  }
}
