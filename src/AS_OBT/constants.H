#ifndef CONSTANTS_H
#define CONSTANTS_H

extern "C" {
#include "AS_MSG_pmesg.h" 
}

//
//  Constants used for overlap based trimming.
//


//  Mostly just to catch for when we don't include the correct assembler
//  header.
#ifdef AS_FRAG_MIN_LEN
#if AS_FRAG_MIN_LEN != 64
#error AS_FRAG_MIN_LEN changed!
#endif
#endif

#ifndef AS_FRAG_MIN_LEN
#warning AS_FRAG_MIN_LEN not defined!
#define AS_FRAG_MIN_LEN 64
#endif

//  sort-overlaps.C
#define OBT_FAR5PRIME        (29)
#define OBT_MIN_ERATE        (2.0 * 100)
#define OBT_MIN_DIFF         (75)

//  qualityTrim.C

//  merge-trimming.C
//
#define OBT_MODE_WIGGLE      (5)

#define OBT_CQ_LENGTH        (100)
#define OBT_CQO_LENGTH       (200)
#define OBT_CQO_OVERLAP      (100)
#define OBT_CQ_SHORT         (5)

#define OBT_MIN_LENGTH       (AS_FRAG_MIN_LEN)
#define OBT_QLT_CLOSE_5      (10)  // 5,6  use 5'mode, use 5'min>1
#define OBT_QLT_FAR_5        (50)  // 11   use min5'
#define OBT_QLT_MODE3        (150) // 9    use 3'mode
#define OBT_QLT_CLOSE_MAXM3  (30)  // 14   use max>1 close to max
#define OBT_QLT_CLOSE_MAX3   (100) // 12   use max3'


#endif  //  CONSTANTS_H
