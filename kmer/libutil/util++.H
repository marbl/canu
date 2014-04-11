#ifndef UTIL_PLUS_PLUS_H
#define UTIL_PLUS_PLUS_H

#include "util.h"

//  These are all inlined, and C doesn't want to listen to that, so
//  they're here.
//
#include "bitOperations.h"
#include "bitPacking.h"
#include "endianess.H"

//  Various methods for encoding numbers into a bitstream.
//
//  Still missing:
//    minimal binary
//    golomb (actually rice, since power of two)
//    teuhola exponential golomb
//
//  And a nice way of getting parameters to those (and generalizedUnary)
//
#include "unaryEncoding.h"
#include "generalizedUnaryEncoding.h"
#include "eliasGammaEncoding.h"
#include "eliasDeltaEncoding.h"
#include "fibonacciEncoding.h"

//  Lists?
#include "uint32List.H"

//  Now the good stuff!
//
#include "speedCounter.H"
//#include "bzipBuffer.H"
#include "readBuffer.H"
#include "splitToWords.H"
#include "bitPackedArray.H"
#include "bitPackedFile.H"
#include "recordFile.H"
#include "intervalList.H"
#include "bigQueue.H"
#include "sweatShop.H"
#include "logMsg.H"

#endif  //  UTIL_PLUS_PLUS_H

