#ifndef LIBBRITYPES_H
#define LIBBRITYPES_H

#if defined(__alpha) || defined(_AIX)
#define TRUE64BIT
#endif

//
//  Useful types.
//
//  *MASK() is defined for only unsigned types
//

#ifdef TRUE64BIT

//  Generic 64-bit platform

typedef unsigned long       u64bit;
typedef unsigned int        u32bit;
typedef signed long         s64bit;
typedef signed int          s32bit;
typedef unsigned short      u16bit;
typedef signed short        s16bit;
typedef unsigned char       u8bit;
typedef signed char         s8bit;

#define  u64bitZERO      (0x0000000000000000LU)
#define  u64bitONE       (0x0000000000000001LU)
#define  u64bitMAX       (0xffffffffffffffffLU)
#define  u64bitMASK(X)   ((~u64bitZERO) >> (64 - (X)))
#define  u64bitFMT       "%lu"
#define  u64bitHEX       "0x%016lx"
#define  s64bitFMT       "%ld"

#define  u32bitZERO      (0x00000000U)
#define  u32bitONE       (0x00000001U)
#define  u32bitMAX       (0xffffffffU)
#define  u32bitMASK(X)   ((~u32bitZERO) >> (32 - (X)))
#define  u32bitFMT       "%u"
#define  u32bitHEX       "0x%08x"
#define  s32bitFMT       "%d"

#define  u16bitZERO      (0x0000)
#define  u16bitONE       (0x0001)
#define  u16bitMAX       (0xffff)
#define  u16bitMASK(X)   ((~u16bitZERO) >> (16 - (X)))

#define  u8bitZERO       (0x00)
#define  u8bitONE        (0x01)
#define  u8bitMAX        (0xff)
#define  u8bitMASK(X)    ((~u8bitZERO) >> (8 - (X)))

#else

//  Generic 32-bit platform

typedef unsigned long long  u64bit;
typedef unsigned long       u32bit;
typedef signed long long    s64bit;
typedef signed long         s32bit;
typedef unsigned short      u16bit;
typedef signed short        s16bit;
typedef unsigned char       u8bit;
typedef signed char         s8bit;

#define  u64bitZERO      (0x0000000000000000LLU)
#define  u64bitONE       (0x0000000000000001LLU)
#define  u64bitMAX       (0xffffffffffffffffLLU)
#define  u64bitMASK(X)   ((~u64bitZERO) >> (64 - (X)))
#define  u64bitFMT       "%llu"
#define  u64bitHEX       "0x%016llx"
#define  s64bitFMT       "%lld"

#define  u32bitZERO      (0x00000000LU)
#define  u32bitONE       (0x00000001LU)
#define  u32bitMAX       (0xffffffffLU)
#define  u32bitMASK(X)   ((~u32bitZERO) >> (32 - (X)))
#define  u32bitFMT       "%lu"
#define  u32bitHEX       "0x%08lx"
#define  s32bitFMT       "%ld"

#define  u16bitZERO      (0x0000)
#define  u16bitONE       (0x0001)
#define  u16bitMAX       (0xffff)
#define  u16bitMASK(X)   ((~u16bitZERO) >> (16 - (X)))

#define  u8bitZERO       (0x00)
#define  u8bitONE        (0x01)
#define  u8bitMAX        (0xff)
#define  u8bitMASK(X)    ((~u8bitZERO) >> (8 - (X)))

#endif  //  TRUE64BIT

#endif  //  LIBBRITYPES_H
