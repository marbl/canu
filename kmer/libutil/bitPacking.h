#ifndef BRI_BITPACKING_H
#define BRI_BITPACKING_H

//  Routines used for stuffing bits into a word array.
//
//


//  Returns 'siz' bits from the stream based at 'ptr' and currently at
//  location 'pos'.  The position of the stream is not changed.
//
u64bit
getDecodedValue(u64bit *ptr,
                u64bit  pos,
                u64bit  siz);


//  Copies the lowest 'siz' bits in 'val' to the stream based at 'ptr'
//  and currently at 'pos'.  The position of the stream is not
//  changed.
//
void
setDecodedValue(u64bit *ptr,
                u64bit  pos,
                u64bit  siz,
                u64bit  val);


//  Like getDecodedValue() but will pre/post increment/decrement the
//  value stored in the stream before in addition to returning the
//  value.
//
//  preIncrementDecodedValue(ptr, pos, siz) === x = getDecodedValue(ptr, pos, siz) + 1;
//                                              setDecodedValue(ptr, pos, siz, x);
//
//  preDecrementDecodedValue(ptr, pos, siz) === x = getDecodedValue(ptr, pos, siz) - 1;
//                                              setDecodedValue(ptr, pos, siz, x);
//
//  postIncrementDecodedValue(ptr, pos, siz) === x = getDecodedValue(ptr, pos, siz);
//                                               setDecodedValue(ptr, pos, siz, x + 1);
//
//  postDecrementDecodedValue(ptr, pos, siz) === x = getDecodedValue(ptr, pos, siz);
//                                               setDecodedValue(ptr, pos, siz, x - 1);
//
u64bit
preIncrementDecodedValue(u64bit *ptr,
                         u64bit  pos,
                         u64bit  siz);

u64bit
preDecrementDecodedValue(u64bit *ptr,
                         u64bit  pos,
                         u64bit  siz);

u64bit
postIncrementDecodedValue(u64bit *ptr,
                          u64bit  pos,
                          u64bit  siz);

u64bit
postDecrementDecodedValue(u64bit *ptr,
                          u64bit  pos,
                          u64bit  siz);





//  NOTE:
//
//  I assume the addresses of bits in two consectuve words would be:
//
//    [ 3 2 1 0 ][ 3 2 1 0 ]
//
//  However, setDecodedValue() works as:
//
//    bitstream = [xxxx][xx00][0000]
//    value = 0cba (three bits have information. the fourth is empty)
//
//    setDecodedValue(bitstream, value, 3);
//
//    bitstream = [xxxx][xxcb][a000]
//
//  So, if you're reading the bitstream bit-by-bit, you get the values
//  high-order first.
//



inline
u64bit
getDecodedValue(u64bit *ptr,
                u64bit  pos,
                u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);
  }

  ret &= u64bitMASK(siz);

  return(ret);
}



inline
void
setDecodedValue(u64bit *ptr,
                u64bit  pos,
                u64bit  siz,
                u64bit  val) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1

  val &= u64bitMASK(siz);

  if (b1 >= siz) {
    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= val << (b1-siz);
  } else {
    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (val & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (val & (u64bitMASK(b2))) << (64-b2);
  }
}



inline
u64bit
preIncrementDecodedValue(u64bit *ptr,
                         u64bit  pos,
                         u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret  = 0;

  if (b1 >= siz) {
    ret  = ptr[wrd] >> (b1 - siz);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  return(ret);
}



inline
u64bit
preDecrementDecodedValue(u64bit *ptr,
                         u64bit  pos,
                         u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  return(ret);
}



inline
u64bit
postIncrementDecodedValue(u64bit *ptr,
                          u64bit  pos,
                          u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  ret--;
  ret &= u64bitMASK(siz);

  return(ret);
}





inline
u64bit
postDecrementDecodedValue(u64bit *ptr,
                          u64bit  pos,
                          u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  ret++;
  ret &= u64bitMASK(siz);

  return(ret);
}



#endif  //  BRI_BITPACKING_H
