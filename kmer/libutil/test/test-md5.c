#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

//  Performs the md5 test suite using libbri.  MD5 itself is tested in
//  external/md5.
//
//  Appendix 5 of RFC 1321;
//
//  MD5 test suite:
//  MD5 ("") = d41d8cd98f00b204e9800998ecf8427e
//  MD5 ("a") = 0cc175b9c0f1b6a831c399e269772661
//  MD5 ("abc") = 900150983cd24fb0d6963f7d28e17f72
//  MD5 ("message digest") = f96b697d7cb7938d525a2f31aaf161d0
//  MD5 ("abcdefghijklmnopqrstuvwxyz") = c3fcd3d76192e4007dfb496cca67e13b
//  MD5 ("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789") = d174ab98d277d9f5a5611c2c9f419d9f
//  MD5 ("12345678901234567890123456789012345678901234567890123456789012345678901234567890") = 57edf4a22be3c955ac49da2e2107b67a
//

int
testit(char *str, char *ans) {
  md5_s  m;
  char   r[33];
  int    ret = 0;

  md5_toascii(md5_string(&m, str, strlen(str)), r);
  ret = strcmp(r, ans);
  if (ret)
    printf("ERROR: expect %s, got %s for %s\n", ans, r, str);
  return(ret == 0);
}

int
main(int argc, char **argv) {
  int ret = 7;
  
  ret -= testit("", "d41d8cd98f00b204e9800998ecf8427e");
  ret -= testit("a", "0cc175b9c0f1b6a831c399e269772661");
  ret -= testit("abc", "900150983cd24fb0d6963f7d28e17f72");
  ret -= testit("message digest", "f96b697d7cb7938d525a2f31aaf161d0");
  ret -= testit("abcdefghijklmnopqrstuvwxyz", "c3fcd3d76192e4007dfb496cca67e13b");
  ret -= testit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789", "d174ab98d277d9f5a5611c2c9f419d9f");
  ret -= testit("12345678901234567890123456789012345678901234567890123456789012345678901234567890", "57edf4a22be3c955ac49da2e2107b67a");
  exit(ret);
}
