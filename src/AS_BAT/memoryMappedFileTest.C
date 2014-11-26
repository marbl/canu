#include "memoryMappedFile.H"

int
main(int argc, char **argv) {

  memoryMappedFileRW<uint32>  *write = new memoryMappedFileRW<uint32>("testWrite");

  for (uint32 ii=0; ii<1000000000; ii++)
    (*write)[ii] = ii;

  delete write;

  exit(0);
}
