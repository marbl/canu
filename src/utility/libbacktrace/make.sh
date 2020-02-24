#!/bin/sh

OPTS="-W -Wall -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wold-style-definition -Wmissing-format-attribute -Wcast-qual -g -O2"

gcc48 -I. -funwind-tables $OPTS -c atomic.c     -o atomic.o
gcc48 -I. -funwind-tables $OPTS -c dwarf.c      -o dwarf.o
gcc48 -I. -funwind-tables $OPTS -c fileline.c   -o fileline.o
gcc48 -I. -funwind-tables $OPTS -c posix.c      -o posix.o
gcc48 -I. -funwind-tables $OPTS -c print.c      -o print.o
gcc48 -I. -funwind-tables $OPTS -c sort.c       -o sort.o
gcc48 -I. -funwind-tables $OPTS -c state.c      -o state.o
gcc48 -I. -funwind-tables $OPTS -c backtrace.c  -o backtrace.o
gcc48 -I. -funwind-tables $OPTS -c simple.c     -o simple.o
gcc48 -I. -funwind-tables $OPTS -c elf.c        -o elf.o
gcc48 -I. -funwind-tables $OPTS -c mmapio.c     -o mmapio.o
gcc48 -I. -funwind-tables $OPTS -c mmap.c       -o mmap.o

#/bin/sh ./libtool --tag=CC   --mode=link gcc48 -funwind-tables -frandom-seed=libbacktrace.la -W -Wall -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wold-style-definition -Wmissing-format-attribute -Wcast-qual  -g -O2   -o libbacktrace.la -rpath /work/canu-stack-trace/src/AS_UTL/libbacktrace/../../../FreeBSD-amd64/lib atomic.lo dwarf.lo fileline.lo posix.lo print.lo sort.lo state.lo backtrace.lo simple.lo elf.lo mmapio.lo mmap.lo

ar cru libbacktrace.a  atomic.o dwarf.o fileline.o posix.o print.o sort.o state.o backtrace.o simple.o elf.o mmapio.o mmap.o
ranlib libbacktrace.a

