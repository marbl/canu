PROG    = 
INCLUDE = -I../libbri
LIBS    = -L../libbri -lbri -lm

SRCS    = sim4polish-comparisons.c \
          sim4polish-copy.c \
          sim4polish-deleteexon.c \
          sim4polish-destroypolish.c \
          sim4polish-exons.c \
          sim4polish-polishtostring.c \
          sim4polish-printpolish.c \
          sim4polish-printpolishnormalized.c \
          sim4polish-read.c \
          sim4polish-stringtopolish.c \
          sim4polish.h \
          sim4polishBuilder.H \
          sim4polishBuilder.C
OBJS    = sim4polish-comparisons.o \
          sim4polish-copy.o \
          sim4polish-deleteexon.o \
          sim4polish-destroypolish.o \
          sim4polish-exons.o \
          sim4polish-polishtostring.o \
          sim4polish-printpolish.o \
          sim4polish-printpolishnormalized.o \
          sim4polish-read.o \
          sim4polish-stringtopolish.o \
          sim4polishBuilder.o

all: libsim4polish.a

install:

libsim4polish.a: $(OBJS) buildinfo-libsim4polish.o
	rm -f libsim4polish.a
	ar $(AROPTS) libsim4polish.a buildinfo-libsim4polish.o $(OBJS)

buildinfo-libsim4polish.c: $(SRCS)
	../util/buildinfo.pl libsim4polish "$(CCOPTS)" $(SRCS)

clean:
	rm -f $(PROG) *.o *.a buildinfo-*

include ../Makefile-configuration
