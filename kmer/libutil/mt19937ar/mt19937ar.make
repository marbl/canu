# -*- makefile -*-

LIBBRI/     :=$(call MakePath,$/../../libbri/)

$/.C_SRCS     := $/mt19937ar.c
$/.C_LIBS     := $/libmt19937ar.a

$/.CLEAN      :=$/*.o

$(eval $/%.d $/%.o:  CFLAGS+= -I${LIBBRI/})

$/libmt19937ar.a: $/mt19937ar.o


