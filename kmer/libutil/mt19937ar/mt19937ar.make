# -*- makefile -*-

LIBBRI/     :=$(call MakePath,$/../../libbri/)

#  Just like in md5lib, $/ is not this directory when building from
# BRI/.  This is a gross hack to fix it, but it also breaks building
# in this directory.
#
MTDIR = ${EXTERNAL/}mt19937ar


$/.C_EXES     := $/mt19937ar-test
$/.C_SRCS     := $/mt19937ar.c $/test.c
$/.C_LIBS     := $/libmt19937ar.a

$/.CLEAN      := $/*.o $/test.c
$/.REAL-CLEAN := $/*.o $/test.c $/mt19937ar-test

$/libmt19937ar.a: $/mt19937ar.o $/test.o

$/mt19937ar-test: $/mt19937ar.o $/mt19937ar-test.o

$/test.c: $/mt19937ar-test
	@-${MTDIR}/mt19937ar-test | diff - ${MTDIR}/mt19937ar.out > ${MTDIR}/diffs 2>&1
	@if test -s ${MTDIR}/diffs ; then echo 'MT19937: TEST FAILED'; else echo 'MT19937: Test Passed'; fi
	@touch ${MTDIR}/test.c
	@${MTDIR}/mt19937ar-test | diff - ${MTDIR}/mt19937ar.out > ${MTDIR}/diffs 2>&1

$(eval $/%.d $/%.o:  CFLAGS+= -I${LIBBRI/})



