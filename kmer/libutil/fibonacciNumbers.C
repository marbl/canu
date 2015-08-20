
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2004-APR-27 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-APR-11
 *      are Copyright 2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "util.h"

//
//  Argh, 64-bit guys use LU as their modifier, but 32-bit guys use LLU.
//

#ifdef TRUE64BIT
#define _(VAL) VAL ## LU
#else
#define _(VAL) VAL ## LLU
#endif

uint32
fibonacciValuesLen = 92;

uint64
fibonacciValues[92] = { _(1),
                        _(2),
                        _(3),
                        _(5),
                        _(8),
                        _(13),
                        _(21),
                        _(34),
                        _(55),
                        _(89),
                        _(144),
                        _(233),
                        _(377),
                        _(610),
                        _(987),
                        _(1597),
                        _(2584),
                        _(4181),
                        _(6765),
                        _(10946),
                        _(17711),
                        _(28657),
                        _(46368),
                        _(75025),
                        _(121393),
                        _(196418),
                        _(317811),
                        _(514229),
                        _(832040),
                        _(1346269),
                        _(2178309),
                        _(3524578),
                        _(5702887),
                        _(9227465),
                        _(14930352),
                        _(24157817),
                        _(39088169),
                        _(63245986),
                        _(102334155),
                        _(165580141),
                        _(267914296),
                        _(433494437),
                        _(701408733),
                        _(1134903170),
                        _(1836311903),
                        _(2971215073),
                        _(4807526976),
                        _(7778742049),
                        _(12586269025),
                        _(20365011074),
                        _(32951280099),
                        _(53316291173),
                        _(86267571272),
                        _(139583862445),
                        _(225851433717),
                        _(365435296162),
                        _(591286729879),
                        _(956722026041),
                        _(1548008755920),
                        _(2504730781961),
                        _(4052739537881),
                        _(6557470319842),
                        _(10610209857723),
                        _(17167680177565),
                        _(27777890035288),
                        _(44945570212853),
                        _(72723460248141),
                        _(117669030460994),
                        _(190392490709135),
                        _(308061521170129),
                        _(498454011879264),
                        _(806515533049393),
                        _(1304969544928657),
                        _(2111485077978050),
                        _(3416454622906707),
                        _(5527939700884757),
                        _(8944394323791464),
                        _(14472334024676221),
                        _(23416728348467685),
                        _(37889062373143906),
                        _(61305790721611591),
                        _(99194853094755497),
                        _(160500643816367088),
                        _(259695496911122585),
                        _(420196140727489673),
                        _(679891637638612258),
                        _(1100087778366101931),
                        _(1779979416004714189),
                        _(2880067194370816120),
                        _(4660046610375530309),
                        _(7540113804746346429),
                        _(12200160415121876738) };
