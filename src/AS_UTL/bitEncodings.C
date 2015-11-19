
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
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "bitEncodings.H"

const
uint32
fibonacciValuesLen = 92;

const
uint64
fibonacciValues[92] = { 1LLU,
                        2LLU,
                        3LLU,
                        5LLU,
                        8LLU,
                        13LLU,
                        21LLU,
                        34LLU,
                        55LLU,
                        89LLU,
                        144LLU,
                        233LLU,
                        377LLU,
                        610LLU,
                        987LLU,
                        1597LLU,
                        2584LLU,
                        4181LLU,
                        6765LLU,
                        10946LLU,
                        17711LLU,
                        28657LLU,
                        46368LLU,
                        75025LLU,
                        121393LLU,
                        196418LLU,
                        317811LLU,
                        514229LLU,
                        832040LLU,
                        1346269LLU,
                        2178309LLU,
                        3524578LLU,
                        5702887LLU,
                        9227465LLU,
                        14930352LLU,
                        24157817LLU,
                        39088169LLU,
                        63245986LLU,
                        102334155LLU,
                        165580141LLU,
                        267914296LLU,
                        433494437LLU,
                        701408733LLU,
                        1134903170LLU,
                        1836311903LLU,
                        2971215073LLU,
                        4807526976LLU,
                        7778742049LLU,
                        12586269025LLU,
                        20365011074LLU,
                        32951280099LLU,
                        53316291173LLU,
                        86267571272LLU,
                        139583862445LLU,
                        225851433717LLU,
                        365435296162LLU,
                        591286729879LLU,
                        956722026041LLU,
                        1548008755920LLU,
                        2504730781961LLU,
                        4052739537881LLU,
                        6557470319842LLU,
                        10610209857723LLU,
                        17167680177565LLU,
                        27777890035288LLU,
                        44945570212853LLU,
                        72723460248141LLU,
                        117669030460994LLU,
                        190392490709135LLU,
                        308061521170129LLU,
                        498454011879264LLU,
                        806515533049393LLU,
                        1304969544928657LLU,
                        2111485077978050LLU,
                        3416454622906707LLU,
                        5527939700884757LLU,
                        8944394323791464LLU,
                        14472334024676221LLU,
                        23416728348467685LLU,
                        37889062373143906LLU,
                        61305790721611591LLU,
                        99194853094755497LLU,
                        160500643816367088LLU,
                        259695496911122585LLU,
                        420196140727489673LLU,
                        679891637638612258LLU,
                        1100087778366101931LLU,
                        1779979416004714189LLU,
                        2880067194370816120LLU,
                        4660046610375530309LLU,
                        7540113804746346429LLU,
                        12200160415121876738LLU };
