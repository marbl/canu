
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
 *    Sergey Koren beginning on 2015-DEC-28
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <cassert>
#include "Alignment.H"

using namespace dagcon;

Alignment::Alignment() :
    start(0),
    end(0),
    frgid(""),
    qstr(""),
    tstr("") { }


Alignment normalizeGaps(Alignment& aln) {
    size_t qlen = aln.qstr.length(), tlen = aln.tstr.length();
    assert(qlen == tlen);
    std::string qNorm, tNorm;

    // convert mismatches to indels
    for (size_t i=0; i < qlen; i++) {
        char qb = aln.qstr[i], tb = aln.tstr[i];
        if (qb != tb && qb != '-' && tb != '-') {
            qNorm += '-';
            qNorm += qb;
            tNorm += tb;
            tNorm += '-';
        } else {
            qNorm += qb;
            tNorm += tb;
        }
    }

    // update lengths
    qlen = qNorm.length();
    tlen = tNorm.length();

    // push gaps to the right, but not past the end
    for (size_t i=0; i < qlen-1; i++) {
        // pushing target gaps
        if (tNorm[i] == '-') {
            size_t j = i;
            while (true) {
                char c = tNorm[++j];
                if (c != '-' || j > qlen - 1) {
                    if (c == qNorm[i]) {
                        tNorm[i] = c;
                        tNorm[j] = '-';
                    }
                    break;
                }
            }
        }

        // pushing query gaps
        if (qNorm[i] == '-') {
            size_t j = i;
            while (true) {
                char c = qNorm[++j];
                if (c != '-' || j > tlen - 1) {
                    if (c == tNorm[i]) {
                        qNorm[i] = c;
                        qNorm[j] = '-';
                    }
                    break;
                }
            }
        }
    }

    // generate the final, normalized alignment strings
    Alignment finalNorm;
    finalNorm.frgid = aln.frgid;
    finalNorm.start = aln.start;
    for (size_t i=0; i < qlen; i++) {
        if (qNorm[i] != '-' || tNorm[i] != '-') {
            finalNorm.qstr += qNorm[i];
            finalNorm.tstr += tNorm[i];
        }
    }

    return finalNorm;
}
