/*
Copyright (C) 2011, Battelle National Biodefense Institute (BNBI);
all rights reserved. Authored by: Sergey Koren

This Software was prepared for the Department of Homeland Security
(DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
part of contract HSHQDC-07-C-00020 to manage and operate the National
Biodefense Analysis and Countermeasures Center (NBACC), a Federally
Funded Research and Development Center.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the Battelle National Biodefense Institute nor
  the names of its contributors may be used to endorse or promote
  products derived from this software without specific prior written
  permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using namespace std;

static const char *rcsid_AS_PBR_UTIL_C = "$Id: AS_PBR_util.cc,v 1.3 2012-08-22 14:41:00 skoren Exp $";

#include "AS_PBR_util.hh"

#include "AS_PER_encodeSequenceQuality.h"

uint32 loadSequence(gkStore *fs, map<AS_IID, uint8> &readsToPrint, map<AS_IID, char*> &frgToEnc) {
  gkFragment  fr;
  uint32 counter = 0;

  //fprintf(stderr, "Loading fragment seq\n");
  // figure out which libraries we want to use and store fragment info
  for (map<AS_IID, uint8>::const_iterator i = readsToPrint.begin(); i != readsToPrint.end(); i++) {
	 if (i->second != 0) {
		fs->gkStore_getFragment(i->first, &fr, GKFRAGMENT_QLT);
		int32 len = fr.gkFragment_getClearRegionLength();

		char *seq = fr.gkFragment_getSequence() + fr.gkFragment_getClearRegionBegin();
		char *qlt = fr.gkFragment_getQuality()  + fr.gkFragment_getClearRegionBegin();
		seq[len] = 0;
		qlt[len] = 0;
		char *enc = new char[len + 1];
		encodeSequenceQuality(enc, seq, qlt);
		frgToEnc[fr.gkFragment_getReadIID()] = enc;

		counter++;
	 }
  }
  return counter;
}

void convertOverlapToPosition(const OVSoverlap& olap, SeqInterval &pos, SeqInterval &bClr, uint32 alen, uint32 blen, bool forB) {
	if (olap.dat.ovl.type == AS_OVS_TYPE_OVL) {
		if (forB) {
			bClr.bgn = 0;
			bClr.end = alen;
			pos.bgn = -olap.dat.ovl.a_hang;
			pos.end = blen - olap.dat.ovl.b_hang;

			if (olap.dat.ovl.flipped) {
			      pos.bgn = blen - pos.bgn;
			      pos.end = blen - pos.end;
			}
		}
		else {
			bClr.bgn = 0;
			bClr.end = blen;
			pos.bgn = olap.dat.ovl.a_hang;
		    pos.end = alen + olap.dat.ovl.b_hang;

		   if (olap.dat.ovl.flipped) {
			  uint32 x = pos.end;
			  pos.end = pos.bgn;
			  pos.bgn = x;
		   }
		}
	} else if (olap.dat.ovl.type == AS_OVS_TYPE_OBT) {
		if (forB) {
			fprintf(stderr, "Error, OBT is not currently supported in this mode\n");
			exit(1);
		}
	   pos.bgn = olap.dat.obt.a_beg;
	   pos.end = olap.dat.obt.a_end;

	   if (!olap.dat.obt.fwd) {
		  pos.bgn = olap.dat.obt.a_end;
		  pos.end = olap.dat.obt.a_beg;
	   }
	   uint32 bend = olap.dat.obt.b_end_hi << 9 | olap.dat.obt.b_end_lo;
	   bClr.bgn = MIN(olap.dat.obt.b_beg, bend);
	   bClr.end = MAX(olap.dat.obt.b_beg, bend);
	}

	// update end points if necessary
	uint32 len = MAX(pos.bgn, pos.end) - MIN(pos.bgn, pos.end);
	if (forB) {
        if (len > alen) {
           if (pos.bgn > pos.end) {
              pos.bgn = pos.end + alen;
           } else {
              pos.end = pos.bgn + alen;
          }
        }
	}
	else {
	    if (len > blen) {
           if (pos.bgn > pos.end) {
              pos.bgn = pos.end + blen;
           } else {
              pos.end = pos.bgn + blen;
          }
	    }
	}
}
