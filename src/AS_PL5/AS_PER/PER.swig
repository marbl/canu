/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Guillaume Marcais <guillaume@marcais.net>
 *************************************************************************/

%module "AS::PER"
%include "typemaps.i"

%{
#include "AS_MSG/AS_MSG_pmesg.h"
#include "AS_PER/AS_PER_gkpStore.h"
%}

/* Load basic int types
 */
#define __extension__
%import stdint.h
/* Define __cpluplus, otherwise AS_global.h will define the type bool
 * as an int, on which swig choaks.
 */
#define __cplusplus
%import AS_global.h
#undef __cplusplus

#define __IID AS_IID
#define __UID AS_UID

%apply unsigned int *OUTPUT { uint32 *start, uint32 *end };
%apply unsigned int *OUTPUT { __IID *iid, uint32 *type };

/* Lot of stuff copied from AS_PER_gkpStore.h.  Should we load it
 * instead, even if it defines way more stuff than what we want?
 */

// Stripped down version of the structs
typedef struct fragRecord {} fragRecord;
typedef struct GateKeeperStore {} GateKeeperStore;
typedef struct {
  char           comment[AS_PER_COMMENT_LEN];
  double         mean;
  double         stddev;
} GateKeeperLibraryRecord;

%extend fragRecord {
  ~fragRecord() { safe_free($self); }
  void ClearRegionUntrim(uint32 *start, uint32 *end) {
    getFragRecordClearRegion($self, start, end, AS_READ_CLEAR_UNTRIM);
  }
  void ClearRegion(uint32 *start, uint32 *end) {
    getFragRecordClearRegion($self, start, end, AS_READ_CLEAR_LATEST);
  }
};
%extend GateKeeperStore {
  GateKeeperStore(const char *path) {
    return openGateKeeperStore(path, 0);
  }
  ~GateKeeperStore() { closeGateKeeperStore($self); }

  //  int32 NumBatches() { return getNumGateKeeperBatches($self); };
  int32 NumLibraries() { return getNumGateKeeperLibraries($self); };
  int32 NumFragments() { return getNumGateKeeperFragments($self); };

  fragRecord *getFrag(__IID iid, int32 flags = 0) {
    fragRecord *fr = new_fragRecord();
    getFrag($self, iid, fr, flags);
    return fr;
  }
  GateKeeperLibraryRecord *getLibrary(int libiid) {
    return getGateKeeperLibrary($self, libiid);
  }

  __IID UIDtoIID(__UID uid) {
    uint32 type;
    return getGatekeeperUIDtoIID($self, uid, &type);
  }
  void UIDtoIID_type(__UID uid, __IID *iid, uint32 *type) {
    *iid = getGatekeeperUIDtoIID($self, uid, type);
  }
};

#define FRAG_SEQ 0x01
#define FRAG_QLT 0x02
#define FRAG_HPS 0x04
#define FRAG_SRC 0x08
#define FRAG_ALL 0x0f

#define IID_UNK     0
#define IID_BAT     1
#define IID_FRG     2
#define IID_LIB     3
%perlcode %{
  %IIDType = ( $IID_UNK => "UNK",
               $IID_BAT => "BAT",
               $IID_FRG => "FRG",
               $IID_LIB => "LIB",
               );
%}

__UID   getFragRecordUID(fragRecord *fr);
__IID   getFragRecordIID(fragRecord *fr);
__IID   getFragRecordMateIID(fragRecord *fr);
__IID   getFragRecordLibraryIID(fragRecord *fr);
int     getFragRecordIsDeleted(fragRecord *fr);
int     getFragRecordIsNonRandom(fragRecord *fr);
int     getFragRecordSequenceLength(fragRecord *fr);
int     getFragRecordQualityLength(fragRecord *fr);
int     getFragRecordHPSLength(fragRecord *fr);
int     getFragRecordSourceLength(fragRecord *fr);
char    *getFragRecordSequence(fragRecord *fr);
char    *getFragRecordQuality(fragRecord *fr);
char    *getFragRecordHPS(fragRecord *fr);
char    *getFragRecordSource(fragRecord *fr);

#define FRAGMETH(name)   *swig_ ## name ## _get = *::AS::PER::getFragRecord ## name
%perlcode {
  package AS::PER::GateKeeperStore;
  *open = *new;

  package AS::PER::fragRecord;
  FRAGMETH(UID);
  FRAGMETH(IID);
  FRAGMETH(MateIID);
  FRAGMETH(LibraryIID);
  FRAGMETH(IsDeleted);
  FRAGMETH(IsNonRandom);
  FRAGMETH(SequenceLength);
  FRAGMETH(QualityLength);
  FRAGMETH(HPSLength);
  FRAGMETH(SourceLength);
  FRAGMETH(Sequence);
  FRAGMETH(Quality);
  FRAGMETH(HPS);
  FRAGMETH(Source);
}

