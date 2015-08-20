
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
 *  This file is derived from:
 *
 *    src/AS_CNS/MultiAlignment_CNS.h
 *    src/AS_CNS/MultiAlignment_CNS_private.H
 *    src/AS_CNS/MultiAlignment_CNS_private.h
 *    src/utgcns/libcns/MultiAlignment_CNS_private.H
 *
 *  Modifications by:
 *
 *    Gennady Denisov from 2005-MAY-23 to 2007-OCT-25
 *      are Copyright 2005-2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUL-08 to 2013-AUG-01
 *      are Copyright 2005-2009,2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2006-FEB-13 to 2008-FEB-13
 *      are Copyright 2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-JAN-28 to 2009-SEP-25
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2011-SEP-21
 *      are Copyright 2011 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz beginning on 2014-NOV-17
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

//
//  Strongly typed integer indices
//
//  operator<() is used to compare the ID against the length of the array holding the data

class abSeqID {
public:
  abSeqID()           { _id = UINT32_MAX; };
  abSeqID(uint32 id)  { _id = id;         };
  ~abSeqID()          {                   };

  //operator size_t () const { return(_id); };
  abSeqID  &operator++(void) { ++_id;  return(*this); };  //  Prefix ++.

  void   set(size_t id)                { _id = id; };
  uint32 get(void)               const { return(_id); };

  bool   isInvalid(void)         const { return(_id == UINT32_MAX); };
  bool   isValid(void)           const { return(_id != UINT32_MAX); };

  bool   operator!=(abSeqID that) const  { return(_id != that._id); };
  bool   operator==(abSeqID that) const  { return(_id == that._id); };
  //bool   operator< (uint32  that) const  { return(_id < that); };
private:
  uint32   _id;
};


class abBaseID {
public:
  abBaseID()  { _id = UINT32_MAX; };
  ~abBaseID() {                   };

  //operator size_t () const { return(_id); };
  abBaseID  &operator+(uint32 x) { _id += x;  return(*this); };
  abBaseID  &operator-(uint32 x) { _id -= x;  return(*this); };

  void   set(size_t id)         { _id = id; };
  uint32 get(void)        const { return(_id); };

  bool   isInvalid(void)  const { return(_id == UINT32_MAX); };
  bool   isValid(void)    const { return(_id != UINT32_MAX); };

  bool   operator!=(abBaseID that) const  { return(_id != that._id); };
  bool   operator==(abBaseID that) const  { return(_id == that._id); };
  //bool   operator< (uint32  that) const  { return(_id < that); };
private:
  uint32   _id;
};


class abBeadID {
public:
  abBeadID()  { _id = UINT32_MAX; };
  ~abBeadID() {                   };

  //operator const size_t () { return(_id); };
  abBeadID  &operator+(uint32 x) { _id += x;  return(*this); };
  abBeadID  &operator-(uint32 x) { _id -= x;  return(*this); };

  void   set(size_t id)         { _id = id; };
  uint32 get(void)        const { return(_id); };

  bool   isInvalid(void)  const { return(_id == UINT32_MAX); };
  bool   isValid(void)    const { return(_id != UINT32_MAX); };

  bool   operator!=(abBeadID that) const  { return(_id != that._id); };
  bool   operator==(abBeadID that) const  { return(_id == that._id); };
  //bool   operator< (uint32  that) const  { return(_id < that); };
private:
  uint32   _id;
};


class abColID {
public:
  abColID()  { _id = UINT32_MAX; };
  ~abColID() {                   };

  //operator size_t () const { return(_id); };

  void   set(size_t id)         { _id = id; };
  uint32 get(void)        const { return(_id); };

  bool   isInvalid(void)  const { return(_id == UINT32_MAX); };
  bool   isValid(void)    const { return(_id != UINT32_MAX); };

  bool   operator!=(abColID that) const  { return(_id != that._id); };
  bool   operator==(abColID that) const  { return(_id == that._id); };
  //bool   operator< (uint32  that) const  { return(_id < that); };
private:
  uint32   _id;
};


class abMultiAlignID {
public:
  abMultiAlignID()  { _id = UINT32_MAX; };
  ~abMultiAlignID() {                   };

  //operator size_t () const { return(_id); };

  void   set(size_t id)         { _id = id; };
  uint32 get(void)        const { return(_id); };

  bool   isInvalid(void)  const { return(_id == UINT32_MAX); };
  bool   isValid(void)    const { return(_id != UINT32_MAX); };

  bool   operator!=(abMultiAlignID that) const  { return(_id != that._id); };
  bool   operator==(abMultiAlignID that) const  { return(_id == that._id); };
  //bool   operator< (uint32  that) const  { return(_id < that); };
private:
  uint32   _id;
};

