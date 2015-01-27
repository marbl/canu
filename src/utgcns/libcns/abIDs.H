
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

