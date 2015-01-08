
//
//  Strongly typed integer indices
//
//  operator<() is used to compare the ID against the length of the array holding the data

class abSeqID {
public:
  abSeqID()  { _id = UINT32_MAX; };
  ~abSeqID() {                   };

  //operator size_t () const { return(_id); };

  void   set(size_t id)   { _id = id; };
  uint32 get(void)        { return(_id); };

  bool   isInvalid(void)  { return(_id == UINT32_MAX); };
  bool   isValid(void)    { return(_id != UINT32_MAX); };

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

  void   set(size_t id)   { _id = id; };
  uint32 get(void)        { return(_id); };

  bool   isInvalid(void)  { return(_id == UINT32_MAX); };
  bool   isValid(void)    { return(_id != UINT32_MAX); };

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

  void   set(size_t id)   { _id = id; };
  uint32 get(void)        { return(_id); };

  bool   isInvalid(void)  { return(_id == UINT32_MAX); };
  bool   isValid(void)    { return(_id != UINT32_MAX); };

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

  void   set(size_t id)   { _id = id; };
  uint32 get(void)        { return(_id); };

  bool   isInvalid(void)  { return(_id == UINT32_MAX); };
  bool   isValid(void)    { return(_id != UINT32_MAX); };

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

  void   set(size_t id)   { _id = id; };
  uint32 get(void)        { return(_id); };

  bool   isInvalid(void)  { return(_id == UINT32_MAX); };
  bool   isValid(void)    { return(_id != UINT32_MAX); };

  bool   operator!=(abMultiAlignID that) const  { return(_id != that._id); };
  bool   operator==(abMultiAlignID that) const  { return(_id == that._id); };
  //bool   operator< (uint32  that) const  { return(_id < that); };
private:
  uint32   _id;
};

