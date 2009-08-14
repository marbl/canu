////////////////////////////////////////
//
//  Placement constraint for reads.
//
//  When adding new fields, please update:
//    AS_PER/gkPlacement.H      -- this file
//
//
class gkPlacement {
public:
  gkPlacement() {
  };
  ~gkPlacement() {
  };

private:
public:
  AS_IID         frag;
  AS_IID         bound1;
  AS_IID         bound2;
};

