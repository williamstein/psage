class zz_pEXRat {
 private:

 public:
  zz_pEXRat();
  zz_pEXRat(const zz_PEXRat&);
  void operator=(const zz_PEXRat&);

};

// elliptic surface over F_q

class ellsurf_pEInfoT {
 private:
  ellsurf_pEInfoT();
  ellsurf_pEInfoT(const ellsurf_pEInfoT&);
  void operator=(const ellsurf_pEInfoT&);

 public: 
  class jTable {
  public:
    long ref_count;

    jTable(long q);
    ~jTable();

    zz_pEX *table; // needs to be rational functions
  
};

class ellfiber_pE {
 private:

 public:

};
