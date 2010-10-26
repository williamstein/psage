#ifndef ELL_SURFACE_H
#define ELL_SURFACE_H

#include "helper.h"
#include "lzz_pEratX.h"

// meant for elliptic surfaces in char != 2,3

class ell_surfaceInfoT {
  public:
    class affine_model {
      public:
	// divisors for additive reduction
	zz_pEX	A, I_star, II, II_star, III, III_star, IV, IV_star;

	// divisors for multiplicative reduction
	zz_pEX	M_sp, M_ns;

      private:
	void init();

	// cardinality of scalar field
	long    q;

	void minimize();

      public:
	affine_model();

	void init(const zz_pEratX& a4, const zz_pEratX& a6);
	void init(const zz_pEX& a4, const zz_pEX& a6);

	// j-invariant
	zz_pEratX  j;

	// coefficients
	zz_pEX  a4, a6;

	// discriminant
	zz_pEX  disc;

	// contribution to sign of functional equation
	int	epsilon;

	// calculates Kodaira type about pi
	void kodaira(const zz_pEX& pi);
    };

  public:
    long    ref_count;	// for garbage collection
    long    q;

    ell_surfaceInfoT(const zz_pEratX& a4, const zz_pEratX& a6);

    ~ell_surfaceInfoT();

    affine_model  finite_model, infinite_model;
    int     sign, deg_L;

    zz_pEX finite_A();
    int    infinite_A();

    zz_pEX finite_M();
    int    infinite_M();
};

typedef ell_surfaceInfoT  *ell_surfaceInfoPtr;

extern ell_surfaceInfoPtr   ell_surfaceInfo;

class ell_surfaceContext {
  private:
    ell_surfaceInfoT   *ptr;

  public:
    void save();
    void restore() const;

    ell_surfaceContext() { ptr = NULL; }
    ell_surfaceContext(const zz_pEratX& a4, const zz_pEratX& a6);

    ell_surfaceContext(const ell_surfaceContext&);

    ell_surfaceContext& operator=(const ell_surfaceContext&);

    ~ell_surfaceContext();
};

class ell_surface {
  private:

  public:
    // set curve to E/F_q : y^2 = x^3 + a4*x + a6
    static void init(const zz_pEX& a4, const zz_pEX& a6);
    static void init(const zz_pEratX& a4, const zz_pEratX& a6);
    static void init(const zz_pEratX& a1, const zz_pEratX& a2,
		     const zz_pEratX& a3, const zz_pEratX& a4, 
		     const zz_pEratX& a5);

    static ell_surfaceInfoPtr getSurfaceInfo() { return ell_surfaceInfo; }
};

void get_an(ZZ_pEX& a4, ZZ_pEX& a6);
void get_j_invariant(ZZ_pEX& j_num, ZZ_pEX& j_den);

void get_reduction(ZZ_pEX **divisors, int finite);
void get_disc(ZZ_pEX& disc, int finite);

#endif // ELL_SURFACE_H

