//
// (c) copyright 2006 Salman Butt and Chris Hall
//
// C++ elliptic curve class
//
// this code offers basic routines for elliptic curves over a finite field.
// currently it is limited to curves with a model of the form
// y^2=x^3+a4*x+a6, which excludes characteristic two and some curves in
// characteristic three.  because the latter two characteristics present
// other difficulties for the L-function code, we felt it was ok to make
// this exclusion.
//
// our class is modeled after NTL classes such as the finite field class
// ZZ_pE.  in particular, there is a global elliptic curve which is assumed
// to be defined over NTL's global finite field, and one must use context
// loading and unloading in order to switch between curves and finite
// fields.  on the one hand, the SAGE people assure us that this has proven
// to be a source of many bugs because one must be very careful about
// switching to the appropriate context.  on the other hand, we do not want
// to rewrite our code using a library other than NTL, so we will have to
// deal with this issue regardless of whether or not our own classes use
// the model.
//
// this code is *not* meant to provide full elliptic curve functionality,
// rather merely what we need to perform point counting.  this means that
// we have code for computing P+Q and [m]P for points P,Q and an integer m,
// but these routines are only mildly optimized.  our experience is that
// for a fixed finite field, one really only needs to use this code to
// compute the Euler factors for a (semi-)universal family over the j-line,
// and then to regard any other family as a combined pullback and quadratic
// twist.  the point is that the Euler factors for the latter kind of curve
// can be computed without any elliptic curve arithmetic given the Euler
// factors for the former, and in particular, one only needs the elliptic
// curve code when getting started over a new finite field.
//

#include <assert.h>
#include <string.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

#include "ell.h"
#include "helper.h"
#include "lzz_pEExtra.h"

#define INV_TABLE  1	// precompute a table of inverses in F_q^*?
#define X_ORDER    1	// order points by x-coordinates first? (then y?)

///////////////////////////////////////////////////////////////////////////
// elliptic curve parameter class
//   - we keep track of the coeffs a4,a6, the order of the finite field,
//   and a table of inverses for F_q^*
///////////////////////////////////////////////////////////////////////////

// constructor

ell_pEInfoT::ell_pEInfoT(const zz_pE& new_a4, const zz_pE& new_a6)
{
    ref_count = 1;

    a4 = new_a4;
    a6 = new_a6;
    q = to_long(zz_pE::cardinality());
}

// destructor

ell_pEInfoT::~ell_pEInfoT()
{
}

// copied from NTL

static void CopyPointer(ell_pEInfoPtr& dst, ell_pEInfoPtr src)
{
    if (src == dst) return;

    if (dst) {
        dst->ref_count--;

        if (dst->ref_count < 0) 
            Error("internal error: negative ell_pEContext ref_count");

        if (dst->ref_count == 0) delete dst;
   }

    if (src) {
        if (src->ref_count == NTL_MAX_LONG) 
            Error("internal error: ell_pEContext ref_count overflow");

        src->ref_count++;
    }

    dst = src;
}

////////////////////////////////////////////////////////////
// context switching class, shamelessly copied from NTL
////////////////////////////////////////////////////////////

ell_pEInfoPtr   ell_pEInfo = NULL;	// info for current elliptic curve

ell_pEContext::ell_pEContext(const zz_pE& a4, const zz_pE& a6)
{
    ptr = new ell_pEInfoT(a4, a6);
}

ell_pEContext::ell_pEContext(const ell_pEContext& a)
{
    ptr = NULL;
    CopyPointer(ptr, a.ptr);
}

ell_pEContext& ell_pEContext::operator=(const ell_pEContext& a)
{
    CopyPointer(ptr, a.ptr);
    return *this;
}

ell_pEContext::~ell_pEContext()
{
    CopyPointer(ptr, NULL);
}

void ell_pEContext::save()
{
    CopyPointer(ptr, ell_pEInfo);
}

void ell_pEContext::restore() const
{
    CopyPointer(ell_pEInfo, ptr);
}

//////////////////////////////

void ell_pE::init(zz_pE& a4, zz_pE& a6)
{
    ell_pEContext c(a4, a6);
    c.restore();
}

//////////////////////////////

ellpoint_pE::ellpoint_pE()
{
    x = 0;
    y = 0;
    identity_f = 1;
}

ellpoint_pE::ellpoint_pE(const zz_pE& new_x, const zz_pE& new_y)
{
    x = new_x;
    y = new_y;
    identity_f = 0;
}

// P + Q = R

void add(ellpoint_pE& R, const ellpoint_pE& P, const ellpoint_pE& Q)
{
    // temporary variables for calculation
    static zz_pE  *Rx = NULL, *Ry = NULL, *Rz = NULL, *Ru = NULL;
    static zz_pE  *lambda = NULL;

    static long   q = -1;

    // reuse point from previous invocation if possible
    if (ell_pEInfo->q != q) {
	q = ell_pEInfo->q;
	if (lambda != NULL) {
	    delete lambda;
	    delete Rx;
	    delete Ry;
	    delete Rz;
	    delete Ru;
	}
	// NB: the only way the following are deallocated is if they
	//     are replaced during a later invocation (see previous
	//     lines), so technically this is a tiny memory leak
	lambda = new zz_pE;
	Rx = new zz_pE;
	Ry = new zz_pE;
	Rz = new zz_pE;
	Ru = new zz_pE;
    }

    // check if either is identity
    if (P.identity_f) { R = Q; return; }
    if (Q.identity_f) { R = P; return; }
    
    // calculate slope between lines
    //  - special cases arise when P=+/-Q
    if (P.x == Q.x) {
	if (P.y != Q.y || P.y == 0) {
	    R.set();
	    return;
	}

	// point-doubling
	if (zz_pEExtraInfo->inv_table != NULL) {
	    // lambda  = (3*x(P)^2 + a4)/(2*y(P))
	    sqr(*Rx, P.x);
		add(*Ry, *Rx, *Rx);
		add(*Rz, *Rx, *Ry);
		add(*Rx, *Rz, ell_pEInfo->a4);
	    add(*lambda, P.y, P.y);
		*Ru  = zz_pEExtraInfo->inv_table->table[to_ulong(*lambda)];
	    mul(*lambda, *Rx, *Ru);
	} else {
	    *lambda = (3*sqr(P.x) + ell_pEInfo->a4) / (2*P.y);
	}
    } else {
	if (zz_pEExtraInfo->inv_table != NULL) {
	    // lambda  = (y(P)-y(Q))/(x(P)-x(Q))
	    sub(*Rz, P.x, Q.x);
		*Ru  = zz_pEExtraInfo->inv_table->table[to_ulong(*Rz)];
	    sub(*Rz, P.y, Q.y);
	    mul(*lambda, *Ru, *Rz);
	} else {
	    *lambda = (P.y - Q.y) / (P.x - Q.x);
	}
    }

    // lambda = slope of line between P,Q(,and -R)

    // x(R) = lambda^2 - x(P) - x(Q)
    sqr(*Ru, *lambda);
	sub(*Rz, *Ru, P.x);
	sub(*Rx, *Rz, Q.x);

    // y(R) = lambda*(x(P) - x(R)) - y(P) = -(lambda*(x(R)-x(P))+y(P))
    sub(*Ru, P.x, *Rx);
	mul(*Rz, *Ru, *lambda);
	sub(*Ry, *Rz, P.y);

    // set coordinates of R according to above
    R.set(*Rx, *Ry);
}

void neg(ellpoint_pE& Q, const ellpoint_pE& P)
{
    Q.set(P.x, -P.y);
}

// P - Q = R

void sub(ellpoint_pE& R, const ellpoint_pE& P, const ellpoint_pE& Q)
{
    static ellpoint_pE   T;
    neg(T, Q);
    add(R, P, T);
}

// addition/subtraction
ellpoint_pE operator+(const ellpoint_pE& P, const ellpoint_pE& Q)
   { ellpoint_pE x; add(x, P, Q); return x; }

ellpoint_pE operator-(const ellpoint_pE& P, const ellpoint_pE& Q)
   { ellpoint_pE x; sub(x, P, Q); return x; }

// negation
ellpoint_pE operator-(const ellpoint_pE& P)
   { ellpoint_pE x; x = P; if (x.identity_f == 0) x.y = -x.y; return x; }

ellpoint_pE& operator+=(ellpoint_pE& P, const ellpoint_pE& Q)
   { add(P, P, Q); return P; }

ellpoint_pE& operator-=(ellpoint_pE& P, const ellpoint_pE& Q)
   { sub(P, P, Q); return P; }

ellpoint_pE operator*(int m, const ellpoint_pE& P)
{
    static ellpoint_pE  *result = NULL, *Q = NULL;
    static long    q = -1;
    long    i;

    // reallocate as infrequently as possible
    if (ell_pEInfo->q != q) {
	q = ell_pEInfo->q;
	if (result != NULL) {
	    delete result;
	    delete Q;
	}
	// NB: the only way the following are deallocated is if they
	//     are replaced during a later invocation (see previous
	//     lines), so technically this is a tiny memory leak
	result = new ellpoint_pE;
	Q      = new ellpoint_pE;
    } else
	result->set();

    if (m == 0) return *result;
    if (m < 0)  return (-m)*P;

    for (i = 0, *Q = P; (1<<i) <= m; i++) {
	if (m & (1<<i))
	    *result += *Q;
	*Q += *Q;
    }

    return *result;
}

int is_zero(const ellpoint_pE& P) { return P.identity_f; }

// *slow* but safe routine for computing order of a point

long ellpoint_pE::order()
{
    ellpoint_pE  Q;
    long    o;

    // compute multiples of self until find 0
    Q = *this;
    for (o = 1; !is_zero(Q); o++)
        Q += *this;

    return o;
}

//////////////////////////////

unsigned long* to_ulong(ellpoint_pE& P)
{
    static unsigned long  result[2];
    if (is_zero(P)) return NULL;  // identity
    result[0] = to_ulong(P.x);
    result[1] = to_ulong(P.y);
    return result;
}

//////////////////////////////

long operator==(const ellpoint_pE& P, const ellpoint_pE& Q)
{
    if (is_zero(P) && is_zero(Q))  return 1;
    if (is_zero(P) || is_zero(Q))  return 0;
    return (P.x == Q.x) && (P.y == Q.y);
}

long operator<( ellpoint_pE& P, ellpoint_pE& Q)
{
    // O < anything else
    if (is_zero(P))   return 1;
    if (is_zero(Q))   return 0;

#if X_ORDER
    if (P.x < Q.x)    return 1;
    if (P.x == Q.x && P.y < Q.y)  return 1;
#else
    if (P.y < Q.y)    return 1;
    if (P.y == Q.y && P.x < Q.x)  return 1;
#endif
    return 0;
}

long operator<=(ellpoint_pE& P, ellpoint_pE& Q)
    { return (P < Q) || (P == Q); }

long operator>( ellpoint_pE& P, ellpoint_pE& Q)
    { return (Q <= P); }

long operator>=(ellpoint_pE& P, ellpoint_pE& Q)
    { return (Q < P); }

ostream& operator<<(ostream& s, const ellpoint_pE& P)
{
    if (P.identity_f)
	s << "0";
    else
	s << "(" << P.x << ", " << P.y << ")";
    return s;
}

//////////////////////////////
//
// routines for calculating group order

///////////////////////////////////////////////////////////////////////////

// binary tree code
//  - used for storing points in bsgs algorithm

struct node {
    unsigned long point[2];
    long	i;		// logarithm with respect to some point
    struct node *left;
    struct node *right;
};

typedef struct node node_t;
typedef struct node *node_p;

// compare two nodes

int cmp(unsigned long *a, unsigned long *b) {
    if (a[0] == b[0] && a[1] == b[1]) return 0;
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return  1;
    if (a[1] < b[1]) return -1;
    return 1;
}

// look to see if node lies in tree

static long lookup(struct node* node, unsigned long *target) {
    if (node == NULL)
        return(-1);
    else {
        int	r = cmp(node->point, target);
        if (r == 0)
	    return(node->i);
        else {
            if (r > 0)
		return(lookup(node->left, target));
            else
		return(lookup(node->right, target));
        }
    }
}

// add node to tree

void insert(node_p head, node_p node)
{
    node->left = node->right = NULL;

    if (node == head)
	return;

    int	r = cmp(node->point, head->point);

    if (r<0) {
        if (head->left == NULL)
            head->left = node;
        else
            insert(head->left, node);
    } else {
        if (head->right == NULL)
            head->right = node;
        else
            insert(head->right, node);
    }
}

///////////////////////////////////////////////////////////////////////////
//
// baby-step giant-step code

// - used to find an annihilator of a point on an elliptic curve
//
// input:
//   P : point on elliptic curve
//   M : number of baby steps
//   N : anchor point for giant step
//
// sketch:
//
// (1) compute [i]P for i=1..M-1 and return i if find [i]P=O
//
// (2) compute [N-M^2+j*M]P for j=0..M-1 and return (N-M^2+j*M-i) if
//     [i]P = [N-M^2+j*M]P
//
// (3) return -1
//
// - if terminates in (1), will return smallest possible i
// - if terminates in (2), will return smallest possible j
// - a return value of -1 implies the search failed
// - if called with N=M^2+1 and returns r>0, then r=order(P)

static node_t	*bsgs_annihilate_nodes = NULL;

long bsgs_annihilate(const ellpoint_pE& P, const long M, const long N)
{
    static ellpoint_pE  *R = NULL;
    static long		q = -1, max_M = -1;
    node_t	*nodes = bsgs_annihilate_nodes;
    long	i;

    // create workspace for storing multiples of P, if necessary
    if (nodes == NULL || ell_pEInfo->q != q || max_M < M) {
	// either need more workspace or have switched finite fields

	// clean up previous workspace, if necessary
	if (nodes != NULL) {
	    delete nodes;
	    delete R;
	}

	// NB: one can clean up the following allocation for nodes using
	//     ell_cleanup(), but the last copy of R used will always
	//     stick around

	// allocate new workspace
        nodes = bsgs_annihilate_nodes = new(node_t[M]);
	R = new ellpoint_pE;
	assert(nodes != NULL && R != NULL);

	// keep track of values current workspace depends on
        q = ell_pEInfo->q;
        max_M = M;
    }

    // baby steps: compute [i]P for i=1..M-1
    //  - obviously [0]P=O, hence the reason we don't consider i=0
    static unsigned long *point;
    *R = P;
    for (long i=1; i < M; i++, *R += P) {
	// claim: R = [i]P

	// convert R to convenient format for storing in a binary tree
        point = to_ulong(*R);
	if (point == NULL)
	    return i;	// [i]*P = O

	// add R to our binary tree
	//  - not on tree already because otherwise R=[i]P=[j]P and
	//    hence [i-j]P=O, which we would have seen
        memcpy(nodes[i-1].point, point, sizeof(unsigned long)*2);
        nodes[i-1].i = i;
        insert(nodes, nodes+(i-1));

	// NB: point points to statically allocated memory hence
	//     doesn't need to be freed
    }

    // claim: R=[M]P 

    // compute base point for giant steps
    ellpoint_pE Q = (N-((M-1)*M))*P;
    if (Q == *R)
	return N-(M*M);	// claim: order(P) | (N-M^2+M)-M = N^2-M

    // giant steps: compute Q + [M*i]P = for i=1..M
    for (long j = M-1; j >= 0; j--, Q += *R) {
	// claim: Q = [N-j*M]P
        if (is_zero(Q))
	    return N-j*M;

	// check to see if point was computed as a baby step
        if ((i=lookup(nodes, to_ulong(Q))) != -1)
            return N-j*M-i;	// [N-j*M]P = [i]P

	// don't do the last addition since we don't need it
        if (j == 0)
	    break;
    }

    // failed
    return -1;
}

// baby-step giant-step search for m in [0,M^2-1] so [m]P = Q
//
// input:
//   P : base point on elliptic curve
//   Q : point on to compute log_P for
//   M : number of baby steps
//
// sketch:
//
// (1) compute Q-[i]P for i=0..M-1 and return i if find [i]P=Q
//
// (2) compute [j*M]P for j=1..M-1 and return (i+j*M) if Q-[i]P=[j*M]P
//
// (3) return -1
//
// - if terminates in (1), will return smallest possible i
// - if terminates in (2), will return smallest possible i+j*M
// - a return value of -1 implies the search failed

static node_t   *bsgs_log_nodes = NULL;

long bsgs_log(const ellpoint_pE& P, const ellpoint_pE& Q, const long M)
{
    static ellpoint_pE	*R = NULL, *S = NULL;
    static long		q = -1, max_M = -1;
    node_t	*nodes = bsgs_log_nodes;
    long	i;

    // create workspace for storing multiples of P, if necessary
    if(nodes == NULL || ell_pEInfo->q != q || max_M < M) {
	// either need more workspace or have switched finite fields

	// clean up previous workspace, if necessary
	if (nodes != NULL)
	    delete nodes;
	if (R != NULL)
	    delete R;
	if (S != NULL)
	    delete S;

	// NB: one can clean up the following allocation for nodes using
	//     ell_cleanup(), but the last copies of R,S used will always
	//     stick around

	// allocate new workspace
        bsgs_log_nodes = nodes = new(node_t[M]);
	R = new ellpoint_pE;
	S = new ellpoint_pE;
	assert(nodes != NULL && R != NULL && S != NULL);

	// keep track of values current workspace depends on
        q = ell_pEInfo->q;
        max_M = M;
    }

    // baby steps: compute Q-[i]P for i=0,...,M-1
    static unsigned long *point;
    *R = Q;
    for(long i = 0; i < M; i++, *R -= P) {
	// claim: R = Q - [i]P

	// convert R to convenient format for storing in a binary tree
        point = to_ulong(*R);
        if (point == NULL)
	    return i;	// Q = [i]P

	// add R to our binary tree
	//  - if on tree already, then R=Q-[i]P=Q-[j]P and [i-j]P=O, which
	//    implies we won't find k so [k]P=Q because necessarily k<i-j
	//    and we already noticed [k]P!=Q for all k<i
	//  - rather that noticing we're in previous case, we just keep
	//    running since we expect to call this routine only when the
	//    order(P)>M, hence it can't occur
	// ==> we can assume [i]P!=[j]P for j<i when necessary
        memcpy(nodes[i].point, point, sizeof(unsigned long)*2);
        nodes[i].i = i;
        insert(nodes, nodes+i);

	// NB: point points to statically allocated memory hence
	//     doesn't need to be freed
    }

    // giant steps: compute S = [j*M]P for j=1..M-1
    *R = M*P;
    *S = *R;
    for (long j = 1; j < M; j++, *S += *R) {
	// S = [j*M]P

	// if [j*M]P==0, then Q does not lie in <P> because we already
	// checked whether [k]P==Q for 0<k<j*M
	if (is_zero(*S))
	    return -1;		// [j*M]P = 0

	i = lookup(nodes, to_ulong(*S));
        if (i != -1)
            return i+j*M;	// [j*M]P = Q - [i]P

	// don't need last addition, so skip it
        if (j == M-1 )
	    break;
    }

    // failed
    return -1;
}

///////////////////////////////////////////////////////////////////////////

// compute group size of current elliptic curve
//  - memoryless, so every invocation does calculation from scratch!!

long ell_pE::order()
{
    static zz_pE    x, y, y_sqr, z;
    static ellpoint_pE	P, Q, R, S;
    long	a_tau_q = 0, q = ell_pEInfo->q;

    // for really small finite fields, do stupid brute-force count
    if (q <= 36) {
        x = 0;
        do {
	    // y^ = (x^2 + a4)*x + a6;
            sqr(z, x);
		add(y_sqr, z, ell_pEInfo->a4);
		mul(z, y_sqr, x);
		add(y_sqr, z, ell_pEInfo->a6);

	    // if y^2=0, then trace contribution is 0
            if (y_sqr == 0)
		continue;

	    // trace contribution is -chi(x^3+a4*x+a6)
	    //  - note, this is trace on H^1, hence the minus sign
	    if (zz_pEExtraInfo->legendre_char(y_sqr) == 1)
		a_tau_q += -1;
	    else
		a_tau_q += +1;
        } while(inc(x) == NO_CARRY);

        return (q+1-a_tau_q);
    }

    // for larger fields, generate one or more points so can determine
    // order of group as the unique integer in the Hasse-Weil interval
    // [q+1-2*sqrt(q),q+1+2*sqrt(q)] satisfying certain properties.
    // 
    // - in general it suffices to find P so there is a unique multiple of
    //   of |<P>|
    // - might need to find Q so there's a unique multiple of |<P,Q>|

    long    long_two_sqrt_q = (long)floor(2.0*sqrt((double)q));
    while ((long_two_sqrt_q+1)*(long_two_sqrt_q+1) <= 2*q)
	long_two_sqrt_q++;

    // integral form of H-W interval
    long    long_l, long_r;
    long_l = q+1-long_two_sqrt_q;
    long_r = q+1+long_two_sqrt_q;

    long    HW_int = long_r-long_l;
    long    M = (long)ceil(sqrt(HW_int+1));
    long    sqrt_l = (long)ceil(sqrt(long_l));

    long    order_P = -1, order_Q, sqrt_order_P;
    unsigned long ul;

    x = 0; y = 0; z = 0; y_sqr = 0;
    do {
	// step 1: see if there is R in E(F_q) so x(R)=x

	// y^2 = (x^2+a4)*x + a6
	//  - a priori, y lives in F_{q^2} and we only know y^2
        sqr(z, x);
	    add(y_sqr, z, ell_pEInfo->a4);
	    mul(y, y_sqr, x);
	    add(y_sqr, y, ell_pEInfo->a6);

	// check to see if y lives in F_q
	ul = to_ulong(y_sqr);
	if (y_sqr != 0 && zz_pEExtraInfo->legendre_char(y_sqr) == 0)
	    continue;	// nope, so no such point

	// R=(x,y) is on our curve
	y = (y_sqr == 0) ? y_sqr : zz_pEExtraInfo->square_root(y_sqr);
        R.set(x, y);
        assert(R.on_curve());

	// step 2: check to see if order(R) <= r-l <= M^2
        long order_R = bsgs_annihilate(R, M, M*M+1);
        if (order_R == -1 || order_R > HW_int) {
	    // order(R) > r-l so there is a unique multiple in HW interval
            order_P = bsgs_annihilate(R, M, long_l + M*M);
            assert(order_P > 0 && order_P <= long_r);
            return order_P;
        }

	// step 3: if is first time we made it this far, note the point
        if (order_P == -1) {
            P = R; order_P = order_R;
            order_Q = -1;
	    continue;	// continue looking for more points
        }

	// step 4: if new point has higher order than recorded point, swap
        if (order_R > order_P) { // && order_Q == -1) {
            Q = P; order_Q = order_P;
            P = R; order_P = order_R;
	} else {
	    Q = R; order_Q = order_R;
	}
	sqrt_order_P = (long)ceil(sqrt(order_P));

	// step 5: check if order(P) >= sqrt(l)
	// - if P has maximal order, this will be satisfied
        if (order_P < sqrt_l)
	    continue; // no, so continue looking for R which works

	// step 6: determine which multiples of order_P lie in HW interval
	//
	//
	//  - is at least one; corresponds to |G|

	// step 6.a: find smallest and largest multiples
        long   min_d = (long)floor(long_l/order_P) - 1, max_d;
        while ((q+1-order_P*min_d)*(q+1-order_P*min_d) > 4*q)
            min_d++;
	max_d = min_d;
        while ((q+1-order_P*max_d)*(q+1-order_P*max_d) <= 4*q)
            max_d++;
	max_d--;

	//  claim: if order_P >= sqrt(q+1-2*sqrt(q)) and q>36, then
	//    #multiples <= 1 + (4*sqrt(q)+1)/order_P < 6
	assert(max_d - min_d < 5);

	// step 6.b: check if had exactly one multiple
        if (min_d == max_d)
	    return min_d*order_P;

	// prop: q>34 and |<P,Q>| >= q+1-2*sqrt(q) ==> |<P,Q>|/order_P > 4
	//
	// pf:
	//  - let n = max_d-min_d
	//
	//  ==> n*order(P) = (max_d-min_d)*order(P) <= 4*sqrt(q)
	//  ==> min_d >= (q+1-2*sqrt(q))/order(P)
	//            >= n*(q+1-2*sqrt(q))/(4*sqrt(q))
	//
	//  - q > 34, so min_d >= n*(q+1-2*sqrt(q))/(4*sqrt(q)) > n
	//
	//  - gcd(d1,d2) <= n for all d1<d2 in [min_d,max_d]
	//
	//  ==> suffices to show there is a unique d in [min_d,max_d] so
	//    d*Q lies in <P>
	//
	//  - let i = order(Q mod P) = |<P,Q>|/|<P>|
	//
	//  - if |<P,Q>| lies in H-W, then i >= min_d > n
	//
	//  - if d1*Q,d2*Q are both in <P> for d1<=d2 in [min_d,max_d],
	//    then i|gcd(d1,d2), hence d1=d2 or |<P,Q>|<q+1-2*sqrt(q)

	// step 7: check if m*Q lies in <P> for m=1,...,4
	//  - prop implies we can find P,Q so this doesn't happen
	int	num_d = max_d - min_d + 1;
	if (bsgs_log(P, Q, sqrt_order_P) != -1)
	    continue;
	if (num_d > 2 && order_Q % 2 == 0 && bsgs_log(P, 2*Q, sqrt_order_P) != -1)
	    continue;
	if (num_d > 3 && order_Q % 3 == 0 && bsgs_log(P, 3*Q, sqrt_order_P) != -1)
	    continue;
	if (num_d > 4 && order_Q % 4 == 0 && bsgs_log(P, 4*Q, sqrt_order_P) != -1)
	    continue;

	// step 8: determine if there is unique d in [min_d,max_d] so
	//   d*Q lies in <P>
	long    d, order_G = -1;
	num_d = 0;
	for (d = min_d; d <= max_d; d++) {
	    long    log_P = bsgs_log(P, GCD(d, order_Q)*Q, sqrt_order_P);
	    if (log_P != -1) {
		num_d++;
		order_G = d*order_P;
	    }
	}
	assert(num_d == 1);

        return order_G;
    } while(inc(x) == NO_CARRY /*&& test < MAX_POINT_TEST*/);

    // should never get here!
    assert(0);

    return -1;
}

void ell_cleanup()
{
    if (bsgs_annihilate_nodes != NULL) {
	delete bsgs_annihilate_nodes;
	bsgs_annihilate_nodes = NULL;
    }
    if (bsgs_log_nodes != NULL) {
	delete bsgs_log_nodes;
	bsgs_log_nodes = NULL;
    }
}

//////////////////////////////

#ifdef MAIN
int main(int argc, char **argv)
{
    long    p = 101, d = 1, q, i;

    if (argc == 2)
	d = atoi(argv[1]);

    init_NTL_ff(p, d, 0, 0, 0);
    q = to_ulong(zz_pE::cardinality());
    assert(q < (1L << 28));

    // precompute list of square roots
    zz_pE   sqr_roots[q], x, x2;
    do {
	x2 = x*x;
	sqr_roots[to_ulong(x2)] = x;
    } while (inc(x) == NO_CARRY);

#if 0
    {
	zz_pE   X;

	X = p-1;
	inc(X);
	cout << "X = " << X << "\n";
	cout << "X^2 = " << (X^2) << "\n";
    }
#endif

    zz_pE   one;

    // initialize the curve y^2 = x^3 + x + 1
    one = 1;
    ell_pE::init(one /*a4*/, one /*a6*/);

    zz_pE   y, y2;

    // brute force points on curve (mostly in pairs) by trying every x
    // and using model y^2 = x^3 + x + 1
    do {
	// try to compute square root of x^3 + x + 1
	y2 = (sqr(x) + ell_pEInfo->a4)*x + ell_pEInfo->a6;
	y  = sqr_roots[to_ulong(y2)];

	// zero entry in table with y2 = 0 implies no square root
	if (y2 != 0 && y == 0)
	    continue; // not a square

	cout << x << " " << y << " ";

	// initialize point on current curve
	ellpoint_pE  P(x,y);

	assert(P.on_curve());
	cout << P << " ";

	long o = P.order();
	cout << o << "\n";

	// make sure that [o]P = O
	ellpoint_pE Q;
	Q = o*P;
	assert(is_zero(Q));

	// extra sanity check.  make sure [i]P != O for 0 < i < o
	long    i;
	for (i = 1; i < o; i++)
	    assert(!is_zero(i*P));
    } while (inc(x) == NO_CARRY);

    return 0;
}
#endif // MAIN
