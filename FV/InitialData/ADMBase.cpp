#include "ADMBase.h"
#include "InitialData.h"

#include "CCZ4Solver_FV_Variables.h"
// or:   "CCZ4Solver_ADERDG_Variables.h", as you prefer. Then adopt the aliases below.
namespace idx = CCZ4::CCZ4Solver_FV_Variables::shortcuts;
constexpr int nVars = CCZ4::AbstractCCZ4Solver_FV::NumberOfVariables;

#include "kernels/aderdg/generic/c/computeGradients.cpph" // computeGradQ

#define TDIM DIMENSIONS
#include "SVEC/tensish.cpph"

#include <cmath>
#include <exception>


using namespace CCZ4_InitialData;
using namespace tensish;

typedef tensish::metric::generic<square::const_shadow_D> our_metric;

void CCZ4_InitialData::ADMBase::zero() {
	alpha = 0;
	for(int i=0; i<3; i++) {
		beta_lo[i] = 0;
		for(int j=0; j<3; j++) {
			gamma_lo[i][j] = 0;
			Kextr_lo[i][j] = 0;
		}
	}
}

void CCZ4_InitialData::ADMBase::one() {
	// Minkowski flat space
	alpha = 1.0;
	for(int i=0; i<3; i++) {
		beta_lo[i] = 0;
		for(int j=0; j<3; j++) {
			gamma_lo[i][j] = sym::delta(i,j);
			Kextr_lo[i][j] = 0; //sym::delta(i,j);
		}
	}
}

void CCZ4_InitialData::ADMBase::copyToStateVector(double* Q) const {
	// Export C matrices in Fortran row-major ordering	
	Q[idx::G + 0] = gamma_lo[0][0];
	Q[idx::G + 1] = gamma_lo[0][1];
	Q[idx::G + 2] = gamma_lo[0][2];
	Q[idx::G + 3] = gamma_lo[1][1];
	Q[idx::G + 4] = gamma_lo[1][2];
	Q[idx::G + 5] = gamma_lo[2][2];
	
	Q[idx::K + 0] = Kextr_lo[0][0];
	Q[idx::K + 1] = Kextr_lo[0][1];
	Q[idx::K + 2] = Kextr_lo[0][2];
	Q[idx::K + 3] = Kextr_lo[1][1];
	Q[idx::K + 4] = Kextr_lo[1][2];
	Q[idx::K + 5] = Kextr_lo[2][2];

	Q[idx::lapse] = alpha;
	
	Q[idx::shift + 0] = beta_lo[0]; 
	Q[idx::shift + 1] = beta_lo[1];
	Q[idx::shift + 2] = beta_lo[2];
}

void CCZ4_InitialData::ADMBase::copyFromStateVector(const double* const Q) {
	// Import C matrices from Fortran row-major ordering
	gamma_lo[0][0] = Q[idx::G + 0];
	gamma_lo[0][1] = Q[idx::G + 1];
	gamma_lo[0][2] = Q[idx::G + 2];
	gamma_lo[1][1] = Q[idx::G + 3];
	gamma_lo[1][2] = Q[idx::G + 4];
	gamma_lo[2][2] = Q[idx::G + 5];
	
	gamma_lo[1][0] = gamma_lo[0][1];
	gamma_lo[2][0] = gamma_lo[0][2];
	gamma_lo[2][1] = gamma_lo[1][2];

	Kextr_lo[0][0] = Q[idx::K + 0];
	Kextr_lo[0][1] = Q[idx::K + 1];
	Kextr_lo[0][2] = Q[idx::K + 2];
	Kextr_lo[1][1] = Q[idx::K + 3];
	Kextr_lo[1][2] = Q[idx::K + 4];
	Kextr_lo[2][2] = Q[idx::K + 5];
	
	Kextr_lo[1][0] = Kextr_lo[0][1];
	Kextr_lo[2][0] = Kextr_lo[0][2];
	Kextr_lo[2][1] = Kextr_lo[1][2];

	alpha = Q[idx::lapse];
	
	beta_lo[0] = Q[idx::shift + 0];
	beta_lo[1] = Q[idx::shift + 1];
	beta_lo[2] = Q[idx::shift + 2];
}

void CCZ4_InitialData::ADMBase::importSOCCZ4(double* Q) const {
	const our_metric gam(&gamma_lo[0][0]);
	
	// conformal factor
	//% phi = 1/cbrt(sqrt(det)) = det^{1/2}^{1/3}^{-1} = det^{-1/6}
	double phi = std::pow(gam.det, -1./6.); //1./std::cbrt(gam.sqdet);
	double phi2 = phi*phi;
	//% psi4 = det^{1/3} = 1/phi**2
	double psi4 = std::cbrt(gam.det);
	
	// logs are evolved of these two scalar fields
	double logW = std::log(phi);
	double logLapse = std::log(alpha);

	double traceK = 0; CONTRACT2(i,j) traceK += gam.up(i,j) * Kextr_lo[i][j];
	
	for(int i=0; i<nVars; i++) Q[i] = 0; // safety, could also set NaNs.
	copyToStateVector(Q); // this does all F<->C ordering for us
	
	DOSYM(ij) Q[idx::G + ij] *= phi2;
	DOSYM(ij) Q[idx::K + ij]  = phi2*Q[idx::K + ij] - 1./3. * traceK * Q[idx::G + ij];

	Q[idx::traceK] = traceK;
	Q[idx::lapse] = logLapse;
	Q[idx::phi] = logW;
	
	// idx::Z is Gtilde and is set with the help of the gradients.
	Q[idx::theta] = 0;
	DFOR(i) Q[idx::Z + i] = 0;
	Q[idx::K0] = 0;
}

void CCZ4_InitialData::adjustAuxilliaries(double *Q, const double* const gradQ[3]) {
	// gradQ are the gradients of the "conserved" state vector (logLapse, logW)
	// and not the "primitive" one, therefore the conversions from pure derivatives
	// to auxillaries differ from Dumbsers code.
	
	//% A_i = (partial_i lapse) / lapse = \partial_i \log(lapse)
	DFOR(i) Q[idx::dLapse + i] = gradQ[i][idx::lapse];
	DFOR(i) Q[idx::P + i] = gradQ[i][idx::phi];

	//% B^i_k = partial_k beta^i
	DFOR(i) Q[idx::dxShift + i] = gradQ[0][idx::shift + i];
	DFOR(i) Q[idx::dyShift + i] = gradQ[1][idx::shift + i];
	DFOR(i) Q[idx::dzShift + i] = gradQ[2][idx::shift + i];
	
	//% D_{ijk} = 0.5 partial_k gammat_ij
	DOSYM(ij) Q[idx::dxG + ij] = 0.5 * gradQ[0][idx::G + ij];
	DOSYM(ij) Q[idx::dyG + ij] = 0.5 * gradQ[1][idx::G + ij];
	DOSYM(ij) Q[idx::dzG + ij] = 0.5 * gradQ[2][idx::G + ij];
	
	// Without the ADMBase quantities accessible, we need to recover not only the
	// the conformal factor (cons2prim) and the metric inverse, but also at the first
	// place the regular 3-metric from the conformal one :-(
	/*
	// Only neccessary for adopting Dumbsers equation
	ADMBase state;
	state.copyFromStateVector(Q);
	const double phi = std::exp(Q[idx::phi]); // Boxing of phi would be applied here (if important).
	const double psi4 = 1./(phi*phi);
	DFOR(i) DFOR(j) state.gamma_lo[i][j] *= psi4; // restore regular 3-metric from conformal one
	const our_metric g(&state.gamma_lo[0][0]);
	*/
	
	// Same thing but with conformal metric: Recomputing \tilde\gamma^{ij}
	ADMBase state2;
	state2.copyFromStateVector(Q);
	const our_metric gammat(&state2.gamma_lo[0][0]);
	
	// The way to avoid this double work would be an extended ADMBase storing
	// phi, psi4, g_up, g_lo. However, come on, this is only reading in InitialData.
	
	// Access DD_{kij} as DD(k)(i,j)
	const tensish::generic::stored<square::const_shadow<sym::row_major, DIMENSIONS>, DIMENSIONS>
		DD { Q+idx::dxG, Q+idx::dyG, Q+idx::dzG };
	const vec::const_shadow_D PP(Q+idx::P);

	// Gamma_d^i = gamma^{jk} Christoffel^i_{jk}
	vec::shadow_D Gtilde(Q+idx::Z);
	Gtilde = 0;

	// We cannot use exactly the equation from Dumbsers code since there, he uses \partial_l gamma_{jk}
	// but we only have \partial_l \tilde\gamma_{jk}. Therefore, instead...
	/*
	// Fortran equivalent:
	// Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  )
	DFOR(i) CONTRACT3(j,k,l)
	Gtilde(i) += psi4*g.up(i,j)*g.up(k,l)*(2*DD(l)(j,k) + 2*PP(l)*g.lo(j,k));
	*/
	// ... we adopt an equation from our paper (Dumbser+2017)
	DFOR(i) CONTRACT3(j,k,l) Gtilde(i) += gammat.up(i,j)*gammat.up(k,l)*2*DD(l)(j,k);
	
	// For debugging
	/*
	printf("Cphi=%.20e\n", phi);
	DFOR(i) DFOR(j) printf("Cgup[%d,%d]=%.20e\n", i, j, g.up(i,j));
	DFOR(i) DFOR(j) printf("Cglo[%d,%d]=%.20e\n", i, j, g.lo(i,j));
	DFOR(k) DFOR(i) DFOR(j) printf("CDD[%d,%d,%d]=%.20e\n", k, i, j, DD(k)(i,j));
	DFOR(i) printf("CPP[%d]=%.20e\n", i, PP(i));
	DFOR(i) printf("CGtilde[%d]=%.20e\n", i, Gtilde(i));
	
	if(std::getenv("ABORT_C")) std::abort();
	//std::abort(); // debug me!
	*/
}


void CCZ4_InitialData::ADMBase::importStateVector(const double* const Q) {
	copyFromStateVector(Q); // this does all F<->C ordering for us
	alpha = std::exp(alpha);
	double phi = std::exp(Q[idx::phi]);
	double phi2 = phi*phi;
	DFOR(i) DFOR(j) gamma_lo[i][j] /= phi2; // conformal -> non-conformal
	
	// TODO: Check wether this is right:
	DFOR(i) DFOR(j) Kextr_lo[i][j] = Kextr_lo[i][j] / phi2 + 1./3.* gamma_lo[i][j] * Q[idx::traceK];
}

#ifdef COMPILES	
#endif
