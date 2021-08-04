#include "tarch/la/Vector.h"
#include "PDE.h"
#include "Parameters/mexa.h"
#include "InitialData/InitialData.h"
#include <cstdlib> // std::abort
#include <sstream>

// use a file-global _log which does not pollute the global object namespace
namespace CCZ4Fortran_local { tarch::logging::Log _log("CCZ4Fortran"); }
using namespace CCZ4Fortran_local;

void CCZ4Fortran::tEquations::setByParameters(const mexa::mexafile& para) {
	#define SET(var) var = para[#var].as_double(); logInfo("tEquations::setByParameters(mexa)", "Setting " #var " = " << var);
	SET(k1);
	SET(k2);
	SET(k3);
	SET(eta);
	SET(itau);
	SET(f);
	SET(g);
	SET(xi);
	SET(e);
	SET(c);
	SET(mu);
	SET(ds);
	SET(sk);
	SET(bs);
	SET(LapseType);
}

std::string CCZ4Fortran::tEquations::toString() const {
	#define PRINT(var) ret << #var << ":" << var << ",";

	std::stringstream ret;
	ret << "{";
	PRINT(k1);
	PRINT(k2);
	PRINT(k3);
	PRINT(eta);
	PRINT(itau);
	PRINT(f);
	PRINT(g);
	PRINT(xi);
	PRINT(e);
	PRINT(c);
	PRINT(mu);
	PRINT(ds);
	PRINT(sk);
	PRINT(bs);
	PRINT(LapseType);
	ret << "}";
	return ret.str();
}

// Fortran functions:
extern "C" {

// PDE.f90
void pdeflux_(double* Fx, double* Fy, double* Fz, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void pdencp_(double* BgradQ, const double* const Q, const double* const gradQ);
void pdencpvector_(double* BgradQ, const double* const  Q,  const double* const Qx, const double* const Qy, const double* const Qz);
void pdesource_(double* S, const double* const Q);
void pdefusedsrcncp_(double* S, const double* const Q, const double* const gradQ);
void pdefusedsrcncpvector_(double* const S, const double* Q, const double* Qx, const double* Qy, const double* Qz);
void pdematrixb_(double* Bn, const double* Q, double* nv);

// c2p for log(alpha) and log(phi) conversion
void pdecons2prim_(double *V, const double* const Q, int* iErr);
void pdeprim2cons_(double *Q, const double* const V);

// ADMConstraints.f90
void admconstraints_(double* constraints, const double* const Q, const double* const gradQ);

// AdjustSolution.f90
void enforceccz4constraints_(double* Q);

// Psi4.f90
void calcpsi4_(const double* const Q, const double* const gradQ, double* psi0, double* psi4, double* phi0, double* phi1, double* phi2);

// Init.f90 for typesdef
void eqninit_();

// for calling abort from Fortran, always helpful:
void exahype_abort_() { std::abort(); }

// Reference to the global module storage from Fortran
extern CCZ4Fortran::tEquations typesDef_eqn;

}/* extern "C" */

CCZ4Fortran::GlobalPDEParameters::GlobalPDEParameters() : parameters(&typesDef_eqn) {}

void CCZ4Fortran::GlobalPDEParameters::setByParameters(const mexa::mexafile& mp) {
	std::string ccz4_pde_key = "ccz4";
	mexa::mexafile pde_param = mp(ccz4_pde_key);
	parameters->setByParameters(pde_param);
	logInfo("GlobalPDEParameters::setByParameters(mexa)", "CCZ4 parameters = " << parameters->toString());
}

CCZ4Fortran::GlobalPDEParameters& CCZ4Fortran::GlobalPDEParameters::getInstance() {
	static CCZ4Fortran::GlobalPDEParameters* me = nullptr;
	if(!me) me = new CCZ4Fortran::GlobalPDEParameters();
	return *me;
}


void CCZ4Fortran::PDE::nonConservativeProduct(const double* const Q, const double* const gradQ, double* BgradQ) { pdencp_(BgradQ, Q, gradQ); }
//void CCZ4Fortran::PDE::nonConservativeProductVector(const double* Q, const double* Qx, const double* Qy, const double* Qz, double* BgradQ) { pdencpvector_(BgradQ, Q, Qx, Qy, Qz); }
void CCZ4Fortran::PDE::algebraicSource(const double* const Q, double* S) { pdesource_(S, Q); }
void CCZ4Fortran::PDE::fusedSource(const double* const Q, const double* const gradQ, double* S) { pdefusedsrcncp_(S, Q, gradQ); }
//void CCZ4Fortran::PDE::fusedSourceVector(double* Q, double* Qx, double* Qy, double* Qz, double* S) { pdefusedsrcncpvector_(S, Q, Qx, Qy, Qz); }
void CCZ4Fortran::PDE::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
	double nv[3] = {0.};
	nv[dIndex] = 1;
	pdeeigenvalues_(lambda, Q, nv);
}

struct Cons2Prim {
	Cons2Prim(double* const V, const double* const Q, bool crash_on_failure=true);
	bool failed; ///< Stores whether the C2P failed or not
};
void Prim2Cons(double* const Q, const double* const V);
	
void CCZ4Fortran::Cons2Prim(double* const V, const double* const Q) {
	int iErr;
	//pdecons2prim_(V, Q, &iErr);
	// The CCZ4 C2P can never fail. It's not "real" conserved variables. I don't like
	// calling it conserved here, either.
}
void CCZ4Fortran::Prim2Cons(double* const Q, const double* const V) {
	//pdeprim2cons_(Q, V);
}

void CCZ4Fortran::DeriveADMConstraints(double* constraints, const double* const Q, const double* const gradQ) {
	admconstraints_(constraints, Q, gradQ);
}

void CCZ4Fortran::EnforceCCZ4Constraints(double *Q) {
	enforceccz4constraints_(Q);
}

void CCZ4Fortran::AdjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
	using namespace CCZ4Fortran;
	if(tarch::la::equals(t,0.0) ) {
		CCZ4_InitialData::InitialData(x,t,Q);
	} else {
		EnforceCCZ4Constraints(Q);
	}
}

void CCZ4Fortran::CalcPsi4(const double* const Q, const double* const gradQ, double* psi0, double* psi4, double* phi0, double* phi1, double* phi2) {
	double psi0_[2], psi4_[2], phi0_[2], phi1_[2], phi2_[2];
	
	// Allow optional return values (set nullptr if you are not interested in these quantities)
	if(psi0==nullptr) psi0 = psi0_;
	if(psi4==nullptr) psi4 = psi4_;
	if(phi0==nullptr) phi0 = phi0_;
	if(phi1==nullptr) phi1 = phi1_;
	if(phi2==nullptr) phi2 = phi2_;
	
	calcpsi4_(Q,gradQ,psi0,psi4,phi0,phi1,phi2);
}
