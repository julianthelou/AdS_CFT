#ifndef __EXAHYPE_CCZ4_PDE_FORTRAN__
#define __EXAHYPE_CCZ4_PDE_FORTRAN__

// forward declaration instead of #include "Parameters/mexa.h"
namespace mexa { class mexafile; }

namespace CCZ4Fortran {

	/**
	* The CCZ4 PDE system runtime parameters.
	* This structure has to been kept in sync with the Fortran user type in typesDef.f90
	**/
	struct tEquations {
		double k1;        // CCZ4 damping parameter k1 
		double k2;        // CCZ4 damping parameter k2  
		double k3;        // CCZ4 damping parameter k3 
		double eta;       // CCZ4 damping parameter for the PDE of b^i in the gamma driver 
		double itau;      // inverse relaxation time for the relaxed approach in order to enforce the unit determinant of the conformal metric and the trace-free matrix A_ij
		double f;         // set f=0.75 or f=1.0 for the gamma driver. Typical BSSNOK value: f=0.75. Set f=0 to avoid shift evolution, i.e. to have d/dt beta^i=0.  
		double g;         // not used at the moment: reserved for the generalized harmonic shift   
		double xi;        // set to zero to switch off the gamma driver and to avoid shift evolution (i.e. to have d/dt b^i = 0) 
		double e;         // cleaning speed e>=1 for the Hamiltonian constraint. Typical value for CCZ4 is e=1. However, e>1 gives better constraints and better hyperbolicity//  
		double c;         // set c=0 to remove the algebraic source terms of the type -2*Theta 
		double mu;        // mu=1 adds second order ordering constraints. Has shown to be important for the gamma driver, but too large values can also become a problem... 
		double ds;        // set this value always to ds=1, unless you know what you are doing. It allows to increase the cleaning speed for the momentum constraints, but in CCZ4 this does not seem to do what it should do...  
		double sk;        // setting sk=0 removes the contribution of the shift also in the auxiliary variables. If you want to evolve the shift, set sk=1. 
		double bs;        // set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
		int LapseType;    // LapseType = 0 is harmonic lapse, LapseType = 1 is 1+log slicing.
		
		void setByParameters(const mexa::mexafile& parameters);
		std::string toString() const;
	};

	struct GlobalPDEParameters {
		tEquations* const parameters;
		GlobalPDEParameters();

		/// Tries to parameters from mexa. In case of problems, raises exceptions.
		void setByParameters(const mexa::mexafile& parameters);
	
		/// Gives a singleton instance
		static GlobalPDEParameters& getInstance();
	};

	// it is pointless to have this being a class as the methods only mask Fortran functions.
	namespace PDE {
		void nonConservativeProduct(const double* const Q, const double* const gradQ, double* BgradQ);
		void nonConservativeProductVector(const double* const Q, const double* const Qx, const double* const Qy, const double* const Qz, double* BgradQ);
		void algebraicSource(const double* const Q, double* S);
		void fusedSource(const double* const Q, const double* const gradQ, double* S);
		void fusedSourceVector(double*   Q,  double*   Qx,  double*   Qy,  double*   Qz, double*   S);
		void eigenvalues(const double* const Q, const int dIndex, double* lambda);
	}

	void Cons2Prim(double* const V, const double* const Q);
	void Prim2Cons(double* const Q, const double* const V);
	void EnforceCCZ4Constraints(double *Q);
	void DeriveADMConstraints(double* constraints, const double* const Q, const double* const gradQ);
	void AdjustPointSolution(const double* const x,const double t,const double dt,double* Q);
	void CalcPsi4(const double* const Q, const double* const gradQ, double* psi0, double* psi4, double* phi0, double* phi1, double* phi2);
	
	
	// Sequentialization which really goes down for each tensor to the individiual position.
	// In contrast, in CCZ4Solver_{ADERDG,FV}_Variables.h, everything is a 1d vector.
	namespace StateVectorSequentialization {
		constexpr int g11    =  0;
		constexpr int g12    =  1;
		constexpr int g13    =  2;
		constexpr int g22    =  3;
		constexpr int g23    =  4;
		constexpr int g33    =  5;
		constexpr int A11    =  6;
		constexpr int A12    =  7;
		constexpr int A13    =  8;
		constexpr int A22    =  9;
		constexpr int A23    = 10;
		constexpr int A33    = 11;
		constexpr int Theta  = 12;
		constexpr int G1     = 13;
		constexpr int G2     = 14;
		constexpr int G3     = 15;
		constexpr int lapse  = 16;
		constexpr int shift1 = 17;
		constexpr int shift2 = 18;
		constexpr int shift3 = 19;
		constexpr int b1     = 20;
		constexpr int b2     = 21;
		constexpr int b3     = 22;
		constexpr int A1     = 23;
		constexpr int A2     = 24;
		constexpr int A3     = 25;
		constexpr int B11    = 26;
		constexpr int B21    = 27;
		constexpr int B31    = 28;
		constexpr int B12    = 29;
		constexpr int B22    = 30;
		constexpr int B32    = 31;
		constexpr int B13    = 32;
		constexpr int B23    = 33;
		constexpr int B33    = 34;
		constexpr int D111   = 35;
		constexpr int D112   = 36;
		constexpr int D113   = 37;
		constexpr int D122   = 38;
		constexpr int D123   = 39;
		constexpr int D133   = 40;
		constexpr int D211   = 41;
		constexpr int D212   = 42;
		constexpr int D213   = 43;
		constexpr int D222   = 44;
		constexpr int D223   = 45;
		constexpr int D233   = 46;
		constexpr int D311   = 47;
		constexpr int D312   = 48;
		constexpr int D313   = 49;
		constexpr int D322   = 50;
		constexpr int D323   = 51;
		constexpr int D333   = 52;
		constexpr int K      = 53;
		constexpr int phi    = 54;
		constexpr int P1     = 55;
		constexpr int P2     = 56;
		constexpr int P3     = 57;
		constexpr int K0     = 58;
		constexpr int domain = 59;
		constexpr int posx   = 60;
		constexpr int posy   = 61;
		constexpr int posz   = 62;
	} // ns StateVectorSequentialization
} // ns CCZ4Fortran

#endif /* __EXAHYPE_CCZ4_PDE_FORTRAN__ */
