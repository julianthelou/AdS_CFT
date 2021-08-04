#ifndef __EXAHYPE_ASTROPHYSICS_CCZ4_SVENS_SOMMERFELD_BOUNDARY_CONDITION_CODE__
#define __EXAHYPE_ASTROPHYSICS_CCZ4_SVENS_SOMMERFELD_BOUNDARY_CONDITION_CODE__

#include "AbstractCCZ4Solver_ADERDG.h"
#include "tarch/la/Vector.h"

namespace CCZ4_RadBC {
	using BaseSolverType = CCZ4::AbstractCCZ4Solver_ADERDG;

	/**
	 * Describe the PDE to be solved at the boundary in an ExaHyPE-Solver manner.
	 * The PDE is very simple:
	 * 
	 *   partial_t f + v r^i / r partial_i f = - v (f-f_0) / r
	 * 
	 * This PDE is purely non-conservative (has zero conserved flux), so it has a NCP and
	 * an algebraic Source. See also the EinsteinDoc directory for further information.
	 * 
	 * In our code, we call f_0=Q0 and v=v0. Obviously, we need the radius as an additional information.
	 * 
	 * Note that the Q0 are the field values at infinity which have to be initialized for all fields
	 * at startup, similar to NewRad initialization in Cactus.
	 * 
	 * The routines in this class read like in the solver API in ExaHyPE.
	 * 
	 * Important: Note the subtle difference:  Here, the radius (r) are passed, as the CCZ4
	 *     state vector, as written in the literature, does not hold them (by intention).
	 **/
	struct RadiativeWavePDE {
		// The standard constexpr's needed by the kernels
		static constexpr int NumberOfVariables = BaseSolverType::NumberOfVariables;
		static constexpr int NumberOfParameters = BaseSolverType::NumberOfParameters;
		static constexpr int Order = BaseSolverType::Order;
		
		/// The constructor sets the default FO-CCZ4 values for v0, Q0.
		RadiativeWavePDE();
		
		// Since all other Fortran PDEs are global, use this for singleton access.
		// (You are not forced to do so)
		static RadiativeWavePDE& getSingleton();
		
		double v0[NumberOfVariables] = {0}; ///< wave speeds
		double Q0[NumberOfVariables] = {0}; ///< field values at infinity
		
		/// A schmankerl (syntactic sugar) for our Cactus friends, reads like
		///   newrad(cctkGH, varName, rhsName, var0, v0) in Cactus newRad thorn.
		void newrad(int pos, double infval, double vel) { Q0[pos] = infval; v0[pos] = vel; }
		
		/// This method is a bit awkward: It actually computes the BgradQ, but calls it flux,
		///   because it is used as such in the boundary treatment. Don't be confused: The
		///   conserved flux is zero.
		void boundaryFluxContribution(const double* const x, const double* const Q, double** F) const;
		
		/// Will just return a copy of v0.
		void eigenvalues(const double* const x, const double* const Q, double* Lambda) const;
		
		/// In ExaHyPE slang, this is called the fusedSource.
		void rightHandSide(const double* const x, const double* const Q, const double* const gradQ, double* S) const;

		/// The regular NCP as in ExaHyPE
		void nonConservativeProduct(const double* const x, const double* const Q, const double* const gradQ, double* BgradQ) const;
		
		/// The regular algebraicSource as in ExaHyPE
		void algebraicSource(const double* const x, const double* const Q, double* S) const;
	};
	
	/**
	 * A variant of the radiative BC PDE which is more comfortable with CCZ4:
	 *   - Uses no fluxes (but instead the NCP only)
	 *   - Provides a distinct algebraic Source
	 **/
	
	/**
	 * Diffusive Interface Method: Coupled FO-CCZ4 and BohrSommerfeld
	 * An indicator field is used to tell where the physical domain is (where to solve CCZ4) and where
	 * the outer domain is (where to solve Bohr Sommerfeld).
	 **/
	namespace DIM_CCZ4 {
		void getPosition(const double* const Q, double* x);
		void mix(const double* const Q, const double* const CCZ4_Solution, const double* const BohrSommerfeld_Solution, double* Mixed_Solution);
		
		void nonConservativeProduct(const double* const Q, const double* const gradQ, double* BgradQ);
		void algebraicSource(const double* const Q, double* S);
		void fusedSource(const double* const Q, const double* const gradQ, double* S);
		void eigenvalues(const double* const Q, const int dIndex, double* lambda);
	}
	
	/**
	 * A call to this function replaces the call to
	 * 
	 *   kernels::aderdg::generic::c::boundaryConditions<CCZ4Solver_ADERDG>(*static_cast<CCZ4Solver_ADERDG*>(this),fluxOut,stateOut,fluxIn,stateIn,cellCentre,cellSize,t,dt,faceIndex,direction);
	 *
	 * in the AbstractCCZ4Solver.cpp, in the function
	 *   void CCZ4::AbstractCCZ4Solver_ADERDG::boundaryConditions(....)
	 *
	 * What is does is to evaluate the BC for the whole patch using the
	 * Integrator below instead of computing BC point by point.
	 * 
	 * Notes to keep in mind:
	 * 
	 * 1) there is *No* call to the CCZ4 PDE system here but instead to the
	 *    RadiativeWavePDE.
	 * 2) As we want to provide BC for the ADER-DG scheme, we need to provide
	 *    the state at the boundary surface (QOut) and the fluxes at the boundary
	 *    surface (FOut).
	 * 
	 **/
	struct RadiativeBoundaryConditions {
		using SolverType = RadiativeWavePDE; // no need for templates here.
		SolverType pde;
		
		/**
		 * Slightly modified signature from the kernel call:
		 * 
		 *  - luh: Volumetric cell solution DOF
		 *  - then comes the rest
		 **/
		void boundaryConditions(
			double* luh, // volumetric solution!
			double* FOut,
			double* QOut,
			const double* const FIn,
			const double* const QIn,
			const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
			const tarch::la::Vector<DIMENSIONS,double>& cellSize,
			const double t,const double dt,
			const int faceIndex,
			const int normalNonZero);
	};
} // namespace CCZ4_BC

#endif /* __EXAHYPE_ASTROPHYSICS_CCZ4_SVENS_SOMMERFELD_BOUNDARY_CONDITION_CODE__ */
