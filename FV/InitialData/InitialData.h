#ifndef __EXAHYPE_CCZ4_INITIAL_DATA__
#define __EXAHYPE_CCZ4_INITIAL_DATA__

#include "AbstractCCZ4Solver_ADERDG.h"
	
#include <stdio.h>
#include "Parameters/mexa.h"

namespace CCZ4_InitialData {
	struct InitialDataCode; // forward definition
	class GlobalInitialData; // forward definition
	struct ADMBase; // ADMBase.h
}

namespace CCZ4_InitialData {
/**
 * The Initial Data code interface.
 * 
 * This structure is similar to the Initial Data interface in the GRMHD application
 * and all CCZ4-specific or Gradient-based-specific details have been moved to the
 * GradientBasedInitialDataCode class.
 **/
class InitialDataCode {
public:
	virtual void readParameters(const mexa::mexafile& parameters) {} ///< by default do nothing.
	virtual void prepare() {} ///< Run a ID code. Parameters should be set before.

	/**
	 * This is the interface to an actual initial data code. They won't set
	 * BSSNOK, SO-CCZ4 or FO-CCZ4 quantities, all this has to be computed in a post-step.
	 **/
	virtual void Interpolate(const double* const x, double t, ADMBase& out) = 0;
	
	/**
	 * ID codes can post-process the computed state vector by using this method.
	 * It is called after all auxilliary parameters were computed, on every point.
	 * This is suitable for setting the Diffusive Interface.
	 * By default, does nothing.
	 **/
	virtual void PostProcessStateVector(const double* const x, double t, double* Q) {}
	
	/**
	 * This is the user frontend for the classical adjustPointSolution function
	 * which spills out a FO-CCZ4 state vector. It computes the auxilliary
	 * quantities by means of a high order FD stencil around x which requires
	 * several evaluations of Interpolate().
	 **/
	virtual void adjustPointSolution(const double* const x,const double t,double* Q);
	
	/**
	 * This is the user frontend for the patch wise evaluation where a whole cell
	 * (all DOF in a patch) is filled with FO-CCZ4 state vectors. This is done
	 * internally by looping over the grid, computing the SO-CCZ4 quantities and
	 * then computing the DG derivatives for the FO-CCZ4 quantities. In sum,
	 * this yields a faster evaluation (less calls to Interpolate) and
	 * potentially more consistent data (real DG derivatives).
	 **/
	void adjustPatchSolution(
		const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
		const tarch::la::Vector<DIMENSIONS, double>& dx,
		const double t, const double dt, double* luh);
	
	/**
	 * Set the "inside" default for the diffusive interface indicator, as well
	 * as store the coordinates in the state vector.
	 **/
	void setExaHyPEspecificQuantities(const double* const x, double* Q);
	
	/// This is the actual initial data registration, it is implemented by us.
	/// TODO: This should at sometime move to another class. Is not really good here.
	static InitialDataCode* getInstanceByName(std::string name);
	
private:
	/**
	* Internal, but here for documentation:
	* Do the a determination of gradQ by Interpolation with finite differences
	**/
	void FD_Interpolate(const double* const x, const double t, double **gradQ);
};

// a shorthand for pointSolutions.
void InitialData(const double* x, double t, double* Q);

/**
 * Global initial data, a singleton for the time when there is only (logically)
 * a single solver instance in the application. That means especially when there
 * is a single set of Limiting/FD/DG solvers, the singleton ensures that the ID
 * are initialized only once and not for every of the two FD/DG solvers.
 * That typically means the earlier initialized solver also initializes the
 * ID and the later solver just uses it.
 **/
class GlobalInitialData {
	tarch::logging::Log _log;
	bool alreadySetParameters; ///< internal flag to avoid double setting parameters
	bool alreadyPrepared; ///< internal flag to avoid calling preparation routine several times
public:
	/// The pointer is by default null, meaning no ID can be used. The setIdByName method
	/// sets this id pointer, finally.
	InitialDataCode* id;
	std::string name; ///< set by setIdByName
	
	GlobalInitialData(InitialDataCode* _id, std::string _name)
		: _log("GlobalInitialData"), alreadySetParameters(false), alreadyPrepared(false), id(_id), name(_name) {}
	
	/// Returns true in case of success
	bool setIdByName(const std::string& name);
	bool alreadySetByName() { return (id != nullptr); }
	
	/// Tries to read ID from parameters. In case of problems, raises exceptions.
	GlobalInitialData& setByParameters(const mexa::mexafile& parameters);
	
	/// Gives a singleton instance
	static GlobalInitialData& getInstance();
	
	// as a shorthand
	static InitialDataCode& getInitialDataCode() { return *(getInstance().id); }
	
	// Run the initial data code.
	GlobalInitialData& prepare();
	
	// Ensure the initial data code runned already.
	static bool ensureGlobalIDareSet(bool doCrash=true);
};

} // namespace CCZ4_InitialData

#include "ToyID.h"
#include "TwoPunctures_Binding.h"
#include "UIUC_Binding.h"

#endif /* __EXAHYPE_CCZ4_INITIAL_DATA__ */
