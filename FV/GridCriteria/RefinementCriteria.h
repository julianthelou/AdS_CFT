#ifndef __CCZ4_REFINEMENT_CRITERIA__
#define __CCZ4_REFINEMENT_CRITERIA__

#include "AbstractCCZ4Solver_ADERDG.h"
#include "Parameters/mexa.h"

#define REFINEMENT_CONTROL_SIGNATURE \
	const double* luh, \
	const tarch::la::Vector<DIMENSIONS,double>& center, \
	const tarch::la::Vector<DIMENSIONS,double>& dx, \
	double t, \
	const int level
#define REFINEMENT_CONTROL_CALL \
	luh, center, dx, t, level
#define REFINEMENT_CONTROL_RETURN \
	exahype::solvers::Solver::RefinementControl

struct RefinementCriterionCode {
	virtual REFINEMENT_CONTROL_RETURN refinementCriterion(REFINEMENT_CONTROL_SIGNATURE) const { std::abort(); /* pls implement */ }
	virtual void readParameters(const mexa::mexafile& parameters) { /* by default do nothing. */ }

	/// This is the actual initial data registration
	static RefinementCriterionCode* getInstanceByName(std::string name);
};

// here are some codes:

struct NoRefinement : public RefinementCriterionCode {
	REFINEMENT_CONTROL_RETURN refinementCriterion(REFINEMENT_CONTROL_SIGNATURE) const override;
};

struct GeometricSphereRefinement : public RefinementCriterionCode {
	double radius; ///< Radius of the sphere where to limit, default -1 (= no limiting)
	
	GeometricSphereRefinement() : radius(-1) {}
	
	REFINEMENT_CONTROL_RETURN refinementCriterion(REFINEMENT_CONTROL_SIGNATURE) const override;
	void readParameters(const mexa::mexafile& parameters) override;
};

class GlobalRefinementCriterion {
	tarch::logging::Log _log;
public:
	/// The pointer is by default null, meaning no code can be used. The setCodeByName method
	/// sets the refiner pointer, finally.
	RefinementCriterionCode* refiner;
	std::string name; ///< set by setIdByName
	
	GlobalRefinementCriterion(RefinementCriterionCode* _refiner, std::string _name) :
		_log("GlobalRefinementCriterion"), refiner(_refiner), name(_name) {}
	
	/// Returns true in case of success
	bool setIdByName(const std::string& name);
	
	/// Tries to read ID from parameters. In case of problems, raises exceptions.
	void setByParameters(const mexa::mexafile& parameters);
	
	/// Gives a singleton instance
	static GlobalRefinementCriterion& getInstance();
	
	// as a shorthand
	static RefinementCriterionCode& getCode() { return *(getInstance().refiner); }
	
	// as a shorthand
	static REFINEMENT_CONTROL_RETURN refinementCriterion(REFINEMENT_CONTROL_SIGNATURE);
};


#endif /* __CCZ4_REFINEMENT_CRITERIA__ */
