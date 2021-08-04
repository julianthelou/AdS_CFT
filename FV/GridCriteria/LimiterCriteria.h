#ifndef __CCZ4_LIMITER_CRITERIA__
#define __CCZ4_LIMITER_CRITERIA__

#include "AbstractCCZ4Solver_ADERDG.h"
#include "Parameters/mexa.h"


#define IS_PHYSICALLY_ADMISSIBLE_SIGNATURE \
  const double* const solution, \
  const double* const observablesMin,const double* const observablesMax, \
  const bool wasTroubledInPreviousTimeStep, \
  const tarch::la::Vector<DIMENSIONS,double>& center, \
  const tarch::la::Vector<DIMENSIONS,double>& dx, \
  const double t

#define IS_PHYSICALLY_ADMISSIBLE_CALL \
	solution, observablesMin, observablesMax, wasTroubledInPreviousTimeStep, center, dx, t


struct LimitingCriterionCode {
	virtual bool isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) const { std::abort(); /* pls implement */ }
	virtual void readParameters(const mexa::mexafile& parameters) { /* by default do nothing. */ }

	/// This is the actual initial data registration
	static LimitingCriterionCode* getInstanceByName(std::string name);
};

// here are some codes:

struct NoLimiting : public LimitingCriterionCode {
	bool isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) const override {
		return true;
	}
};

struct DynamicLapseLimiting : public LimitingCriterionCode {
	double valid_alp_min; ///< valid minimum lapse, default -1  (= all valid)
	double valid_alp_max; ///< valid maximum lapse, default 1e100 (= all valid)
	
	DynamicLapseLimiting() : valid_alp_min(-1), valid_alp_max(1e100) {}
	
	bool isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) const override;
	void readParameters(const mexa::mexafile& parameters) override;
};

struct GeometricSphereLimiting : public LimitingCriterionCode {
	double radius; ///< Radius of the sphere where to limit, default -1 (= no limiting)
	
	GeometricSphereLimiting() : radius(-1) {}
	
	bool isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) const override;
	void readParameters(const mexa::mexafile& parameters) override;
};

struct GeometricBallLimiting : public LimitingCriterionCode {
	double radius; ///< Radius of the ball (default -1 = no limiting)
	bool limit_inside; ///< Whether to limit inside or outside (default inside)
	
	GeometricBallLimiting() : radius(-1), limit_inside(true) {}
	
	bool isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) const override;
	void readParameters(const mexa::mexafile& parameters) override;
};


class GlobalLimitingCriterion {
	tarch::logging::Log _log;
public:
	/// The pointer is by default null, meaning no code can be used. The setCodeByName method
	/// sets the limiter pointer, finally.
	LimitingCriterionCode* limiter;
	std::string name; ///< set by setIdByName
	
	GlobalLimitingCriterion(LimitingCriterionCode* _limiter, std::string _name) :
		_log("GlobalLimitingCriterion"), limiter(_limiter), name(_name) {}
	
	/// Returns true in case of success
	bool setIdByName(const std::string& name);
	
	/// Tries to read ID from parameters. In case of problems, raises exceptions.
	void setByParameters(const mexa::mexafile& parameters);
	
	/// Gives a singleton instance
	static GlobalLimitingCriterion& getInstance();
	
	// as a shorthand
	static LimitingCriterionCode& getCode() { return *(getInstance().limiter); }
	
	// as a shorthand
	static bool isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE); 
};


#endif /* __CCZ4_LIMITER_CRITERIA__ */
