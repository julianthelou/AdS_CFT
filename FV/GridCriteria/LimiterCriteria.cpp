#include "GridCriteria/LimiterCriteria.h"
#include "tarch/logging/Log.h"

#include <algorithm> // transform
#include <cctype> // tolower

/************** ADMINISTRATIVE STUFF   ***********************/

/// String to lowercase, inplace.
namespace LimiterTools {
	
void toLower(std::string& data) {
	std::transform(data.begin(), data.end(), data.begin(), ::tolower);
}

} // IDtools

using namespace LimiterTools;

LimitingCriterionCode* LimitingCriterionCode::getInstanceByName(std::string name) {
	toLower(name);
	if(name == "dynamical_lapse") 	return new DynamicLapseLimiting();
	if(name == "geometric_sphere")	return new GeometricSphereLimiting();
	if(name == "geometric_ball") 	return new GeometricBallLimiting();
	return nullptr;
}

GlobalLimitingCriterion& GlobalLimitingCriterion::getInstance() {
	static GlobalLimitingCriterion* me = nullptr;
	if(!me) me = new GlobalLimitingCriterion(nullptr, "null");
	return *me;
}

bool GlobalLimitingCriterion::setIdByName(const std::string& _name) {
	static tarch::logging::Log _log("GlobalLimitingCriterion");
	name = _name;
	if(limiter) logWarning("setIdByName()", "Have already set global limtier criterion, now overwriting with "<< name << ".");
	
	limiter = LimitingCriterionCode::getInstanceByName(name);
	if(limiter) {
		logInfo("setIdByName()", "Successfully loaded "<< name << " Limiter criterion code.");
		return true;
	} else {
		logError("setIdByName()", "Requested Limiter criterion code '"<< name << "' not known.");
		return false;
	}
}

void GlobalLimitingCriterion::setByParameters(const mexa::mexafile& parameters) {
	std::string seckey = "limiter";
	std::string namekey = "criterion";
	mexa::mexafile param = parameters(seckey);
	if(!param.contains(namekey)) {
		logError("setByParameters()", "For setting up the Limiter criterion, I need a section " << seckey  << " in the parameters, as well as the key " << seckey << "/" << namekey << " to be set to a valid identifier. Instead, I got these parameters: " << parameters.toString());
		std::abort();
	}
	bool couldSetId = setIdByName(param(namekey).get_string());
	if(!couldSetId) {
		logError("setByParameters()", "Could not setup Limiter criterion code from the given parameters:" << param.toString());
		std::abort();
	} else {
		getCode().readParameters(param);
	}
}


bool GlobalLimitingCriterion::isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) {
	return getCode().isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_CALL);
}

/************** PARAMETER READING FOR LIMITING CRITERIA   ***********************/


void DynamicLapseLimiting::readParameters(const mexa::mexafile& para) {
	// TODO: check whether given or not, use default values instead of unphysical values like -1
	valid_alp_min = para["valid_alp_min"].as_double();
	valid_alp_max = para["valid_alp_max"].as_double();
	
	static tarch::logging::Log _log("DynamicLapseLimiting");
	logInfo("readParameters()", "Limiting in regions where " << valid_alp_min << " < lapse < " << valid_alp_max);
}

void GeometricSphereLimiting::readParameters(const mexa::mexafile& para) {
	radius = para["radius"].as_double();
	
	static tarch::logging::Log _log("GeometricSphereLimiting");
	logInfo("readParameters()", "Limiting at the shell with radius=" << radius);
}


void GeometricBallLimiting::readParameters(const mexa::mexafile& para) {
	static tarch::logging::Log _log("GeometricBallLimiting");
	
	radius = para["radius"].as_double();
	
	std::string where = para["where"].as_string();
	toLower(where);
	if(where == "inside") limit_inside = true;
	else if(where == "outside") limit_inside = false;
	else {
		logError("readParameters()", "Valid values for where are 'inside' and 'outside'. Keeping default.");
	}
	
	logInfo("readParameters()", "Limiting " << (limit_inside ? "within" : "outside of") << " a ball with radius=" << radius);
}

/************** THE ACTUAL LIMITING CRITERIA   ***********************/

bool DynamicLapseLimiting::isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) const {
	constexpr int NumberOfObservables = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfDMPObservables;
	if(NumberOfObservables < 1) {
		return false; // no observables!
	}
	const double lapse = observablesMin[0];
	return (valid_alp_min < lapse) && (lapse < valid_alp_max);
}

bool GeometricSphereLimiting::isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) const {
	// lower left, uppe right radius of cell
	double l = tarch::la::norm2(center - dx/2.0);
	double r = tarch::la::norm2(center + dx/2.0);
	bool isAdmissible = (l > radius) || (r <= radius);
	//printf("Cell has l=%f,r=%f => isAdmissible=%s\n", l, r, isAdmissible?"true":"false");
	return isAdmissible;
}

bool GeometricBallLimiting::isPhysicallyAdmissible(IS_PHYSICALLY_ADMISSIBLE_SIGNATURE) const {
	// lower left, uppe right radius of cell
	double l = tarch::la::norm2(center - dx/2.0);
	double r = tarch::la::norm2(center + dx/2.0);
	bool isAdmissible;
	
	if(limit_inside) {
		isAdmissible = r > radius;
	} else {
		isAdmissible = l < radius;
	}
	
	return isAdmissible;
}
