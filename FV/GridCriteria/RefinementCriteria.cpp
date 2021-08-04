#include "GridCriteria/RefinementCriteria.h"
#include "tarch/logging/Log.h"

#include <algorithm> // transform
#include <cctype> // tolower

/************** ADMINISTRATIVE STUFF   ***********************/

/// String to lowercase, inplace.
namespace RefinerTools {
	
void toLower(std::string& data) {
	std::transform(data.begin(), data.end(), data.begin(), ::tolower);
}

} // IDtools

using namespace RefinerTools;

RefinementCriterionCode* RefinementCriterionCode::getInstanceByName(std::string name) {
	toLower(name);
	if(name == "no_refinement") 	return new NoRefinement();
	if(name == "geometric_sphere")	return new GeometricSphereRefinement();
	return nullptr;
}

GlobalRefinementCriterion& GlobalRefinementCriterion::getInstance() {
	static GlobalRefinementCriterion* me = nullptr;
	if(!me) me = new GlobalRefinementCriterion(nullptr, "null");
	return *me;
}

bool GlobalRefinementCriterion::setIdByName(const std::string& _name) {
	static tarch::logging::Log _log("GlobalRefinementCriterion");
	name = _name;
	if(refiner) logWarning("setIdByName()", "Have already set global refinement criterion, now overwriting with "<< name << ".");
	
	refiner = RefinementCriterionCode::getInstanceByName(name);
	if(refiner) {
		logInfo("setIdByName()", "Successfully loaded "<< name << " Refinement criterion code.");
		return true;
	} else {
		logError("setIdByName()", "Requested Refinement criterion code '"<< name << "' not known.");
		return false;
	}
}

void GlobalRefinementCriterion::setByParameters(const mexa::mexafile& parameters) {
	std::string seckey = "refinement";
	std::string namekey = "criterion";
	mexa::mexafile param = parameters(seckey);
	if(!param.contains(namekey)) {
		logError("setByParameters()", "For setting up the Refinement criterion, I need a section " << seckey  << " in the parameters, as well as the key " << seckey << "/" << namekey << " to be set to a valid identifier. Instead, I got these parameters: " << parameters.toString());
		std::abort();
	}
	bool couldSetId = setIdByName(param(namekey).get_string());
	if(!couldSetId) {
		logError("setByParameters()", "Could not setup Refinement criterion code from the given parameters:" << param.toString());
		std::abort();
	} else {
		getCode().readParameters(param);
	}
}


REFINEMENT_CONTROL_RETURN GlobalRefinementCriterion::refinementCriterion(REFINEMENT_CONTROL_SIGNATURE) {
	return getCode().refinementCriterion(REFINEMENT_CONTROL_CALL);
}

/************** PARAMETER READING FOR Refinement CRITERIA   ***********************/


void GeometricSphereRefinement::readParameters(const mexa::mexafile& para) {
	radius = para["radius"].as_double();
	
	static tarch::logging::Log _log("GeometricSphereRefinement");
	logInfo("readParameters()", "Refinement at the shell with radius=" << radius);
}


/************** THE ACTUAL Refinement CRITERIA   ***********************/

REFINEMENT_CONTROL_RETURN NoRefinement::refinementCriterion(REFINEMENT_CONTROL_SIGNATURE) const {
	return exahype::solvers::Solver::RefinementControl::Keep;
}

REFINEMENT_CONTROL_RETURN GeometricSphereRefinement::refinementCriterion(REFINEMENT_CONTROL_SIGNATURE) const {
  // lower left, uppe right radius of cell
  double l = tarch::la::norm2(center - dx/2.0);
  double r = tarch::la::norm2(center + dx/2.0);
  constexpr double radius_shell = 0.5; // limit around this shell
  bool isAdmissible = (l > radius_shell) || (r <= radius_shell);

  if (isAdmissible) {
    return exahype::solvers::Solver::RefinementControl::Keep;
  } else {
    return exahype::solvers::Solver::RefinementControl::Refine;
  }
  //printf("Cell has l=%f,r=%f => isAdmissible=%s\n", l, r, isAdmissible?"true":"false");
}
