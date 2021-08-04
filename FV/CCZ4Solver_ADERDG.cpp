#include "CCZ4Solver_ADERDG.h"

#include "CCZ4Solver_ADERDG_Variables.h"

#include "kernels/GaussLegendreBasis.h"
#include "kernels/KernelUtils.h" // matrix indexing
#include <cstring> // memset
#include "Fortran/PDE.h"
//#include "RadiativeBoundaryConditions.h"
#include "InitialData/InitialData.h"
#include "GridCriteria/LimiterCriteria.h"
#include "GridCriteria/RefinementCriteria.h"
#include "Parameters/mexa.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/Lock.h"


// #include <fenv.h>

// This doesn't seem to be needed anymore
//#include "exahype/disableOptimization.h" // we experience compiler bugs sometimes.

const int nDim = DIMENSIONS;
const int nVar = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
const int order = CCZ4::AbstractCCZ4Solver_ADERDG::Order;
const int basisSize = order+1;

tarch::logging::Log CCZ4::CCZ4Solver_ADERDG::_log( "CCZ4::CCZ4Solver_ADERDG" );

using namespace CCZ4_InitialData;

namespace PDE = CCZ4Fortran::PDE; // Solving CCZ4 only
//namespace PDE = CCZ4_RadBC::DIM_CCZ4  // Solving the coupled CCZ4-BohrSommerfeld system


void CCZ4::CCZ4Solver_ADERDG::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
	mexa::mexafile mf = mexa::fromSpecfile(constants.getAllAsOrderedMap(), constants.toString());
	
	// Debugging
	logInfo("init(...)", "Mexa configuration: \n" << mf.toString());
	logInfo("init(...)", "ID configuration: \n" << mf("initialdata").toString());
	logInfo("init(...)", "ID NAME: '" << mf.get("initialdata/name").get_string() << "'\n");
	
	static tarch::multicore::BooleanSemaphore initializationSemaphoreDG;
	tarch::multicore::Lock lock(initializationSemaphoreDG);
	// PDE Runtime Parameters
	CCZ4Fortran::GlobalPDEParameters::getInstance().setByParameters(mf);
	std::cout << "Set Parameters" << std::endl;
	
	// Initial Data
	GlobalInitialData::getInstance().setByParameters(mf).prepare();
	lock.free();
	// Grid Criteria
	//GlobalLimitingCriterion::getInstance().setByParameters(mf);
	//GlobalRefinementCriterion::getInstance().setByParameters(mf);
	
	std::cout << "Finished init" << std::endl;
}


bool CCZ4::CCZ4Solver_ADERDG::isPhysicallyAdmissible(
      const double* const solution,
      const double* const observablesMin,const double* const observablesMax,
      const bool wasTroubledInPreviousTimeStep,
      const tarch::la::Vector<DIMENSIONS,double>& center,
      const tarch::la::Vector<DIMENSIONS,double>& dx,
      const double t) const {

      if(center[0] < 6.0 && center[0] > -6.0 && fabs(center[1]) < 2.5 && fabs(center[2]) < 2.5)
      //if(sqrt(center[0]*center[0]+center[1]*center[1]+center[2]*center[2]) < 0.7)
      	return false;
      return true;//GlobalLimitingCriterion::isPhysicallyAdmissible(solution, observablesMin, observablesMax, wasTroubledInPreviousTimeStep, center, dx, t);
}

void CCZ4::CCZ4Solver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* const Q) {
	CCZ4Fortran::AdjustPointSolution(x,t,dt,Q); //sets initial conditions and enforces ADM constraints
	//if(tarch::la::equals(t,0.0))
	//	GlobalInitialData::getInstance().setByParameters(mf).free();
}

void CCZ4::CCZ4Solver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,const double* const fluxIn,const double* const stateIn,const double* const gradStateIn,double* const fluxOut,double* const stateOut) {
	double Qgp[nVar];

	// Impose exact boundary conditions
	std::memset(stateOut, 0, nVar * sizeof(double));
	std::memset(fluxOut, 0, nVar * sizeof(double));

	//double F[nDim][nVar];
	//kernels::idx2 F_idx(nDim, nVar);

	const int basisSize = order+1;

	// Integrate solution in gauss points (Qgp) in time
	for(int i=0; i < basisSize; i++)  { // i == time
	const double weight = kernels::legendre::weights[order][i];
	const double xi = kernels::legendre::nodes[order][i];
	double ti = t + xi * dt;

	InitialData(x,ti,Qgp);
	for(int m=0; m < nVar; m++) {
		stateOut[m] += weight * Qgp[m];
	}
	}
	//for(int i=0; i< nVar; i++)
	//	stateOut[i] = stateIn[i]; 

}


exahype::solvers::Solver::RefinementControl CCZ4::CCZ4Solver_ADERDG::refinementCriterion(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
      if(center[0] < 6.5 && center[0] > -6.5 && fabs(center[1]) < 3.5 && fabs(center[2]) < 3.5)
      	return exahype::solvers::Solver::RefinementControl::Refine; //GlobalRefinementCriterion::refinementCriterion(luh, center, dx, t, level);
      return exahype::solvers::Solver::RefinementControl::Keep; //GlobalRefinementCriterion::refinementCriterion(luh, center, dx, t, level);
}


void CCZ4::CCZ4Solver_ADERDG::eigenvalues(const double* const Q,const int d,double* const lambda) {
	PDE::eigenvalues(Q, d, lambda);
}

void CCZ4::CCZ4Solver_ADERDG::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
	PDE::algebraicSource(Q, S);
}

void  CCZ4::CCZ4Solver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
	PDE::nonConservativeProduct(Q, gradQ, BgradQ);
}

void CCZ4::CCZ4Solver_ADERDG::fusedSource(const double* const restrict Q, const double* const restrict gradQ, double* const restrict S){
	PDE::fusedSource(Q, gradQ, S);
}

//void CCZ4::CCZ4Solver_ADERDG::fusedSource_vect(double* * restrict Q,  double*  *  *  restrict gradQ, double*  *  restrict S, const int size) {
  //printf("size = %d\n", size);
  //PDE::fusedSourceVector(Q[0], gradQ[0][0], gradQ[1][0], gradQ[2][0], S[0]);
  //printf("Calling S\n");
  // for (int i=0; i<size; i++) {
  // for (int j=0; j<CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables; j++) {
    // printf("S[%d] = %.16f\n", i*4+j, S[j][i]);
  // }
  // }
/*
   for (int i=0; i<size; i++) {
    double QScal[CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables];
    double gradQScal[CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables*3];
    double SScal[CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables];
    for(int j=0; j<CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables; j++) {
      QScal[j] = Q[j][i];
      gradQScal[j] = gradQ[0][j][i];
      gradQScal[j+CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables] = gradQ[1][j][i];
      gradQScal[j+CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables*2] = gradQ[2][j][i];
    }
    fusedSource(QScal, gradQScal, SScal);
    for(int j=0; j<CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables; j++) {
      S[j][i] = SScal[j];
    }
  } 
  */
//}

//void  CCZ4::CCZ4Solver_ADERDG::nonConservativeProduct_vect(const double* const * const Q, const double* const * const * const gradQ, double** const BgradQ, const int size) {
   /*
   for (int i=0; i<size; i++) {
    double QScal[CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables];
    double gradQScal[CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables*3];
    double ncpScal[CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables];
    for(int j=0; j<CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables; j++) {
      QScal[j] = Q[j][i];
      gradQScal[j] = gradQ[0][j][i];
      gradQScal[j+CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables] = gradQ[1][j][i];
      gradQScal[j+CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables*2] = gradQ[2][j][i];
    }
    nonConservativeProduct(QScal, gradQScal, ncpScal);
    for(int j=0; j<CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables; j++) {
      BgradQ[j][i] = ncpScal[j];
    }
  } 
  */
  //PDE::nonConservativeProductVector(Q[0], gradQ[0][0], gradQ[1][0], gradQ[2][0], BgradQ[0]);
//}
