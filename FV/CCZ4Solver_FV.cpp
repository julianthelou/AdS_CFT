#include "CCZ4Solver_FV.h"

#include "CCZ4Solver_FV_Variables.h"

#include <cstdlib> // exit()
#include "Fortran/PDE.h"
#include "InitialData/InitialData.h"
#include "Parameters/mexa.h"


tarch::logging::Log CCZ4::CCZ4Solver_FV::_log( "CCZ4::CCZ4Solver_FV" );

// enable nan tracker
//#include <fenv.h>

//using namespace CCZ4Fortran;
using namespace CCZ4_InitialData;

namespace PDE = CCZ4Fortran::PDE; // Solving CCZ4 only
//namespace PDE = CCZ4_RadBC::DIM_CCZ4  // Solving the coupled CCZ4-BohrSommerfeld system


void CCZ4::CCZ4Solver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  //mexa::mexafile mf = mexa::fromSpecfile(constants.getAllAsOrderedMap(), constants.toString());
  
  //static tarch::multicore::BooleanSemaphore initializationSemaphoreDG;
  //tarch::multicore::Lock lock(initializationSemaphoreDG);
	 
  // PDE Runtime Parameters
  //CCZ4Fortran::GlobalPDEParameters::getInstance().setByParameters(mf);
  
  // Initial Data
  //GlobalInitialData::getInstance().setByParameters(mf).prepare();
  //lock.free();
  // Boundary conditions for CCZ4: Currently always exact.
  //GlobalBoundaryConditions::getInstance().initializeDG(this).readParameters(mf("boundaries"));
}

void CCZ4::CCZ4Solver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  CCZ4Fortran::AdjustPointSolution(x,t,dt,Q);
  //if(tarch::la::equals(t,0.0))
  //	  GlobalInitialData::getInstance().setByParameters(mf).free();
}

void CCZ4::CCZ4Solver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  PDE::eigenvalues(Q, dIndex, lambda);
}

void CCZ4::CCZ4Solver_FV::boundaryValues(const double* const x,  const double t,const double dt, const int faceIndex, const int d, const double* const stateInside, double* const stateOutside) {
  // exact BC without integration
  InitialData(x, t, stateOutside);
  //for(int i=0; i< nVar; i++)
  //	stateOutside[i] = stateInside[i]; 
}

// if you want to use fluxes...
/*
void CCZ4::CCZ4Solver_FV::flux(const double* const Q,double** const F) {
  pdeflux_(F[0], F[1], DIMENSIONS==3 ? F[2] : nullptr, Q);
}
*/

//void CCZ4::CCZ4Solver_FV::fusedSource(const double* const restrict Q, const double* const restrict gradQ, double* const restrict S){
//	PDE::fusedSource(Q, gradQ, S);
//}
void CCZ4::CCZ4Solver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  PDE::algebraicSource(Q, S);
}

void  CCZ4::CCZ4Solver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  PDE::nonConservativeProduct(Q, gradQ, BgradQ);
}
