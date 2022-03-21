#include "FOCCZ4Solver.h"

#include "FOCCZ4Solver_Variables.h"

#include "kernels/finitevolumes/musclhancock/c/musclhancock.h"
#include "InitialData.h"
#include "PDE.h"

#include "AdS.h"

tarch::logging::Log FOCCZ4::FOCCZ4Solver::_log( "FOCCZ4::FOCCZ4Solver" );

void FOCCZ4::FOCCZ4Solver::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver.h".
    std::string reference = constants.getValueAsString("reference");
    const int length=reference.length();
    logInfo("init(...)","Reference setup:"<<reference);

    initparameters_(&length,&reference[0]);
}

void FOCCZ4::FOCCZ4Solver::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
 const int nVars = FOCCZ4::FOCCZ4Solver::NumberOfVariables;
 
  if (tarch::la::equals(t,0.0)) {
    int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
    double cms = exahype::solvers::Solver::getCoarsestMeshSize();
    const int order = FOCCZ4::FOCCZ4Solver::PatchSize;
    std::fill_n(Q,96,0.0);

    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);

    initialdata_(x_3, &t, Q);
  }
  else
  {
	  enforceccz4constraints_(Q); 
  }

}

void FOCCZ4::FOCCZ4Solver::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
 
  double nv[3] = {0.};
  nv[dIndex] = 1;
  // pdeeigenvalues_(lambda, Q, nv);

  for(int i = 0; i< 59 ; i++){
   	lambda[i]=1;
   // assert(std::isfinite(lambda[i]));
  }
 
 
}

void FOCCZ4::FOCCZ4Solver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* const stateOutside) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

    const int nVars = FOCCZ4::FOCCZ4Solver::NumberOfVariables;

	double Qgp[nVars];

	double ti = t + 0.5 * dt;
	// Compute the outer state according to the initial condition
  double x_3[3];
  x_3[2]=0;
  std::copy_n(&x[0],DIMENSIONS,&x_3[0]);

  if (isAdS3_){
  // r = sqrt(x**2 + y**2)
  //double r = sqrt(x[0]**2 + x[1]**2)
  // rho = r/(1+r)
  double rho = AdSrho(x[0]);
  // q=1-rho
  double q = AdSq(rho);
  // a=q^2 + rho^2
  double a = AdSa(AdSq(rho), AdSrho(x[0]));
  double Da = DAdSa(AdSrho(x[0]));
  
  // Setting all statevariables wich are zero to -stateInside
  for(int m=0; m < nVars; m++){
  	stateOutside[m]=-stateInside[m];
  }
  
  // Overwriting the components which are not zero
  // conformally decomposed spatial metric
  // at the boundary q=0, rho=1, a=1 and Da=2 
   	stateOutside[0] = 1;
  	stateOutside[3] = 1;
  	
  //	stateOutside[0] = 1/(sqrt(a)*rho);
  //	stateOutside[3] = sqrt(a)*rho;
  	
  // conformally traceless-part of the extrinsic curvature

  // Theta
  	
  // Gamma hat
	stateOutside[13] = - 2;
  	//stateOutside[13] = - 0.5*rho*Da/sqrt(a) - sqrt(a);
  	
  // lapse function
    	stateOutside[16] = 0;
    	
  //	stateOutside[16] = log(1-q) + 0.5*log(1+q*q/(rho*rho)) - log(q);
  	
  // shift vector
  	
  // Gamma driver condition
  	
  // Auxiliary variable A_i
    	stateOutside[23] = 1/0.01;   	
  //	stateOutside[23] = rho/(a*q);
  	
  // Auxiliary variable B_i^j
  	
  // Auxiliary variable D_ijk
    	stateOutside[35] = -1;
  	stateOutside[38] = 1;
  	
  //	stateOutside[35] = -1/(2*pow(a,3/2)) + q/(2*rho*pow(a,3/2)) - 1/(2*sqrt(a)*rho*rho);
  //	stateOutside[38] = rho*rho/(2*sqrt(a)) - q*rho/(2*sqrt(a)) + 0.5*sqrt(a);
  	
  // Trace of the extrinsic curvature
  	
  // conformal factor
    	stateOutside[54] = - 100;
    	
  //	stateOutside[54] = log(q) + 0.25*log(a) - 0.5*log(rho);
  	
  // Auxiliary variable P_i
    	stateOutside[55] = -1/0.01;

  //	stateOutside[55] = -1/q + rho/(2*a) - q/(2*a) - 1/(2*rho);
  	
  // Trace of the extrinsic curvature at t=0
  }
  else {
    initialdata_(x_3, &ti, Qgp);
    for(int m=0; m < nVars; m++) {
          stateOutside[m] = Qgp[m];
    }
  }
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void FOCCZ4::FOCCZ4Solver::flux(const double* const Q,double** const F) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  const int nVar = FOCCZ4::FOCCZ4Solver::NumberOfVariables;
  if(DIMENSIONS == 2){
    double F_3[nVar];
    pdeflux_(F[0], F[1],F_3, Q);
  }else{
    pdeflux_(F[0], F[1],F[2], Q);
  }
  
}




//You can either implement this method or modify fusedSource
void FOCCZ4::FOCCZ4Solver::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  pdesource_(S, Q);

}

void  FOCCZ4::FOCCZ4Solver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  // Tip: You find documentation for this method in header file "FOCCZ4::FOCCZ4Solver.h".
  // Tip: See header file "FOCCZ4::AbstractFOCCZ4Solver.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
 pdencp_(BgradQ, Q, gradQ);
 
}

double FOCCZ4::FOCCZ4Solver::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, const double* gradQL, const double* gradQR, const double* cellSize, int direction) {


    //return kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, FOCCZ4Solver>(*static_cast<FOCCZ4Solver*>(this), fL,fR,qL,qR,gradQL, gradQR, cellSize, direction);
	constexpr int numberOfVariables = AbstractFOCCZ4Solver::NumberOfVariables;

	//printf("SONO QUI IN riemannSolver");
	/* HLLEM */
	
	//const int numberOfVariables = GRMHDb::AbstractGRMHDbSolver::NumberOfVariables;
	const int numberOfParameters = FOCCZ4::AbstractFOCCZ4Solver::NumberOfParameters;
	const int numberOfData = numberOfVariables + numberOfParameters;
	const int order = 0;  // for finite volume we use one single d.o.f., i.e. the cell average.
	const int basisSize = order + 1;
	// Compute the average variables and parameters from the left and the right
	double QavL[numberOfData] = { 0.0 }; // ~(numberOfVariables+numberOfParameters)
	double QavR[numberOfData] = { 0.0 }; // ~(numberOfVariables+numberOfParameters)
#ifdef CCZ4GRHD 	
    double lambda;
     hllemfluxfv_(&lambda, fL, fR, qL, qR, &direction);
#else
    double lambda = kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, FOCCZ4Solver>(*static_cast<FOCCZ4Solver*>(this), fL, fR, qL, qR, gradQL, gradQR, cellSize, direction);
#endif
	
	//std::cout << lambda << std::endl;
	//double1 lambda = 10.0;
	return lambda; 
}

