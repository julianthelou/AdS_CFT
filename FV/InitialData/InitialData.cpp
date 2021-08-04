#include "InitialData/InitialData.h"
#include "InitialData/ADMBase.h"

#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"

#include "kernels/KernelUtils.h" // idx classes
#include "kernels/GaussLegendreBasis.h" // gaussLegendreNodes
#include "kernels/aderdg/generic/c/computeGradients.cpph"

#include <algorithm> // transform
#include <cctype> // tolower
#include <exception>

using namespace CCZ4_InitialData;

/******************************************************************************\
|*                                                                            *|
|*                                                                            *|
|*                 Global Initial Data registration                           *|
|*                 ================================                           *|
|*                                                                            *|
|*     If you want your initial data, write an InitialDataCode function       *|
|*     and register it afterwards here.                                       *|
|*                                                                            *|
\******************************************************************************/

InitialDataCode* CCZ4_InitialData::InitialDataCode::getInstanceByName(std::string idname) {
	std::transform(idname.begin(), idname.end(), idname.begin(), ::tolower); // toLower(idname)
	
	if(idname == "gaugewave")       return new GaugeWaveID();
	if(idname == "twopunctures")    return new ImportedTwoPunctures();
	if(idname == "uiuc")            return new ImportedUIUC_BlackHole();
	
	//if(idname == "gaugewave") 	return new FortranPointwiseID(typesDef_ICType.CCZ4GaugeWave);
	//if(idname == "puncture")	return new FortranPointwiseID(typesDef_ICType.CCZ4Puncture);
	//if(idname == "twopunctures") 	return new FortranPointwiseID(typesDef_ICType.CCZ4TwoPunctures);
	//if(idname == "grmhdaccretion") 	return new FortranPointwiseID(typesDef_ICType.CCZ4GRMHDAccretion);
	//if(idname == "kerr2d") 		return new FortranPointwiseID(typesDef_ICType.CCZ4Kerr2D);
	return nullptr;
}

/******************************************************************************\
|*                                                                            *|
|*                                                                            *|
|*     Administrative Initial Data loading stuff                              *|
|*     ========================================                               *|
|*                                                                            *|
|*                                                                            *|
\******************************************************************************/


GlobalInitialData& CCZ4_InitialData::GlobalInitialData::getInstance() {
	static GlobalInitialData* me = nullptr;
	if(!me) me = new GlobalInitialData(nullptr, "null");
	return *me;
}

bool CCZ4_InitialData::GlobalInitialData::setIdByName(const std::string& _name) {
	static tarch::logging::Log _log("GlobalInitialData");
	name = _name;
	
	if(alreadySetByName()) {
		// This must not be an error but is good practise: Having a DG and FV solver which both
		// try to initialize the global initial data just to be safe.
		logWarning("setIdByName(str)", "Have already set global initial data, will not overwrite it.");
	}
	
	id = CCZ4_InitialData::InitialDataCode::getInstanceByName(name);
	if(id) {
		logInfo("setIdByName(str)", "Successfully loaded "<< name << " Initial Data.");
		return true;
	} else {
		logError("setIdByName(str)", "Requested Initial Data '"<< name << "' not known.");
		return false;
	}
}

GlobalInitialData& CCZ4_InitialData::GlobalInitialData::setByParameters(const mexa::mexafile& parameters) {
	if(alreadySetParameters) {
		// This must not be an error but is good practise: Having a DG and FV solver which both
		// try to initialize the global initial data just to be safe.
		logWarning("setByParameters(mexa)", "Have already set global initial data by parameters, will not override it.");
		return *this;
	}
	
	std::string idseckey = "initialdata";
	std::string idnamekey = "name";
	mexa::mexafile idparam = parameters(idseckey);
	if(!idparam.contains(idnamekey)) {
		logError("setByParameters(mexa)", "For setting up the initila data, I need a section " << idseckey  << " in the parameters, as well as the key " << idseckey << "/" << idnamekey << " to be set to a valid initial data function. Instead, I got these parameters: " << parameters.toString());
		std::abort();
	}
	bool couldSetId = setIdByName(idparam(idnamekey).get_string());
	if(!couldSetId) {
		logError("setByParameters(mexa)", "Could not create Initial Data. Cannot solve an initial value problem without initial data.");
		std::abort();
	} else {
		id->readParameters(idparam);
	}
	alreadySetParameters = true;
	return *this;
}

GlobalInitialData& CCZ4_InitialData::GlobalInitialData::prepare() {
	if(!alreadyPrepared) {
		id->prepare();
	}
	alreadyPrepared = true;
	return *this;
}

bool CCZ4_InitialData::GlobalInitialData::ensureGlobalIDareSet(bool doCrash) {
	bool isSet = (getInstance().id != nullptr);
	if(doCrash && !isSet) {
		static tarch::logging::Log _log("GlobalInitialData");
		logError("InitialData()", "Cannot access InitialData because no initial Data has been defined yet.");
		std::abort();
	}
	return isSet;
}


// DEBUGGING 
#include "CCZ4Solver_ADERDG_Variables.h"
// or:   "CCZ4Solver_FV_Variables.h", as you prefer. Then adopt the aliases below.
namespace idx = CCZ4::CCZ4Solver_ADERDG_Variables::shortcuts;

void compare(double*Q, double*F, int i) {
	if(std::abs(Q[i] - F[i]) > 1e-12) {
		printf("C[%d]=%e, F[%d]=%e\n",i,Q[i],i,F[i]);
	}
}

/* Debugging: Fortran for comparison: Declaration */
// extern "C" void initialfield_(double* Q, const double* const x, const double* const t);

void CCZ4_InitialData::InitialData(const double* x, double t, double* Q) {
	GlobalInitialData::ensureGlobalIDareSet();
	GlobalInitialData::getInitialDataCode().adjustPointSolution(x,t,Q);
	
	// Debugging
	/*
	// Compare the data with Fortran! :-)
	constexpr int nVars = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
	double F[nVars];
	
	initialfield_(F, x, &t);
	// Make this test:
	// GlobalInitialData::getInitialDataCode().adjustPointSolution(x,t,F); // identical output!
	
	// Debug only conformal factor:
	for(int i=0; i<nVars; i++) {
		//if(i== idx::phi)
		//	compare(Q,F,i);
		double diff = std::abs(Q[i] - F[i]);
		if(diff > 1e-10) {
			printf("at (x;t)=(%f,%f,%f;%f), ", x[0],x[1],x[2],t);
			printf("C[%d]=%e; F[%d]=%e; diff[%d]=%e\n",i,Q[i],i,F[i],i,diff);
		}
	}
	
	//std::abort();
	*/
}

/******************************************************************************\
|*                                                                            *|
|*                                                                            *|
|*     Pointwise and Patchwise Initial Data Reader for FO-CCZ4                *|
|*     =======================================================                *|
|*                                                                            *|
|*     The following functions are the main workhorses, they do all           *|
|*     the gradient computation. However, all PDE-related work                *|
|*     (conversion from and to ADM, BSSNOK, SO/FO-CCZ4) is delegated          *|
|*     to functions in ADMBase.h/cpp.                                         *|
|*                                                                            *|
\******************************************************************************/

void CCZ4_InitialData::InitialDataCode::adjustPointSolution(const double* const x,const double t,double* Q) {
	constexpr int nVars = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
	
	// Zeroth, for safety initialize everything to zero.
	for(int i=0; i<nVars; i++) Q[i] = 0;
	
	setExaHyPEspecificQuantities(x, Q);
	
	// First, set the incomplete ID at the point
	ADMBase adm;
	Interpolate(x, t, adm);
	adm.importSOCCZ4(Q);

	// Second, compute gradients around point
	double Qx[nVars], Qy[nVars], Qz[nVars];
	double *gradQ[DIMENSIONS] = { Qx,Qy,Qz };
	FD_Interpolate(x, t, gradQ);
	
	// Third, compute FO-CCZ4 quantities based on the gradients
	adjustAuxilliaries(Q, gradQ);
	
	// Fourth, let the user have a final look
	PostProcessStateVector(x,t,Q);
}

void CCZ4_InitialData::InitialDataCode::adjustPatchSolution(
  const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
  const tarch::la::Vector<DIMENSIONS, double>& dx,
  const double t, const double dt, double* luh) {
	assert(DIMENSIONS == 3);
	
	// Attention: This function is so far untested because I have not looked into the
	// adjustPatchSolution for a long time.
	
	using namespace kernels; // idx4, idx5
	using namespace tarch::la;
	typedef Vector<DIMENSIONS,double> dvec;
	constexpr int numberOfVariables = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
	constexpr int order = CCZ4::AbstractCCZ4Solver_ADERDG::Order;
	constexpr int basisSize = order + 1;
	constexpr int basisSize3 = basisSize*basisSize*basisSize;

	// First, set the incomplete ID on all DOF
	idx4 idx_luh(basisSize, basisSize, basisSize, numberOfVariables);
	for(int iz = 0; iz < basisSize; iz++) {
		dvec q;
		q[2] = kernels::legendre::nodes[order][iz];
		for(int iy = 0; iy < basisSize; iy++) {
			q[1] = kernels::legendre::nodes[order][iy];
			for(int ix = 0; ix < basisSize; ix++) {	
				q[0] = kernels::legendre::nodes[order][ix];
				
				dvec r  = cellCentre + dx * ( q - 0.5 );
				double *Q  = &luh[idx_luh(iz,iy,ix,0)];
				for(int i=0; i<numberOfVariables; i++) Q[i] = 0;
				
				setExaHyPEspecificQuantities(r.data(), Q);
				
				ADMBase adm;
				Interpolate(r.data(), t, adm);
				adm.importSOCCZ4(Q);
			} // x
		} // y
	} // z
	
	// Second, for the time being, just compute *all* gradients.
	// No more call to the ID.
	double gradQ[basisSize3 * DIMENSIONS * numberOfVariables];
        kernels::aderdg::generic::c::computeGradQ<CCZ4::AbstractCCZ4Solver_ADERDG>(gradQ, luh, dx);
	
	// Third, compute the FO-CCZ4 quantities based on the gradients
	idx5 idx_gradQ(basisSize, basisSize, basisSize, DIMENSIONS, numberOfVariables);
	for(int iz = 0; iz < basisSize; iz++) {
		for(int iy = 0; iy < basisSize; iy++) {
			for(int ix = 0; ix < basisSize; ix++) {
				// gradients in x, y, z direction
				double *local_gradQ[DIMENSIONS];
				for(int d=0; d < 3; d++) local_gradQ[d] = gradQ + idx_gradQ(iz,iy,ix,d,0);
				double *Q  = luh + idx_luh(iz,iy,ix,0);
				adjustAuxilliaries(Q, local_gradQ);
			} // x
		} // y
	} // z
	
	// Fourth, let the user have a final look
	for(int iz = 0; iz < basisSize; iz++) {
		dvec q;
		q[2] = kernels::legendre::nodes[order][iz];
		for(int iy = 0; iy < basisSize; iy++) {
			q[1] = kernels::legendre::nodes[order][iy];
			for(int ix = 0; ix < basisSize; ix++) {	
				q[0] = kernels::legendre::nodes[order][ix];
				dvec r  = cellCentre + dx * ( q - 0.5 );
				double *Q  = &luh[idx_luh(iz,iy,ix,0)];
				
				PostProcessStateVector(r.data(), t, Q);
			} // x
		} // y
	} // z
	
}


void CCZ4_InitialData::InitialDataCode::FD_Interpolate(const double* const x, const double t, double **gradQ) {
	constexpr int nVars = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
	using namespace tarch::la;
	typedef Vector<DIMENSIONS,double> dvec;
	const dvec xc(x[0],x[1],x[2]); // ugly workaround for missing tarch::la::vector(const double* const) constructor.
	
	constexpr double epsilon = 1e-7, eps4 = 1e-4;
	ADMBase Ap1,Am1,Ap2,Am2;
	double Qp1[nVars],Qm1[nVars],Qp2[nVars],Qm2[nVars];
	
	// Metric derivative computed with a fourth order central finite difference 
	// as done by Michael D
	for(int d=0; d<DIMENSIONS; d++) {
		dvec xp1(xc), xm1(xc), xp2(xc), xm2(xc);
		xp1(d) += eps4;
		xm1(d) -= eps4;
		xp2(d) += 2*eps4;
		xm2(d) -= 2*eps4;
		
		Interpolate(xp1.data(), t, Ap1); Ap1.importSOCCZ4(Qp1); 
		Interpolate(xm1.data(), t, Am1); Am1.importSOCCZ4(Qm1);
		Interpolate(xp2.data(), t, Ap2); Ap2.importSOCCZ4(Qp2);
		Interpolate(xm2.data(), t, Am2); Am2.importSOCCZ4(Qm2);
		
		for(int i=0; i<nVars; i++) {
			gradQ[d][i] = ( 8.0*Qp1[i] - 8.0*Qm1[i]  + Qm2[i]   - Qp2[i]  )/(12.0*eps4);
		}
	}
}

void CCZ4_InitialData::InitialDataCode::setExaHyPEspecificQuantities(const double* const x, double* Q) {
	Q[idx::domain] = 1.0; // inside

	// Store coordinates in state vector.
	for(int i=0; i<DIMENSIONS; i++)
		Q[idx::pos + i] = x[i];
}

/*
double CCZ4_InitialData::getDIMradius(const mexa::mexafile& params, double* Q) {
	tarch::logging::Log _log("CCZ4_InitialData");

	if(params.contains("dim-radius")) {
	} else {
		logInfo("getDIMradius(mexa)", "No Radius given.");
	}
}
*/
