#include "InitialData.h"
#include "ADMBase.h"
#include "Fortran/PDE.h"

namespace idx = CCZ4Fortran::StateVectorSequentialization;

using namespace CCZ4_InitialData;


inline double my_norm2(const double* const x) {
	double r(0);
	for(int d=0; d<DIMENSIONS; d++) r+= x[d]*x[d];
	return std::sqrt(r);
}

/* GUARD to compile only if the TwoPunctures code is available */
#if !defined(TWOPUNCTURES_AVAILABLE)

#include <exception>

void twopunctures_not_available() {
         std::cout << ("The TwoPunctures Initial Data code is not available. Please recompile with -DTWO_PUNCTURES.") << std::endl;
	exit(-1);
}

ImportedTwoPunctures::ImportedTwoPunctures() { twopunctures_not_available(); }
void ImportedTwoPunctures::readParameters(const mexa::mexafile& parameters) { twopunctures_not_available(); }
void ImportedTwoPunctures::prepare() { twopunctures_not_available(); }
void ImportedTwoPunctures::Interpolate(const double* const x, double t, ADMBase& out) { twopunctures_not_available(); }
void ImportedTwoPunctures::PostProcessStateVector(const double* const x, double t, double* Q) { twopunctures_not_available(); }

#else

// TP Black Holes Initial Data code
#include "libtwopunctures/TwoPunctures.h"

TP::TwoPunctures *cheater = nullptr;

// Tie TwoPunctures output to Tarch Logging
struct ExaHyPE_TwoPunctures_Logger : public TP::logger {
	tarch::logging::Log _log;
	std::string state; ///< a convenient name for a method or so
	ExaHyPE_TwoPunctures_Logger() : _log("CCZ4::ImportedTwoPunctures"), state("") {}

	void log(const std::string& msg) override { logInfo(state.c_str(), msg); }
	void error(const std::string& msg) override { logError(state.c_str(), msg); }
	void info(const std::string& msg) override { logInfo(state.c_str(), msg); }
	void warn(const std::string& msg) override { logWarning(state.c_str(), msg); }
};

ImportedTwoPunctures::ImportedTwoPunctures() {
	tp = new TP::TwoPunctures(); // Pointer for lesser compil dependency
	cheater = tp;
	logger = new ExaHyPE_TwoPunctures_Logger;
	tp->log = logger;
}

void ImportedTwoPunctures::readParameters(const mexa::mexafile& para) {
	logger->state = "readParameters(mexa)";
	tp->log->info("Reading parameters");
	// These neat macros allow a 1:1 mapping of (m)exahype parameters
	// to TwoPunctures paramteres.
	
	#define PAR(varname, typemethod) tp->varname = para[#varname].typemethod();
	#define VEC(varname, typemethod) for(int i=0;i<3;i++) tp->varname[i] = para[#varname].vec().typemethod()[i];
	// example typesmethods are: as_double, get_string, get_bool, get_double, etc...
	
	PAR(par_b             , as_double);
	PAR(target_M_plus     , as_double);
	VEC(par_P_plus        , as_double);
	VEC(par_S_plus        , as_double);
	VEC(center_offset     , as_double);
	
	PAR(target_M_minus    , as_double);
	VEC(par_P_minus       , as_double);
	VEC(par_S_minus       , as_double);
	
	PAR(grid_setup_method , get_string);
	PAR(TP_epsilon        , get_double);
	
	// etc., could list here all parameters in TP_Parameters.h
	
	tp->PrintParameters();
}

void ImportedTwoPunctures::PostProcessStateVector(const double* const x, double t, double* Q) {
	double radius = my_norm2(x);

	double r0 = 10; // Domain radius
	double w = 0.1; // Interface width
	
	//Q[idx::domain] = (radius < r0) ? 1 : 0;   // hard cut
	
	Q[idx::domain] = 1./(1 + std::exp( (radius-r0)/w )); // Logistic curve
}

void ImportedTwoPunctures::prepare() {
	logger->state = "prepare()";
	tp->log->info("Preparation run started");
	tp->Run();
	logger->state = "Interpolate(x,t,ADMBase out)";
}

void ImportedTwoPunctures::Interpolate(const double* const x, double t, ADMBase& out) {
	constexpr int numberOfVariables = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
	double Q[numberOfVariables] = {0.0};

	// TwoPunctures sets us only gij, kij, alpha, but the interface is written
	// in a way that it works on the full CCZ4 state vector. This was from old
	// times when I thought it makes sense to deal with the full state vector
	// even when it is only set partially.
	tp->Interpolate(x, Q);
	
	// Therefore, we copy now a bit idiotic around.
	out.copyFromStateVector(Q);
	//printf("Queried ID for C at x=[%.20e,%.20e,%.20e], t=%e\n", x[0],x[1],x[2],t);
}

/*
// This is for testing Fortran.
extern "C" void twopunctures_interpolate(double* x, double* Q) {
	if(cheater == nullptr) std::abort();
	cheater->Interpolate(x,Q);
	//printf("Queried ID for F at x=[%.20e,%.20e,%.20e]\n", x[0],x[1],x[2]);
}
*/


#endif /* TWOPUNCTURES_AVAILABLE */
