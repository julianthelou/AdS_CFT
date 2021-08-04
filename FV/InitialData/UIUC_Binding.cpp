#include "InitialData.h"
#include "ADMBase.h"

using namespace CCZ4_InitialData;

/*
 * Since the UIUC code has no external dependencies, it is copied into the
 * ExaHyPE Astrophysics repository and we can safely enable it here.
 */
#define UIUC_BlackHole_AVAILABLE
#include "libuiuc_bh/UIUCInitialData.h"

#if !defined(UIUC_BlackHole_AVAILABLE)

#include <exception>

void uiuc_not_available() {
	throw std::exception("The UIUC Initial Data code is not available. Please recompile with -DUIUC_BlackHole_AVAILABLE.");
}

ImportedUIUC_BlackHole::ImportedTwoPunctures() { uiuc_not_available(); }
void ImportedUIUC_BlackHole::readParameters(const mexa::mexafile& parameters) { uiuc_not_available(); }
void ImportedUIUC_BlackHole::prepare() { uiuc_not_available(); }
void ImportedUIUC_BlackHole::Interpolate(const double* const x, double t, ADMBase& out) { uiuc_not_available(); }

#else

// No include here as it is... in general... above.

ImportedUIUC_BlackHole::ImportedUIUC_BlackHole() {
	uiuc = new UIUC_BlackHole(); // Pointer for lesser compile dependency
}

void ImportedUIUC_BlackHole::readParameters(const mexa::mexafile& para) {
	//tp->log->info("Reading parameters");
	
	#define PAR(varname, typemethod) uiuc->varname = para[#varname].typemethod();
	#define VEC(varname, typemethod) for(int i=0;i<3;i++) uiuc->varname[i] = para[#varname].vec().typemethod()[i];
	// example typesmethods are: as_double, get_string, get_bool, get_double, etc...
	
	// these are all parameters.
	// We could offer a way of falling back to the standard if not any parameter was provided.
	
	PAR(lapse_type, get_int);
	PAR(shift_type, get_int);
	PAR(BH_gauge_choice, get_int);
	VEC(center_offset, as_double);
	PAR(avoid_puncture, get_bool);
	PAR(par_M, as_double);
	PAR(par_chi, as_double);
	VEC(par_P, as_double);
	PAR(beta_r_interior, as_double)
	PAR(beta_r_x0, as_double);
	PAR(beta_r_w, as_double);
	
	tarch::logging::Log _log("ImportedUIUC_BlackHole");
	logInfo("readParameters(mexa)", uiuc->toString());
}

void ImportedUIUC_BlackHole::Interpolate(const double* const x, double t, ADMBase& out) {
	double N[4];
	double pos[3] = {x[0],x[1],x[2]}; // look for a better way
	
	uiuc->InitialData(pos, out.gamma_lo, out.Kextr_lo, N);
	
	out.alpha = N[0];
	for(int i=0; i<3; i++) out.beta_lo[i] = N[i];
}

#endif /* UIUC_BlackHole_AVAILABLE */
