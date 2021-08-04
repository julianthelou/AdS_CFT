#ifndef __EXAHYPE_CCZ4_ID_TWOPUNCTURES_BINDING__
#define __EXAHYPE_CCZ4_ID_TWOPUNCTURES_BINDING__

// C interface to TwoPunctures:
namespace TP { class TwoPunctures; } // Forward declaration
struct ExaHyPE_TwoPunctures_Logger;  // another forward declaration

#include "InitialData.h"

class ImportedTwoPunctures : public CCZ4_InitialData::InitialDataCode {
public:
	TP::TwoPunctures* tp;
	ExaHyPE_TwoPunctures_Logger* logger;
	
	ImportedTwoPunctures();

	/// Accept runtime parameters
	void readParameters(const mexa::mexafile& parameters) override;

	/// Find a consistent solution of Einsteins Equations on a spectral domain
	void prepare() override;

	void Interpolate(const double* const x, double t, CCZ4_InitialData::ADMBase& out) override;
	
	void PostProcessStateVector(const double* const x, double t, double* Q) override;
};

#endif /* __EXAHYPE_CCZ4_ID_TWOPUNCTURES_BINDING__ */
