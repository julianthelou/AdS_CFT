#ifndef __EXAHYPE_CCZ4_UIUC_BLACH_HOLES_ID_CODE_BINDING__
#define __EXAHYPE_CCZ4_UIUC_BLACH_HOLES_ID_CODE_BINDING__

// UIUCInitialData.h, the external code.
struct UIUC_BlackHole; // Forward declaration

#include "InitialData.h"

class ImportedUIUC_BlackHole : public CCZ4_InitialData::InitialDataCode {
public:
	UIUC_BlackHole* uiuc;
	
	ImportedUIUC_BlackHole();

	/// Accept runtime parameters
	void readParameters(const mexa::mexafile& parameters) override;

	// This is an exact code, no preparation neccessary.
	// void prepare() override;

	void Interpolate(const double* const x, double t, CCZ4_InitialData::ADMBase& out) override;
};

#endif /* __EXAHYPE_CCZ4_UIUC_BLACH_HOLES_ID_CODE_BINDING__ */
