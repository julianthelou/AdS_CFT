#ifndef __EXAHYPE_CCZ4_TOY_ID_CODE_BINDING__
#define __EXAHYPE_CCZ4_TOY_ID_CODE_BINDING__

#include "InitialData.h"

struct GaugeWaveID : public CCZ4_InitialData::InitialDataCode {
  
  void Interpolate(const double* const x, double t, CCZ4_InitialData::ADMBase& out) override;
  
  void adjustPointSolution(const double* const x,const double t,double* Q) override;
};

#endif /* __EXAHYPE_CCZ4_TOY_ID_CODE_BINDING__ */
