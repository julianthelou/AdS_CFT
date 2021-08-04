#include "ToyID.h"
#include "ADMBase.h"

#include <cmath>
using namespace std;
constexpr double pi = M_PI;

constexpr double ICA = 0.1; ///< Amplitude of the wave

void GaugeWaveID::Interpolate(const double* const x, double t, CCZ4_InitialData::ADMBase& adm) {
	double w = 2.0*pi*(x[0]-t);

	// mind the negative sign in contrast to eq 4.4 in https://arxiv.org/pdf/gr-qc/0305023.pdf
	// here, we follow Dumbsers InitialField.f90 GaugeWave.
	double H = 1.0 - ICA*sin(w);
	double Kxx = -pi*ICA*cos(w)/sqrt(1-ICA*sin(w));

	adm.one();
	adm.alpha = sqrt(H);
	adm.gamma_lo[0][0] = H;
	adm.Kextr_lo[0][0] = Kxx;
	
	// to be sure, because adm.one() is probably wrong.
	adm.Kextr_lo[1][1] = 0;
	adm.Kextr_lo[2][2] = 0;
	
	//printf("Queried GW for C at x=[%.20e,%.20e,%.20e], t=%e,\n", x[0],x[1],x[2],t);
}

void GaugeWaveID::adjustPointSolution(const double* const x,const double t,double* Q) {
  constexpr int nVars = CCZ4::AbstractCCZ4Solver_ADERDG::NumberOfVariables;
  double HH     = 1.0 - ICA*sin( 2.0*pi*( x[0] - t));
  double dxH    = -2.0*pi*ICA*cos( 2.0 * pi*(x[0] - t));
  double dxphi  = - pow(HH,(-7.0/6.0))*dxH/6.0;
  double phi    = pow(( 1.0 / HH),(1.0/6.0));
  double Kxx    = - pi*ICA*cos( 2.0 * pi*(x[0] - t))/sqrt( 1.0 - ICA*sin( 2.0*pi*( x[0] - t))  );
  double traceK = Kxx/HH;
  
  memset(Q, 0, nVars*sizeof(double));
  Q[0]  = phi*phi*HH                         ;
  Q[3]  = phi*phi                            ;
  Q[5]  = phi*phi                            ;
  Q[6]  = phi*phi*(Kxx - 1.0/3.0*traceK*HH ) ;
  Q[9] =  phi*phi*(0.0 - 1.0/3.0*traceK*1.0) ;
  Q[11] = phi*phi*(0.0 - 1.0/3.0*traceK*1.0)  ;
  Q[16] = log(sqrt(HH));
  Q[13] = 2.0/(3.0*pow(HH,(5.0/3.0)))*dxH        ;
  Q[23] = 1.0/(2.0*HH)*dxH               ;
  Q[35] = pow(HH,(-1.0/3.0))*dxH/3.0         ;
  Q[38] = phi*dxphi                     ;
  Q[40] = phi*dxphi                    ;
  Q[53] = traceK;
  Q[54] = log(phi);
  Q[55] = dxphi/phi;
}
// A debugging Fortran callout function
/*
extern "C" void debugme_interpolate(double* x, double* Q) {
	// versatile debugging point
	if(std::getenv("ABORT_F")) std::abort();
	printf("Queried GW for F at x=[%.20e,%.20e,%.20e],\n", x[0],x[1],x[2]);
}
*/
