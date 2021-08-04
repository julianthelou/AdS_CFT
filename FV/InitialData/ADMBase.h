#ifndef __EXAHYPE_CCZ4_ADMBASE__
#define __EXAHYPE_CCZ4_ADMBASE__

/*
 * ADMBase.h contains functions to import ADM quantities to the FO-CCZ4 language.
 * This is a pure C(++) variant and replaces the Fortran stuff which does similar
 * work. Actually, there are other variants of similar code available, for instance
 * within the Antelope code.
 */

namespace CCZ4_InitialData {

struct InitialDataCode; // Forward declaration from InitialData.h

/**
 * This structure holds the Einstein field information at a given space-time point
 * in the 3+1 split with alpha, beta, gamma and extrinsic curvature.
 *
 * For convention, we store the full matrices here without taking care of their
 * symmetries. Tensor libraries such as Tensish offer options to switch between
 * those.
 **/
struct ADMBase {
	double alpha;
	double beta_lo[3];
	double gamma_lo[3][3];
	double Kextr_lo[3][3];
	
	/// Set all quantities to 0
	void zero();
	
	/// Set all quantities to flat space (Minkowski)
	void one();
	
	/**
	 * Copy the ADM quantities to their respective positions in a
	 * FO-CCZ4 state vector.
	 * This will do no conversion, just copies
	 * the values. That means the outcome is not a proper state vector
	 * and adoptions such as conformal factors, trace freedom, have
	 * to be applied afterwards. The same is true with lapse.
	 **/
	void copyToStateVector(double* Q) const;
	
	/**
	 * In contrast to toStateVector, this will do the full conversion, i.e.
 	 * it computes
	 *   - log(lapse)
	 *   - the conformal metric
	 *   - the traceless extrinsic curvature
	 *   - and all directly accessible SO-CCZ4 quantities: traceK, logW
	 * This is fully local, i.e. computes no derivatives. The state vector
	 * is then just filled with zeros at that place.
	 * 
	 * See the helper functions below for the subsequent steps.
	 **/
	void importSOCCZ4(double *Q) const; // todo refactoring: should be named exportSOCCZ4
 
	/**
	 * Reads the ADM quantities from a FO-CCZ4 state vector. As in the
	 * copyToStateVector function, this will do no conversion and assumes
	 * the state vector to hold the real metric (not the conformal one),
	 * the real Kextr (not a traceless one) and the real lapse (not
	 * exp(lapse)).
	 **/
	void copyFromStateVector(const double* const Q);
	
	// This is a slim ADM Quantity recoverer
	void importStateVector(const double* const Q);
};

/**
 * Fixes the FO-CCZ4 Auxilliaries in Q based on the gradients of the state vector
 * given in gradQ.
 **/
void adjustAuxilliaries(double *Q, const double* const gradQ[3]);


} // namespace CCZ4_InitialData


#endif /* __EXAHYPE_CCZ4_ADMBASE__ */
