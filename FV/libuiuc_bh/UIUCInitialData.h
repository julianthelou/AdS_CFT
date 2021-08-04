#ifndef __UIUC_INITIAL_DATA_POINTWISE_STANDALONE__
#define __UIUC_INITIAL_DATA_POINTWISE_STANDALONE__

#include <string>

/**
 * This POD collects the physical arguments for creating the black hole.
 *
 * These are analytic Initial Data, there is no preparation and thus
 * there are no intermediate output data on some grid.
 * 
 * 
 **/
struct UIUC_BlackHole {
	// Choose the initial lapse from parameter initial_lapse_type
	int lapse_type = 0;
	// 0 : stationary
	// 1 : unit
	// 2 : trumpet
	// 3 : iso_schwarzschild
	// 4 : puncture_r2

	// Choose the initial shift from parameter initial_shift_type
	int shift_type = 1;
	// 0 : zero
	// 1 : stationary
	// 2 : radial Gaussian

	int BH_gauge_choice = 0;
	// 0 : UIUC Liu - Shapiro - Etienne 2009
	// 1 : Quasi-isotropic


	//For a single BH = location of the puncture
	double center_offset[3] {
		0., // x offset
		0., // y offset
		0.  // z offset
	};

	//limit r to r_min to not evaluate on the puncture
	bool avoid_puncture = true;

	// BH Mass in M_sun
	double par_M = 2.0;

	//dimensionless spin J/M^2
	double par_chi = 0.0;


	//For a single BH = momentum of the puncture
	double par_P[3] {
		0.4, // x offset
		0.4, // y offset
		0.0  // z offset
	};

	//"Constant value of radial shift near puncture."
	double beta_r_interior = 0;

	// "Radial shift falloff radius."
	double beta_r_x0  = 0.25;

	//"Radial shift falloff width."
	double beta_r_w = -1.;

	// Calculate the 3D (spatial) Lorentz transformation and its inverse
	void Lorentz_Transformation_3D(double LT[3][3], double iLT[3][3]) const;

	/**
	 * University of Illinois Urbana-Champaign Kerr initial data.
	 * UIUC definitions appear in Liu et al PRD 80 121503(R) (2009).
	 * 
	 * @param x: X coordinate where to ask for ID, puncture at origin.
	 * @param y: Y coordinate
	 * @param z: Z coordinate 
	 * @param g: 3-metric (output)
	 * @param K: 3-extrinsic curvature (output)
	 * @param N: lapse and shift (output)
	 **/
	void UIUC_ID(double x, double y, double z, double g[3][3], double K[3][3], double N[4]) const;

	/**
	 * Quasi-isotropic Kerr initial data
	 * QI definitions appear in Hannam et al CQG 24 S15 (2007)
	 *
	 * @param x: X coordinate where to ask for ID, puncture at origin.
	 * @param y: Y coordinate
	 * @param z: Z coordinate
	 * @param g: 3-metric (output)
	 * @param K: 3-extrinsic curvature (output)
	 * @param N: lapse and shift (output)
	 **/
	void QI_ID(double x, double y, double z, double g[3][3], double K[3][3], double N[4]) const;

	/**
	 * Give either UIUC_ID or QI_ID based on the choice of BH_gauge_choice.
	 * Also position the puncture at center_offset.
	 **/
	void InitialData(double pos[3], double g[3][3], double K[3][3], double N[4]) const;
	
	/// Reports about the parameters
	std::string toString() const;
	
	/// Just for convenience, show the content of the fields
	static std::string outputToString(double g[3][3], double K[3][3], double N[4]);
};



#endif /* __UIUC_INITIAL_DATA_POINTWISE_STANDALONE__ */
