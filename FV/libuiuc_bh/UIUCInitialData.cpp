// This is the "real-standalone" C++ version of the UIUC Initial Data
// It has no external dependencies, it is just plain C++.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>
#include <sstream>

#include "UIUCInitialData.h"

#define ERF(X, X0,W) (0.5*(erf( ( (X) - (X0) )/ (W) ) + 1.0))

// Calculate the 3D (spatial) Lorentz transformation and its inverse
void UIUC_BlackHole::Lorentz_Transformation_3D(double LT[3][3], double iLT[3][3]) const
{

  double KD[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; // 3D Kronecker delta
  double Pmag = sqrt(par_P[0] * par_P[0] + par_P[1] * par_P[1] + par_P[2] * par_P[2]); // Momentum magnitude
  double v[3]; // Three-velocity

  for(int i = 0; i < 3; i++)
    {
      v[i] = par_P[i] / sqrt(par_M * par_M + Pmag * Pmag);
    }

  double W = pow(1.0 - v[0] * v[0] - v[1] * v[1] - v[2] * v[2], -0.5); // Lorentz factor

  for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  LT[i][j] = KD[i][j] + W * W / (1.0 + W) * v[i] * v[j]; // Lorentz transformation
	  iLT[i][j] = KD[i][j] - W / (1.0 + W) * v[i] * v[j]; // Inverse Lorentz transformation
	}
    }

  return;
}



// University of Illinois Urbana-Champaign Kerr initial data.
// UIUC definitions appear in Liu et al PRD 80 121503(R) (2009).
void UIUC_BlackHole::UIUC_ID(double x, double y, double z, double g[3][3], double K[3][3], double N[4]) const
{

  double epsilon = 1.0e-10;

  double M = par_M; // Mass parameter
  double a = M * par_chi; // Spin per unit mass, in terms of dimensionless spin parameter par_chi
  double rm = M - sqrt(M * M - a * a); // UIUC inner horizon, equation appears in text after Eq (1) in Liu
  double rp = M + sqrt(M * M - a * a); // UIUC outer horizon, equation appears in text after Eq (1) in Liu

  // Calculate Lorentz transformation matrix
  double LT[3][3], iLT[3][3];
  Lorentz_Transformation_3D(LT, iLT);

  double xf[3] = {x, y, z}; // Field point
  double xS[3]; // Stationary coordinates

  for(int i = 0; i < 3; i++)
    {
      xS[i] = 0;

      for(int j = 0; j < 3; j++)
	{
	  xS[i] += iLT[i][j] * xf[j]; // Transform field point from boosted to stationary frame
	}
    }

  double rUIUC = sqrt(xS[0] * xS[0] + xS[1] * xS[1] + xS[2] * xS[2]); // UIUC radial coordinate
  double r_small = 0.001;
  
  // If we're too close to the coordinate singularity at the puncture, shift away
  if(rUIUC < r_small && avoid_puncture)
    {
      rUIUC = r_small;
    }

  double th = acos(xS[2] / rUIUC); // UIUC latitude coordinate
  double ph = atan2(xS[1], xS[0]); // UIUC longitude coordinate
  double Sth = sin(th), Cth = cos(th), Tth = tan(th), Sph = sin(ph), Cph = cos(ph);
  double rBL = rUIUC * pow(1 + 0.25 * rp / rUIUC, 2); // Boyer-Lindquist radial coordinate, Liu Eq (11)
  double SIG = rBL * rBL + pow(a * Cth, 2); // Boyer-Lindquist "Sigma", equation appears in text after Eq (2) in Liu
  double DEL = rBL * rBL - 2 * M * rBL + a * a; // Boyer-Lindquist "Delta", equation appears in text after Eq (2) in Liu
  double A = pow(rBL * rBL + a * a, 2) - DEL * pow(a * Sth, 2); // equation appears in text after Eq (2) in Liu

  // Physical spatial metric in spherical basis, Liu Eq (13)
  double grr = SIG * pow(rUIUC + 0.25 * rp, 2) * pow(rUIUC, -3) / (rBL - rm);
  double gthth = SIG;
  double gphph = A * Sth * Sth / SIG;

  // Physical extrinsic curvature in spherical basis, Liu Eqs (14) and (15)
  double Krph = M * a * Sth * Sth * (3 * pow(rBL, 4) + 2 * pow(a * rBL, 2) - pow(a, 4) - (rBL * rBL - a * a) * pow(a * Sth, 2)) * (1 + 0.25 * rp / rUIUC) / (SIG * sqrt(A * SIG * rUIUC * (rBL - rm)));
  double Kthph = - 2 * pow(a, 3) * M * rBL * Cth * pow(Sth, 3) * (rUIUC - 0.25 * rp) * sqrt(rBL - rm) / (SIG * sqrt(A * SIG * rUIUC));

  double g0[3][3], K0[3][3]; // Stationary spatial metric and extrinsic curvature

  // Near the z-axis
  if(rUIUC * Sth < epsilon)
    {
      g0[0][0] = (a * a + pow(rp + 4 * fabs(xS[2]), 4) / (256 * xS[2] * xS[2])) / (xS[2] * xS[2]);
      g0[0][1] = 0;
      g0[0][2] = 0;
      g0[1][1] = g0[0][0];
      g0[1][2] = 0;
      g0[2][2] = -pow(rp + 4 * fabs(xS[2]), 2) * g0[0][0] / (a * a - 2 * (M * rp + 8 * xS[2] * xS[2]) + 8 * fabs(xS[2]) * (M - 3 * sqrt(M * M - a * a)));
      
      g0[1][0] = g0[0][1];
      g0[2][0] = g0[0][2];
      g0[2][1] = g0[1][2];

      int i, j;

      for(i = 0; i < 3; i++)
	{
	  for(j = 0; j < 3; j++)
	    {
	      K0[i][j] = 0;
	    }
	}
    }
  else
    {
      // ADMBase physical spatial metric (Cartesian basis)
      g0[0][0] = (gthth * pow(Cph * Cth, 2)) / pow(rUIUC, 2) + (gphph * pow(Sph, 2)) / (pow(rUIUC, 2) * pow(Sth, 2)) + grr * pow(Cph, 2) * pow(Sth, 2);
      g0[0][1] = (gthth * Cph * pow(Cth, 2) * Sph) / pow(rUIUC, 2) - (gphph * Cph * Sph) / (pow(rUIUC, 2) * pow(Sth, 2)) + grr * Cph * Sph * pow(Sth, 2);
      g0[0][2] = grr * Cph * Cth * Sth - (gthth * Cph * Cth * Sth) / pow(rUIUC, 2);
      g0[1][1] = (gphph * pow(Cph, 2)) / (pow(rUIUC, 2) * pow(Sth, 2)) + (gthth * pow(Cth, 2) * pow(Sph, 2)) / pow(rUIUC, 2) + grr * pow(Sph, 2) * pow(Sth, 2);
      g0[1][2] = grr * Cth * Sph * Sth - (gthth * Cth * Sph * Sth) / pow(rUIUC, 2);
      g0[2][2] = grr * pow(Cth, 2) + (gthth * pow(Sth, 2)) / pow(rUIUC, 2);
      
      g0[1][0] = g0[0][1];
      g0[2][0] = g0[0][2];
      g0[2][1] = g0[1][2];

      // ADMBase physical extrinsic curvature (Cartesian basis)
      K0[0][0] = (-2 * Krph * Cph * Sph) / rUIUC - (2 * Kthph * Cph / Tth * Sph) / pow(rUIUC, 2);
      K0[0][1] = (Krph * pow(Cph, 2)) / rUIUC + (Kthph * pow(Cph, 2) / Tth) / pow(rUIUC, 2) - (Krph * pow(Sph, 2)) / rUIUC - (Kthph / Tth * pow(Sph, 2)) / pow(rUIUC, 2);
      K0[0][2] = (Kthph * Sph) / pow(rUIUC, 2) - (Krph / Tth * Sph) / rUIUC;
      K0[1][1] = (2 * Krph * Cph * Sph) / rUIUC + (2 * Kthph * Cph / Tth * Sph) / pow(rUIUC, 2);
      K0[1][2] = -((Kthph * Cph) / pow(rUIUC, 2)) + (Krph * Cph / Tth) / rUIUC;
      K0[2][2] = 0;

      K0[1][0] = K0[0][1];
      K0[2][0] = K0[0][2];
      K0[2][1] = K0[1][2];
    }

  // Lorentz transformation
  for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  g[i][j] = 0;
	  K[i][j] = 0;

	  for(int k = 0; k < 3; k++)
	    {
	      for(int l = 0; l < 3; l++)
		{
		  g[i][j] += LT[i][k] * LT[j][l] * g0[k][l];
		  K[i][j] += LT[i][k] * LT[j][l] * K0[k][l];
		}
	    }
	}
    }

  // Determinant of spatial metric
  double detg = -g[0][2] * g[0][2] * g[1][1] + 2 * g[0][1] * g[0][2] * g[1][2] - g[0][0] * g[1][2] * g[1][2] - g[0][1] * g[0][1] * g[2][2] + g[0][0] * g[1][1] * g[2][2];
  double psi = pow(detg, 1.0/12.0);

  // Conformal factor, equation appears in last paragraph of Sec. IIB in Liu

  double alpha;
  double beta_x;
  double beta_y;
  double beta_z;

  if(lapse_type == 0)
    {
      // stationary
      alpha = sqrt(DEL * SIG / A); // Lapse function, Liu Eq (6)
    }
  else if(lapse_type == 1)
    {
      // unit
      alpha = 1.0;
    }
  else if(lapse_type == 2)
    {
      // trumpet
      alpha = 1.0 / (2.0 * psi - 1.0);
    }
  else if(lapse_type == 3)
    {
      // iso_schw
      alpha = 2.0 / (1.0 + pow(psi, 4));
    }
  else if(lapse_type == 4)
    {
      // puncture_r2
      alpha = pow(psi, -2);
    }

  if(shift_type == 0)
    {
      // zero
      beta_x = 0;
      beta_y = 0;
      beta_z = 0;
    }
  else if(shift_type == 1)
    {
      // stationary
      double beta_phi = - 2 * M * a * rBL / A; // Shift vector in spherical basis, Liu Eq (7)
      beta_x = - beta_phi * rUIUC * Sth * Sph;
      beta_y =   beta_phi * rUIUC * Sth * Cph;
      beta_z = 0; 
    }
  else if(shift_type == 2)
    {
      double beta_r = beta_r_interior * exp(-(rUIUC-beta_r_x0)*(rUIUC-beta_r_x0)/(beta_r_w*beta_r_w)); //ERF(rUIUC, beta_r_x0, beta_r_w);
      beta_x = beta_r * Sth * Cph;
      beta_y = beta_r * Sth * Sph;
      beta_z = beta_r * Cth;
    }

  double N0[3] = {beta_x, beta_y, beta_z};

  // Lorentz transformation
  for(int i = 0; i < 3; i++)
    {
      N[i+1] = 0;

      for(int j = 0; j < 3; j++)
	{
	  N[i+1] += LT[i][j] * N0[j];
	}
    }

  // Lapse and shift
  N[0] = alpha;
  //N[1] = beta_x;
  //N[2] = beta_y;
  //N[3] = beta_z;

  return;
}

// Quasi-isotropic Kerr initial data
// QI definitions appear in Hannam et al CQG 24 S15 (2007)
void UIUC_BlackHole::QI_ID(double x, double y, double z, double g[3][3], double K[3][3], double N[4]) const
{

  double epsilon = 1.0e-10;

  double M = par_M; // Mass parameter
  double a = M * par_chi; // Spin per unit mass, in terms of dimensionless spin parameter par_chi

  // Calculate Lorentz transformation matrix
  double LT[3][3], iLT[3][3];
  Lorentz_Transformation_3D(LT, iLT);

  double xf[3] = {x, y, z}; // Field point
  double xS[3]; // Stationary coordinates

  for(int i = 0; i < 3; i++)
    {
      xS[i] = 0;

      for(int j = 0; j < 3; j++)
	{
	  xS[i] += iLT[i][j] * xf[j]; // Transform field point from boosted to stationary frame
	}
    }

  double rQI = sqrt(xS[0] * xS[0] + xS[1] * xS[1] + xS[2] * xS[2]); // Quasi-isotropic radius
  double r_small = 0.001;
  
  // If we're too close to the coordinate singularity at the puncture, shift away
  if(rQI < r_small && avoid_puncture)
    {
      rQI = r_small;
    }

  double th = acos(xS[2] / rQI); // Quasi-isotropic latitude
  double ph = atan2(xS[1], xS[0]); // Quasi-isotropic longitude
  
  double rBL = rQI + M + 0.25 * (M * M - a * a) / rQI; // Boyer-Lindquist radius, Eq (22) in Hannam et al.
  double SIG = rBL * rBL + a * a * cos(th) * cos(th); // Eq (21) in Hannam et al.
  double sig = 2 * M * rBL / SIG; // Eq (23) in Hannam et al.
  double h = (1 + sig) / (SIG * rQI * rQI); // Eq (20) in Hannam et al.
  double psi = pow(SIG / (rQI * rQI), 0.25); // Quasi-isotropic conformal factor, Eq (18) in Hannam et al.

  double DEL = rBL * rBL - 2 * M * rBL + a * a;
  double A = pow(rBL * rBL + a * a, 2) - DEL * a * a * sin(th) * sin(th);

  double eq2 = SIG / (rBL * rBL + a * a * (1 + sig * sin(th) * sin(th))); // Eq (26) in Hannam et al.
  double HE = sqrt(eq2) * a * M * ((rBL * rBL - a * a) * SIG + 2 * rBL * rBL * (rBL * rBL + a * a)) / (SIG * SIG); // Eq (27) in Hannam et al.
  double HF = 0.5 * sqrt(eq2) * a * a * a * M * rBL * (M * M - a * a - 4 * rQI * rQI) * cos(th) * sin(th) * sin(th) / (rQI * SIG * SIG); // Eq (28) in Hannam et al.
  double Arph = HE * sin(th) * sin(th) / (rQI * rQI); // Eq (24) in Hannam et al.
  double Athph = HF * sin(th) / rQI; // Eq (25) in Hannam et al.

  // Kronecker delta
  double delta[3][3] = {{1, 0, 0},
			   {0, 1, 0},
			   {0, 0, 1}};

  // Eq (17) in Hannam et al.
  double v[3][3] = {{ xS[1] * xS[1], -xS[0] * xS[1], 0},
		       {-xS[0] * xS[1],  xS[0] * xS[0], 0},
		       {             0,              0, 0}};

  int i, j;
  double g0[3][3], K0[3][3]; // Stationary spatial metric and extrinsic curvature

  // Evaluate spatial metric
  for(i = 0; i < 3; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  g0[i][j] = pow(psi, 4) * (delta[i][j] + a * a * h * v[i][j]); // Eq (16) in Hannam et al.
	}
    }

  // Evaluate extrinsic curvature
  // Near the spin axis
  if(rQI * sin(th) < epsilon)
    {
      for(i = 0; i < 3; i++)
	{
	  for(j = 0; j < 3; j++)
	    {
	      K0[i][j] = 0;
	    }
	}
    }
  else
    {
      K0[0][0] = pow(psi, -2) * (-2 * cos(ph) * (Arph * rQI + Athph * (1.0 / tan(th))) * sin(ph)) / (rQI * rQI);
      K0[0][1] = pow(psi, -2) * (cos(2 * ph) * (Arph * rQI + Athph * (1.0 / tan(th)))) / (rQI * rQI);
      K0[0][2] = pow(psi, -2) * ((Athph - Arph * rQI * (1.0 / tan(th))) * sin(ph)) / (rQI * rQI);
      K0[1][1] = pow(psi, -2) * (2 * cos(ph) * (Arph * rQI + Athph * (1.0 / tan(th))) * sin(ph)) / (rQI * rQI);
      K0[1][2] = pow(psi, -2) * (cos(ph) * (-Athph + Arph * rQI * (1.0 / tan(th)))) / (rQI * rQI);
      K0[2][2] = 0;

      K0[1][0] = K0[0][1];
      K0[2][0] = K0[0][2];
      K0[2][1] = K0[1][2];
    }

  // Lorentz transformation
  for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  g[i][j] = 0;
	  K[i][j] = 0;

	  for(int k = 0; k < 3; k++)
	    {
	      for(int l = 0; l < 3; l++)
		{
		  g[i][j] += LT[k][i] * LT[l][j] * g0[k][l];
		  K[i][j] += LT[k][i] * LT[l][j] * K0[k][l];
		}
	    }
	}
    }

  double alpha;
  double beta_x;
  double beta_y;
  double beta_z;

  if(lapse_type == 0)
    {
      // stationary
      alpha = pow(1 + 2 * M * rBL * (rBL * rBL + a * a) / (DEL * SIG), -0.5);
    }
  else if(lapse_type == 1)
    {
      // unit
      alpha = 1.0;
    }
  else if(lapse_type == 2)
    {
      // trumpet
      alpha = 1.0 / (2.0 * psi - 1.0);
    }
  else if(lapse_type == 3)
    {
      // iso_schw
      alpha = 2.0 / (1.0 + pow(psi, 4));
    }
  else if(lapse_type == 4)
    {
      // puncture_r2
      alpha = pow(psi, -2);
    }

  if(shift_type == 0)
    {
      // zero
      beta_x = 0;
      beta_y = 0;
      beta_z = 0;
    }
  else if(shift_type == 1)
    {
      // stationary
      double beta_phi = - 2 * M * a * rBL / A; // Shift vector in spherical basis, Liu Eq (7)
      beta_x = - beta_phi * rQI * sin(th) * sin(ph);
      beta_y =   beta_phi * rQI * sin(th) * cos(ph);
      beta_z = 0; 
    }
  else if(shift_type == 2)
    {
      double beta_r = beta_r_interior * exp(-(rQI-beta_r_x0)*(rQI-beta_r_x0)/(beta_r_w*beta_r_w)); //ERF(rQI, beta_r_x0, beta_r_w);
      beta_x = beta_r * sin(th) * cos(ph);
      beta_y = beta_r * sin(th) * sin(ph);
      beta_z = beta_r * cos(th);
    }

  double N0[3] = {beta_x, beta_y, beta_z};

  // Lorentz transformation
  for(int i = 0; i < 3; i++)
    {
      N[i+1] = 0;

      for(int j = 0; j < 3; j++)
	{
	  N[i+1] += LT[i][j] * N0[j];
	}
    }

  // Lapse and shift
  N[0] = alpha;
  //N[1] = beta_x;
  //N[2] = beta_y;
  //N[3] = beta_z;

  return;
}


// Calculate physical spatial metric and extrinsic curvature of Kerr BH.
void UIUC_BlackHole::InitialData(double pos[3], double g[3][3], double K[3][3], double N[4]) const {
	// Cartesian coordinates with puncture offset
	double pp[3];
	for(int i=0; i<3; i++) pp[i] = pos[i] - center_offset[i];


	// Pick a gauge
	if(BH_gauge_choice == 0) {
		UIUC_ID(pp[0], pp[1], pp[2], g, K, N);
	} else if(BH_gauge_choice == 1) {
		QI_ID(pp[0], pp[1], pp[2], g, K, N);
	} else {
		std::abort();
	}
}

// string output helpers

#define SS(x) (static_cast<std::ostringstream&>(( std::ostringstream() << x )).str())
#define INDENT "    "
#define NL SS(", " << std::endl << INDENT)
#define DESC(var) SS(#var << " = " << var)
#define MEANS(var, num, text) SS((var == num ? SS(DESC(var) << " = " text) : ""))
#define VEC3(var) SS(#var << " = " << "[" << var[0] << ","	<< var[1] << "," << var[2] << "]")
#define VEC4(var) SS(#var << " = " << "[" << var[0] << ","	<< var[1] << "," << var[2] << "," << var[3] << "]")
#define MAT(var) SS(#var << " = " << "[" << NL << VEC3(var[0]) << NL << VEC3(var[1]) << NL << VEC3(var[2]) << NL <<  "]")

std::string UIUC_BlackHole::toString() const {
	std::stringstream out;
	out << "UIUC InitialData(" << NL;
	
	out
	<< MEANS(lapse_type, 0, "stationary")
	<< MEANS(lapse_type, 1, "unit")
	<< MEANS(lapse_type, 2, "trumpet")
	<< MEANS(lapse_type, 3, "iso_schwarzschild")
	<< MEANS(lapse_type, 4, "puncture")
	<< NL;

	out
	<< MEANS(shift_type, 0, "zero")
	<< MEANS(shift_type, 1, "stationary")
	<< MEANS(shift_type, 2, "radial Gaussian")
	<< NL;
	
	out
	<< MEANS(BH_gauge_choice, 0, "UIUC Liu - Shapiro - Etienne 2009")
	<< MEANS(BH_gauge_choice, 1, "Quasi-isotropic")
	<< NL;

	out << "Puncture location (for a single BH) = " << VEC3(center_offset) << NL;
	out << DESC(avoid_puncture) << NL;
	out << "BH Mass in M_sun: " << DESC(par_M) << NL;
	out << "Dimensionless spin J/M^2: " << DESC(par_chi) << NL;
	out << "For a single BH = momentum of the puncture = " << VEC3(par_P) << NL;
	out << "Constant value of radial shift near puncture: " << DESC(beta_r_interior) << NL;
	out << "Radial shift falloff radius: " << DESC(beta_r_x0) << NL;
	out << "Radial shift falloff width: " << DESC(beta_r_w) << NL;
	
	out << ")"; // closing bracket of InitialData()

	return out.str();
}

std::string UIUC_BlackHole::outputToString(double g[3][3], double K[3][3], double N[4]) {
	std::stringstream out;
	
	out << "3-Metric: " << MAT(g) << "\n";
	out << "3-Extrinsic Curvature: " << MAT(K) << "\n";
	out << "(Lapse,Shift): " << VEC4(N) << "\n";
	
	return out.str();
}

