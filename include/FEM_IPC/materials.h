#ifndef MATERIALS_H
#define MATERIALS_H

#include "utils.h" 


// Struct of material
struct Material
{
	double density = 1000; // particle density
	double E = 1.0E4; // Young's modulus
	double nu = 0.3; //Poisson ratio

	double mu = E / (2.0 * (1.0 + nu)); // lame parameter mu / shear modulus  ::::only for isotropic material
	double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); // lame parameter lambda  ::::only for isotropic material
	double K = 2.0 / 3.0 * mu + lambda; // bulk modulus  ::::only for isotropic material


	// local damage field parameters
	double thetaf = 480; // the reference stress for conventional phase field method
	double Gf = 3;
	double lch = sqrt(2) * 1.0;
	double HsBar = thetaf * thetaf / 2.0 / E / Gf;
	double Hs = HsBar * lch / (1.0 - HsBar * lch);

	double thetaF = 480; // the reference stress for no-weakening phase field


	// return mapping stress threshold(only for bending stress)
	double bendingStressThreshold = 1.0E6;


	void updateDenpendecies()
	{
		mu = E / (2.0 * (1.0 + nu)); // lame parameter mu / shear modulus  ::::only for isotropic material
		lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); // lame parameter lambda  ::::only for isotropic material
		K = 2.0 / 3.0 * mu + lambda; // bulk modulus  ::::only for isotropic material

		HsBar = thetaf * thetaf / 2.0 / E / Gf;
		Hs = HsBar * lch / (1.0 - HsBar * lch);
	}

};

#endif
