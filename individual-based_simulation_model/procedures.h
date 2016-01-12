/*
 * Simulation program for Kubisch et al. (2016) -
 * "The downward spiral: eco-evolutionary feedback loops lead to the
 *  emergence of 'elastic' ranges."
 *
 * procedure simplifications
 *
 * Dr. Alexander Kubisch
 * University of Hohenheim
 * Institute for Landscape and Plant Ecology
 * Stuttgart, Germany
 * akubisch@posteo.de
 * 2016
 * © MIT License
 */

const gsl_rng *gBaseRand;

//________________________________________________________________________________________
//------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);

	srand(randSeed);
	unsigned long r = random();
	gsl_rng_set(gBaseRand, r);
}

//________________________________________________________________________________________
//-------------------------------------------------------------------------Simplifications

//-------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

//---------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

//-----------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}

//---------------------------------------------------------------Simplify Lognormal Random

double lognorm(double zeta, double sigma)
{
	double var;											//variance of resulting randoms
	double mu,s_d;										//mean and sd of lognormal distr.
														//to be calculated by mean and
														//sigma of resulting randoms
	var = sigma*sigma;
	s_d = sqrt(log((var/(2*zeta))+1));
	mu = log(zeta)-0.5*(s_d*s_d);
	return gsl_ran_lognormal(gBaseRand,mu,s_d);
}
