/*
 * Simulation program for Kubisch et al. (2016) -
 * "The downward spiral: eco-evolutionary feedback loops lead to the
 *  emergence of 'elastic' ranges."
 *
 * class definitions
 *
 * Dr. Alexander Kubisch
 * University of Hohenheim
 * Institute for Landscape and Plant Ecology
 * Stuttgart, Germany
 * akubisch@posteo.de
 * 2016
 * © MIT License
 */

const int 	YMAX	= 50;				//max. y-dimension of the world
const int 	XMAX 	= 200;				//max. x-dimension of the world
const float VARIANCE	= 0.2;			//variance used mutations

const int   BURN_IN_TIME = 2500;		//how long is the burn-in period?
const int   BURN_IN_REGION = 10;        //how large the area?

class tInd // class for individuals
{
public:
	double  ld;						//one locus with 1 dispersal allele
};

class tPatch
{
public:
	vector<tInd> indiv;			//population
	vector<tInd> immis;			//immigrants
	float		mu;						//dispersal mortality
	float		K;						//carrying capacity
	float		lambda_0;				//per capita growth rate
	float		ext;					//extinction rate
	float		sigma;					//environmental stochasticity
};

class tParameters
{
public:
	int 	capacity;					//habitat capacity
	double	lambda_0;					//per capita growth rate
	double	disp_mort;					//dispersal mortality
	double	sigma;						//environmental fluctuations
	double	ext_prob;					//extinction rate
	double	mut_prob;					//mutation probability
	int		tmax;						//max. no. of generations to be simulated
	bool    burnin;                     //is currently the burn-in running?
	double	allee;						//Allee effect strength
};
