/*
 * Simulation program for Kubisch et al. (2015) - 
 * "The downward spiral: Eco-evolutionary feedback loops lead to the
 *  emergence of 'elastic' ranges."
 * 
 * class definitions
 * 
 * Copyright (C) 2015 Dr. Alexander Kubisch
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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