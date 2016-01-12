/*
 * Simulation program for Kubisch et al. (2015) - 
 * "The downward spiral: Eco-evolutionary feedback loops lead to the
 *  emergence of 'elastic' ranges."
 * 
 * main source file
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

#include <iostream>
#include <cstdlib>							//standard C library
#include <ctime>							//access system time library
#include <fstream>							//file streaming library
#include <string>							//string library included
#include <sstream>							//string streaming for reading numbers from
//infiles
#include <vector>
#include <cmath>							//standard math library

#include <gsl/gsl_rng.h>					//gsl random number generator
#include <gsl/gsl_randist.h>				//gsl random distributions
#include <gsl/gsl_statistics.h>				//gsl some statistical methods
#include <gsl/gsl_statistics_int.h>			//gsl some integer statistical methods
#include <gsl/gsl_sort_int.h>				//gsl sort integer arrays
#include <gsl/gsl_math.h>					//provides additional standard math functions
#include <gsl/gsl_histogram.h>                          //provides histogram functions needed for relatedness

using namespace std;

#include "procedures.h"						//procedure simplifications
#include "classes.h"						//class definitions

//Variables

vector<vector<tPatch> > world;				//simulated landscape
tParameters	par;							//parameters

fstream 	parinfile("data/para.in", ios::in);		//parameter infile
fstream		outfile("data/output.out", ios::out);	//results outfile

///////////////////////////////////////////////////////////////////////////////////////
//---------------------------------SET PARAMETERS------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

void set_parameters() {					//read in standard simulation parameters

	string buffer;
	istringstream is;
	getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> par.capacity;					//habitat capacity
	getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> par.lambda_0;					//per capita growth rate
	getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> par.sigma;					//environmental stochasticity
	getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> par.disp_mort;				//dispersal mortality
	getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> par.ext_prob;					//extinction probability
	getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> par.mut_prob;					//mutation probability
	getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> par.tmax;						//simulated time
	getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> par.allee;	//Allee effect strength
}

///////////////////////////////////////////////////////////////////////////////////////
//-----------------------------INITIALIZATION FUNCTIONS------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

void apply_gradient() { // this function applies the according gradients, here exemplarily the dispersal mortality gradient
	for (int x=0; x<BURN_IN_REGION; ++x) { // no gradient within the burn-in region
		for (int y=0; y<YMAX; ++y) {
			world[x][y].mu = par.disp_mort;
			world[x][y].K = par.capacity;
			world[x][y].lambda_0 = par.lambda_0;
			world[x][y].ext = par.ext_prob;
			world[x][y].sigma = par.sigma;
		}
	}
	for (int x=BURN_IN_REGION; x<XMAX+BURN_IN_REGION; ++x) { // afterwards it starts
		for (int y=0; y<YMAX; ++y) {
			world[x][y].mu = 0.2 + (0.8/double(XMAX))*double(x-BURN_IN_REGION);
			world[x][y].K = par.capacity;
			world[x][y].lambda_0 = par.lambda_0;
			world[x][y].ext = par.ext_prob;
			world[x][y].sigma = par.sigma;
		}
	}
}

void init_individual(tInd *newind) { // a new individual is created
	newind->ld = ran(); // emigration probability is random
}

void init_world() { // the landscape is initialized

	tInd newind;

	world.resize(XMAX+BURN_IN_REGION);
	for (int x=0; x < XMAX+BURN_IN_REGION; x++) {
		world[x].resize(YMAX);
	}

	apply_gradient();

	for (int x = 0; x < XMAX+BURN_IN_REGION; x++)
		for (int y = 0; y < YMAX; y++) {

			world[x][y].indiv.clear();

			if (x<BURN_IN_REGION) { // only in the burn-in region K individuals are created in each population
				for (int i = 0; i < par.capacity; i++) {
					init_individual(&newind);
					world[x][y].indiv.push_back(newind);
				}
			}

		}

}

///////////////////////////////////////////////////////////////////////////////////////
//-------------------------------SIMULATION FUNCTIONS--------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________
//-Dispersal---------------------------------------------------------------------------

double disp_prob(tInd ind) { // calculate dispersal probability

	double trait = ind.ld;
	return(trait);
}

void find_target_patch(int xs, int ys, int *xt, int *yt) { // determine dispersal target patch
	double dist;				//drawn dispersal distance
	double phi;

	dist = 1;						// this distance can be drawn from distributions - kernel implementation

	phi = 2 * M_PI * ran();			//direction is uniformly drawn (0°-360°)

	double xstart = double(xs) + ran() - 0.5; // area-to-area dispersal to avoid artifacts
	double ystart = double(ys) + ran() - 0.5;		

	int xtt = round(xstart + dist*cos(phi)); // actual target determination
	int ytt = round(ystart + dist*sin(phi));

	if (par.burnin) { // restricted dispersal during burn-in
		if (xtt>=BURN_IN_REGION) {
			xtt = xtt-BURN_IN_REGION;
		}
		if (xtt<0) {
			xtt = BURN_IN_REGION+xtt;
		}
	}
	
	if (ytt>=YMAX) { // wrapping in y-direction (needs to be adjusted for kernel implementation)
		ytt = ytt - YMAX;
	}
	if (ytt<0) {
		ytt = BURN_IN_REGION + ytt;
	}

	*xt = xtt;
	*yt = ytt;
}

void dispersal(int x, int y, int xt, int yt, unsigned int i) { //actual dispersal

	if ((xt>=0)&&(xt<XMAX+BURN_IN_REGION)) {
		double mort = (world[x][y].mu+world[xt][yt].mu)/double(2); // calculating effective mortality
		if (ran()>mort) {					//a tube with absorbing x-conditions
			world[xt][yt].immis.push_back(world[x][y].indiv.at(i));	//immigration
		}
	}
}

void dispersal_loop(int t) { // here it is defined, whether individuals disperse
	int xt,yt;

	for (int x = 0; x < XMAX+BURN_IN_REGION; x++)
		for (int y = 0; y < YMAX; y++) {

			int pop = world[x][y].indiv.size();

			if (pop>0) {	//dispersal only, when there are beetles in the patch

				for (unsigned int i = 0; i < world[x][y].indiv.size(); i++) {
					double d = disp_prob(world[x][y].indiv.at(i));

					if (ran()<d) {
						find_target_patch(x,y,&xt,&yt);
						dispersal(x,y,xt,yt,i);		
					} else {
						world[x][y].immis.push_back(world[x][y].indiv.at(i));
					}

				}

			}

		}

	for (int x = 0; x < XMAX+BURN_IN_REGION; x++) // in the end all individuals are replaced by the immigrants
		for (int y = 0; y < YMAX; y++) {
			world[x][y].indiv 	= world[x][y].immis;
			world[x][y].immis.clear();
		}

}

//_____________________________________________________________________________________
//-Reproduction------------------------------------------------------------------------

void genetics(tInd *newborn, tInd mother) { // genetic inheritance of dispersal

	newborn->ld = mother.ld;

	if (ran()<par.mut_prob) {
		newborn->ld += gauss(VARIANCE);
	}
}

void reproduction(int x, int y, double lambda, double survival) { //actual reproduction

	vector<tInd> offspring;
	tInd newborn;

	offspring.clear();

	for (unsigned int i = 0; i < world[x][y].indiv.size(); i++) { // each individual reproduces
		int kids = poisson(lambda*survival);	//children numbers are poisson distributed

		for (int child = 0; child < kids; child++) {

			genetics(&newborn,world[x][y].indiv.at(i)); // inheritance

			offspring.push_back(newborn);

		}
	}

	world[x][y].indiv = offspring;
}

void reproduction_loop() { // density-regulation etc.

	for (int x = 0; x < XMAX+BURN_IN_REGION; x++)
		for (int y = 0; y < YMAX; y++) {

			// implementing environmental stochasticity
			double lambda = lognorm(world[x][y].lambda_0,world[x][y].sigma);

			double dens = double(world[x][y].indiv.size())/double(par.capacity);

			// logistic growth according to the Beverton & Holt-model
			double a = (world[x][y].lambda_0-double(1))/world[x][y].K;
			double b = (dens*dens)/(par.allee*par.allee+dens*dens);
			double survival = b/(1+a*double(world[x][y].indiv.size()));

			if (world[x][y].indiv.size()>0) {
				reproduction(x,y,lambda,survival);
			} else {
				world[x][y].indiv.clear();
			}

		}
}

//_____________________________________________________________________________________
//-Extinction--------------------------------------------------------------------------

void extinction() { // random catastrophic extinction
	for (int x = 0; x < XMAX+BURN_IN_REGION; ++x)
		for (int y = 0; y < YMAX; y++) {
			if (ran()<world[x][y].ext) {
				world[x][y].indiv.clear();
			}
		}
}

///////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------ANALYSIS--------------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

void analyze(int t) { // some simple analysis of the absolute range limit position R

	int R = 0;
	int occ = 0;
	
	// we determine the range margin

	for (int x=0; x<XMAX+BURN_IN_REGION; ++x) {
		occ = 0;
		for (int y=0; y<YMAX; ++y) {
			if (world[x][y].indiv.size()>0) {occ++;}
		}
		if (occ==0) {
			R = x-1;
			break;
		}
	}

	outfile << t << "  " << R << endl;

}

///////////////////////////////////////////////////////////////////////////////////////
//----------------------------------MAIN FUNCTIONS-----------------------------------//
///////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________
//-Initialize--------------------------------------------------------------------------

void initialize() {

	outfile << "t	R" << endl;

	init_world();

}


void simulate() {
	
	par.burnin = true;

	for (int t = 0; t < par.tmax+BURN_IN_TIME; t++) {
		
		if (t == BURN_IN_TIME) {par.burnin = false;}

		if (par.burnin) {cout << "BURNIN! t: " << t << endl;} else {cout << "t: " << t << endl;}

		dispersal_loop(t);
		reproduction_loop();
		extinction();

		if (t%10==0) {analyze(t);}
	}
}

int main() {

	cout.setf(ios_base::scientific);		//exponential data output

	specify_rng(time(NULL));				//initialize random number generator
	
	set_parameters();
	
	initialize();
	
	simulate();

	outfile.close();

	return 0;
}

