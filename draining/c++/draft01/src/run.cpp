/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_bessel.h>
#include "../include/solve.h"
#include "../include/write.h"

using namespace std;

void init(int, double, double *, double *);

int main(){
	// system parameters
	double B  = 1000.0;	// Bond number
	double U0 = 0.0;		// left boundary value
	double U1 = 0.001;		// right boundary value
	
	double p[3];
	p[0] = B ;
	p[1] = U0;
	p[2] = U1;

	// time and space discretization
	int N = 21;
	int J = 50;
	double xmax = 5.0;
	double tmax = 1.0;
	double dx = xmax/double(J);
	double dt = tmax/double(N);

	// write parameters
	string dir = "../output";	// output directory
	string fn = "data";				// output file name
	int nw = 20;								// number of timesteps to write

	// solution vector
	double u0[J+2];
	vector<double> u;

	// initialize
	init(J, dx, p, u0);

	// evolve system in time
	solve_tevol(N, J, dt, dx, p, u0, u);

	// write to file
	write(dir, fn, nw, N, J, dt, dx, u.data());

	return 0;
}

void init(int J, double dx, double *p, double *u0){
	int j;
	double x;
	double U0 = p[1];
	double U1 = p[2];
	
	// initialize H
	for (j = 0; j < J+2; j++){
		x = (double(j) - 0.5)*dx;
		if (x < 1.0) u0[j] = U0 + 1.0 - 4.0*(x - 0.5)*(x - 0.5);
		else         u0[j] = U1;
	}
}

