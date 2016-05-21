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
	double H = 1.0;
	double X = 1000;

	double visc = 1.0;
	double tens = 0.0;
	double hinf = 0.0000001;
	double angl = 0.0;
	double vlme = X*X*(M_PI*hinf + 1.45352*H);
	
	double p[7];
	p[0] = visc;
	p[1] = tens;
	p[2] = hinf;
	p[3] = angl;
	p[4] = vlme;
	p[5] = H   ;
	p[6] = X   ;

	// time and space discretization
	int N = 21;
	int J = 50;
	double xmax = 1.0;
	double tmax = 1.0;
	double dx = xmax/double(J);
	double dt = tmax/double(N);

	// write parameters
	string dir = "../output";	// output directory
	string fn = "data";				// output file name
	int nw = N;								// number of timesteps to write

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
	double hinf = p[2];
	double H    = p[5];
	double X    = p[6];
	
	// initialize film thickness
	for (j = 0; j < J+2; j++){
		x = (double(j) + 0.5)*dx;
		u0[j] = hinf + H*gsl_sf_cos(M_PI*x/(2.0));
	}
}

