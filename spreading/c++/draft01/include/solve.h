/* TIME INTEGRATION
 *  Time evolve the system for the solution vector u. Time-advance using
 *  a centered difference scheme.
 *
 * REFERENCES
 *  Lopez et al, J Coll Int Sci (1976) - gravitational spreading
 *  Moriarty, Tuck, and Schwartz, Phys Fluids A (1991) - gravitational thinning w/surface tension
 *  von Rosenberg, Methods for the Numerical Solution of PDEs (1969)
 *  
 * PARAMETERS
 *  J		[input]			number of grid points
 *  p		[input]			parameters
 *  u		[input]			solution vector
 */

/* The following problem is adapted from Moriarty, Tuck, and Schwartz, Physics
 * of Fluids A (1991). The differential system is
 *   u_t = - q_x,
 *     q = (u^3/3)*(1 + (1/B)*u_xxx)
 * where u is the film thickness, q is the flux, B is the Bond number.
 * 
 * NOTE: The independent variable u is shifted 1/2 step from the nodal points,
 *       so that u_{j+1/2,n} = u(x_{j}, t_{n}). The vector u has J+2 elements,
 *       and there are J+1 nodal points x_{j} where j = 0, 1, ..., J-1, J.
 *       The elements u_{0} and u_{J+1} lie outside the domain.
 */

#ifndef SOLVE_H
#define SOLVE_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "./lalg.h"

using namespace std;

/* PROTOTYPES */
void solve_tevol (int, int, double, double, double *, double *, vector<double> &);
void solve_tstep (int, double, double, double *, double *, double *);
void solve_coeffs(int, double, double, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *);
void solve_radius(int, double, double *, double *, double &);

/* IMPLEMENTATIONS */

// time evolve u
void solve_tevol(int N, int J, double dt, double dx, double *p, double *u0, vector<double> &u){
	int n, j;
	double v0[J+1], v1[J+1];
	
	// initialize
	for (j = 0; j < J+1; j++){
		v0[j] = u0[j];
		v1[j] = u0[j];
	}

	// time evolution
	for (n = 0; n < N+1; n++){	
		cout << "ts = " << n << " / " << N << ", dt/dx = " << dt/dx << endl;
		for (j = 0; j < J+1; j++)
			u.push_back(v0[j]);

		solve_tstep(J, dt, dx, p, v0, v1);

		for (j = 0; j < J+1; j++)
			v0[j] = v1[j];
	}
}

// advance u from t_{n} to t_{n+1}
void solve_tstep(int J, double dt, double dx, double *p, double *u0, double *u1){
	int i, j;
	int M = J - 1; // system size
	double v0[J+1], v1[J+1];
	double w[M];
	double a[M], b[M], c[M], d[M], e[M], f[M];

	// parameters 
	double visc = p[0];		// viscosity (scaled by body force)
	double tens = p[1];		// surface tension (scaled by body force)
	double hinf = p[2];		// right boundary value
	double angl = p[3];		// right slope
	double vlme = p[4];		// zeroth moment of u (volume)

	// initialize
	double r, R;
	for (j = 0; j < J+1; j++)
		v0[j] = u0[j];

	// calculate radius
	solve_radius(J, dx, p, v0, r);
	R = angl*r*dx/2.0;
	
	// prediction: project v0 onto intermediate time point to get v1
	solve_coeffs(J, 0.5*dt, dx, p, v0, v0, a, b, c, d, e, f);
	lalg_pent   (M, a, b, c, d, e, f, w);
	
	for (i = 0; i < M; i++)
		v1[i] = w[i];
	v1[M  ] = hinf + R;
	v1[M+1] = hinf - R;

	// recalculate radius and update boundary values
	solve_radius(J, dx, p, v1, r);
	R = angl*r*dx/2.0;
	v1[M  ] = hinf + R;
	v1[M+1] = hinf - R;
	
	// correction: correct v1 by taking the full time step
	solve_coeffs(J,     dt, dx, p, v0, v1, a, b, c, d, e, f);
	lalg_pent   (M, a, b, c, d, e, f, w);
	
//	// no predictor-corrector: just time advance using previous coefficients
//	solve_coeffs(J,     dt, dx, p, v0, v0, a, b, c, d, e, f);
//	lalg_pent   (M, a, b, c, d, e, f, w);
	
	for (i = 0; i < M; i++)
		u1[i] = w[i];
	u1[M  ] = hinf + R;
	u1[M+1] = hinf - R;
	
	// recalculate radius and update boundary values
	solve_radius(J, dx, p, u1, r);
	R = angl*r*dx/2.0;
	u1[M  ] = hinf + R;
	u1[M+1] = hinf - R;
}

// compute pentadiagonal coefficients
/* u0[j], v[j], j = 0, 1, ..., J
 * BCs:	u0[ -1] = u0[0], v[-1] = v[0]
 *			u0[J-1] = v[J-1] = hinf + R
 *      u0[J  ] = v[J  ] = hinf - R
 */
void solve_coeffs(int J, double dt, double dx, double *p, double *u0, double *v,
                 double *a, double *b, double *c, double *d, double *e, double *f){
	int i, j;
	int M = J - 1; // system size
	double dx2  = dx*dx;

	// parameters 
	double visc = p[0];		// viscosity (scaled by body force)
	double tens = p[1];		// surface tension (scaled by body force)
	double hinf = p[2];		// right boundary value
	double angl = p[3];		// right slope
	double vlme = p[4];		// zeroth moment of u (volume)

	// calculate radius of drop
	double r, r2;
	solve_radius(J, dx, p, v, r);
	
	// scalar coefficients
	double T, S, R;
	T = 2.0*dt/(r2*dx2);
	S = 2.0*tens/(r2*dx2);
	R = angl*r*dx/2.0;
	r2 = r*r;

	// diffusion coefficients
	double k1[M], k2[M];
	k1[0] = 0.0;
	k2[0] = pow(0.5*(v[0  ] + v[1  ]),3.0)/(3.0*visc);
	for (j = 1; j < J-1; j++){
		k1[j] =  j   *pow(0.5*(v[j-1] + v[j  ]),3.0)/(3.0*visc);
		k2[j] = (j+1)*pow(0.5*(v[j  ] + v[j+1]),3.0)/(3.0*visc);
	}
	k1[J-1] = (J-1)*pow(0.5*(v[J-2] + hinf + R),3.0)/(3.0*visc);
	k2[J-1] =  J   *pow(              hinf     ,3.0)/(3.0*visc);

	// pentadiagonal coefficients
	double cfj;
	for (j = 0; j < J-1; j++){
		cfj = -dt/(double(2*j+1)*r2*dx2);

		a[j] =        - S*(double(  j-1)/double(2*j-1)) *k1[j]                                                ;
		b[j] =   (1.0 + S*(double(3*j+1)/double(2*j+1)))*k1[j]        + S*(double(  j  )/double(2*j+1)) *k2[j];
		c[j] = - (1.0 + S*(double(3*j-1)/double(2*j-1)))*k1[j] - (1.0 + S*(double(3*j+4)/double(2*j+3)))*k2[j];
		d[j] =          S*(double(  j+1)/double(2*j+1)) *k1[j] + (1.0 + S*(double(3*j+2)/double(2*j+1)))*k2[j];
		e[j] =                                                        - S*(double(  j+2)/double(2*j+3)) *k2[j];
		
		a[j] *= cfj; 
		b[j] *= cfj; 
		c[j] *= cfj; 
		d[j] *= cfj; 
		e[j] *= cfj; 
		
		f[j] = 0.0;
	}

	// apply boundary conditions
	j = 0;
	c[j] += b[j];
	b[j]  = 0.0;
	a[j]  = 0.0;
	f[j] += -((c[j] - 1.0)*u0[j] + d[j]*u0[j+1] + e[j]*u0[j+2]);

	j = 1;
	b[j] += a[j];
	a[j]  = 0.0;
	f[j] += -(b[j]*u0[j-1] + (c[j] - 1.0)*u0[j] + d[j]*u0[j+1] + e[j]*u0[j+2]);

	j = J-3;
	f[j] += -e[j]*(2.0*hinf + R);
	e[j]  = 0.0;
	f[j] += -(a[j]*u0[j-2] + b[j]*u0[j-1] + (c[j] - 1.0)*u0[j] + d[j]*u0[j+1]);

	j = J-2;
	f[j] += -d[j]*(2.0*hinf + R) - e[j]*(2.0*hinf - R);
	d[j]  = 0.0;
	e[j]  = 0.0;
	f[j] += -(a[j]*u0[j-2] + b[j]*u0[j-1] + (c[j] - 1.0)*u0[j]);

	// finish assembling the right-hand side
	for (j = 2; j < J-3; j++)
		f[j] += -(a[j]*u0[j-2] + b[j]*u0[j-1] + (c[j] - 1.0)*u0[j] + d[j]*u0[j+1] + e[j]*u0[j+2]);

	// finish assembling the main diagonal
	for (j = 0; j < J-1; j++)
		c[j] += 1.0;

	// FOR DEBUGGING
	//printf("\n");
	//for (i = 0; i < M; i++){
	////	printf("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", a[i], b[i], c[i], d[i], e[i], f[i]);
	//	printf("%.4f ", c[i]);
	//}
	//printf("\n");
}

void solve_radius(int J, double dx, double *p, double *u, double &r){
	int j;
	double dx2 = dx*dx;
	double sum = 0.0  ;

	// parameters 
	double visc = p[0];		// viscosity (scaled by body force)
	double tens = p[1];		// surface tension (scaled by body force)
	double hinf = p[2];		// right boundary value
	double angl = p[3];		// right slope
	double vlme = p[4];		// zeroth moment of u (volume)

	// compute sum
	for (j = 1; j < J-1; j++)
		sum += j*(u[j-1] + u[j]);
	sum += J*hinf;

	r  = sqrt(vlme/(M_PI*dx2*sum));
}


#endif
