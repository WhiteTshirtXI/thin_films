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

/* IMPLEMENTATIONS */

// time evolve u
void solve_tevol(int N, int J, double dt, double dx, double *p, double *u0, vector<double> &u){
	int n, j;
	double v0[J+1], v1[J+1];
	
	// initialize
	for (j = 0; j < J+2; j++){
		v0[j] = u0[j];
		v1[j] = u0[j];
	}

	// time evolution
	for (n = 0; n < N+1; n++){	
		cout << "ts = " << n << " / " << N << ", dt/dx = " << dt/dx << endl;
		for (j = 0; j < J+2; j++)
			u.push_back(v0[j]);

		solve_tstep(J, dt, dx, p, v0, v1);

		for (j = 0; j < J+2; j++)
			v0[j] = v1[j];
	}
}


// advance u from t_{n} to t_{n+1}
void solve_tstep(int J, double dt, double dx, double *p, double *u0, double *u1){
	int i, j;
	int R = J - 2; // system size
	double v0[R], v1[R];
	double a[R], b[R], c[R], d[R], e[R], f[R];

	double U0 = p[1];
	double U1 = p[2];

	// initialize
	for (j = 1; j < J; j++){
		i = j-1;
		v0[i] = u0[j];
	}
	
	// prediction: project v0 onto intermediate time point to get v1
	solve_coeffs(J, 0.5*dt, dx, p, v0, v0, a, b, c, d, e, f);
	lalg_pent   (R, a, b, c, d, e, f, v1);

	// correction: correct v1 by taking the full time step
	solve_coeffs(J,     dt, dx, p, v0, v1, a, b, c, d, e, f);
	lalg_pent   (R, a, b, c, d, e, f, v1);
	
	// no predictor-corrector: just time advance using previous coefficients
	solve_coeffs(J,     dt, dx, p, v0, v0, a, b, c, d, e, f);
	lalg_pent   (R, a, b, c, d, e, f, v1);
	
	// copy solution
	for (j = 1; j < J; j++){
		i = j-1;
		u1[j] = v1[i];
	}
	u1[0  ] = -u1[1];
	u1[J  ] =  U1;
	u1[J+1] =  u1[J]; // extra entry b/c u is shifted 1/2 step
}

// compute pentadiagonal coefficients
void solve_coeffs(int J, double dt, double dx, double *p, double *u0, double *v,
                 double *a, double *b, double *c, double *d, double *e, double *f){
	int i, j;
	int R = J - 2; // system size

	// parameters
	double Bo = p[0];		// Bond number
	double U0 = p[1];		// left boundary value
	double U1 = p[2];		// right boundary value
	
	// setup
	double cf1 = dt/dx;
	double cf2 = 1.0/(2.0*Bo*dx*dx*dx);
	double g[R], g1, g2;
	double A, B, C, D, E, F;
	for (i = 0; i < R-1; i++)
		g[i] = (cf1/3.0)*pow((v[i] + v[i+1])/2.0,3);
	g[R-1] = (cf1/3.0)*pow((v[R-1] + U1)/2.0,3);

	// compute coefficients for the particular problem
	/*---------------*/
	j = 1  ; i = j-1;
	g1   = 0.0;
	g2   = cf2*g[i  ];
	
	C    =            4.0*g2 ;
	D    = -          3.0*g2 ;
	E    =                g2 ;
	F    = g[i]              ;

	a[i] = 0.0;
	b[i] = 0.0;
	c[i] = 1.0 + C;
	d[i] = D;
	e[i] = E;
	f[i] = (1.0 - C)*u0[j] - D*u0[j+1] - E*u0[j+2] - F;

	/*---------------*/
	j = 2  ; i = j-1;
	g1 = cf2*g[i-1];
	g2 = cf2*g[i  ];

	B    = -(4.0*g1 +     g2);
	C    =  (3.0*g1 + 3.0*g2);
	D    = -(    g1 + 3.0*g2);
	E    =                g2 ;
	F    = g[i] - g[i-1]     ;

	a[i] = 0.0;
	b[i] = B;
	c[i] = 1.0 + C;
	d[i] = D;
	e[i] = E;
	f[i] = -B*u0[j-1] + (1.0 - C)*u0[j] - D*u0[j+1] - E*u0[j+2] - F;

	/*---------------*/
	for (j = 3; j < J-2; j++){
		i = j-1;
		g1 = cf2*g[i-1];
		g2 = cf2*g[i  ];

		A    =       g1          ;
		B    = -(3.0*g1 +     g2);
		C    =  (3.0*g1 + 3.0*g2);
		D    = -(    g1 + 3.0*g2);
		E    =                g2 ;
		F    = g[i] - g[i-1]     ;

		a[i] = A;
		b[i] = B;
		c[i] = 1.0 + C;
		d[i] = D;
		e[i] = E;
		f[i] = -A*u0[j-2] - B*u0[j-1] + (1.0 - C)*u0[j] - D*u0[j+1] - E*u0[j+2] - F;
	}
	
	/*---------------*/
	j = J-2; i = j-1;
	g1 = cf2*g[i-1];
	g2 = cf2*g[i  ];

	A    =       g1          ;
	B    = -(3.0*g1 +     g2);
	C    =  (3.0*g1 + 3.0*g2);
	D    = -(    g1 + 3.0*g2);
	F    = g[i] - g[i-1] + 2.0*U1;

	a[i] = A;
	b[i] = B;
	c[i] = 1.0 + C;
	d[i] = D;
	e[i] = 0.0;
	f[i] = -A*u0[j-2] - B*u0[j-1] + (1.0 - C)*u0[j] - D*u0[j+1] - F;

	/*---------------*/
	j = J-1; i = j-1;
	g1 = cf2*g[i-1];
	g2 = cf2*g[i  ];

	A    =       g1          ;
	B    = -(3.0*g1 +     g2);
	C    =  (3.0*g1 + 3.0*g2);
	F    = g[i] - g[i-1] + 2.0*U1 - 2.0*(g1 + 2.0*g2);

	a[i] = A;
	b[i] = B;
	c[i] = 1.0 + C;
	d[i] = 0.0;
	e[i] = 0.0;
	f[i] = -A*u0[j-2] - B*u0[j-1] + (1.0 - C)*u0[j] - F;

	// FOR DEBUGGING
	//printf("\n");
	//for (i = 0; i < R; i++){
	////	printf("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", a[i], b[i], c[i], d[i], e[i], f[i]);
	//	printf("%.4f ", c[i]);
	//}
	//printf("\n");
}


#endif
