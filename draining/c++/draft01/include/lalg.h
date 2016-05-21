/* LINEAR ALGEBRA
 *  Solve linear system A*u = b.
 *
 * REFERENCES
 *  von Rosenberg, Methods for the Numerical Solution of PDEs (1969) - Appendix
 *  
 * PARAMETERS
 *  R  						[input]   system size
 *  a,b,c,d,e,f		[input]		coefficients
 *  u  						[output]	solution vector
 */

#ifndef LALG_H
#define LALG_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>

/* PROTOTYPES */
void lalg_tri (int, double *, double *, double *, double *, double *);
void lalg_pent(int, double *, double *, double *, double *, double *, double *, double *);

/* IMPLEMENTATIONS */

/* tridiagonal algorithm
 *  a(i)*u(i-1) + b(i)*u(i) + c(i)*u(i+1) = d(i)
 *  for  i = 0, 1, ..., R-1
 *  with a(0) = c(R-1) = 0
 */
void lalg_tri(int R, double *a, double *b, double *c, double *d, double *u){
	int i;
	double beta[R], gamm[R];
	
	// compute coefficients
	beta[0] = b[0]     ;
	gamm[0] = d[0]/b[0];

	for (i = 1; i < R; i++){
		beta[i] =  b[i-1] - a[i]*c   [i-1] /beta[i-1];
		gamm[i] = (d[i  ] - a[i]*gamm[i-1])/beta[i  ];
	}

	// compute solution vector
	u[R-1] = gamm[R-1];

	for (i = R-2; i > -1; i--){
		u[i] = gamm[i] - c[i]*u[i+1]/beta[i];
	}
}

/* pentadiagonal algorithm
 *  a(i)*u(i-2) + b(i)*u(i-1) + c(i)*u(i) + d(i)*u(i+1) + e(i)*u(i+2) = f(i)
 *  for  i = 0, 1, ..., R-1
 *  with a(0) = b(0) = a(1) = e(R-2) = d(R-1) = e(R-1) = 0
 */
void lalg_pent(int R, double *a, double *b, double *c, double *d, double *e, double *f, double *u){
	int i;
	double delt[R], lamb[R], gamm[R];
	double beta, mu;
	/* NOTE: beta and mu are used only to compute delt, lamb, and gamm
	 *       and need not be stored after they are computed. The delt, 
	 *       lamb, and gamm must be stored, as they are used in the
	 *       back solution */

	// compute coefficients
	delt[0] = d[0]/c[0];
	lamb[0] = e[0]/c[0];
	gamm[0] = f[0]/c[0];

	mu      =  c[1] - b[1]*delt[1]    ;
	delt[1] = (d[1] - b[1]*lamb[0])/mu;
	lamb[1] =  e[1]                /mu;
	gamm[1] = (f[1] - b[1]*gamm[0])/mu;

	for (i = 2; i < R-2; i++){
		beta    =  b[i] - a[i]*delt[i-2]                     ;
		mu      =  c[i] - beta*delt[i-1] - a[i]*lamb[i-2]    ;
		delt[i] = (d[i] - beta*lamb[i-1]                 )/mu;
		lamb[i] =  e[i]                                   /mu;
		gamm[i] = (f[i] - beta*gamm[i-1] - a[i]*gamm[i-2])/mu;
	}

	beta      =  b[R-2] - a[R-2]*delt[R-4]                       ;
	mu        =  c[R-2] - beta  *delt[R-3] - a[R-2]*lamb[R-4]    ;
	delt[R-2] = (d[R-2] - beta  *lamb[R-3]                   )/mu;
	lamb[R-2] =  e[R-2]                                       /mu;
	gamm[R-2] = (f[R-2] - beta  *gamm[R-3] - a[R-2]*gamm[R-4])/mu;

	beta      =  b[R-1] - a[R-1]*delt[R-3]                       ;
	mu        =  c[R-1] - beta  *delt[R-2] - a[R-1]*lamb[R-3]    ;
	delt[R-1] =  0.0                                             ;
	lamb[R-1] =  0.0                                             ;
	gamm[R-1] = (f[R-1] - beta  *gamm[R-2] - a[R-1]*gamm[R-3])/mu;

	// compute solution vector
	u[R-1] = gamm[R-1]                   ;
	u[R-2] = gamm[R-2] - delt[R-2]*u[R-1];

	for (i = R-3; i > -1; i--){
		u[i] = gamm[i] - delt[i]*u[i+1] - lamb[i]*u[i+2];
	}
}

void lalg_bitri(){
}

void lalg_tritri(){
}

#endif
