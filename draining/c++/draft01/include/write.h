/* WRITE FILE
 *  Write data file containing the abscissa and solution vectors.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 */

#ifndef WRITE_H
#define WRITE_H

/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/* PROTOTYPES */
void writeFile(string, string, int, int, int, double, double, double *);
void write    (string, string, int, int, int, double, double, double *);

/* IMPLEMENTATIONS */

// Write output to file.
void writeFile(string dir, string fn, int ctr, int n, int J, double dt, double dx, double *uw){
	int j;
	int width = 12;
	int tacc = 4;
	int xacc = 4;
	int uacc = 4;
	int fill = 6;
	string fmt = ".txt";
	string col1 = "x";
	string col2 = "u";
	ofstream file;
	ostringstream time, ts, params, header;
	string line, path;

	// setup: assign strings and stringstreams
	time << fixed << setprecision(tacc) << double(n)*dt;
	ts << setfill('0') << setw(fill) << ctr;
	params << "t = " << time.str();
	cout << params.str() << endl;
	header << setw(width) << col1 << setw(width) << col2;
	path = dir + "/" + fn + ts.str() + fmt;

	// write to file
	file.open(path.c_str());
	file << params.str() << endl;
	file << header.str() << endl;
	for (j = 0; j < J+1; j++){
		ostringstream xx, uu;
		xx << fixed << setprecision(xacc) << setw(width) << (double(j) - 0.5)*dx;
		uu << fixed << setprecision(uacc) << setw(width) << uw[j];
		line = xx.str() + uu.str();

		file << line << endl;
	}
	file.close();

	cout << "Output written to " << path << "." << endl;
}

void write(string dir, string fn, int nw, int N, int J, double dt, double dx, double *u){
	int n, j;
	double uw[J+2];
	int m = N/nw;
	int ctr = 0;

	for (n = 0; n < N+1; n++){
		if (n % m == 0){
			ctr++;
			for (j = 0; j < J+2; j++)
				uw[j] = u[n*(J+2) + j];
			writeFile(dir, fn, ctr, n, J, dt, dx, uw);
		}
		else
			continue;
	}
}


#endif
