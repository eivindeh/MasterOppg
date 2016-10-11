#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;


mat getA(int N, vec r) {
	mat A = zeros<mat>(N,N);
	for (int i = 1; i < N-1; i++){
		A(i,i-1) = r(i-1);
		A(i,i)   = -r(i-1)-r(i);
		A(i,i+1) = r(i);
	}
	A(0,0) = -r(0)-r(1);
	A(0,1) = r(1);
	A(N-1,N-1) = -r(N-1)-r(N-2);
	A(N-1,N-2) = r(N-2);
	return A;
}

vec getB(int N, double P0, double PN, vec r){
	vec B = zeros<vec>(N);
	B(0) = -P0*r(0);
	B(N-1) = -PN*r(N-1);
	return B;
}

int main() {
	double P0 = 0;
	double PN = 1;
	int N = 9;
	arma_rng::set_seed_random();
	vec r(N,fill::randu);
	r = pow(r+0.5,4);
	mat A = getA(N,r);
	cout << A <<endl;
	
	vec B = getB(N,P0,PN,r);
	//cout << B <<endl;
	
	vec p = solve(A,B);
	//cout << p <<endl;
	
	//p.save("p.txt",raw_ascii);
	
	return 0;
}
