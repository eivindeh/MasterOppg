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

mat getR(int N) {
	int k = 0;
	mat R = zeros<mat>(pow(N,2),pow(N,2));
	for (int j = 0; j < N-1; j++){
		for (int i = 0; i < N; i++){
			k = i+(j)*N;
			cout<<k<<endl;
			R(k+1,k) = pow((double)rand()/RAND_MAX+0.5,4);
			R(k,k+1) = R(k+1,k);
			R(k+N,k) = pow((double)rand()/RAND_MAX+0.5,4);
			R(k,k+N) = R(k+N,k);
		}
		R(k+1,k) = 0;
		R(k,k+1) = 0;
		R(k,k-N+1) = pow((double)rand()/RAND_MAX+0.5,4);
		R(k-N+1,k) = R(k,k-N+1);
		//R(k+1,k) = R(k-(N-1)+1,k-(N-1));
		//R(k,k+1) = R(k+1,k);
	}
	for (int i = 1; i<N; i++){
		R(k+i,k+i+1) = pow((double)rand()/RAND_MAX+0.5,4);
		R(k+i+1,k+i) = R(k+i,k+i+1);
	}
	R(pow(N,2)-1,pow(N,2)-3) =  pow((double)rand()/RAND_MAX+0.5,4);
	R(pow(N,2)-3,pow(N,2)-1) = R(pow(N,2)-1,pow(N,2)-3);
	return R;
	
}



mat getA(int N, mat R) {
	vec v = ones<vec>(pow(N,2));
	R = R-diagmat(R*v);
	v(span(N,pow(N,2)-N-1)) = zeros<vec>(pow(N,2)-2*N);
	//cout<<v<<endl;
	return R-diagmat(v);
}

vec getB(int N){
	vec B = zeros<vec>(pow(N,2));
	B(span(0,N-1)).ones();
	return -B;
}

mat vec2mat(vec p, int N){
	mat pMat = zeros<mat>(N,N);
	for (int j = 0; j < N; j++){
		for (int i = 0; i < N; i++){
			pMat(j,i) = p(i+j*N);
		}
	}
	return pMat;
}

int main() {
	srand(time(NULL));
	double P0 = 0;
	double PN = 1;
	int N = 20;
	//arma_rng::set_seed_random();
	//vec r(N,fill::randu);
	//r = pow(r+0.5,4);
	mat R = getR(N);
	//cout<<R<<endl;
	
	mat A = getA(N,R);
	//cout << A <<endl;
	
	vec B = getB(N);
	//cout << B <<endl;
	
	vec p = solve(A,B);
	//cout << p <<endl;
	
	mat pMat = vec2mat(p,N);
	//cout <<pMat<<endl;
	pMat.save("p.txt",raw_ascii);
	
	return 0;
}
