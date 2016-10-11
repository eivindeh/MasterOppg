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

mat setOffDiag(mat M,int N,int l,vec v){
	int N2 = pow(N,2);
	for(int i = 0; i<N-l; i++){
		M(i,i+l) = v(i);
		M(i+l,i) = v(i);
	}
	return M;
}

mat getR1(int N){
	mat R = zeros<mat>(N,N);
	vec v(N-1,fill::randu);
	v = pow(v+0.5,4);
	R = setOffDiag(R,N,1,v);
	R(N-1,0) = pow(((double)rand()/RAND_MAX+0.5),4);
	R(0,N-1) = R(N-1,0);
	//cout<<"R1 = "<<endl<<R<<endl;
	return R;
}

mat getR2(int N){
	mat R = zeros<mat>(N,N);
	vec v(N,fill::randu);
	v = pow(v+0.5,4);
	for (int i = 0;i<N;i = i+2){
		v(i) = 0;
	}
	R = diagmat(v);
	//cout<<"R2 = "<<endl<<R<<endl;
	return R;
}

mat getR3(int N){
	mat R = zeros<mat>(N,N);
	vec v(N,fill::randu);
	v = pow(v+0.5,4);
	for (int i = 1;i<N;i = i+2){
		v(i) = 0;
	}
	R = diagmat(v);
	cout<<"R3 = "<<endl<<R<<endl;
	return R;
}

mat getR(int N,int M) {
	mat RN;
	int K = N*M;
	mat R = zeros<mat>(K,K);
	for (int i = 0; i<M; i++){
		R(span(N*i,i*N+N-1),span(N*i,N*i+N-1)) = getR1(N);
	}
	for (int i = 0; i<M-1; i = i+2){
		R(span(N*i+N,N*i+2*N-1),span(N*i,N*i+N-1)) = getR3(N);
		R(span(N*i,N*i+N-1),span(N*i+N,N*i+2*N-1)) = R(span(N*i+N,N*i+2*N-1),span(N*i,N*i+N-1));
	}
	for (int i = 1; i<M-1; i = i+2){
		R(span(N*i+N,N*i+2*N-1),span(N*i,N*i+N-1)) = getR2(N);
		R(span(N*i,N*i+N-1),span(N*i+N,N*i+2*N-1)) = R(span(N*i+N,N*i+2*N-1),span(N*i,N*i+N-1));
	}
	
	R(span(0,N-1),span(N*M-N,N*M-1)) = getR2(N);
	R(span(N*M-N,N*M-1),span(0,N-1)) = R(span(0,N-1),span(N*M-N,N*M-1));
	
	//v1 = pow(v1+0.5,4);
	//v2 = pow(v2+0.5,4);
	//v3 = pow(v3+0.5,4);
	//v4 = pow(v4+0.5,4);
	
	return R;
}



mat getA(int N,int M, mat R) {
	int N2 = pow(N,2);
	mat A = zeros(N*M+1,N*M);
	vec v = ones<vec>(N*M);
	A(span(0,N*M-1),span(0,N*M-1)) = R-diagmat(R*v);
	A(N*M,0) = 1;
	return A;
}

vec getB(int N,int M,mat R){
	vec B = zeros<vec>(N*M+1);
	for(int i = 1; i<N; i = i+2){
		B(i) = -1*R(i,i+N*M-N);
		B(i+N*M-N) = 1*R(i,i+N*(M-1));
	}
	cout<<B<<endl;
	return B;
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

mat getQ(mat R,vec p, int N,int M){
	mat Q = zeros(N*M,N*M);
	for (int j = 0; j < N*M; j++){
		for (int i = 0; i < N*M; i++){
			if (R(j,i) != 0){
				if(j<N && i>N*(M-1)){
				Q(j,i) = (p(i) - p(j)+1)*R(j,i);
				}
				else if(i<N && j>N*(M-1)){
				Q(j,i) = (p(i) - p(j)-1)*R(j,i);
				}
				else{
				Q(j,i) = (p(i) - p(j))*R(j,i);
				}
			}
		}
	}
	return Q;
}

int main() {
	srand(time(NULL));
	double P0 = 0;
	double PN = 1;
	int N = 6;
	int M = 4;
	arma_rng::set_seed_random();
	//vec r(N,fill::randu);
	//r = pow(r+0.5,4);
	mat R = getR(N,M);
	//cout<<R<<endl;
	R.save("R.txt",raw_ascii);
	mat A = getA(N,M,R);
	//cout << A <<endl;
	//cout << pinv(A) <<endl;
	
	vec B = getB(N,M,R);
	cout << B <<endl;
	vec p = pinv(A)*B;
	
	cout << p <<endl;
	//A.save("A.txt",raw_ascii);
	//B.save("B.txt",raw_ascii);
	//mat pMat = vec2mat(p,N);
	mat qMat = getQ(R,p,N,M);
	//cout<<qMat<<endl;
	//cout<<p<<endl;
	//cout <<pMat<<endl;
	qMat.save("Q.txt",raw_ascii);
	R.save("R.txt",raw_ascii);
	
	return 0;
}
