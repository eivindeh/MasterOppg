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
	for(int i = 0; i<N2-l; i++){
		M(i,i+l) = v(i);
		M(i+l,i) = v(i);
	}
	return M;
}


mat getR(int N) {
	int k = 0;
	int N2 = pow(N,2);
	mat R = zeros<mat>(pow(N,2),pow(N,2));
	vec v1(N2-1,fill::randu); //horizontal connections
	vec v2(N2-N+1,fill::randu); // horizontal boundery connections
	vec v3(N2-N,fill::randu); // Vertical connections
	vec v4(N,fill::randu); // vertical boundary connections
	
	v1 = pow(v1+0.5,4);
	v2 = pow(v2+0.5,4);
	v3 = pow(v3+0.5,4);
	v4 = pow(v4+0.5,4);
	//set horizontal elements not connected to zero
	for (int i = 0; i<N2-1; i++){
		if( !((i+1)%N)){
			v1(i) = 0;
		}
	}
	
	for (int i = 0; i<N2-N; i++){
		if(i%N){
			v2(i) = 0;
		}
	}
	cout <<v1<<endl<<v2<<endl<<v3<<endl<<v4<<endl;
	R = setOffDiag(R,N,1,v1);
	R = setOffDiag(R,N,N-1,v2);
	R = setOffDiag(R,N,N,v3);
	R = setOffDiag(R,N,N2-N,v4);
	return R;
	
}



mat getA(int N, mat R) {
	int N2 = pow(N,2);
	vec v = ones<vec>(N2);
	R = R-diagmat(R*v);
	mat A = zeros<mat>(N2+2*N+1,N2+2*N);
	
	A(span(0,N2-1),span(N,N2+N-1)) = R;
	A(span(0,N-1),span(N2+N,N2+2*N-1)) = A(span(0,N-1),span(N2,N2+N-1));
	A(span(N2-N,N2-1),span(0,N-1)) = A(span(N2-N,N2-1),span(N,2*N-1)); 
	A(span(0,N-1),span(N2,N2+N-1)) = zeros(N,N);
	A(span(N2-N,N2-1),span(N,2*N-1)) = zeros(N,N);
	
	A(span(N2,N2+N-1),span(N,2*N-1))= eye(N,N);
	A(span(N2,N2+N-1),span(N2,N2+N-1)) = -eye(N,N);
	A(span(N2+N,N2+2*N-1),span(0,N-1)) = eye(N,N);
	A(span(N2+N,N2+2*N-1),span(N,2*N-1)) = -eye(N,N);
	A(span(N2+N,N2+2*N-1),span(N2,N2+N-1)) = -eye(N,N);
	A(span(N2+N,N2+2*N-1),span(N2+N,N2+2*N-1)) = eye(N,N);
	A(N2+2*N,N2) = 1;
	return A;
}

vec getB(int N){
	vec B = zeros<vec>(pow(N,2)+2*N+1);
	B(span(pow(N,2),pow(N,2)+N-1)).ones();
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

mat getQ(mat R,vec p, int N){
	int N2 = pow(N,2);
	mat Q = zeros(N2,N2);
	for (int j = 0; j < N2; j++){
		for (int i = 0; i < N2; i++){
			if (R(j,i) != 0){
				if(i < N && j >= pow(N,2)-N){ 
					Q(j,i) = (p(i+N) - p(2*N+j))*R(j,i);
				}
				else if (j < N && i >= pow(N,2)-N){
					Q(j,i) = (p(i+N)-p(j))*R(j,i);
					//cout<<"i = " << i <<"   j = "<< j <<endl;
					//cout <<"p(i) ="<<p(i+N)<<"   p(j) = " <<p(j)<<endl;
				}
				else{
					Q(j,i) = (p(i+N)-p(j+N))*R(j,i);
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
	int N = 15;
	int N2 = pow(N,2);
	//arma_rng::set_seed_random();
	//vec r(N,fill::randu);
	//r = pow(r+0.5,4);
	mat R = getR(N);
	//cout<<R<<endl;
	
	mat A = getA(N,R);
	//cout << A <<endl;
	
	vec B = getB(N);
	//cout << B <<endl;
	
	vec p_full = solve(A,B);
	
	vec p = p_full(span(N,N2+N-1));
	
	cout << p <<endl;
	//A.save("A.txt",raw_ascii);
	//B.save("B.txt",raw_ascii);
	mat pMat = vec2mat(p,N);
	mat qMat = getQ(R, p_full,N);
	cout<<qMat<<endl;
	cout<<p_full<<endl;
	//cout <<pMat<<endl;
	qMat.save("Q.txt",raw_ascii);
	R.save("R.txt",raw_ascii);
	cout<<R;
	
	return 0;
}
