#include <armadillo>
#include "node.h"
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

void node::setConnections(Pipe pipes[],node nodes[],int N,int M){

	int n = (index-1)%(N);
	int m = floor((index-1)/N);

	leftPipe  = &pipes[(m*((N+M-1))+n)+(n==0)*(N)];
	leftNode  = &nodes[index-1+(n==0)*(N)];
	
	rightPipe = &pipes[m*((N+M-1))+n+1];
	rightNode = &nodes[index+1-(n==(N-1))*(N)];
	
	double above = m*((N+M-1))+N+(double)(n+1+((m+1)%(2)))/2;
	double below = m*((N+M-1))-(double)(M-n+((m+1)%(2)))/2+(m==0)*(N+M-1)*M;

	if (ceil(above) == above){
		horizontalPipe = &pipes[(int)above];
		horizontalNode = &nodes[index+N-(m==(M-1))*N*M];
	}
	else{
		horizontalPipe = &pipes[(int)below];
		horizontalNode = &nodes[index-N+(m==0)*N*M];
	}
}

void node::fillAandB(mat *A, mat *B){
	cout<<"node is "<<index<<endl<<"horizontal is "<<horizontalNode->index<<"lef is "<<leftNode->index<<endl<<"right is " <<rightNode->index<<endl;
	(*A)(index-1,index-1) = -pow(rightPipe->radius,4)-pow(leftPipe->radius,4)-pow(horizontalPipe->radius,4);
	(*A)(index-1,horizontalNode->index-1) = pow(horizontalPipe->radius,4);
	(*A)(index-1,leftNode->index-1) = pow(leftPipe->radius,4);
	(*A)(index-1,rightNode->index-1) = pow(rightPipe->radius,4);
}
