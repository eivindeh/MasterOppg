#pragma once
#include "pipe.h"
#include <armadillo>
using namespace std;
using namespace arma;

class Pipe;

class node{
	public:
	int index;
	double pressure;
	bool upper;
	int n;
	int m;
	double x;
	double y;
	Pipe* leftPipe;
	Pipe* rightPipe;
	Pipe* horizontalPipe;
	node* leftNode;
	node* rightNode;
	node* horizontalNode;
	double nodeFluxTL;
	double nodeFluxNW;
	double nodeFrac;
	void setConnections(Pipe pipes[],node nodes[],int N,int M);
	void fillAandB(mat *A, vec *B);
	void getFlow();
	void updateNodeFrac();
	void bubbleAdjust(Pipe* pipes);
};
