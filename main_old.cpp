#include "pipe.h"
#include "node.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <armadillo>
#include <time.h>
//#include "draw.cpp"
#define maxBubbles 15
#define ERR 1.0e-14

using namespace std;
using namespace arma;

class pipe;
class node;

int M = 4*3;
int N = 6*3;
int xLen = N;
int yLen = M;
int nNodes = N*M;
int nPipes = N*M*3/2;
double dt = 1;
int k = 0;
double totalT = 0;
mat bubbleStartOut = zeros(nPipes*100,15);
mat bubbleStopOut = zeros(nPipes*100,15);
vec TotFlow = zeros(2000);

node* nodes = new node[nNodes];
Pipe* pipes = new Pipe[nPipes];

void initializeNodesAndPipes(){
	for (int i = 0; i<nPipes; i++){
		pipes[i].index = i;
	}
	for (int i = 0; i<nNodes; i++){
		nodes[i].index = i;
		nodes[i].setConnections(pipes,nodes,N,M);
	}
	for(int i = 0; i<nPipes*1/4; i++){
		pipes[i].bubbleStop[0] = 1.0;
		pipes[i].bubbleStart[0] = 0.0;
		pipes[i].nBubbles = 1;
	}
	for(int i = 0; i<nPipes; i++){
		pipes[i].calcPcandMobility();
		pipes[i].getFlow();
	}
}

/*void getPressure(){
	mat A = zeros(nNodes,nNodes);
	vec B = zeros(nNodes);
	for (int i = 0; i<nNodes; i++){	
		nodes[i].fillAandB(&A,&B);
	}
	vec p = solve(A,B);

	for (int i = 0; i<nNodes; i++){
		nodes[i].pressure = p(i);
	}
}*/

void getFlow(){
	for (int i = 0; i<nPipes; i++){
		pipes[i].getFlow();
	}
}

void updateBubbles(){
		for (int i = 0; i<nNodes; i++){
			nodes[i].nodeFrac = 0;
		}
		for (int i = 0; i<nPipes; i++){
			pipes[i].moveBubble(dt);
		}
		for (int i = 0; i<nPipes; i++){
			pipes[i].killBubbles();
		}
		for (int i = 0; i<nNodes; i++){
			nodes[i].updateNodeFrac();
		}
		for (int i = 0; i<nPipes; i++){
			pipes[i].distributeBubbles();
		}
		for (int i = 0; i<nNodes; i++){
			nodes[i].bubbleAdjust();
		}
}

void setOutput(int i, int K){
	bool swappedd;
	if(i%K==0){
		k++;
		for (int i = 0; i<nPipes;i++){
			//cout<<i<endl;
			swappedd = false;
			if((!pipes[i].swapped && pipes[i].nodeN->horizontalPipe->index != pipes[i].index) || (pipes[i].swapped && pipes[i].nodeN->horizontalPipe->index == pipes[i].index)){
				swappedd = true;
				pipes[i].Swap();
			}
			for(int j = 0; j<15; j++){
				bubbleStartOut(i+N*M/2*3*(k-1),j) = pipes[i].bubbleStart[j];
				bubbleStopOut(i+N*M/2*3*(k-1),j) = pipes[i].bubbleStop[j];
			}
			if(swappedd){
				pipes[i].Swap();
			}
		}
	}
}

void solvePressure(){
	double p[nNodes];
	double r[nNodes];
	double ap[nNodes];
	double rpn = 0.0;
	node c_node;
	Pipe link;
	for (int i = 0; i < nNodes; i++){
		p[i] = 0;
		for(int j = 0; j<3; j++){
		if(j == 0){
        		link = *nodes[i].leftPipe;
        		c_node = *nodes[i].leftNode;}
        	else if (j == 1){
        		link = *nodes[i].rightPipe;
        		c_node = *nodes[i].rightNode;}
        	else{
        		link = *nodes[i].horizontalPipe;
        		c_node = *nodes[i].horizontalNode;}
			
			p[i] -= link.mobility*((c_node.pressure - nodes[i].pressure)+ (-link.Pc + link.boundary)*(nodes[i].upper? 1: -1) );
			}
			//p[i] = (nodes[i].leftPipe->flow + nodes[i].rightPipe->flow + nodes[i].horizontalPipe->flow)*(nodes[i].upper? -1: 1);
		
		r[i] = p[i];
		rpn += r[i]*r[i];
	}
	
	int v = 0;
	while (rpn >= ERR){
		v++;
		double rps = rpn;
		double am = 0.0;
		for(int i = 0; i < nNodes; i++){
			double api = nodes[i].leftPipe->mobility * (p[nodes[i].leftNode->index] - p[i]) +
			             nodes[i].rightPipe->mobility * (p[nodes[i].rightNode->index] - p[i]) +
			             nodes[i].horizontalPipe->mobility * (p[nodes[i].horizontalNode->index] - p[i]);
			ap[i] = api;
			am   += p[i]*api;
		}
		am = rps/am;
		rpn = 0.0;
		for(int i = 0; i<nNodes; i++){
			nodes[i].pressure += am * p[i];
			r[i] -= am * ap[i];
			rpn += r[i] * r[i];
		}
		
		double bm = rpn / rps;
		for (int i = 0; i < nNodes; i++){
			p[i] = r[i] + bm*p[i];
		}
		//cout<<" numerical check: "<<rpn<<endl;
		//for(int i = 0; i < nNodes; i++){
			//cout<<ap[i];
		//}
	}
}

void debug(int i){
	cout<<endl<<" Iteration "<<i<<endl;
	//cout<<"nBubbles "<<pipes[22].nBubbles<<endl;
	//cout<<"unswapped: ";
	/*for(int j = 0; j<pipes[26].nBubbles; j++){
	cout<<" bStart "<<pipes[26].bubbleStart[j]<<" bStop "<<pipes[26].bubbleStop[j]<<endl;
	}
	pipes[26].Swap();*/
	
	for(int j = 0; j<pipes[22].nBubbles; j++){
	cout<<" bStart "<<pipes[22].bubbleStart[j]<<" bStop "<<pipes[22].bubbleStop[j]<<endl;
	}
	cout<<" swapped? "<<pipes[22].swapped<<endl;
	cout<<"flow: "<<pipes[22].flow<<endl;
	/*cout<<"swapping... "<<endl;
	pipes[22].Swap();
	for(int j = 0; j<pipes[22].nBubbles; j++){
	cout<<" bStart "<<pipes[22].bubbleStart[j]<<" bStop "<<pipes[22].bubbleStop[j]<<endl;
	}
	cout<<" swapped? "<<pipes[22].swapped<<endl;
	pipes[22].Swap();*/
	//cout<<"Nodeindex: "<<i<<"frac: "<<nodes[21].nodeFrac<<endl;
	//cout<<"NodeNI "<<pipes[34].nodeN->index<<" NodePI " <<pipes[34].nodeP->index<<endl;
	//cout<<"NodeN: "<<pipes[34].nodeN->pressure<<" NodeP: "<<pipes[34].nodeP->pressure<< " Flow "<<pipes[34].flow<< " Swapped? "<<pipes[34].swapped<<endl; 
	for (int j = 0; j<nNodes; j++){
		/*cout<<" Flow: "<<nodes[j].rightPipe->flow<<" Pthis: "<<nodes[j].pressure<<" P horz "<<nodes[j].rightNode->pressure;
		cout<<" NodeFrac "<<nodes[j].nodeFrac;
		cout<<" Pc: "<<nodes[j].rightPipe->Pc<<" Mob "<<nodes[j].rightPipe->mobility;
		cout<<"Pressure: "<<nodes[j].pressure<<endl;
		cout<<" swap? "<< nodes[j].rightPipe->swapped;
		*/
		//cout<<" Pressure: "<<nodes[j].pressure<<"   totalFlow "<<nodes[j].leftPipe->flow+nodes[j].rightPipe->flow+nodes[j].horizontalPipe->flow<<endl;
	}
	for(int j = 0; j<nPipes;j++){
		//cout<<" Pipeindex "<<pipes[j].index<<endl;
		//cout<<" bubleStart "<<pipes[j].bubbleStart[0];
		//cout<<" bubbleStop " <<pipes[j].bubbleStop[0];
		//cout<<endl<<" AltNodeIndex "<< pipes[22].nodeP->index<<endl;
	}
}

/*void timestep(){
	totalT += dt;
	
	for(int i = 0; i<nPipes; i++){
		pipes[i].calcPcandMobility();
	}
	solvePressure();
	getFlow();
	//debug(i);
	updateBubbles();
}*/

void meashureFlow(int i){
	double totalFlow = 0;
	for (int i = N; i<3*N/2; i++){
		totalFlow += (double)pipes[i].lenTL;
	}
	TotFlow(i) = totalFlow;
}
int main(){
	mat A = zeros(nNodes,nNodes);
	vec B = zeros(nNodes);
	initializeNodesAndPipes();	
	
	cout<<endl;
	for (int i = 0; i<2000; i++){
		for(int i = 0; i<nPipes; i++){
			pipes[i].calcPcandMobility();
		}
		solvePressure();
		getFlow();
		//debug(i);
		updateBubbles();
		
		meashureFlow(i);
		//setOutput(i,10);	
		cout<<"iteration: "<<i<<endl;
	}
	
	mat Output = zeros(5,N*M/2*3);
	for (int i = 0; i<N*M/2*3; i++){
		Output(0,i) = pipes[i].flow;
		Output(1,i) = pipes[i].x1;
		Output(2,i) = pipes[i].y1;
		Output(3,i) = pipes[i].x2;
		Output(4,i) = pipes[i].y2;
	}
	cout<<"her?"<<endl;
	/*bubbleStartOut.save("bStart.txt",raw_ascii);
	bubbleStopOut.save("bStop.txt",raw_ascii);
	Output.save("Output.txt",raw_ascii);*/
	TotFlow.save("TotalFlow.txt",raw_ascii);
	//cout<<"herda? "<<endl;
	delete[] nodes;
	delete[] pipes;
return 0;
}
