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
double Pi = 4.0*atan(1.0);
int M = 40;
int N = 40;
int yLen = M;
double stens = 3.0;         // surface tension: units = dyn/mm == 10 mN/m
double linkRadMin = 0.1;
double linkRadMax = 0.4;
int nNodes = N*M;
int nPipes = N*M*3/2;
double dt = 1;
int k = 0;
double saturationNW = 0.5;
double totalT = 0;
mat bubbleStartOut = zeros(nPipes*104,15);
mat bubbleStopOut = zeros(nPipes*104,15);
vec TotFlow = zeros(1000);
vec Flow = zeros(1000*10);
vec avgFlow =zeros(10);
double pressure = 1 * ((sin(Pi / 6) + 1) / cos(Pi / 6)*yLen * 2 * stens / (linkRadMin + linkRadMax)*0.5);
double poreVolume = N*M*3/2*0.25*0.25*3.1415;
double saturation;

node* nodes; 
Pipe* pipes;

void initializeNodesAndPipes(){
	nodes = new node[nNodes];
	pipes = new Pipe[nPipes];
	for (int i = 0; i<nPipes; i++){
		pipes[i].index = i;
	}
	for (int i = 0; i<nNodes; i++){
		nodes[i].index = i;
		nodes[i].setConnections(pipes,nodes,N,M);
	}
	for(int i = 0; i<nPipes*saturation; i++){
		pipes[i].bubbleStop[0] = 1.0;
		pipes[i].bubbleStart[0] = 0.0;
		pipes[i].nBubbles = 1;
	}
	
	for(int i = 0; i<nPipes; i++){
		pipes[i].calcPcandMobility();
		pipes[i].getFlow(pressure);
	}
}

void getFlow(){
	for (int i = 0; i<nPipes; i++){
		pipes[i].getFlow(pressure);
	}
}


double calcDeltaT() {
    double velMax = 0.0;
    int velMaxPos = -1;
    for (int i = 0; i < nPipes; i++) {
        double vel = pipes[i].flow / pipes[i].area; // volume-normalized velocity
	cout<<"vel: "<<vel<<endl;
        if (saturationNW == 0.0 || saturationNW == 1.0) {
            if (velMax < fabs(vel)) {
                velMax = fabs(vel);
                velMaxPos = i;
            }
        } else {
            if (velMax < fabs(vel) && pipes[i].nBubbles > 0) {
                velMax = fabs(vel);
                velMaxPos = i;
            }
        }
    }
    if (velMaxPos == -1) {
        printf("No maximum velocity found.\n"); fflush(stdout);
        exit(1);
    }
    return 0.10 / velMax * 1.0;
}

double totSaturation(){
	double totSat = 0;
	for (int i = 0; i < nPipes; i++){
		totSat += pipes[i].calcLinkSaturation()*pipes[i].area;
	}
}

void updateBubbles(){
	double prevSat;
	double curSat;
		for (int i = 0; i<nNodes; i++){
			nodes[i].nodeFrac = 0;
		}
		dt = calcDeltaT();
		prevSat = totSaturation();
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
		curSat = totSaturation();
		cout<<"satDiff "<<curSat-prevSat<<endl;
		for (int i = 0; i<nNodes; i++){
			nodes[i].bubbleAdjust(pipes);
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
			
			p[i] -= link.mobility*((c_node.pressure - nodes[i].pressure)+ (-link.Pc + pressure*link.boundary)*(nodes[i].upper? 1: -1) );
			}
		
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
	}
}

void debug(int i){
	cout<<endl<<" Iteration "<<i<<endl;
	int TotBubbles = 0;
	for (int i = 0; i<nPipes; i++){
		TotBubbles += pipes[i].nBubbles;
	}
	cout<<"Total Bubbles: "<<TotBubbles<<endl;
	
	for(int j = 0; j<pipes[20].nBubbles;j++){
		cout<<" bubleStart "<<pipes[20].bubbleStart[j];
		cout<<" bubbleStop " <<pipes[20].bubbleStop[j];
		cout<<endl;
	}
	
}

void meashureFlow(int i){
	double totalFlow = 0;
	double totalFlowNW = 0;
	for (int j = N; j<3*N/2; j++){
		totalFlow += pipes[j].lenTL*pipes[j].area*(pipes[j].swapped?-1:1);
		totalFlowNW += pipes[j].lenNW*pipes[j].area;
	}
	TotFlow(i) = totalFlow;
}

int main(){
	int iter = 0;
	for(saturation = 0.5; saturation < 1; saturation += 0.1){
	initializeNodesAndPipes();	
	
	cout<<endl;
	double TotQ = 0;
	int i = 0;
	while (i < 1000){
		for(int j = 0; j<nPipes; j++){
			pipes[j].calcPcandMobility();
		}
		solvePressure();
		getFlow();
		//debug(i);
		updateBubbles();
			
		meashureFlow(i);
//		setOutput(i,10);	
		TotQ += TotFlow(i);
		i++;
		cout<<"percent done: "<<(TotQ/2/poreVolume)*100.0<<endl;
		cout<<"totalSaturation: "<<totSaturation()/poreVolume<<endl;
	}
	
	/*mat Output = zeros(5,N*M/2*3);
	for (int i = 0; i<N*M/2*3; i++){
		Output(0,i) = pipes[i].flow;
		Output(1,i) = pipes[i].x1;
		Output(2,i) = pipes[i].y1;
		Output(3,i) = pipes[i].x2;
		Output(4,i) = pipes[i].y2;
	}*/
	bubbleStartOut.save("bStart.txt",raw_ascii);
	bubbleStopOut.save("bStop.txt",raw_ascii);
	//Output.save("Output.txt",raw_ascii);
	TotFlow.save("TotalFlow.txt",raw_ascii);
	for(int i = 500; i<1000; i++){
		avgFlow(iter) += TotFlow(i);
	}
	avgFlow(iter) = avgFlow(iter)/1000.0;
	/*for( int i = 0; i<nPipes; i++){
		delete[] pipes[i].bubbleStart;
		delete[] pipes[i].bubbleStop;
	}*/
	Flow(span((iter)*1000,(iter+1)*1000-1)) = TotFlow;
	iter++;
	cout<<iter<<endl;
	delete[] nodes;
	delete[] pipes;
	}
	
	avgFlow.save("AvgFlow.txt",raw_ascii);
	Flow.save("Flow.txt",raw_ascii);
	//cout<<"herda? "<<endl;
return 0;
}
