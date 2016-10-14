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

int M = 20;
int N = 20;
int yLen = M;
int nNodes = N*M;
int nPipes = N*M*3/2;
int k = 0;

double Pi = 4.0*atan(1.0);
double stens = 3.0;         // surface tension: units = dyn/mm == 10 mN/m
double linkRadMin = 0.1;
double linkRadMax = 0.4;
double dt = 1.0;
double saturationNW = 0.5;
double totalT = 0;
double pressure = 0.7 * yLen * 2 * stens / ((linkRadMin + linkRadMax)*0.5);//*((sin(Pi / 6.0) + 1) / cos(Pi / 6.0);
double poreVolume = N*M*3.0/2.0*0.25*0.25*3.1415;
double saturation;

int iterations = 10000;
int numPoints = 21;

mat bubbleStartOut = zeros(nPipes*100,15);
mat bubbleStopOut = zeros(nPipes*100,15);
mat Ainv = zeros(nNodes,nNodes);
mat A = zeros(nNodes,nNodes);
vec B = zeros(nNodes);
vec TotFlow = zeros(iterations);
vec Flow = zeros(iterations*numPoints);
vec avgFlow =zeros(numPoints);
vec avgFrac = zeros(numPoints);
vec FracFlow = zeros(iterations);
vec TotFracFlow = zeros(iterations*numPoints);
node* nodes = new  node[nNodes]; 
Pipe* pipes = new Pipe[nPipes];

void initializeNodesAndPipes();
void getFlow();
void updateBubbles();
void setOutput(int i, int K);
void solvePressure();
void debug(int i);
void meashureFlow(int i);
double calcDeltaT();
double totSaturation();


void initializeNodesAndPipes(){
	//nodes = new node[nNodes];
	//pipes = new Pipe[nPipes];
	for (int i = 0; i<nPipes; i++){
		pipes[i].index = i;
		pipes[i].flow = 1;
		pipes[i].nBubbles = 0;
		pipes[i].swapped = false;
	}
	for (int i = 0; i<nNodes; i++){
		nodes[i].index = i;
		nodes[i].setConnections(pipes,nodes,N,M);
		nodes[i].pressure = 0;
	}
	cout<<"nubpipes "<<(int)round(nPipes*saturation)<<endl;
	for(int i = 0; i<(int)(nPipes*saturation); i++){
		pipes[i].bubbleStop[0] = 1.0;
		pipes[i].bubbleStart[0] = 0.0;
		pipes[i].nBubbles = 1;
	}
	for(int i = 0; i<nPipes; i++){
		pipes[i].calcPcandMobility();
		pipes[i].getFlow(pressure);
	}
}

void getPoreVolume(){
	poreVolume = 0;
	for(int i = 0; i<nPipes; i++){
		poreVolume += pipes[i].area;
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
        if (saturation == 0.0 || saturation == 1.0) {
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
	return totSat;
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
		//cout<<"satDiff "<<curSat-prevSat<<endl;
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

void ldlt(){
	//A = zeros(nNodes,nNodes);
	A(0,0) = 1;
	for(int i = 1; i<nNodes; i++){
		A(i,i) 				= - nodes[i].rightPipe->mobility - nodes[i].leftPipe->mobility - nodes[i].horizontalPipe->mobility;
		A(i,nodes[i].horizontalNode->index) 	= nodes[i].horizontalPipe->mobility;
		A(i,nodes[i].leftNode->index) 		=  nodes[i].leftPipe->mobility;
		A(i,nodes[i].rightNode->index) 		=  nodes[i].rightPipe->mobility;
	}
	Ainv = inv(A);
	//cout<<A<<endl<<Ainv<<endl;
}

void ldlt_solve(){
	int n,m;
	double b[nNodes];
	B(0) = 0;
	for(int i = 1; i<nNodes; i++){
		B(i) = 0;
		for(int j = 0; j < 3; j++){
			if(j == 0){
				m = nodes[i].leftPipe->index;
			}else if (j == 1){
				m = nodes[i].rightPipe->index;
			}
			else{
				m = nodes[i].horizontalPipe->index;
			}
			
			B(i) -= pipes[m].mobility*(-pipes[m].Pc+pressure*pipes[m].boundary)*(nodes[i].upper? 1: -1);//pluss/minus?
			//cout<<"mob: "<<pipes[m].mobility<<endl;
		}
	}

	vec P = Ainv*B;
	//cout<<"P: "<<P<<endl<<endl;
	//cout<<"B"<<endl<<B<<endl<<endl;
	//cout<<solve(A,B)<<endl;
	for(int i = 0; i<nNodes; i++){
		nodes[i].pressure = P(i);
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
	//cout<<"pressure 20: "<<nodes[20].pressure;
	//cout<<"   flow 30: " <<pipes[30].flow<<endl;
	
}

void meashureFlow(int i){
	double totalFlow = 0.0;
	double totalFlowNW = 0.0;
	for (int j = N; j<3*N/2; j++){
		totalFlow += pipes[j].lenTL*pipes[j].area*(pipes[j].swapped?-1:1);
		totalFlowNW += pipes[j].lenNW*pipes[j].area*(pipes[j].swapped?-1:1);
	}
	TotFlow(i) = totalFlow/dt;
	FracFlow(i) = totalFlowNW/totalFlow;
	if(totalFlowNW>totalFlow+10e-12){
		/*cout<<"KRISEE! "<<"totFlow "<<totalFlow<<"NWflow "<<totalFlowNW<<endl;
		for(int j = N; j<3*N/2; j++){
			if(pipes[j].lenTL*pipes[j].area*(pipes[j].swapped?-1:1) < pipes[j].lenNW*pipes[j].area*(pipes[j].swapped?-1:1)-10e-12){
			cout<<pipes[j].lenTL*pipes[j].area*(pipes[j].swapped?-1:1)<<"  "<<pipes[j].lenNW*pipes[j].area*(pipes[j].swapped?-1:1)<<endl;
			}
		}*/
	}
}



int main(){
	int iter = 0;
	for(saturation = 0.0; saturation <= 1.01; saturation += 1/((double)numPoints-1)){
	cout<<"iteration: " << iter<<endl;
	cout<<"saturation: "<<saturation<<endl;
	initializeNodesAndPipes();	
	getPoreVolume();
	cout<<"real Saturation: "<<totSaturation()/poreVolume<<endl;
	//debug(i);
	cout<<endl;
	double TotQ = 0;
	int i = 0;
	ldlt();
	while (i < iterations){
		for(int j = 0; j<nPipes; j++){
			pipes[j].calcPcandMobility();
		}
		
		ldlt_solve();
		//solvePressure();
		//debug(i);
		getFlow();
		updateBubbles();
			
		meashureFlow(i);
		//setOutput(i,10);	
		TotQ += TotFlow(i);
		i++;
		//cout<<"percent done: "<<(TotQ/2.0/poreVolume)*100.0<<endl;
		//cout<<"totalSaturation: "<<totSaturation()/poreVolume<<endl;
	}
	
	/*mat Output = zeros(5,N*M/2*3);
	for (int i = 0; i<N*M/2*3; i++){
		Output(0,i) = pipes[i].flow;
		Output(1,i) = pipes[i].x1;
		Output(2,i) = pipes[i].y1;
		Output(3,i) = pipes[i].x2;
		Output(4,i) = pipes[i].y2;
	}
	bubbleStartOut.save("bStart.txt",raw_ascii);
	bubbleStopOut.save("bStop.txt",raw_ascii);*/
	//Output.save("Output.txt",raw_ascii);
	//TotFlow.save("TotalFlow.txt",raw_ascii);
	for(int i = iterations/2; i<iterations; i++){
		avgFlow(iter) += TotFlow(i);
		avgFrac(iter) += FracFlow(i);
	}
	avgFlow(iter) = avgFlow(iter)/1000.0/((double)N/2);
	avgFrac(iter) = avgFrac(iter)/1000.0;
	cout<<"avgFlow: " <<avgFlow<<" avgFrac: "<<avgFrac<<endl;
	/*for( int i = 0; i<nPipes; i++){
		delete[] pipes[i].bubbleStart;
		delete[] pipes[i].bubbleStop;
	}*/
	Flow(span((iter)*iterations,(iter+1)*iterations-1)) = TotFlow;
	TotFracFlow(span((iter)*iterations,(iter+1)*iterations-1)) = FracFlow;
	iter++;
	//delete[] nodes;
	//delete[] pipes;
	}
	
	avgFlow.save("AvgFlow.txt",raw_ascii);
	Flow.save("Flow.txt",raw_ascii);
	avgFrac.save("AvgFrac.txt",raw_ascii);
	TotFracFlow.save("TotFrac.txt",raw_ascii);
	//cout<<"herda? "<<endl;
return 0;
}
