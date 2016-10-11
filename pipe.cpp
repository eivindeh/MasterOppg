#include "pipe.h"
#include <armadillo>
#define maxBubbles 15
#define SIGN (x < 0) ? -1 : (x > 0)
#define clamp(x,a,b) (x < a ? a : x > b ? b : x)
#define unlerp(t,a,b)  ( ((t) - (a)) / (double) ((b) - (a)) ) // inverse linear interpolation
#define TOL 1.0e-12
#define maxBubbles 15
using namespace std;
using namespace arma;

double PI = 4.0*atan(1.0);
double sTens = 3.0;         // surface tension: units = dyn/mm == 10 mN/m
double muNON = 0.01;        // viscosity: units = Po == 10 Pa*s
double muWET = 0.01;

Pipe::Pipe(){
	bubbleStart = new double[maxBubbles];
	bubbleStop  = new double[maxBubbles];
	for ( int i = 0; i<maxBubbles; i++){
		bubbleStart[i] = 0;
		bubbleStop[i] = 0;
	}
	boundary = 0;
	radius = ((double)rand()/RAND_MAX*3.0+1.0)*0.1;
	area = PI*pow(radius,2);
	//mobility = pow(radius,4);
}

void Pipe::setCoordinates(double x, double y){
	if (x1 ==-1){
		x1 = x;
	}else{
		x2 = x;
	}
	if(y1 == -1){
		y1 = y;
	}else{
		y2 = y;
	}
}

void Pipe::moveBubble(double dt){
	//cout<<"firstMove "<<nodeP->index<<endl;
	lenTL = abs(flow)/area*dt;
	lenNW = 0.0;
	for (int bubble = 0; bubble< nBubbles; bubble++){
		bubbleStart[bubble] += lenTL;
		bubbleStop[bubble] += lenTL; 
	}
	//cout<<"secondMove "<<nodeP->index<<endl;
	double start 	= 1.0 - lenTL;
        double end 	= 1.0;
        
        for (int i = 0; i < nBubbles; i++) {
            lenNW += clamp(bubbleStop[i], start, end) - clamp(bubbleStart[i], start, end);
        }
        //cout<<" Pipe "<<index<<endl;
        //cout<< " NodeP "<<nodeP->nodeFluxTL<<endl; 
        //cout<<" pIndex "<<" nIndex "<<nodeP->index<<endl;
        nodeP->nodeFluxTL += lenTL*area;
        nodeP->nodeFluxNW += lenNW*area;	
}

void Pipe::distributeBubbles(){
	double dx1 = bubbleStart[0] - lenTL;                  // size of first W  bubble
        double dx2 = bubbleStop[0]  - bubbleStart[0];    // size of first NW bubble
     	if(nBubbles > 0){
     		if(dx1 < TOL){
     			if(dx2 < radius){
     				bubbleStart[0] -= lenTL;
     				bubbleStop[0]  -= (1.0 - nodeN->nodeFrac)*lenTL;
     			}
     			else{
     				bubbleStart[0] -= nodeN->nodeFrac*lenTL;
     			}
     		}
     		else{
     			if(dx1 < radius){
     				bubbleStart[0] -= nodeN->nodeFrac*lenTL;
     			}
     			else{
     				for(int i = maxBubbles -1 ; i > 0; i--){
     					bubbleStart[i] = bubbleStart[i-1];
     					bubbleStop[i] = bubbleStop[i-1];	
     				}
     				nBubbles++;
     				bubbleStart[0] = 0.0;
     				bubbleStop[0]  = nodeN->nodeFrac*lenTL;
     			}
		}
	}
	else{
		bubbleStart[0] = 0.0;
		bubbleStop[0]  = nodeN->nodeFrac*lenTL;
		nBubbles++;
	}
	for (int i = 0; i<nBubbles; i++){
		bubbleStart[i] = clamp(bubbleStart[i],0,1);
		bubbleStop[i] = clamp(bubbleStop[i],0,1);
	}
	if(index == 2){
	//cout<<" bubbleStop after distribute "<<nodeN->nodeFrac<<endl;}
	}
}


void Pipe::Swap(){
	//Swap the direction of pipe and bubbles in it so that flow always goes from 0 to 1.
	node* tempNode = nodeP;
	nodeP = nodeN;
	nodeN = tempNode;
	double temp;
	//Pc = -Pc;
	for(int i = 0; i < ceil((double)nBubbles/2.0); i++){
		temp = bubbleStart[nBubbles-i-1];
		bubbleStart[nBubbles-i-1] =1.0 - bubbleStop[i];
		bubbleStop[i] =1.0 - temp;
		
		if(nBubbles-i-1 != i){
		temp = bubbleStop[nBubbles-i-1];
		bubbleStop[nBubbles-i-1] =1.0- bubbleStart[i];
		bubbleStart[i] =1.0 - temp;}
	}
	

	swapped = !swapped;
}

void Pipe::killBubbles(){
	for (int i = nBubbles-1; i >= 0; i--){
		if(abs(bubbleStop[i]-bubbleStart[i]) < 1.0e-12){
			for(int j = i; j < nBubbles; j++){
				bubbleStart[j] = bubbleStart[j+1];
				bubbleStop[j] = bubbleStop[j+1];
			}
			nBubbles--;
		}
	}
}

void Pipe::getFlow(double pressure){
	double temp_flow = mobility*(nodeP->pressure - nodeN->pressure + (swapped ? -1 : 1)*((boundary? pressure : 0) - Pc));
	
	flow = temp_flow * (swapped ? -1 : 1);
	if(temp_flow < 0){
		//switch direction
		//flow = -flow;
		Swap();
	}
}

void Pipe::calcPcandMobility(){
	double zeroZoneSize = 0.1;
	Pc = 0;
	for (int i = 0; i<nBubbles; i++){
		double str = clamp(unlerp(bubbleStart[i], zeroZoneSize, 1.0 - zeroZoneSize), 0.0, 1.0);
            	double end = clamp(unlerp(bubbleStop[i], zeroZoneSize, 1.0 - zeroZoneSize), 0.0, 1.0);
            	Pc += ((1.0 - cos(2.0*PI*str)) - (1.0 - cos(2.0*PI*end)));
	}
	Pc *= 2.0*sTens/radius;
	Pc *= (swapped?-1:1);
            	//cout<<Pc<<"    ";
	saturation = calcLinkSaturation();
        mobility = PI * pow(radius, 4) / (8.0 * (muNON * saturation + muWET * (1.0 - saturation)) * length);
}

double Pipe::calcLinkSaturation() {
    double sat = 0.0;
    for (int i = 0; i < nBubbles; i++) {
        sat += bubbleStop[i] - bubbleStart[i];
    }
    return sat;
}

Pipe::~Pipe(){
	//delete[] bubbleStart;
	//delete[] bubbleStop;
}
