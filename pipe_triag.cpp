#include "pipe.h"
#include <armadillo>
#define maxBubbles 15
#define SIGN (x < 0) ? -1 : (x > 0)
#define clamp(x,a,b) (x < a ? a : x > b ? b : x)
#define unlerp(t,a,b)  ( ((t) - (a)) / (double) ((b) - (a)) ) // inverse linear interpolation
#define r(x) (2*radius/(1 -0.99*cos(2*PI*(x))))  
#define alpha(x) atan((-2*PI*0.99*radius*sin(2*PI*(x))/(pow((1-0.99*cos(2*PI*(x))),2))))
#define r2(x) (radius+cos(2*PI*x))  
#define r3(x) (radius+(0.5-sqrt(0.25-pow((x-0.5),2))))
#define alpha2(x) (atan(-2*PI*sin(2*PI*x)))
#define alpha3(x) (atan(-(1-2*x)/sqrt((1-x)*x)/2))
#define TOL 1.0e-12
#define maxBubbles 15
using namespace std;
using namespace arma;


double PI = 4.0*atan(1.0);
double sTens = 7.0; //3.0;         // surface tension: units = dyn/mm == 10 mN/m
double muNON = 0.01;        // viscosity: units = Po == 10 Pa*s
double muWET = 0.01;
bool seeded = true;
vec cap = {18.3785,18.1614, 17.5342,16.5598,15.3204,13.9183,12.4758,11.1307,10.0274,9.2997,9.0450};

RadiusFunc radiusFunc = UNREAL;
TubeType tubeType = TRIAG1;



Pipe::Pipe(){
	bubbleStart = new double[maxBubbles];
	bubbleStop  = new double[maxBubbles];
	for ( int i = 0; i<maxBubbles; i++){
		bubbleStart[i] = 0;
		bubbleStop[i] = 0;
	}
	if(!seeded){
	srand (time(NULL));
	seeded = true;
	}
	boundary = 0;
	flow = 1;
	//G = 1;
	if(tubeType == CIRC){
		radius = ((double)rand()/RAND_MAX*3.0+1.0)*0.1;
		area = PI*pow(radius,2);
	}
	else if(tubeType == TRIAG1){
		G = 0.04811;
		radius = ((double)rand()/RAND_MAX*3.0+1.0)*0.1*pow(10.0/3.0*PI*G,0.25);
		area = radius*radius/(G*4);
	}
	else{
		G = 0.03944;
		radius = ((double)rand()/RAND_MAX*3.0+1.0)*0.1*pow(10.0/3.0*PI*G,0.25);
		area = radius*radius/(G*4);
	}
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
	double start 	= 1.0;
        double end 	= 1.0 +lenTL;
        
        for (int i = 0; i < nBubbles; i++) {
            lenNW += clamp(bubbleStop[i], start, end) - clamp(bubbleStart[i], start, end);
        }
        if(lenTL+10e-12<lenNW &&nBubbles == 1){
        cout<<"Alert   "<<"nbubbles "<<nBubbles<<endl;
        cout<<"lenTL: "<<lenTL<<"lenNW: "<<lenNW<<endl;
        }
        nodeP->nodeFluxTL += lenTL*area;
        nodeP->nodeFluxNW += lenNW*area;	
}

void Pipe::distributeBubbles(){
	double dx1 = bubbleStart[0] - lenTL;                  // size of first W  bubble
        double dx2 = bubbleStop[0]  - bubbleStart[0];    // size of first NW bubble
     	if(nBubbles > 0){
     		if(dx1 < TOL){
     			if(dx2 < radius){
     				//nW first
     				bubbleStart[0] -= lenTL;
     				bubbleStop[0]  -= (1.0 - nodeN->nodeFrac)*lenTL;
     				//cout<<"a";
     			}
     			else{
     				//W first
     				bubbleStart[0] -= nodeN->nodeFrac*lenTL;
     				//cout<<"b";
     			}
     		}
     		else{
     			if(dx1 < radius){
     				// putW first
     				bubbleStart[0] -= nodeN->nodeFrac*lenTL;
     				//cout<<"c";
     			}
     			else{
     				//put nW first
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
		Swap();
	}
}

void Pipe::calcPcandMobility(){
	double zeroZoneSize = 0.1;
	double muEff;
	Pc = 0;
	for (int i = 0; i<nBubbles; i++){
		double str = clamp(unlerp(bubbleStart[i], zeroZoneSize, 1.0 - zeroZoneSize), 0.0, 1.0);
            	double end = clamp(unlerp(bubbleStop[i], zeroZoneSize, 1.0 - zeroZoneSize), 0.0, 1.0);
            	
            	if(tubeType == CIRC){
            		muEff = 2.0*sTens;
            	}else{
            		muEff = (1.0 +2.0*sqrt(PI*G))*sTens;
            	}
            	
            	if(radiusFunc == OLD){
            	Pc += cos(theta)*((1.0 - cos(2.0*PI*str)) - (1.0 - cos(2.0*PI*end)))*muEff/radius;}
            	else if(radiusFunc == UNREAL){
            		Pc += muEff*(1/r(str)*cos(theta-alpha(str))-1/r(end)*cos(-theta-alpha(end)));
            	}
            	else if(radiusFunc == REAL){	
            		
         		//Pc += muEff*(1/r2(str)*cos(theta-alpha2(str))-1/r2(end)*cos(-theta-alpha2(end)));
         		Pc += muEff*(1/r3(str)*cos(theta-alpha3(str))*((str==1 || str == 0) ? 0 : 1)-1/r3(end)*cos(-theta-alpha3(end))*((end==1 || end == 0) ? 0 : 1));
         		//Pc += (str == 0 ? 0 : str == 1? 0: cap(i)*(1+(0.2*-radius)/0.2*(1-0.57)))+(end = 0 ? 0 : end = 1? 0: cap(i)*(1+(0.2*-radius)/0.2*(1-0.57)));
    			/*if(true){
    			cout<<"Pc 1: " <<muEff*(1/r2(str)*cos(theta-alpha2(str))-1/r2(end)*cos(-theta-alpha2(end)))<<"Pc 2: "<< muEff*(1/r3(str)*cos(theta-alpha3(str))*((str==1 || str == 0) ? 0 : 1)-1/r3(end)*cos(-theta-alpha3(end))*((end==1 || end == 0) ? 0 : 1))<<"Pc 3: "<<(str == 0 ? 0 : str == 1? 0: cap(i)*(1+(0.2*-radius)/0.2*(1-0.57)))+(end = 0 ? 0 : end = 1? 0: cap(i)*(1+(0.2*-radius)/0.2*(1-0.57)))<<endl;
    			}*/
            	}
            	//if(Pc != 0)
            	//cout<<"Pc: "<<Pc;
	}
	Pc *= (swapped?-1:1);
	saturation = calcLinkSaturation();
	
        if(tubeType == CIRC){
		mobility = PI * pow(radius, 4) / (8.0 * (muNON * saturation + muWET * (1.0 - saturation)) * length);
	}
	else{
		mobility = 3.0/5.0*G*area*area / ((muNON * saturation + muWET * (1.0 - saturation)) * length);
		
	}
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
