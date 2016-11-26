#include "pipe.h"
#include <armadillo>
#define maxBubbles 15
#define SIGN (x < 0) ? -1 : (x > 0)
#define clamp(x,a,b) (x < a ? a : x > b ? b : x)
#define unlerp(t,a,b)  ( ((t) - (a)) / (double) ((b) - (a)) ) // inverse linear interpolation
#define r(x) (2*radius/(1 -0.99*cos(2*PI*(x))))  
#define alpha(x) atan((-2*PI*0.99*radius*sin(2*PI*(x))/(pow((1-0.99*cos(2*PI*(x))),2))))
#define r2(x) (radius+cos(2*PI*x))  
#define r3(X) (radius-0.0482+(0.5-sqrt(0.25-pow((X-0.5),2))))
#define alpha2(x) (atan(-2*PI*sin(2*PI*x)))
#define alpha3(X) (atan(-(1.0-2.0*X)/sqrt((1.0-X)*X)/2.0))
#define TOL 1.0e-12
#define maxBubbles 15
using namespace std;
using namespace arma;


double PI = 4.0*atan(1.0);
double sTens = 3.0; //3.0;         // surface tension: units = dyn/mm == 10 mN/m
double muNON = 0.01;        // viscosity: units = Po == 10 Pa*s
double muWET = 0.01;
bool seeded = true;
double r_m = sqrt(3)*(0.5+0.25-0.0482)-1.0;
vec cap = {18.3785,18.1614, 17.5342,16.5598,15.3204,13.9183,12.4758,11.1307,10.0274,9.2997,9.0450};

RadiusFunc radiusFunc = REAL;
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
	else if(nodeN->nodeFrac > TOL){
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


double Pipe::PcStart(double pos,double theta, double muEff){
	double Pc_s;
	double rTreshold;
	//double rMinisci = r3(pos)/cos(theta-alpha3(pos));
	if(pos > x_t_start && pos < x_t_stop){
		Pc_s = muEff*(1/r3(pos)*cos(theta-alpha3(pos)))*((pos < TOL) ? 0 : 1);
	}else if(pos < x_t_start){
		Pc_s = muEff/r3(x_t_start)*cos(theta-alpha3(x_t_start))*((pos < TOL) ? 0 : 1);
		}
	else{
		Pc_s = muEff/r3(x_t_stop)*cos(theta-alpha3(x_t_stop))*((pos < TOL) ? 0 : 1);
		//if(Pc_s != 0){
		//cout<<Pc_s<<"   "<<radius<<"   "<<x_t_stop<<endl;}
	}
	return Pc_s;
}

double Pipe::PcStop(double pos,double theta, double muEff){
	theta = -theta;
	double Pc_s;
	double rTreshold;
	//double rMinisci = r3(pos)/cos(theta-alpha3(pos));
	if(pos > 1.0-x_t_stop && pos < 1.0-x_t_start){
		Pc_s = muEff*(1/r3(pos)*cos(theta-alpha3(pos)))*((pos > 1-TOL) ? 0 : 1);
		if(Pc_s != 0){
			//cout<<"A"<<endl;
		}
	}else if(pos < 1.0-x_t_stop){
		Pc_s = muEff/r3(1.0-x_t_stop)*cos(theta-(atan(-(1.0-2.0*(1.0-x_t_stop))/sqrt((1.0-(1.0-x_t_stop))*(1.0-x_t_stop))/2.0)))*((pos > 1 - TOL) ? 0 : 1);
		if(Pc_s != 0){
			//cout<<"B"<<(1.0-x_t_stop)<<alpha3(1.0-x_t_stop)<<endl;
			//cout<< (atan(-(1.0-2.0*(1.0-x_t_stop))/sqrt((1.0-(1.0-x_t_stop))*(1.0-x_t_stop))/2.0));
		}
		}
	else{
		Pc_s = muEff/r3(1.0-x_t_start)*cos(theta-alpha3(1.0-x_t_start))*((pos > 1 - TOL) ? 0 : 1);
		if(Pc_s != 0){
			//cout<<"C"<<endl;
		}
	}
	theta = -theta;

	return Pc_s;
	
}

void Pipe::calcPcandMobility(){
	double zeroZoneSize = 0.0;
	double muEff;
	Pc = 0;
	for (int i = 0; i<nBubbles; i++){
		double str = clamp(unlerp(bubbleStart[i], zeroZoneSize, 1.0 - zeroZoneSize), 0.0, 1.0);
            	double end = clamp(unlerp(bubbleStop[i], zeroZoneSize, 1.0 - zeroZoneSize), 0.0, 1.0);
           	//double str = bubbleStart[i];
           	//double end = bubbleStop[i];
           
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
         		//Pc += muEff*(1/r3(str)*cos(theta-alpha3(str))*((str < TOL) ? 0 : 1)-1/r3(end)*cos(-theta-alpha3(end))*((end > 1-TOL) ? 0 : 1));
         		Pc += PcStart(str,theta,muEff) - PcStop(end,theta,muEff);
         		//Pc += muEff*Pc_max*((str < 0.1) ? 0 : 1)-muEff*Pc_max*((end > 0.9)? 0 : 1);
         		//Pc += muEff*(((bubbleStart[i] > TOL ? 1 : (nodeN->connectingBubbles ? 0 : 1) )*1/r3(str)*cos(theta-alpha3(str)))-(bubbleStop[i] < 1 - TOL ? 1 : (nodeP->connectingBubbles ? 0 : 1) )*1/r3(end)*cos(-theta-alpha3(end)));
         		//Pc += (str == 0 ? 0 : str == 1? 0: cap(i)*(1+(0.2*-radius)/0.2*(1-0.57)))+(end = 0 ? 0 : end = 1? 0: cap(i)*(1+(0.2*-radius)/0.2*(1-0.57)));
    			//cout<<muEff*Pc_max*((str < 0.3) ? 0 : 1)<<"  "<<muEff*Pc_max*((end > 0.7) ? 0 : 1)<<endl;
    		
            	}
            	if(Pc != 0){
            		//cout<<nodeN->connectingBubbles<<nodeP->connectingBubbles;
            		//cout<<str<<"     "<<end<<endl;
            	}
   
		if(Pc != 0){
         	//cout<<PcStart(str,theta,muEff)<<"   "<<PcStop(end,theta,muEff)<<endl;
         }
	}
	Pc *= (swapped?1:-1);
	
	saturation = calcLinkSaturation();
	
        if(tubeType == CIRC){
		mobility = PI * pow(radius, 4) / (8.0 * (muNON * saturation + muWET * (1.0 - saturation)) * length);
	}
	else{
		mobility = 3.0/5.0*G*area*area / ((muNON * saturation + muWET * (1.0 - saturation)) * length);
		
	}
}



void Pipe::getCapBoundaries(){
	double N = 1000.0;
	double l_c;
	double Alpha;
	x_t_start = 0.0;
	x_t_stop = 1.0;
	for (double x = 1.0/2.0 ;x < 1.0; x += 1.0/N){
		Alpha = (atan(-(1.0-2.0*x)/sqrt((1.0-x)*x)/2.0));
		l_c = r3(x)*(1-sin(theta-Alpha))/cos(theta-Alpha);
		//cout<<"l_c "<<l_c<<endl;
		if(x+l_c > 1+r_m){
			x_t_stop = x;
			break;
		}	
	} 
	for (double x = 1.0/2.0 ;x > 0.0; x -= 1.0/N){
		Alpha = (atan(-(1.0-2.0*x)/sqrt((1.0-x)*x)/2.0));
		l_c = r3(x)*(1-sin(theta-Alpha))/cos(theta-Alpha);
		if(x+l_c < -r_m){
			x_t_start = x;
			break;
		}	
	} 
}

void Pipe::getPcMax(){
	double N = 1000.0;
	Pc_max = -1000.0;
	for(double x = 0.0; x<1.0; x += 1.0/N){
		if(Pc_max<(1/r3(x)*cos(theta-alpha3(x)))){
			Pc_max = (1/r3(x)*cos(theta-alpha3(x)));
		}
	}
	//cout<<Pc_max<<endl;
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
