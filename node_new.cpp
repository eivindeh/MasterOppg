#include "node.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <armadillo>
#include <time.h>
#define SIGN(x) (x < 0) ? -1 : (x > 0)
#define TOL 1.0e-12

using namespace std;
using namespace arma;

void node::setConnections(Pipe pipes[],node nodes[],int N,int M){

	n = (index)%(N);
	m = floor((index)/(N));

	x = (n+1)*sqrt(3)/2;
	y = ((m%2)?(n%2)*0.5:(1-(n%2))*0.5)+1.5*m;
	
	leftPipe  = &pipes[(m*N*3/2+n)-1+(n==0)*(N)];
	leftNode  = &nodes[index-1+(n==0)*(N)];
	
	rightPipe = &pipes[m*N*3/2+n];
	rightNode = &nodes[index+1-(n==(N-1))*(N)];
	
	double above = m*N*3/2+N+(double)(n+1+((m+1)%(2)))/2-1;
	double below = m*N*3/2-(double)(N-n+((m+1)%(2)))/2+(m==0)*M*N*3/2;
	if (ceil(above) == above){
		horizontalPipe = &pipes[(int)above];
		horizontalNode = &nodes[index+N];
		upper = true;
		if (m == (M-1)){
			horizontalNode = &nodes[index+N-(m==(M-1))*N*M];
			horizontalPipe->boundary = 1; 
		}
	}
	else{
		horizontalPipe = &pipes[(int)below];
		horizontalNode = &nodes[index-N];
		upper = false;
		if(m==0){
			horizontalNode = &nodes[index-N+N*M];
			horizontalPipe->boundary = 1;
		}
	}
	
	if(upper){
		leftPipe->nodeN = this;
		rightPipe->nodeN = this;
		horizontalPipe->nodeN = this;
	}else{
		leftPipe->nodeP = this;
		rightPipe->nodeP = this;
		horizontalPipe->nodeP = this;
	}
	
	leftPipe->setCoordinates(x,y);
	rightPipe->setCoordinates(x,y);
	horizontalPipe->setCoordinates(x,y);
}

void node::updateNodeFrac(){
	nodeFrac = nodeFluxNW/(nodeFluxTL+10e-50);
	nodeFluxNW = 0.0;
	nodeFluxTL = 0.0;
}

void node::bubbleAdjust(Pipe* pipes) {
	double dx1s[3] = { 1.0, 1.0, 1.0 }; // sizes of first W  bubble in each neighouring tube
        double dx2s[3] = { 0.0, 0.0, 0.0 }; // sizes of first NW bubble
        int typeFirst[3] = { 0,0,0 }; // 1 == NW, 0 == W
        int links[3];
        int link;
        int linksId[3] = { 0,0,0 };
        int ctr = 0;
        bool isSmallW = false;
        bool isSmallNW = false;
        for(int i = 0; i < 3; i++){
        	if(i == 0){
        		link = leftPipe->index;
        	}else if (i == 1){
        		link = rightPipe->index;
        	}
        	else{
        		link = horizontalPipe->index;
        	}
        	if (pipes[link].nBubbles == 0){
        		//cout<<"skjerdette? "<<endl;
        		continue;
        	}
        	int bubble = (pipes[link].nodeN->index == index ? 0 : pipes[link].nBubbles-1);
        	dx1s[i] = (pipes[link].nodeN->index == index ? pipes[link].bubbleStart[bubble] : 1.0 - pipes[link].bubbleStop[bubble]); //switched
        	dx2s[i] = pipes[link].bubbleStop[bubble] - pipes[link].bubbleStart[bubble];
        	
        	if(dx1s[i] < TOL){
        		if(dx2s[i] < pipes[link].radius && dx2s[i]> TOL){
        			typeFirst[ctr] = 1;
        			links[ctr] = link;
        			linksId[ctr++] = i;
        			isSmallNW = true;
        		}
        		}
        	else{
        		if (dx1s[i] <pipes[link].radius) {
        			links[ctr] = link;
        			linksId[ctr++] = i;
        			isSmallW = true;
        		}
        	}
        }
        
        
            if (isSmallW && isSmallNW) {
            double sumW = 0.0;
            double sumNW = 0.0;
            for (int i = 0; i < ctr; i++) {
                if (typeFirst[i]) {
                    sumNW += dx2s[linksId[i]]*pipes[links[i]].area;
                } else {
                    sumW += dx1s[linksId[i]]*pipes[links[i]].area;
                }
            }
            double rem = 0.0;
            double add = 0.0;
            double minVol = min(sumNW, sumW);

            for (int i = 0; i < ctr; i++) {
                int link = links[i];
                if (typeFirst[i]) {
                    // remove NW from tubes with small NW
                    if (pipes[link].nodeP->index == index) {
                        pipes[link].bubbleStart[pipes[link].nBubbles-1] += dx2s[linksId[i]]*minVol/sumNW;
                        rem +=  (dx2s[linksId[i]]*minVol/sumNW)*pipes[link].area;
                    } else {
                        pipes[link].bubbleStop[0] -= dx2s[linksId[i]]*minVol/sumNW;
                        rem +=  (dx2s[linksId[i]]*minVol/sumNW)*pipes[link].area;
                    }                     
                } else {
                    // add NW to tubes with small W
                    if (pipes[link].nodeP->index == index) {
                        pipes[link].bubbleStop[pipes[link].nBubbles-1] += dx1s[linksId[i]] * minVol/sumW;
                        add +=  (dx1s[linksId[i]] * minVol/sumW)*pipes[link].area;
                    } else {
                        pipes[link].bubbleStart[0] -= dx1s[linksId[i]] * minVol/sumW;
                        add +=  (dx1s[linksId[i]] * minVol/sumW)*pipes[link].area;
                    }     
                }
            }
        } else if (isSmallNW && ctr == 1) {
           // double lenTL = pipes[links[0]].flow * deltaT / pipes[links[0]].area;
            int bubble = (pipes[links[0]].nodeP->index == index ? pipes[links[0]].nBubbles - 1 : 0);
            if (nodeFrac < TOL) {
                pipes[links[0]].bubbleStart[bubble] += pipes[links[0]].lenTL;
                pipes[links[0]].bubbleStop[bubble] += pipes[links[0]].lenTL;
            }
        } else if (ctr == 1 && isSmallW) {
            // Add option to move small W bubbles also
            // but maybe not with full linkQ?, due to less surface area?
            /*
                For example:
                    double alpha = 0.5;
                    double lenTL = linkQ[links[0]]*deltaT/linkVolume[links[0]];
                    int bubble = (node == linkNodeUp[links[0]] ? linkNumBubbles[links[0]]-1 : 0);
                    if (nodeFrac[node] > 1.0 - TOL) {
                        bubbleStart[links[0]][bubble] += alpha*lenTL;
                        bubbleStop[links[0]][bubble] += alpha*lenTL;
                    }
            */
        } else if (isSmallW xor isSmallNW) {
            // if two small bubbles, handle?
            // NOTE: If typeFirst[0] == 1 && typeFirst[1] == 1, move to one side?
        }
        
        
        
/*
    	
    	int linkNW = (typeFirst[0] ? links[0] : links[1]);
        int linkW = (typeFirst[0] ? links[1] : links[0]);
        
        double volNW = dx2s[(typeFirst[0] ? linksId[0] : linksId[1])] * pipes[linkNW].area;
        double volW  = dx1s[(typeFirst[0] ? linksId[1] : linksId[0])] * pipes[linkW].area;
        
        if (pipes[linkNW].nodeP->index == index)
        	pipes[linkNW].bubbleStart[pipes[linkNW].nBubbles-1] += min(volNW,volW)/pipes[linkNW].area;
        else
        	pipes[linkNW].bubbleStop[0] -= min(volNW,volW) / pipes[linkNW].area;
        
        if(pipes[linkW].nodeP->index == index) //switched to P from N.
        	pipes[linkW].bubbleStop[pipes[linkW].nBubbles-1] += min(volNW,volW)/pipes[linkW].area;
        else
        	pipes[linkW].bubbleStart[0] -= min(volNW,volW) / pipes[linkW].area;
        
        }*/
        //delete[] links;
}
/*void node::fillAandB(mat *A, vec *B){
	if(index == 0){
		(*A)(index,index) = 1;
	}
	else{
		(*A)(index,index) 			= -rightPipe->mobility - leftPipe->mobility - horizontalPipe->mobility;
		(*A)(index,horizontalNode->index) 	= horizontalPipe->mobility;
		(*A)(index,leftNode->index) 		= leftPipe->mobility;
		(*A)(index,rightNode->index) 		= rightPipe->mobility;
		(*B)(index) 				= horizontalPipe->mobility*(horizontalPipe->boundary*(upper?-1:1)+ (upper?-1:1)*3.0) + 0*(
		rightPipe->mobility*rightPipe->Pc + leftPipe->mobility*leftPipe->Pc + horizontalPipe->mobility*horizontalPipe->Pc);
	}
}*/



/*void node::getFlow(){
	double leftFlow = -leftPipe->mobility*(leftNode->pressure-pressure)*(upper ? -1 : 1);
	double rightFlow = rightPipe->mobility*(pressure-rightNode->pressure)*(upper ? -1 : 1);
	double horizontalFlow = horizontalPipe->mobility*(pressure - horizontalNode->pressure + (upper ? -horizontalPipe->boundary : horizontalPipe->boundary))*(upper ? -1 : 1);
	if(index ==1){
		cout<<endl<<endl<<"pipeI " <<rightFlow<<" pipeJ "<<rightPipe->flow<<endl<<endl;
	}
	
	if(SIGN(leftPipe->flow) != SIGN(leftFlow)){
		leftPipe->Swap(); //???
	}
	if(SIGN(rightPipe->flow) != SIGN(rightFlow)){
		rightPipe->Swap(); //???
	}
	if(SIGN(horizontalPipe->flow) != SIGN(horizontalFlow)){
		horizontalPipe->Swap(); //???
	}
	leftPipe->flow = leftFlow;
	rightPipe->flow = rightFlow;
	horizontalPipe->flow = horizontalFlow;
	//cout<<"horizontal flow "<< horizontalPipe->flow<<endl;
}*/
