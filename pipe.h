#pragma once
#include <armadillo>
#include "node.h"
#define maxBubbles 15
using namespace std;
using namespace arma;

class node;

class Pipe{
	public:
	bool swapped = false;
	int index;
	int nBubbles=0;
	double x1 = -1;
	double x2 = -1;
	double y1 = -1;
	double y2 = -1;
	double radius;
	double Pc;
	double area;
	double length = 1;
	double lenTL;
	double lenNW;
	double saturation;
	bool boundary;
	double flow = 1; //OBS!!
	double* bubbleStart;
	double* bubbleStop;
	double mobility;
	node* nodeP;
	node* nodeN;
	void setCoordinates(double x, double y);
	void moveBubble(double dt);
	void distributeBubbles();
	void Swap();
	void killBubbles();
	void getFlow(double pressure);
	double calcLinkSaturation();
	void calcPcandMobility();
	Pipe();
	~Pipe();
};
