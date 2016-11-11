#pragma once
#include <armadillo>
#include "node.h"
#define maxBubbles 15
using namespace std;
using namespace arma;


class node;

enum RadiusFunc {OLD, UNREAL, REAL};
enum TubeType {CIRC, TRIAG1, TRIAG2};

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
	double G;
	double saturation;
	bool boundary;
	double flow = 1; //OBS!!
	double* bubbleStart;
	double* bubbleStop;
	double mobility;
	double theta = 0;
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
	//RadiusFunc radiusFunc;
	//TubeType tubeType;

	Pipe();
	~Pipe();
};
