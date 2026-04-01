#pragma once
#ifndef POINT2D_H
#define POINT2D_H
#include"Parameter.h"

struct Point2D
{
	REAL x;
	REAL y;
};
struct Position2D
{
	int x;
	int y;
};
struct Force
{
	REAL x;
	REAL y;
};

struct Node2D
{
	int dup;
	//double dup;
	double rho;
	Point2D u;
	double p;
};

struct Commpair
{
	int localLatt;
	int buffLatt;
};

struct Commpair_adv
{
	int localLatt;
	int buffLatt;
	int fainfo;
	int soninfo;
};

#endif