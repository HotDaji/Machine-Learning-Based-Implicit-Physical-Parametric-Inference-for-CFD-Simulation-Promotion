#pragma once
#ifndef MESH_H
#define MESH_H
#define _USE_MATH_DEFINES
#include<vector>
#include <cmath> // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
#include <cstdlib> // standard library
#include <string>
#include <cstdio>
#include <algorithm>
#include <unordered_map>//for mangement based on Xflow
#include<mpi.h>
#include"NodeAndLattice.h"
#include <torch/script.h>//libtorch
#include <pthread.h>//绑定torch到目标cpu的函数
#include <ATen/Parallel.h>  // 包含设置线程数的正确头文件
struct PointData {//"VelocityModule","StaticPressure","Vort","TotalPressure","Turbulence","Velocity:0","Velocity:1","Velocity:2","Vertex:0","Vertex:1","Vertex:2"
    double VelocityModule;
    double StaticPress;
   
    double TotalPressure;
	double viscosity;
    double ux;
    double uy;
    double uz;
    double vx;
    double vy;
    double vz;
    };
class Mesh{
public:
	Mesh();
	~Mesh();
	vector<Particlenode> IBParticle;
	vector<vector<Lattice>> LatMesh;
	int LatmeshSizeNx;
	int LatmeshSizeNy;
	Point2D center;
	double radius;
	const double eps = 1.0e-9;
    double startTau;
	Point2D startU;
	int MyrankID;//本对象运行的cpu核心id
	int MyCoreNum;//总共的运行cpu核心数量
	int stepH,stepL;
	int startAvg=0;;
	double StartCs;
	int errorflag = 0;
	int ParallelSize=0;
	
public:
#ifdef D2Q9
	void InitialLatMesh(int Nx_P, int Ny, double tau, Point2D StartVel, double Rho , int RankID , int CoreNum, int PartitionRemainderFlag);
	void InitialIBParticle(int NumNodes, double Radius, double ParticlePositionX, double ParticlePositionY , Point2D ParticleVel);
	void Streaming(int RankID, int CoreNum);
	void BoundaryCondition(int RankID, int CoreNum);
	void OneLBMStep(int RankID, int CoreNum,REAL Cs,REAL Pr,int t);//
	void CalParticleForce(int Radius);//
	void UpdateVelocity(int RankID, int CoreNum);
	void InterpolateVel(int RankID, int CoreNum, int Nx_P);//
	void SpreadForce(int RankID, int CoreNum, int Nx_P);
	void CaculateLESTau(int i, int j,REAL Cs, REAL Pr,int RankID,int CoreNum);//通过S模型更新tau
	REAL DiracDeltaD2Q9(REAL x, REAL y);
	void UpdateBoundary();
	double ComputeDistance( double x0, double y0, double x1, double y1 );
#endif // D2Q9
public:
	void JudgeStyle(Point2D CylinderCenter1, REAL Radius, int RankID , int CoreNum);
	void InitialFile();
	void WriteDatFluid(int t , int RankID , int NX , int CoreNum);
	void MergeDatFluid(int t, int CoreNum);
	void WriteDatParticle(int t);
	REAL DiracDelta(REAL x);
	REAL GetDis2D(REAL x1, REAL y1, REAL x2, REAL y2);
	void OutputCd(int RankID,int CoreNum, int time , double Radius);
	void OutputCp(int RankID,int CoreNum, int time , double Radius);
	void MergeCp(int t, int CoreNum);
	void ComputeCp(int startcp,int t);
	void ComputeAvgVelocity(int t);
	void WriteAvgu(int t , int RankID , int NX , int CoreNum);
	void MergeAvgu(int t, int CoreNum);
	void SetNormalDis();
	void ErrorCatch(int i , int j);
	
//------------math--------------------//
	double GetDuxDx(int i , int j);
	double GetDuxDy(int i , int j);
	double GetDuyDx(int i , int j);
	double GetDuyDy(int i , int j);
	void GetSt(int t);
	double ComputeDynamicCs(int i, int j);
	void FilterLBM(int i, int j);//滤波操作
	bool isAvailable(int i , int j);

	vector<PointData> ReadXflowdata(int t);
	void ManageDataFromXflow(int t);
	void XflowBoundary(int t);
	void XflowABS(int t);
	
	void FindAndWriteTau(int i, int j,int t);
	void ChangeTauFromXflow(int t);
	//void CalInletMSE(int t);
	int NNmodel_Type01(int t);
	int NNmodel_Type05(int t);
	int NNmodel_Type03(int t);
	int NNmodel_Type07(int t);
	int NNmodel_Type00(int t);
	int NNmodel_Type02(int t);
	int NNmodel_Type04(int t);
	int NNmodel_Type06(int t);
	
};
#endif